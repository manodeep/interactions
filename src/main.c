#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
/* #include <signal.h> */

/* For file open and other file handling routines*/
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "macros.h"
#include "utils.h"

#include "interactions.h"
#include "tree_conversion.h"
#include "progressbar.h"

int ThisTask=0, NTasks=1;
int exit_status = EXIT_FAILURE;

void bye()
{
#ifdef MPI
    MPI_Finalize();
#endif
}


int main(int argc, char **argv)
{
    
#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTasks);
#else
    ThisTask = 0;
    NTasks = 1;
#endif

    if(argc <= 3) {
        fprintf(stderr,"\n  usage: `%s' <scale-factors> <output filebase> <input file1> ...\n", argv[0]);
        fprintf(stderr,"<scale-factors>   \t : string, filename for a file with a list of scale-factors\n");
        fprintf(stderr,"<output filebase> \t : string, output filename. Each input file will have a corresponding output file\n");
        fprintf(stderr,"<input file1>     \t : string, input files to be processed\n\n");
        bye();
        return EXIT_FAILURE;
    }
    atexit(bye);
    
    struct lh_vtree allnodes;
    /* Initialize the statically allocated pointer array.
       This way, on malloc or function failure, free(allnodes->tree[i]) functions correctly
       within the goto 'fail' label. with functions */
    {
        FILE *fp=my_fopen(argv[1],"r");
        if(fp == NULL) {
            return EXIT_FAILURE;
        }
        for(int32_t i=0;i<MAX_SNAPSHOTS;i++) {
            allnodes.tree[i] = NULL;
            double a;
            fscanf(fp, "%lf",&a);
            allnodes.redshifts[i] = 1.0/a - 1.0;//convert scale-factor to redshift
        }
        fclose(fp);
    }
    

    
    for(int iarg=3 + ThisTask; iarg < argc; iarg += NTasks){
        char outputfile[MAXLEN];
        char filename[MAXLEN];
        int status = my_snprintf(outputfile, MAXLEN, "%s.%d", argv[2], iarg-3);
        if(status == -1) {
            fprintf(stderr,"Error: (On ThisTask = %d) Could not write name of output file into variable. Increase `MAXLEN' in 'macros.h' and recompile\n",
                    ThisTask);
            perror(NULL);
            ABORT(EXIT_FAILURE);
        }
        status = my_snprintf(filename, MAXLEN, "%s",argv[iarg]);
        if(status == -1) {
            fprintf(stderr,"Error: (On ThisTask = %d) Could not write name of output file into variable. Increase `MAXLEN' in 'macros.h' and recompile\n",
                    ThisTask);
            perror(NULL);
            ABORT(EXIT_FAILURE);
        }
        fprintf(stderr,"Reading file `%s'\n", filename);

        FILE *fp = my_fopen(outputfile, "w");
        int32_t ntrees=0;
        int32_t totnhalos = 0;
        int32_t *nhalos_per_tree=NULL;
        size_t offset = 0;
        
#if SINGLE_TREE_ATATIME

        status = read_file_headers_lhalotree(filename, &ntrees, &totnhalos, &nhalos_per_tree);
        int32_t max_nhalos_in_tree = 100;
        for(int32_t i=0;i<ntrees;i++) {
            max_nhalos_in_tree = max_nhalos_in_tree < nhalos_per_tree[i]  ? nhalos_per_tree[i]:max_nhalos_in_tree;
        }
        struct lhalotree *tree = my_malloc(sizeof(*tree), max_nhalos_in_tree);
        int fd = open(filename,  O_RDONLY, O_NONBLOCK);
        if(fd < 0) {
            fprintf(stderr,"Error in opening input file\n");
            perror(NULL);
            exit_status = fd;
            goto fail;
        }
#else
        //Load in the entire tree in one shot
        struct lhalotree *all_trees = read_entire_lhalotree(argv[iarg], &ntrees, &totnhalos, &nhalos_per_tree);
#endif
        int interrupted=0;
        init_my_progressbar(ntrees, &interrupted);
        
        for(int32_t itree=0;itree<ntrees;itree++) {
            my_progressbar(itree, &interrupted);
            const int32_t nhalos = nhalos_per_tree[itree];

#ifdef SINGLE_TREE_ATATIME
            off_t file_offset = sizeof(int32_t) /* ntrees */
                + sizeof(int32_t)                           /* totnhalos */
                + sizeof(int32_t)*ntrees   /* nhalos per tree */
                + offset*sizeof(struct lhalotree);
        
            status = pread_single_lhalotree_with_offset(fd, tree, nhalos, file_offset);
            if(status != EXIT_SUCCESS) {
                fprintf(stderr,"Error: (On ThisTask = %d) could not read treenum %d from disk. expected nhalos = %d at offset (nhalos) = %zu offset(bytes) = %zu\n",
                        ThisTask, itree, nhalos, offset, file_offset);
                goto fail;
            }

#else
            struct lhalotree *tree = &(all_trees[offset]);
#endif

            
            /* Validate tree */
            for(int64_t i=0;i<nhalos;i++) {
                if(tree[i].NextProgenitor == -1) continue;
                XRETURN(tree[i].NextProgenitor >=0 && tree[i].NextProgenitor < nhalos, EXIT_FAILURE,
                        "halo = %"PRId64", NextProgenitor = %d must be within [0, %d)\n", i, tree[i].NextProgenitor, nhalos);
            }
            
            
            /* Un-fix flybys. This step might be redundant since the standard lhalotree output
               does not contain flybys. lhalotree files generated by https://github.com/manodeep/ConvertCTrees
               has negative MostBoundID within the group to denote a flyby.
            */

#ifdef FIX_FLYBYS
            
            /*(Un)Fix flybys */
            for(int32_t i=0;i<nhalos;i++){
                if(tree[i].MostBoundID < 0) {
                    /* This is a flyby FOF halo */
                    tree[i].MostBoundID = -tree[i].MostBoundID;
                    int32_t curr_fof_index = tree[i].FirstHaloInFOFgroup;
                    tree[i].FirstHaloInFOFgroup = i;
                    /* All of the original satellites of this FOF would point to the stitched FOF as the host FOF.
                       Move them back and make them point to this FOF */
                    int32_t next = tree[i].NextHaloInFOFgroup;
                    while(next != -1) {
                        tree[next].FirstHaloInFOFgroup = i;
                        if(next == tree[next].NextHaloInFOFgroup) {
                            fprintf(stderr,"WARNING: Loop detected next = %d tree[next].NextHaloInFOFgroup = %d i=%d\n",
                                    next, tree[next].NextHaloInFOFgroup, i);
                            tree[next].NextHaloInFOFgroup = -1;
                        } else {
                            next = tree[next].NextHaloInFOFgroup;
                        }
                    }

                    /* Some other (sub)halo must point to this i'th halo as the "NextHaloInFOFgroup" -> fix that */
                    int32_t curr = curr_fof_index;
                    while(curr != -1 && tree[curr].NextHaloInFOFgroup != i) {
                        curr = tree[curr].NextHaloInFOFgroup;
                    }
                    if(curr == -1) {
                        fprintf(stderr,"Error: (On ThisTask = %d) While fixing flyby halo, could not locate the previous subhalo pointing to this FOF. "
                                "Treenum = %d (file = `%s') forced FOF = %d (flyby) halo index=%d \n",
                                ThisTask, itree, filename, curr_fof_index, i);
                        fprintf(stderr,"Please file bug report with ConvertCTrees code here: https://github.com/manodeep/ConvertCTrees/issues \n");
                        fprintf(stderr,"Here is some debug info for the bug report\n");
                        fprintf(stderr,"*************************************\n");
                        curr = curr_fof_index;
                        fprintf(stderr,"\n## Properties of previous Host FOF halo \n\n");
                        while(curr != -1) {
                            fprintf(stderr,"Host halo index = %9d len = %8d MboundID = %12lld\n",curr, tree[curr].Len, tree[curr].MostBoundID);
                            curr = tree[curr].NextHaloInFOFgroup;
                        }

                        fprintf(stderr,"\n## Properties of current host FOF halo \n\n");
                        curr = i;
                        while(curr != -1) {
                            fprintf(stderr,"Host halo index = %9d len = %8d MboundID = %12lld\n",curr, tree[curr].Len, tree[curr].MostBoundID);
                            curr = tree[curr].NextHaloInFOFgroup;
                        }
                        fprintf(stderr,"*************************************\n");
                    } else {
                        tree[curr].NextHaloInFOFgroup = -1;
                    }
                }
            }
#endif //FIX_FLYBYS            

#if 0
            status = write_interactions_to_stream(tree, nhalos,fp);
            if(status != EXIT_SUCCESS) {
                fprintf(stderr,"Error: (On ThisTask = %d) could not write interactions to file\n", ThisTask);
                goto fail;
            }
#endif

            status = convert_tree_to_hinge(tree, nhalos, LHALOTREE, &allnodes);
            if(status != EXIT_SUCCESS) {
                goto fail;
            }
            
            /*assign haloids*/
            status = assign_haloid(&allnodes);
            if(status != EXIT_SUCCESS) {
                goto fail;
            }

            
            /* The trees are now in the locally-horizontal vertical tree format -> record interactions */
            int64_t nhalos_so_far=0;
            for(int isnap=allnodes.max_snap;isnap>=allnodes.min_snap;isnap--) {
                struct node_data *BaseNode = allnodes.tree[isnap];//Get the starting halo at this snapshot
                const int64_t nhalos_snap = allnodes.Nhalos_per_snap[isnap];
                nhalos_so_far += nhalos_snap;
                for(int64_t igroup=0;igroup<nhalos_snap;igroup++) {
                    struct node_data *halo = &(BaseNode[igroup]);//This the unique halo given by (isnap, igroup)
                    assert(halo->snapshot == isnap && halo->halonum == igroup);
                    assert(halo->FofHalo != NULL);
                    fprintf(fp,"%3d  %12"PRId64"  %2d  %12"PRId64"  %12.4e  %8d  %12.6e %2d\n",
                            isnap,    igroup, halo->isFOF, halo->haloid, halo->Mtot*1e10,
                            halo->FofHalo->halonum, halo->FofHalo->Mtot*1e10, halo->FofHalo->isFOF);
                }

            }
            assert(nhalos_so_far == nhalos);
            
            
            
            /* Free up the locally horizontal vertical tree */
            free_lh_vtree(&allnodes);
            offset += nhalos;
        }//loop over ntrees in a single file
        finish_myprogressbar(&interrupted);
    fail:
#ifdef SINGLE_TREE_ATATIME            
        close(fd);
        free(tree);
#else
        free(all_trees);
#endif            
        free(nhalos_per_tree);
        fclose(fp);
        
        if(status != EXIT_SUCCESS) {
            ABORT(status);
        }
    }//loop over all files


    exit_status = 0;
    return 0;
}


#if 0   //#if 0 -> Code completely commented out since the trees are now converted into HINGE style decomposition into (snapshot, nhalos) + pointers for navigation         

            //Record every interaction that begins with a new subhalo appearing
            for(int32_t i=0;i<nhalos;i++) {

                /* Various kinds of interactions where a new subhalo appears within a FOF -> as defined by HINGE. The interactions
                   themselves are defined based on what happens in the future and noted at the snapshot where the new subhalo appeared.
                   
                   0      -> merger, subhalo dissolves within a FOF halo. The standard LCDM merger idea we have
                   2      -> merger, subhalo remains intact. All subhalos at final snapshot *must* have a tag=2 merger
                   (3, 5) -> flyby, new FOF halo appears where the progenitor was a subhalo and even further in the past the progenitor was a FOF.
                             Characterized by the transitio of FOF -> subhalo -> FOF for the secondary. 3 is the case where the initial
                             relative v1 \cdot v2 was negative, i.e., approaching (FOF) halos. 5 is the case where the the
                             initial "vdotv" is >= 0. 5's are faster interactions typically -> so should be less damaging to the structure

                   4      -> split, a new subhalo appears inside the primary FOF and becomes an independent FOF later
                   1      -> the primary FOF halo fell into a different FOF halo -> which means the primary FOF is now a subhalo
                             so the original case, of tracing the relation between a subhalo and the host FOF halo, can not be
                             followed any more. 
                  -2      -> swap?, the primary FOF does not have a progenitor but the subhalo does. Somehow, this new FOF
                             got enough mass in one time-step to be considered the primary FOF, while the secondary halo,
                             even though it has a progenitor, is considered a subhalo at this timestep.
                             

                   While seemingly simple (a new subhalo appears -> follow into the future with host FOF halo and tag), the
                   tagging is quite complex, since halos can skip snapshots. Ensuring same snapshot for both the future-descendent-of-current-FOF and the
                   future-descendant-of-current-subhalo is incredibly difficult in the general case.

                   An easier way might would be to look at the end-product (so the interaction tag is immediately known) and then go backwards in snapshot.
                   The reason this is simpler is because FirstProg *always* tracks the main branch (or in HINGE parlance, preserves haloid -> because that is
                   how haloid is assigned in HINGE, i.e., haloid is preserved in HINGE along FirstProg(==BigChild in HINGE) *by construction*).

                   Simply following halos back should work. In theory. MS 11th Nov, 2016
                   
                 */

                int32_t firstprog = tree[i].FirstProgenitor;
                int32_t is_fof = tree[i].FirstHaloInFOFgroup == i;
                int32_t desc = tree[i].Descendant;

                //Does this halo have a progenitor?
                if(firstprog == -1) {
                    //So if a new FOF halo appears, then do nothing
                    if(is_fof) {
                        //new FOF halo without progenitor
                        continue;//simply that a new FOF halo has appeared --> *NOT* an interaction
                    } else {
                        //new subhalo without progenitor
                        if(desc == -1) {
                            fprintf(stderr,"New subhalo (halo num = %d) without (both) progenitor, descendant -> probably spurious halo\n",i);
                            continue;//new subhalo without valid progenitor *and* descendant. 
                        } 
                    }
                } 

                //halos that reach here should satisfy: (firstprog == -1 && desc != -1 and is_fof==0) || (firstprog != -1)
                if(is_fof == 1) {
                    //fof halos that reach here *must* have firstprog >=0
                    if( firstprog < 0)  {
                        fprintf(stderr,"Error: (On ThisTask = %d) - Bug in logic. all FOF halos that reach this line should have a valid progenitor. However "
                                "halo i=%d in tree=%d (filename = '%s') has firstprog = %d", ThisTask, i, itree, filename, firstprog);
                        status = EXIT_FAILURE;
                        goto fail;
                    }
                    status = write_fof_halo_interactions_to_stream(tree, nhalos, i, LHALOTREE, fp);
                } else {
                    status = write_subhalo_interactions_to_stream(tree, nhalos, i, LHALOTREE, fp);
                }

                if(status != EXIT_SUCCESS) {
                    goto fail;
                }
                
            }//loop over nhalos
#endif//END of #if 0 -> Code completely commented out since the trees are now converted into HINGE style decomposition into (snapshot, nhalos) + pointers for navigation

