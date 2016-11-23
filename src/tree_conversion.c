#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>

#include "utils/macros.h"
#include "interactions.h"
#include "tree_conversion.h"

#include "utils/utils.h"
#include "utils/sglib.h"

int convert_tree_to_hinge_lhalotree(struct lhalotree *tree, const int64_t nhalos, struct lh_vtree *allnodes);
int convert_tree_to_hinge_ctree(struct ctree *tree, const int64_t nhalos, struct lh_vtree *allnodes);
int arrange_lhalotree_into_snapshot_halo(struct lhalotree *tree, const int64_t nhalos, struct lh_vtree *allnodes, uint64_t *new_index);
uint64_t set_snapshot_and_halo_index(const uint32_t snap, const uint32_t halonum);
void get_snapshot_and_halo_index(const int64_t new_index, int32_t *snap, int32_t *halonum);

void free_lh_vtree(struct lh_vtree *allnodes)
{
    for(int i=allnodes->min_snap;i<=allnodes->max_snap;i++) {
        free(allnodes->tree[i]);
        allnodes->tree[i] = NULL;
        allnodes->Nhalos_per_snap[i] = 0;
    }
}

int convert_tree_to_hinge(void *tree, const int64_t nhalos, const tree_kind tree_type, struct lh_vtree *allnodes)
{
    switch(tree_type) {
    case(LHALOTREE): return convert_tree_to_hinge_lhalotree((struct lhalotree *) tree, nhalos, allnodes);
    /* case(CTREE): return convert_tree_to_hinge_ctree((struct ctree *) tree, nhalos, allnodes); */
    default: return EXIT_FAILURE;
    }

    return EXIT_FAILURE;//unreachable
}

double behroozi_fx(const double x, const double alpha, const double delta, const double gamma)
{
  if(x == 0.0) {
    return -log10(2.0) + delta* pow(log10(2.0), gamma)/(1.0 + exp(1.0));
  } else {
    double first_term = -log10(pow(10.0,alpha*x) + 1.0);
    double second_term_numerator = pow(log10(1.0 + exp(x)), gamma);
    double second_term_denom   = 1.0 + exp(pow(10.0,-x));
    double second_term = delta * second_term_numerator/second_term_denom;
    return first_term + second_term;
  }
}  


double assign_stellar_mass_from_mvir(struct node_data * const thisnode,int model)
{
  //Taken from data compiled by Stewart, K arxiv:1109.3207v1 Table 1
  double alpha,beta,gamma,m,M1,M2;//
  double mstar=0.0;
  float z = thisnode->redshift;
  const double ActualMassUnits = 1e10;
  double Mvir = thisnode->InfallMass*ActualMassUnits;

  const char modelnames[][MAXLEN] = {"Conroy & Wechsler (2009)","Moster et al (2010)","Behroozi et al (2010)","Behroozi et al (2013)"};
  const int nmodels = sizeof(modelnames)/MAXLEN;
  const float MinZ_for_Models[] = {2.0,2.0,2.0,8.0};

  if(model < nmodels) {
    if(z <= MinZ_for_Models[model]) {
        switch(model) {
        case 0: //Conroy & Wechsler (2009)
            {
                M1 = exp10(0.056*z*z + 0.068*z + 9.5);
                M2 = exp10(0.320*z*z + 0.018*z + 11.2);
                alpha = 0.021*pow(z,4.86) + 3.39;
                beta  = 0.085*z + 0.36;
                
                mstar = M1*pow(Mvir,alpha)*pow(M2,-beta)*pow(0.5*(M2+Mvir),  beta-alpha);
                break;
            }
        case 1://Moster (2010)
            {
                m  = 0.0282*pow(1.0+z, -0.72);
                M1 = exp10( 11.884*pow(1.0+z, 0.019));
                beta  = 0.17*z + 1.06;
                gamma = 0.556*pow(1.0+z, -0.26);
                
                mstar = 2.0*Mvir*m/( pow(Mvir/M1, -beta) + pow(Mvir/M1,gamma) );
                break;
            }
        case 2://Behroozi (2010)
            {
                M1 = exp10( 0.03500*z*z - 0.19200*z + 10.199);
                M2 = exp10( 0.00509*z*z + 0.00299*z + 11.824);
                alpha =       -0.20760*z*z + 0.75200*z + 2.423;
                beta  =        0.12000*z*z - 0.09940*z + 0.206;  
                
                mstar = M1*pow(Mvir,alpha)*pow(M2,-beta)*pow(0.5*(M2+Mvir),  beta-alpha);
                break;
            }
        case 3://Behroozi (2013)
            {
                double scale_factor = 1.0/(1.0 + (double) z);
                double nu = exp(-4.0*scale_factor*scale_factor);
                double log10_epsilon = -1.777 + (-0.006*(scale_factor-1.0) + (-0.0)*z ) * nu +
                    (-0.119 * (scale_factor-1.0));
                
                double log10M1 = 11.514 + (-1.793*(scale_factor-1.0) + (-0.251)*z ) * nu ;
                alpha = -1.412 + (0.731*(scale_factor-1.0))*nu;
                double delta = 3.508 + (2.608 *(scale_factor-1.0) + (-0.043)*z )*nu;
                gamma = 0.316 + (1.319 *(scale_factor-1.0) + (0.279 )*z )*nu; 
                
                double first_term  = log10_epsilon + log10M1;
                double log10Mh_over_M1 = log10(Mvir)-log10M1;
                double second_term = behroozi_fx(log10Mh_over_M1,alpha,delta,gamma);
                double third_term  = behroozi_fx(0.0, alpha,delta,gamma);
                
                mstar = exp10(first_term + second_term - third_term);
                break;
            }
        default:
            {
                fprintf(stderr,"Mvir-Mstar model = %d not implemented\n The options for assigning stellar mass as a function of Mvir are :\n",model);
                for(int i=0;i<nmodels;i++)
                    fprintf(stderr,"%s  [%d]\n",modelnames[i],i);
                exit(EXIT_FAILURE);
            }
        }
        
      if(mstar <= 0.0 || mstar >= Mvir) {
	fprintf(stderr,"mstar has an unphysical value (with model = %s [option %d]). Mvir = %lf at z = %f with mstar = %lf\n",modelnames[model],model,Mvir,z,mstar);
	fprintf(stderr,"exiting..\n");
	exit(EXIT_FAILURE);
      }
    } else {
      mstar = 0.0;
    }
  } else {
    fprintf(stderr,"Mvir-Mstar model = %d not implemented\n The options for assigning stellar mass as a function of Mvir are :\n",model);
    for(int i=0;i<nmodels;i++)
      fprintf(stderr,"%s  [%d]\n",modelnames[i],i);
    exit(EXIT_FAILURE);
  }
  return mstar/ActualMassUnits;
}



int assign_haloid(struct lh_vtree *allnodes)
{
  int64_t haloid = 0;
  const int model = 3;//0 - Conroy & Wechsler (2009), 1 - Moster (2010), 2 - Behroozi (2010), 3 - Behroozi (2013)
  
  /* reset all haloid's (in case the tree has been modified) */
  for(int isnapshot=allnodes->min_snap;isnapshot<=allnodes->max_snap;isnapshot++) {
      struct node_data *base = allnodes->tree[isnapshot];
      const int64_t nhalos = allnodes->Nhalos_per_snap[isnapshot];
      for(int64_t igroup=0;igroup<nhalos;igroup++) {
          struct node_data *thisnode = &base[igroup];
          thisnode->haloid = -1;
          
          /* Make sure that thisnode->Nchild is correct */
          if(thisnode->BigChild != NULL) {
              thisnode->Nchild = 1;
              struct node_data *tmp_node = thisnode->BigChild->Sibling;
              while(tmp_node != NULL) {
                  thisnode->Nchild++;
                  tmp_node = tmp_node->Sibling;
              }
          }
      }
  }

  for(int isnapshot=allnodes->max_snap;isnapshot>=allnodes->min_snap;isnapshot--) {
      struct node_data *base = allnodes->tree[isnapshot];
      const int64_t nhalos = allnodes->Nhalos_per_snap[isnapshot];
      for(int64_t igroup=0;igroup<nhalos;igroup++) {
          struct node_data *thisnode = &base[igroup];
          if (thisnode->haloid < 0) {
              float destructionz = thisnode->redshift;
              float formationz = -1.0;
              thisnode->haloid = haloid;
              
              while(thisnode->BigChild != NULL) {
                  thisnode = thisnode->BigChild;
                  if (thisnode->haloid >= 0) {
                      fprintf(stderr,"\n This should not have happened.. found a haloid while assigning haloids -- exiting \n");
                      fprintf(stderr,"snapshot = %d this haloid = %"PRId64"  this halonum = %d this parent haloid  = %"PRId64"  this parent halonum = %d snapshot = %d\n",
                              thisnode->snapshot,thisnode->haloid,thisnode->halonum,thisnode->Parent->haloid,thisnode->Parent->halonum,thisnode->Parent->snapshot);
                      fprintf(stderr,"thisnode->parent->Nchild = %d thisnode->parent->bigchild->haloid = %"PRId64" at snapshot = %d \n",
                              thisnode->Parent->Nchild,thisnode->Parent->BigChild->haloid,thisnode->Parent->BigChild->snapshot);
                      
                      return EXIT_FAILURE;
                  } else {
                      thisnode->haloid = haloid;
                      thisnode->DestructionRedshift = destructionz;
                      formationz = thisnode->redshift;
                  }
              }

              double infallmass = thisnode->Mtot;
              int32_t infallsnap    = thisnode->snapshot;
              while(thisnode != NULL) {
                  if(thisnode->isFOF == 1) {
                      infallmass = thisnode->Mtot;
                      infallsnap    = thisnode->snapshot;
                  }
                  
                  if(thisnode->haloid == haloid) {
                      thisnode->FormationRedshift = formationz;
                      thisnode->InfallMass = infallmass;
                      thisnode->InfallSnapshot = infallsnap;
                      thisnode->Mstar = assign_stellar_mass_from_mvir(thisnode,model);
                      
                      //ensure that the stellar mass does not reduce. 
                      if(thisnode->BigChild != NULL && thisnode->BigChild->Mstar > thisnode->Mstar)
                          thisnode->Mstar = thisnode->BigChild->Mstar;
                      
                  } else {
                      break;
                  }
                  thisnode = thisnode->Parent;
              }
              haloid++;
          }
      }
  }
  
  return EXIT_SUCCESS;
}

uint64_t set_snapshot_and_halo_index(const uint32_t snap, const uint32_t halonum)
{
    return (((uint64_t) snap) << 32) |  halonum;
}


void get_snapshot_and_halo_index(const int64_t new_index, int32_t *snap, int32_t *halonum)
{
    *snap = new_index >> 32;//select the top 32 bits
    *halonum = (new_index << 32) >> 32;//select the bottom 32 bits. 
}



int convert_tree_to_hinge_lhalotree(struct lhalotree *tree, const int64_t nhalos, struct lh_vtree *allnodes)
{
    /* Validate tree */
    for(int64_t i=0;i<nhalos;i++) {
        if(tree[i].NextProgenitor == -1) continue;
        XRETURN(tree[i].NextProgenitor >=0 && tree[i].NextProgenitor < nhalos, EXIT_FAILURE,
                "halo = %"PRId64", NextProgenitor = %d must be within [0, %"PRId64")\n", i, tree[i].NextProgenitor, nhalos);
    }


    /* This sorts and then also fixes the mergertree indices within this new sorted order
       The last argument is to setup tests within the repo and should be always 0 (as in no need to test the sorting routines here)
     */
    int status = sort_lhalotree_in_snapshot_and_fof_groups(tree, nhalos, 0);
    if(status != EXIT_SUCCESS) {
        goto fail;
    }

    /* Validate tree */
    for(int64_t i=0;i<nhalos;i++) {
        if(tree[i].NextProgenitor == -1) continue;
        XRETURN(tree[i].NextProgenitor >=0 && tree[i].NextProgenitor < nhalos, EXIT_FAILURE,
                "halo = %"PRId64", NextProgenitor = %d must be within [0, %"PRId64")\n", i, tree[i].NextProgenitor, nhalos);
    }


    
    /* Reset the number of halos per snapshot */
    for(int32_t i=0;i<MAX_SNAPSHOTS;i++) {
        allnodes->Nhalos_per_snap[i] = 0;
        free(allnodes->tree[i]);
        allnodes->tree[i] = NULL;
    }

    int32_t min_snap=MAX_SNAPSHOTS, max_snap=-1;
    for(int64_t i=0;i<nhalos;i++) {
        struct lhalotree *this_halo = &(tree[i]);
        if(this_halo->SnapNum < min_snap) {
            min_snap = this_halo->SnapNum;
        }
        
        if(this_halo->SnapNum > max_snap) {
            max_snap = this_halo->SnapNum;
        }
    }
    if(min_snap < 0 || max_snap < 0) {
        fprintf(stderr,"Error: Min/max. snapshot number=(%d,%d) must be greater than 0\n", min_snap, max_snap);
        return EXIT_FAILURE;
    }
    if(max_snap >= MAX_SNAPSHOTS) {
        fprintf(stderr,"Error: Max. snapshot number=%d must be <= %d. Please change `MAX_SNAPSHOTS` in the `tree_utils.h` "
                "(warning generated in `%s`)\n", max_snap, MAX_SNAPSHOTS, __FILE__);
        return EXIT_FAILURE;
    }

    allnodes->min_snap = min_snap;
    allnodes->max_snap = max_snap;
    allnodes->nsnapshots = max_snap + 1;
    
    /* Create an array with the index mappings */
    /* Even though nhalos is 64 bit integer, LHaloTree can only handle int32_t nhalos (by design).
       However, I need two indices now -> one for snapshot and another for index with allnodes->tree[isnap]
       -> Using a 64 bit integer: top 32 bits are for snapshot while bottom are for the index in allnodes->tree[isnap]
    */
    uint64_t *new_index = my_malloc(sizeof(*new_index), nhalos);
            
    /* Now figure out the mapping between the Lhalo Vertical Tree within the new 2_D (snapshot, halos) arrangement */
    /* Also, figure out how many halos per snapshot */
    for(int64_t i=0;i<nhalos;i++) {
        struct lhalotree *this_halo = &(tree[i]);
        int32_t snap = this_halo->SnapNum;
        if(snap < allnodes->min_snap || snap > allnodes->max_snap) {
            fprintf(stderr,"Error: snap = %d must be within [%d, %d]", snap, allnodes->min_snap, allnodes->max_snap);
            status = EXIT_FAILURE;
            goto fail;
        }

        //new_index[i] = (((int64_t) snap) << 32) +  allnodes->Nhalos_per_snap[snap];
        //Upper 32 bits are set by snap while lower 32 bits are nhalos (for a total of 64 bits in new_index)
        new_index[i] = set_snapshot_and_halo_index((uint32_t) snap, (uint32_t) allnodes->Nhalos_per_snap[snap]);

        allnodes->Nhalos_per_snap[snap]++;
        if(allnodes->Nhalos_per_snap[snap] > INT_MAX) {
            fprintf(stderr,"Error: Indexing scheme will not work since at snap = %d Nhalos = %"PRId64" is larger than INT_MAX\n",
                    snap, allnodes->Nhalos_per_snap[snap]);
            status = EXIT_FAILURE;
            goto fail;
        }
        

        {
            /* Test that the algorithm is correct */
            int test_snap, test_idx;
            get_snapshot_and_halo_index(new_index[i], &test_snap, &test_idx);
            assert(test_snap == snap && test_idx != allnodes->Nhalos_per_snap[snap]);
        }
    }


    /* Allocate memory for each snapshot */    
    for(int32_t i=allnodes->min_snap;i<=allnodes->max_snap;i++) {
        const int64_t nhalos_snap = allnodes->Nhalos_per_snap[i];
        if(nhalos_snap > 0) {
            allnodes->tree[i] = my_malloc(sizeof(struct node_data), nhalos_snap);
            if(allnodes->tree[i] == NULL) {
                fprintf(stderr,"Error: malloc failed for allocating tree at snapshot = %d with nhalos = %"PRId64"\n",i, nhalos_snap);
                status = EXIT_FAILURE;
                goto fail;
            }
        } 
    }
    
    /* Memory has been allocated everything is set to copy the lhalotree into the new (snapshot, nhalos_per_snap) struct */
    status = arrange_lhalotree_into_snapshot_halo(tree, nhalos, allnodes, new_index);
    if(status != EXIT_SUCCESS) {
        /* fprintf(stderr,"Error: Failed to fix mergertree indices into the new (snapshot, nhalos) order\n"); */
        goto fail;
    }
    free(new_index);
    return status;
    
 fail:
    fprintf(stderr,"Failed\n");
    free_lh_vtree(allnodes);

    return status;
    
}



int arrange_lhalotree_into_snapshot_halo(struct lhalotree *tree, const int64_t nhalos, struct lh_vtree *allnodes, uint64_t *new_index)
{
#define CHECK_FIELD_VALUE(FIELD) {              \
        const int32_t index = this_halo->FIELD; \
        if(index == -1 || (index >= 0 && index < nhalos)) { \
            /* all good */                                              \
        } else {                                                        \
            fprintf(stderr, "For field = `"#FIELD"` index = %d must be either -1 or within [0, %"PRId64") \n", index, nhalos); \
            assert(0);                                                  \
        }                                                               \
    }
    
#define CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(FIELD) {        \
        const int32_t index = this_halo->FIELD;                         \
        if(index >= 0) {                                                \
            assert(index < nhalos);                                     \
            int32_t dst_snap, dst_halo_idx;                             \
            get_snapshot_and_halo_index(new_index[index], &dst_snap, &dst_halo_idx); \
            if(dst_snap >= allnodes->min_snap && dst_snap <= allnodes->max_snap) { \
                struct node_data *dst = &((allnodes->tree[dst_snap])[dst_halo_idx]); \
                src->FIELD = dst;                                       \
            } else {                                                    \
                fprintf(stderr,"snap = %d is not within [min_snap, max_snap] = [%d, %d] \n", dst_snap,  allnodes->min_snap, allnodes->max_snap); \
                assert(0);                                              \
            }                                                           \
        } else {                                                        \
            src->FIELD = NULL;                                          \
        }                                                               \
    }
    
    for(int64_t i=0;i<nhalos;i++) {
        struct lhalotree *this_halo=&(tree[i]);
        int32_t snap, ihalo;
        get_snapshot_and_halo_index(new_index[i], &snap, &ihalo);
        XRETURN(snap <= allnodes->max_snap, EXIT_FAILURE, "snap = %d must be less than max_snap = %d\n", snap, allnodes->max_snap);
        XRETURN(ihalo <= allnodes->Nhalos_per_snap[snap], EXIT_FAILURE,
                "ihalo = %d must be less than Nhalos_per_snap[%d] = %"PRId64"\n",ihalo, snap,  allnodes->Nhalos_per_snap[snap]);
        struct node_data *src = &((allnodes->tree[snap])[ihalo]);
        src->lht = *this_halo;//Copy the entire Lhalotree struct (NOT pointer. all of struct lhalotree gets copied)
        src->snapshot = this_halo->SnapNum;
        src->halonum = ihalo;

        CHECK_FIELD_VALUE(FirstProgenitor);
        CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(FirstProgenitor);

        CHECK_FIELD_VALUE(NextProgenitor);
        CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(NextProgenitor);        

        CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(Descendant);

        CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(FirstHaloInFOFgroup);

        /* CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX(NextHaloInFOFgroup); */


        /* Populate and initialize the standard fields */
        src->Mtot = this_halo->Mvir;//Mvir
        src->npart = this_halo->Len;
        src->InfallMass = -1.0;
        src->Mstar = -1.0;//MS - added on 6th Oct, 2011. Assigned using SAM models
        src->InfallSnapshot = -1;
        src->FormationRedshift = -1.0;
        src->DestructionRedshift = -1.0;
        src->RedshiftofLastMerger = -1.0;
        src->isFOF = (this_halo->FirstHaloInFOFgroup == i) ? 1:0;
        src->redshift = allnodes->redshifts[this_halo->SnapNum];
    }
    
#undef CONVERT_LINEAR_INDEX_TO_SNAPSHOT_AND_HALO_INDEX    

    return EXIT_SUCCESS;
}

