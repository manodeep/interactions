#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_SNAPSHOTS  (200)

#include "io_trees.h"


struct node_data
{
    union{
        struct lhalotree lht;
        struct ctree ctr;
    };

    /* Structure hierarchy pointers. These are pointers to halos at the same snapshot */
    union{
        struct node_data *FofHalo;/* FOF container at this snapshot */
        struct node_data *FirstHaloInFOFgroup;/* Identical variable type -> only present for Lhalotree users*/
    };
    /* Pointer from HINGE and for Ctrees. Hence, no union for LHaloTree users. Same as FOFHalo in most cases (identical to FOFHalo for LHalotree) */
    struct node_data *ContainerHalo; /* container halo; for a FOF/sub-halo, container is the Fof; a sub-subhalo will have a container as a subhalo and so on*/


    /* Mergertree pointers -> point to halos potentially at different redshifts. */
    union{
        struct node_data *Sibling;/* Pointer to the next halo that shares the same descendent, arranged in decreasing mass order. Might not be at same redshift as currnet halo */
        struct node_data *NextProgenitor;/* Identical variable type -> only present for Lhalotree users*/
    };
    
    union{
        struct node_data *Parent; /* The halo this is going to go into in a future snapshot. Depends on the parent matching algorithm - mostly at the immediately next snapshot */
        struct node_data *Descendant;/* Identical variable type -> only present for Lhalotree users*/
    };
    
    union{
        struct node_data *BigChild; /* Progenitor at a previous snapshot. Similar to Parent (mostly at the immediately previous snapshot but not necessarily so ) */
        struct node_data *FirstProgenitor;/* Identical variable type -> only present for Lhalotree users */
    };

    int64_t haloid;
    double Mtot;
    double Mstar;//MS - added on 6th Oct, 2011. Assigned using SAM models
    double InfallMass;
    
    int32_t halonum;
    int32_t snapshot;
    float redshift;
    int32_t Nchild;
    int32_t npart;

    int32_t isFOF;
    int32_t InfallSnapshot;
    float FormationRedshift;
    float DestructionRedshift;
    float RedshiftofLastMerger;
};
    
    
/* Locally horizontal vertical tree can store vertical store but still process them in snapshot order */
/* This can be forced to behave as a fully horizontal tree simply by storing all halos (from all trees) together in snapshot order */
struct lh_vtree
{
    int64_t Nhalos_per_snap[MAX_SNAPSHOTS];//the number of halos per snapshot (should be directly accessible by snapshot number, i.e., allocated for nsnapshots )
    /* the starting halo pointer per snapshot, for all snapshots.
       To be accessed as tree[isnapshot][ihalo] (where, 0 <= ihalo < Ngroups[isnapshot] and 0 <= isnapshot <= max_snap) */
    struct node_data *tree[MAX_SNAPSHOTS];//declaring as an array of pointers to reduce chances of malloc failure.

    double redshifts[MAX_SNAPSHOTS];//the redshift corresponding to each snapshot -> to be accessed directly by snapshot number
    
    int nsnapshots; // nsnapshots := max_snap + 1 -> here for convenience and to limit possible confusion.
    int min_snap;//the first snapshot -> inclusive
    int max_snap;//the last snapshot -> inclusive. Must be less than MAX_SNAPSHOTS
};

    extern int convert_tree_to_hinge(void *tree, const int64_t nhalos, const tree_kind tree_type, struct lh_vtree *allnodes);
    extern void free_lh_vtree(struct lh_vtree *allnodes);

    /* Tree-walking utility*/
    struct node_data * walk_tree(struct node_data *start);
    extern int assign_haloid(struct lh_vtree *allnodes);
    
#ifdef __cplusplus
}
#endif
