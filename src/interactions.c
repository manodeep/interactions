#include <stdio.h>
#include <stdlib.h>

#include "interactions.h"
#include "tree_conversion.h"
#include "io_trees.h"


/* Private functions with actual implementations */
/* For subhalos*/
int write_subhalo_interactions_to_stream_ctrees(struct ctree *tree, const int64_t nhalos, const int64_t index, FILE *fp);
int write_subhalo_interactions_to_stream_lhalotree(struct lhalotree *tree, const int64_t nhalos, const int64_t index, FILE *fp);

/* For FOF halos */
int write_fof_halo_interactions_to_stream_ctrees(struct ctree *tree, const int64_t nhalos, const int64_t index, FILE *fp);
int write_fof_halo_interactions_to_stream_lhalotree(struct lhalotree *tree, const int64_t nhalos, const int64_t index, FILE *fp);


/* Public-facing API */
int write_subhalo_interactions_to_stream(void *tree, const int64_t nhalos, const int64_t index, const tree_kind tree_type, FILE *fp)
{
    switch(tree_type) {
    case(LHALOTREE): return write_subhalo_interactions_to_stream_lhalotree((struct lhalotree *) tree, nhalos, index, fp);
    case(CTREE): return write_subhalo_interactions_to_stream_ctrees((struct ctree *) tree, nhalos, index, fp);
    default: return EXIT_FAILURE;
    }

    return EXIT_FAILURE;//unreachable
}


int write_fof_halo_interactions_to_stream(void *tree, const int64_t nhalos, const int64_t index, const tree_kind tree_type, FILE *fp)
{

    switch(tree_type) {
    case(LHALOTREE): return write_fof_halo_interactions_to_stream_lhalotree((struct lhalotree *) tree, nhalos, index, fp);
    case(CTREE): return write_fof_halo_interactions_to_stream_ctrees((struct ctree *) tree, nhalos, index, fp);
    default: return EXIT_FAILURE;
    }

    return EXIT_FAILURE;//unreachable
}


/* Lhalotree implementation for tagging interactions where the current halo is a FOF halo */
int write_fof_halo_interactions_to_stream_lhalotree(struct lhalotree *tree, const int64_t nhalos, const int64_t index, FILE *fp)
{




    return EXIT_SUCCESS;
}

int write_subhalo_interactions_to_stream_lhalotree(struct lhalotree *tree, const int64_t nhalos, const int64_t index, FILE *fp)
{

    return EXIT_SUCCESS;
}

int write_fof_halo_interactions_to_stream_ctrees(struct ctree *tree, const int64_t nhalos, const int64_t index, FILE *fp)
{

    return EXIT_SUCCESS;
}

int write_subhalo_interactions_to_stream_ctrees(struct ctree *tree, const int64_t nhalos, const int64_t index, FILE *fp)
{

    return EXIT_SUCCESS;
}




struct node_data * walk_tree(struct node_data *start)
{
    struct node_data *tmp = start;

    if( tmp->BigChild != NULL) {
        return tmp->BigChild;
    }   else {
        if(tmp->Sibling !=NULL) {
            return tmp->Sibling;
        } else {

            while( (tmp->Sibling==NULL) && (tmp->Parent !=NULL))
                tmp = tmp->Parent;

            tmp = tmp->Sibling;
        }
    }
    return tmp;

}

