#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>//defines int64_t datatype -> *exactly* 8 bytes int

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    DEFAULT = -42, /* present to make the enum a signed variety */
    LHALOTREE = 0, /* Input comes from LHALOTREE output */
    CTREE = 1,     /* Input comes from consistent-tree generated output */
    NTREETYPES     /* Number of different mergertree inputs supported */
} tree_kind;
    
    int write_subhalo_interactions_to_stream(void *tree, const int64_t nhalos, const int64_t index, const tree_kind tree_type, FILE *fp);
    int write_fof_halo_interactions_to_stream(void *tree, const int64_t nhalos, const int64_t index, const tree_kind tree_type, FILE *fp);
    
#ifdef __cplusplus
}
#endif
    
