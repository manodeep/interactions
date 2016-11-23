#pragma once

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdio.h>
#include <inttypes.h>
    
#ifndef EXTRA_HALO_INFO
#define EXTRA_HALO_INFO
#endif
    
    
struct ctree
{
    float scale;
    int64_t id, num_prog, phantom, pid, upid, mmp;
    int64_t breadth_first_id, depth_first_id, tree_root_id, orig_halo_id, next_coprogenitor_depthfirst_id, last_progenitor_depthfirst_id, tidal_id;
    int32_t snap_num;
    union {
        int32_t desc;
        int32_t Parent;//points to the future
    };

    union {
        int32_t prog;
        int32_t BigChild;//points to the past -> most massive progenitor

    };

    union {
        int32_t next_coprog;
        int32_t Sibling;
    };

    //Indices at the same level
    union {
        int32_t parent;
        int32_t LeastMassiveHost;
    };
    union {
        int32_t uparent;
        int32_t MostMassiveHost;
    };

    float mvir, orig_mvir, rvir, rs, vrms, scale_of_last_MM,
        vmax, pos[3], vel[3], J[3], spin, tidal_force;
    EXTRA_HALO_INFO
};

#ifdef __cplusplus
}
#endif






                             
