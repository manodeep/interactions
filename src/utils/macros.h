#pragma once
#include <stdio.h>
#include <assert.h>

#define MAXLEN  (1024)

#define ADD_DIFF_TIME(t0,t1)            ((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))
#define REALTIME_ELAPSED_NS(t0, t1)     ((t1.tv_sec - t0.tv_sec)*1000000000.0 + (t1.tv_nsec - t0.tv_nsec))
    
#define ALIGNMENT                32

#define STRINGIFY(x)   #x
#define STR(x) STRINGIFY(x)


#define ABORT(sigterm)                                                  \
    do {                                                                \
        printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
        perror(NULL);                                                   \
    } while(0)

#ifdef NDEBUG
#define XPRINT(EXP, ...)                                do{} while(0)
#else
#define XPRINT(EXP, ...)                                                \
    do { if (!(EXP)) {                                                  \
            printf("Warning in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
        }                                                               \
    } while (0)
#endif

#ifdef NDEBUG
#define XRETURN(EXP, VAL, ...)                                do{} while(0)
#else
#define XRETURN(EXP, VAL, ...)                                           \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
             return VAL;                                                \
         }                                                              \
     } while (0)
#endif


#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                           \
     do { if (!(EXP)) {                                                 \
             fprintf(stderr,"Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
             fprintf(stderr,__VA_ARGS__);                               \
             fprintf(stderr,ANSI_COLOR_BLUE "Hopefully, input validation. Otherwise, bug in code: please email Manodeep Sinha <manodeep@gmail.com>"ANSI_COLOR_RESET"\n"); \
             assert(EXP);                                               \
         }                                                              \
     } while (0)
#endif


/* Macro Constants */
//Just to output some colors
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_RESET   "\x1b[0m"
