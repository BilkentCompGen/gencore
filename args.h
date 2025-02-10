#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>

#define MAGIC_LCP_FA_CONSTANT 2.20  // the constant reduction of cores is 2.33 but to be 
                                    // safe, it is selected lower than that

#define MAGIC_LCP_FQ_CONSTANT 2.00  // the constant reduction of cores is 1.5 but to be 
                                    // more efficient, it is selected higher than that

typedef enum {
    INFO,
    WARN,
    ERROR
} LogLevel;

typedef enum {
    FA,
    FQ,
    LOAD
} program_mode;

typedef enum {
    SET,
    VECTOR
} sim_calculation_type;

typedef uint64_t simple_core; // first 32 bits are ulabel, last 32 is length of the core

struct pargs {
    program_mode mode;
    int thread_number;
    char *prefix;
    int number_of_genomes;
};

struct gargs {
    int apply_filter;
    uint32_t min_cc;
    uint32_t max_cc;
    char *inFileName;
    char *shortName;
    char *outFileName;
    uint64_t cores_len;
    simple_core *cores;
    double total_len;
    // other
    sim_calculation_type sct;
    int lcp_level;
    int write_lcpt; // 1: true, 0: false
    int verbose;  // 1: true, 0: false
};

#endif