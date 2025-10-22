#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>
#include <stdio.h>

#define DEFAULT_FA_MIN_CC 0
#define DEFAULT_FQ_MIN_CC 32
#define DEFAULT_FA_MAX_CC UINT32_MAX
#define DEFAULT_FQ_MAX_CC UINT32_MAX
#define DEFAULT_SIM_CALC_MODE SET
#define DEFAULT_LCP_LEVEL 5
#define DEFAULT_THREAD_NUMBER 8
#define DEFAULT_VERBOSE 0
#define DEFAULT_WRITE_LCP_CORES 0
#define DEFAULT_PREFIX "gc"
#define DEFAULT_COMPRESSION_RATIO 4
#define MAGIC_LCP_FA_CONSTANT 2.20  // the constant reduction of cores is 2.33 but to be 
                                    // safe, it is selected lower than that
#define MAGIC_LCP_FQ_CONSTANT 2.00  // the constant reduction of cores is 1.5 but to be 
                                    // more efficient, it is selected higher than that
#define INITIAL_SEQUENCE_SIZE 300000000
#define SEPERATOR '$'

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

typedef struct {
    char *fasta;
    char *name;
    int length;
    int seq_idx;
    simple_core *cores;
    uint64_t core_count;
    double exec_time; // exec. time per thread (including LCP processing)
} seq_t;

typedef struct {
    seq_t *seqs;
    int seq_count;
    int capacity;
    long total_seq_len;
    int lcp_level;
    int verbose;
} fa_thread_t;

typedef struct {
    double lcp;
    double merging;
    double sorting;
    double filtering;
} time_stats_t;

typedef struct {
    program_mode mode;
    int n_threads;
    char *prefix;
    int n_genomes;
    sim_calculation_type sct;
    int lcp_level;
    int write_lcpt; // 1: true, 0: false
    int verbose; // 1: true, 0: false
} p_args_t;

typedef struct {
    int apply_filter;
    uint32_t min_cc;
    uint32_t max_cc;
    char *inFileName;
    char *shortName;
    char *outFileName;
    simple_core *cores;
    uint64_t core_count;
    double total_len;
    // other
    sim_calculation_type sct;
    int lcp_level;
    int write_lcpt; // 1: true, 0: false
    int verbose;  // 1: true, 0: false
    time_stats_t time_stats;
} g_args_t;

typedef struct {
    uint64_t value;
    size_t array_index;
    size_t element_index;
} heap_node;

typedef struct {
    heap_node *data;
    size_t size;
    size_t capacity;
} min_heap;


#endif