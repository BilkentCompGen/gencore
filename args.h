#ifndef ARGS_H
#define ARGS_H

#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include <stdatomic.h>

/* =============================================
 * Defaults / constants
 * ============================================= */

#define DEFAULT_FA_MIN_CC               0u
#define DEFAULT_FQ_MIN_CC               32u
#define DEFAULT_FA_MAX_CC               UINT32_MAX
#define DEFAULT_FQ_MAX_CC               UINT32_MAX
#define DEFAULT_SIM_CALC_MODE           SIM_CALC_SET
#define DEFAULT_LCP_LEVEL               5
#define DEFAULT_CORE_OUTSPAN            0
#define DEFAULT_THREAD_NUMBER           8
#define DEFAULT_FQ_READER_NUMBER        1
#define DEFAULT_VERBOSE                 0
#define DEFAULT_WRITE_LCP_CORES         0
#define DEFAULT_PREFIX                  "gc"
#define DEFAULT_COMPRESSION_RATIO       4

/* =============================================
 * Historical estimation constants.
 * ============================================= */

#define MAGIC_LCP_FA_CONSTANT           2.20    // the constant reduction of cores is 2.33 but to be 
                                            // safe, it is selected lower than that
#define MAGIC_LCP_FQ_CONSTANT           2.20    // the constant reduction of cores is 1.5 but to be 
                                            // more efficient, it is selected higher than that
#define INITIAL_SEQUENCE_SIZE           300000000u

/* =============================================
 * FASTQ batching.
 * ============================================= */

#define DEFAULT_FASTQ_BATCH_SIZE        (1u * 1024u * 1024u)
#define DEFAULT_FASTQ_QUEUE_CAPACITY    128
#define SEPARATOR                       '$'
#define FQ_PARALLEL_MIN_CORE_CAP        1024
#define FQ_WORKERS_PER_READER           4

/* =============================================
 * Legacy worker-state constants.
 * ============================================= */

#define THREAD_EXIT_SIGNAL             -2
#define BUFFER_NOT_INITIALIZED         -1
#define THREAD_BUSY                    0
#define THREAD_AVAILABLE               1
#define SHORT_WAIT_TIME                4096
#define MODERATE_WAIT_TIME             32768

/* =============================================
 * Enums
 * ============================================= */

typedef enum {
    INFO,
    WARN,
    ERROR
} log_level_t;

typedef enum {
    PROGRAM_MODE_FA,
    PROGRAM_MODE_FQ,
    PROGRAM_MODE_LOAD
} program_mode_t;

typedef enum {
    SIM_CALC_SET,
    SIM_CALC_VECTOR
} sim_calculation_type_t;

typedef enum {
    FASTQ_INPUT_DIR,
    FASTQ_INPUT_FILE
} fastq_input_type_t;

/* =============================================
 * General
 * ============================================= */

typedef struct {
    double running;
    double idle;
    double lcp;
    double merging;
    double sorting;
    double filtering;
} time_stats_t;

typedef uint64_t simple_core; // first 32 bits are ulabel, last 32 is length of the core

typedef struct {
    simple_core *cores;
    uint64_t count;
    uint64_t capacity;
    uint64_t total_core_len;
    uint64_t total_sequence_len;
} core_result_t;

/* =============================================
 * Configuration structs
 * ============================================= */

 typedef struct {
    char *prefix;
    int n_threads;
    int n_readers;
    int n_genomes;
    int lcp_level;
    int core_span;
    int write_lcpt; // 1: true, 0: false
    int verbose; // 1: true, 0: false
    program_mode_t mode;
    sim_calculation_type_t sct;
} p_args_t;

typedef struct {
    char *inFileName;
    char *shortName;
    char *outFileName;
    int apply_filter;
    int lcp_level;
    int write_lcpt; // 1: true, 0: false
    int verbose;  // 1: true, 0: false
    uint32_t min_cc;
    uint32_t max_cc;
    // other
    sim_calculation_type_t sct;
    core_result_t result;
    time_stats_t time_stats;
} g_args_t;

/* =============================================
 * Helper structs (FA)
 * ============================================= */

typedef struct {
    char *fasta;
    char *name;
    int length;
    int seq_idx;
    double exec_time; // exec. time per thread (including LCP processing)
    core_result_t result;
} seq_t;

typedef struct {
    seq_t *seqs;
    int seq_count;
    int capacity;
    int lcp_level;
    int core_span;
    int verbose;
#if NUMA_AVAILABLE
    int numa_node_id;
#endif
    long total_seq_len;
} fa_thread_t;

/* =============================================
 * Helper structs (FQ)
 * ============================================= */

typedef struct {
    int n_readers;
    int n_workers;
    size_t batch_size;
    int queue_capacity;
    int lcp_level;
    uint64_t estimated_total_core_capacity;
    int verbose;
} fq_parallel_options_t;

typedef struct {
    simple_core *cores;
    uint64_t count;
    uint64_t capacity;
    time_stats_t time_stats;
} fq_parallel_result_t;

typedef struct {
    char *data;
    size_t len;
    int file_id;
} fq_batch_t;

typedef struct {
    fq_batch_t *items;
    int capacity;
    int head;
    int tail;
    int count;
    int closed;

    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} fq_batch_queue_t;

typedef struct {
    fq_batch_queue_t *queue;
    int worker_id;
    int lcp_level;

    simple_core *cores;
    uint64_t count;
    uint64_t capacity;

    uint64_t total_core_len;
    uint64_t total_sequence_len;

    time_stats_t time_stats;
} fq_worker_t;

typedef struct {
    const char **files;
    int n_files;
    atomic_int next_file_id;

    atomic_int active_readers;
    atomic_int active_workers;

    fq_batch_queue_t *queue;
    size_t batch_size;
    int verbose;
} fq_reader_pool_t;

typedef struct {
    int reader_id;
    fq_reader_pool_t *pool;

    fq_worker_t *worker;
} fq_reader_t;

/* =============================================
 * Heap
 * ============================================= */

typedef struct {
    uint64_t value;
    uint64_t element_index;
    size_t array_index;
} heap_node;

typedef struct {
    heap_node *data;
    size_t size;
    size_t capacity;
} min_heap;

#endif