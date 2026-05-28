#ifndef RFASTQ_H
#define RFASTQ_H

#include "args.h"
#include "tqueue.h"
#include "lps.h"
#include "utils.h"
#include <dirent.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <pthread.h>
#include <time.h>


/**
 * @brief Reads multiple FASTQ files concurrently using a thread pool.
 * 
 * This function processes a collection of FASTQ files by spawning multiple
 * threads based on the provided thread arguments (`targs`) and program
 * settings (`pargs`). It manages the number of active threads and ensures
 * that each thread reads a FASTQ file using the `read_fastq` function, which
 * operates on the given thread arguments and shared program settings.
 * 
 * @param genome_args  Reference to a vector of `gargs` structures representing
 *                     the arguments specific to each genome.
 * @param program_args Constant reference to a `pargs` structure representing
 *                     the global program arguments.
 * 
 * @note Parallelization occurs at the FASTQ file level: multiple main threads
 *       read separate FASTQ files and distribute read batches to worker threads
 *       for processing.
 */
void read_fastqs(g_args_t *genome_args, p_args_t *program_args);

/**
 * @brief Reads and processes all FASTQ/FASTQ.GZ files for a genome directory.
 * 
 * This function processes all FASTQ and FASTQ.GZ files located in the directory
 * associated with the given genome (`genome_args`). The genome arguments are
 * already initialized and specify the target genome’s directory and metadata.
 * 
 * Each FASTQ file within the directory is read sequentially, but the reading
 * and processing of reads within each file are parallelized using multiple
 * threads for improved performance.
 * 
 * @param genome_args  Pointer to an initialized `g_args_t` structure representing
 *                     the genome and its associated FASTQ directory.
 * @param program_args Pointer to a `p_args_t` structure containing global program
 *                     parameters and runtime configuration.
 * 
 * @note The function handles both uncompressed (*.fastq) and compressed (*.fastq.gz)
 *       input files. File-level processing is sequential, while read-level operations
 *       within each file are executed concurrently.
 */
void process_dir_fastq(g_args_t *genome_args, p_args_t *program_args);

/**
 * @brief Reads and processes a single FASTQ or FASTQ.GZ file for a genome.
 * 
 * This function reads sequencing reads from a FASTQ or FASTQ.GZ file specified
 * in the genome arguments (`genome_args`). It decodes, parses, and processes
 * the reads in batches using multiple threads to maximize throughput.
 * 
 * @param genome_args  Pointer to a `g_args_t` structure representing the genome
 *                     currently being processed, including the FASTQ file path
 *                     and metadata.
 * @param program_args Pointer to a `p_args_t` structure containing global program
 *                     parameters and runtime configuration.
 * 
 * @note This function is typically invoked by higher-level routines such as
 *       `process_dir_fastq()` or parallel thread workers to handle individual
 *       FASTQ files. It supports both compressed (*.fastq.gz) and uncompressed
 *       (*.fastq) input formats.
 */
void read_fastq(g_args_t *genome_args, p_args_t *program_args);

#endif