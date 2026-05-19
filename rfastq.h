#ifndef RFASTQ_H
#define RFASTQ_H

#include "args.h"
#include "utils.h"
#include "tpool.h"
#include "lps.h"
#include <htslib/kseq.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <time.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>

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
 * @param fq_type      A mode to determine if the genome consists of a single 
 *                     FASTQ file or not. It helps to determine if the function
 *                     is called from `read_fastqs` or `process_dir_fastq`.
 * 
 * @note The function handles both uncompressed (*.fastq) and compressed (*.fastq.gz)
 *       input files. File-level processing is sequential, while read-level operations
 *       within each file are executed concurrently.
 */
void process_dir_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type_t fq_type);

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
 * @param fq_type      A mode to determine if the genome consists of a single 
 *                     FASTQ file or not. It helps to determine if the function
 *                     is called from `read_fastqs` or `process_dir_fastq`.
 * 
 * @note This function is typically invoked by higher-level routines such as
 *       `process_dir_fastq()` or parallel thread workers to handle individual
 *       FASTQ files. It supports both compressed (*.fastq.gz) and uncompressed
 *       (*.fastq) input formats.
 */
void read_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type_t fq_type);

/**
 * @brief Thread entry function that processes a batch of FASTQ reads.
 * 
 * This function serves as a worker thread routine, executed by each thread
 * created to process assigned batches of FASTQ reads. The main thread divides
 * the input data into batches and passes a thread argument structure to each
 * worker for processing.
 * 
 * @param args Pointer to a thread argument structure containing the read batch,
 *             genome context, and shared program parameters.
 * 
 * @note Designed to be used as a thread function (e.g., passed to pthread_create()).
 *       Each worker thread processes its batch independently to enable
 *       parallelized read processing.
 */
void process_reads(void *args);

/**
 * @brief Reads multiple FASTQ files concurrently using a pool of threads.
 * 
 * This function processes a collection of FASTQ files by spawning threads 
 * based on the thread arguments (`targs`) and program settings (`pargs`). It 
 * manages the number of concurrent threads and ensures that each thread reads 
 * a FASTQ file using the `read_fastq` function, which operates on the provided 
 * thread arguments and shared program settings. 
 * 
 * @param genome_args A reference to a vector of `gargs` structures 
 *        representing the arguments specific to each genome.
 * @param program_args A constant reference to a `pargs` structure 
 *        representing the global program arguments.
 */
void read_fastqs_trivial(g_args_t *genome_args, p_args_t *program_args);

/**
 * @brief Processes a genome files to extract LCP cores using multiple threads.
 *
 * This function reads genomic sequences from a files found in specified directory and
 * computes LCP cores for the sequences at a given LCP level and aggregates these cores 
 * into a shared array. Then, it combined the cores of each file into single array. 
 * The function tracks the total number of reads processed and their combined length. 
 * It ensures efficient and thread-safe handling of genomic data, leveraging parallel 
 * processing to enhance performance.
 *
 * @param args A reference to the `gargs` structure that contains the genome-specific 
 *        arguments, including the input FASTQ file name, the output data structures.
 */
void process_dir_fastq_trivial(void *arg);

/**
 * @brief Processes a genome file to extract LCP cores using multiple threads.
 *
 * This function reads genomic sequences from a specified file and computes LCP cores 
 * for the sequences at a given LCP level and aggregates these cores into a shared array. 
 * The function tracks the total number of reads processed and their combined length. 
 * It ensures efficient and thread-safe handling of genomic data, leveraging parallel 
 * processing to enhance performance.
 *
 * @param args A reference to the `gargs` structure that contains the genome-specific 
 *        arguments, including the input FASTQ file name, the output data structures.
 */
void read_fastq_trivial(void *arg);

/**
 * @brief Processes a DNA sequence for both forward and reverse complement strands.
 *
 * This function processes a given DNA sequence in both its forward and reverse 
 * complement forms. It initializes the `lps` structure for both strands, deepens 
 * the `lps` structure based on the specified LCP level, and saves the processed 
 * result if the `write_lcpt` flag is set. After processing both strands, the 
 * allocated memory for the `lps` structures is freed.
 *
 * @param sequence A pointer to the DNA sequence to be processed.
 * @param seq_size The length of the DNA sequence.
 * @param capacity The pointer to the capacity value of the cores array.
 * @param genome_args Pointer to the genome arguments structure, which 
 *        contains settings such as the LCP level and whether to save results.
 */
void process_reads_trivial(char *sequence, uint64_t seq_size, uint64_t *capacity, g_args_t *genome_args);

#endif