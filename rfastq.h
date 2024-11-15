#ifndef RFASTQ_H
#define RFASTQ_H

#include <cstdint>
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "args.h"
#include "lps.h"
#include "similarity_metrics.h"
#include "helper.h"
#include "utils/GzFile.hpp"
#include "fileio.h"


/**
 * @brief Reads multiple FASTQ files concurrently using a pool of threads.
 * 
 * This function processes a collection of FASTQ files by spawning threads 
 * based on the thread arguments (`targs`) and program settings (`pargs`). It 
 * manages the number of concurrent threads and ensures that each thread reads 
 * a FASTQ file using the `read_fastq` function, which operates on the provided 
 * thread arguments and shared program settings. 
 * 
 * @param thread_arguments A reference to a vector of `targs` structures 
 *        representing the arguments specific to each thread.
 * @param program_arguments A constant reference to a `pargs` structure 
 *        representing the global program arguments.
 * 
 * @details 
 * - The function uses a pool of threads, limiting the number of active threads 
 *   based on `program_arguments.threadNumber`. 
 * - Threads are launched until the maximum thread limit is reached, and once a 
 *   thread finishes its work, it is joined and removed from the active thread pool.
 * - The function continues to launch new threads until all FASTQ files in 
 *   `thread_arguments` have been processed.
 * - Once all threads are launched, the function ensures that all threads are 
 *   joined before exiting, preventing any orphan threads from continuing execution.
 */
void read_fastqs( std::vector<struct targs>& thread_arguments, const struct pargs& program_arguments );


/**
 * @brief Processes a genome file to extract LCP cores using multiple threads.
 *
 * This function reads genomic sequences from a specified file and distributes the processing
 * tasks among several worker threads. Each thread computes LCP cores for the sequences at a
 * given LCP level and aggregates these cores into a shared vector. The function tracks the
 * total number of reads processed and their combined length. It ensures efficient and
 * thread-safe handling of genomic data, leveraging parallel processing to enhance performance.
 *
 * @param filename Path to the input file containing genomic sequences.
 * @param infile Reference to an open GzFile object for reading the input file.
 * @param lcp_level The depth of LCP analysis for extracting cores from sequences.
 * @param cores A reference to a shared vector where extracted LCP cores are aggregated.
 * @param read_count Reference to a variable that will hold the total number of processed reads.
 * @param total_length Reference to a variable that will hold the combined length of all reads.
 * @param thread_number The number of worker threads to use for processing.
 */
void read_fastq( struct targs& arguments, const struct pargs program_arguments );


#endif