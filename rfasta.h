#ifndef RFASTA_H
#define RFASTA_H

#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include "args.h"
#include "logging.h"
#include "helper.h"
#include "lps.h"
#include "fileio.h"


/**
 * @brief Reads multiple FASTA files concurrently using a pool of threads.
 * 
 * This function processes a collection of FASTA files by spawning threads 
 * based on the thread arguments (`targs`) and program settings (`pargs`). It 
 * manages the number of concurrent threads and ensures that each thread reads 
 * a FASTA file using the `read_fasta` function, which operates on the provided 
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
 * - The function continues to launch new threads until all FASTA files in 
 *   `thread_arguments` have been processed.
 * - Once all threads are launched, the function ensures that all threads are 
 *   joined before exiting, preventing any orphan threads from continuing execution.
 */
void read_fastas( std::vector<struct targs>& thread_arguments, const struct pargs& program_arguments );

/**
 * @brief Reads a FASTA file and processes its sequences using the LCP (Locally Consistent Parsing) method.
 * 
 * This function is responsible for reading a FASTA file specified in the `thread_arguments`, 
 * processing each sequence using LCP techniques, and storing the results. The function manages 
 * logging for verbose output, tracks the size of processed sequences, and handles thread-safe 
 * operations, as it is designed to be run in a multithreaded environment.
 * 
 * @param thread_arguments A reference to the `targs` structure that contains the thread-specific 
 *        arguments, including the input FASTA file name, the output data structures (`cores`, 
 *        `counts`), and size tracking.
 * @param program_arguments A constant reference to the `pargs` structure that contains the 
 *        program-wide settings, such as the LCP depth level (`lcpLevel`), verbosity, and whether 
 *        to write LCP cores to file.
 * 
 * @details
 * - The function opens the FASTA file specified in `thread_arguments.inFileName` and processes 
 *   each chromosome or sequence individually.
 * - For each sequence, an `lps` (locally parsed string) object is created, and its depth is increased using 
 *   `deepen()`. The sequence is then stored in a vector for further processing.
 * - If the `verbose` flag in `program_arguments` is set, the function logs detailed information about each 
 *   sequence, including its ID and size.
 * - Once all sequences are processed, the function optionally saves the LCP cores to a file if the `writeLcpt` 
 *   flag is enabled.
 * - Finally, the LCP cores are flattened and processed to generate a signature and initialize sets and counts.
 * - The function ensures proper memory management by cleaning up dynamically allocated `lps` objects.
 * 
 * @note This function is designed to be run in a multithreaded environment, with each thread handling a 
 *       separate FASTA file or sequence set. The function logs the initiation and completion of each thread 
 *       using the thread ID for traceability.
 * 
 * @see flatten(), generateSignature(), initializeSetAndCounts(), save(), log()
 */
void read_fasta( struct targs& thread_arguments, const struct pargs& program_arguments );

#endif