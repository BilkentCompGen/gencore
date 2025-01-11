#ifndef RLOAD_H
#define RLOAD_H

#include "args.h"
#include "utils.h"
#include "tpool.h"
#include "lps.h"
#include <stdint.h>

/**
 * @brief Reads core data from multiple files using multithreading.
 * 
 * This function spawns multiple threads to concurrently read core data from files specified in the 
 * `genome_arguments` structure. The number of threads spawned is controlled by the `thread_number` 
 * parameter in the `program_arguments` structure. Once the reading tasks are completed, the threads are 
 * joined and cleaned up.
 * 
 * @param genome_arguments A reference to a array of `gargs` structures, where each element contains 
 *        file information (e.g., input file names) and is passed to the respective threads for reading.
 * @param program_arguments A reference to the `pargs` structure, which contains general program settings, 
 *        including the number of threads to spawn (`thread_number`).
 */
void read_lcpts(struct gargs *genome_arguments, struct pargs *program_arguments);

/**
 * @brief Reads LCP cores from a file and processes them.
 * 
 * This function loads LCP cores from a specified file and extracts their 
 * hashes for further processing. It manages memory by deleting loaded LCP core objects after 
 * extracting their hashes and setting the results into the provided `genome_arguments`.
 * 
 * @param args A reference to the `gargs` structure that contains the genome-specific 
 *        arguments, including the input cores file name, the output data structures.
 */
void read_lcpt(void *arg);

#endif