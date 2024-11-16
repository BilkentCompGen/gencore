#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include <string>
#include <thread>
#include <fstream>
#include <iostream>
#include <sstream>
#include "args.h"
#include "lps.h"
#include "helper.h"
#include "logging.h"



/**
 * @brief Saves the current state of the lcp::lps object to the specified output file.
 * 
 * This function writes a `false` value (indicating work is not complete),
 * followed by the serialized data of the lps object, to a binary file.
 *
 * @param out Lps object will be stored using this oftream.
 * @param str Pointer to the lcp::lps object that is being saved.
 */
void save( std::ofstream& out, const lcp::lps* str );


/**
 * @brief Marks the thread's processing as done by writing a boolean flag and size to the specified output file.
 * 
 * This function writes a `true` value to indicate the completion of the thread's work,
 * followed by the size of the data processed, to a binary file.
 *
 * @param out Notifying that writing is completed will be done using this oftream.
 */
void done( std::ofstream& out );


/**
 * @brief Reads core data from multiple files using multithreading.
 * 
 * This function spawns multiple threads to concurrently read core data from files specified in the 
 * `thread_arguments` structure. The number of threads spawned is controlled by the `threadNumber` 
 * parameter in the `program_arguments` structure. Once the reading tasks are completed, the threads are 
 * joined and cleaned up.
 * 
 * @param thread_arguments A reference to a vector of `targs` structures, where each element contains 
 *        file information (e.g., input file names) and is passed to the respective threads for reading.
 * @param program_arguments A reference to the `pargs` structure, which contains general program settings, 
 *        including the number of threads to spawn (`threadNumber`).
 */
void read_lcpts( std::vector<struct targs>& thread_arguments, struct pargs& program_arguments );


/**
 * @brief Reads LCP cores from a file and processes them.
 * 
 * This function loads LCP (Longest Common Prefix) cores from a specified file and extracts their 
 * hashes for further processing. It manages memory by deleting loaded LCP core objects after 
 * extracting their hashes and setting the results into the provided `thread_arguments`.
 * 
 * @param thread_arguments A reference to a `targs` structure containing file-specific data (e.g., input file name) 
 *        and will be updated with the extracted LCP cores and their counts.
 * @param program_arguments A reference to a `pargs` structure containing program-wide settings needed 
 *        for loading and processing the LCP cores.
 */
void read_lcpt( struct targs& thread_arguments, struct pargs& program_arguments );


#endif
