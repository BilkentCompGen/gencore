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
 * @brief Marks the thread's processing as done by writing a boolean flag and size to the specified output file.
 * 
 * This function writes a `true` value to indicate the completion of the thread's work,
 * followed by the size of the data processed, to a binary file.
 *
 * @param thread_arguments A struct containing the output file name and size of data processed.
 */
void done( const struct targs& thread_arguments );


/**
 * @brief Saves the current state of the lcp::lps object to the specified output file.
 * 
 * This function writes a `false` value (indicating work is not complete),
 * followed by the serialized data of the lps object, to a binary file.
 *
 * @param thread_arguments A struct containing the output file name and other thread-related parameters.
 * @param str Pointer to the lcp::lps object that is being saved.
 */
void save( const struct targs& thread_arguments, const lcp::lps* str );

/**
 * @brief Saves the LCP cores to a binary file.
 * 
 * This function writes the contents of the provided LCP cores to a binary file specified 
 * by the `outFileName` in the `thread_arguments` structure. It also records metadata about the 
 * number of cores and the total size of the genome being processed.
 * 
 * @param thread_arguments A constant reference to the `targs` structure, which contains file output 
 *        information (such as the output file name) and the total size of the processed genome.
 * @param cores A constant reference to a vector of pointers to `lcp::lps` objects, representing
 *        the LCP cores that will be saved to the file.
 */
void save( const struct targs& thread_arguments, const std::vector<lcp::lps*>& cores );

/**
 * @brief Loads LCP cores from a binary file into a vector.
 * 
 * This function reads the specified binary file to load LCP core data into a provided vector of 
 * `lcp::lps` pointers. It extracts the total number of cores and the associated genome size from 
 * the file, which are then stored in the respective fields of the `thread_arguments`.
 * 
 * @param thread_arguments A reference to a `targs` structure that contains the input file name and
 *        will be updated with the total genome size after loading the cores.
 * @param program_arguments A reference to a `pargs` structure that contains program-wide settings, 
 *        including the level of LCP depth to be applied to each loaded core.
 * @param cores A reference to a vector of pointers to `lcp::lps` objects, where the loaded LCP cores 
 *        will be stored.
 */
void load( struct targs& arguments, struct pargs& program_arguments, std::vector<lcp::lps*>& cores );

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
void read_from_file( struct targs& thread_arguments, struct pargs& program_arguments );

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
void read_cores( std::vector<struct targs>& thread_arguments, struct pargs& program_arguments );


#endif
