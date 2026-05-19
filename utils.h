#ifndef UTILS_H
#define UTILS_H

#include "args.h"
#include "lps.h"
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>


// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: ENUM Stringification
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Stringify the program mode enum `program_mode`.
 * 
 * @param m program mode
 * @return the string version of given parameter
 */
const char* mode2str(program_mode m);

/**
 * @brief Stringify the similarity calculation mode enum `sim_calculation_type`.
 * 
 * @param m sim calculation mode
 * @return the string version of given parameter
 */
const char* sct2str(sim_calculation_type m);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Similarity score calculation functions
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Calculates the intersection and union sizes of LCP core sets from two thread arguments.
 * 
 * This function computes the intersection and union sizes between the `cores` vectors 
 * in two thread-specific argument structures (`argument1` and `argument2`). It performs the calculations 
 * based on a set-based mode.
 * 
 * @param argument1 A constant reference to the `gargs` structure representing the first set of LCP cores 
 *        and counts for comparison.
 * @param argument2 A constant reference to the `gargs` structure representing the second set of LCP cores 
 *        and counts for comparison.
 * @param interSize A reference to the variable where the computed size of the intersection will be stored.
 * @param unionSize A reference to the variable where the computed size of the union will be stored.
 */
void calcUISize(const g_args_t *argument1, const g_args_t *argument2, uint64_t *interSize, uint64_t *unionSize);

/**
 * @brief Calculates the Jaccard similarity between two genomes.
 *
 * Computes the Jaccard similarity metric based on the intersection and union sizes
 * of hashed LCP cores from two genomes. This metric provides a measure of similarity
 * in terms of shared genomic features.
 *
 * @param interSize The size of the intersection between the two sets of cores.
 * @param unionSize The size of the union between the two sets of cores.
 * @return The Jaccard similarity coefficient as a double.
 */
double calcJaccardSim(uint64_t interSize, uint64_t unionSize);

/**
 * @brief Calculates the Dice similarity between two genomes.
 *
 * Computes the Dice similarity metric based on the intersection size of hashed LCP cores
 * and the sizes of individual core sets from two genomes. This metric is another measure
 * of similarity focusing on shared genomic features.
 *
 * @param interSize The size of the intersection between the two sets of cores.
 * @param size1 A size of the first set's core count.
 * @param size2 A size of the second set's core count.
 * @return The Dice similarity coefficient as a double.
 */
double calcDiceSim(uint64_t interSize, uint64_t size1, uint64_t size2);

/**
 * @brief Calculates the Hamming distance from Dice similarity.
 * 
 * This function computes the (approx) Hamming distance using the formula:
 * Hamming Distance = 1 - (Jaccard Similarity)^(1/avgLen)
 * 
 * @param jaccardSim The Jaccard similarity value (range: 0 to 1).
 * @param avgLen   The average length (e.g., k-mer size).
 * @return The computed Hamming distance.
 */
double calcHammDist(double jaccardSim, double avgLen);

/**
 * @brief Calculates the Hamming distance from Dice similarity.
 * 
 * This function computes the (approx) Hamming distance using the formula:
 * Evolutionary Distance = - 1 / kmerSize * ln( 2 * Jaccard Similarity / (1 +  Jaccard Similarity))
 * 
 * @param jaccardSim The Jaccard similarity value (range: 0 to 1).
 * @param avgLen   The average length (e.g., k-mer size, average lcp core size).
 * @return The computed evolutionary distance.
 */
double calcEvolDist(double jaccardSim, double avgLen);

/**
 * @brief Applies the Jukes-Cantor correction to a Hamming distance.
 * 
 * This function converts a raw Hamming distance into a phylogenetic
 * distance using the Jukes-Cantor correction formula:
 * JC Distance = -(3/4) * log(1 - (4/3) * Hamming Distance)
 * 
 * @param hammingDist The raw Hamming distance (range: 0 to 1).
 * @return The Jukes-Cantor corrected phylogenetic distance.
 */
double calcJukesCantorCor(double hammingDist);

/**
 * @brief Computes and writes distance matrices for genome comparisons.
 *
 * This function calculates pairwise distance matrices (Dice, Jaccard, and 
 * Jukes-Cantor) for the given genomes using their cores and lengths. It writes 
 * the resulting matrices to output files with filenames based on the program 
 * prefix, genome type, and LCP level.
 *
 * @param genome_args Pointer to the array of genome arguments (`gargs`) 
 *                         containing genome data for comparison.
 * @param program_args Pointer to the program arguments (`pargs`) 
 *                          containing program-wide parameters such as the 
 *                          number of genomes and file prefix.
 */
void calcDistances(const g_args_t *genome_args, const p_args_t* program_args);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: LCP cores related functions
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Compares two simple_core values.
 *
 * This function compares two values of type simple_core (uint64_t).
 * It is typically used as a comparison callback for sorting or searching
 * functions such as qsort() or bsearch().
 *
 * @param a Pointer to the first simple_core value.
 * @param b Pointer to the second simple_core value.
 * @return Negative if *a < *b, zero if equal, positive if *a > *b.
 */
int compare_simple_core(const void *a, const void *b);

/**
 * @brief Sorts the provided vector of hash values in ascending order.
 *
 * This function modifies the input vector `hash_values` by sorting it in-place
 * The resulting vector will contain the same values arranged in ascending order.
 *
 * It also records the time statistics for sorting and filtering.
 * 
 * @param args A reference to a vector of `gargs` structures
 *        representing the arguments specific to each genome which is needed for cores.
 */
void genSign(void *args);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: File I/O operations
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Saves the current state of the lc`lps` object to the specified output file.
 * 
 * This function writes a `false` value (indicating work is not complete),
 * followed by the serialized data of the lps object, to a binary file.
 *
 * @param out The output file pointer to save processed results.
 * @param str Pointer to the `lps` object that is being saved.
 */
void save(FILE *out, struct lps *str);

/**
 * @brief Marks the thread's processing as done by writing a boolean flag and size to 
 * the specified output file.
 * 
 * This function writes a `true` value to indicate the completion of the thread's work.
 *
 * @param out The output file pointer to save processed results.
 */
void done(FILE *out);

/**
 * @brief Makes an estimation about LCP cores that will be computed for give file.
 * 
 * This function returns number of cores that will be generated with LCP for all reads
 * in given fastq file.
 *
 * @param filename The input fastq filename.
 * @param lcp_level The lcp level that reads will be processed.
 * @return The expected number of LCP cores.
 */
uint64_t est_core_fq(const char *filename, int lcp_level);

/**
 * @brief Helper funct that checks the give filename is fq/fastq/fq.gz/fastq.gz
 * 
 * @param str File name to be checked
 * @return 1 if it fastq file, 0 if not
 */
int ends_with_fq(const char *str);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Logging
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Logs a formatted message with a timestamp and log level.
 *
 * This function prints a log message prefixed with the current timestamp and 
 * the specified log level (INFO, WARN, or ERROR). It uses a `printf`-style 
 * format string and additional arguments for the message content.
 *
 * @param level The log level (INFO, WARN, or ERROR) to categorize the log message.
 * @param format A `printf`-style format string for the log message.
 * @param ... Additional arguments for the format string.
 * @return Always returns 1 upon completion.
 */
int log1(LogLevel level, const char *format, ...);

/**
 * @brief Logs a formatted message with a timestamp and log level in thread-safe 
 * manner.
 *
 * This function prints a log message prefixed with the current timestamp and 
 * the specified log level (INFO, WARN, or ERROR). It uses a `printf`-style 
 * format string and additional arguments for the message content. It locks the
 * mutex while performing printing.
 *
 * @param level The log level (INFO, WARN, or ERROR) to categorize the log message.
 * @param mutex The mutex locked to be used for safe printing.
 * @param format A `printf`-style format string for the log message.
 * @param ... Additional arguments for the format string.
 * @return Always returns 1 upon completion.
 */
int log3(LogLevel level, pthread_mutex_t *mutex, const char *format, ...);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Cleanup
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Frees allocated memory for genome and program arguments.
 *
 * This function releases the memory allocated for the `cores` array in each 
 * genome argument within the `genome_args` structure. It also resets the 
 * length of the `cores` array (`core_count`) to 0 to avoid dangling references. 
 * Finally, it frees the entire `genome_args` array.
 *
 * @param genome_args Pointer to the array of genome arguments (`gargs`) 
 *                         to be freed.
 * @param program_args Pointer to the program arguments (`pargs`) that 
 *                          contains the number of genomes.
 */
void free_args(g_args_t * genome_args, p_args_t * program_args);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Heap Operations
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Merges multiple sorted arrays using a fixed-size min-heap.
 *
 * This function efficiently merges multiple sorted arrays of simple_core
 * elements into a single sorted array. A fixed-size min-heap, whose size
 * equals the number of arrays, is used to maintain the smallest current
 * elements across all arrays.
 * 
 * It also records merging time and logging purposes.
 *
 * @param cores         A 2D array containing the sorted input arrays.
 * @param sizes         An array containing the size of each input array.
 * @param file_count    The number of input arrays.
 * @param genome_args   Pointer to genome-level configuration and thresholds.
 *
 * @note The function applies core count thresholds specified in genome_args
 *       (e.g., minimum and maximum allowed values) when merging.
 */
void merge_sorted_arrays(simple_core **cores, uint64_t *sizes, uint64_t file_count, g_args_t *genome_args);

/**
 * @brief Merges sorted arrays produced by multiple threads.
 *
 * This function merges sorted simple_core arrays generated by parallel
 * worker threads into a single output array.
 *
 * @param args   Array of thread argument structures containing metadata
 *               and pointers to thread-local arrays.
 * @param n_args Number of thread argument structures.
 * @param cores  Pointer to the output array that will store merged results.
 *
 * @return The total number of merged simple_core elements.
 */
uint64_t merge_thread_arrays(fqw_args_t *args, int n_args, simple_core **cores);

uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed);

#endif
