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
void calcUISize(const struct gargs *argument1, const struct gargs *argument2, uint64_t *interSize, uint64_t *unionSize);

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
 * @param avgLen   The average length (e.g., k-mer size).
 * @return The computed evolutionary distance.
 */
double calcEvolDist(double jaccardSim, double avgLen);

/**
 * @brief Applies the Jukes-Cantor correction to a Hamming distance.
 * 
 * This function converts a raw Hamming distance into a phylogenetic
 * distance using the Jukes-Cantor correction formula:
 * JC Distance = -(3/4) * log(1 - (3/4) * Hamming Distance)
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
 * @param genome_arguments Pointer to the array of genome arguments (`gargs`) 
 *                         containing genome data for comparison.
 * @param program_arguments Pointer to the program arguments (`pargs`) 
 *                          containing program-wide parameters such as the 
 *                          number of genomes and file prefix.
 */
void calcDistances(const struct gargs *genome_arguments, const struct pargs* program_arguments);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: LCP cores related functions
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Sorts the provided vector of hash values in ascending order.
 *
 * This function modifies the input vector `hash_values` by sorting it in-place
 * The resulting vector will contain the same values arranged in ascending order.
 *
 * @param genome_arguments A reference to a vector of `gargs` structures
 *        representing the arguments specific to each genome which is needed for cores.
 * @param apply_filter A boolena integer that decides whether to apply filter or not.
 */
void genSign(struct gargs *genome_arguments, int apply_filter);

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

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Cleanup
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Frees allocated memory for genome and program arguments.
 *
 * This function releases the memory allocated for the `cores` array in each 
 * genome argument within the `genome_arguments` structure. It also resets the 
 * length of the `cores` array (`cores_len`) to 0 to avoid dangling references. 
 * Finally, it frees the entire `genome_arguments` array.
 *
 * @param genome_arguments Pointer to the array of genome arguments (`gargs`) 
 *                         to be freed.
 * @param program_arguments Pointer to the program arguments (`pargs`) that 
 *                          contains the number of genomes.
 */
void free_args(struct gargs * genome_arguments, struct pargs * program_arguments);

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Heap Operations
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

/**
 * @brief Merges sorted arrays efficiently.
 *
 * This function merges the sorted arrays into single array efficiently using
 * min-heap. The heap is set to be fixed-sized, to the number of arrays to be exact.
 *
 * @param cores The 2D array that stores sorted arrays
 * @param size The array that stores sizes of the each array
 * @param file_count The number of array
 * @param min_cc A minimum core count threshold so the core will be stored.
 * @param max_cc A maximum core count threshold so the core will not be descared.
 * @param out_total_size The total number of element resulted in merging.
 * @param out_total_len The total length of the LCP cores.
 * @return Merged array
 */
uint64_t *merge_sorted_arrays(simple_core **cores, uint64_t *sizes, size_t file_count, uint32_t min_cc, uint32_t max_cc, uint64_t *out_total_size, double *out_total_len);

#endif
