#ifndef SIMILARITY_METRICS_H
#define SIMILARITY_METRICS_H

#include <cstdint>
#include <vector>
#include <iterator>
#include <algorithm>
#include "logging.h"
#include "args.h"


/**
 * @brief Calculates the intersection and union sizes of LCP core sets from two thread arguments.
 * 
 * This function computes the intersection and union sizes between the `cores` and `counts` vectors 
 * in two thread-specific argument structures (`argument1` and `argument2`). It performs the calculations 
 * based on whether the program is operating in a set-based mode or a vector-based mode, as specified 
 * by the `program_arguments.type`.
 * 
 * @param argument1 A constant reference to the `targs` structure representing the first set of LCP cores 
 *        and counts for comparison.
 * @param argument2 A constant reference to the `targs` structure representing the second set of LCP cores 
 *        and counts for comparison.
 * @param program_arguments A constant reference to the `pargs` structure, which contains program-wide settings, 
 *        including the type of operation (set-based or vector-based).
 * @param interSize A reference to the variable where the computed size of the intersection will be stored.
 * @param unionSize A reference to the variable where the computed size of the union will be stored.
 * 
 * @details
 * - The function uses two iterators (`core1` and `core2`) to traverse the `cores` vectors of `argument1` 
 *   and `argument2` respectively, and similarly uses `count1` and `count2` for the corresponding `counts` vectors.
 * - If the LCP core in `argument1` is smaller than the core in `argument2`, the iterator for `argument1` advances, 
 *   and vice versa. When the cores match, the intersection size is updated, depending on the program type:
 *   - If `SET` mode is enabled, the intersection size is incremented by one.
 *   - If `VECTOR` mode is enabled, the intersection size is incremented by the minimum of the counts.
 * - The union size is calculated similarly, incrementing based on the maximum count when in `VECTOR` mode.
 * - After one vector is exhausted, any remaining elements in the other vector are added to the union size.
 * 
 * @note
 * - In `SET` mode, the function treats the `cores` vectors as sets and only counts unique cores for the 
 *   intersection and union.
 * - In `VECTOR` mode, the `counts` vectors are considered, and the union and intersection sizes are determined 
 *   by comparing the counts associated with each core.
 */
void calculateIntersectionAndUnionSizes( const struct targs& argument1, const struct targs& argument2, const struct pargs& program_arguments, size_t& interSize, size_t& unionSize );

/**
 * @brief Calculates the Jaccard similarity between two genomes.
 *
 * Computes the Jaccard similarity metric based on the intersection and union sizes
 * of hashed LCP cores from two genomes. This metric provides a measure of similarity
 * in terms of shared genomic features.
 *
 * @param intersectionSize The size of the intersection between the two sets of cores.
 * @param unionSize The size of the union between the two sets of cores.
 * @return The Jaccard similarity coefficient as a double.
 */
double calculateJaccardSimilarity( size_t intersectionSize, size_t unionSize );

/**
 * @brief Calculates the Dice similarity between two genomes.
 *
 * Computes the Dice similarity metric based on the intersection size of hashed LCP cores
 * and the sizes of individual core sets from two genomes. This metric is another measure
 * of similarity focusing on shared genomic features.
 *
 * @param intersectionSize The size of the intersection between the two sets of cores.
 * @param argument1 A constant reference to the `targs` structure representing the first set of LCP cores 
 *        and counts for comparison.
 * @param argument2 A constant reference to the `targs` structure representing the second set of LCP cores 
 *        and counts for comparison.
 * @param program_arguments A constant reference to the `pargs` structure, which contains program-wide settings, 
 *        including the type of operation (set-based or vector-based).
 * @return The Dice similarity coefficient as a double.
 */
double calculateDiceSimilarity( const size_t intersectionSize, const struct targs& argument1, const struct targs& argument2 );

#endif