#ifndef HELPER_H
#define HELPER_H

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include "logging.h"
#include "args.h"
#include "lps.h"

#ifndef BUFFERSIZE
#define BUFFERSIZE      100000
#endif


/**
 * @brief Generates the reverse complement of a DNA sequence.
 *
 * This function takes a DNA sequence as input, reverses it, and then replaces each
 * nucleotide with its complement (A <-> T, C <-> G). This operation is commonly used
 * in bioinformatics for DNA sequence analysis, especially when working with palindromic
 * sequences or when preparing to align sequence reads from both strands.
 *
 * @param sequence The DNA sequence to be reversed and complemented.
 * @return Returns a new string containing the reverse complement of the input sequence.
 */
bool reverseComplement( std::string& sequence );

/**
 * @brief Flattens a collection of locally parsed strings into a single vector of core labels.
 * 
 * This function processes a vector of pointers to `lcp::lps` (Locally Consisted Parsing) structures
 * and extracts the labels from their associated cores, appending them to the output vector `lcp_cores`.
 * The resulting flattened vector contains the labels of all `core` elements, maintaining the same 
 * order as they appear in the original structures.
 * 
 * @param strs A reference to a vector of pointers to `lcp::lps` structures, each containing a collection of cores.
 * @param lcp_cores A reference to a vector of 32-bit unsigned integers that will be populated with core labels.
 *                  This vector is modified in-place, and its capacity is pre-allocated based on the total number 
 *                  of cores across all locally parsed strings.
 */
void flatten(std::vector<lcp::lps*>& strs, std::vector<uint32_t>& lcp_cores);

/**
 * @brief Sorts the provided vector of hash values in ascending order.
 * 
 * This function modifies the input vector `hash_values` by sorting it in-place
 * using the standard library's `std::sort` algorithm. The resulting vector 
 * will contain the same values arranged in ascending order.
 * 
 * @param thread_arguments A reference to a vector of `targs` structures 
 *        representing the arguments specific to each thread which is needed for cores.
 * @param program_arguments A constant reference to a `pargs` structure 
 *        representing the global program arguments.
 */
void generateSignature( struct targs& thread_arguments, const struct pargs& program_arguments );

#endif