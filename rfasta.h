#ifndef RFASTA_H
#define RFASTA_H

#include "args.h"
#include "utils.h"
#include "tpool.h"
#include "lps.h"
#include <stdint.h>

#define INITIAL_SEQUENCE_SIZE 300000000

/**
 * @brief Reads multiple FASTA files concurrently using a pool of threads.
 * 
 * This function processes a collection of FASTA files by spawning threads 
 * based on the genome arguments (`gargs`) and program settings (`pargs`). It 
 * manages the number of concurrent threads and ensures that each thread reads 
 * a FASTA file using the `read_fasta` function, which operates on the provided 
 * thread arguments and shared program settings. 
 * 
 * @param genome_arguments A reference to a array of `gargs` structures 
 *        representing the arguments specific to each genome.
 * @param program_arguments A constant reference to a `pargs` structure 
 *        representing the global program arguments.
 */
void read_fastas(struct gargs *genome_arguments, struct pargs *program_arguments);

/**
 * @brief Reads a FASTA file and processes its sequences using the LCP (Locally 
 * Consistent Parsing) method.
 * 
 * This function is responsible for reading a FASTA file specified in the `genome_arguments`, 
 * processing each sequence using LCP technique, and storing the results. The function manages 
 * logging for verbose output, tracks the size of processed sequences, and handles thread-safe 
 * operations, as it is designed to be run in a multithreaded environment.
 * 
 * @param args A reference to the `gargs` structure that contains the genome-specific 
 *        arguments, including the input FASTA file name, the output data structures.
 */
void read_fasta(void *arg);

/**
 * @brief Processes a DNA sequence with LCP technique and extracts cores.
 *
 * This function processes a given DNA sequence and initializes the `lps` structure 
 * for as strand, deepens the `lps` structure based on the specified LCP level, 
 * and saves the processed result if the `write_lcpt` flag is set.
 *
 * @param sequence A pointer to the DNA sequence to be processed.
 * @param seq_size The length of the DNA sequence.
 * @param capacity The pointer to the capacity value of the cores array.
 * @param genome_arguments Pointer to the genome arguments structure, which 
 *        contains settings such as the LCP level and whether to save results.
 * @param out The output file pointer to save the processed results.
 */
void process_chrom(char *sequence, size_t seq_size, uint64_t *capacity, struct gargs *genome_arguments, FILE *out);

#endif