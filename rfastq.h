#ifndef RFASTQ_H
#define RFASTQ_H

#include "args.h"
#include "utils.h"
#include "tpool.h"
#include "lps.h"
#include <htslib/kseq.h>
#include <sys/stat.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Reads multiple FASTQ files concurrently using a pool of threads.
 * 
 * This function processes a collection of FASTQ files by spawning threads 
 * based on the thread arguments (`targs`) and program settings (`pargs`). It 
 * manages the number of concurrent threads and ensures that each thread reads 
 * a FASTQ file using the `read_fastq` function, which operates on the provided 
 * thread arguments and shared program settings. 
 * 
 * @param genome_arguments A reference to a vector of `gargs` structures 
 *        representing the arguments specific to each genome.
 * @param program_arguments A constant reference to a `pargs` structure 
 *        representing the global program arguments.
 */
void read_fastqs(struct gargs *genome_arguments, struct pargs *program_arguments);

/**
 * @brief Processes a genome file to extract LCP cores using multiple threads.
 *
 * This function reads genomic sequences from a specified file and computes LCP cores 
 * for the sequences at a given LCP level and aggregates these cores into a shared array. 
 * The function tracks the total number of reads processed and their combined length. 
 * It ensures efficient and thread-safe handling of genomic data, leveraging parallel 
 * processing to enhance performance.
 *
 * @param args A reference to the `gargs` structure that contains the genome-specific 
 *        arguments, including the input FASTQ file name, the output data structures.
 */
void read_fastq(void *arg);

/**
 * @brief Processes a DNA sequence for both forward and reverse complement strands.
 *
 * This function processes a given DNA sequence in both its forward and reverse 
 * complement forms. It initializes the `lps` structure for both strands, deepens 
 * the `lps` structure based on the specified LCP level, and saves the processed 
 * result if the `write_lcpt` flag is set. After processing both strands, the 
 * allocated memory for the `lps` structures is freed.
 *
 * @param sequence A pointer to the DNA sequence to be processed.
 * @param seq_size The length of the DNA sequence.
 * @param capacity The pointer to the capacity value of the cores array.
 * @param genome_arguments Pointer to the genome arguments structure, which 
 *        contains settings such as the LCP level and whether to save results.
 * @param out The output file pointer to save the processed results.
 */
void process_read(char *sequence, size_t seq_size, uint64_t *capacity, struct gargs *genome_arguments, FILE *out);

#endif