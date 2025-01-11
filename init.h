#ifndef INIT_H
#define INIT_H

#include "args.h"
#include "utils.h" // logging
#include <stdio.h>
#include <errno.h> // errno
#include <limits.h> // UINT32_MAX
#include <getopt.h>

#ifndef THREAD_NUMBER
#define THREAD_NUMBER 8
#endif

#ifndef VERBOSE 
#define VERBOSE false
#endif

#ifndef PREFIX 
#define PREFIX "gc"
#endif

#ifndef COMPRESSION_RATIO
#define COMPRESSION_RATIO 4
#endif

/**
 * @brief Parses command-line arguments.
 *
 * This function processes the command-line arguments passed to
 * the program and populates the genome_arguments and program_arguments
 * structures based on the provided options. It validates the input
 * and handles any errors, potentially displaying usage information.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @param genome_arguments A pointer to array of gargs structures to be filled
 *                        with genome-specific arguments based on
 *                        parsed command-line options.
 * @param program_arguments A pointer to pargs structure that stores global
 *                         program settings and parameters.
 *
 * @note This function may exit the program if the provided arguments
 *       are invalid or if there are missing required options.
 */
void parse(int argc, char **argv, struct gargs **genome_arguments, struct pargs *program_arguments);

#endif