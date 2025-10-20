#ifndef INIT_H
#define INIT_H

#include "args.h"
#include "utils.h" // logging
#include <stdio.h>
#include <errno.h> // errno
#include <limits.h> // UINT32_MAX
#include <getopt.h>

/**
 * @brief Parses command-line arguments.
 *
 * This function processes the command-line arguments passed to
 * the program and populates the genome_args and program_args
 * structures based on the provided options. It validates the input
 * and handles any errors, potentially displaying usage information.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @param genome_args A pointer to array of gargs structures to be filled
 *                        with genome-specific arguments based on
 *                        parsed command-line options.
 * @param program_args A pointer to pargs structure that stores global
 *                         program settings and parameters.
 *
 * @note This function may exit the program if the provided arguments
 *       are invalid or if there are missing required options.
 */
void parse(int argc, char **argv, g_args_t **genome_args, p_args_t *program_args);

#endif