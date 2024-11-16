#ifndef INIT_H
#define INIT_H

#include <sstream>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include "args.h"
#include "program_mode.h"
#include "logging.h"

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
 * @brief Displays the usage information for the program.
 *
 * This function prints the command-line usage instructions,
 * including available options and their descriptions. It is
 * typically called when the user requests help or when the
 * input arguments are invalid.
 */
void printUsage();

/**
 * @brief Parses command-line arguments.
 *
 * This function processes the command-line arguments passed to
 * the program and populates the thread_arguments and program_arguments
 * structures based on the provided options. It validates the input
 * and handles any errors, potentially displaying usage information.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @param thread_arguments A vector of targs structures to be filled
 *                        with thread-specific arguments based on
 *                        parsed command-line options.
 * @param program_arguments A pargs structure that stores global
 *                         program settings and parameters.
 *
 * @note This function may exit the program if the provided arguments
 *       are invalid or if there are missing required options.
 */
void parse( int argc, char **argv, std::vector<struct targs>& thread_arguments, struct pargs& program_arguments);

#endif