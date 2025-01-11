#include "init.h"

void printUsage() {
    printf("Usage: ./gencore [PROGRAM] [OPTIONS]\n\n");
    printf("[PROGRAM]: \n");
    printf("\tfa:   Processing assembled genomes.\n");
    printf("\tfq:   Processing genomes' reads.\n");
    printf("\tld:   Processing precomputed cores.\n");
}

void printFaUsage() {
    printf("Usage: ./gencore fa [OPTIONS]\n\n");
    printf("Options:\n");
    printf("\t-i [filename]   The file contains filenames of genomes.\n\n");
    printf("\t-l [num]        Lcp-level. [Default: 4]\n\n");
    printf("\t-t [num]        Number of threads. [Default: 8]\n\n");
    printf("\t--min-cc [num]  Minimum frequency (core count) for a core. [Default: 1]\n\n");
    printf("\t--max-cc [num]  Maximum frequency (core count) for a core. [Default: UINT64_MAX]\n\n");
    printf("\t[--set|--vec]   Distances based or set or vector of cores. [Default: set]\n\n");
    printf("\t-o [filename]   Store cores.\n\n");
    printf("\t-p [prefix]     Prefix for the results. [Default: gc]\n\n");
    printf("\t-s [filename]   Set short names of input files. Default is first 10 characters of input file names.\n\n");
    printf("\t-v              Verbose. [Default: false]\n\n");
}

void printFqUsage() {
    printf("Usage: ./gencore fa [OPTIONS]\n\n");
    printf("Options:\n");
    printf("\t-i [filename]   The file contains filenames of genomes.\n\n");
    printf("\t-l [num]        Lcp-level. [Default: 4]\n\n");
    printf("\t-t [num]        Number of threads. [Default: 8]\n\n");
    printf("\t--min-cc [num]  Minimum frequency (core count) for a core. [Default: 15]\n\n");
    printf("\t--max-cc [num]  Maximum frequency (core count) for a core. [Default: 256]\n\n");
    printf("\t[--set|--vec]   Distances based or set or vector of cores. [Default: set]\n\n");
    printf("\t-o [filename]   Store cores.\n\n");
    printf("\t-p [prefix]     Prefix for the results. [Default: gc]\n\n");
    printf("\t-s [filename]   Set short names of input files. Default is first 10 characters of input file names.\n\n");
    printf("\t-v              Verbose. [Default: false]\n\n");
}

void printUsage2(program_mode mode) {
    switch(mode) {
    case FA:
        printFaUsage();
        break;
    case FQ:
        printFqUsage();
        break;
    default:
        break;
    }
}

int get_line_count(const char *filename) {
    int line_count = 0;

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        log1(ERROR, "Could not open file: %s", filename);
        return -1;
    }

    char buffer[1024];

    while (fgets(buffer, sizeof(buffer), file)) {
        line_count++;
    }
    
    fclose(file);
    return line_count;
}

int read_line_uint32(FILE *file, char buffer[1024], uint32_t *val) {

    if (fgets(buffer, 1024, file)) {

        char *endptr;
        errno = 0;
        unsigned long value = strtoul(buffer, &endptr, 10); // base 10 conversion

        if (errno == ERANGE || value > UINT32_MAX) {
            log1(ERROR, "Value out of range for uint32_t.");
            exit(EXIT_FAILURE);
        }
        if (endptr == buffer || *endptr != '\0') {
            log1(ERROR, "Invalid numeric string: %s", buffer);
            exit(EXIT_FAILURE);
        }

        *val = (uint32_t)value;
        return 1;
    }
    return 0;
}

int read_line(FILE *file, char buffer[1024], char **result) {

    if (fgets(buffer, 1024, file)) {
 
        size_t len = strlen(buffer);

        if (len > 0 && buffer[len - 1] == '\n') {
            buffer[len - 1] = '\0';
        }

        *result = strdup(buffer); // might not work on different systems

        if (*result == NULL) {
            log1(ERROR, "Memory allocation failed.");
            return -1;
        }

        return 0;
    }

    // if we reach here, the line_num was out of range
    log1(ERROR, "Out of range in file.");
    return -1;
}

void free_targs(struct gargs **genome_arguments, struct pargs *program_arguments) {
    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        // clean inFileName
        if ((*genome_arguments)[i].inFileName != NULL)
            free((*genome_arguments)[i].inFileName);
        (*genome_arguments)[i].inFileName = NULL;
        // clean shortName  
        if ((*genome_arguments)[i].shortName != NULL)
            free((*genome_arguments)[i].shortName);
        (*genome_arguments)[i].shortName = NULL;
        // clean outFileName
        if ((*genome_arguments)[i].outFileName != NULL)
            free((*genome_arguments)[i].outFileName);
        (*genome_arguments)[i].outFileName = NULL;
    }
}

void parse(int argc, char **argv, struct gargs **genome_arguments, struct pargs *program_arguments) {

    if (argc < 2) {
        printUsage();
        exit(1);
    }

    uint32_t min_cc;
    uint32_t max_cc;

    if (strcmp(argv[1], "fa") == 0) {
        program_arguments->mode = FA;
        min_cc = 0;
        max_cc = UINT32_MAX;
    } else if (strcmp(argv[1], "fq") == 0) {
        program_arguments->mode = FQ;
        min_cc = 15;
        max_cc = 256;
    } else if (strcmp(argv[1], "ld") == 0) {
        program_arguments->mode = LOAD;
        min_cc = 0;
        max_cc = UINT32_MAX;
    } else {
        log1(ERROR, "Invalid program mode '%s'", argv[1]);
        printUsage2(program_arguments->mode);
        exit(EXIT_FAILURE);
    }

    // set program arguments with their default values
    program_arguments->thread_number = 8;
    program_arguments->prefix = "gc";
    program_arguments->number_of_genomes = 0;

    struct option long_options[] = {
        {"min-cc", required_argument, NULL, 1},
        {"min-cc-file", required_argument, NULL, 2},
        {"max-cc", required_argument, NULL, 3},
        {"max-cc-file", required_argument, NULL, 4},
        {"set", no_argument, NULL, 5},
        {"vec", no_argument, NULL, 6},
        {NULL, 0, NULL, 0}
    };

    char *filename_min_cc = NULL;
    char *filename_max_cc = NULL;
    char *filename_inputs = NULL;
    char *filename_names = NULL;
    char *filename_outputs = NULL;
    sim_calculation_type sct = SET;
    int lcp_level = 4;
    int write_lcpt = 0;
    int verbose = 0;

    int opt;
    int long_index;
    char *endptr;

    // Parsing options
    while ((opt = getopt_long(argc, argv, "i:l:t:o:p:s:v", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'i':
                filename_inputs = optarg;
                break;
            case 'l':
                lcp_level = atoi(optarg);
                break;
            case 't':
                program_arguments->thread_number = atoi(optarg);
                break;
            case 'o':
                filename_outputs = optarg;
                write_lcpt = 1;
                break;
            case 'p':
                program_arguments->prefix = optarg;
                break;
            case 's':
                filename_names = optarg;
                break;
            case 'v':
                verbose = 1;
                break;
            case 1: // --min-cc
                min_cc = (uint32_t)strtol(optarg, &endptr, 10);
                break;
            case 2: // --min-cc-file
                filename_min_cc = optarg;
                break;
            case 3: // --max-cc
                max_cc = (uint32_t)strtol(optarg, &endptr, 10);
                break;
            case 4: // --max-cc-file
                filename_max_cc = optarg;
                break;
            case 5: // --set
                sct = SET;
                break;
            case 6: // --vec
                sct = VECTOR;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    if (filename_inputs == NULL) {
        log1(ERROR, "Please provide input function");
        printUsage2(program_arguments->mode);
        exit(EXIT_FAILURE);
    }

    program_arguments->number_of_genomes = get_line_count(filename_inputs);

    if (program_arguments->number_of_genomes == -1) {
        exit(EXIT_FAILURE);
    }

    (*genome_arguments) = (struct gargs*)malloc(program_arguments->number_of_genomes * sizeof(struct gargs));
    if ((*genome_arguments) == NULL) {
        log1(ERROR, "Memory allocation failed for genome arguments.");
        exit(EXIT_FAILURE);
    }

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        (*genome_arguments)[i].min_cc = min_cc;
        (*genome_arguments)[i].max_cc = max_cc;
        (*genome_arguments)[i].inFileName = NULL;
        (*genome_arguments)[i].shortName = NULL;
        (*genome_arguments)[i].outFileName = NULL;
        (*genome_arguments)[i].cores_len = 0;
        (*genome_arguments)[i].cores = NULL;
        (*genome_arguments)[i].total_len = 0.0;
        (*genome_arguments)[i].sct = sct;
        (*genome_arguments)[i].lcp_level = lcp_level;
        (*genome_arguments)[i].write_lcpt = write_lcpt;
        (*genome_arguments)[i].verbose = verbose;
    }

    program_arguments->thread_number = program_arguments->thread_number < program_arguments->number_of_genomes ? program_arguments->thread_number : program_arguments->number_of_genomes;

    // check filename_inputs
    if (filename_inputs != NULL) {

        FILE *file = fopen(filename_inputs, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_inputs);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            if (read_line(file, buffer, &((*genome_arguments)[i].inFileName)) == -1) {
                free_targs(genome_arguments, program_arguments);
                free(*genome_arguments);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }

        fclose(file);
    }
    
    // check filename_min_cc
    if (filename_min_cc != NULL) {

        FILE *file = fopen(filename_min_cc, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_min_cc);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            if (read_line_uint32(file, buffer, &((*genome_arguments)[i].min_cc)) == -1) {
                free_targs(genome_arguments, program_arguments);
                free(*genome_arguments);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }

        fclose(file);
    }
    
    // check filename_max_cc
    if (filename_max_cc != NULL) {

        FILE *file = fopen(filename_max_cc, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_max_cc);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            if (read_line_uint32(file, buffer, &((*genome_arguments)[i].max_cc)) == -1) {
                free_targs(genome_arguments, program_arguments);
                free(*genome_arguments);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }

        fclose(file);
    }

    // check filename_names
    if (filename_names != NULL) {

        FILE *file = fopen(filename_names, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_names);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            if (read_line(file, buffer, &((*genome_arguments)[i].shortName)) == -1) {
                free_targs(genome_arguments, program_arguments);
                free(*genome_arguments);
                fclose(file);
                exit(EXIT_FAILURE);
            }
            (*genome_arguments)[i].shortName[10] = '\0';
        }

        fclose(file);
    } else {
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            (*genome_arguments)[i].shortName = strdup((*genome_arguments)[i].inFileName);
            (*genome_arguments)[i].shortName[10] = '\0';
        }
    }

    // check filename_outputs
    if (filename_outputs != NULL) {

        FILE *file = fopen(filename_outputs, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_outputs);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            if (read_line(file, buffer, &((*genome_arguments)[i].outFileName)) == -1) {
                free_targs(genome_arguments, program_arguments);
                free(*genome_arguments);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }

        fclose(file);
    }

    // log parameters
    if (strcmp(argv[1], "fa") == 0) {
        log1(INFO, "Program mode: FA");
    } else if (strcmp(argv[1], "fq") == 0) {
        log1(INFO, "Program mode: FQ");
    } else if (strcmp(argv[1], "bam") == 0) {
        log1(INFO, "Program mode: BAM");
    } else if (strcmp(argv[1], "ld") == 0) {
        log1(INFO, "Program mode: LOAD");
    }

    log1(INFO, "Thread number: %d", program_arguments->thread_number);
    log1(INFO, "Prefix: %s", program_arguments->prefix);
    log1(INFO, "LCP level: %d", (*genome_arguments)[0].lcp_level);
    log1(INFO, "Distance calculation mode: %s", ((*genome_arguments)[0].sct == SET ? "set" : "vector"));

    if ((*genome_arguments)[0].write_lcpt) { 
        log1(INFO, "Program will write cores to files.");
    }

    if ((*genome_arguments)[0].verbose) {
        for (int i=0; i<program_arguments->number_of_genomes; i++) {
            log1(INFO, "inFileName: %s, shortName: %s, outFileName: %s, min-cc %ld, max-cc: %ld", (*genome_arguments)[i].inFileName, (*genome_arguments)[i].shortName, (*genome_arguments)[i].outFileName, (*genome_arguments)[i].min_cc, (*genome_arguments)[i].max_cc);
        }
    }
}