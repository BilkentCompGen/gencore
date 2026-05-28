#include "init.h"

void printUsage() {
    printf("Usage: ./gencore [PROGRAM] [OPTIONS]\n");
    printf("[PROGRAM]: \n");
    printf("\tfa:   Processing assembled genomes.\n");
    printf("\tfq:   Processing genomes' reads.\n");
    printf("\tld:   Processing precomputed cores.\n");
}

void printFaUsage() {
    printf("Usage: ./gencore fa [OPTIONS]\n");
    printf("Options:\n");
    printf("\t-i [filename]   The file contains filenames of genomes.\n");
    printf("\t-l [num]        Lcp-level. [Default: %d]\n", DEFAULT_LCP_LEVEL);
    printf("\t-t [num]        Number of threads. [Default: %d]\n", DEFAULT_THREAD_NUMBER);
    printf("\t--min-cc [num]  Minimum frequency (core count) for a core. [Default: %d]\n", DEFAULT_FA_MIN_CC);
    printf("\t--max-cc [num]  Maximum frequency (core count) for a core. [Default: %d]\n", DEFAULT_FA_MAX_CC);
    printf("\t[--set|--vec]   Distances based or set or vector of cores. [Default: %s]\n", sct2str(DEFAULT_SIM_CALC_MODE));
    printf("\t-o [filename]   Store cores.\n");
    printf("\t-p [prefix]     Prefix for the results. [Default: %s]\n", DEFAULT_PREFIX);
    printf("\t-s [filename]   Set short names of input files. Default is first 10 characters of input file names.\n");
    printf("\t-v              Verbose. [Default: %d]\n\n", DEFAULT_VERBOSE);
}

void printFqUsage() {
    printf("Usage: ./gencore fq [OPTIONS]\n");
    printf("Options:\n");
    printf("\t-i [filename]   The file contains filenames of genomes.\n");
    printf("\t-l [num]        Lcp-level. [Default: %d]\n", DEFAULT_LCP_LEVEL);
    printf("\t-t [num]        Number of threads. [Default: %d]\n", DEFAULT_THREAD_NUMBER);
    printf("\t-r [num]        Number of reader threads. [Default: %d]\n", DEFAULT_FQ_READER_NUMBER);
    printf("\t--min-cc [num]  Minimum frequency (core count) for a core. [Default: %d]\n", DEFAULT_FQ_MIN_CC);
    printf("\t--max-cc [num]  Maximum frequency (core count) for a core. [Default: %d]\n", DEFAULT_FA_MAX_CC);
    printf("\t[--set|--vec]   Distances based or set or vector of cores. [Default: %s]\n", sct2str(DEFAULT_SIM_CALC_MODE));
    printf("\t-o [filename]   Store cores.\n");
    printf("\t-p [prefix]     Prefix for the results. [Default: %s]\n", DEFAULT_PREFIX);
    printf("\t-s [filename]   Set short names of input files. Default is first 10 characters of input file names.\n");
    printf("\t-v              Verbose. [Default: %d]\n\n", DEFAULT_VERBOSE);
}

void printUsage2(program_mode_t mode) {
    switch(mode) {
    case PROGRAM_MODE_FA:
        printFaUsage();
        break;
    case PROGRAM_MODE_FQ:
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
 
        uint64_t len = strlen(buffer);

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

void free_targs(g_args_t **genome_args, p_args_t *program_args) {
    for (int i = 0; i < program_args->n_genomes; i++) {
        // clean inFileName
        if ((*genome_args)[i].inFileName != NULL)
            free((*genome_args)[i].inFileName);
        (*genome_args)[i].inFileName = NULL;
        // clean shortName  
        if ((*genome_args)[i].shortName != NULL)
            free((*genome_args)[i].shortName);
        (*genome_args)[i].shortName = NULL;
        // clean outFileName
        if ((*genome_args)[i].outFileName != NULL)
            free((*genome_args)[i].outFileName);
        (*genome_args)[i].outFileName = NULL;
    }
}

void parse(int argc, char **argv, g_args_t **genome_args, p_args_t *program_args) {

    if (argc < 2) {
        printUsage();
        exit(1);
    }

    int apply_filter;
    uint32_t min_cc;
    uint32_t max_cc;

    if (strcmp(argv[1], "fa") == 0) {
        program_args->mode = PROGRAM_MODE_FA;
        min_cc = DEFAULT_FA_MIN_CC;
        max_cc = DEFAULT_FA_MAX_CC;
        apply_filter = 0;
    } else if (strcmp(argv[1], "fq") == 0) {
        program_args->mode = PROGRAM_MODE_FQ;
        min_cc = DEFAULT_FQ_MIN_CC;
        max_cc = DEFAULT_FQ_MAX_CC;
        apply_filter = 1;
    } else if (strcmp(argv[1], "ld") == 0) {
        program_args->mode = PROGRAM_MODE_LOAD;
        min_cc = 0;
        max_cc = UINT32_MAX;
        apply_filter = 0;
    } else {
        log1(ERROR, "Invalid program mode '%s'", argv[1]);
        printUsage();
        exit(EXIT_FAILURE);
    }

    // set program arguments with their default values
    program_args->n_threads = DEFAULT_THREAD_NUMBER;
    program_args->n_readers = DEFAULT_FQ_READER_NUMBER;
    program_args->prefix = DEFAULT_PREFIX;
    program_args->n_genomes = 0;

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
    sim_calculation_type_t sct = DEFAULT_SIM_CALC_MODE;
    int lcp_level = DEFAULT_LCP_LEVEL;
    int core_span = DEFAULT_CORE_OUTSPAN;
    int write_lcpt = DEFAULT_WRITE_LCP_CORES;
    int verbose = DEFAULT_VERBOSE;

    int opt;
    int long_index;
    char *endptr;

    // Parsing options
    while ((opt = getopt_long(argc, argv, "i:l:e:t:r:o:p:s:v", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'i':
                filename_inputs = optarg;
                break;
            case 'l':
                lcp_level = atoi(optarg);
                break;
            case 'e':
                core_span = atoi(optarg);
                break;
            case 't':
                program_args->n_threads = atoi(optarg);
                break;
            case 'r':
                program_args->n_readers = atoi(optarg);
                break;
            case 'o':
                filename_outputs = optarg;
                write_lcpt = 1;
                break;
            case 'p':
                program_args->prefix = optarg;
                break;
            case 's':
                filename_names = optarg;
                break;
            case 'v':
                verbose = 1;
                break;
            case 1: // --min-cc
                min_cc = (uint32_t)strtol(optarg, &endptr, 10);
                apply_filter = 1;
                break;
            case 2: // --min-cc-file
                filename_min_cc = optarg;
                apply_filter = 1;
                break;
            case 3: // --max-cc
                max_cc = (uint32_t)strtol(optarg, &endptr, 10);
                apply_filter = 1;
                break;
            case 4: // --max-cc-file
                filename_max_cc = optarg;
                apply_filter = 1;
                break;
            case 5: // --set
                sct = SIM_CALC_SET;
                break;
            case 6: // --vec
                sct = SIM_CALC_SET;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    if (filename_inputs == NULL) {
        log1(ERROR, "Please provide input files.");
        printUsage2(program_args->mode);
        exit(EXIT_FAILURE);
    }

    program_args->n_genomes = get_line_count(filename_inputs);
    program_args->sct = sct;
    program_args->lcp_level = lcp_level;
    program_args->core_span = core_span;
    program_args->write_lcpt = write_lcpt;
    program_args->verbose = verbose;

    if (program_args->n_genomes == -1) {
        exit(EXIT_FAILURE);
    }

    (*genome_args) = (g_args_t *)malloc(sizeof(g_args_t) * program_args->n_genomes);
    if ((*genome_args) == NULL) {
        log1(ERROR, "Memory allocation failed for genome arguments.");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < program_args->n_genomes; i++) {
        (*genome_args)[i].inFileName = NULL;
        (*genome_args)[i].shortName = NULL;
        (*genome_args)[i].outFileName = NULL;
        (*genome_args)[i].apply_filter = apply_filter;
        (*genome_args)[i].lcp_level = lcp_level;
        (*genome_args)[i].write_lcpt = write_lcpt;
        (*genome_args)[i].verbose = verbose;
        (*genome_args)[i].min_cc = min_cc;
        (*genome_args)[i].max_cc = max_cc;
        (*genome_args)[i].sct = sct;
        memset(&((*genome_args)[i].result), 0, sizeof(core_result_t));
        memset(&((*genome_args)[i].time_stats), 0, sizeof(time_stats_t));
    }

    // program_args->n_threads = program_args->n_threads < program_args->n_genomes ? program_args->n_threads : program_args->n_genomes;

    // check filename_inputs
    if (filename_inputs != NULL) {

        FILE *file = fopen(filename_inputs, "r");
        if (file == NULL) {
            log1(ERROR, "Could not open file: %s", filename_inputs);
            exit(EXIT_FAILURE);
        }

        char buffer[1024];
 
        for (int i = 0; i < program_args->n_genomes; i++) {
            if (read_line(file, buffer, &((*genome_args)[i].inFileName)) == -1) {
                free_targs(genome_args, program_args);
                free(*genome_args);
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
 
        for (int i = 0; i < program_args->n_genomes; i++) {
            if (read_line_uint32(file, buffer, &((*genome_args)[i].min_cc)) == -1) {
                free_targs(genome_args, program_args);
                free(*genome_args);
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
 
        for (int i = 0; i < program_args->n_genomes; i++) {
            if (read_line_uint32(file, buffer, &((*genome_args)[i].max_cc)) == -1) {
                free_targs(genome_args, program_args);
                free(*genome_args);
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
 
        for (int i = 0; i < program_args->n_genomes; i++) {
            if (read_line(file, buffer, &((*genome_args)[i].shortName)) == -1) {
                free_targs(genome_args, program_args);
                free(*genome_args);
                fclose(file);
                exit(EXIT_FAILURE);
            }
            if (strlen((*genome_args)[i].shortName) > 10)
                (*genome_args)[i].shortName[10] = '\0';
        }

        fclose(file);
    } else {
        for (int i = 0; i < program_args->n_genomes; i++) {
            (*genome_args)[i].shortName = strdup((*genome_args)[i].inFileName);
            if (strlen((*genome_args)[i].shortName) > 10)
                (*genome_args)[i].shortName[10] = '\0';
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
 
        for (int i = 0; i < program_args->n_genomes; i++) {
            if (read_line(file, buffer, &((*genome_args)[i].outFileName)) == -1) {
                free_targs(genome_args, program_args);
                free(*genome_args);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }

        fclose(file);
    }

    // log parameters
    log1(INFO, "%s, threads: %d, LCP: %d, span: %d, mode: %s, prefix: %s", mode2str(program_args->mode), program_args->n_threads, program_args->lcp_level, program_args->core_span, sct2str(program_args->sct), program_args->prefix);

    if ((*genome_args)[0].write_lcpt) { 
        log1(INFO, "Program will write cores to files.");
    }

    if ((*genome_args)[0].verbose) {
        for (int i = 0; i < program_args->n_genomes; i++) {
            if ((*genome_args)[i].apply_filter) {
                log1(INFO, "in: %s, short: %s, out: %s, min-cc: %ld, max-cc: %ld", (*genome_args)[i].inFileName, (*genome_args)[i].shortName, (*genome_args)[i].outFileName, (*genome_args)[i].min_cc, (*genome_args)[i].max_cc);
            } else {
                log1(INFO, "in: %s, short: %s, out: %s", (*genome_args)[i].inFileName, (*genome_args)[i].shortName, (*genome_args)[i].outFileName);
            }
        }
    }
}