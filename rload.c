#include "rload.h"

pthread_mutex_t console_mutex_rload;

void read_lcpts(g_args_t *genome_args, p_args_t *program_args) {

    struct tpool *tm;

    tm = tpool_create(program_args->n_threads);

    for (int i = 0; i < program_args->n_genomes; i++) {
        tpool_add_work(tm, read_lcpt, genome_args+i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);

    calcDistancesFrequencySweep(genome_args, program_args);
}

void read_lcpt(void *arg) {

    g_args_t *genome_args = (g_args_t *)arg;

    if (genome_args->verbose) {
        pthread_mutex_lock(&console_mutex_rload);
        log1(INFO, "Thread ID: %ld started processing %s", pthread_self(), genome_args->inFileName);
        pthread_mutex_unlock(&console_mutex_rload);
    }

    // open binary file
    FILE *in = fopen(genome_args->inFileName, "rb");

    if (in == NULL) {
        log1(ERROR, "Error opening file %s", genome_args->inFileName);
        return;
    }

    uint64_t core_count = 0;

    if (fread(&core_count, sizeof(uint64_t), 1, in) != 1) {
        log1(ERROR, "Error reading core count from file %s", genome_args->inFileName);
        fclose(in);
        return;
    }

    genome_args->result.count = core_count;
    genome_args->result.cores = (simple_core *)malloc(sizeof(simple_core) * core_count);

    if (genome_args->result.cores == NULL && core_count > 0) {
        log1(ERROR, "Memory allocation failed");
        fclose(in);
        return;
    }

    size_t read_count = fread(
        genome_args->result.cores,
        sizeof(simple_core),
        core_count,
        in
    );

    if (read_count != core_count) {
        log1(ERROR, "Expected %lu cores, but read %zu", core_count, read_count);

        free(genome_args->result.cores);
        genome_args->result.cores = NULL;
        genome_args->result.count = 0;

        fclose(in);
        return;
    }

    fclose(in);

    // sort
    time_t start, end;
    time(&start);

    sort_u64_radix(genome_args->result.cores, genome_args->result.count);

    time(&end);
    genome_args->time_stats.sorting = difftime(end, start);

    // log ending of processing fasta
    if (genome_args->verbose) {
        log1(INFO, "Thread ID: %ld ended processing %s, size: %ld", pthread_self(), genome_args->inFileName, genome_args->result.count);
    }
}

void calcDistancesFrequencySweep(g_args_t *genome_args, p_args_t *program_args) {
    for (uint32_t min_cc = 0; min_cc <= 200; min_cc += 4) {

        log1(INFO, "Running distance calculation for min_cc=%u", min_cc);

        g_args_t *filtered_args = calloc(program_args->n_genomes, sizeof(g_args_t));

        if (filtered_args == NULL) {
            log1(ERROR, "Memory allocation failed for filtered genome args");
            return;
        }

        int failed = 0;

        for (int i = 0; i < program_args->n_genomes; i++) {
            filtered_args[i] = genome_args[i];

            uint32_t max_cc = genome_args[i].max_cc;

            if (!genome_args[i].apply_filter) {
                max_cc = UINT32_MAX;
            }

            if (build_filtered_result(&genome_args[i], &filtered_args[i], min_cc, max_cc) != 0) {
                failed = 1;
                break;
            }
        }

        if (genome_args->verbose) {
            for (int i = 0; i < program_args->n_genomes; i++) {   
                log1(INFO, "Processing ended for %s, size: %ld", genome_args[i].inFileName, filtered_args->result.count);
            }
        }

        if (!failed) {
            char *temp_prefix = program_args->prefix;
            char filename_buffer[256];
            if (snprintf(filename_buffer, 256, "%s.cc%u", temp_prefix, min_cc) < 0) {
                log1(ERROR, "Filename buffer for dice overflow.");
                exit(EXIT_FAILURE);
            }
            program_args->prefix = temp_prefix;

            calcDistances(filtered_args, program_args);
            
            program_args->prefix = temp_prefix;
        }

        for (int i = 0; i < program_args->n_genomes; i++) {
            free(filtered_args[i].result.cores);
            filtered_args[i].result.cores = NULL;
            filtered_args[i].result.count = 0;
        }

        free(filtered_args);

        if (failed) {
            log1(ERROR, "Stopping frequency sweep because filtering failed");
            return;
        }
    }
}