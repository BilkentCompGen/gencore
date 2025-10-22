#include "rload.h"

pthread_mutex_t console_mutex_rload;

void read_lcpts(g_args_t *genome_args, p_args_t *program_args) {

    struct tpool *tm;

    tm = tpool_create(program_args->n_threads);

    for (int i=0; i<program_args->n_genomes; i++) {
        tpool_add_work(tm, read_lcpt, genome_args+i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void read_lcpt(void *arg) {

    g_args_t *genome_args = (g_args_t *)arg;

    if (genome_args->verbose) {
        pthread_mutex_lock(&console_mutex_rload);
        log1(INFO, "Thread ID: %ld started processing %s", pthread_self(), genome_args->inFileName);
        pthread_mutex_unlock(&console_mutex_rload);
    }

    // open fasta file
    FILE *in = fopen(genome_args->inFileName, "r");

    if (in == NULL) {
        log1(ERROR, "Error opening file %s", genome_args->inFileName);
        return;
    }

    // create file for writing cores
    FILE *out = NULL;

    if (genome_args->write_lcpt) {
        out = fopen(genome_args->inFileName, "wb");
        if (out == NULL) {
            log1(ERROR, "Error opening file for saving into file %s", genome_args->outFileName);
            return;
        }
    }

    // in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));

    // while(!isDone) {
    //     lcp::lps *str = new lcp::lps(in);
    //     str->deepen(program_args.lcpLevel);
        
    //     append(str, thread_arguments.cores);

    //     delete str;

    //     in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));
    // } 
    
    fclose(in);

    // end writing cores to file if user specified to do so
    if (genome_args->write_lcpt) {
        done(out);
        fclose(out);
    }

    // sort and filter the cores
    genSign(genome_args);

    // log ending of processing fasta
    if (genome_args->verbose) {
        log1(INFO, "Thread ID: %s ended processing %s, size: %ld", pthread_self(), genome_args->inFileName, genome_args->core_count);
    }
}