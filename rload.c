#include "rload.h"

pthread_mutex_t console_mutex_rload;

void read_lcpts(struct gargs *genome_arguments, struct pargs *program_arguments) {

    struct tpool *tm;

    tm = tpool_create(program_arguments->thread_number);

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        tpool_add_work(tm, read_lcpt, genome_arguments+i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void read_lcpt(void *arg) {

    struct gargs *genome_arguments = (struct gargs *)arg;

    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rload);
        log1(INFO, "Thread ID: %ld started processing %s", pthread_self(), genome_arguments->inFileName);
        pthread_mutex_unlock(&console_mutex_rload);
    }

    // open fasta file
    FILE *in = fopen(genome_arguments->inFileName, "r");

    if (in == NULL) {
        log1(ERROR, "Error opening file %s", genome_arguments->inFileName);
        return;
    }

    // create file for writing cores
    FILE *out = NULL;

    if (genome_arguments->write_lcpt) {
        out = fopen(genome_arguments->inFileName, "wb");
        if (out == NULL) {
            log1(ERROR, "Error opening file for saving into file %s", genome_arguments->outFileName);
            return;
        }
    }

    // in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));

    // while(!isDone) {
    //     lcp::lps *str = new lcp::lps(in);
    //     str->deepen(program_arguments.lcpLevel);
        
    //     append(str, thread_arguments.cores);

    //     delete str;

    //     in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));
    // } 
    
    fclose(in);

    // end writing cores to file if user specified to do so
    if (genome_arguments->write_lcpt) {
        done(out);
        fclose(out);
    }

    // sort and filter the cores
    genSign(genome_arguments, genome_arguments->sct);

    // log ending of processing fasta
    if (genome_arguments->verbose) {
        log1(INFO, "Thread ID: %s ended processing %s, size: %ld", pthread_self(), genome_arguments->inFileName, genome_arguments->cores_len);
    }
}