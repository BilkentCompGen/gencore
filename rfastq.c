#include "rfastq.h"

pthread_mutex_t console_mutex_rfastq;

KSEQ_INIT(gzFile, gzread)

void read_fastqs(struct gargs *genome_arguments, struct pargs *program_arguments) {
    
    struct tpool *tm;

    tm = tpool_create(program_arguments->thread_number);

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        tpool_add_work(tm, read_fastq, genome_arguments+i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void read_fastq(void *arg) {

    struct gargs *genome_arguments = (struct gargs *)arg;

    gzFile in = gzopen(genome_arguments->inFileName, "r");
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

    // estimate core counts and make array allocation
    struct stat st;
    if (stat(genome_arguments->inFileName, &st) != 0) {
        log1(ERROR, "Error getting file size of %s", genome_arguments->inFileName);
        return;
    }
    uint64_t file_size = st.st_size;

    uint64_t estimated_uncompressed_size = file_size; // default assumption
    if (strstr(genome_arguments->inFileName, ".gz")) {
        // estimate uncompressed size using typical compression ratio (~4:1 for FASTQ)
        estimated_uncompressed_size = file_size * 4;
    }

    uint64_t estimated_bp_count = estimated_uncompressed_size / 2;
    uint64_t estimated_core_size = (int)(estimated_bp_count / pow(MAGIC_LCP_FQ_CONSTANT, genome_arguments->lcp_level));
    genome_arguments->cores = (simple_core*)malloc(estimated_core_size * sizeof(simple_core));

    if (genome_arguments->cores == NULL) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread ID: %ld couldn't allocate memory of size %ld for cores", pthread_self(), estimated_core_size);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }

    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread ID: %ld, in: %s, cc: %ld", pthread_self(), genome_arguments->inFileName, estimated_core_size);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }

    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        process_read(seq->seq.s, seq->seq.l, &estimated_core_size, genome_arguments, out);
    }

    kseq_destroy(seq);
    gzclose(in);

    // end writing cores to file if user specified to do so
    if (genome_arguments->write_lcpt) {
        done(out);
        fclose(out);
    }

    // log ending of reading fasta
    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread ID: %ld ended reading %s, size: %ld", pthread_self(), genome_arguments->inFileName, genome_arguments->cores_len);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }

    // sort and filter the cores
    genSign(genome_arguments, genome_arguments->sct);

    // log ending of processing fasta
    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread ID: %ld ended processing %s, size: %ld", pthread_self(), genome_arguments->inFileName, genome_arguments->cores_len);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }
}

void process_read(char *sequence, size_t seq_size, uint64_t *capacity, struct gargs *genome_arguments, FILE *out) {

    uint64_t cap = *capacity;

    // process forward
    struct lps str_fwd;
    init_lps(&str_fwd, sequence, seq_size);
    lps_deepen(&str_fwd, genome_arguments->lcp_level);

    if (genome_arguments->write_lcpt) {
        save(out, &str_fwd);
    }

    uint64_t len = genome_arguments->cores_len;

    if (cap <= len+str_fwd.size) {
        cap = cap * 1.5;
        simple_core* temp = (simple_core*)realloc(genome_arguments->cores, cap);
        if (temp == NULL) {
            log1(ERROR, "Couldn't increase cores array size.");
            return;
        }
        genome_arguments->cores = temp;
    }

    simple_core *cores = genome_arguments->cores;

    for (int i=0; i<str_fwd.size; i++) {
        cores[len] = ((uint64_t)str_fwd.cores[i].label << 32) + (str_fwd.cores[i].end-str_fwd.cores[i].start);
        len++;
    }
    
    free_lps(&str_fwd);

    // process reverse complement
    struct lps str_rev;
    init_lps2(&str_rev, sequence, seq_size);
    lps_deepen(&str_rev, genome_arguments->lcp_level);

    if (genome_arguments->write_lcpt) {
        save(out, &str_rev);
    }

    if (*capacity <= len+str_rev.size) {
        *capacity = *capacity * 1.5;
        simple_core* temp = (simple_core*)realloc(genome_arguments->cores, *capacity);
        if (temp == NULL) {
            log1(ERROR, "Couldn't increase cores array size.");
            return;
        }
        genome_arguments->cores = temp;
    }

    cores = genome_arguments->cores;

    for (int i=0; i<str_rev.size; i++) {
        cores[len] = ((uint64_t)str_rev.cores[i].label << 32) + (str_rev.cores[i].end-str_rev.cores[i].start);
        len++;
    }

    genome_arguments->cores_len = len;
    *capacity = cap;

    free_lps(&str_rev);
}