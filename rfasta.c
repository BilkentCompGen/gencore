#include "rfasta.h"

pthread_mutex_t console_mutex_rfasta;

void read_fastas(struct gargs *genome_arguments, struct pargs *program_arguments) {
    
    struct tpool *tm;

    tm = tpool_create(program_arguments->thread_number);

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        tpool_add_work(tm, read_fasta, genome_arguments+i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void read_fasta(void *arg) {

    time_t start, point1, point2;
    time(&start);

    struct gargs *genome_arguments = (struct gargs *)arg;

    // open fasta file
    FILE *in = fopen(genome_arguments->inFileName, "r");
    if (in == NULL) {
        log1(ERROR, "Error opening file %s", genome_arguments->inFileName);
        return;
    }

    fseek(in, 0, SEEK_END); // seek to end of file
    uint64_t size = ftell(in); // get current file pointer
    fseek(in, 0, SEEK_SET); // seek back to beginning of fill

    uint64_t estimated_core_size = (int)(size / pow(MAGIC_LCP_FA_CONSTANT, genome_arguments->lcp_level));
    
    genome_arguments->cores = (simple_core*)malloc(estimated_core_size * sizeof(simple_core));
    if (genome_arguments->cores == NULL) {
        pthread_mutex_lock(&console_mutex_rfasta);
        log1(INFO, "Thread ID: %ld couldn't allocate memory of size %ld for cores", pthread_self(), size);
        pthread_mutex_unlock(&console_mutex_rfasta);
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

    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfasta);
        log1(INFO, "Thread - in: %s, cc: %ld", genome_arguments->inFileName, estimated_core_size);
        pthread_mutex_unlock(&console_mutex_rfasta);
    }

    // read file
    char *sequence = (char*)malloc(INITIAL_SEQUENCE_SIZE);
    if (!sequence) {
        pthread_mutex_lock(&console_mutex_rfasta);
        log1(ERROR, "Memory allocation failed for sequence buffer\n");
        exit(EXIT_FAILURE);
    }
    uint64_t sequence_size = 0;
    uint64_t sequence_capacity = INITIAL_SEQUENCE_SIZE;

    char line[1024];

    while (fgets(line, sizeof(line), in)) {
        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            if (sequence_size != 0) {
                process_chrom(sequence, sequence_size, &estimated_core_size, genome_arguments, out);
                sequence_size = 0;
            }
        } else {
            size_t line_len = strlen(line);

            if (sequence_size + line_len >= sequence_capacity) {
                sequence_capacity = (size_t)(sequence_capacity * 1.5);
                sequence = realloc(sequence, sequence_capacity);
                if (!sequence) {
                    log1(ERROR, "Memory reallocation failed\n");
                    exit(EXIT_FAILURE);
                }
            }

            memcpy(sequence + sequence_size, line, line_len);
            sequence_size += line_len;
        }
    }

    if (sequence_size != 0) {
        process_chrom(sequence, sequence_size, &estimated_core_size, genome_arguments, out);
    }

    free(sequence);
    fclose(in);

    // end writing cores to file if user specified to do so
    if (genome_arguments->write_lcpt) {
        done(out);
        fclose(out);
    }

    time(&point1);
    double diff1 = difftime(point1, start);

    // sort and filter the cores
    genSign(genome_arguments, genome_arguments->apply_filter);

    time(&point2);
    double diff2 = difftime(point2, point1);

    // log ending of processing fasta
    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfasta);
        log1(INFO, "Thread - %s, LCP [%d:%d:%d], gen-sign [%d:%d:%d]", genome_arguments->inFileName, (int)(diff1/3600), (int)(diff1/60), (int)(diff1)%60, (int)(diff2/3600), (int)(diff2/60), (int)(diff2)%60);
        pthread_mutex_unlock(&console_mutex_rfasta);
    }
}

void process_chrom(char *sequence, size_t seq_size, uint64_t *capacity, struct gargs *genome_arguments, FILE *out) {

    int valid_chars[256] = {0};
    valid_chars['A'] = valid_chars['C'] = valid_chars['T'] = valid_chars['G'] = 1;
    valid_chars['a'] = valid_chars['c'] = valid_chars['t'] = valid_chars['g'] = 1;
    uint64_t index = 0;

    while (index < seq_size) {
        
        while (index < seq_size && !valid_chars[(unsigned char)sequence[index]]) index++;

        if (index == seq_size) break;

        uint64_t end = index;
        
        while (end < seq_size && valid_chars[(unsigned char)sequence[end]]) end++;

        struct lps str;
        init_lps_offset(&str, sequence+index, end-index, index);
        lps_deepen(&str, genome_arguments->lcp_level);

        if (str.size) {
            
            if (genome_arguments->write_lcpt) save(out, &str);
        
            uint64_t len = genome_arguments->cores_len;

            if (*capacity <= len+str.size) {
                *capacity = *capacity * 1.5;
                simple_core* temp = (simple_core*)realloc(genome_arguments->cores, *capacity);
                if (temp == NULL) {
                    log1(ERROR, "Couldn't increase cores array size.");
                    return;
                }
                genome_arguments->cores = temp;
            }
        
            simple_core *cores = genome_arguments->cores;
        
            for (int i=0; i<str.size; i++) {
                cores[len] = ((uint64_t)str.cores[i].label << 32) + (str.cores[i].end-str.cores[i].start);
                len++;
            }
        
            genome_arguments->cores_len = len;        
        }

        index = end;
        free_lps(&str); 
    }
}