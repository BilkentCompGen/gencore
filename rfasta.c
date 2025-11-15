#include "rfasta.h"

pthread_mutex_t console_mutex_rfasta;

static int cmp_seq_desc(const void *a, const void *b) {
    const seq_t *sa = *(const seq_t **)a;
    const seq_t *sb = *(const seq_t **)b;
    return sb->length - sa->length;  // descending by length
}

#if NUMA_AVAILABLE
void bind_to_node(int node_id) {
    struct bitmask *bm = numa_allocate_nodemask();
    numa_bitmask_setbit(bm, node_id);
    numa_bind(bm);
    numa_free_nodemask(bm);
}
#endif

fa_thread_t *init_threads(p_args_t *program_args) {

#if NUMA_AVAILABLE
    int num_nodes = numa_num_configured_nodes();
#endif

    fa_thread_t *threads = calloc(program_args->n_threads, sizeof(fa_thread_t));
    for (int i = 0; i < program_args->n_threads; i++) {
        threads[i].capacity = 16;
        threads[i].seqs = calloc(threads[i].capacity, sizeof(seq_t));
        threads[i].seq_count = 0;
        threads[i].total_seq_len = 0;
        threads[i].lcp_level = program_args->lcp_level;
        threads[i].verbose = program_args->verbose;
#if NUMA_AVAILABLE
        threads[i].numa_node_id = i % num_nodes;
#endif
    }
    return threads;
}

void assign_seq2thd(fa_thread_t *threads, int n_threads, seq_t *s) {
    // find thread with smallest total_len
    int min_idx = 0;
    long min_len = threads[0].total_seq_len;
    for (int i = 1; i < n_threads; i++) {
        if (threads[i].total_seq_len < min_len) {
            min_len = threads[i].total_seq_len;
            min_idx = i;
        }
    }

    fa_thread_t *t = &threads[min_idx];
    if (t->seq_count == t->capacity) {
        t->capacity *= 2;
        t->seqs = realloc(t->seqs, sizeof(seq_t) * t->capacity);
    }
    t->seqs[t->seq_count].fasta = strdup(s->fasta);
    t->seqs[t->seq_count].name = strdup(s->name);
    t->seqs[t->seq_count].length = s->length;
    t->seqs[t->seq_count].seq_idx = s->seq_idx;
    t->seqs[t->seq_count].cores = NULL;
    t->seqs[t->seq_count++].core_count = 0;
    t->total_seq_len += s->length;
}

void process_chr(char *sequence, seq_t *seq, fa_thread_t *thread_args) {

    int valid_chars[256] = {0};
    valid_chars['A'] = valid_chars['C'] = valid_chars['T'] = valid_chars['G'] = 1;
    valid_chars['a'] = valid_chars['c'] = valid_chars['t'] = valid_chars['g'] = 1;
    uint64_t index = 0;
    uint64_t seq_length = seq->length;
    
    uint64_t est_core_size = (int)(seq->length / pow(MAGIC_LCP_FA_CONSTANT, thread_args->lcp_level));
    seq->cores = (simple_core *)malloc(sizeof(simple_core) * est_core_size);
    if (!seq->cores) {
        log3(ERROR, &console_mutex_rfasta, "Thread %ld, couldn't allocate cores array size.", pthread_self());
        return;
    }
    
    while (index < seq_length) {
        
        while (index < seq_length && !valid_chars[(unsigned char)sequence[index]]) index++;

        if (index == seq_length) break;

        uint64_t end = index;
        
        while (end < seq_length && valid_chars[(unsigned char)sequence[end]]) end++;

        struct lps str;
        init_lps_offset(&str, sequence+index, end-index, index);
        lps_deepen(&str, thread_args->lcp_level);

        if (str.size) {
                    
            uint64_t len = seq->core_count;

            if (est_core_size <= len + str.size) {
                est_core_size = est_core_size * 1.5;
                simple_core *temp = (simple_core *)realloc(seq->cores, sizeof(simple_core) * est_core_size);
                if (temp == NULL) {
                    log3(ERROR, &console_mutex_rfasta, "Thread %ld, couldn't increase cores array size.", pthread_self());
                    return;
                }
                seq->cores = temp;
            }
        
            simple_core *cores = seq->cores;
        
            for (int i = 0; i < str.size; i++) {
                cores[len] = ((uint64_t)str.cores[i].label << 32) + (str.cores[i].end-str.cores[i].start);
                len++;
            }
        
            seq->core_count = len;        
        }

        index = end;
        free_lps(&str); 
    }
}

void thread_process_seqs(void *arg) {

    fa_thread_t *thread_args = (fa_thread_t *)arg;

#if NUMA_AVAILABLE
    bind_to_node(thread_args->numa_node_id);
#endif

    for (int i = 0; i < thread_args->seq_count; i++) {
        
        time_t seq_time_start, seq_time_end;
        time(&seq_time_start);

        faidx_t *fai = fai_load(thread_args->seqs[i].fasta);

        if (!fai) {
            log3(ERROR, &console_mutex_rfasta, "Thread %ld, error loading FASTA index for %s", pthread_self(), thread_args->seqs[i].fasta);
            return;
        }

        int seq_len;
        char *seq = fai_fetch(fai, thread_args->seqs[i].name, &seq_len);
        if (!seq) {
            log3(ERROR, &console_mutex_rfasta, "Thread %ld, region %s not found in %s", pthread_self(), thread_args->seqs[i].name, thread_args->seqs[i].fasta);
            fai_destroy(fai);
            return;
        }

        process_chr(seq, thread_args->seqs + i, thread_args);

        free(seq);
        fai_destroy(fai);

        time(&seq_time_end);
        thread_args->seqs[i].exec_time = difftime(seq_time_end, seq_time_start);
    }
}

void read_fastas(g_args_t *genome_args, p_args_t *program_args) {
    
    // initialize threads arguments
    fa_thread_t *threads = init_threads(program_args);

    int n_genomes = program_args->n_genomes;
    int n_threads = program_args->n_threads;

    for (int i = 0; i < n_genomes; i++) {

        faidx_t *fai = fai_load(genome_args[i].inFileName);
        if (!fai) {
            log3(ERROR, &console_mutex_rfasta, "Failed to read FAI for %s", genome_args[i].inFileName);
            return;
        }

        int nseqs = faidx_nseq(fai);
        seq_t **seqs = malloc(sizeof(seq_t *) * nseqs);
        for (int j = 0; j < nseqs; j++) {
            const char *name = faidx_iseq(fai, j);
            int len = faidx_seq_len(fai, name);
            seqs[j] = malloc(sizeof(seq_t));
            seqs[j]->fasta = strdup(genome_args[i].inFileName);
            seqs[j]->name = strdup(name);
            seqs[j]->length = len;
            seqs[j]->seq_idx = i;
            seqs[j]->exec_time = 0;
        }
        fai_destroy(fai);

        // sort by descending sequence length
        qsort(seqs, nseqs, sizeof(seq_t *), cmp_seq_desc);

        // greedy allocation
        for (int j = 0; j < nseqs; j++) {
            assign_seq2thd(threads, n_threads, seqs[j]);
        }

        // clean up
        for (int j = 0; j < nseqs; j++) {
            free(seqs[j]->fasta);
            free(seqs[j]->name);
            free(seqs[j]);
        }
        free(seqs);
    }

    // create pool and assign jobs
    struct tpool *tm;

    tm = tpool_create(n_threads);

    for (int i = 0; i < n_threads; i++) {
        tpool_add_work(tm, thread_process_seqs, threads+i);
    }

    tpool_wait(tm);
    tpool_destroy(tm);

    for (int i = 0; i < n_genomes; i++) {
        for (int j = 0; j < n_threads; j++) {
            for (int k = 0; k < threads[j].seq_count; k++) {
                if (strcmp(threads[j].seqs[k].fasta, genome_args[i].inFileName) == 0) {
                    genome_args[i].core_count += threads[j].seqs[k].core_count;
                    genome_args[i].time_stats.lcp += threads[j].seqs[k].exec_time;
                }
            }
        }
    }

    for (int i = 0; i < n_genomes; i++) {

        time_t start, end;
        time(&start);

        genome_args[i].cores = (simple_core *)malloc(sizeof(simple_core) * genome_args[i].core_count);
        if (!genome_args[i].cores) {
            log1(ERROR, "Couldn't allocate array.");
            return;
        }

        int idx = 0;
        for (int j = 0; j < n_threads; j++) {
            for (int k = 0; k < threads[j].seq_count; k++) {
                if (strcmp(threads[j].seqs[k].fasta, genome_args[i].inFileName) == 0 && threads[j].seqs[k].core_count) {
                    memcpy(genome_args[i].cores + idx, threads[j].seqs[k].cores, sizeof(simple_core) * threads[j].seqs[k].core_count);
                    idx += threads[j].seqs[k].core_count;

                    free(threads[j].seqs[k].cores);
                    threads[j].seqs[k].core_count = 0;
                    threads[j].seqs[k].cores = NULL;
                }
            }
        }

        time(&end);
        genome_args[i].time_stats.merging = difftime(end, start);
    }

    // generate signitures in parallel
    int n_threads_gen_sign = (n_threads < n_genomes ? n_threads : n_genomes);
    tm = tpool_create(n_threads_gen_sign);

    for (int i = 0; i < n_genomes; i++) {
        tpool_add_work(tm, genSign, genome_args + i);
    }

    tpool_wait(tm);
    tpool_destroy(tm);

    // print summary and stats
    if (genome_args->verbose) {
        for (int i = 0; i < n_genomes; i++) {
            log1(INFO, "Name: %s, cc: %ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args[i].shortName,
                genome_args[i].core_count,
                (int)genome_args[i].time_stats.lcp/3600, ((int)genome_args[i].time_stats.lcp%3600)/60, (int)genome_args[i].time_stats.lcp%60,
                (int)genome_args[i].time_stats.merging/3600, ((int)genome_args[i].time_stats.merging%3600)/60, (int)genome_args[i].time_stats.merging%60,
                (int)genome_args[i].time_stats.sorting/3600, ((int)genome_args[i].time_stats.sorting%3600)/60, (int)genome_args[i].time_stats.sorting%60,
                (int)genome_args[i].time_stats.filtering/3600, ((int)genome_args[i].time_stats.filtering%3600)/60, (int)genome_args[i].time_stats.filtering%60
            );
        }
    }

    // clean up
    for (int j = 0; j < program_args->n_threads; j++) {
        for (int k = 0; k < threads[j].seq_count; k++) {
            free(threads[j].seqs[k].fasta);
            free(threads[j].seqs[k].name);
        }
        free(threads[j].seqs);
        threads[j].seq_count = 0;
    }
    free(threads);
}

void read_fastas_trivial(g_args_t *genome_args, p_args_t *program_args) {
    
    struct tpool *tm;

    tm = tpool_create(program_args->n_threads);

    for (int i = 0; i < program_args->n_genomes; i++) {
        tpool_add_work(tm, read_fasta, genome_args + i);
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void read_fasta(void *arg) {

    time_t start, end;
    time(&start);

    g_args_t *genome_args = (g_args_t *)arg;

    // open fasta file
    FILE *in = fopen(genome_args->inFileName, "r");
    if (in == NULL) {
        log3(ERROR, &console_mutex_rfasta, "Thread %ld, error opening file %s", pthread_self(), genome_args->inFileName);
        return;
    }

    fseek(in, 0, SEEK_END); // seek to end of file
    uint64_t size = ftell(in); // get current file pointer
    fseek(in, 0, SEEK_SET); // seek back to beginning of fill

    uint64_t estimated_core_size = (int)(size / pow(MAGIC_LCP_FA_CONSTANT, genome_args->lcp_level));
    
    genome_args->cores = (simple_core *)malloc(sizeof(simple_core) * estimated_core_size);
    if (genome_args->cores == NULL) {
        log3(ERROR, &console_mutex_rfasta, "Thread %ld couldn't allocate memory of size %ld for cores", pthread_self(), size);
        return;
    }

    // create file for writing cores
    FILE *out = NULL;
    if (genome_args->write_lcpt) {
        out = fopen(genome_args->inFileName, "wb");
        if (out == NULL) {
            log3(ERROR, &console_mutex_rfasta, "Thread %ld, Error opening file for saving into file %s", pthread_self(), genome_args->outFileName);
            return;
        }
    }

    // read file
    char *sequence = (char *)malloc(INITIAL_SEQUENCE_SIZE);
    if (!sequence) {
        log3(ERROR, &console_mutex_rfasta, "Thread %ld, memory allocation failed for sequence buffer.", pthread_self());
        return;
    }
    uint64_t sequence_size = 0;
    uint64_t sequence_capacity = INITIAL_SEQUENCE_SIZE;

    char line[1024];

    while (fgets(line, sizeof(line), in)) {
        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            if (sequence_size != 0) {
                process_chrom(sequence, sequence_size, &estimated_core_size, genome_args, out);
                sequence_size = 0;
            }
        } else {
            uint64_t line_len = strlen(line);

            if (sequence_size + line_len >= sequence_capacity) {
                sequence_capacity = (uint64_t)(sequence_capacity * 1.5);
                sequence = realloc(sequence, sequence_capacity);
                if (!sequence) {
                    log3(ERROR, &console_mutex_rfasta, "Thread %ld, memory reallocation failed.", pthread_self());
                    return;
                }
            }

            memcpy(sequence + sequence_size, line, line_len);
            sequence_size += line_len;
        }
    }

    if (sequence_size != 0) {
        process_chrom(sequence, sequence_size, &estimated_core_size, genome_args, out);
    }

    free(sequence);
    fclose(in);

    // end writing cores to file if user specified to do so
    if (genome_args->write_lcpt) {
        done(out);
        fclose(out);
    }

    time(&end);
    genome_args->time_stats.lcp = difftime(end, start);

    // sort and filter the cores
    genSign(genome_args);

    // log ending of processing fasta
    if (genome_args->verbose) {
        log3(INFO, &console_mutex_rfasta, "Thread - name: %s, cc: %ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
            genome_args->shortName,
            genome_args->core_count,
            (int)genome_args->time_stats.lcp/3600, ((int)genome_args->time_stats.lcp%3600)/60, (int)genome_args->time_stats.lcp%60,
            (int)genome_args->time_stats.merging/3600, ((int)genome_args->time_stats.merging%3600)/60, (int)genome_args->time_stats.merging%60,
            (int)genome_args->time_stats.sorting/3600, ((int)genome_args->time_stats.sorting%3600)/60, (int)genome_args->time_stats.sorting%60,
            (int)genome_args->time_stats.filtering/3600, ((int)genome_args->time_stats.filtering%3600)/60, (int)genome_args->time_stats.filtering%60
        );
    }
}

void process_chrom(char *sequence, uint64_t seq_size, uint64_t *capacity, g_args_t *genome_args, FILE *out) {

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
        lps_deepen(&str, genome_args->lcp_level);

        if (str.size) {
            
            if (genome_args->write_lcpt) save(out, &str);
        
            uint64_t core_count = genome_args->core_count;

            if (*capacity <= core_count + str.size) {
                *capacity = *capacity * 1.5;
                simple_core *temp = (simple_core *)realloc(genome_args->cores, sizeof(simple_core) * (*capacity));
                if (temp == NULL) {
                    log3(ERROR, &console_mutex_rfasta, "Couldn't increase cores array size.");
                    return;
                }
                genome_args->cores = temp;
            }
        
            simple_core *cores = genome_args->cores;
        
            for (int i = 0; i < str.size; i++) {
                cores[core_count] = ((uint64_t)str.cores[i].label << 32) + (str.cores[i].end-str.cores[i].start);
                core_count++;
            }
        
            genome_args->core_count = core_count;        
        }

        index = end;
        free_lps(&str); 
    }
}