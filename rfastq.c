#include "rfastq.h"

pthread_mutex_t console_mutex_rfastq;

KSEQ_INIT(gzFile, gzread)

void read_fastqs(g_args_t *genome_args, p_args_t *program_args) {
    
    struct tpool *tm;

    tm = tpool_create(program_args->n_threads);

    for (int i=0; i<program_args->n_genomes; i++) {
        struct stat s;
        if (lstat(genome_args->inFileName, &s) == 0) {
            if (S_ISDIR(s.st_mode)) { // directory
                tpool_add_work(tm, process_dir_fastq, genome_args+i);
            } else if (S_ISREG(s.st_mode)) { // file
                tpool_add_work(tm, read_fastq, genome_args+i);
            } else if (S_ISLNK(s.st_mode)) {
                // symbolic link
            } else {
                // something else
            }
        } else {
            //error
            log1(ERROR, "Couldn't process %s.", genome_args->inFileName);
        }
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void process_dir_fastq(void *arg) {

    g_args_t *genome_args = (g_args_t *)arg;

    DIR *dir;
    struct dirent *entry;
    dir = opendir(genome_args->inFileName);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", genome_args->inFileName);
        return;
    }

    char *dir_name = genome_args->inFileName;

    int file_count = 0;
    struct stat file_info;
    char full_path[1024];
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 && ends_with_fq(entry->d_name)) {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);
            if (stat(full_path, &file_info) == 0) {
                if (S_ISREG(file_info.st_mode)) {
                    file_count++;
                }
            }
        }
    }

    closedir(dir);

    dir = opendir(dir_name);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", dir_name);
        return;
    }

    int temp_filter = genome_args->apply_filter;
    sim_calculation_type temp_mode = genome_args->sct;

    simple_core **cores = (simple_core **)malloc(sizeof(simple_core *) * file_count);
    memset(cores, 0, file_count * sizeof(simple_core *));
    uint64_t *sizes = (uint64_t *)malloc(sizeof(uint64_t) * file_count);
    int index = 0;
    genome_args->apply_filter = 0;
    genome_args->sct = VECTOR;
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 && ends_with_fq(entry->d_name)) {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);
            if (stat(full_path, &file_info) == 0) {
                if (S_ISREG(file_info.st_mode)) {
                    genome_args->inFileName = full_path;
                    read_fastq((void*)genome_args);
                    cores[index] = genome_args->cores;
                    sizes[index++] = genome_args->core_count;
                }
            }
        }
    }
    closedir(dir);

    double total_len = 0;
    uint64_t core_count = 0;
    simple_core *new_cores = merge_sorted_arrays(cores, sizes, file_count, genome_args->min_cc, genome_args->max_cc, &core_count, &total_len);
    genome_args->inFileName = dir_name;
    genome_args->cores = new_cores;
    genome_args->core_count = core_count;
    genome_args->total_len = total_len;
    genome_args->apply_filter = temp_filter;
    genome_args->sct = temp_mode;
}

void read_fastq(void *arg) {

    time_t start, point1, point2;
    time(&start);

    g_args_t *genome_args = (g_args_t *)arg;
    uint64_t estimated_core_size = est_core_fq(genome_args->inFileName, genome_args->lcp_level);
    if (!estimated_core_size) {
        log3(ERROR, &console_mutex_rfastq, "Couldn't calculate core size for %s", genome_args->inFileName);
        return;
    }
    
    genome_args->cores = (simple_core *)malloc(sizeof(simple_core) * estimated_core_size);
    if (genome_args->cores == NULL) {
        log3(ERROR, &console_mutex_rfastq, "Thread coultdn't allocate memory - in: %s", estimated_core_size);
        return;
    }
    genome_args->core_count = 0;
    genome_args->total_len = 0;

    gzFile in = gzopen(genome_args->inFileName, "r");
    if (in == NULL) {
        log3(ERROR, &console_mutex_rfastq, "Error opening file %s", genome_args->inFileName);
        return;
    }

    // kseq_t *seq = kseq_init(in);

    // while (kseq_read(seq) >= 0) {
    //     process_read(seq->seq.s, seq->seq.l, &estimated_core_size, genome_args);
    // }

    char *buffer = malloc(INITIAL_SEQUENCE_SIZE);
    if (!buffer) { 
        log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
        return; 
    }

    size_t buffer_len = 0;
    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        if (buffer_len + seq->seq.l >= INITIAL_SEQUENCE_SIZE) {
            // batch is full, process and reset
            process_reads(buffer, buffer_len, &estimated_core_size, genome_args);
            buffer_len = 0;
        }

        memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
        buffer_len += seq->seq.l;
        // add delimiter if to know read boundaries
        buffer[buffer_len++] = SEPERATOR; 
    }

    if (buffer_len) {
        process_reads(buffer, buffer_len, &estimated_core_size, genome_args);
    }

    kseq_destroy(seq);
    gzclose(in);

    time(&point1);
    double diff1 = difftime(point1, start);

    // sort and filter the cores
    genSign(genome_args);

    time(&point2);
    double diff2 = difftime(point2, point1);

    // log ending of processing fastq
    if (genome_args->verbose) {
        log3(INFO, &console_mutex_rfastq, "Thread - %s, LCP [%d:%d:%d], gen-sign [%d:%d:%d], cc: %lu/%lu", genome_args->inFileName, (int)(diff1/3600), (int)(((int)(diff1)%3600)/60), (int)(diff1)%60, (int)(diff2/3600), (int)(((int)(diff2)%3600)/60), (int)(diff2)%60, genome_args->core_count, estimated_core_size);
    }
}

void process_reads(char *sequence, size_t seq_size, uint64_t *capacity, g_args_t *genome_args) {

    uint64_t cap = *capacity;
    size_t start = 0;

    while (start < seq_size) {
        size_t end = start;

        while (end < seq_size && sequence[end] != SEPERATOR) end++;

        // process forward
        struct lps str_fwd;
        init_lps(&str_fwd, sequence+start, end - start);
        lps_deepen(&str_fwd, genome_args->lcp_level);

        uint64_t core_count = genome_args->core_count;

        if (cap <= core_count + str_fwd.size) {
            cap = cap * 1.5;
            simple_core *temp = (simple_core *)realloc(genome_args->cores, sizeof(simple_core) * cap);
            if (temp == NULL) {
                log3(ERROR, &console_mutex_rfastq, "Couldn't increase cores array size.");
                return;
            }
            genome_args->cores = temp;
        }

        simple_core *cores = genome_args->cores;

        for (int i = 0; i < str_fwd.size; i++) {
            cores[core_count] = ((uint64_t)str_fwd.cores[i].label << 32) + (str_fwd.cores[i].end-str_fwd.cores[i].start);
            core_count++;
        }
        
        free_lps(&str_fwd);

        // process reverse complement
        struct lps str_rev;
        init_lps2(&str_rev, sequence+start, end-start);
        lps_deepen(&str_rev, genome_args->lcp_level);

        if (cap <= core_count + str_rev.size) {
            cap = cap * 1.5;
            simple_core *temp = (simple_core *)realloc(genome_args->cores, sizeof(simple_core) * cap);
            if (temp == NULL) {
                log3(ERROR, &console_mutex_rfastq, "Couldn't increase cores array size.");
                return;
            }
            genome_args->cores = temp;
        }

        cores = genome_args->cores;

        for (int i = 0; i < str_rev.size; i++) {
            cores[core_count] = ((uint64_t)str_rev.cores[i].label << 32) + (str_rev.cores[i].end-str_rev.cores[i].start);
            core_count++;
        }
        
        free_lps(&str_rev);

        genome_args->core_count = core_count;
        *capacity = cap;

        start = end + 1;
    }
}