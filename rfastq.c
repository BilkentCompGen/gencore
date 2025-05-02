#include "rfastq.h"

pthread_mutex_t console_mutex_rfastq;

KSEQ_INIT(gzFile, gzread)

void read_fastqs(struct gargs *genome_arguments, struct pargs *program_arguments) {
    
    struct tpool *tm;

    tm = tpool_create(program_arguments->thread_number);

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        struct stat s;
        if (lstat(genome_arguments->inFileName, &s) == 0) {
            if (S_ISDIR(s.st_mode)) { // directory
                tpool_add_work(tm, process_dir_fastq, genome_arguments+i);
            } else if (S_ISREG(s.st_mode)) { // file
                tpool_add_work(tm, read_fastq, genome_arguments+i);
            } else if (S_ISLNK(s.st_mode)) {
                // symbolic link
            } else {
                // something else
            }
        } else {
            //error
            log1(ERROR, "Couldn't process %s.", genome_arguments->inFileName);
        }
    }

    tpool_wait(tm);

    tpool_destroy(tm);
}

void process_dir_fastq(void *arg) {

    struct gargs *genome_arguments = (struct gargs *)arg;

    DIR *dir;
    struct dirent *entry;
    dir = opendir(genome_arguments->inFileName);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", genome_arguments->inFileName);
        return;
    }

    char *dir_name = genome_arguments->inFileName;

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

    int temp_filter = genome_arguments->apply_filter;
    sim_calculation_type temp_mode = genome_arguments->sct;

    simple_core **cores = (simple_core **)malloc(file_count * sizeof(simple_core *));
    memset(cores, 0, file_count * sizeof(simple_core *));
    uint64_t *sizes = (uint64_t *)malloc(file_count * sizeof(uint64_t));
    int index = 0;
    genome_arguments->apply_filter = 0;
    genome_arguments->sct = VECTOR;
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 && ends_with_fq(entry->d_name)) {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);
            if (stat(full_path, &file_info) == 0) {
                if (S_ISREG(file_info.st_mode)) {
                    genome_arguments->inFileName = full_path;
                    read_fastq((void*)genome_arguments);
                    cores[index] = genome_arguments->cores;
                    sizes[index++] = genome_arguments->cores_len;
                }
            }
        }
    }
    closedir(dir);

    double total_len = 0;
    uint64_t cores_len = 0;
    simple_core *new_cores = merge_sorted_arrays(cores, sizes, file_count, genome_arguments->min_cc, genome_arguments->max_cc, &cores_len, &total_len);
    genome_arguments->inFileName = dir_name;
    genome_arguments->cores = new_cores;
    genome_arguments->cores_len = cores_len;
    genome_arguments->total_len = total_len;
    genome_arguments->apply_filter = temp_filter;
    genome_arguments->sct = temp_mode;
}

void read_fastq(void *arg) {

    time_t start, point1, point2;
    time(&start);

    struct gargs *genome_arguments = (struct gargs *)arg;
    uint64_t estimated_core_size = est_core_fq(genome_arguments->inFileName, genome_arguments->lcp_level);
    if (!estimated_core_size) {
        return;
    }
    
    genome_arguments->cores = (simple_core*)malloc(estimated_core_size * sizeof(simple_core));
    if (genome_arguments->cores == NULL) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread coultdn't allocate memory - in: %s", estimated_core_size);
        pthread_mutex_unlock(&console_mutex_rfastq);
        return;
    }
    genome_arguments->cores_len = 0;
    genome_arguments->total_len = 0;

    gzFile in = gzopen(genome_arguments->inFileName, "r");
    if (in == NULL) {
        log1(ERROR, "Error opening file %s", genome_arguments->inFileName);
        return;
    }

    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread - in: %s, cc: %ld", genome_arguments->inFileName, estimated_core_size);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }

    // kseq_t *seq = kseq_init(in);

    // while (kseq_read(seq) >= 0) {
    //     process_read(seq->seq.s, seq->seq.l, &estimated_core_size, genome_arguments);
    // }

    char *buffer = malloc(INITIAL_SEQUENCE_SIZE);
    if (!buffer) { 
        log1(ERROR, "Malloc failed."); 
        return; 
    }

    size_t buffer_len = 0;
    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        if (buffer_len + seq->seq.l >= INITIAL_SEQUENCE_SIZE) {
            // batch is full, process and reset
            process_reads(buffer, buffer_len, &estimated_core_size, genome_arguments);
            buffer_len = 0;
        }

        memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
        buffer_len += seq->seq.l;
        // add delimiter if to know read boundaries
        buffer[buffer_len++] = SEPERATOR; 
    }

    if (buffer_len) {
        process_reads(buffer, buffer_len, &estimated_core_size, genome_arguments);
    }

    kseq_destroy(seq);
    gzclose(in);

    time(&point1);
    double diff1 = difftime(point1, start);

    // sort and filter the cores
    genSign(genome_arguments, genome_arguments->apply_filter);

    time(&point2);
    double diff2 = difftime(point2, point1);

    // log ending of processing fastq
    if (genome_arguments->verbose) {
        pthread_mutex_lock(&console_mutex_rfastq);
        log1(INFO, "Thread - %s, LCP [%d:%d:%d], gen-sign [%d:%d:%d], cc: %lu", genome_arguments->inFileName, (int)(diff1/3600), (int)(diff1/60), (int)(diff1)%60, (int)(diff2/3600), (int)(diff2/60), (int)(diff2)%60, genome_arguments->cores_len);
        pthread_mutex_unlock(&console_mutex_rfastq);
    }
}

void process_reads(char *sequence, size_t seq_size, uint64_t *capacity, struct gargs *genome_arguments) {

    uint64_t cap = *capacity;
    size_t start = 0;

    while (start < seq_size) {
        size_t end = start;

        while (end < seq_size && sequence[end] != SEPERATOR) end++;

        // process forward
        struct lps str_fwd;
        init_lps(&str_fwd, sequence+start, end - start);
        lps_deepen(&str_fwd, genome_arguments->lcp_level);

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
        init_lps2(&str_rev, sequence+start, end-start);
        lps_deepen(&str_rev, genome_arguments->lcp_level);

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
        
        free_lps(&str_rev);

        genome_arguments->cores_len = len;
        *capacity = cap;

        start = end + 1;
    }
}