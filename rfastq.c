#include "rfastq.h"

pthread_mutex_t console_mutex_rfastq;

KSEQ_INIT(gzFile, gzread)

static inline void cpu_relax(void) {
#if defined(__x86_64__)
    __asm__ __volatile__("pause" ::: "memory");
#elif defined(__aarch64__)
    __asm__ __volatile__("yield" ::: "memory");
#else
    __asm__ __volatile__("" ::: "memory");
#endif
}

#define GET_FILENAME_INDEX(path, skip_index) {                          \
    const char *p = (path);                                             \
    const char *unix_sep = strrchr(p, '/');                             \
    const char *win_sep = strrchr(p, '\\');                             \
    const char *last_sep = (unix_sep > win_sep) ? unix_sep : win_sep;   \
                                                                        \
    /* If a separator was found, return its index + 1 (the start of the filename). */ \
    /* Otherwise, return 0 (the file is the whole string). */           \
    skip_index = (last_sep) ? (int)(last_sep - p + 1) : 0;              \
}

int wait_for_available(atomic_int *flags, int n_threads) {
    while (1) {
        for (int i = 0; i < n_threads; i++) {
            if (atomic_load(&flags[i]) == THREAD_AVAILABLE)
                return i;
        }
        for (int i = 0; i < SHORT_WAIT_TIME; i++) cpu_relax();
    }
}

void read_fastqs(g_args_t *genome_args, p_args_t *program_args) {
    for (int i = 0; i < program_args->n_genomes; i++) {
        struct stat s;
        if (lstat(genome_args->inFileName, &s) == 0) {
            if (S_ISDIR(s.st_mode)) { // directory
                process_dir_fastq(genome_args+i, program_args, FQT_DIR);
            } else if (S_ISREG(s.st_mode)) { // file
                read_fastq(genome_args+i, program_args, FQT_FILE);
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

    // log exec time stats
    if (program_args->verbose) {
        for (int i = 0; i < program_args->n_genomes; i++) {
            log1(INFO, "Main - name: %s, cc: %ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args[i].shortName,
                genome_args[i].core_count,
                (int)genome_args[i].time_stats.lcp/3600, ((int)genome_args[i].time_stats.lcp%3600)/60, (int)genome_args[i].time_stats.lcp%60,
                (int)genome_args[i].time_stats.merging/3600, ((int)genome_args[i].time_stats.merging%3600)/60, (int)genome_args[i].time_stats.merging%60,
                (int)genome_args[i].time_stats.sorting/3600, ((int)genome_args[i].time_stats.sorting%3600)/60, (int)genome_args[i].time_stats.sorting%60,
                (int)genome_args[i].time_stats.filtering/3600, ((int)genome_args[i].time_stats.filtering%3600)/60, (int)genome_args[i].time_stats.filtering%60
            );
        }
    }
}

void process_dir_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type fq_type) {

    DIR *dir;
    struct dirent *entry;
    dir = opendir(genome_args->inFileName);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", genome_args->inFileName);
        return;
    }

    char *dir_name = genome_args->inFileName;

    // count fastq file number under given directory
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

    // initialize lcp cores array (for each file)
    simple_core **cores = (simple_core **)malloc(sizeof(simple_core *) * file_count);
    uint64_t *sizes = (uint64_t *)malloc(sizeof(uint64_t) * file_count);
    
    dir = opendir(dir_name);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", dir_name);
        return;
    }

    if (program_args->verbose) {
        log1(INFO, "Processing %s...", genome_args->shortName);
    }

    int index = 0;
    uint64_t pre_filtering_total_core_count = 0;
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 && ends_with_fq(entry->d_name)) {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);
            if (stat(full_path, &file_info) == 0) {
                if (S_ISREG(file_info.st_mode)) {
                    genome_args->inFileName = full_path;
                    read_fastq(genome_args, program_args, fq_type);
                    cores[index] = genome_args->cores;
                    pre_filtering_total_core_count += genome_args->core_count;
                    sizes[index++] = genome_args->core_count;
                    
                }
            }
        }
    }
    closedir(dir);

    genome_args->inFileName = dir_name;

    // merge lcp cores to single array
    merge_sorted_arrays(cores, sizes, file_count, genome_args);

    // log ending of processing fastq
    if (genome_args->verbose) {
        log1(INFO, "Summary - name: %s, cc: %ld/%ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
            genome_args->shortName,
            genome_args->core_count,
            pre_filtering_total_core_count,
            (int)genome_args->time_stats.lcp/3600, ((int)genome_args->time_stats.lcp%3600)/60, (int)genome_args->time_stats.lcp%60,
            (int)genome_args->time_stats.merging/3600, ((int)genome_args->time_stats.merging%3600)/60, (int)genome_args->time_stats.merging%60,
            (int)genome_args->time_stats.sorting/3600, ((int)genome_args->time_stats.sorting%3600)/60, (int)genome_args->time_stats.sorting%60,
            (int)genome_args->time_stats.filtering/3600, ((int)genome_args->time_stats.filtering%3600)/60, (int)genome_args->time_stats.filtering%60
        );
    }
}

void read_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type fq_type) {

    uint64_t estimated_core_size = est_core_fq(genome_args->inFileName, genome_args->lcp_level);
    if (!estimated_core_size) {
        log3(ERROR, &console_mutex_rfastq, "Couldn't calculate core size for %s", genome_args->inFileName);
        return;
    }

    gzFile in = gzopen(genome_args->inFileName, "r");
    if (in == NULL) {
        log3(ERROR, &console_mutex_rfastq, "Error opening file %s", genome_args->inFileName);
        return;
    }

    // init pool for cummunication and job assignment
    atomic_int *available_buffers = (atomic_int *)malloc(sizeof(atomic_int) * program_args->n_threads);
    if (!available_buffers) {
        log3(ERROR, &console_mutex_rfastq, "Error allocating thread buffer array.");
        return;
    }
    for (int i = 0; i < program_args->n_threads; i++) {
        atomic_init(available_buffers + i, BUFFER_NOT_INITIALIZED);
    }

    struct tpool *tm;
    tm = tpool_create(program_args->n_threads);

    fqw_args_t *args = malloc(sizeof(fqw_args_t) * program_args->n_threads);
    if (!args) {
        log3(ERROR, &console_mutex_rfastq, "Error allocating thread arguments.");
        return;
    }

    // init threads
    for (int i = 0; i < program_args->n_threads; i++) {
        args[i].available = available_buffers + i;
        args[i].buffer_len = -1;
        args[i].lcp_level = program_args->lcp_level;
        args[i].core_count = 0;
        args[i].estimated_core_count = estimated_core_size / program_args->n_threads;
        args[i].time_stats = (time_stats_t){0, 0, 0, 0};
        tpool_add_work(tm, process_reads, args + i);
    }

    char *buffer = (char *)malloc(FASTQ_WORKER_BUFFER_SIZE);
    if (!buffer) { 
        log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
        return; 
    }
    uint64_t buffer_len = 0;
    
    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        if (buffer_len + seq->seq.l >= FASTQ_WORKER_BUFFER_SIZE) {
            // batch is full, assign to available thread and reset
            int idx = wait_for_available(available_buffers, program_args->n_threads);

            args[idx].buffer = buffer;
            args[idx].buffer_len = buffer_len;
            atomic_store(available_buffers + idx, THREAD_BUSY);

            buffer = (char *)malloc(FASTQ_WORKER_BUFFER_SIZE);
            if (!buffer) { 
                log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
                return; 
            }
            buffer_len = 0;
        }

        memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
        buffer_len += seq->seq.l;
        // add delimiter if to know read boundaries
        buffer[buffer_len++] = SEPERATOR; 
    }

    kseq_destroy(seq);
    gzclose(in);

    // assign last batch
    if (buffer_len) {
        int idx = wait_for_available(available_buffers, program_args->n_threads);

        args[idx].buffer = buffer;
        args[idx].buffer_len = buffer_len;
        atomic_store(available_buffers + idx, THREAD_BUSY);
    } else {
        free(buffer);
    }

    for (int i = 0; i < program_args->n_threads; i++) {
        int idx = wait_for_available(available_buffers, program_args->n_threads);

        atomic_store(available_buffers + idx, THREAD_EXIT_SIGNAL);
    }

    tpool_wait(tm);
    tpool_destroy(tm);
    free(available_buffers);

    // merge cores stored in thread arguments
    simple_core *cores;
    time_t start, end;
    time(&start);
    genome_args->core_count = merge_thread_arrays(args, program_args->n_threads, &cores);
    time(&end);
    genome_args->time_stats.merging += difftime(end, start);

    // assign main core array to arguments to be passed back
    genome_args->cores = cores;

    // final merging
    estimated_core_size = 0;
    time_stats_t temp = (time_stats_t){0, 0, 0, 0};
    for (int i = 0; i < program_args->n_threads; i++) {
        // add time stats
        temp.lcp += args[i].time_stats.lcp;
        temp.merging += args[i].time_stats.merging;
        temp.sorting += args[i].time_stats.sorting;
        temp.filtering += args[i].time_stats.filtering;
        // sum total capacities
        estimated_core_size += args[i].estimated_core_count;
    }
        
    genome_args->time_stats.lcp += temp.lcp;
    genome_args->time_stats.merging += temp.merging;
    genome_args->time_stats.sorting += temp.sorting;
    genome_args->time_stats.filtering += temp.filtering;

    // log ending of processing fastq
    if (genome_args->verbose) {
        if (fq_type == FQT_DIR) {
            int skip_index;
            GET_FILENAME_INDEX(genome_args->inFileName, skip_index);
            log1(INFO, "Main - fq: %s, cc: %ld/%ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args->inFileName + skip_index,
                genome_args->core_count,
                estimated_core_size,
                (int)temp.lcp/3600, ((int)temp.lcp%3600)/60, (int)temp.lcp%60,
                (int)temp.merging/3600, ((int)temp.merging%3600)/60, (int)temp.merging%60,
                (int)temp.sorting/3600, ((int)temp.sorting%3600)/60, (int)temp.sorting%60,
                (int)temp.filtering/3600, ((int)temp.filtering%3600)/60, (int)temp.filtering%60
            );
        } else {
            // TODO: genSign (only filter)
            log1(INFO, "Summary - name: %s, cc: %ld/%ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args->shortName,
                genome_args->core_count,
                estimated_core_size,
                (int)genome_args->time_stats.lcp/3600, ((int)genome_args->time_stats.lcp%3600)/60, (int)genome_args->time_stats.lcp%60,
                (int)genome_args->time_stats.merging/3600, ((int)genome_args->time_stats.merging%3600)/60, (int)genome_args->time_stats.merging%60,
                (int)genome_args->time_stats.sorting/3600, ((int)genome_args->time_stats.sorting%3600)/60, (int)genome_args->time_stats.sorting%60,
                (int)genome_args->time_stats.filtering/3600, ((int)genome_args->time_stats.filtering%3600)/60, (int)genome_args->time_stats.filtering%60
            );
        }
    }

    // cleanup
    free(args);
}

void process_reads(void *args) {

    fqw_args_t *fastq_worker_args = (fqw_args_t *)args;

    atomic_int *available = fastq_worker_args->available;
    int lcp_level = fastq_worker_args->lcp_level;
    double lcp_exec_time = 0;
    uint64_t core_count = 0;
    uint64_t core_capacity = fastq_worker_args->estimated_core_count;
    simple_core *cores = (simple_core *)malloc(sizeof(simple_core) * core_capacity);
    if (!cores) {
        log3(ERROR, &console_mutex_rfastq, "Couldn't allocated core array.");
        return;
    }
    
    atomic_store(available, THREAD_AVAILABLE);

    while ( atomic_load(available) != THREAD_EXIT_SIGNAL) {

        if ( atomic_load(available) == THREAD_BUSY ) {

            char *sequence = fastq_worker_args->buffer;
            int buffer_len = fastq_worker_args->buffer_len;

            int start = 0;

            while (start < buffer_len) {
                int end = start;

                while (end < buffer_len && sequence[end] != SEPERATOR) end++;
                
                // process forward
                time_t start_time, end_time;
                time(&start_time);
                struct lps str_fwd;
                init_lps(&str_fwd, sequence + start, end - start);
                lps_deepen(&str_fwd, lcp_level);
                time(&end_time);
                lcp_exec_time += difftime(end_time, start_time);

                if (core_capacity <= core_count + str_fwd.size) {
                    core_capacity = core_capacity * 1.5;
                    simple_core *temp = (simple_core *)realloc(cores, sizeof(simple_core) * core_capacity);
                    if (temp == NULL) {
                        log3(ERROR, &console_mutex_rfastq, "Couldn't increase cores array size.");
                        return;
                    }
                    cores = temp;
                }

                for (int i = 0; i < str_fwd.size; i++) {
                    cores[core_count] = ((uint64_t)str_fwd.cores[i].label << 32) + (str_fwd.cores[i].end - str_fwd.cores[i].start);
                    core_count++;
                }
                
                free_lps(&str_fwd);

                // process reverse complement
                time(&start_time);
                struct lps str_rev;
                init_lps2(&str_rev, sequence + start, end - start);
                lps_deepen(&str_rev, lcp_level);
                time(&end_time);
                lcp_exec_time += difftime(end_time, start_time);

                if (core_capacity <= core_count + str_rev.size) {
                    core_capacity = core_capacity * 1.5;
                    simple_core *temp = (simple_core *)realloc(cores, sizeof(simple_core) * core_capacity);
                    if (temp == NULL) {
                        log3(ERROR, &console_mutex_rfastq, "Couldn't increase cores array size.");
                        return;
                    }
                    cores = temp;
                }

                for (int i = 0; i < str_rev.size; i++) {
                    cores[core_count] = ((uint64_t)str_rev.cores[i].label << 32) + (str_rev.cores[i].end - str_rev.cores[i].start);
                    core_count++;
                }
                
                free_lps(&str_rev);

                start = end + 1;
            }
            
            free(fastq_worker_args->buffer);
            fastq_worker_args->buffer_len = -1;

            atomic_store(available, THREAD_AVAILABLE);
        } else {
            for (int i = 0; i < SHORT_WAIT_TIME; i++) cpu_relax();
        }
    }

    // sort cores for convenience
    time_t start, end;
    time(&start);
    qsort(cores, core_count, sizeof(simple_core), compare_simple_core);
    time(&end);

    fastq_worker_args->cores = cores;
    fastq_worker_args->core_count = core_count;
    fastq_worker_args->estimated_core_count = core_capacity;
    fastq_worker_args->time_stats.lcp = lcp_exec_time;
    fastq_worker_args->time_stats.sorting = difftime(end, start);
}

void read_fastqs_trivial(g_args_t *genome_args, p_args_t *program_args) {
    
    struct tpool *tm;

    tm = tpool_create(program_args->n_threads);

    for (int i = 0; i < program_args->n_genomes; i++) {
        struct stat s;
        if (lstat(genome_args->inFileName, &s) == 0) {
            if (S_ISDIR(s.st_mode)) { // directory
                tpool_add_work(tm, process_dir_fastq_trivial, genome_args+i);
            } else if (S_ISREG(s.st_mode)) { // file
                tpool_add_work(tm, read_fastq_trivial, genome_args+i);
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

void process_dir_fastq_trivial(void *arg) {

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
                    read_fastq_trivial((void*)genome_args);
                    cores[index] = genome_args->cores;
                    sizes[index++] = genome_args->core_count;
                }
            }
        }
    }
    closedir(dir);

    merge_sorted_arrays(cores, sizes, file_count, genome_args);
    genome_args->inFileName = dir_name;
    genome_args->apply_filter = temp_filter;
    genome_args->sct = temp_mode;

    if (genome_args->verbose) {
        log1(INFO, "Summary - name: %s, cc: %ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
            genome_args->shortName,
            genome_args->core_count,
            (int)genome_args->time_stats.lcp/3600, ((int)genome_args->time_stats.lcp%3600)/60, (int)genome_args->time_stats.lcp%60,
            (int)genome_args->time_stats.merging/3600, ((int)genome_args->time_stats.merging%3600)/60, (int)genome_args->time_stats.merging%60,
            (int)genome_args->time_stats.sorting/3600, ((int)genome_args->time_stats.sorting%3600)/60, (int)genome_args->time_stats.sorting%60,
            (int)genome_args->time_stats.filtering/3600, ((int)genome_args->time_stats.filtering%3600)/60, (int)genome_args->time_stats.filtering%60
        );
    }
}

void read_fastq_trivial(void *arg) {

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

    char *buffer = malloc(INITIAL_SEQUENCE_SIZE);
    if (!buffer) { 
        log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
        return; 
    }

    uint64_t buffer_len = 0;
    kseq_t *seq = kseq_init(in);

    while (kseq_read(seq) >= 0) {
        if (buffer_len + seq->seq.l >= INITIAL_SEQUENCE_SIZE) {
            // batch is full, process and reset
            process_reads_trivial(buffer, buffer_len, &estimated_core_size, genome_args);
            buffer_len = 0;
        }

        memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
        buffer_len += seq->seq.l;
        // add delimiter if to know read boundaries
        buffer[buffer_len++] = SEPERATOR; 
    }

    if (buffer_len) {
        process_reads_trivial(buffer, buffer_len, &estimated_core_size, genome_args);
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

void process_reads_trivial(char *sequence, uint64_t seq_size, uint64_t *capacity, g_args_t *genome_args) {

    uint64_t cap = *capacity;
    uint64_t start = 0;

    while (start < seq_size) {
        uint64_t end = start;

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