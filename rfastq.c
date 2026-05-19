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
                process_dir_fastq(genome_args+i, program_args, FASTQ_INPUT_DIR);
            } else if (S_ISREG(s.st_mode)) { // file
                read_fastq(genome_args+i, program_args, FASTQ_INPUT_FILE);
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
                genome_args[i].result.count,
                (int)genome_args[i].time_stats.lcp/3600, ((int)genome_args[i].time_stats.lcp%3600)/60, (int)genome_args[i].time_stats.lcp%60,
                (int)genome_args[i].time_stats.merging/3600, ((int)genome_args[i].time_stats.merging%3600)/60, (int)genome_args[i].time_stats.merging%60,
                (int)genome_args[i].time_stats.sorting/3600, ((int)genome_args[i].time_stats.sorting%3600)/60, (int)genome_args[i].time_stats.sorting%60,
                (int)genome_args[i].time_stats.filtering/3600, ((int)genome_args[i].time_stats.filtering%3600)/60, (int)genome_args[i].time_stats.filtering%60
            );
        }
    }
}

void process_dir_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type_t fq_type) {

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

    genome_args->result.total_core_len = 0;
    genome_args->result.total_sequence_len = 0;

    int index = 0;
    uint64_t pre_filtering_total_core_count = 0;
    
    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 && ends_with_fq(entry->d_name)) {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);
            if (stat(full_path, &file_info) == 0) {
                if (S_ISREG(file_info.st_mode)) {
                    genome_args->inFileName = full_path;
                    read_fastq(genome_args, program_args, fq_type);
                    cores[index] = genome_args->result.cores;
                    pre_filtering_total_core_count += genome_args->result.count;
                    sizes[index++] = genome_args->result.count;
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
            genome_args->result.count,
            pre_filtering_total_core_count,
            (int)genome_args->time_stats.lcp/3600, ((int)genome_args->time_stats.lcp%3600)/60, (int)genome_args->time_stats.lcp%60,
            (int)genome_args->time_stats.merging/3600, ((int)genome_args->time_stats.merging%3600)/60, (int)genome_args->time_stats.merging%60,
            (int)genome_args->time_stats.sorting/3600, ((int)genome_args->time_stats.sorting%3600)/60, (int)genome_args->time_stats.sorting%60,
            (int)genome_args->time_stats.filtering/3600, ((int)genome_args->time_stats.filtering%3600)/60, (int)genome_args->time_stats.filtering%60
        );
    }
}

void read_fastq(g_args_t *genome_args, p_args_t *program_args, fastq_input_type_t fq_type) {

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
        memset(&(args[i].result), 0, sizeof(core_result_t));
        args[i].result.capacity = estimated_core_size / program_args->n_threads;
        memset(&(args[i].time_stats), 0, sizeof(time_stats_t));
        tpool_add_work(tm, process_reads, args + i);
    }

    char *buffer = (char *)malloc(DEFAULT_FASTQ_BATCH_SIZE);
    if (!buffer) { 
        log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
        return; 
    }
    uint64_t buffer_len = 0;
    
    kseq_t *seq = kseq_init(in);

    double total_idle_time = 0;
    time_t running_start, running_end;
    time(&running_start);

    while (kseq_read(seq) >= 0) {
        if (buffer_len + seq->seq.l >= DEFAULT_FASTQ_BATCH_SIZE) {
            // batch is full, assign to available thread and reset
            time_t idle_start, idle_end;
            time(&idle_start);
            int idx = wait_for_available(available_buffers, program_args->n_threads);
            time(&idle_end);
            total_idle_time += difftime(idle_end, idle_start);

            args[idx].buffer = buffer;
            args[idx].buffer_len = buffer_len;
            atomic_store(available_buffers + idx, THREAD_BUSY);

            buffer = (char *)malloc(DEFAULT_FASTQ_BATCH_SIZE);
            if (!buffer) { 
                log3(ERROR, &console_mutex_rfastq, "Malloc failed."); 
                return; 
            }
            buffer_len = 0;
        }

        memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
        buffer_len += seq->seq.l;
        // add delimiter if to know read boundaries
        buffer[buffer_len++] = SEPARATOR; 
    }
    time(&running_end);

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
    genome_args->result.count = merge_thread_arrays(args, program_args->n_threads, &cores);
    time(&end);
    genome_args->time_stats.merging += difftime(end, start);

    // assign main core array to arguments to be passed back
    genome_args->result.cores = cores;

    // final merging
    uint64_t total_core_cap = 0;
    time_stats_t temp = (time_stats_t){0, 0, 0, 0, 0, 0};
    for (int i = 0; i < program_args->n_threads; i++) {
        // add time stats
        temp.running += args[i].time_stats.running;
        temp.idle += args[i].time_stats.idle;
        temp.lcp += args[i].time_stats.lcp;
        temp.merging += args[i].time_stats.merging;
        temp.sorting += args[i].time_stats.sorting;
        temp.filtering += args[i].time_stats.filtering;
        // sum total capacities
        total_core_cap += args[i].result.capacity;
        genome_args->result.total_core_len += args[i].result.total_core_len;
        genome_args->result.total_sequence_len += args[i].result.total_sequence_len;
    }
    
    genome_args->time_stats.running += temp.running;
    genome_args->time_stats.idle += temp.idle;
    genome_args->time_stats.lcp += temp.lcp;
    genome_args->time_stats.merging += temp.merging;
    genome_args->time_stats.sorting += temp.sorting;
    genome_args->time_stats.filtering += temp.filtering;

    // log ending of processing fastq
    if (genome_args->verbose) {
        if (fq_type == FASTQ_INPUT_DIR) {
            int skip_index;
            GET_FILENAME_INDEX(genome_args->inFileName, skip_index);
            log1(INFO, "Main - fq: %s, cc: %ld/%ld, idle: %.2f/%.2f, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args->inFileName + skip_index,
                genome_args->result.count,
                total_core_cap,
                total_idle_time / difftime(running_end, running_start),
                (temp.idle) / (temp.idle + temp.running) * 100,
                (int)temp.lcp/3600, ((int)temp.lcp%3600)/60, (int)temp.lcp%60,
                (int)temp.merging/3600, ((int)temp.merging%3600)/60, (int)temp.merging%60,
                (int)temp.sorting/3600, ((int)temp.sorting%3600)/60, (int)temp.sorting%60,
                (int)temp.filtering/3600, ((int)temp.filtering%3600)/60, (int)temp.filtering%60
            );
        } else {
            // TODO: genSign (only filter)
            log1(INFO, "Summary - name: %s, cc: %ld/%ld, LCP [%d:%02d:%02d], merge: [%d:%02d:%02d], sort: [%d:%02d:%02d], filter: [%d:%02d:%02d]", 
                genome_args->shortName,
                genome_args->result.count,
                total_core_cap,
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
    double total_idle_time = 0;
    double total_running_time = 0;
    uint64_t core_count = 0;
    uint64_t total_core_len = 0;
    uint64_t total_genome_len = 0;
    uint64_t core_capacity = fastq_worker_args->result.capacity;
    simple_core *cores = (simple_core *)malloc(sizeof(simple_core) * core_capacity);
    if (!cores) {
        log3(ERROR, &console_mutex_rfastq, "Couldn't allocated core array.");
        return;
    }
    
    atomic_store(available, THREAD_AVAILABLE);

    while ( atomic_load(available) != THREAD_EXIT_SIGNAL) {

        if ( atomic_load(available) == THREAD_BUSY ) {

            time_t running_start, running_end;
            time(&running_start);

            char *sequence = fastq_worker_args->buffer;
            int buffer_len = fastq_worker_args->buffer_len;

            int start = 0;

            while (start < buffer_len) {
                int end = start;

                while (end < buffer_len && sequence[end] != SEPARATOR) end++;
                
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

                if (str_fwd.size) {
                    for (int i = 0; i < str_fwd.size; i++) {
                        cores[core_count] = ((uint64_t)str_fwd.cores[i].label << 32) + (str_fwd.cores[i].end - str_fwd.cores[i].start);
                        core_count++;
                        total_core_len += (str_fwd.cores[i].end - str_fwd.cores[i].start);
                    }
                    total_genome_len = str_fwd.cores[str_fwd.size-1].end - str_fwd.cores[0].start;
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

                if (str_rev.size) {
                    for (int i = 0; i < str_rev.size; i++) {
                        cores[core_count] = ((uint64_t)str_rev.cores[i].label << 32) + (str_rev.cores[i].end - str_rev.cores[i].start);
                        core_count++;
                        total_core_len += (str_rev.cores[i].end - str_rev.cores[i].start);
                    }
                    total_genome_len = str_rev.cores[str_rev.size-1].end - str_rev.cores[0].start;
                }
                
                free_lps(&str_rev);

                start = end + 1;
            }
            
            free(fastq_worker_args->buffer);
            fastq_worker_args->buffer_len = -1;

            atomic_store(available, THREAD_AVAILABLE);
            
            time(&running_end);
            total_running_time += difftime(running_end, running_start);
        } else {
            time_t idle_start, idle_end;
            time(&idle_start);
            for (int i = 0; i < SHORT_WAIT_TIME; i++) cpu_relax();
            time(&idle_end);
            total_idle_time += difftime(idle_end, idle_start);
        }
    }

    // sort cores for convenience
    time_t start, end;
    time(&start);
    qsort(cores, core_count, sizeof(simple_core), compare_simple_core);
    time(&end);

    fastq_worker_args->result.cores = cores;
    fastq_worker_args->result.count = core_count;
    fastq_worker_args->result.capacity = core_capacity;
    fastq_worker_args->result.total_core_len = total_core_len;
    fastq_worker_args->result.total_sequence_len = total_genome_len;
    fastq_worker_args->time_stats.idle = total_idle_time;
    fastq_worker_args->time_stats.running = total_running_time;
    fastq_worker_args->time_stats.lcp = lcp_exec_time;
    fastq_worker_args->time_stats.sorting = difftime(end, start);
}