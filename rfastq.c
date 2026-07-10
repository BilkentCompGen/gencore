#include "rfastq.h"

pthread_mutex_t console_mutex_rfastq;

KSEQ_INIT(gzFile, gzread)

static inline double fq_seconds_now(void) {
    struct timespec ts;
#if defined(CLOCK_MONOTONIC)
    clock_gettime(CLOCK_MONOTONIC, &ts);
#else
    timespec_get(&ts, TIME_UTC);
#endif
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0;
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

static void fq_reader_done_producing(fq_reader_pool_t *rp) {
    /*
     * atomic_fetch_sub returns the old value.
     * If old value was 1, this was the last active reader.
     */
    if (atomic_fetch_sub(&rp->active_readers, 1) == 1) {
        fq_queue_close(rp->queue);
    }
}

static int fq_try_enter_worker(fq_reader_pool_t *rp) {
    int readers = atomic_load(&rp->active_readers);

    /* Once all readers are done, preserve draining behaviour. */
    if (readers <= 0) {
        atomic_fetch_add(&rp->active_workers, 1);
        return 1;
    }

    int workers = atomic_load(&rp->active_workers);
    int limit = readers * FQ_WORKERS_PER_READER;

    while (workers < limit) {
        if (atomic_compare_exchange_weak(&rp->active_workers, &workers, workers + 1)) {
            return 1;
        }
    }

    return 0;
}

static int fq_worker_add_lps_cores(fq_worker_t *worker, struct lps *str) {
    if (!str || str->size <= 0) return 0;

    if (fq_worker_reserve(worker, (uint64_t)str->size) != 0) return -1;

    for (int i = 0; i < str->size; i++) {
        uint64_t core_len = (uint64_t)(str->cores[i].end - str->cores[i].start);
        worker->cores[worker->count++] = ((uint64_t)str->cores[i].label << 32) + core_len;
    }

    return 0;
}

static int fq_process_one_read(fq_worker_t *worker, char *seq, int len) {
    double time0, time1;

    // process forward strand
    struct lps str_fwd;
    time0 = fq_seconds_now();
    init_lps(&str_fwd, seq, len);
    lps_deepen(&str_fwd, worker->lcp_level);
    time1 = fq_seconds_now();
    worker->time_stats.lcp += time1 - time0;

    if (fq_worker_add_lps_cores(worker, &str_fwd) != 0) {
        free_lps(&str_fwd);
        return -1;
    }
    free_lps(&str_fwd);

    // process reverse strand
    struct lps str_rev;
    time0 = fq_seconds_now();
    init_lps2(&str_rev, seq, len);
    lps_deepen(&str_rev, worker->lcp_level);
    time1 = fq_seconds_now();
    worker->time_stats.lcp += time1 - time0;

    if (fq_worker_add_lps_cores(worker, &str_rev) != 0) {
        free_lps(&str_rev);
        return -1;
    }
    free_lps(&str_rev);

    return 0;
}

static int fq_process_one_batch(fq_worker_t *worker, char *data, size_t len) {
    size_t start = 0;

    while (start < len) {
        size_t end = start;
        while (end < len && data[end] != SEPARATOR) end++;

        if (end > start) {
            size_t read_len = end - start;
            if (read_len > (size_t)INT32_MAX) return -1;

            if (fq_process_one_read(worker, data + start, (int)read_len) != 0) {
                return -1;
            }
        }

        start = end + 1;
    }

    return 0;
}

static void *fq_worker_thread_main(void *ptr) {
    fq_worker_t *worker = (fq_worker_t *)ptr;
    fq_batch_t batch;

    pthread_setname_np(pthread_self(), "fq-worker");

    while (1) {
        double idle0 = fq_seconds_now();
        int ok = fq_queue_pop(worker->queue, &batch);
        double idle1 = fq_seconds_now();

        worker->time_stats.idle += idle1 - idle0;

        if (ok == 0) break;

        double run0 = fq_seconds_now();

        if (fq_process_one_batch(worker, batch.data, batch.len) != 0) {
            // On worker error, close the queue so other threads can stop
            fq_queue_close(worker->queue);
        }

        free(batch.data);

        double run1 = fq_seconds_now();
        worker->time_stats.running += run1 - run0;
    }

    double sort0 = fq_seconds_now();
    if (worker->cores && worker->count > 1) {
        sort_u64_radix(worker->cores, worker->count);
    }
    double sort1 = fq_seconds_now();

    worker->time_stats.sorting += sort1 - sort0;
    return NULL;
}

static void *fq_reader_thread_main(void *ptr) {
    fq_reader_t *reader = (fq_reader_t *)ptr;
    fq_reader_pool_t *rp = reader->pool;

    pthread_setname_np(pthread_self(), "fq-reader");

    while (1) {
        int file_id = atomic_fetch_add(&rp->next_file_id, 1);
        if (file_id >= rp->n_files) {
            fq_reader_done_producing(rp);

            if (!fq_try_enter_worker(rp)) {
                return NULL;
            }
            
            pthread_setname_np(pthread_self(), "fq-rd-worker");
            return fq_worker_thread_main(reader->worker);
        }

        const char *path = rp->files[file_id];
        gzFile in = gzopen(path, "r");
        if (!in) {
            if (rp->verbose) {
                fprintf(stderr, "[fq_parallel] could not open FASTQ: %s\n", path);
            }
            continue;
        }

        kseq_t *seq = kseq_init(in);
        if (!seq) {
            gzclose(in);
            continue;
        }

        char *buffer = (char *)malloc(rp->batch_size);
        if (!buffer) {
            kseq_destroy(seq);
            gzclose(in);
            fq_queue_close(rp->queue);
            return NULL;
        }

        size_t buffer_len = 0;

        while (kseq_read(seq) >= 0) {
            size_t need = seq->seq.l + 1;

            /*
             * Extremely long read: allocate a dedicated oversized batch.
             */
            if (need > rp->batch_size) {
                if (buffer_len > 0) {
                    fq_batch_t b = { buffer, buffer_len, file_id };
                    if (fq_queue_push(rp->queue, b) != 0) {
                        free(buffer);
                        kseq_destroy(seq);
                        gzclose(in);
                        return NULL;
                    }

                    buffer = (char *)malloc(rp->batch_size);
                    if (!buffer) {
                        kseq_destroy(seq);
                        gzclose(in);
                        fq_queue_close(rp->queue);
                        return NULL;
                    }
                    buffer_len = 0;
                }

                char *large = (char *)malloc(need);
                if (!large) {
                    free(buffer);
                    kseq_destroy(seq);
                    gzclose(in);
                    fq_queue_close(rp->queue);
                    return NULL;
                }

                memcpy(large, seq->seq.s, seq->seq.l);
                large[seq->seq.l] = SEPARATOR;

                fq_batch_t b = { large, need, file_id };
                if (fq_queue_push(rp->queue, b) != 0) {
                    free(large);
                    free(buffer);
                    kseq_destroy(seq);
                    gzclose(in);
                    return NULL;
                }

                continue;
            }

            if (buffer_len + need > rp->batch_size) {
                fq_batch_t b = { buffer, buffer_len, file_id };
                if (fq_queue_push(rp->queue, b) != 0) {
                    free(buffer);
                    kseq_destroy(seq);
                    gzclose(in);
                    return NULL;
                }

                buffer = (char *)malloc(rp->batch_size);
                if (!buffer) {
                    kseq_destroy(seq);
                    gzclose(in);
                    fq_queue_close(rp->queue);
                    return NULL;
                }
                buffer_len = 0;
            }

            memcpy(buffer + buffer_len, seq->seq.s, seq->seq.l);
            buffer_len += seq->seq.l;
            buffer[buffer_len++] = SEPARATOR;
        }

        if (buffer_len > 0) {
            fq_batch_t b = { buffer, buffer_len, file_id };
            if (fq_queue_push(rp->queue, b) != 0) {
                free(buffer);
                kseq_destroy(seq);
                gzclose(in);
                return NULL;
            }
        } else {
            free(buffer);
        }

        kseq_destroy(seq);
        gzclose(in);

        if (rp->verbose) {
            log1(INFO, "Processed - name: %s", path);
        }
    }

    return NULL;
}

static uint64_t fq_emit_core_run(simple_core value, uint64_t freq, simple_core *result, uint64_t result_index, int apply_filter, uint32_t min_cc, uint32_t max_cc, int keep_duplicates) {
    int passes = 1;

    if (apply_filter) {
        passes = ((uint64_t)min_cc <= freq && freq <= (uint64_t)max_cc);
    }

    if (!passes) {
        return result_index;
    }

    if (keep_duplicates) {
        for (uint64_t i = 0; i < freq; ++i) {
            result[result_index++] = value;
        }
    } else {
        result[result_index++] = value;
    }

    return result_index;
}

static uint64_t merge_filter_thread_arrays(fq_worker_t *workers, int n_workers, simple_core **cores, g_args_t *genome_args) {
    if (!workers || n_workers <= 0 || !cores || !genome_args) {
        if (cores) {
            *cores = NULL;
        }
        return 0;
    }

    double filter0 = fq_seconds_now();

    uint64_t total_size = 0;

    for (int i = 0; i < n_workers; ++i) {
        total_size += workers[i].count;
    }

    if (total_size == 0) {
        *cores = NULL;
        return 0;
    }

    if (genome_args->write_lcpt) {
        FILE *file = fopen(genome_args->outFileName, "wb");

        if (file != NULL) {
            if (fwrite(&total_size, sizeof(uint64_t), 1, file) == 1) {
                for (int i = 0; i < n_workers; ++i) {
                    if (workers[i].count) {
                        size_t written = fwrite(workers[i].cores, sizeof(simple_core), workers[i].count, file);
                        if (written != workers[i].count) {
                            perror("Failed to write cores");
                            fclose(file);
                        }
                    }
                }
            } else {
                perror("Failed to write total_size");
            }
            fclose(file);
        } else {
            perror("Failed to open file");
        }
    }

    min_heap heap = {
        .data = malloc(sizeof(heap_node) * n_workers),
        .size = 0,
        .capacity = n_workers
    };

    if (!heap.data) {
        *cores = NULL;
        return 0;
    }

    const int keep_duplicates = (genome_args->sct == SIM_CALC_VECTOR);

    simple_core *result;
    if (keep_duplicates) {
        result = (simple_core *)malloc(sizeof(simple_core) * total_size);
    } else {
        result = (simple_core *)malloc(sizeof(simple_core) * total_size / genome_args->min_cc);
    }

    if (!result) {
        free(heap.data);
        *cores = NULL;
        return 0;
    }

    for (int i = 0; i < n_workers; ++i) {
        if (workers[i].count > 0) {
            heap_push(&heap, (heap_node){
                .value = workers[i].cores[0],
                .element_index = 0,
                .array_index = i
            });
        }
    }

    uint64_t result_index = 0;
    uint64_t freq = 0;
    simple_core current = 0;
    int have_current = 0;

    while (heap.size > 0) {
        heap_node min = heap_pop(&heap);
        simple_core value = min.value;

        if (!have_current) {
            current = value;
            freq = 1;
            have_current = 1;
        } else if (value == current) {
            freq++;
        } else {
            result_index = fq_emit_core_run(current, freq, result, result_index, genome_args->apply_filter, genome_args->min_cc, genome_args->max_cc, keep_duplicates);
            current = value;
            freq = 1;
        }

        uint64_t next_idx = min.element_index + 1;

        if (next_idx < workers[min.array_index].count) {
            heap_push(&heap, (heap_node){
                .value = workers[min.array_index].cores[next_idx],
                .element_index = next_idx,
                .array_index = min.array_index
            });
        }
    }

    if (have_current) {
        result_index = fq_emit_core_run(current, freq, result, result_index, genome_args->apply_filter, genome_args->min_cc, genome_args->max_cc, keep_duplicates);
    }

    free(heap.data);

    if (result_index == 0) {
        free(result);
        *cores = NULL;
        return 0;
    }

    simple_core *shrunk = (simple_core *)realloc(result, sizeof(simple_core) * result_index);

    if (shrunk) {
        result = shrunk;
    }

    *cores = result;

    uint64_t total_sequence_len = 0;

    for (uint64_t i = 0; i < result_index; i++) {
        total_sequence_len += result[i] & 0xFFFFFFFF;
    }

    genome_args->result.total_core_len = total_sequence_len;
    genome_args->result.total_sequence_len = total_sequence_len;

    double filter1 = fq_seconds_now();
    genome_args->time_stats.filtering += filter1 - filter0;

    return result_index;
}

static int fq_merge_workers(fq_worker_t *workers, int n_workers, fq_parallel_result_t *out, g_args_t *genome_args) {
    
    if (!workers || n_workers <= 0 || !out || !genome_args) {
        return -1;
    }

    uint64_t raw_count = 0;

    for (int i = 0; i < n_workers; i++) {
        raw_count += workers[i].count;
        out->capacity += workers[i].capacity;
    }

    // merge cores stored in thread arguments
    double merge0 = fq_seconds_now();
    uint64_t final_count = merge_filter_thread_arrays(workers, n_workers, &(out->cores), genome_args);
    double merge1 = fq_seconds_now();

    out->time_stats.merging += merge1 - merge0;

    out->count = final_count;
    out->capacity = final_count;

    for (int i = 0; i < n_workers; i++) {
        out->time_stats.running += workers[i].time_stats.running;
        out->time_stats.idle += workers[i].time_stats.idle;
        out->time_stats.lcp += workers[i].time_stats.lcp;
        out->time_stats.sorting += workers[i].time_stats.sorting;
    }

    return 0;
}

void fq_normalize_options(fq_parallel_options_t *dst, const fq_parallel_options_t *src) {
    memset(dst, 0, sizeof(*dst));
    if (src) *dst = *src;

    if (dst->n_readers <= 0) dst->n_readers = 1;
    if (dst->n_workers <= 0) dst->n_workers = 1;
    if (dst->batch_size == 0) dst->batch_size = DEFAULT_FASTQ_BATCH_SIZE;
    if (dst->queue_capacity <= 0) dst->queue_capacity = DEFAULT_FASTQ_QUEUE_CAPACITY;
    if (dst->lcp_level <= 0) dst->lcp_level = DEFAULT_LCP_LEVEL;
}

int fq_parallel_process_files(const char **files, int n_files, const fq_parallel_options_t *opt_in, fq_parallel_result_t *out, g_args_t *genome_args) {
    if (!files || n_files <= 0 || !out) return -1;

    memset(out, 0, sizeof(fq_parallel_result_t));

    fq_parallel_options_t opt;
    fq_normalize_options(&opt, opt_in);

    if (opt.n_readers > n_files) opt.n_readers = n_files;

    fq_batch_queue_t queue;
    if (fq_queue_init(&queue, opt.queue_capacity) != 0) {
        return -1;
    }

    pthread_t *reader_threads = NULL;
    pthread_t *worker_threads = NULL;
    fq_worker_t *workers = NULL;
    fq_reader_t *readers = NULL;

    int rc = -1;

    reader_threads = (pthread_t *)calloc((size_t)opt.n_readers, sizeof(pthread_t));
    worker_threads = (pthread_t *)calloc((size_t)opt.n_workers, sizeof(pthread_t));
    workers = (fq_worker_t *)calloc((size_t)(opt.n_readers+opt.n_workers), sizeof(fq_worker_t));
    readers = (fq_reader_t *)calloc((size_t)(opt.n_readers), sizeof(fq_reader_t));

    if (!reader_threads || !worker_threads || !workers) goto cleanup;

    uint64_t per_worker_cap = (opt.estimated_total_core_capacity > 0) ? (opt.estimated_total_core_capacity / (uint64_t)opt.n_workers + 1) : FQ_PARALLEL_MIN_CORE_CAP;
    uint64_t per_reader_cap = (per_worker_cap * (n_files - (n_files % opt.n_readers))) / (opt.n_readers + opt.n_workers);

    if (per_worker_cap < FQ_PARALLEL_MIN_CORE_CAP) {
        per_worker_cap = FQ_PARALLEL_MIN_CORE_CAP;
    }
    
    if (per_reader_cap < FQ_PARALLEL_MIN_CORE_CAP) {
        per_reader_cap = FQ_PARALLEL_MIN_CORE_CAP;
    }

    for (int i = 0; i < (opt.n_workers); i++) {
        workers[i].queue = &queue;
        workers[i].worker_id = i;
        workers[i].lcp_level = opt.lcp_level;
        workers[i].capacity = per_worker_cap;
        workers[i].cores = (simple_core *)malloc(sizeof(simple_core) * per_worker_cap);
        if (!workers[i].cores) goto cleanup;
    }

    for (int i = opt.n_workers; i < (opt.n_readers+opt.n_workers); i++) {
        workers[i].queue = &queue;
        workers[i].worker_id = i;
        workers[i].lcp_level = opt.lcp_level;
        workers[i].capacity = per_reader_cap;
        workers[i].cores = (simple_core *)malloc(sizeof(simple_core) * per_reader_cap);
        if (!workers[i].cores) goto cleanup;
    }

    fq_reader_pool_t reader_pool;
    memset(&reader_pool, 0, sizeof(fq_reader_pool_t));
    reader_pool.files = files;
    reader_pool.n_files = n_files;
    reader_pool.queue = &queue;
    reader_pool.batch_size = opt.batch_size;
    reader_pool.verbose = opt.verbose;
    atomic_init(&reader_pool.next_file_id, 0);
    atomic_init(&reader_pool.active_readers, opt.n_readers);
    atomic_init(&reader_pool.active_workers, opt.n_workers);

    for (int i = 0; i < (opt.n_readers); i++) {
        readers[i].reader_id = i;
        readers[i].pool = &reader_pool;
        readers[i].worker = workers + opt.n_workers + i;
    }

    for (int i = 0; i < opt.n_workers; i++) {
        if (pthread_create(&worker_threads[i], NULL, fq_worker_thread_main, &workers[i]) != 0) {
            fq_queue_close(&queue);
            goto cleanup_join;
        }
    }

    for (int i = 0; i < opt.n_readers; i++) {
        if (pthread_create(&reader_threads[i], NULL, fq_reader_thread_main, &readers[i]) != 0) {
            fq_queue_close(&queue);
            goto cleanup_join;
        }
    }

    for (int i = 0; i < opt.n_readers; i++) {
        if (reader_threads[i]) pthread_join(reader_threads[i], NULL);
    }

    fq_queue_close(&queue);

    for (int i = 0; i < opt.n_workers; i++) {
        if (worker_threads[i]) pthread_join(worker_threads[i], NULL);
    }

    if (fq_merge_workers(workers, opt.n_workers + opt.n_readers, out, genome_args) != 0) goto cleanup;

    rc = 0;
    goto cleanup;

cleanup_join:
    for (int i = 0; i < opt.n_readers; i++) {
        if (reader_threads[i]) pthread_join(reader_threads[i], NULL);
    }

    fq_queue_close(&queue);

    for (int i = 0; i < opt.n_workers; i++) {
        if (worker_threads[i]) pthread_join(worker_threads[i], NULL);
    }

cleanup:
    if (workers) {
        for (int i = 0; i < opt.n_workers + opt.n_readers; i++) {
            free(workers[i].cores);
        }
    }

    if (readers) {
        free(readers);
    }

    free(workers);
    free(reader_threads);
    free(worker_threads);
    fq_queue_destroy(&queue);

    if (rc != 0) {
        fq_parallel_result_destroy(out);
    }

    return rc;
}

void read_fastqs(g_args_t *genome_args, p_args_t *program_args) {
    for (int i = 0; i < program_args->n_genomes; i++) {
        struct stat s;
        if (lstat(genome_args->inFileName, &s) == 0) {
            if (S_ISDIR(s.st_mode)) { // directory
                process_dir_fastq(genome_args+i, program_args);
            } else if (S_ISREG(s.st_mode)) { // file
                read_fastq(genome_args+i, program_args);
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

void process_dir_fastq(g_args_t *genome_args, p_args_t *program_args) {

    DIR *dir = opendir(genome_args->inFileName);
    if (dir == NULL) {
        log1(ERROR, "Couldn't open dir %s.", genome_args->inFileName);
        return;
    }

    char *dir_name = genome_args->inFileName;

    // count fastq file number under given directory
    int capacity = 128;
    int file_count = 0;
    char **files = malloc(sizeof(char *) * capacity);

    struct dirent *entry;
    struct stat file_info;
    char full_path[1024];

    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 ||
            strcmp(entry->d_name, "..") == 0 ||
            !ends_with_fq(entry->d_name)) {
            continue;
        }
        
        snprintf(full_path, sizeof(full_path), "%s/%s", dir_name, entry->d_name);

        if (stat(full_path, &file_info) == 0 && S_ISREG(file_info.st_mode)) {
            if (file_count == capacity) {
                capacity *= 2;
                char **tmp = realloc(files, sizeof(char *) * capacity);
                if (!tmp) {
                    log1(ERROR, "Could not grow FASTQ file list.");
                    closedir(dir);
                    return;
                }
                files = tmp;
            }

            files[file_count++] = strdup(full_path);
        }
    }

    closedir(dir);

    if (file_count == 0) {
        log1(WARN, "No FASTQ files found under %s.", dir_name);
        free(files);
        return;
    }

    uint64_t estimated_total_core_size = 0;
    for (int i = 0; i < file_count; i++) {
        estimated_total_core_size += est_core_fq(files[i], genome_args->lcp_level);
    }

    log1(INFO, "Total estimated %lu", estimated_total_core_size);

    fq_parallel_options_t opt = {
        .n_readers = program_args->n_readers,
        .n_workers = program_args->n_threads - program_args->n_readers,
        .batch_size = 4 * 1024 * 1024,
        .queue_capacity = (program_args->n_threads - program_args->n_readers) * 12,
        .lcp_level = program_args->lcp_level,
        .estimated_total_core_capacity = estimated_total_core_size,
        .verbose = program_args->verbose
    };

    fq_parallel_result_t result;

    if (fq_parallel_process_files((const char **)files, file_count, &opt, &result, genome_args) != 0) {
        log1(ERROR, "Parallel directory FASTQ processing failed for %s.", dir_name);

        for (int i = 0; i < file_count; i++) free(files[i]);
        free(files);
        return;
    }

    genome_args->result.cores = result.cores;
    genome_args->result.count = result.count;

    genome_args->time_stats.running = result.time_stats.running;
    genome_args->time_stats.idle = result.time_stats.idle;
    genome_args->time_stats.lcp = result.time_stats.lcp;
    genome_args->time_stats.merging = result.time_stats.merging;
    genome_args->time_stats.sorting = result.time_stats.sorting;

    // // filter based on fequency
    // genSign(genome_args);

    if (program_args->verbose) {
        log1(INFO,
             "Summary - name: %s, files: %d, readers: %d, workers: %d, cc: %ld/%ld, worker idle %.2f%%",
             genome_args->shortName,
             file_count,
             program_args->n_readers,
             (program_args->n_threads - program_args->n_readers),
             genome_args->result.count,
             result.capacity,
             result.time_stats.idle /
                 (result.time_stats.idle + result.time_stats.running + 1e-9) * 100.0);
    }

    for (int i = 0; i < file_count; i++) free(files[i]);
    free(files);
}

void read_fastq(g_args_t *genome_args, p_args_t *program_args) {

    uint64_t estimated_core_size = est_core_fq(genome_args->inFileName, genome_args->lcp_level);
    if (!estimated_core_size) {
        log3(ERROR, &console_mutex_rfastq, "Couldn't calculate core size for %s", genome_args->inFileName);
        return;
    }

    const char *files[1];
    files[0] = genome_args->inFileName;

    fq_parallel_options_t opt = {
        .n_readers = 1,  // keep 1 for one normal .fastq.gz
        .n_workers = program_args->n_threads,
        .batch_size = 4 * 1024 * 1024,
        .queue_capacity = program_args->n_threads * 4,
        .lcp_level = program_args->lcp_level,
        .estimated_total_core_capacity = estimated_core_size,
        .verbose = program_args->verbose
    };

    fq_parallel_result_t result;

    if (fq_parallel_process_files(files, 1, &opt, &result, genome_args) != 0) {
        log3(ERROR, &console_mutex_rfastq,
             "Parallel FASTQ processing failed for %s", genome_args->inFileName);
        return;
    }

    genome_args->result.cores = result.cores;
    genome_args->result.count = result.count;

    genome_args->time_stats.running = result.time_stats.running;
    genome_args->time_stats.idle = result.time_stats.idle;
    genome_args->time_stats.lcp = result.time_stats.lcp;
    genome_args->time_stats.merging = result.time_stats.merging;
    genome_args->time_stats.sorting = result.time_stats.sorting;

    if (genome_args->verbose) {
        log1(INFO,
             "Summary - name: %s, cc: %ld/%ld, LCP %.2f, sort %.2f, worker idle %.2f%%",
             genome_args->shortName,
             genome_args->result.count,
             result.capacity,
             result.time_stats.lcp,
             result.time_stats.sorting,
             result.time_stats.idle /
                 (result.time_stats.idle + result.time_stats.running + 1e-9) * 100.0);
    }
}