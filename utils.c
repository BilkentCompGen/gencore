#include "utils.h"

void calcUISize(const struct gargs *argument1, const struct gargs *argument2, uint64_t *interSize, uint64_t *unionSize) {
    
    uint64_t is = 0;
    uint64_t us = 0;
    uint64_t size1 = argument1->cores_len;
    uint64_t size2 = argument2->cores_len;
    uint64_t index1 = 0;
    uint64_t index2 = 0;

    const simple_core *cores1 = argument1->cores;
    const simple_core *cores2 = argument2->cores;

    while (index1 < size1 && index2 < size2) {
        us++;
        
        if (cores1[index1] == cores2[index2]) {
            is++;
            index1++;
            index2++;
        } else if (cores1[index1] < cores2[index2]) {
            index1++;
        } else {
            index2++;
        }
    }

    us += (size1-index1);
    us += (size2-index2);

    *interSize = is;
    *unionSize = us;
}

double calcJaccardSim(uint64_t interSize, uint64_t unionSize) {
    return (double)interSize / (double)unionSize;
}

double calcDiceSim(uint64_t interSize, uint64_t size1, uint64_t size2) {
    return 2 * (double)interSize / ((double)size1+(double)size2);
}

double calcHammDist(double jaccardSim, double avgLen) {
    return 1 - pow(jaccardSim, 1.0/avgLen);
}

double calcJukesCantorCor(double hammingDist) {
    return - 3.0/4.0 * log(1 - hammingDist * 4.0/3.0);
}

void calcDistances(const struct gargs *genome_arguments, const struct pargs* program_arguments) {

    log1(INFO, "Calculating distance matrices...");

    // Initialize similarity matrices
    double jaccard[program_arguments->number_of_genomes][program_arguments->number_of_genomes];
    double dice[program_arguments->number_of_genomes][program_arguments->number_of_genomes];
    double jukes_cantor[program_arguments->number_of_genomes][program_arguments->number_of_genomes];
    
    for (int i = 0; i < program_arguments->number_of_genomes; i++) {
        jaccard[i][i] = 0.0;
        dice[i][i] = 0.0;
        jukes_cantor[i][i]= 0.0;
    }

    // Compute similarity scores
    for (int i=0; i<program_arguments->number_of_genomes; i++) { 
        for (int j=i+1; j<program_arguments->number_of_genomes; j++) {
            
            size_t interSize, unionSize;
            calcUISize(&(genome_arguments[i]), &(genome_arguments[j]), &interSize, &unionSize);

            double diceDist = 1.0 - calcDiceSim(interSize, genome_arguments[i].cores_len, genome_arguments[j].cores_len);
            double jaccardDist = 1.0 - calcJaccardSim(interSize, unionSize);
            
            double avg_len = (genome_arguments[i].total_len+genome_arguments[j].total_len)/(genome_arguments[i].cores_len+genome_arguments[j].cores_len);
            double jukesCantorDist = calcHammDist(1.0 - jaccardDist, avg_len);
            jukesCantorDist = calcJukesCantorCor(jukesCantorDist);

            dice[i][j] = diceDist;
            jaccard[i][j] = jaccardDist;
            jukes_cantor[i][j] = jukesCantorDist;

            // set values to transposed locations
            dice[j][i] = diceDist;
            jaccard[j][i] = jaccardDist;
            jukes_cantor[j][i] = jukesCantorDist;
        }
    }

    log1(INFO, "Writing distance matrices to files...");

    int lcp_level = genome_arguments[0].lcp_level;
    
    // Write outputs to files
    char *program_type = genome_arguments[0].sct == SET ? "set" : "vec";
    FILE *dice_out, *jaccard_out, *jukes_cantor_out;
    char filename_buffer[256];
    if (snprintf(filename_buffer, 256, "%s.%s.%s%d.phy", program_arguments->prefix, program_type, "dice.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for dice overflow.");
        exit(EXIT_FAILURE);
    }
    dice_out = fopen(filename_buffer, "w");

    if (snprintf(filename_buffer, 256, "%s.%s.%s%d.phy", program_arguments->prefix, program_type, "jaccard.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for jaccard overflow.");
        exit(EXIT_FAILURE);
    }
    jaccard_out = fopen(filename_buffer, "w");

    if (snprintf(filename_buffer, 256, "%s.%s.%s%d.phy", program_arguments->prefix, program_type, "jc.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for jc overflow.");
        exit(EXIT_FAILURE);
    }
    jukes_cantor_out = fopen(filename_buffer, "w");

    // write dice
    if (dice_out) {
        fprintf(dice_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(dice_out, "%-10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(dice_out, " %-.15f", dice[i][j]);
            }
            fprintf(dice_out, "\n");
        }
        fclose(dice_out);
    }

    // write jaccard
    if (jaccard_out) {
        fprintf(jaccard_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(jaccard_out, "%-10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(jaccard_out, " %-.15f", jaccard[i][j]);
            }
            fprintf(jaccard_out, "\n");
        }
        fclose(jaccard_out);
    }

    // write jukes cantor
    if (jukes_cantor_out) {
        fprintf(jukes_cantor_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(jukes_cantor_out, "%-10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(jukes_cantor_out, " %-.15f", jukes_cantor[i][j]);
            }
            fprintf(jukes_cantor_out, "\n");
        }
        fclose(jukes_cantor_out);
    }
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: LCP cores related functions
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

void quicksort(simple_core *array, int low, int high) {
    if (low < high) {
        simple_core pivot = array[high];
        int i = low - 1;

        for (int j = low; j < high; j++) {
            if (array[j] < pivot) {
                i++;
                simple_core temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        simple_core temp = array[i+1];
        array[i+1] = array[high];
        array[high] = temp;

        quicksort(array, low, i);
        quicksort(array, i + 2, high);
    }
}

void genSign(struct gargs *genome_arguments, int apply_filter) {

    simple_core *cores = genome_arguments->cores;
    uint64_t len = genome_arguments->cores_len;
    double total_len = genome_arguments->total_len;

    quicksort(cores, 0, len);
    
    if (apply_filter) {
        uint32_t min_cc = genome_arguments->min_cc;
        uint32_t max_cc = genome_arguments->max_cc;
        uint64_t index = 0;
        uint64_t i = 0;

        while (i<len) {
            uint64_t freq = 1;

            for (uint64_t j=i+1; j<len && cores[i]==cores[j]; j++, freq++);

            if (min_cc<=freq && freq<=max_cc) {
                memcpy(&(cores[index]), &(cores[i]), freq * sizeof(simple_core));
                index += freq;
            }

            i += freq;
        }
        genome_arguments->cores_len = index;
        len = index;
    }
    
    if (genome_arguments->sct == VECTOR) {
        return;
    }

    uint64_t index = 0;
    uint64_t i=1;
    total_len += cores[0] & 0xFFFFFFFF;

    while (i<len) {
        if (cores[index] != cores[i]) {
            index++;
            cores[index] = cores[i];
            total_len += cores[i] & 0xFFFFFFFF;
        }
        i++;
    }

    if (index) {
        simple_core *new_cores = (simple_core *)realloc(cores, index * sizeof(simple_core));
        if (new_cores) {
            genome_arguments->cores = new_cores;
        } else {
            free(cores);
            genome_arguments->cores = NULL;
            index = 0;
        }
    } else {
        free(cores);
        genome_arguments->cores = NULL;
    }

    genome_arguments->cores_len = index; 
    genome_arguments->total_len = total_len;
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: File I/O operations
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

void save(FILE *out, struct lps *str) {
    int isDone = 0;

    // notify that there is an output to be written
    fwrite(&isDone, sizeof(int), 1, out);

    // write lps object
    write_lps(str, out);
}

void done(FILE *out) {
    int isDone = 1;

    // notify that there will be no output after this
    fwrite(&isDone, sizeof(int), 1, out);
}

uint64_t est_core_fq(const char *filename, int lcp_level) {

    struct stat st;

    if (stat(filename, &st) != 0) {
        log1(ERROR, "Error getting file size of %s", filename);
        return 0;
    }
    uint64_t file_size = st.st_size;

    uint64_t estimated_uncompressed_size = file_size; // default assumption
    if (strstr(filename, ".gz")) {
        // estimate uncompressed size using typical compression ratio (~4:1 for FASTQ)
        estimated_uncompressed_size = file_size * 4;
    }

    uint64_t estimated_bp_count = estimated_uncompressed_size / 2;
    return (uint64_t)(estimated_bp_count / pow(MAGIC_LCP_FQ_CONSTANT, lcp_level));
}

int ends_with_fq(const char *str) {
    size_t str_len = strlen(str);

    return  (str_len >= 3 && strcmp(str + str_len - 3, ".fq") == 0) || 
            (str_len >= 6 && strcmp(str + str_len - 6, ".fastq") == 0) || 
            (str_len >= 6 && strcmp(str + str_len - 6, ".fq.gz") == 0) || 
            (str_len >= 9 && strcmp(str + str_len - 9, ".fastq.gz") == 0);
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Logging
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

int log1(LogLevel level, const char *format, ...) {
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    printf("[%02d-%02d-%04d %02d:%02d:%02d] ", local->tm_mday, local->tm_mon + 1, local->tm_year + 1900,
           local->tm_hour, local->tm_min, local->tm_sec);

    switch (level) {
        case INFO:
            printf("[INFO] ");
            break;
        case WARN:
            printf("[WARN] ");
            break;
        case ERROR:
            printf("[ERROR] ");
            break;
    }

    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);

    printf("\n");

    return 1;
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Cleanup
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

void free_args(struct gargs * genome_arguments, struct pargs * program_arguments) {

    for (int i=0; i<program_arguments->number_of_genomes; i++) {
        if (genome_arguments[i].cores_len) {
            if (genome_arguments[i].cores_len)
                free(genome_arguments[i].cores);
            genome_arguments[i].cores_len = 0;
        }
    }

    free(genome_arguments);
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// MARK: Heap Operations
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

void heap_swap(heap_node *a, heap_node *b) {
    heap_node temp = *a;
    *a = *b;
    *b = temp;
}

void heapify_up(min_heap *heap, size_t idx) {
    while (idx > 0) {
        size_t parent = (idx - 1) / 2;
        if (heap->data[parent].value <= heap->data[idx].value)
            break;
        heap_swap(&heap->data[parent], &heap->data[idx]);
        idx = parent;
    }
}

void heapify_down(min_heap *heap, size_t idx) {
    while (2 * idx + 1 < heap->size) {
        size_t left = 2 * idx + 1;
        size_t right = 2 * idx + 2;
        size_t smallest = idx;

        if (left < heap->size && heap->data[left].value < heap->data[smallest].value)
            smallest = left;
        if (right < heap->size && heap->data[right].value < heap->data[smallest].value)
            smallest = right;
        if (smallest == idx) break;

        heap_swap(&heap->data[idx], &heap->data[smallest]);
        idx = smallest;
    }
}

void heap_push(min_heap *heap, heap_node node) {
    if (heap->size == heap->capacity) {
        heap->capacity *= 2;
        heap->data = realloc(heap->data, heap->capacity * sizeof(heap_node));
    }
    heap->data[heap->size++] = node;
    heapify_up(heap, heap->size - 1);
}

heap_node heap_pop(min_heap *heap) {
    heap_node min = heap->data[0];
    heap->data[0] = heap->data[--heap->size];
    heapify_down(heap, 0);
    return min;
}

uint64_t *merge_sorted_arrays(simple_core **cores, uint64_t *sizes, size_t file_count, uint32_t min_cc, uint32_t max_cc, uint64_t *out_total_size, double *out_total_len) {
    min_heap heap = (min_heap){malloc(file_count * sizeof(heap_node)), 0, file_count};
    uint64_t total_size = 0;
    for (size_t i = 0; i < file_count; ++i)
        total_size += sizes[i];

    uint64_t *result = malloc(total_size * sizeof(uint64_t));
    size_t result_index = 0;

    for (size_t i = 0; i < file_count; ++i) {
        if (sizes[i] > 0) {
            heap_push(&heap, (heap_node){cores[i][0], i, 0});
        }
    }

    while (heap.size > 0) {
        heap_node min = heap_pop(&heap);
        result[result_index++] = min.value;

        size_t next_idx = min.element_index + 1;
        if (next_idx < sizes[min.array_index]) {
            heap_push(&heap, (heap_node){
                cores[min.array_index][next_idx],
                min.array_index,
                next_idx
            });
        }
    }

    free(heap.data);

    // apply filtering
    uint64_t index = 0;
    uint64_t i = 0;

    while (i<total_size) {
        uint64_t freq = 1;

        for (uint64_t j=i+1; j<total_size && result[i]==result[j]; j++, freq++);

        if (min_cc<=freq && freq<=max_cc) {
            memcpy(&(result[index]), &(result[i]), freq * sizeof(simple_core));
            index += freq;
        }

        i += freq;
    }
    total_size = index;

    // remove duplicates
    if (index) {
        index = 0;
        i=1;
        double total_len = result[0] & 0xFFFFFFFF;
    
        while (i<total_size) {
            if (result[index] != result[i]) {
                index++;
                result[index] = result[i];
                total_len += result[i] & 0xFFFFFFFF;
            }
            i++;
        }
        *out_total_len = total_len;
        *out_total_size = index;
    } else {
        *out_total_len = 0;
        *out_total_size = 0;
    }

    for (size_t i = 0; i < file_count; i++) {
        if (cores[i] != NULL) free(cores[i]);
    }
    free(cores);
    free(sizes);
    return result;
}
