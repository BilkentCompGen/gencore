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

double calcHammDist(double diceSimm, double avgLen) {
    return 1 - pow(diceSimm, 1.0/avgLen);
}

double calcJukesCantorCor(double hammingDist) {
    return -3.0/4.0 * log(1-hammingDist * 3.0/4.0);
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

            double diceSim = 1.0 - calcDiceSim(interSize, genome_arguments[i].cores_len, genome_arguments[j].cores_len);
            double jaccardSim = 1.0 - calcJaccardSim(interSize, unionSize);
            double avg_len = (genome_arguments[i].total_len+genome_arguments[j].total_len)/(genome_arguments[i].cores_len+genome_arguments[j].cores_len);
            double jukesCantorSim = calcHammDist(diceSim, avg_len);
            jukesCantorSim = calcJukesCantorCor(jukesCantorSim);

            dice[i][j] = diceSim;
            jaccard[i][j] = jaccardSim;
            jukes_cantor[i][j] = jukesCantorSim;

            // set values to transposed locations
            dice[j][i] = diceSim;
            jaccard[j][i] = jaccardSim;
            jukes_cantor[j][i] = jukesCantorSim;
        }
    }

    log1(INFO, "Writing distance matrices to files...");

    int lcp_level = genome_arguments[0].lcp_level;
    
    // Write outputs to files
    char *program_type = genome_arguments[0].sct == SET ? "set" : "vec";
    FILE *dice_out, *jaccard_out, *jukes_cantor_out;
    char filename_buffer[256];
    if (snprintf(filename_buffer, 256, "%s.%s.%s%03d.phy", program_arguments->prefix, program_type, "dice.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for dice overflow.");
        exit(EXIT_FAILURE);
    }
    dice_out = fopen(filename_buffer, "w");

    if (snprintf(filename_buffer, 256, "%s.%s.%s%03d.phy", program_arguments->prefix, program_type, "jaccard.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for jaccard overflow.");
        exit(EXIT_FAILURE);
    }
    jaccard_out = fopen(filename_buffer, "w");

    if (snprintf(filename_buffer, 256, "%s.%s.%s%03d.phy", program_arguments->prefix, program_type, "jc.lvl", lcp_level) < 0) {
        log1(ERROR, "Filename buffer for jc overflow.");
        exit(EXIT_FAILURE);
    }
    jukes_cantor_out = fopen(filename_buffer, "w");

    // write dice
    if (dice_out) {
        fprintf(dice_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(dice_out, "%10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(dice_out, " %.15f", dice[i][j]);
            }
            fprintf(dice_out, "\n");
        }
        fclose(dice_out);
    }

    // write jaccard
    if (jaccard_out) {
        fprintf(jaccard_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(jaccard_out, "%10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(jaccard_out, " %.15f", jaccard[i][j]);
            }
            fprintf(jaccard_out, "\n");
        }
        fclose(jaccard_out);
    }

    // write jukes cantor
    if (jukes_cantor_out) {
        fprintf(jukes_cantor_out, "%d\n", program_arguments->number_of_genomes);

        for (int i = 0; i < program_arguments->number_of_genomes; i++) {
            fprintf(jukes_cantor_out, "%10s", genome_arguments[i].shortName);

            for (int j = 0; j < program_arguments->number_of_genomes; j++) {
                fprintf(jukes_cantor_out, " %.15f", jukes_cantor[i][j]);
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

void genSign(struct gargs *genome_arguments, sim_calculation_type mode) {

    simple_core *cores = genome_arguments->cores;
    uint64_t len = genome_arguments->cores_len;
    double totalLen = genome_arguments->total_len;

    quicksort(cores, 0, len);

    if (mode == VECTOR) {
        return;
    }

    uint64_t index = 0;
    for (uint64_t i=1; i<len; i++) {
        if (cores[index] != cores[i]) {
            index++;
            cores[index] = cores[i];
            totalLen += cores[i] & 0xFFFFFFFF;
        }
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
    genome_arguments->total_len = totalLen;
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
            free(genome_arguments[i].cores);
            genome_arguments[i].cores_len = 0;
        }
    }

    free(genome_arguments);
}