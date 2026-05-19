#include "args.h"
#include "init.h"
#include "utils.h"
#include "rfasta.h"
#include "rfastq.h"
#include "rload.h"

int main(int argc, char **argv) {

    // parse and initialize arguments
    g_args_t *genome_args;
    p_args_t program_args;

    parse(argc, argv, &genome_args, &program_args);

    // initialize coefficient arrays
    LCP_INIT();

    // process files program
    switch (program_args.mode) {
    case PROGRAM_MODE_FA:
        read_fastas(genome_args, &program_args);
        break;
    case PROGRAM_MODE_FQ:
        read_fastqs(genome_args, &program_args);
        break;
    case PROGRAM_MODE_LOAD:
        read_lcpts(genome_args, &program_args);
        break;
    default:
        log1(ERROR, "Invalid program mode provided. It should not happen.");
        exit(1);
    }
    
    // calculate distances and store them in files
    calcDistances(genome_args, &program_args);

    // cleanup
    free_args(genome_args, &program_args);

    return 0;
}
