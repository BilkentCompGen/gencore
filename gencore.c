#include "args.h"
#include "init.h"
#include "utils.h"
#include "rfasta.h"
#include "rfastq.h"
#include "rload.h"

int main(int argc, char **argv) {

    // parse and initialize arguments
    struct gargs *genome_arguments;
    struct pargs program_arguments;

    parse(argc, argv, &genome_arguments, &program_arguments);

    // initialize coefficient arrays
    LCP_INIT();

    // process files program
    switch (program_arguments.mode) {
    case FA:
        read_fastas(genome_arguments, &program_arguments);
        break;
    case FQ:
        read_fastqs(genome_arguments, &program_arguments);
        break;
    case LOAD:
        read_lcpts(genome_arguments, &program_arguments);
        break;
    default:
        log1(ERROR, "Invalid program mode provided. It should not happen.");
        exit(1);
    }
    
    // calculate distances and store them in files
    calcDistances(genome_arguments, &program_arguments);

    // cleanup
    free_args(genome_arguments, &program_arguments);

    return 0;
}
