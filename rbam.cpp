#include "rbam.h"


void read_bam( struct targs& arguments, const struct pargs& program_arguments ) {
    
    chtslib c(arguments.inFileName.c_str());

    bam1_t* aln = bam_init1();

    while( c.read(aln) > 0 ) {
        printf("%s\n", bam_get_qname(aln));
    }

    bam_destroy1(aln);

    program_arguments.verbose && log(INFO, "Successfully tested bam reading functionatly.");
}