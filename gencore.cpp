#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "args.h"
#include "init.h"
#include "rfasta.h"
#include "rfastq.h"
#include "rbam.h"
#include "similarity_metrics.h"


int main(int argc, char **argv) {

    // Parse and initialize arguments
    std::vector<struct targs> thread_arguments;
    struct pargs program_arguments;

    parse(argc, argv, thread_arguments, program_arguments);

    const size_t numGenomes = thread_arguments.size();

    // Initialize coefficient arrays
    lcp::init_coefficients( program_arguments.verbose );

    // Process files program
    if ( program_arguments.readCores ) {
        read_cores( thread_arguments, program_arguments );
    } else {
        switch ( program_arguments.mode ) {
        case FA:
            read_fastas( thread_arguments, program_arguments );
            break;
        case FQ:
            for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
                read_fastq( *it, program_arguments );
            }
            break;
        case BAM:
            for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
                read_bam( *it, program_arguments );
            }
            break;
        default:
            throw std::invalid_argument("Invalid program mode provided");
            exit(1);
        }
    }

    log(INFO, "Calculating distance matrices...");

    // Initialize similarity matrices
    double jaccard[numGenomes][numGenomes];
    double dice[numGenomes][numGenomes];
    
    for (size_t i = 0; i < numGenomes; i++ ) {
        jaccard[i][i] = 1;
        dice[i][i] = 1;
    }

    // Compute similarity scores
    for ( std::vector<struct targs>::iterator it1 = thread_arguments.begin(); it1 < thread_arguments.end(); it1++ ) { 
        for ( std::vector<struct targs>::iterator it2 = it1+1; it2 < thread_arguments.end(); it2++ ) {
            
            size_t interSize, unionSize;
            calculateIntersectionAndUnionSizes( *it1, *it2, program_arguments, interSize, unionSize );

            double jaccard_similarity = calculateJaccardSimilarity( interSize, unionSize );
            double dice_similarity = calculateDiceSimilarity( interSize, *it1, *it2 );

            dice[it1 - thread_arguments.begin()][it2 - thread_arguments.begin()] = dice_similarity;
            jaccard[it1 - thread_arguments.begin()][it2 - thread_arguments.begin()] = jaccard_similarity;

            // set values to transposed locations
            dice[it2 - thread_arguments.begin()][it1 - thread_arguments.begin()] = dice_similarity;
            jaccard[it2 - thread_arguments.begin()][it1 - thread_arguments.begin()] = jaccard_similarity;
        }
    }

    log(INFO, "Writing distance matrices to files...");
    
    // Write outputs to files
    std::string program_type = program_arguments.type == SET ? ".set" : ".vec";
    std::fstream dice_out, jaccard_out, distance_out;
    dice_out.open( program_arguments.prefix + program_type + ".dice.lvl" + std::to_string(program_arguments.lcpLevel) + ".phy", std::ios::out );
    jaccard_out.open( program_arguments.prefix + program_type + ".jaccard.lvl" + std::to_string(program_arguments.lcpLevel) + ".phy", std::ios::out );

    if ( dice_out.is_open() ) {  

        dice_out << numGenomes << std::endl;
        dice_out << std::setprecision(15);
        
        for(size_t i = 0; i < numGenomes; i++ ) {
            
            dice_out << thread_arguments[i].shortName;  

            for(size_t j = 0; j < numGenomes; j++ ) {
                dice_out << std::fixed << std::setprecision(15) << " " << 1-dice[i][j];
            }
            dice_out << std::endl;
        }
        dice_out.close();
    }
    
    if ( jaccard_out.is_open() ) {  

        jaccard_out << numGenomes << std::endl;
        jaccard_out << std::setprecision(15);

        for(size_t i = 0; i < numGenomes; i++ ) {
            
            jaccard_out << thread_arguments[i].shortName;  

            for(size_t j = 0; j < numGenomes; j++ ) {
                jaccard_out << std::fixed << std::setprecision(15) << " " << 1-jaccard[i][j];
            }
            jaccard_out << std::endl;
        }
        jaccard_out.close();
    }

    return 0;
};
