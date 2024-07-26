#include "rfasta.h"


void read_fasta( args& arguments ) {
    
    std::fstream file;
    file.open( arguments.infilename, std::ios::in );

    int chrom_index = 0;
    std::string line;

    std::vector<lcp::lps*> strs;

    // read file
    if ( file.is_open() ) {  

        std::string sequence, id;
        sequence.reserve(250000000);
        
        while (getline(file, line)) {

            if (line[0] == '>') {

                chrom_index++;

                // process previous chromosome before moving into new one
                if (sequence.size() != 0) {
                
                    lcp::lps* str = new lcp::lps(sequence);
                    
                    for ( size_t i = 1; i < arguments.lcpLevel; i++ ) {
                        str->deepen();
                    }

                    std::cout << "Processing is done for " << id << std::endl;
                    std::cout << "Length of the processed sequence: " << sequence.size() << std::endl;
                    std::cout << "Number of cores for lcp level " << arguments.lcpLevel << ": " << str->cores.size() << std::endl;
                    std::cout << std::endl;
                    
                    strs.push_back(str);
                    
                    sequence.clear();
                }
                
                id = line.substr(1);
                std::cout << "Processing started for " << id << std::endl;
                continue;
                    
            }
            else if (line[0] != '>'){
                sequence += line;
            }
        }

        file.close();
    }
    
    flatten(strs, arguments.lcp_cores);

    if ( arguments.writeCores ) {
        save( arguments.outfilename, strs );
    }

    generateMinhashSignature(arguments.lcp_cores);
};
