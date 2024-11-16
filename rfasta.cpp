#include "rfasta.h"


void read_fastas( std::vector<struct targs>& all_thread_arguments, const struct pargs& program_arguments ) {
    
    std::vector<std::thread> threads;
    std::vector<struct targs>::iterator current_argument = all_thread_arguments.begin();

    while (current_argument < all_thread_arguments.end()) {

        while (threads.size() < program_arguments.threadNumber && current_argument < all_thread_arguments.end()) {
            threads.emplace_back(read_fasta, std::ref(*current_argument), std::ref(program_arguments));
            current_argument++;
        }

        for ( std::vector<std::thread>::iterator it = threads.begin(); it != threads.end(); ) {
            if (it->joinable()) {
                it->join();
                it = threads.erase(it);
            } else {
                it++;
            }
        }
    }

    for (std::vector<std::thread>::iterator it = threads.begin(); it != threads.end(); it++ ) {
        it->join();
    }
};


void read_fasta( struct targs& thread_arguments, const struct pargs& program_arguments ) {
    
    // get thread id
    std::ostringstream ss;
    ss << std::this_thread::get_id();

    // log initiation of reading fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s started processing %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // open fasta file
    std::fstream file;
    file.open( thread_arguments.inFileName, std::ios::in );

    // create file for writing cores
    std::ofstream out; 

    if ( program_arguments.writeLcpt ) {
        out.open(thread_arguments.outFileName, std::ios::binary);

        if ( !out ) {
            log(ERROR, "Error opening file for saving into file %s", thread_arguments.outFileName.c_str());
            return;
        }
    }

    // read file
    if ( file.is_open() ) {
        
        std::string sequence, id, line;
        sequence.reserve(250000000);
        
        while (getline(file, line)) {

            if (line[0] == '>') {

                // process previous chromosome before moving into new one
                if (sequence.size() != 0) {
                
                    lcp::lps* str = new lcp::lps(sequence);
                    str->deepen( program_arguments.lcpLevel );
                    str->get_labels( thread_arguments.cores );

                    if ( program_arguments.writeLcpt ) {
                        save( out, str );
                    }

                    delete str;                  

                    // cleanup
                    sequence.clear();
                }
                
                // get new chromosome's id
                id = line.substr(1);
            }
            else if (line[0] != '>'){
                sequence += line;
            }
        }

        // process last chromosome set into sequence string
        if ( sequence.size() != 0 ) {

            lcp::lps* str = new lcp::lps(sequence);
            str->deepen( program_arguments.lcpLevel );
            str->get_labels( thread_arguments.cores );

            if ( program_arguments.writeLcpt ) {
                save( out, str );
            }

            delete str;
            
            // cleanup
            sequence.clear();
        }
        
        file.close();
    } else {
        log(ERROR, "Could not be able to open %s", thread_arguments.inFileName.c_str());
        exit(1);
    }

    // end writing cores to file if user specified to do so
    if ( program_arguments.writeLcpt ) {
        done( out );
        out.close();
    }

    // sort and filter the cores
    generateSignature( thread_arguments, program_arguments );

    // log ending of processing fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s ended processing %s, size: %ld", ss.str().c_str(), thread_arguments.inFileName.c_str(), thread_arguments.cores.size());

};
