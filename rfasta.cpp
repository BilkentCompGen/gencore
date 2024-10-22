#include "rfasta.h"


void read_fastas( std::vector<struct targs>& thread_arguments, const struct pargs& program_arguments ) {
    
    std::vector<std::thread> threads;
    std::vector<struct targs>::iterator current_argument = thread_arguments.begin();

    while (current_argument < thread_arguments.end()) {

        while (threads.size() < program_arguments.threadNumber && current_argument < thread_arguments.end()) {
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

    thread_arguments.size = 0;

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

                    if ( program_arguments.writeCores ) {
                        save( thread_arguments, str );
                    }

                    delete str;

                    // increment processed sequence size
                    thread_arguments.size += sequence.size();                    

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

            if ( program_arguments.writeCores ) {
                save( thread_arguments, str );
            }

            delete str;

            // increment processed sequence size
            thread_arguments.size += sequence.size();

            // cleanup
            sequence.clear();
        }
        
        file.close();
    } else {
        log(ERROR, "Could not be able to open %s", thread_arguments.inFileName.c_str());
        exit(1);
    }

    // log ending of processing fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s ended processing %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // end writing cores to file if user specified to do so
    if ( program_arguments.writeCores ) {
        done( thread_arguments );
    }

    // set lcp cores and counts to arguments
    generateSignature( thread_arguments, program_arguments );
};
