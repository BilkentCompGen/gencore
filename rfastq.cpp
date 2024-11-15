#include "rfastq.h"


void read_fastqs( std::vector<struct targs>& thread_arguments, const struct pargs& program_arguments ) {
    
    std::vector<std::thread> threads;
    std::vector<struct targs>::iterator current_argument = thread_arguments.begin();

    while (current_argument < thread_arguments.end()) {

        while (threads.size() < program_arguments.threadNumber && current_argument < thread_arguments.end()) {
            threads.emplace_back(read_fastq, std::ref(*current_argument), std::ref(program_arguments));
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


void read_fastq( struct targs& thread_arguments, const struct pargs program_arguments ) {

    // get thread id
    std::ostringstream ss;
    ss << std::this_thread::get_id();
    
    GzFile infile( thread_arguments.inFileName.c_str(), "rb" );

    program_arguments.verbose && log(INFO, "Thread ID: %s started processing %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    if ( infile.is_open() ) {

        // variables
        char buffer[BUFFERSIZE];

        // read file
        std::string sequence;
            
        sequence.reserve(BUFFERSIZE);

        while (true) {

            if ( infile.gets(buffer) == Z_NULL) {
                // End of file or an error
                if ( ! infile.eof() ) {
                    std::cerr << "Error reading file." << std::endl;
                }
                break;
            }

            infile.gets(buffer);
            sequence.append(buffer);

            lcp::lps *str = new lcp::lps(sequence);
            str->deepen(program_arguments.lcpLevel);
            
            if ( str->cores != nullptr ) {
                for ( std::vector<lcp::core>::iterator it = str->cores->begin(); it != str->cores->end(); it++ ) {
                    thread_arguments.cores.push_back( (it)->label );
                }
            }

            if ( program_arguments.writeCores ) {
                save( thread_arguments, str );
            }

            delete str;
            str = nullptr;

            str = new lcp::lps(sequence, false, true);
            str->deepen(program_arguments.lcpLevel);

            if ( str->cores != nullptr ) { 
                for ( std::vector<lcp::core>::iterator it = str->cores->begin(); it != str->cores->end(); it++ ) {
                    thread_arguments.cores.push_back( (it)->label );
                }
            }

            if ( program_arguments.writeCores ) {
                save( thread_arguments, str );
            }

            delete str;
            str = nullptr;
            
            infile.gets(buffer);
            infile.gets(buffer);

            sequence.clear();
        }
    } else {
        log(ERROR, "Could not be able to open %s", thread_arguments.inFileName.c_str());
        exit(1);
    }

    program_arguments.verbose && log(INFO, "Thread ID: %s ended processing %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // end writing cores to file if user specified to do so
    if ( program_arguments.writeCores ) {
        done( thread_arguments );
    }

    generateSignature( thread_arguments, program_arguments );
};