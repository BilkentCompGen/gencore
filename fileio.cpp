#include "fileio.h"


void save( std::ofstream& out, const lcp::lps* str ) {

    bool isDone = false;

    // notify that there is an output to be written
    out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
    
    // write lps object
    str->write(out);
};


void done( std::ofstream& out ) {

    bool isDone = true;

    // notify that there will be no output after this
    out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
};


void read_lcpts( std::vector<struct targs>& thread_arguments, struct pargs& program_arguments ) {

    std::vector<std::thread> threads;
    std::vector<struct targs>::iterator current_argument = thread_arguments.begin();

    while (current_argument < thread_arguments.end()) {

        while (threads.size() < program_arguments.threadNumber && current_argument < thread_arguments.end()) {
            threads.emplace_back(read_lcpt, std::ref(*current_argument), std::ref(program_arguments));
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


void read_lcpt( struct targs& thread_arguments, struct pargs& program_arguments ) {

    // get thread id
    std::ostringstream ss;
    ss << std::this_thread::get_id();

    // log initiation of reading fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s started loading %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // load lcp cores
    thread_arguments.cores.reserve( 10000000 ); // it is defalt and I don't like it
    bool isDone;

    // open file
    std::ifstream in( thread_arguments.inFileName, std::ios::binary );

    if ( !in ) {
        log(ERROR, "Error opening file for reading %s", thread_arguments.inFileName.c_str());
        return;
    }

    in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone) );

    while( !isDone ) {
        lcp::lps *str = new lcp::lps(in);
        str->deepen( program_arguments.lcpLevel );
        
        append(str, thread_arguments.cores);

        delete str;

        in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));
    } 
    
    in.close();

    // set lcp cores and counts to arguments
    generateSignature( thread_arguments, program_arguments );

    // log ending of processing fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s ended processing %s", ss.str().c_str(), thread_arguments.inFileName.c_str());
};
