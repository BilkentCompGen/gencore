#include "fileio.h"


void done(  const struct targs& thread_arguments ) {
    std::ofstream out(thread_arguments.outFileName, std::ios::binary | std::ios::app);

    if (!out) {
        log(ERROR, "Error opening file for writing %s", thread_arguments.inFileName.c_str());
        return;
    }

    bool isDone = true;

    out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
    out.write(reinterpret_cast<const char*>(&(thread_arguments.size)), sizeof(thread_arguments.size));
    
    out.close();
};


void save( const struct targs& thread_arguments, const lcp::lps* str ) {
    
    std::ofstream out(thread_arguments.outFileName, std::ios::binary | std::ios::app);

    if (!out) {
        log(ERROR, "Error opening file for writing %s", thread_arguments.inFileName.c_str());
        return;
    }

    bool isDone = false;

    out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
    str->write(out);

    out.close();
};


void save( const struct targs& thread_arguments, const std::vector<lcp::lps*>& cores ) {
    std::ofstream out(thread_arguments.outFileName, std::ios::binary);
    if (!out) {
        log(ERROR, "Error opening file for writing %s", thread_arguments.inFileName.c_str());
        return;
    }

    log(INFO, "Saving cores to file %s", thread_arguments.outFileName.c_str());

    size_t isDone = false;

    for ( std::vector<lcp::lps*>::const_iterator it = cores.begin(); it != cores.end(); it++ ) {
        out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
        (*it)->write(out);
    }

    isDone = true;

    out.write(reinterpret_cast<const char*>(&isDone), sizeof(isDone));
    out.write(reinterpret_cast<const char*>(&(thread_arguments.size)), sizeof(thread_arguments.size));

    out.close();
};


void load( struct targs& thread_arguments, struct pargs& program_arguments, std::vector<lcp::lps*>& cores ) {
    std::ifstream in(thread_arguments.inFileName, std::ios::binary);

    if (!in) {
        log(ERROR, "Error opening file for reading %s", thread_arguments.inFileName.c_str());
        return;
    }
    
    bool isDone;

    in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));

    while( !isDone ) {
        lcp::lps *str = new lcp::lps(in);
        str->deepen(program_arguments.lcpLevel);

        cores.push_back(str);

        in.read(reinterpret_cast<char*>(&isDone), sizeof(isDone));
    } 

    in.read(reinterpret_cast<char*>(&(thread_arguments.size)), sizeof(thread_arguments.size));

    in.close();
};


void read_from_file( struct targs& thread_arguments, struct pargs& program_arguments ) {

    // get thread id
    std::ostringstream ss;
    ss << std::this_thread::get_id();

    // log initiation of reading fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s started loading %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // load lcp cores
    std::vector<lcp::lps*> strs;
    load( thread_arguments, program_arguments, strs );

    // log ending of processing fasta
    program_arguments.verbose && log(INFO, "Thread ID: %s ended loading %s", ss.str().c_str(), thread_arguments.inFileName.c_str());

    // get lcp core hashes
    flatten(strs, thread_arguments.cores );

    // delete lcp cores
    for ( std::vector<lcp::lps*>::iterator it = strs.begin(); it != strs.end(); it++ ) {
        delete (*it);
    }
    strs.clear();

    // set lcp cores and counts to arguments
    generateSignature( thread_arguments, program_arguments );
};


void read_cores( std::vector<struct targs>& thread_arguments, struct pargs& program_arguments ) {

    std::vector<std::thread> threads;
    std::vector<struct targs>::iterator current_argument = thread_arguments.begin();

    while (current_argument < thread_arguments.end()) {

        while (threads.size() < program_arguments.threadNumber && current_argument < thread_arguments.end()) {
            threads.emplace_back(read_from_file, std::ref(*current_argument), std::ref(program_arguments));
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
