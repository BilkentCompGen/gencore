#include "helper.h"


bool reverseComplement( std::string& sequence ) {
    
    // reverse the sequence
    std::reverse( sequence.begin(), sequence.end() );
    
    // replace each nucleotide with its complement
    std::transform( sequence.begin(), sequence.end(), sequence.begin(),
                   [](char nucleotide) -> char {
                       switch (nucleotide) {
                           case 'A': return 'T';
                           case 'T': return 'A';
                           case 'G': return 'C';
                           case 'C': return 'G';
                           default: return nucleotide;
                       }
                   });
    
    return true;
};


void flatten(std::vector<lcp::lps*>& strs, std::vector<uint32_t>& lcp_cores) {
    
    size_t size = 0;
    
    for(std::vector<lcp::lps*>::iterator it_str = strs.begin(); it_str != strs.end(); it_str++) {
        size += (*it_str)->cores->size();
    }

    lcp_cores.reserve(size);

    for(std::vector<lcp::lps*>::iterator it_str = strs.begin(); it_str != strs.end(); it_str++) {
        for ( std::vector<lcp::core*>::iterator it_lcp = (*it_str)->cores->begin(); it_lcp != (*it_str)->cores->end(); it_lcp++ ) {
            lcp_cores.push_back( (*it_lcp)->label );
        }
    }
};


void generateSignature( struct targs& thread_arguments, const struct pargs& program_arguments ) {
    
    std::sort(thread_arguments.cores.begin(), thread_arguments.cores.end());

    if ( program_arguments.type == SET ) {
        size_t size = 1;
        for( std::vector<uint32_t>::const_iterator it = thread_arguments.cores.begin() + 1; it < thread_arguments.cores.end(); it++ ) {
            *(it-1) != *it && size++;
        }
        thread_arguments.core_size = size;
    } else {
        thread_arguments.core_size = thread_arguments.cores.size();
    }    
};
