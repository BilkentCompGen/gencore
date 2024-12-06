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
        if ( (*it_str)->cores != nullptr ) {
            size += (*it_str)->cores->size();
        }
    }

    lcp_cores.reserve(size);

    for(std::vector<lcp::lps*>::iterator it_str = strs.begin(); it_str != strs.end(); it_str++) {
        if ( (*it_str)->cores != nullptr ) {
            for ( std::vector<lcp::core>::iterator it_lcp = (*it_str)->cores->begin(); it_lcp != (*it_str)->cores->end(); it_lcp++ ) {
                lcp_cores.push_back( (it_lcp)->label );
            }
        }
    }
};


void append(lcp::lps *str, std::vector<uint32_t>& lcp_cores) {
    if ( str->cores != nullptr ) {
        for ( std::vector<lcp::core>::iterator it_lcp = str->cores->begin(); it_lcp != str->cores->end(); it_lcp++ ) {
            lcp_cores.push_back( it_lcp->label );
        }
    }
};


void generateSignature( struct targs& thread_arguments, const struct pargs& program_arguments ) {
    
    std::sort(thread_arguments.cores.begin(), thread_arguments.cores.end());
    
    size_t count = 1;
    std::vector<uint32_t>::iterator write_index = thread_arguments.cores.begin(), temp;

    for ( std::vector<uint32_t>::const_iterator it = thread_arguments.cores.begin() + 1; it < thread_arguments.cores.end(); it++ ) {
        if ( *it == *(it - 1) ) {
            count++;
        } else {
            if ( count >= thread_arguments.min_cc && count <= thread_arguments.max_cc ) {
                while ( count-- ) {
                    *write_index++ = *(it - 1);
                }
            }
            count = 1;
        }
    }

    thread_arguments.cores.resize(write_index - thread_arguments.cores.begin());

    if ( program_arguments.type == SET ) {
        count = 1;
        for( std::vector<uint32_t>::const_iterator it = thread_arguments.cores.begin() + 1; it < thread_arguments.cores.end(); it++ ) {
            *(it-1) != *it && count++;
        }
        thread_arguments.core_count = count;
    } else {
        thread_arguments.core_count = thread_arguments.cores.size();
    }    
};
