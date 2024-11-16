#include "similarity_metrics.h"

#include <iostream>

void calculateIntersectionAndUnionSizes( const struct targs& argument1, const struct targs& argument2, const struct pargs& program_arguments, size_t& intersectionSize, size_t& unionSize ) {
    
    std::vector<uint32_t>::const_iterator it_core1 = argument1.cores.begin(), it_core2 = argument2.cores.begin();
    std::vector<uint32_t>::const_iterator it_end1 = argument1.cores.end(), it_end2 = argument2.cores.end();
    intersectionSize = 0;
    unionSize = 0;
    
    bool countDistinct = ( program_arguments.type == SET );

    while ( it_core1 < it_end1 && it_core2 < it_end2 ) {
        unionSize++;
        
        if ( *it_core1 == *it_core2 ) {
            intersectionSize++;
            it_core1++;
            it_core2++;

            if (countDistinct) {
                while ( it_core1 < it_end1 && *(it_core1-1) == *(it_core1)) it_core1++;
                while ( it_core2 < it_end2 && *(it_core2-1) == *(it_core2)) it_core2++;
            }
        } else if ( *it_core1 < *it_core2 ) {
            it_core1++;
            
            if (countDistinct) {
                while ( it_core1 < it_end1 && *(it_core1-1) == *(it_core1)) it_core1++;
            } else {
                it_core1++;
            }
        } else {
            it_core2++;

            if (countDistinct) {
                while ( it_core2 < it_end2 && *(it_core2-1) == *(it_core2)) it_core2++;
            }
        }
    }

    while ( it_core1 < it_end1 ) {
        unionSize++;
        it_core1++;

        if (countDistinct) {
            while ( it_core1 < it_end1 && *(it_core1-1) == *(it_core1)) it_core1++;
        }
    }

    while ( it_core2 < it_end2 ) {
        unionSize++;
        it_core2++;

        if (countDistinct) {
            while ( it_core2 < it_end2 && *(it_core2-1) == *(it_core2)) it_core2++;
        }
    }
};


double calculateJaccardSimilarity( size_t intersectionSize, size_t unionSize ) {
    return static_cast<double>(intersectionSize) / static_cast<double>(unionSize);
};


double calculateDiceSimilarity( const size_t intersectionSize, const struct targs& argument1, const struct targs& argument2 ) {
    return 2 * static_cast<double>(intersectionSize) / ( static_cast<double>(argument1.core_count) + static_cast<double>(argument2.core_count) );
};