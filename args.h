#ifndef ARGS_H
#define ARGS_H

#include <cstdint>
#include <string>
#include <vector>
#include "program_mode.h"


struct pargs {
    program_mode mode;
    data_type type;
    bool readLcpt;
    bool writeLcpt;
    std::string prefix;
    size_t threadNumber;
    int lcpLevel;
    bool verbose;
    uint32_t min_cc;
    uint32_t max_cc;
};


struct targs {
    std::string inFileName;
    std::string outFileName;
    std::string shortName;
    std::vector<uint32_t> cores;
    size_t estimated_core_count;
    size_t core_count;
};


#endif