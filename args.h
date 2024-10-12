#ifndef ARGS_H
#define ARGS_H

#include <cstdint>
#include <string>
#include <vector>
#include "program_mode.h"


struct pargs {
    program_mode mode;
    data_type type;
    bool readCores;
    bool writeCores;
    std::string prefix;
    size_t threadNumber;
    int lcpLevel;
    bool verbose;
};


struct targs {
    std::string inFileName;
    std::string outFileName;
    std::string shortName;
    std::vector<uint32_t> cores;
    std::vector<size_t> counts;
    size_t size;
};


#endif