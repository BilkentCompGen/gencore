#include "init.h"


void printUsage() {
    std::cout << "Usage: ./gencore [OPTIONS] [FILES]" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -r [filename]   Read stored cores from specified files (read mode)" << std::endl << std::endl;
    std::cout << "  [fa|fq|bam]     Execute program with specified files in the given format" << std::endl;
    std::cout << "                  Supported formats: [ fa | fq.gz | bam ]" << std::endl << std::endl;
    std::cout << "  -f [filename]   Execute program with a file that contains input/output file names" << std::endl << std::endl;
    std::cout << "  -l [num]        Set lcp-level. [Default: 4]" << std::endl << std::endl;
    std::cout << "  -t [num]        Set number of threads. [Default: 8]" << std::endl << std::endl;
    std::cout << "  --min-cc [num]  Set the minimum frequency (core count) for a core. [Default: 1]" << std::endl << std::endl;
    std::cout << "  --max-cc [num]  Set the maximum frequency (core count) for a core. [Default: 4294967295]" << std::endl << std::endl;
    std::cout << "  [--set|--vec]   Set program to calculate distances based or set or vector of cores. [Default: set]" << std::endl << std::endl;
    std::cout << "  -w [filename]   Store cores processed from input files." << std::endl << std::endl;
    std::cout << "  -p [prefix]     Prefix for the output of the similarity matrices results. [Default: gc]" << std::endl << std::endl;
    std::cout << "  -s [filename]   Set short names of input files. Default is first 10 characters of input file names." << std::endl << std::endl;
    std::cout << "  -v              Verbose. [Default: false]" << std::endl << std::endl;
};


void parse( int argc, char **argv, std::vector<struct targs>& all_thread_arguments, struct pargs& program_arguments ) {

    if ( argc < 2 ) {
        printUsage();
        exit(1);
    }

    // set program arguments with their default values
    program_arguments.mode = FA;
    program_arguments.type = SET;
    program_arguments.readLcpt = false;
    program_arguments.writeLcpt = false;
    program_arguments.prefix = PREFIX;
    program_arguments.threadNumber = THREAD_NUMBER;
    program_arguments.lcpLevel = 7;
    program_arguments.verbose = false;

    int index = 1;
    
    // ------------------------------------------------------------------
    // Read `read mode` 
    // ------------------------------------------------------------------
    if ( strcmp(argv[index], "-r") == 0 ) {
        program_arguments.readLcpt = true;

        // move next argument
        index++;
    }

    // ------------------------------------------------------------------
    // Read `program mode` 
    // ------------------------------------------------------------------
    if ( !program_arguments.readLcpt ) {

        // validate if following next argument exists
        if ( index >= argc ) {
            log(ERROR, "Missing program mode.");
            exit(1);
        }

        // set program mode
        if ( strcmp(argv[index], "fa") == 0 ) {
            program_arguments.mode = FA;
        } else if ( strcmp(argv[index], "fq") == 0 ) {
            program_arguments.mode = FQ;
        } else if ( strcmp(argv[index], "bam") == 0 ) {
            program_arguments.mode = BAM;
        } else {
            log(ERROR, "Invalid mode provided.");
            exit(1);
        }

        // move next argument
        index++;
    }

    // validate if input file names exists (should be given in next the argument)
    if ( index >= argc ) {
        log(ERROR, "Missing file names.");
        exit(1);
    }

    // ------------------------------------------------------------------
    // Read `input file names` 
    // ------------------------------------------------------------------
    if ( strcmp(argv[index], "-f") == 0 ) {     // if file names are given in txt file
        
        // move next argument, skip `-f`
        index++;

        // validate if following next argument exists
        if ( index >= argc ) {
            log(ERROR, "Missing file name.");
            exit(1);
        }

        // get txt file name and read its content
        std::string filename(argv[index]), line;
        std::fstream file;
        file.open( filename, std::ios::in );

        if ( file.is_open() ) {  
            while ( getline( file, line ) ) {
                struct targs args;
                args.inFileName = line;
                all_thread_arguments.push_back(args);
            }
        } else {
            log(ERROR, "Couldn't open %s", filename.c_str());
            exit(1);
        }

        file.close();
        
        // validate if at least 2 files are provided
        if ( all_thread_arguments.size() < 2 ) {
            log(ERROR, "There should be at least 2 files in %s", filename.c_str());
            exit(1);
        }
    } else {                                    // if file names are given in comma seperated format
        std::stringstream ss(argv[index]);
        std::string filename;
        
        // parse given file names w.r.t comma and set to thread arguments' inFileName
        while ( std::getline(ss, filename, ',') ) {
            struct targs args;
            args.inFileName = filename;
            all_thread_arguments.push_back(args);
        }

        // validate if at least 2 files are provided
        if ( all_thread_arguments.size() < 2 ) {
            log(ERROR, "There should be at least 2 files in %s, separated by commas.", argv[index]);
            exit(1);
        }
    }
    // move next argument
    index++;

    // Set short names' default values (first 10 characters of input file names)
    // If the file name is less than 10 characters, fill with space.
    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
        std::string name = it->inFileName;

        if ( name.size() < 10 ) {
            name.append(10 - name.size(), ' ');
        } else if ( name.size() > 10 ) {
            name = name.substr(0, 10);
        }

        it->shortName = name;
        it->min_cc = 1;
        it->max_cc = 4294967295;
    }

    // Read rest of the arguments if provided
    while ( index < argc ) {
        // ------------------------------------------------------------------
        // Read `thread number` 
        // ------------------------------------------------------------------
        if( strcmp(argv[index], "-t") == 0 ) {

            // move next argument, skip `-t`
            index++;

            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing value for thread number.");
                exit(1);
            }

            // get thread number and validate it
            try {
                if ( std::stoi(argv[index]) < 1 ) {
                    throw std::invalid_argument("Invalid thread number provided.");
                }
                program_arguments.threadNumber = std::stoi(argv[index]);            
            } catch (const std::invalid_argument& e ) {
                log(ERROR, "Invalid thread number provided.");
                exit(1);
            }

            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Read `data type` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "--set" ) == 0 || strcmp(argv[index], "--vec") == 0 ) {
            program_arguments.type = strcmp(argv[index], "--set" ) == 0 ? SET : VECTOR;
            
            // move next argument
            index++;               
        } 
        // ------------------------------------------------------------------
        // Read `LCP level` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "-l") == 0 ) {
            
            // move next argument, skip `-l`
            index++;
            
            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing value for LCP level.");
                exit(1);
            }

            // get LCP level and validate it
            try {
                if ( std::stoi(argv[index]) <= 0 ) {
                    throw std::invalid_argument("Invalid LCP level");
                }
                program_arguments.lcpLevel = std::stoi(argv[index]);
            } catch ( const std::invalid_argument& e ) {
                log(ERROR, "Invalid LCP level provided.");
                exit(1);
            }

            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Read `output file names` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "-w") == 0 ) {

            // move next argument, skip `-w`
            index++;

            // set program argument's writeLcpt mode to `true`
            program_arguments.writeLcpt = true;
            
            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing output files names.");
                exit(1);
            }
            
            if ( strcmp(argv[index], "-f") == 0 ) {         // if output file names are provided in txt file
                
                // move next argument, skip `-f`
                index++;
                
                // validate if following next argument exists
                if ( index >= argc ) {
                    log(ERROR, "Missing file name.");
                    exit(1);
                }

                // read file names from file
                std::string filename(argv[index]);        
                std::fstream file;
                file.open( filename, std::ios::in );
                
                if ( file.is_open() ) {  
                    try {
                        for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                            if ( !std::getline(file, it->outFileName) ) {
                                throw std::invalid_argument("Missing output file name.");
                            }
                        }
                    } catch ( const std::invalid_argument& e ) {
                        log(ERROR, "Number of input file should match output file names count.");
                        exit(1);
                    }
                } else {
                    log(ERROR, "Couldn't open %s", filename.c_str());
                    exit(1);
                }

                file.close();
            } else {                                        // if output file names are provided in comma seperated format
                try {
                    std::stringstream ss(argv[index]);
                    // parse given file names w.r.t comma and set to thread arguments' outFileName
                    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                        if ( !std::getline(ss, it->outFileName, ',') ) {
                            throw std::invalid_argument("Failed to parse output file name.");
                        }
                    }
                } catch ( const std::invalid_argument& e ) {
                    log(ERROR, "Number of input file should match output file names count.");
                    exit(1);
                }
            }
            
            // move next argument
            index++;
        }
        // ------------------------------------------------------------------
        // Read `prefix` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "-p") == 0 ) { 

            // move next argument, skip `-p`           
            index++;

            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing prefix.");
                exit(1);
            }

            program_arguments.prefix = argv[index];
            
            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Read `short names` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "-s") == 0 ) { 

            // move next argument, skip `-s`
            index++;
                        
            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing output files names.");
                exit(1);
            }

            if ( strcmp(argv[index], "-f") == 0 ) {     // if file names are given in txt file
                
                // move next argument, skip `-f`
                index++;

                // validate if following next argument exists
                if ( index >= argc ) {
                    log(ERROR, "Missing file name.");
                    exit(1);
                }

                // get txt file name and read its content
                std::string filename(argv[index]), line;
                std::fstream file;
                file.open( filename, std::ios::in );

                if ( file.is_open() ) {  
                    try {
                        for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                            if ( !std::getline(file, it->shortName) ) {
                                throw std::invalid_argument("Missing short name name.");
                            }

                            // make its size 10
                            if ( it->shortName.size() < 10 ) {
                                it->shortName.append(10 - it->shortName.size(), ' ');
                            } else if ( it->shortName.size() > 10 ) {
                                it->shortName = it->shortName.substr(0, 10);
                            }
                        }
                    } catch ( const std::invalid_argument& e) {
                        log(ERROR, "Number of short names should match input files' count.");
                        exit(1);
                    }
                } else {
                    log(ERROR, "Couldn't open %s", filename.c_str());
                    exit(1);
                }

                file.close();
            } else {                                    // if file names are given in comma seperated format
                try {
                    std::stringstream ss(argv[index]);
                    // parse given file names w.r.t comma and set to thread arguments' shortName
                    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                        if ( !std::getline(ss, it->shortName, ',') ) {
                            throw std::invalid_argument("Failed to parse short name.");
                        }

                        // make its size 10
                        if ( it->shortName.size() < 10 ) {
                            it->shortName.append(10 - it->shortName.size(), ' ');
                        } else if ( it->shortName.size() > 10 ) {
                            it->shortName = it->shortName.substr(0, 10);
                        }
                    }
                } catch ( const std::invalid_argument& e ) {
                    log(ERROR, "Number of short names should match input files' count.");
                    exit(1);
                }
            }

            // move next argument
            index++;
        }
        // ------------------------------------------------------------------
        // Read `min_cc` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "--min-cc") == 0 ) {
            
            // move next argument, skip `--min-cc`
            index++;
            
            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing value for minimum core count value.");
                exit(1);
            }

            if ( strcmp(argv[index], "-f") == 0 ) {         // if minimum core counts are provided in txt file
                
                // move next argument, skip `-f`
                index++;
                
                // validate if following next argument exists
                if ( index >= argc ) {
                    log(ERROR, "Missing file name.");
                    exit(1);
                }

                // read file names from file
                std::string filename(argv[index]);        
                std::fstream file;
                file.open( filename, std::ios::in );
                std::string min_cc; 
                
                if ( file.is_open() ) {  
                    try {
                        for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                            if ( !std::getline(file, min_cc) ) {
                                throw std::invalid_argument("Missing value.");
                            }

                            // get minimum core count and validate it
                            try {
                                if ( std::stoi(min_cc) <= 0 ) {
                                    throw std::invalid_argument("Invalid minimum core count value.");
                                }
                                it->min_cc = std::stoi(min_cc);
                            } catch ( const std::invalid_argument& e ) {
                                log(ERROR, "Invalid minimum core count provided.");
                                exit(1);
                            }
                        }
                    } catch ( const std::invalid_argument& e ) {
                        log(ERROR, "Number of input file should match output file names count.");
                        exit(1);
                    }
                } else {
                    log(ERROR, "Couldn't open %s", filename.c_str());
                    exit(1);
                }

                file.close();
            } else {                                        // if minimum core counts are provided as same for all
                // get minimum core count and validate it
                try {
                    if ( std::stoi(argv[index]) <= 0 ) {
                        throw std::invalid_argument("Invalid minimum core count value.");
                    }
                    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                        it->min_cc = std::stoi(argv[index]);
                    }
                } catch ( const std::invalid_argument& e ) {
                    log(ERROR, "Invalid minimum core count provided.");
                    exit(1);
                }
            }

            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Read `max_cc` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "--max-cc") == 0 ) {
            
            // move next argument, skip `--max-cc`
            index++;
            
            // validate if following next argument exists
            if ( index >= argc ) {
                log(ERROR, "Missing value for maximum core count value.");
                exit(1);
            }

            if ( strcmp(argv[index], "-f") == 0 ) {         // if minimum core counts are provided in txt file
                
                // move next argument, skip `-f`
                index++;
                
                // validate if following next argument exists
                if ( index >= argc ) {
                    log(ERROR, "Missing file name.");
                    exit(1);
                }

                // read file names from file
                std::string filename(argv[index]);        
                std::fstream file;
                file.open( filename, std::ios::in );
                std::string max_cc; 
                
                if ( file.is_open() ) {  
                    try {
                        for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                            if ( !std::getline(file, max_cc) ) {
                                throw std::invalid_argument("Missing value.");
                            }

                            // get minimum core count and validate it
                            try {
                                if ( std::stoi(max_cc) <= 0 ) {
                                    throw std::invalid_argument("Invalid maximum core count value.");
                                }
                                it->max_cc = std::stoi(max_cc);
                            } catch ( const std::invalid_argument& e ) {
                                log(ERROR, "Invalid maximum core count provided.");
                                exit(1);
                            }
                        }
                    } catch ( const std::invalid_argument& e ) {
                        log(ERROR, "Number of input file should match output file names count.");
                        exit(1);
                    }
                } else {
                    log(ERROR, "Couldn't open %s", filename.c_str());
                    exit(1);
                }

                file.close();
            } else {                                        // if minimum core counts are provided as same for all
                // get minimum core count and validate it
                try {
                    if ( std::stoi(argv[index]) <= 0 ) {
                        throw std::invalid_argument("Invalid maximum core count value.");
                    }
                    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
                        it->max_cc = std::stoi(argv[index]);
                    }
                } catch ( const std::invalid_argument& e ) {
                    log(ERROR, "Invalid maximum core count provided.");
                    exit(1);
                }
            }

            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Read `verbose` 
        // ------------------------------------------------------------------
        else if( strcmp(argv[index], "-v") == 0 ) {
            program_arguments.verbose = true;
            
            // move next argument
            index++;
        } 
        // ------------------------------------------------------------------
        // Handle INVALID OPTION 
        // ------------------------------------------------------------------
        else {
            log(ERROR, "Invalid option provided: %s", argv[index]);
            exit(1);
        }
    }

    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
        

        if ( !program_arguments.readLcpt ) {
            switch ( program_arguments.mode ) {
                case FA:
                    it->estimated_core_count = 0;
                    break;
                case FQ:
                    it->estimated_core_count = 0;
                    struct stat stat_buf;
                    if (stat(it->inFileName.c_str(), &stat_buf) == 0) {
                        size_t estimated_genome_size = static_cast<size_t>(stat_buf.st_size * COMPRESSION_RATIO) / 2;
                        it->estimated_core_count = static_cast<size_t>(estimated_genome_size * std::pow(0.45, program_arguments.lcpLevel));
                    }
                    break;
                case BAM:
                    it->estimated_core_count = 0;
                    break;
            }
        }
    }

    // Log parameters
    if( program_arguments.readLcpt ) { 
        log(INFO, "Reading cores from file.");
    } else {
        log(INFO, "Program mode: %s", ( program_arguments.mode == FA ? "FA" : ( program_arguments.mode == FQ ? "FQ" : "BAM" ) ) );
    }

    log(INFO, "Distance calculation mode: %s", ( program_arguments.type == SET ? "set" : "vector" ) );
    log(INFO, "Thread number: %d", program_arguments.threadNumber);
    log(INFO, "LCP level: %d", program_arguments.lcpLevel);
    log(INFO, "Prefix: %s", program_arguments.prefix.c_str());

    if( program_arguments.writeLcpt ) { 
        log(INFO, "Program will write cores to files.");
    }

    for ( std::vector<struct targs>::iterator it = all_thread_arguments.begin(); it < all_thread_arguments.end(); it++ ) {
        program_arguments.verbose && log(INFO, "estimated count: %ld, inFileName: %s, shortName: %s, outFileName: %s, min-cc %ld, max-cc: %ld", it->estimated_core_count, it->inFileName.c_str(), it->shortName.c_str(), it->outFileName.c_str(), it->min_cc, it->max_cc);
    }
};
