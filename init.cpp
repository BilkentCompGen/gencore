#include "init.h"


void printUsage() {
    std::cout << "Usage: ./gencore [OPTIONS] [FILES]" << std::endl << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -r              Read cores from specified files (read mode)" << std::endl;
    std::cout << "                  Usage: ./gencore -r file1.cores,file2.cores" << std::endl;
    std::cout << "                  Usage: ./gencore -r -f files.txt" << std::endl << std::endl;
    std::cout << "  [fa|fq|bam]     Execute program with specified files in the given format" << std::endl;
    std::cout << "                  Supported formats: [ fa | fq.gz | bam ]" << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa" << std::endl;
    std::cout << "                         ./gencore fq reads1.fq.gz,reads2.fq.gz" << std::endl;
    std::cout << "                         ./gencore bam aln1.bam,aln2.bam" << std::endl << std::endl;
    std::cout << "  -f [filename]   Execute program with a file that contains input/output file names" << std::endl;
    std::cout << "                  Usage: ./gencore fa -f files.txt" << std::endl << std::endl;
    std::cout << "  -l [level]      Set lcp-level. [Default: 4]" << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa -l 4" << std::endl << std::endl;
    std::cout << "  -t [number]     Set number of threads. [Default: 8]" << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa -t 2" << std::endl << std::endl;
    std::cout << "  [--set|--vec]   Set program to calculate distances based or set or vector of cores. [Default: set]" << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa --set" << std::endl << std::endl;
    std::cout << "  -w [filenames]  Store cores processed from input files." << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa -w -f files.txt" << std::endl << std::endl;
    std::cout << "                         ./gencore fa ref1.fa,ref2.fa -w ref1.cores,ref2.cores" << std::endl;
    std::cout << "  -p [prefix]     Prefix for the output of the similarity matrices results. [Default: gc]" << std::endl;
    std::cout << "                  Usage ./gencore fa -i infiles.txt -o outfiles.txt -p primates" << std::endl << std::endl;
    std::cout << "  -s [shortnames] Set short names of input files. Default is first 10 characters of input file names." << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa -s -f files.txt" << std::endl << std::endl;
    std::cout << "                         ./gencore fa ref1.fa,ref2.fa -s ref1,ref2" << std::endl;
    std::cout << "  -v              Verbose. [Default: false]" << std::endl;
    std::cout << "                  Usage: ./gencore fa ref1.fa,ref2.fa -v" << std::endl << std::endl;
};


void parse( int argc, char **argv, std::vector<struct targs>& thread_arguments, struct pargs& program_arguments ) {

    if ( argc < 2 ) {
        printUsage();
        exit(1);
    }

    // set program arguments with their default values
    program_arguments.mode = FA;
    program_arguments.type = SET;
    program_arguments.readCores = false;
    program_arguments.writeCores = false;
    program_arguments.prefix = PREFIX;
    program_arguments.threadNumber = THREAD_NUMBER;
    program_arguments.lcpLevel = 7;
    program_arguments.verbose = false;
    program_arguments.min_cc = 1;
    program_arguments.max_cc = 4294967295;

    int index = 1;
    
    // ------------------------------------------------------------------
    // Read `read mode` 
    // ------------------------------------------------------------------
    if ( strcmp(argv[index], "-r") == 0 ) {
        program_arguments.readCores = true;

        // move next argument
        index++;
    }

    // ------------------------------------------------------------------
    // Read `program mode` 
    // ------------------------------------------------------------------
    if ( !program_arguments.readCores ) {

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
                thread_arguments.push_back(args);
            }
        } else {
            log(ERROR, "Couldn't open %s", filename.c_str());
            exit(1);
        }

        file.close();
        
        // validate if at least 2 files are provided
        if ( thread_arguments.size() < 2 ) {
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
            thread_arguments.push_back(args);
        }

        // validate if at least 2 files are provided
        if ( thread_arguments.size() < 2 ) {
            log(ERROR, "There should be at least 2 files in %s, separated by commas.", argv[index]);
            exit(1);
        }
    }
    // move next argument
    index++;

    // Set short names' default values (first 10 characters of input file names)
    // If the file name is less than 10 characters, fill with space.
    for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
        std::string name = it->inFileName;

        if ( name.size() < 10 ) {
            name.append(10 - name.size(), ' ');
        } else if ( name.size() > 10 ) {
            name = name.substr(0, 10);
        }

        it->shortName = name;
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

            // set program argument's writeCores mode to `true`
            program_arguments.writeCores = true;
            
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
                        for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
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
                    for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
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
                        for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
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
                    for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
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

            // get minimum core count and validate it
            try {
                if ( std::stoi(argv[index]) <= 0 ) {
                    throw std::invalid_argument("Invalid minimum core count value.");
                }
                program_arguments.min_cc = std::stoi(argv[index]);
            } catch ( const std::invalid_argument& e ) {
                log(ERROR, "Invalid minimum core count provided.");
                exit(1);
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

            // get maximum core count and validate it
            try {
                if ( std::stoi(argv[index]) <= 0 ) {
                    throw std::invalid_argument("Invalid maximum core count value.");
                }
                program_arguments.max_cc = std::stoi(argv[index]);
            } catch ( const std::invalid_argument& e ) {
                log(ERROR, "Invalid maximum core count provided.");
                exit(1);
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

    // Log parameters
    if( program_arguments.readCores ) { 
        log(INFO, "Reading cores from file.");
    } else {
        log(INFO, "Program mode: %s", ( program_arguments.mode == FA ? "FA" : ( program_arguments.mode == FQ ? "FQ" : "BAM" ) ) );
    }

    log(INFO, "Distance calculation mode: %s", ( program_arguments.type == SET ? "set" : "vector" ) );
    log(INFO, "Thread number: %d", program_arguments.threadNumber);
    log(INFO, "LCP level: %d", program_arguments.lcpLevel);
    log(INFO, "Prefix: %s", program_arguments.prefix.c_str());

    if( program_arguments.writeCores ) { 
        log(INFO, "Program will write cores to files.");
    }

    for ( std::vector<struct targs>::iterator it = thread_arguments.begin(); it < thread_arguments.end(); it++ ) {
        program_arguments.verbose && log(INFO, "inFileName: %s, shortName: %s, outFileName: %s", it->inFileName.c_str(), it->shortName.c_str(), it->outFileName.c_str());
    }
};
