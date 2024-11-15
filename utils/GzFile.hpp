/**
 * @file    GzFile.hpp
 * @brief   Wrapper Class for Handling Gzip-Compressed Files
 *
 * This header file defines the GzFile class, which serves as a wrapper for the zlib
 * library, facilitating the handling of gzip-compressed files in a C++ program. 
 * The class provides methods to open, read from, write to, and close gzip files. 
 * It includes functionalities such as checking for the end of file, reading lines, 
 * rewinding the file, and printing formatted data to the compressed file.
 *
 * The GzFile class is designed to simplify operations on gzip files, making it
 * easier to integrate file compression/decompression in genomic data processing
 * or other applications that require handling of large compressed files.
 *
 * Usage Example:
 *     GzFile gzFile("example.gz", "wb");
 *     gzFile.printf("Hello World\n");
 *     // Other operations...
 *     gzFile.~GzFile();
 */


#ifndef GZFILE_HPP
#define GZFILE_HPP

#include <zlib.h>

#ifndef BUFFERSIZE
#define BUFFERSIZE      500000
#endif


class GzFile {
public:
    GzFile(const char* filename, const char* mode) {
        gzFile_ = gzopen(filename, mode);
    }

    ~GzFile() {
        if (gzFile_) gzclose(gzFile_);
    }

    // Returns true if end of file has been reached
    bool eof() const {
        return gzFile_ && gzeof(gzFile_);
    }

    // Wrapper method for gzgets
    char* gets(char* buffer) {
        return gzFile_ ? gzgets(gzFile_, buffer, BUFFERSIZE) : Z_NULL;
    }

    // Rewind file
    void rewind(){
        gzrewind(gzFile_);
    }

    // Is open
    bool is_open() const {
        return gzFile_ != nullptr;
    }

    // Print to compressed file
    int printf(const char* format, ...) {
        char buffer[BUFFERSIZE];
        va_list args;
        va_start(args, format);
        int ret = vsnprintf(buffer, BUFFERSIZE, format, args);
        va_end(args);
        if (ret > 0 && ret < BUFFERSIZE) {
            gzwrite(gzFile_, buffer, ret);
        }
        return ret;
    }

    // Check if the file is valid
    explicit operator bool() const {
        return gzFile_ != nullptr;
    }

private:
    gzFile gzFile_;
};


#endif