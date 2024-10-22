#ifndef LOGGING_H
#define LOGGING_H

#include <stdio.h>
#include <time.h>
#include <stdarg.h>


enum LogLevel {
    INFO,
    WARN,
    ERROR
};


void timestamp();
bool log(enum LogLevel level, const char *format, ...);


#endif