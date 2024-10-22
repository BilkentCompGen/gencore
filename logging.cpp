#include "logging.h"


void timestamp() {
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    printf("[%02d-%02d-%04d %02d:%02d:%02d] ", local->tm_mday, local->tm_mon + 1, local->tm_year + 1900,
           local->tm_hour, local->tm_min, local->tm_sec);
};


bool log(enum LogLevel level, const char *format, ...) {
    timestamp();

    switch (level) {
        case INFO:
            printf("[INFO] ");
            break;
        case WARN:
            printf("[WARN] ");
            break;
        case ERROR:
            printf("[ERROR] ");
            break;
    }

    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);

    printf("\n");

    return true;
};
