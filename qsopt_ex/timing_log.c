#include <stdio.h>
#include <time.h>
#include <stdarg.h>

#define LOG_FILE "qsopt_timing.log"

// creates the log file
void log_timing(const char *label, double seconds) {
    FILE *fp = fopen(LOG_FILE, "a");
    if (fp) {
        time_t now = time(NULL);
        char *timestamp = ctime(&now);  
        timestamp[strcspn(timestamp, "\n")] = 0;  // strips the extra newline provided by the line above
        fprintf(fp, "[%s] %s took %.6f seconds\n", timestamp, label, seconds);
        fclose(fp);
    }
}

// general log line writer
void log_message(const char *format, ...) {
    FILE *fp = fopen(LOG_FILE, "a");
    if (fp) {
        va_list args;
        va_start(args, format);
        vfprintf(fp, format, args);
        va_end(args);
        fprintf(fp, "\n");
        fclose(fp);
    }
}
