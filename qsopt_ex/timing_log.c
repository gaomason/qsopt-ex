#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>  

// global log filename
char log_filename[256] = "qsopt_timing.log";

// sets the log filename based on problem name
void set_log_file(const char *problem_name) {
    snprintf(log_filename, sizeof(log_filename), "%s.log", problem_name);
}

// creates the log file
void log_timing(const char *label, double seconds) {
    FILE *fp = fopen(log_filename, "a");
    if (fp) {
        time_t now = time(NULL);
        char *timestamp = ctime(&now);
        timestamp[strcspn(timestamp, "\n")] = 0;  // strips the newline created by previous live
        fprintf(fp, "[%s] %s%.10f seconds\n", timestamp, label, seconds);
        fclose(fp);
    }
}

// general log line writer
void log_message(const char *format, ...) {
    FILE *fp = fopen(log_filename, "a");
    if (fp) {
        va_list args;
        va_start(args, format);
        vfprintf(fp, format, args);
        va_end(args);
        fprintf(fp, "\n");
        fclose(fp);
    }
}

// session header
void log_session_header(const char *label) {
    FILE *fp = fopen(log_filename, "a");
    if (fp) {
        time_t now = time(NULL);
        char *timestamp = ctime(&now);
        timestamp[strcspn(timestamp, "\n")] = 0;

        fprintf(fp, "============================================================\n");
        fprintf(fp, "LOG SESSION START — %s\n", timestamp);
        fprintf(fp, "Solving Problem: %s\n", label);

        fclose(fp);
    }
}

// session footer
void log_session_footer(const char *label) {
    FILE *fp = fopen(log_filename, "a");
    if (fp) {
        time_t now = time(NULL);
        char *timestamp = ctime(&now);
        timestamp[strcspn(timestamp, "\n")] = 0;

        fprintf(fp, "Solved Problem: %s\n", label);
        fprintf(fp, "LOG SESSION END — %s\n", timestamp);
        fprintf(fp, "============================================================\n");
        fprintf(fp, "\n");
        fclose(fp);
    }
}