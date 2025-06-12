#ifndef TIMING_LOG_H
#define TIMING_LOG_H

void log_timing(const char *label, double seconds);
void log_message(const char *format, ...);
void log_session_header(const char *label);
void log_session_footer(const char *label);

#endif
