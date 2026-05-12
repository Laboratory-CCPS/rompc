#pragma once

#include <time.h>


#ifdef _MSC_VER
#define CLOCK_MONOTONIC 0

static int clock_gettime(int clock_id, struct timespec* tp) { return 0; }

static double get_duration(struct timespec ts_start, struct timespec ts_end) { return NAN; }

#else
static double get_duration(struct timespec ts_start, struct timespec ts_end) {
    double duration = ts_end.tv_sec - ts_start.tv_sec;

    duration += -1.0 + ((1000000000 - ts_start.tv_nsec) + (ts_end.tv_nsec - 0)) * 1e-9;

    return duration;
}
#endif
