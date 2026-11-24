#ifndef CONFIG_STRUCTURES_TS_
#define CONFIG_STRUCTURES_TS_

#include "config_structures.h"
#include "time_scheme.h"

struct time_scheme_config
{
    double start_time;
    double end_time;
    double T;
    double theta;
    bool try_to_import;
    shared_ptr<time_scheme> scheme;

    time_scheme_config(const time_scheme_config & tsc)
    : start_time(tsc.start_time), end_time(tsc.end_time), T(tsc.T),
    theta(tsc.theta), try_to_import(tsc.try_to_import), scheme(tsc.scheme) {}
};

#endif