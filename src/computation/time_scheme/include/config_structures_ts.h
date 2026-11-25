#ifndef CONFIG_STRUCTURES_TS_
#define CONFIG_STRUCTURES_TS_

#include "common.h"
#include "time_scheme.h"


struct time_scheme_config
{
    double start_time;
    double end_time;

    bool try_to_import = false;

    unique_ptr<time_scheme> scheme;
};

#endif