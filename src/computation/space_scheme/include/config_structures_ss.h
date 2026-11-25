#ifndef CONFIG_STRUCTURES_
#define CONFIG_STRUCTURES_

#include "config_structures.h"
#include "eigen_parts.h"
#include "space_scheme.h"

struct subdomain_system
{
    spmat ii;
    spmat ib;
    spmat bi;
    spmat bb;
    spmat gb_subd;
};

struct space_scheme_config
{
    unique_ptr<space_scheme> scheme;
    bool try_to_import;
};

#endif