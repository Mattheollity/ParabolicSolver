#ifndef CONFIG_STRUCTURES_
#define CONFIG_STRUCTURES_

#include "config_structures.h"
#include "grid.h"
#include "eigen_parts.h"
#include "space_scheme.h"

struct subdomain_system
{
    spmat ii;
    spmat ib;
    spmat bb;
    spmat gb_subd;
};

struct space_scheme_config
{
    unique_ptr<space_scheme> scheme;
    unique_ptr<grid> g;
    bool try_to_import;

    space_scheme_config() :
    scheme(unique_ptr<space_scheme>(nullptr)), g(unique_ptr<grid>(nullptr)), try_to_import(false) {};

    space_scheme_config(space_scheme_config & ssc) :
    scheme(std::move(ssc.scheme)), g(std::move(ssc.g)), try_to_import(ssc.try_to_import) {};
};

#endif