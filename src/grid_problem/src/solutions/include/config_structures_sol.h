#ifndef CONFIG_STRUCTURES_SOL_
#define CONFIG_STRUCTURES_SOL_

#include "common.h"

#include "initial_conditions.h"
#include "boundary_conditions.h"


struct initial_boundary_config
{
    unique_ptr<initial_condition> y0;
    unique_ptr<dirichlet> dirichlet_condition;
};

#endif