#ifndef CONFIG_STRUCTURES_SOL_
#define CONFIG_STRUCTURES_SOL_

#include "config_structures.h"
#include "solutions.h"

struct zero_boundary_config
{
    double (*y0)(const point &);
    unique_ptr<dirichlet> dirichlet_condition;

    zero_boundary_config(zero_boundary_config & zbc) :
    y0(zbc.y0), dirichlet_condition(std::move(zbc.dirichlet_condition)) {};
    zero_boundary_config() : y0(nullptr), dirichlet_condition(nullptr) {};
};

#endif