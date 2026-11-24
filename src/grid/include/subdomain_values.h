#ifndef SUBDOMAIN_VALUES_
#define SUBDOMAIN_VALUES_

#include "common.h"
#include "eigen_parts.h"
#include "config_structures_sol.h"

struct subdomain_coords
{
    vector<point> inner;
    vector<point> local_boundary;
    vector<point> global_boundary;
};

struct subdomain_values
{
    const zero_boundary_config & zbc;
    const unique_ptr<subdomain_coords> sc;

    vect inner;
    vect local_boundary;
    vect global_boundary;

    subdomain_values(const zero_boundary_config & zbc_, const subdomain_values & sbd) :
    zbc(zbc_), sc(std::make_unique<subdomain_coords>(&sbd)) {};
};

#endif