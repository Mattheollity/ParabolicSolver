#ifndef GRID_PROBLEM_
#define GRID_PROBLEM_

#include "interactions.h"
#include "eigen_parts.h"
#include "config_structures_sol.h"

struct subdomain_coordinates
{
    vector<point> inner;
    vector<point> local_boundary;
    vector<point> global_boundary;
};

struct subdomain_values
{
    vect inner;
    vect local_boundary;
    vect global_boundary;
};

class grid_problem
{
    zero_boundary_config boundary_conditions;

    vector<subdomain_coordinates> coordinates;
    vector<subdomain_interaction> subdomain_interactions;
    interaction base_boundary_interaction;

    public:

    const vector<subdomain_interaction> & get_subdomain_interactions() const;
    const interaction & get_base_boundary_interaction() const;
    vector<subdomain_values> get_subdomain_values(double t) const;
    vector<subdomain_values> get_initial_values() const;
    vector<subdomain_values> get_boundary_values(double t) const;
};

#endif