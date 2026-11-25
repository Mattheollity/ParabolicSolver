#ifndef GRID_PROBLEM_
#define GRID_PROBLEM_

#include "geometry.h"
#include "eigen_parts.h"
#include "config_structures_sol.h"

struct subdomain_values
{
    vect inner;
    vect local_boundary;
    vect global_boundary;
};



class grid_problem
{
    unique_ptr<geometry> grid_geometry;
    initial_boundary_config boundary_conditions;

    vector<subdomain_coordinates> coordinates;

    public:

    grid_problem(
        unique_ptr<geometry> g,
        initial_boundary_config bc
    );

    const vector<subdomain_interaction> & get_subdomain_interactions() const;
    const interaction & get_base_boundary_interaction() const;
    const vector<int> & get_permutation() const;

    vector<subdomain_values> get_subdomain_values(double t) const;
    vector<subdomain_values> get_initial_values() const;
    vector<subdomain_values> get_boundary_values(double t) const;
};

#endif