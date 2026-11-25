#include "grid_problem.h"

#include <stdexcept>

grid_problem::grid_problem(
    unique_ptr<geometry> g,
    initial_boundary_config bc
) :
grid_geometry(std::move(g)),
boundary_conditions(std::move(bc))
{
    if(!grid_geometry)
    throw std::invalid_argument("Geometry is not defined");

    if(!boundary_conditions.y0)
    throw std::invalid_argument("Initial condition is not defined");

    if(!boundary_conditions.dirichlet_condition)
    throw std::invalid_argument("Boundary condition is not defined");

    coordinates = grid_geometry->get_subdomain_coordinates();
}



const vector<subdomain_interaction> & grid_problem::get_subdomain_interactions() const
{
    return grid_geometry->get_subdomain_interactions();
}



const interaction & grid_problem::get_base_boundary_interaction() const
{
    return grid_geometry->get_boundaries_interaction();
}



const vector<int> & grid_problem::get_permutation() const
{
    return grid_geometry->get_permutation();
}



vector<subdomain_values> grid_problem::get_subdomain_values(double t) const
{
    if(t == 0.0)
    return get_initial_values();

    return get_boundary_values(t);
}



vector<subdomain_values> grid_problem::get_initial_values() const
{
    vector<subdomain_values> values = get_boundary_values(0.0);

    const initial_condition & y0 = *boundary_conditions.y0;

    for(int i=0; i<coordinates.size(); i++)
    {
        const subdomain_coordinates & c = coordinates[i];
        subdomain_values & v = values[i];

        v.inner.resize(c.inner.size());
        v.local_boundary.resize(c.boundary.size());

        for(int j=0; j<c.inner.size(); j++)
        v.inner[j] = y0(*c.inner[j]);

        for(int j=0; j<c.boundary.size(); j++)
        v.local_boundary[j] = y0(*c.boundary[j]);
    }

    return values;
}



vector<subdomain_values> grid_problem::get_boundary_values(double t) const
{
    vector<subdomain_values> values(coordinates.size());

    const dirichlet & bc = *boundary_conditions.dirichlet_condition;

    for(int i=0; i<coordinates.size(); i++)
    {
        const subdomain_coordinates & c = coordinates[i];
        subdomain_values & v = values[i];

        v.global_boundary.resize(c.global_boundary.size());

        for(int j=0; j<c.global_boundary.size(); j++)
        v.global_boundary[j] = bc(*c.global_boundary[j], t);
    }

    return values;
}