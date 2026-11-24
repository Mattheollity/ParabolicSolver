#include "grid_problem.h"

vector<subdomain_values> grid_problem::get_subdomain_values(double t) const
{
    if (t == 0.0)
    return get_initial_values();

    return get_boundary_values(t);
}

vector<subdomain_values> grid_problem::get_initial_values() const
{
    vector<subdomain_values> values = get_boundary_values(0.0);

    for (int i = 0; i < coordinates.size(); ++i)
    {
        values[i].inner.resize(coords.inner.size());
        values[i].local_boundary.resize(coordinates[i].local_boundary.size());

        for (int j = 0; j < coordinates[i].inner.size(); ++j)
        vals.inner[j] = boundary_conditions.y0(coordinates[i].inner[j]);

        for (int j = 0; j < coordinates[i].local_boundary.size(); ++j)
        vals.local_boundary[j] =
        boundary_conditions.y0(coordinates[i].local_boundary[j]);
    }

    return values;
}

vector<subdomain_values> grid_problem::get_boundary_values(double t) const
{
    vector<subdomain_values> values(coordinates.size());

    for(int i=0;i<coordinates.size();i++)
    values[i].global_boundary.resize(coordinates[i].global_boundary);

    for (int i=0;i<coordinates.size();i++)
    for (int j=0;j<coords.size();j++)
    values[i].global_boundary[j] =
    boundary_conditions->dirichlet_condition(coords[j], t);

    return values;
}

const vector<subdomain_interaction> & grid_problem::get_subdomain_interactions() const
{
    return subdomain_interactions;
}

const interaction & grid_problem::get_base_boundary_interaction() const
{
    return base_boundary_interaction;
}