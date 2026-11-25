#ifndef MATRIX_ASSEMBLER_
#define MATRIX_ASSEMBLER_

#include "common.h"
#include "config_structures_ss.h"
#include "space_scheme.h"
#include "eigen_parts.h"
#include "grid_problem.h"


class space_matrix_assembler_io;

class space_matrix_assembler
{
    friend space_matrix_assembler_io;

    const grid_problem & gp;

    list<triplet> base_boundary_system;
    vector<subdomain_system> subdomain_systems;

    void build_subdomain_systems();
    void build_subdomain_system(int i, const subdomain_interaction &);
    void build_base_boundary_system();
    void build_spmat_from_interaction(const interaction &, spmat &);



    public:

    const space_scheme_config & ssc;

    space_matrix_assembler(const grid_problem &, const space_scheme_config &);

    const list<triplet> & get_base_boundary_system() const {return base_boundary_system;}
    const vector<subdomain_system> & get_subdomain_systems() const {return subdomain_systems;}
};

#endif