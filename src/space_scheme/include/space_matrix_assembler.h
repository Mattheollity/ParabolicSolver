#ifndef MATRIX_ASSEMBLER_
#define MATRIX_ASSEMBLER_

#include "common.h"
#include "config_structures.h"
#include "space_scheme.h"
#include <sstream>

struct subdomain_system
{
    spmat ii;
    spmat ib;
    spmat bb;
    spmat gb_subd;
    vector<point> gb_real_coord;
};

class space_matrix_assembler
{
    space_scheme_config ssc;

    spmat permutation;
    list<triplet> base_boundary_system;
    vector<subdomain_system> subdomain_systems;



    void create_from_grid();
    void set_scheme();

    void build_permutation();
    void build_subdomain_systems();
    void build_subdomain_system(int i, const subdomain_interaction &);
    void build_base_boundary_system();
    void build_spmat_from_interaction(const interaction &, spmat &);


    
    public:

    space_matrix_assembler(case_config cc_);
    const spmat & get_permutation() const {return permutation;}
    const list<triplet> & get_base_boundary_system() const {return base_boundary_system;}
    const vector<subdomain_system> & get_subdomain_systems() const {return subdomain_systems;}

    bool import_matrixes(std::istringstream &);
    bool export_matrixes(std::ostringstream &);
};

#endif