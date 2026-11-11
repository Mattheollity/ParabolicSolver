#ifndef MATRIX_ASSEMBLER_
#define MATRIX_ASSEMBLER_

#include "file_io.h"
#include "grid.h"
#include "common.h"
#include "problem_conditions.h"
#include "space_scheme.h"
#include <sstream>

typedef vector<vector<spmat>> block_matrix;

struct subdomain_system
{
    spmat ii;
    spmat ib;
    spmat bb;
    spmat global_boundary;
};

class matrix_assembler
{
    case_config cc;
    grid * g;
    space_scheme scheme;

    spmat permutation;
    list<triplet> base_boundary_system;
    vector<subdomain_system> subdomain_systems;



    void create_from_grid();
    void set_scheme();

    void build_permutation();
    void build_base_boundary_system();
    void build_subdomain_systems();

    void build_spmat_from_interaction(const interaction &, spmat &);

    bool ability_to_import();
    void import_matrixes();
    void export_matrixes();


    
    public:

    matrix_assembler(case_config cc_);
    const spmat & get_permutation();
    const list<triplet> & get_base_boundary_system();
    const subdomain_system & get_subdomain_system(int subdomain);

};

struct space_system
{
    block_matrix system;
    
    spmat right_addition;
    spmat permutation;
    int amount_of_parts;
    vector<int> sizes;
    vector<point> boundary_points;
    approx appr;
    vect boundary_at_time(double t) const;
};

class poisson_system
{
    sparse_system p_system;

    void build_main_diagonal_system(const vector<interaction> &, int i);
    void build_interaction_system(const vector<interaction> &, int i, int j);
    void fill_spmat(spmat &, const interaction &);

    void write_spmat(std::ofstream &, spmat &);
    void load_spmat(std::ifstream &, spmat &);

    void build_right_addition(const grid & g);
    void add_zeros_matrixes();
    void build_permutation_matrix(const grid & g);

    std::stringstream poissonlog;

    public:

    poisson_system(const grid &, const approx &);
    poisson_system(string s) {load_from_file(s);}

    const sparse_system & build_system(const grid &);
    const sparse_system & get_system() const {return p_system;}
    const vector<int> & get_system_sizes() const;
    void write_to_file(string);
    void load_from_file(string);

    vect get_diriclet_boundary(const grid & g);

    string get_poisson_log() const {return poissonlog.str();}
    
};

#endif