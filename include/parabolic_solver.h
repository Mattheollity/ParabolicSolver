#ifndef PARABOLIC_SOLVER_
#define PARABOLIC_SOLVER_

#include "common.h"

#include "configs_io.h"
#include "geometry_io.h"
#include "space_matrix_assembler_io.h"
#include "time_scheme_io.h"

#include "grid_problem.h"
#include "space_matrix_assembler.h"
#include "time_scheme_assembler.h"


class parabolic_solver
{
    string case_dir;

    case_config cfg;

    int rank = 0;
    int process_count = 1;

    unique_ptr<grid_problem> gp;
    unique_ptr<space_matrix_assembler> sma;
    unique_ptr<time_scheme_assembler> tsa;
    unique_ptr<time_scheme_io> tsio;



    bool check_geometry() const;

    bool check_space_matrixes() const;

    int find_last_time_step() const;



    void setup_grid_problem(bool import_geometry);

    void setup_space_scheme(bool import_matrixes);

    void setup_time_scheme();

    void restore_solution(int last_time_step);



    public:

    parabolic_solver(string case_dir_);

    void setup_solver();



    void import_solution(string path);

    void export_solution(string path) const;

    void restart_solver();



    void time_step();

    void solve();
};

#endif