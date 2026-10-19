#ifndef FILE_IO_
#define FILE_IO_

#include "common.h"
#include "poisson.h"
#include "problem_conditions.h"

namespace fs = std::filesystem;

grid_conds load_grid_conds(string grid_file);
time_integr_conds load_tic(string tic_file);
approxs load_approx_conds(string appr_file, const grid & g);

void write_matrix(const spmat & m, string s);
void write_matrix(const mat & m, string s);
void write_matrix_struct(const spmat & m, string s);
void write_matrix_struct(mat, string);
void write_matrix(vect, string);

#endif
