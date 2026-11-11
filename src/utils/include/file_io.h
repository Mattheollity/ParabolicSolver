#ifndef FILE_IO_
#define FILE_IO_

#include "common.h"
#include "problem_conditions.h"

namespace fs = std::filesystem;

struct case_config
{
    string GRID_DIR;
    string SPACE_DIR;
    string TIME_DIR;
};

grid_conds load_grid_conds(string grid_file);

time_integr_conds load_tic(string tic_file);
approx load_approx_conds(string appr_file);

bool file_equality(string FILE1, string FILE2);

void write_matrix(const spmat & m, string s);
void write_matrix(const mat & m, string s);
void write_matrix_struct(const spmat & m, string s);
void write_matrix_struct(mat, string);
void write_matrix(vect, string);

void import_spmat(string FILENAME, spmat &);
void export_spmat(string FILENAME, const spmat &);

#endif
