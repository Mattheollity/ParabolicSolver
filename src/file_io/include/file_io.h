#ifndef FILE_IO_
#define FILE_IO_

#include "common.h"
#include "config_structures.h"
#include <string>

namespace fs = std::filesystem;

пк load_grid_conds(string grid_file);
time_conds load_time_conds(string TIME_DIR);
case_conds load_case_conds(string TIME_DIR);
writter_conds load_writter_conds(string TIME_DIR);


bool file_equality(string FILE1, string FILE2);

void write_matrix(const spmat & m, string s);
void write_matrix(const mat & m, string s);
void write_matrix_struct(const spmat & m, string s);
void write_matrix_struct(mat, string);
void write_matrix(vect, string);

template<typename ARRAY_TRIPLETS>
void import_spmat(string FILENAME, ARRAY_TRIPLETS trp);
template<typename ARRAY_TRIPLETS>
void export_spmat(string FILENAME, ARRAY_TRIPLETS trp);

void import_spmat(string FILENAME, spmat &);
void export_spmat(string FILENAME, const spmat &);

#endif
