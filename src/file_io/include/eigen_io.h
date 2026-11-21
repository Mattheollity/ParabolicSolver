#ifndef EIGEN_IO_
#define EIGEN_IO_

#include "common.h"

bool file_equality(string FILE1, string FILE2);

void write_matrix(const spmat & m, string s);
void write_matrix(const mat & m, string s);
void write_matrix(vect, string);

void write_matrix_struct(const spmat & m, string s);
void write_matrix_struct(mat, string);

template<typename ARRAY_TRIPLETS>
void import_spmat(string FILENAME, ARRAY_TRIPLETS trp);
template<typename ARRAY_TRIPLETS>
void export_spmat(string FILENAME, ARRAY_TRIPLETS trp);

void import_spmat(string FILENAME, spmat &);
void export_spmat(string FILENAME, const spmat &);

#endif
