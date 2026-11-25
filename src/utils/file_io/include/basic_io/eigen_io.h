#ifndef EIGEN_IO_
#define EIGEN_IO_

#include "common.h"
#include "eigen_parts.h"


bool file_equality(string FILE1, string FILE2);



void write_matrix(const spmat & m, string s);

void write_matrix(const mat & m, string s);

void write_matrix(const vect & v, string s);



void write_matrix_struct(const spmat & m, string s);

void write_matrix_struct(const mat & m, string s);



template<typename ARRAY_TRIPLETS>
void import_spmat(string FILENAME, ARRAY_TRIPLETS & trp);

template<typename ARRAY_TRIPLETS>
void export_spmat(string FILENAME, const ARRAY_TRIPLETS & trp);



void import_spmat(string FILENAME, spmat & m);

void export_spmat(string FILENAME, const spmat & m);



#include "eigen_io.tpp"

#endif
