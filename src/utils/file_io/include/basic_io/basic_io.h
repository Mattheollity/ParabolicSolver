#ifndef BASIC_IO_
#define BASIC_IO_

#include "common.h"
#include <filesystem>

namespace fs = std::filesystem;


std::ifstream load_file(string FILE);

void expected_section(string FIELD, std::ifstream & file);



template<typename DATATYPE>
DATATYPE load_field(string FIELD, std::ifstream & file);

int load_int(string FIELD, std::ifstream & file);

bool load_bool(string FIELD, std::ifstream & file);

string load_string(string FIELD, std::ifstream & file);



template<typename DATATYPE>
void build_vector(
    string FIELD,
    vector<DATATYPE> & v,
    std::ifstream & file
);

template<typename DATATYPE>
void build_vector_from_strings(
    string FIELD,
    vector<DATATYPE> & v,
    std::ifstream & file
);



template<typename DATATYPE>
vector<DATATYPE> load_vector(
    string FIELD,
    int amount,
    std::ifstream & file
);

template<typename DATATYPE>
pair<vector<DATATYPE>, vector<DATATYPE>> load_pair_of_vectors(
    string FIELD,
    int amount,
    std::ifstream & file
);



#include "basic_io.tpp"

#endif