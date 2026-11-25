#ifndef SOLUTIONS_IO_
#define SOLUTIONS_IO_

#include "basic_io.h"
#include "config_structures_sol.h"


void build_const_dirichlet(
    unique_ptr<const_dirichlet> & cd,
    std::ifstream & file
);

void write_const_dirichlet(
    const const_dirichlet & cd,
    std::ofstream & file
);

#endif