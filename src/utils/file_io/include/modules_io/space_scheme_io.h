#ifndef SPACE_SCHEME_IO_
#define SPACE_SCHEME_IO_

#include "basic_io.h"
#include "space_scheme.h"
#include "config_structures_ss.h"

void build_anisotrop_laplace(unique_ptr<anisotrop_laplace> &, std::ifstream & file);

#endif