#include "space_scheme_io.h"

void build_anisotrop_laplace(unique_ptr<anisotrop_laplace> & ss, std::ifstream & file)
{
    expected_section("anisotrop_laplace_config", file);
    int dim = load_int("dim", file);
    ss->k = load_vector<double>("anisotrop_coeffs", dim, file);
}