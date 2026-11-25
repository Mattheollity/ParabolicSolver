#include "space_scheme_io.h"


void build_anisotrop_laplace(
    unique_ptr<anisotrop_laplace> & ss,
    std::ifstream & file
)
{
    expected_section("anisotrop_laplace_config", file);

    int dim = load_int("dim", file);

    ss = std::make_unique<anisotrop_laplace>();

    ss->k = load_vector<double>(
        "anisotrop_coeffs", dim, file
    );
}



void write_anisotrop_laplace(
    const anisotrop_laplace & ss,
    std::ofstream & file
)
{
    file<<"anisotrop_laplace_config\n";
    file<<"dim "<<ss.k.size()<<"\n";

    file<<"anisotrop_coeffs ";

    for(double k : ss.k)
    file<<k<<" ";

    file<<"\n";
}