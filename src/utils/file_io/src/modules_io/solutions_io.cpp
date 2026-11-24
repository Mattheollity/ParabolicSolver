#include "solutions_io.h"

void build_const_dirichlet(unique_ptr<const_dirichlet> & cd, std::ifstream & file)
{
    expected_section("const_dirichlet_settings", file);
    cd->C = load_field<double>("C", file);
}