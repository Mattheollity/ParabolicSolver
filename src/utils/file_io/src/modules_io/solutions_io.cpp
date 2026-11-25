#include "solutions_io.h"


void build_const_dirichlet(
    unique_ptr<const_dirichlet> & cd,
    std::ifstream & file
)
{
    expected_section("const_dirichlet_settings", file);

    double C = load_field<double>("C", file);

    cd = std::make_unique<const_dirichlet>(C);
}



void write_const_dirichlet(
    const const_dirichlet & cd,
    std::ofstream & file
)
{
    file<<"const_dirichlet_settings\n";
    file<<"C "<<cd.C<<"\n";
}