#ifndef TIME_SCHEME_
#define TIME_SCHEME_

#include "common.h"
#include "space_matrix_assembler.h"
#include "eigen_parts.h"

struct time_scheme
{
    virtual void build_next_time_step_system(const subdomain_system &, subdomain_system &) const = 0;
    virtual void build_current_time_step_system(const subdomain_system &, subdomain_system &) const = 0;
};

struct weighted_time_scheme : time_scheme
{
    double theta;
    double T;

    virtual void build_next_time_step_system(const subdomain_system &, subdomain_system &) const override;
    virtual void build_current_time_step_system(const subdomain_system &, subdomain_system &) const override;

    weighted_time_scheme(double theta_, double T_) : theta(theta_), T(T_) {};
};



#endif