#ifndef TIME_SCHEME_
#define TIME_SCHEME_

#include "common.h"
#include "config_structures_ss.h"
#include "grid_problem.h"


struct time_scheme
{
    bool static_time_step = true;

    virtual ~time_scheme() = default;

    virtual double compute_local_time_step(
        const subdomain_system & ss,
        const subdomain_values & current
    ) const = 0;

    virtual void build_next_time_step_system(
        const subdomain_system & ss,
        subdomain_system & tss,
        double T
    ) const = 0;

    virtual void build_current_time_step_system(
        const subdomain_system & ss,
        subdomain_system & tss,
        double T
    ) const = 0;

    virtual void build_boundary_time_step_systems(
        const spmat & ss,
        spmat & current_time,
        spmat & next_time,
        double T
    ) const = 0;
};



struct weighted_time_scheme : time_scheme
{
    double theta;
    double T;

    weighted_time_scheme(double theta_, double T_) :
    theta(theta_), T(T_) {};

    double compute_local_time_step(
        const subdomain_system & ss,
        const subdomain_values & current
    ) const override;

    void build_next_time_step_system(
        const subdomain_system & ss,
        subdomain_system & tss,
        double T
    ) const override;

    void build_current_time_step_system(
        const subdomain_system & ss,
        subdomain_system & tss,
        double T
    ) const override;

    void build_boundary_time_step_systems(
        const spmat & ss,
        spmat & current_time,
        spmat & next_time,
        double T
    ) const override;
};

#endif