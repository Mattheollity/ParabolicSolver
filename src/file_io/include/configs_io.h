#ifndef CONFIGS_IO_
#define CONFIGS_IO_

#include "common.h"
#include "config_structures.h"
#include "basic_io.h"
#include "solutions.h"
#include "space_scheme.h"

struct case_config
{
    space_scheme_config ssc;
    time_scheme_config tsc;
    zero_boundary_config zbc;
    grid_config gc;
    computation_config cc;
    writter_config wc;

    case_config
    (
        space_scheme_config ssc_,
        time_scheme_config tsc_,
        zero_boundary_config zbc_,
        grid_config gc_,
        computation_config cc_,
        writter_config wc_
    ) : ssc(ssc_), tsc(tsc_), zbc(zbc_), gc(gc_), cc(cc_), wc(wc_) {};
};

case_config load_case(string CASE_DIR);
space_scheme_config load_space_scheme(string CASE_DIR);
time_scheme_config load_time_scheme(string CASE_DIR);
grid_config load_grid(string CASE_DIR);
zero_boundary_config load_zero_boundary(string CASE_DIR);
computation_config load_computation(string CASE_DIR);
writter_config load_writter(string CASE_DIR);

#endif