#ifndef CONFIGS_IO_
#define CONFIGS_IO_

#include "common.h"
#include "config_structures.h"
#include "solutions_io.h"
#include "space_scheme_io.h"
#include "config_structures_ts.h"


struct case_config
{
    space_scheme_config ssc;
    time_scheme_config tsc;
    initial_boundary_config ibc;
    grid_config gc;
    computation_config cc;
    writter_config wc;
};



case_config load_case(string CASE_DIR);

space_scheme_config load_space_scheme(string CASE_DIR);

time_scheme_config load_time_scheme(string CASE_DIR);

grid_config load_grid(string CASE_DIR);

initial_boundary_config load_initial_boundary(string CASE_DIR);

computation_config load_computation(string CASE_DIR);

writter_config load_writter(string CASE_DIR);



void write_case(string CASE_DIR, const case_config & cfg);

void write_space_scheme(string CASE_DIR, const space_scheme_config & ssc);

void write_time_scheme(string CASE_DIR, const time_scheme_config & tsc);

void write_grid(string CASE_DIR, const grid_config & gc);

void write_initial_boundary(string CASE_DIR, const initial_boundary_config & ibc);

void write_computation(string CASE_DIR, const computation_config & cc);

void write_writter(string CASE_DIR, const writter_config & wc);

#endif