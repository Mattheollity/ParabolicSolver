#ifndef CONFIG_STRUCTES_
#define CONFIG_STRUCTES_

#include "common.h"
#include <memory>

struct case_config;
struct space_scheme_config;
struct space_scheme_request;
struct space_scheme;
struct time_scheme_config;
struct zero_boundary_config;
struct grid_config;
struct computation_config;
struct writter_config;

struct grid_config
{
    int dim;
    pair<point, point> bbox;
    voxel_point N;
    vector<pair<point, point>> subdomains;
    bool try_to_import;
    bool pass_creationg;
};

struct computation_config
{
    int total_memory;
    int total_time;
    int OpenMP_procs;
    int MPI_procs;

    computation_config
    (
        int total_memory_, int total_time_, 
        int OpenMP_procs_, int MPI_procs_
    ) : total_memory(total_memory_), 
        total_time(total_time_), 
        OpenMP_procs(OpenMP_procs_), 
        MPI_procs(MPI_procs_) {};
};

struct writter_config
{
    int write_interval;
    string WRITE_DIR;
    bool paraview_output;

    writter_config(int write_interval_, string WRITE_DIR_, bool paraview_output_) :
    write_interval(write_interval_), WRITE_DIR(WRITE_DIR_), paraview_output(paraview_output_) {};
};

#endif