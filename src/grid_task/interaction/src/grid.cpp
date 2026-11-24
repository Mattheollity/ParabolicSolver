#include "grid.h"

grid::grid(const grid_config & gc_) :
    gc(gc_)
{
    if(gc.pass_creationg) return;
    grid_geometry = std::make_unique<geometry>(gc);
    build_self_subdomain_interactions();
    build_boundaries_interaction();
}