#include "grid.h"
#include "file_io.h"
#include "problem_conditions.h"
#include <cstdio>
#include <omp.h>
#include <vector>

grid::grid(case_config cc_) :
    GRID_DIR(cc.GRID_DIR)
{
    setup_grid_conds();
    bool succed_import = import_grid();
    if(!succed_import) create_grid();
}





