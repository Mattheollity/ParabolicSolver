#include "common.h"
#include "file_io.h"
#include "poisson.h"
#include "problem_conditions.h"
#include "approx.h"
#include "time_scheme.h"
#include <omp.h>

int main()
{
    string gc_file, tic_file, appr_file;
    int ntr;
    cout<<"Grid conds: ";
    cin>>gc_file;
    cout<<"Time integr conds: ";
    cin>>tic_file;
    cout<<"Approx conds: ";
    cin>>appr_file;
    cout<<"OMP THREADS NUM: ";
    cin>>ntr;

    omp_set_num_threads(ntr);

    grid_conds gc = load_grid_conds(gc_file);
    grid g(gc);
    time_integr_conds tic = load_tic(tic_file);
    tic.set_zero_time_data(trip_sin, g.get_inner_points());
    
    approx appr = load_approx_conds(appr_file, g);
    appr.approx_func = approx_func;
    appr.dirichlet = dirich_boundary_conds;

    poisson_system ps(g, appr);
    tsws solver(ps.get_system(), tic);
    solver.solve();

    return 0;
}