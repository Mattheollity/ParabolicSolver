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

void grid::setup_grid_conds()
{
    gc = load_grid_conds(GRID_DIR+"/grid_conds");
    h.resize(gc.dim);
    subdomains.resize(gc.subdomains.size());

    csum = vector<int> (gc.dim, 1);
    for(int d=1;d<gc.dim;d++)
    csum[d]*=csum[d-1]*gc.N[d-1];

    for(int i=0;i<gc.dim;i++)
    h[i] = (gc.bbox.second[i] - gc.bbox.first[i])/(gc.N[i]-1);

    int amount = 1;
    for(int i=0;i<gc.dim;i++) amount*=gc.N[i];
    vertexes.resize(amount);
}

point grid::to_real_coord(voxel_point vp) const
{
    point p(gc.dim);
    for(int i=0; i<gc.dim; i++)
    p[i] = vp[i]*h[i]+gc.bbox.first[i];

    return p;
}

voxel_point grid::to_voxel_coord(point p) const
{
    voxel_point vp(gc.dim);
    for(int i=0; i<gc.dim; i++)
    vp[i] = (p[i]-gc.bbox.first[i])/h[i];

    return vp;
}



void grid::build_permutation_vector()
{
    double t1 = omp_get_wtime();
    permutation.resize(inner.size());
    vector<int> paddings (subdomains.size()+1);

    paddings[0] = 0;
    for(int i=0; i<subdomains.size(); i++)
    paddings[i+1] = paddings[i] + subdomains[i].boundary.size() + subdomains[i].inner.size();

    omp_set_nested(1);
    #pragma omp parallel for shared(paddings, permutation, subdomains, inner)
    for(int i=0; i<subdomains.size(); i++)
    {
        const subdomain & s = subdomains[i];

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
        for(int j=0; j<s.inner.size(); j++)
        permutation[paddings[i]+j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.inner[j]
        ) - inner.begin();

        #pragma omp parallel for shared(paddings, permutation_, subdomains, inner, i)
        for(int j=0; j<s.boundary.size(); j++)
        permutation[paddings[i] + s.inner.size() + j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.boundary[j]
        ) - inner.begin();
    }
    omp_set_nested(0);
    double t2 = omp_get_wtime();

    cout<<"Permutation`s indices was created with time: "<<t2-t1<<endl<<endl;
}
