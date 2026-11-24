#include "geometry.h"

geometry::geometry(const grid_config & gc_) :
    gc(gc_)
{
    setup_geometry();
    if(!gc.try_to_import) create_geometry();
}

void geometry::setup_geometry()
{
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

void geometry::build_permutation_vector()
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

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
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



point geometry::to_real_coord(voxel_point vp) const
{
    point p(gc.dim);
    for(int i=0; i<gc.dim; i++)
    p[i] = vp[i]*h[i]+gc.bbox.first[i];

    return p;
}

vector<point> geometry::to_real_coord(vector<voxel_point> vec_vp) const
{
    vector<point> vec_p(vec_vp.size());
    for(int i=0;i<vec_p.size();i++)
    vec_p[i]=to_real_coord(vec_vp[i]);

    return vec_p;
}

voxel_point geometry::to_voxel_coord(point p) const
{
    voxel_point vp(gc.dim);
    for(int i=0; i<gc.dim; i++)
    vp[i] = (p[i]-gc.bbox.first[i])/h[i];

    return vp;
}

int geometry::find_local_index(int i, const vector<int> & local) const
{
    return std::lower_bound
        (
            local.begin(),
            local.end(),
            i
        ) - local.begin();
}

point geometry::compute_abs_h(const point & a, const point & b) const
{
    point h(a);

    for(int d=0;d<gc.dim;d++)
    h[d]=std::abs(a[d]-b[d]);

    return h;
}

vector<vector<point>> geometry::get_global_boundaries_points() const
{
    vector<vector<point>> global_boundary_points(subdomains.size());
    for(int i=0;i<global_boundary_points.size();i++)
    {
        auto gb_vector_coords = subdomains[i].global_boundary;
        vector<point> gb_vector_points(gb_vector_coords.size());
        for(int j=0;j<gb_vector_points.size();j++)
        gb_vector_points[j] = vertexes[gb_vector_coords[j]].real_coord;
        global_boundary_points[i] = gb_vector_points;
    }

    return global_boundary_points;
}