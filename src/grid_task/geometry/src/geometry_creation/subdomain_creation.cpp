#include "geometry.h"

void geometry::setup_subdomains()
{
    for(int i=0; i<gc.subdomains.size(); i++)
    create_subdomain(i, gc.subdomains[i]);

    find_inner_and_boundary_at_subdomains();
    assert_all_vertexes_has_subdomain();
}

void geometry::create_subdomain(int i, pair<point, point> bbox_)
{
    double t0 = omp_get_wtime();
    cout<<"Creating subdomains..."<<endl;

    subdomains[i] = subdomain();
    subdomains[i].index = i;
    setup_subdomain_N_bbox(subdomains[i], bbox_);
    mark_subdomains_vertexes(subdomains[i]);

    cout<<"End creating subdomains with time: "<<omp_get_wtime()-t0<<endl<<endl;
}

void geometry::setup_subdomain_N_bbox(subdomain & s, pair<point, point> & bbox_)
{
    point & a = bbox_.first;
    point & b = bbox_.second;

    voxel_point va(gc.dim);
    voxel_point vb(gc.dim);

    for(unsigned d=0;d<gc.dim;d++)
    {
        if(a[d]<gc.bbox.first[d]) a[d]=gc.bbox.first[d];
        if(b[d]>gc.bbox.second[d]) b[d]=gc.bbox.second[d];

        double eps = h[d]/2;

        int A = (a[d]-gc.bbox.first[d])/h[d];
        if( a[d] - A*h[d] < eps) va[d] = A;
        else va[d] = A+1;

        int B = (b[d]-gc.bbox.first[d])/h[d];
        if( b[d] - B*h[d] < eps) vb[d] = B;
        else vb[d] = B+1;

        if(va[d]==0) va[d]+=1;
        if(vb[d]==gc.N[d]) vb[d]-=1;
    }

    s.N_bbox.first=va;
    s.N_bbox.second=vb;
}

void geometry::mark_subdomains_vertexes(subdomain & s)
{
    for(dim_iterator di(s.N_bbox); !di.is_overfull(); ++di)
    {
        vertex & v = vertexes[global_index(di.voxel_coord)];

        if(v.subdomain!=-1) throw simple_error("SUBDOMAINS CAN`T HAVE COMMON VERTEXES");
        v.subdomain = s.index;
    } 
}

void geometry::find_inner_and_boundary_at_subdomains()
{
    int inner_vertex_cond = gc.dim*2;

    for(auto v : vertexes)
    if(v.subdomain!=-1)
    {
        subdomain & subd = subdomains[v.subdomain];

        if(inner_vertex_cond == amount_neighbors_at_subdomain(v))
        subd.inner.push_back(v.index);
        else
        subd.boundary.push_back(v.index);
    }
    else
    {
        for(int ni : v.neighbors_indexes)
        if(vertexes[ni].subdomain!=-1)
        {
            int s = vertexes[ni].subdomain;
            subdomains[s].global_boundary.push_back(v.index);
        }
    }
}

int geometry::amount_neighbors_at_subdomain(const vertex & v) const
{
    int amount = 0;
    for(auto n : v.neighbors_indexes)
    {
        auto nv = vertexes[n];
        if(nv.subdomain==v.subdomain) amount++;
    }
    return amount;
}

void geometry::assert_all_vertexes_has_subdomain() const
{
    int inner_vertex_cond = 2*gc.dim;
    #pragma omp parallel for firstprivate(vertexes, subdomains, inner_vertex_cond)
    for(int i=0; i<vertexes.size(); i++)
    if(
        vertexes[i].neighbors_indexes.size()==inner_vertex_cond
        &&
        vertexes[i].subdomain==-1
    ) throw("Any vertex must has subdomain or be boundary vertex");
}