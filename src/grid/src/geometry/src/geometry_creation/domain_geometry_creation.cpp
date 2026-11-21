#include "geometry.h"

void geometry::create_geometry()
{
    fill_domain();
    setup_subdomains();
    build_permutation_vector();
}

void geometry::fill_domain()
{
    double t0, t1;

    t0 = omp_get_wtime();
    cout<<"Start fill_domain"<<endl;

    int i=0;
    for(dim_iterator di({voxel_point(gc.dim,0), gc.N}); !di.is_overfull(); ++di)
    {
        vertexes[i] = {i, di.voxel_coord, to_real_coord(di.voxel_coord)};
        find_neighbors_and_axis(vertexes[i]);
        i++;
    }

    t1 = omp_get_wtime();
    cout<<"End fill_domain with time: "<<t1-t0<<endl;
}

void geometry::find_neighbors_and_axis(vertex & v)
{
    int inner_cond = gc.dim*2;

    auto & vc = v.voxel_coord;
    vector<int> neighbors;
    vector<int> neighbor_axis;

    for(int d=0;d<gc.dim;d++)
    if(vc[d]+1<gc.N[d])
    {
        auto foo = vc;
        foo[d]+=1;
        neighbors.push_back(global_index(foo));
        neighbor_axis.push_back(d+1);
    }
    for(int d=0;d<gc.dim;d++)
    if (vc[d]-1>=0)
    {
        auto foo = vc;
        foo[d]-=1;
        neighbors.push_back(global_index(foo));
        neighbor_axis.push_back(-d-1);
    }

    v.neighbors_indexes=neighbors;
    v.neighbors_dims=neighbor_axis;
    if(neighbors.size()==inner_cond) inner.push_back(v.index);
    else boundary.push_back(v.index);
}

int geometry::global_index(const voxel_point & voxel_coord)
{
    int index=0;
    for(int d=0;d<gc.dim;d++)
    index+=csum[d]*voxel_coord[d];

    return index;
}