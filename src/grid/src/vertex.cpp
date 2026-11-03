#include "vertex.h"

vertex::vertex(int index_, voxel_point vp) :
    index(index_),
    dim(vp.size()),
    voxel_coord(vp),
    subdomain(-1)
{};

vertex & vertex::string_import(string str)
{
    int neigbors_size;

    std::istringstream ss(str);

    ss>>index;
    ss>>subdomain;

    for(int i=0;i<dim;i++) ss>>voxel_coord[i];

    ss>>neigbors_size;
    neighbors_indexes.resize(neigbors_size);
    for(int i=0;i<neigbors_size;i++) ss>>neighbors_indexes[i];
    neighbors_dims.resize(neigbors_size);
    for(int i=0;i<neigbors_size;i++) ss>>neighbors_dims[i];

    return *this;
}

string vertex::string_export() const
{
    int neigbors_size = neighbors_indexes.size();

    std::ostringstream ss;

    ss<<index<<" ";
    ss<<subdomain<<" ";

    for(int i=0;i<dim;i++) ss<<voxel_coord[i]<<" ";
    
    ss<<neigbors_size<<" ";
    for(int i=0;i<neigbors_size;i++) ss<<neighbors_indexes[i]<<" ";
    for(int i=0;i<neigbors_size;i++) ss<<neighbors_dims[i]<<" ";

    return ss.str();
}