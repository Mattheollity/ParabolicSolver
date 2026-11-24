#include "vertex.h"

#include <iomanip>
#include <limits>
#include <stdexcept>

vertex::vertex(int index_, point p) :
    index(index_),
    real_coord(std::move(p))
{}



vertex & vertex::string_import(string str)
{
    std::istringstream ss(str);

    int dim;
    int neighbors_size;

    vertex v;

    if(!(ss>>v.index>>v.subdomain>>dim) || dim<=0)
    throw std::runtime_error("Invalid vertex data");

    v.real_coord.resize(dim);

    for(double & x : v.real_coord)
    if(!(ss>>x)) throw std::runtime_error("Invalid vertex coordinates");

    if(!(ss>>neighbors_size) || neighbors_size<0) throw std::runtime_error("Invalid neighbors size");

    v.neighbors_indexes.resize(neighbors_size);

    for(int & i : v.neighbors_indexes)
    if(!(ss>>i)) throw std::runtime_error("Invalid neighbor index");

    *this = std::move(v);

    return *this;
}



string vertex::string_export() const
{
    std::ostringstream ss;

    ss<<std::setprecision(std::numeric_limits<double>::max_digits10);

    ss<<index<<" ";
    ss<<subdomain<<" ";
    ss<<real_coord.size()<<" ";

    for(double x : real_coord)
    ss<<x<<" ";

    ss<<neighbors_indexes.size()<<" ";

    for(int i : neighbors_indexes)
    ss<<i<<" ";

    return ss.str();
}