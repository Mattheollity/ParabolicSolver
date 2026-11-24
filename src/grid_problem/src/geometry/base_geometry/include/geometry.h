#ifndef GEOMETRY_
#define GEOMETRY_

#include "common.h"

#include "vertex.h"
#include "subdomain.h"

class geometry_io;

class geometry
{
    friend geometry_io;

    protected:

    vector<vertex> vertexes;
    vector<int> inner;
    vector<int> boundary;
    vector<int> permutation;

    vector<subdomain> subdomains;

    geometry() = default;

    void build_permutation_vector();

    public:

    virtual ~geometry() = default;

    virtual void setup_geometry() = 0;

    int find_local_index(int i, const vector<int> &) const;
    point compute_abs_h(const point & a, const point & b) const;

    const vector<int> & get_permutation() const {return permutation;}
};

#endif