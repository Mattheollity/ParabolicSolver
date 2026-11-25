#ifndef SPACE_SCHEME_
#define SPACE_SCHEME_

#include "common.h"
#include "interactions.h"

struct space_scheme
{
    protected:

    point compute_abs_h(const point & base, const point & sub) const;

    public:

    virtual ~space_scheme() = default;

    virtual vector<double> operator()(const vertex_interaction &) const = 0;
};



struct anisotrop_laplace : space_scheme
{
    vector<double> k;

    vector<double> operator()(const vertex_interaction &) const override;
};

#endif