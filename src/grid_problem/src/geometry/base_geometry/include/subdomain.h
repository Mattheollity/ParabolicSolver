#ifndef SUBDOMAIN_
#define SUBDOMAIN_

#include "common.h"

struct subdomain
{
    int index = -1;

    vector<int> inner;
    vector<int> boundary;
    vector<int> global_boundary;
};

struct subdomain_coordinates
{
    vector<const point *> inner;
    vector<const point *> boundary;
    vector<const point *> global_boundary;
};

#endif