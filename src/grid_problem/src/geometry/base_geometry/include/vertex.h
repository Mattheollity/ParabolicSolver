#ifndef VERTEX_
#define VERTEX_

#include "common.h"

struct vertex
{
    int index = -1;
    int subdomain = -1;

    point real_coord;

    vector<int> neighbors_indexes;

    vertex() = default;
    vertex(int index_, point p);
    vertex(string s) {string_import(s);}

    string string_export() const;
    vertex & string_import(string);
};

#endif