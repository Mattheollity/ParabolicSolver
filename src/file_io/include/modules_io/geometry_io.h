#ifndef GEOMETRY_IO_
#define GEOMTRY_IO_

#include "common.h"
#include "basic_io.h"
#include "geometry.h"

class geometry_io
{
    void import_domain(geometry &, string path);
    void import_subdomains(geometry &, string path);
    
    void export_domain(const geometry &, string path) const;
    void export_subdomains(const geometry &, string path) const;

    public:

    void import_geometry(geometry &, string path);
    void export_geometry(const geometry &, string path) const;
};

#endif