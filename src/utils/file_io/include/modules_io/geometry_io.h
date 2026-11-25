#ifndef GEOMETRY_IO_
#define GEOMETRY_IO_

#include "common.h"
#include "basic_io.h"
#include "geometry.h"


class geometry_io
{
    void import_domain(geometry & g, string path);
    void import_subdomains(geometry & g, string path);

    void export_domain(const geometry & g, string path) const;
    void export_subdomains(const geometry & g, string path) const;



    public:

    void import_geometry(geometry & g, string path);

    void export_geometry(const geometry & g, string path) const;
};

#endif