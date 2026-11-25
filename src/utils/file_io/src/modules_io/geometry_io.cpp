#include "geometry_io.h"

#include <iomanip>
#include <limits>
#include <algorithm>


namespace
{
    void write_indexes(
        string FIELD,
        const vector<int> & indexes,
        std::ofstream & file
    )
    {
        file<<FIELD<<"\n";
        file<<"amount "<<indexes.size()<<"\n";

        for(int i : indexes)
        file<<i<<" ";

        file<<"\n";
    }
}



void geometry_io::import_geometry(
    geometry & g,
    string path
)
{
    import_domain(g, path);
    import_subdomains(g, path);

    vector<int> permutation = g.permutation;

    g.build_permutation_vector();

    if(g.permutation!=permutation)
    throw simple_error("Invalid geometry permutation");

    g.subdomain_intr.clear();
    g.boundaries_intr = {};

    g.build_interactions();
}



void geometry_io::import_domain(
    geometry & g,
    string path
)
{
    auto file = load_file(path+"/vertexes");

    build_vector_from_strings(
        "vertexes",
        g.vertexes,
        file
    );

    g.inner.clear();
    g.boundary.clear();

    for(int i=0; i<g.vertexes.size(); i++)
    {
        const vertex & v = g.vertexes[i];

        if(v.index!=i)
        throw simple_error("Invalid vertex index");

        for(int j : v.neighbors_indexes)
        if(j<0 || j>=g.vertexes.size())
        throw simple_error("Invalid neighbor index");

        if(v.subdomain==-1)
        g.boundary.push_back(v.index);
        else
        g.inner.push_back(v.index);
    }

    auto permutation_file = load_file(path+"/permutation");

    build_vector(
        "permutation",
        g.permutation,
        permutation_file
    );
}



void geometry_io::import_subdomains(
    geometry & g,
    string path
)
{
    auto file = load_file(path+"/subdomains");

    int amount = load_int("amount", file);

    if(amount<=0)
    throw simple_error("Invalid subdomains amount");

    g.subdomains.resize(amount);

    for(int i=0; i<amount; i++)
    {
        auto subdomain_file = load_file(path+"/"+to_string(i));

        subdomain & s = g.subdomains[i];

        s.index = load_int("index", subdomain_file);

        if(s.index!=i)
        throw simple_error("Invalid subdomain index");

        build_vector("inner", s.inner, subdomain_file);

        build_vector("boundary", s.boundary, subdomain_file);

        build_vector(
            "global_boundary",
            s.global_boundary,
            subdomain_file
        );
    }

    vector<int> ownership(g.vertexes.size(), 0);

    for(const subdomain & s : g.subdomains)
    {
        for(const auto * indexes : {&s.inner, &s.boundary})
        {
            if(!std::is_sorted(indexes->begin(), indexes->end()))
            throw simple_error("Unsorted subdomain indexes");

            for(int i : *indexes)
            {
                if(i<0 || i>=g.vertexes.size())
                throw simple_error("Invalid subdomain vertex");

                if(g.vertexes[i].subdomain!=s.index)
                throw simple_error("Invalid vertex ownership");

                ownership[i]++;
            }
        }

        if(!std::is_sorted(
            s.global_boundary.begin(),
            s.global_boundary.end()
        ))
        throw simple_error("Unsorted global boundary");

        for(int i : s.global_boundary)
        {
            if(i<0 || i>=g.vertexes.size())
            throw simple_error("Invalid boundary vertex");

            if(g.vertexes[i].subdomain!=-1)
            throw simple_error("Invalid global boundary");
        }
    }

    for(int i : g.inner)
    if(ownership[i]!=1)
    throw simple_error("Invalid subdomain decomposition");
}



void geometry_io::export_geometry(
    const geometry & g,
    string path
) const
{
    fs::create_directories(path);

    export_domain(g, path);
    export_subdomains(g, path);
}



void geometry_io::export_domain(
    const geometry & g,
    string path
) const
{
    std::ofstream file(path+"/vertexes");

    if(!file)
    throw simple_error("Cannot export vertexes");

    file<<"vertexes\n";
    file<<"amount "<<g.vertexes.size()<<"\n";

    for(const vertex & v : g.vertexes)
    file<<v.string_export()<<"\n";



    std::ofstream permutation_file(path+"/permutation");

    if(!permutation_file)
    throw simple_error("Cannot export permutation");

    write_indexes(
        "permutation",
        g.permutation,
        permutation_file
    );
}



void geometry_io::export_subdomains(
    const geometry & g,
    string path
) const
{
    std::ofstream file(path+"/subdomains");

    if(!file)
    throw simple_error("Cannot export subdomains");

    file<<"amount "<<g.subdomains.size()<<"\n";

    for(const subdomain & s : g.subdomains)
    {
        std::ofstream subdomain_file(
            path+"/"+to_string(s.index)
        );

        if(!subdomain_file)
        throw simple_error("Cannot export subdomain");

        subdomain_file<<"index "<<s.index<<"\n";

        write_indexes("inner", s.inner, subdomain_file);

        write_indexes("boundary", s.boundary, subdomain_file);

        write_indexes(
            "global_boundary",
            s.global_boundary,
            subdomain_file
        );
    }
}