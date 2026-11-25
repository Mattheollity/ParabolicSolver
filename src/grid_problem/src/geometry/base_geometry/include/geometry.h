#ifndef GEOMETRY_
#define GEOMETRY_

#include "common.h"

#include "vertex.h"
#include "subdomain.h"
#include "interactions.h"


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

    vector<subdomain_interaction> subdomain_intr;
    interaction boundaries_intr;

    geometry() = default;

    void build_permutation_vector();



    private:

    void build_self_subdomain_interactions();
    void build_subdomain_interaction(int i);

    void build_part_of_subdomain_interaction(
        interaction &,
        const vector<int> & input,
        const vector<int> & output
    );

    vertex_interaction create_vertex_interactions(
        int output_index,
        int vertex_index,
        const vector<int> & input
    ) const;



    void build_boundaries_interaction();
    vector<int> bvio_csum() const;
    void fill_bvio(vector<vertex_interaction> &);



    public:

    virtual ~geometry() = default;

    virtual void setup_geometry() = 0;

    void build_interactions();

    int find_local_index(int i, const vector<int> &) const;
    point compute_abs_h(const point & a, const point & b) const;

    const vector<int> & get_permutation() const {return permutation;}

    const vector<subdomain_interaction> & get_subdomain_interactions() const {return subdomain_intr;}
    const interaction & get_boundaries_interaction() const {return boundaries_intr;}
    vector<subdomain_coordinates> get_subdomain_coordinates() const;
};

#endif