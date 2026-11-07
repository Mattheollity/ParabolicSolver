#include "grid.h"
#include "file_io.h"
#include "problem_conditions.h"
#include <omp.h>
#include <vector>

grid::grid(string GRID_DIR_) :
    GRID_DIR(GRID_DIR_)
{
    setup_grid();
    bool succed_import = import_grid();
    if(!succed_import) create_grid();
}

void grid::setup_grid()
{
    grid_conds gc = load_grid_conds(GRID_DIR+"/grid_conds");
    bbox = gc.bbox;
    dim = gc.N.size();
    N_bbox = gc.N;
    h = vector<double> (dim);
    subdomains = vector<subdomain> (gc.subdomains.size());
}

void grid::create_grid()
{
    double t0;

    cout<<"---------------GRID CREATION--------------"<<endl;
    t0 = omp_get_wtime();

    grid_conds gc = load_grid_conds(GRID_DIR+"/grid_conds");

    for(int i=0;i<dim;i++)
    h[i] = (bbox.second[i] - bbox.first[i])/(N_bbox[i]-1);

    int amount = 1;
    for(int i=0;i<dim;i++) amount*=N_bbox[i];
    vertexes = vector<vertex>(amount);

    fill_domain();

    for(int i=0; i<gc.subdomains.size(); i++)
    create_subdomain(i, gc.subdomains[i]);

    assert_all_vertexes_has_subdomain();

    permutation = build_permutation_vector();

    export_grid();

    cout<<"GRID WAS SUCCESSFULLY CREATED with time: "<<omp_get_wtime()-t0<<endl<<endl;
}

void grid::fill_domain()
{
    double t0, t1;

    t0 = omp_get_wtime();
    cout<<"Start fill_domain"<<endl;

    int i=0;
    for(dim_iterator di({voxel_point(dim,0), N_bbox}); !di.is_overfull(); ++di)
    {
        vertexes[i] = {i, di.voxel_coord};
        find_neighbors_and_axis(vertexes[i]);
        i++;
    }

    t1 = omp_get_wtime();
    cout<<"End fill_domain with time: "<<t1-t0<<endl;
}

void grid::find_neighbors_and_axis(vertex & v)
{
    int inner_cond = dim*2;

    auto & vc = v.voxel_coord;
    vector<int> neighbors;
    vector<int> neighbor_axis;

    for(int d=0;d<dim;d++)
    if(vc[d]+1<N_bbox[d])
    {
        auto foo = vc;
        foo[d]+=1;
        neighbors.push_back(global_index(foo));
        neighbor_axis.push_back(d);
    }
    else if (v.voxel_coord[d]-1>=0)
    {
        auto foo = vc;
        foo[d]-=1;
        neighbors.push_back(global_index(foo));
        neighbor_axis.push_back(-d);
    }

    v.neighbors_indexes=neighbors;
    v.neighbors_dims=neighbor_axis;

    if(neighbors.size()==inner_cond) inner.push_back(v.index);
    else boundary.push_back(v.index);
}

int grid::global_index(const voxel_point & voxel_coord)
{
    int index=0;
    for(int d=0;d<dim;d++)
    index+=csum[d]*voxel_coord[d];

    return index;
}

void grid::create_subdomain(int i, pair<point, point> bbox_)
{
    subdomains[i] = subdomain();
    subdomains[i].index = i;
    setup_subdomain_N_bbox(subdomains[i], bbox_);
    fill_subdomain(subdomains[i]);
}

void grid::setup_subdomain_N_bbox(subdomain & s, pair<point, point> & bbox_)
{
    point & a = bbox_.first;
    point & b = bbox_.second;

    voxel_point va(dim);
    voxel_point vb(dim);

    for(unsigned d=0;d<dim;d++)
    {
        if(a[d]<bbox.first[d]) a[d]=bbox.first[d];
        if(b[d]>bbox.second[d]) b[d]=bbox.second[d];

        double eps = h[d]/2;

        int A = (a[d]-bbox.first[d])/h[d];
        if( a[d] - A*h[d] < eps) va[d] = A;
        else va[d] = A+1;

        int B = (b[d]-bbox.first[d])/h[d];
        if( b[d] - B*h[d] < eps) vb[d] = B;
        else vb[d] = B+1;

        if(va[d]==0) va[d]+=1;
        if(vb[d]==N_bbox[d]) vb[d]-=1;
    }

    s.N_bbox.first=va;
    s.N_bbox.second=vb;
}

void grid::fill_subdomain(subdomain & s)
{
    for(dim_iterator di(s.N_bbox); !di.is_overfull(); ++di)
    {
        vertex & v = vertexes[global_index(di.voxel_coord)];

        if(v.subdomain!=-1) throw simple_error("SUBDOMAINS CAN`T HAVE COMMON VERTEXES");
        v.subdomain = s.index;
    } 
}

void grid::find_inner_boundary_at_subdomains()
{
    int inner_vertex_cond = dim*2;

    for(auto v : vertexes)
    if(inner_vertex_cond == amount_neighbors_at_subdomain(v)) subdomains[v.subdomain].inner.push_back(v.index);
    else subdomains[v.subdomain].boundary.push_back(v.index);
}

int grid::amount_neighbors_at_subdomain(const vertex & v) const
{
    int amount = 0;
    for(auto n : v.neighbors_indexes)
    {
        auto nv = vertexes[n];
        if(nv.subdomain==v.subdomain) amount++;
    }
    return amount;
}

void grid::assert_all_vertexes_has_subdomain() const
{
    int inner_vertex_cond = 2*dim;
    #pragma omp parallel for firstprivate(vertexes, subdomains, inner_vertex_cond)
    for(int i=0; i<vertexes.size(); i++)
    if(
        vertexes[i].neighbors_indexes.size()==inner_vertex_cond
        &&
        vertexes[i].subdomain==-1
    ) throw("Any vertex must has subdomain or be boundary vertex");
}

bool grid::import_grid()
{
    if(fs::exists(GRID_DIR+"/exp_grid_conds"))
    if(file_equality(GRID_DIR+"/grid_conds", GRID_DIR+"/exp_grid_conds"))
    {

        double t0 = omp_get_wtime();
        cout<<"Expected created grid. Try to load..."<<endl;

        bool res1 = import_domain();
        bool res2 = import_subdomains();

        cout<<"Successfully loaded with time: "<<omp_get_wtime()-t0<<endl<<endl;

        return res1 && res2;
    }

    return false;
}

bool grid::import_domain()
{
    string buf;
    int ibuf;

    std::ifstream file(GRID_DIR+"/vertexes");

    file>>buf;
    if(buf=="amount")
    {
        file>>ibuf;
        vertexes = vector<vertex> (ibuf);
        permutation = vector<int> (ibuf);
    }
    else throw simple_error("Expected amount at vertex file");

    for(int i=0;i<vertexes.size();i++)
    {
        file>>buf;
        vertexes[i].string_import(buf);
    }   
    file.close();

    std::ifstream file2(GRID_DIR+"/permutation");

    for(int i=0;i<permutation.size();i++)
    {
        file2>>ibuf;
        permutation[i]=ibuf;
    }

    file2.close();
    return true;
}

bool grid::import_subdomains()
{
    for(int i=0; i<subdomains.size(); i++)
    {
        int ibuf;
        int ibuf2;
        double dbuf;
        string buf;

        subdomain & s = subdomains[i];
        std::ifstream file(GRID_DIR+"/"+to_string(i));

        file>>buf;
        if(buf!="vbbox") throw simple_error("Expected vbbox at "+to_string(i)+" subdomain file");

        s.N_bbox = {voxel_point(dim), voxel_point(dim)};
        for(int d=0;d<dim;d++)
        {
            file>>ibuf;
            s.N_bbox.first[d]=ibuf;
        }
        for(int d=0;d<dim;d++)
        {
            file>>ibuf;
            s.N_bbox.second[d]=ibuf;
        }

        file>>buf;
        if(buf!="neighboring_subdomains") throw simple_error("Expected neighboring_subdomains at "+to_string(i)+" subdomain file");
        
        file>>buf;
        if(buf!="amount") throw simple_error("Expected neighboring_subdomains amount at "+to_string(i)+" subdomain file");
        int foo;
        file>>foo;
        for(int j=0;j<foo;j++)
        {
            file>>ibuf;
            s.neighboring_subdomains.insert(ibuf);
        }

        file>>buf;
        if(buf!="inner_indexes") throw simple_error("Expected inner_indexes at "+to_string(i)+" subdomain file");
        file>>buf;
        if(buf!="amount") throw simple_error("Expected inner_indexes amount at "+to_string(i)+" subdomain file");
        file>>ibuf;
        s.inner = vector<int> (ibuf);
        for(int j=0;j<ibuf;j++) file>>s.inner[j];

        file>>buf;
        if(buf!="boundary_indexes") throw simple_error("Expected boundary_indexes at "+to_string(i)+" subdomain file");
        file>>buf;
        if(buf!="amount") throw simple_error("Expected boundary_indexes amount at "+to_string(i)+" subdomain file");
        file>>ibuf;
        s.boundary = vector<int> (ibuf);
        for(int j=0;j<ibuf;j++) file>>s.boundary[j];

        file.close();
    }

    return true;
}

void grid::export_grid() const
{
    fs::copy_file(GRID_DIR+"/grid_conds", GRID_DIR+"/exp_grid_conds", fs::copy_options::overwrite_existing);
    export_domain();
    export_subdomains();
}

void grid::export_domain() const
{
    std::ofstream file(GRID_DIR+"/vertexes");
    file<<"amount "<<vertexes.size();
    for(vertex v : vertexes) file<<v.string_export()<<endl;
    file.close();

    std::ofstream file2(GRID_DIR+"/permutation");
    for(int i : permutation) file2<<i<<" ";
    file2.close();
}

void grid::export_subdomains() const
{
    for(subdomain s : subdomains)
    {
        std::ofstream file(GRID_DIR+"/"+to_string(s.index));

        file<<"vbbox"<<endl;
        for(auto i : s.N_bbox.first) file<<i<<" ";
        file<<endl;
        for(auto i : s.N_bbox.second) file<<i<<" ";
        file<<endl;

        file<<"neighboring_subdomains"<<endl;
        file<<"amount "<<s.neighboring_subdomains.size()<<endl;
        for(int i : s.neighboring_subdomains) file<<i<<" ";
        file<<endl;

        file<<"inner_indexes"<<endl;
        file<<"amount "<<s.inner.size();
        for(auto i : s.inner) file<<i<<endl;
        file<<"boundary_indexes"<<endl;
        file<<"amount "<<s.boundary.size();
        for(auto i : s.boundary) file<<i<<endl;

        file.close();
    }
}

point grid::to_real_coord(voxel_point vp) const
{
    point p(dim);
    for(int i=0; i<dim; i++)
    p[i] = vp[i]*h[i]+bbox.first[i];

    return p;
}

voxel_point grid::to_voxel_coord(point p) const
{
    voxel_point vp(dim);
    for(int i=0; i<dim; i++)
    vp[i] = (p[i]-bbox.first[i])/h[i];

    return vp;
}

vector<interaction> grid::get_interaction_between(int i, int j) const
{
    if( !(i<subdomains.size() && j<subdomains.size()) ) throw simple_error("index of subdomain out of range");
    if(i==j) return get_self_interaction(i);

    interaction io;

    const vector<int> & output = subdomains[i].boundary;
    const vector<int> & input = subdomains[j].boundary;
    
    interaction_request bb = {
        j,
        i,
        inner,
        output
    };

    return {create_interaction(bb)};
}

vector<interaction> grid::get_self_interaction(int i) const
{
    if( !(i<subdomains.size()) ) throw simple_error("index of subdomain out of range");

    const vector<int> & inner_ = subdomains[i].inner;
    const vector<int> & boundary_ = subdomains[i].boundary;

    interaction_request ii = {
        i,
        i,
        inner_,
        inner_
    };

    interaction_request ib = {
        i,
        i,
        inner_,
        boundary_
    };

    interaction_request bb = {
        i,
        i,
        boundary_,
        boundary_
    };


    return {
        create_interaction(ii),
        create_interaction(ib),
        create_interaction(bb)
    };
}

interaction grid::create_interaction(const interaction_request & ior) const
{
    interaction io;
    io.input_size = ior.input.size();
    io.output_size = ior.output.size();
    io.input_subdomain = ior.input_subdomain;
    io.output_subdomain = ior.output_subdomain;
    io.vertexes_io = vector<vertex_interactions>(io.output_size);

    #pragma omp parallel shared(io, ior)
    #pragma omp for
    for(int i=0; i<io.output_size; i++)
    {
        vertex_interactions_request vir = {
            i,
            ior.output[i],
            ior.input_subdomain,
            ior.input
        };
        io.vertexes_io[i] = create_vertex_interactions(vir);
    }

    return io;
}

vertex_interactions grid::create_vertex_interactions(const vertex_interactions_request & vir) const
{
    const vertex & v = vertexes[vir.vertex_index];
    vertex_interactions vio;
    vio.output_index = vir.output_index;
    vio.output_coord = v.voxel_coord;

    vector<int> input_vertex_indexes;
    for(auto nv : v.neighbors_indexes)
    if(vertexes[nv].subdomain==vir.input_subdomain)
    input_vertex_indexes.push_back(nv);

    for(int i : input_vertex_indexes)
    vio.input_indexes.push_back(
        std::find(vir.input.begin(), vir.input.end(), i)
        -
        vir.input.begin()
    );
    
    for(int i : input_vertex_indexes)
    vio.input_coords.push_back(vertexes[i].voxel_coord);

    return vio;
}

interaction grid::get_inner_boundary_interaction() const
{
    interaction_request ior = {
        0,
        -1,
        inner,
        boundary
    };

    return create_interaction(ior);
}

const vector<point> grid::get_boundary_points() const
{
    vector<point> bp(boundary.size());
    for(int i=0;i<boundary.size();i++)
    bp[i] = to_real_coord(vertexes[boundary[i]].voxel_coord);

    return bp;
}

const vector<point> grid::get_inner_points() const
{
    vector<point> bp(inner.size());
    for(int i=0;i<inner.size();i++)
    bp[i] = to_real_coord(vertexes[inner[i]].voxel_coord);

    return bp;
}

vector<int> grid::build_permutation_vector()
{
    double t1 = omp_get_wtime();
    vector<int> permutation (inner.size());
    vector<int> paddings (subdomains.size()+1);

    paddings[0] = 0;
    for(int i=0; i<subdomains.size(); i++)
    paddings[i+1] = paddings[i] + subdomains[i].boundary.size() + subdomains[i].inner.size();

    omp_set_nested(1);
    #pragma omp parallel for shared(paddings, permutation, subdomains, inner)
    for(int i=0; i<subdomains.size(); i++)
    {
        const subdomain & s = subdomains[i];

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
        for(int j=0; j<s.inner.size(); j++)
        permutation[paddings[i]+j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.inner[j]
        ) - inner.begin();

        #pragma omp parallel for shared(paddings, permutation, subdomains, inner, i)
        for(int j=0; j<s.boundary.size(); j++)
        permutation[paddings[i] + s.inner.size() + j] = std::lower_bound
        (
            inner.begin(),
            inner.end(),
            s.boundary[j]
        ) - inner.begin();
    }
    omp_set_nested(0);
    double t2 = omp_get_wtime();

    cout<<"Permutation`s indices was created with time: "<<t2-t1<<endl<<endl;
    return permutation;
}
