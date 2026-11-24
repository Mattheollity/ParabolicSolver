#include "geometry_io.h"

//template<typename DATATYPE>
//void build_vector_from_strings(vector<DATATYPE> &, std::ifstream & file);

void geometry_io::import_geometry(geometry & geometry_obj, string path)
{
    import_domain(geometry_obj, path);
    import_subdomains(geometry_obj, path);
}

void geometry_io::import_domain(geometry & geometry_obj, string path)
{
    string vertexes_path = path+"/vertexes";
    auto file = load_file(vertexes_path);

    build_vector_from_strings("vertexes", geometry_obj.vertexes, file);
    for(auto v : geometry_obj.vertexes)
    {
        if(v.subdomain==-1) geometry_obj.boundary.push_back(v.index);
        else geometry_obj.inner.push_back(v.index);
    }

    string permutation_path = path+"/vertexes";
    auto file2 = load_file(permutation_path);

    build_vector("permutation", geometry_obj.permutation, file);
}

void geometry_io::import_subdomains(geometry & geometry_obj, string path)
{
    auto & subdomains = geometry_obj.subdomains;
    for(int i=0; i<subdomains.size(); i++)
    {
        auto gc = geometry_obj.gc;

        int ibuf;
        int ibuf2;
        double dbuf;
        string buf;

        subdomain & s = subdomains[i];
        string subdomain_path = path+"/"+to_string(i);
        auto file = load_file(subdomain_path);

        load_pair_of_vectors<voxel_point>("vbbox", gc.dim, file);
        build_vector<int>("inner_indexes", geometry_obj.inner, file);
        build_vector<int>("boundary_indexes", geometry_obj.boundary, file);
    }
}

// void geometry_io::export_grid() const
// {
//     fs::copy_file(GRID_DIR+"/grid_conds", GRID_DIR+"/exp_grid_conds", fs::copy_options::overwrite_existing);
//     export_domain();
//     export_subdomains();
// }

// void geometry_io::export_domain() const
// {
//     std::ofstream file(GRID_DIR+"/vertexes");
//     file<<"amount "<<vertexes.size()<<endl;
//     for(vertex v : vertexes) file<<v.string_export()<<endl;
//     file.close();

//     std::ofstream file2(GRID_DIR+"/permutation");
//     for(int i : permutation) file2<<i<<" ";
//     file2.close();
// }

// void geometry_io::export_subdomains() const
// {
//     for(subdomain s : subdomains)
//     {
//         std::ofstream file(GRID_DIR+"/"+to_string(s.index));

//         file<<"vbbox"<<endl;
//         for(auto i : s.N_bbox.first) file<<i<<" ";
//         file<<endl;
//         for(auto i : s.N_bbox.second) file<<i<<" ";
//         file<<endl;

//         file<<"neighboring_subdomains"<<endl;
//         file<<"amount "<<s.neighboring_subdomains.size()<<endl;
//         for(int i : s.neighboring_subdomains) file<<i<<" ";
//         file<<endl;

//         file<<"inner_indexes"<<endl;
//         file<<"amount "<<s.inner.size()<<endl;
//         for(auto i : s.inner) file<<i<<endl;
//         file<<"boundary_indexes"<<endl;
//         file<<"amount "<<s.boundary.size()<<endl;
//         for(auto i : s.boundary) file<<i<<endl;

//         file.close();
//     }
// }