#include "geometry.h"

bool geometry::import_grid()
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

bool geometry::import_domain()
{
    string buf;
    int ibuf;

    std::ifstream file(GRID_DIR+"/vertexes");

    file>>buf;
    if(buf=="amount")
    {
        file>>ibuf;
        vertexes = vector<vertex> (ibuf);
    }
    else throw simple_error("Expected amount at vertex file");

    for(int i=0;i<vertexes.size();i++)
    {
        getline(file, buf);
        vertexes[i].string_import(buf);

        if(vertexes[i].subdomain==-1) boundary.push_back(i);
        else inner.push_back(i);
    }   
    file.close();

    std::ifstream file2(GRID_DIR+"/permutation");

    permutation = vector<int> (inner.size());

    for(int i=0;i<permutation.size();i++)
    {
        file2>>ibuf;
        permutation[i]=ibuf;
    }

    file2.close();
    return true;
}

bool geometry::import_subdomains()
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

        s.N_bbox = {voxel_point(gc.dim), voxel_point(gc.dim)};
        for(int d=0;d<gc.dim;d++)
        {
            file>>ibuf;
            s.N_bbox.first[d]=ibuf;
        }
        for(int d=0;d<gc.dim;d++)
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

void geometry::export_grid() const
{
    fs::copy_file(GRID_DIR+"/grid_conds", GRID_DIR+"/exp_grid_conds", fs::copy_options::overwrite_existing);
    export_domain();
    export_subdomains();
}

void geometry::export_domain() const
{
    std::ofstream file(GRID_DIR+"/vertexes");
    file<<"amount "<<vertexes.size()<<endl;
    for(vertex v : vertexes) file<<v.string_export()<<endl;
    file.close();

    std::ofstream file2(GRID_DIR+"/permutation");
    for(int i : permutation) file2<<i<<" ";
    file2.close();
}

void geometry::export_subdomains() const
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
        file<<"amount "<<s.inner.size()<<endl;
        for(auto i : s.inner) file<<i<<endl;
        file<<"boundary_indexes"<<endl;
        file<<"amount "<<s.boundary.size()<<endl;
        for(auto i : s.boundary) file<<i<<endl;

        file.close();
    }
}