#include "field_reader.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using std::vector;
using std::string;

namespace
{
    std::ifstream open_file(const fs::path & path)
    {
        std::ifstream file(path);
        if(!file) throw std::runtime_error("Cannot open "+path.string());
        return file;
    }

    void expect(std::istream & file, const string & name)
    {
        string field;
        if(!(file>>field) || field!=name)
        throw std::runtime_error("Expected "+name);
    }

    template<typename T>
    T read(std::istream & file)
    {
        T value;
        if(!(file>>value)) throw std::runtime_error("Invalid file data");
        return value;
    }

    int count(std::istream & file)
    {
        int n = read<int>(file);
        if(n<0) throw std::runtime_error("Negative count");
        return n;
    }

    double real(std::istream & file)
    {
        double x = read<double>(file);
        if(!std::isfinite(x)) throw std::runtime_error("Non-finite value");
        return x;
    }

    vector<int> read_indexes(std::istream & file, const string & name)
    {
        expect(file, name);
        expect(file, "amount");
        vector<int> indexes(count(file));
        for(int & i : indexes) i = read<int>(file);
        return indexes;
    }

    void end_file(std::istream & file)
    {
        file>>std::ws;
        if(!file.eof()) throw std::runtime_error("Unexpected trailing data");
    }

    bool close(double a, double b)
    {
        return std::abs(a-b)<=1e-12*std::max({1.0, std::abs(a), std::abs(b)});
    }
}


field_reader::field_reader(const fs::path & case_dir)
{
    auto wc = open_file(case_dir/"writter");
    expect(wc, "write_interval(steps)");
    count(wc);
    expect(wc, "write_dir");
    directory = read<string>(wc);
    if(directory.is_relative()) directory = case_dir/directory;

    auto bc = open_file(case_dir/"zero_boundary");
    expect(bc, "y0");
    read<string>(bc);
    expect(bc, "boundary_type");
    if(read<string>(bc)!="const_dirichlet")
    throw std::runtime_error("Only const_dirichlet is supported");
    expect(bc, "boundary");
    read<string>(bc);
    expect(bc, "const_dirichlet_settings");
    expect(bc, "C");
    boundary_value = real(bc);

    const fs::path path = case_dir/"geometry";
    auto file = open_file(path/"vertexes");
    expect(file, "vertexes");
    expect(file, "amount");
    int n = count(file);
    if(n==0) throw std::runtime_error("Empty geometry");

    owner.resize(n);
    vtk_index.resize(n);
    vector<std::array<double, 3>> coords(n);
    std::array<vector<double>, 3> axes;
    int dim = 0;

    for(int i=0; i<n; i++)
    {
        if(read<int>(file)!=i) throw std::runtime_error("Invalid vertex index");
        owner[i] = read<int>(file);
        int d = count(file);
        if(i==0) dim = d;
        if(d!=dim || dim<1 || dim>3 || owner[i]<-1)
        throw std::runtime_error("Invalid vertex dimension or owner");

        for(int j=0; j<3; j++)
        {
            coords[i][j] = j<dim ? real(file) : 0.0;
            axes[j].push_back(coords[i][j]);
        }

        int neighbors = count(file);
        for(int j=0; j<neighbors; j++)
        {
            int index = read<int>(file);
            if(index<0 || index>=n) throw std::runtime_error("Invalid neighbor");
        }
    }
    end_file(file);

    std::size_t total = 1;
    for(int d=0; d<3; d++)
    {
        auto & a = axes[d];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        grid.dimensions[d] = a.size();
        grid.origin[d] = a.front();
        grid.spacing[d] = a.size()>1 ? (a.back()-a.front())/(a.size()-1) : 1.0;
        for(std::size_t i=0; i<a.size(); i++)
        if(!close(a[i], grid.origin[d]+i*grid.spacing[d]))
        throw std::runtime_error("VTI requires a uniform Cartesian grid");
        if(total>static_cast<std::size_t>(n)/a.size())
        throw std::runtime_error("Invalid Cartesian grid");
        total *= a.size();
    }
    if(total!=static_cast<std::size_t>(n))
    throw std::runtime_error("Incomplete Cartesian grid");

    vector<bool> occupied(n, false);
    for(int i=0; i<n; i++)
    {
        std::array<std::size_t, 3> p;
        for(int d=0; d<3; d++)
        p[d] = std::lower_bound(axes[d].begin(), axes[d].end(), coords[i][d])-axes[d].begin();
        std::size_t index = p[0]+grid.dimensions[0]*(p[1]+grid.dimensions[1]*p[2]);
        bool physical_boundary = false;
        for(int d=0; d<dim; d++)
        if(p[d]==0 || p[d]+1==grid.dimensions[d]) physical_boundary = true;
        if((owner[i]==-1)!=physical_boundary)
        throw std::runtime_error("Invalid physical boundary ownership");
        if(occupied[index]) throw std::runtime_error("Duplicate grid point");
        occupied[index] = true;
        vtk_index[i] = index;
    }

    auto sd = open_file(path/"subdomains");
    expect(sd, "amount");
    subdomains.resize(count(sd));
    end_file(sd);
    if(subdomains.empty()) throw std::runtime_error("No subdomains");

    vector<bool> owned(n, false);
    for(std::size_t i=0; i<subdomains.size(); i++)
    {
        auto f = open_file(path/std::to_string(i));
        expect(f, "index");
        if(count(f)!=static_cast<int>(i)) throw std::runtime_error("Invalid subdomain index");
        auto & s = subdomains[i];
        s.inner = read_indexes(f, "inner");
        s.boundary = read_indexes(f, "boundary");
        s.global_boundary = read_indexes(f, "global_boundary");
        end_file(f);

        for(const auto * indexes : {&s.inner, &s.boundary})
        for(int j : *indexes)
        {
            if(j<0 || j>=n || owner[j]!=static_cast<int>(i) || owned[j])
            throw std::runtime_error("Invalid owned vertex");
            owned[j] = true;
        }
        for(int j : s.global_boundary)
        if(j<0 || j>=n || owner[j]!=-1)
        throw std::runtime_error("Invalid global boundary vertex");
    }
    for(int i=0; i<n; i++)
    if(owner[i]!=-1 && !owned[i]) throw std::runtime_error("Missing owned vertex");
}


vector<int> field_reader::find_steps() const
{
    vector<int> steps;
    for(const auto & entry : fs::directory_iterator(directory))
    {
        if(!entry.is_directory()) continue;
        string name = entry.path().filename().string();
        int step;
        auto result = std::from_chars(name.data(), name.data()+name.size(), step);
        if(result.ec!=std::errc() || result.ptr!=name.data()+name.size() ||
           step<0 || std::to_string(step)!=name) continue;
        if(fs::is_regular_file(entry.path()/"solution")) steps.push_back(step);
    }
    std::sort(steps.begin(), steps.end());
    if(steps.empty()) throw std::runtime_error("No saved time layers");
    return steps;
}


field_reader::Field field_reader::readField(int step) const
{
    auto file = open_file(directory/std::to_string(step)/"solution");
    Field field;
    expect(file, "time_step");
    field.step = count(file);
    if(field.step!=step) throw std::runtime_error("Time step mismatch");
    expect(file, "time");
    field.time = real(file);
    expect(file, "T");
    field.T = real(file);
    if(field.T<0 || (step>0 && field.T==0)) throw std::runtime_error("Invalid time step");
    expect(file, "amount");
    if(count(file)!=static_cast<int>(subdomains.size()))
    throw std::runtime_error("Subdomain count mismatch");

    field.values.resize(owner.size());
    vector<bool> assigned(owner.size(), false);

    auto read_values = [&](const string & name, const vector<int> & indexes)
    {
        expect(file, name);
        if(count(file)!=static_cast<int>(indexes.size()))
        throw std::runtime_error("Invalid "+name+" vector size");
        for(int j : indexes)
        {
            double value = real(file);
            std::size_t k = vtk_index[j];
            if(assigned[k] && !close(field.values[k], value))
            throw std::runtime_error("Conflicting boundary values");
            field.values[k] = value;
            assigned[k] = true;
        }
    };

    for(std::size_t i=0; i<subdomains.size(); i++)
    {
        expect(file, "subdomain");
        if(count(file)!=static_cast<int>(i)) throw std::runtime_error("Subdomain order mismatch");
        const auto & s = subdomains[i];
        read_values("inner", s.inner);
        read_values("local_boundary", s.boundary);
        read_values("global_boundary", s.global_boundary);
    }
    end_file(file);

    for(std::size_t i=0; i<owner.size(); i++)
    if(!assigned[vtk_index[i]])
    {
        if(owner[i]!=-1) throw std::runtime_error("Missing internal value");
        field.values[vtk_index[i]] = boundary_value;
    }
    return field;
}
