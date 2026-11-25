#pragma once

#include <array>
#include <filesystem>
#include <vector>

class field_reader
{
    struct subdomain
    {
        std::vector<int> inner, boundary, global_boundary;
    };

    std::vector<subdomain> subdomains;
    std::vector<int> owner;
    std::vector<std::size_t> vtk_index;
    std::filesystem::path directory;
    double boundary_value;

    public:

    struct Grid
    {
        std::array<std::size_t, 3> dimensions;
        std::array<double, 3> origin, spacing;
    };

    struct Field
    {
        int step;
        double time;
        double T;
        std::vector<double> values;
    };

    Grid grid;

    field_reader(const std::filesystem::path & case_dir);

    std::vector<int> find_steps() const;

    Field readField(int step) const;
};
