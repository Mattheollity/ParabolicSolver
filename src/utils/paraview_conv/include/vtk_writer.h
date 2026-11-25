#pragma once

#include "field_reader.h"

class VTKWriter
{
    field_reader::Grid grid;
    std::filesystem::path directory;

    public:

    struct Layer
    {
        int step;
        double time;
    };

    VTKWriter(const field_reader::Grid & grid_,
              const std::filesystem::path & directory_);

    void writeVTI(const field_reader::Field & field) const;
    void writePVD(const std::vector<Layer> & layers) const;
};
