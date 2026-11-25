#include "paraview_conv.h"
#include "field_reader.h"
#include "vtk_writer.h"

#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

void ParaViewConverter::convert(const fs::path & case_dir, const fs::path & output_dir)
{
    field_reader reader(case_dir);
    VTKWriter writer(reader.grid, output_dir);
    std::vector<VTKWriter::Layer> layers;

    for(int step : reader.find_steps())
    {
        auto field = reader.readField(step);
        if(!layers.empty() && field.time<=layers.back().time)
        throw std::runtime_error("Non-increasing physical time");
        writer.writeVTI(field);
        layers.push_back({field.step, field.time});
    }
    writer.writePVD(layers);
}


int main(int argc, char ** argv)
{
    if(argc!=2 && argc!=3)
    {
        std::cerr<<"Usage: pvconv <case_directory> [output_directory]\n";
        return 1;
    }

    try
    {
        fs::path case_dir(argv[1]);
        fs::path output_dir = argc==3 ? fs::path(argv[2]) : case_dir/"paraview";
        ParaViewConverter::convert(case_dir, output_dir);
    }
    catch(const std::exception & e)
    {
        std::cerr<<"Conversion failed: "<<e.what()<<"\n";
        return 1;
    }
    return 0;
}
