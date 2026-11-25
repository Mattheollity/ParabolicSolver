#include "vtk_writer.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace fs = std::filesystem;

VTKWriter::VTKWriter(const field_reader::Grid & grid_, const fs::path & directory_) :
grid(grid_), directory(directory_)
{
    fs::create_directories(directory);
}


void VTKWriter::writeVTI(const field_reader::Field & field) const
{
    const auto & n = grid.dimensions;
    if(field.values.size()!=n[0]*n[1]*n[2])
    throw std::runtime_error("Invalid field size");

    std::ofstream out(directory/("step_"+std::to_string(field.step)+".vti"));
    if(!out) throw std::runtime_error("Cannot create VTI file");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);

    auto extent = [&]()
    {
        out<<"0 "<<n[0]-1<<" 0 "<<n[1]-1<<" 0 "<<n[2]-1;
    };

    out<<"<?xml version=\"1.0\"?>\n"
       <<"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
       <<"<ImageData WholeExtent=\"";
    extent();
    out<<"\" Origin=\""<<grid.origin[0]<<" "<<grid.origin[1]<<" "<<grid.origin[2]
       <<"\" Spacing=\""<<grid.spacing[0]<<" "<<grid.spacing[1]<<" "<<grid.spacing[2]<<"\">\n"
       <<"<Piece Extent=\"";
    extent();
    out<<"\">\n<PointData Scalars=\"u\">\n"
       <<"<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";
    for(double value : field.values) out<<value<<"\n";
    out<<"</DataArray>\n</PointData>\n<CellData/>\n</Piece>\n</ImageData>\n</VTKFile>\n";
    out.close();
    if(!out) throw std::runtime_error("Cannot write VTI file");
}


void VTKWriter::writePVD(const std::vector<Layer> & layers) const
{
    std::ofstream out(directory/"solution.pvd");
    if(!out) throw std::runtime_error("Cannot create PVD file");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"<?xml version=\"1.0\"?>\n"
       <<"<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n<Collection>\n";
    for(const Layer & layer : layers)
    out<<"<DataSet timestep=\""<<layer.time<<"\" group=\"\" part=\"0\" file=\"step_"
       <<layer.step<<".vti\"/>\n";
    out<<"</Collection>\n</VTKFile>\n";
    out.close();
    if(!out) throw std::runtime_error("Cannot write PVD file");
}
