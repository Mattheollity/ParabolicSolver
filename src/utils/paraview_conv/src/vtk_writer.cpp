
#include "vtk_writer.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace fs = std::filesystem;

VTKWriter::VTKWriter(const grid_conds& grid, const fs::path& directory)
    : grid_(grid), directory_(directory) {}

void VTKWriter::write(const std::vector<field_reader::Field>& fields) const {
    fs::create_directories(directory_);

    for (const auto& field : fields)
        writeVTI(field);

    writePVD(fields);
}

void VTKWriter::writeVTI(const field_reader::Field& field) const {
    const auto& n = grid_.N;
    const auto& bbox = grid_.bbox;

    if (n.size() != 3 || bbox.first.size() != 3 || bbox.second.size() != 3)
        throw std::runtime_error("Invalid grid dimensions");

    Eigen::Index size = 1;
    for (int count : n) {
        if (count < 2)
            throw std::runtime_error("Invalid grid size");
        size *= count;
    }

    if (field.values.size() != size)
        throw std::runtime_error("Invalid field size");

    const fs::path path = directory_ /
        ("step_" + std::to_string(field.step) + ".vti");

    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot create VTI file");

    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        << "<ImageData WholeExtent=\"0 " << n[0]
        << " 0 " << n[1]
        << " 0 " << n[2] << "\" Origin=\"";

    for (int i = 0; i < 3; ++i) {
        const double h = (bbox.second[i] - bbox.first[i]) / (n[i] - 1);
        out << bbox.first[i] - h / 2 << (i == 2 ? "\"" : " ");
    }

    out << " Spacing=\"";

    for (int i = 0; i < 3; ++i) {
        const double h = (bbox.second[i] - bbox.first[i]) / (n[i] - 1);
        out << h << (i == 2 ? "\">\n" : " ");
    }

    out << "<Piece Extent=\"0 " << n[0]
        << " 0 " << n[1]
        << " 0 " << n[2] << "\">\n"
        << "<CellData Scalars=\"u\">\n"
        << "<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";

    for (Eigen::Index i = 0; i < field.values.size(); ++i)
        out << field.values[i] << '\n';

    out << "</DataArray>\n"
        << "</CellData>\n"
        << "</Piece>\n"
        << "</ImageData>\n"
        << "</VTKFile>\n";
}

void VTKWriter::writePVD(
    const std::vector<field_reader::Field>& fields) const {

    std::ofstream out(directory_ / "solution.pvd");
    if (!out)
        throw std::runtime_error("Cannot create PVD file");

    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        << "<Collection>\n";

    for (const auto& field : fields)
        out << "<DataSet timestep=\"" << field.time
            << "\" group=\"\" part=\"0\" file=\"step_"
            << field.step << ".vti\"/>\n";

    out << "</Collection>\n"
        << "</VTKFile>\n";
}
