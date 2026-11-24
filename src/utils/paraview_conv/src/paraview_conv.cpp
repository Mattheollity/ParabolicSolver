
#include "paraview_conv.h"
#include "field_reader.h"
#include "vtk_writer.h"
#include "file_io.h"

void ParaViewConverter::convert(
    const std::filesystem::path& gridFile,
    const std::filesystem::path& timeFile,
    const std::filesystem::path& fieldsDirectory,
    const std::filesystem::path& outputDirectory) {

    const auto grid = load_grid_conds(gridFile.string());
    const auto time = load_tic(timeFile.string());

    field_reader reader(grid, time, fieldsDirectory);
    VTKWriter writer(grid, outputDirectory);

    writer.write(reader.readFields());
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: converter <grid> <tic> <fields_dir> <output_dir>\n";
        return 1;
    }

    try {
        ParaViewConverter::convert(argv[1], argv[2], argv[3], argv[4]);
    }
    catch (const std::exception& e) {
        std::cerr << "Conversion failed: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
