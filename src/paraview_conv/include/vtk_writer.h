
#pragma once

#include "field_reader.h"
#include <filesystem>

class VTKWriter {
public:
    VTKWriter(const grid_conds& grid,
              const std::filesystem::path& directory);

    void write(const std::vector<field_reader::Field>& fields) const;

private:
    grid_conds grid_;
    std::filesystem::path directory_;

    void writeVTI(const field_reader::Field& field) const;
    void writePVD(const std::vector<field_reader::Field>& fields) const;
};
