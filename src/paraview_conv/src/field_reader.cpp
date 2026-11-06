
#include "field_reader.h"
#include <fstream>
#include <stdexcept>

field_reader::field_reader(const grid_conds& grid,
                         const time_integr_conds& time,
                         const std::filesystem::path& directory)
    : dimensions_(grid.N),
      steps_(time.N),
      writeEach_(time.write_each),
      startTime_(time.time.first),
      dt_(time.T),
      directory_(directory) {

    if (dimensions_.size() != 3 || writeEach_ <= 0 || steps_ <= 0)
        throw std::invalid_argument("Invalid grid or time conditions");

    for (int n : dimensions_)
        if (n < 2)
            throw std::invalid_argument("Invalid grid dimensions");
}

vect field_reader::readField(int step) const {
    std::ifstream file(directory_ / std::to_string(step));

    if (!file)
        throw std::runtime_error("Cannot open field: " + std::to_string(step));

    const auto& n = dimensions_;
    vect values = vect::Zero(n[0] * n[1] * n[2]);

    for (int z = 1; z < n[2] - 1; ++z)
        for (int y = 1; y < n[1] - 1; ++y)
            for (int x = 1; x < n[0] - 1; ++x)
                if (!(file >> values[x + n[0] * (y + n[1] * z)]))
                    throw std::runtime_error("Invalid field: " + std::to_string(step));

    double extra;
    if (file >> extra || !file.eof())
        throw std::runtime_error("Unexpected field data: " + std::to_string(step));

    return values;
}

std::vector<field_reader::Field> field_reader::readFields() const {
    std::vector<Field> fields;

    for (int step = 0; step < steps_; step += writeEach_)
        fields.push_back({
            step,
            startTime_ + step * dt_,
            readField(step)
        });

    return fields;
}
