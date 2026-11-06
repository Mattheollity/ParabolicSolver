
#pragma once

#include "problem_conditions.h"
#include <filesystem>

class field_reader {
public:
    struct Field {
        int step;
        double time;
        vect values;
    };

    field_reader(const grid_conds& grid,
                const time_integr_conds& time,
                const std::filesystem::path& directory);

    std::vector<Field> readFields() const;

private:
    voxel_point dimensions_;
    int steps_;
    int writeEach_;
    double startTime_;
    double dt_;
    std::filesystem::path directory_;

    vect readField(int step) const;
};
