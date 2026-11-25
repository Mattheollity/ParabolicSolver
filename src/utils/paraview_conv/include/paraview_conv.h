#pragma once

#include <filesystem>

class ParaViewConverter
{
    public:

    static void convert(const std::filesystem::path & case_dir,
                        const std::filesystem::path & output_dir);
};
