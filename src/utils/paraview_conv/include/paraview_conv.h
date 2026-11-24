
#pragma once

#include <filesystem>

class ParaViewConverter {
public:
    static void convert(
        const std::filesystem::path& gridFile,
        const std::filesystem::path& timeFile,
        const std::filesystem::path& fieldsDirectory,
        const std::filesystem::path& outputDirectory);
};
