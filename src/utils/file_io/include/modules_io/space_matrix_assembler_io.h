#ifndef SPACE_MATRIX_ASSEMBLER_IO_
#define SPACE_MATRIX_ASSEMBLER_IO_

#include "common.h"
#include "space_matrix_assembler.h"
#include "configs_io.h"
#include "eigen_io.h"

class space_matrix_assembler_io
{
    public:

    void import_matrixes(space_matrix_assembler & sma, string path);
    void export_matrixes(const space_matrix_assembler & sma, string path) const;
};

#endif