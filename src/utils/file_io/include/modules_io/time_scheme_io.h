#ifndef TIME_SCHEME_IO_
#define TIME_SCHEME_IO_

#include "common.h"
#include "config_structures.h"
#include "time_scheme_assembler.h"


struct time_scheme_io
{
    const writter_config & wc;



    time_scheme_io(const writter_config & wc_) :
    wc(wc_) {};



    void import_solution(
        time_scheme_assembler & tsa,
        string path
    ) const;

    void export_solution(
        const time_scheme_assembler & tsa,
        string path
    ) const;

    void write_solution(
        const time_scheme_assembler & tsa
    ) const;
};

#endif