#ifndef TIME_SCHEME_IO_
#define TIME_SCHEME_IO_

#include "config_structures_ts.h"
#include "time_scheme_assembler.h"

struct time_scheme_io
{
    const writter_config & wc;
    int time_step;
    int N;

    time_scheme_io(const writter_config & wc_) : wc(wc_), time_step(0) {};
    void write_solution(const time_scheme_assembler &);
};

#endif