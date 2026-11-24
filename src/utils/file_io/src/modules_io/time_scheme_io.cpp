#include "time_scheme_io.h"

void time_scheme_io::write_solution(const time_scheme_assembler & ts)
{
    N = (ts.tsc.end_time-ts.tsc.start_time)/ts.tsc.T;
    for(;time_step<N;time_step++)
    {
        
    }
}