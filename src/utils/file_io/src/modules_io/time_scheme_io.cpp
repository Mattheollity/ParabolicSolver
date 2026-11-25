#include "time_scheme_io.h"

#include <mpi.h>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <stdexcept>


namespace
{
    template<typename DATATYPE>
    DATATYPE read_value(
        string FIELD,
        std::ifstream & file
    )
    {
        string field;
        DATATYPE value;

        if(!(file>>field>>value) || field!=FIELD)
        throw std::runtime_error("Invalid field "+FIELD);

        return value;
    }



    void write_values(
        string FIELD,
        const vect & v,
        std::ofstream & file
    )
    {
        file<<FIELD<<" "<<v.size()<<"\n";

        for(int i=0; i<v.size(); i++)
        file<<v[i]<<" ";

        file<<"\n";
    }



    vect read_values(
        string FIELD,
        Eigen::Index expected_size,
        std::ifstream & file
    )
    {
        Eigen::Index amount =
            read_value<Eigen::Index>(FIELD, file);

        if(amount!=expected_size)
        throw std::runtime_error("Invalid size of "+FIELD);

        vect v(amount);

        for(Eigen::Index i=0; i<amount; i++)
        if(!(file>>v[i]) || !std::isfinite(v[i]))
        throw std::runtime_error("Invalid value in "+FIELD);

        return v;
    }
}



void time_scheme_io::export_solution(
    const time_scheme_assembler & tsa,
    string path
) const
{
    if(tsa.rank!=0) return;

    std::filesystem::create_directories(path);

    std::ofstream file(path+"/solution");

    if(!file)
    throw std::runtime_error("Cannot export solution");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    file<<"time_step "<<tsa.current_time_step<<"\n";
    file<<"time "<<tsa.current_time<<"\n";
    file<<"T "<<tsa.T<<"\n";

    file<<"amount "<<tsa.current_values.size()<<"\n";

    for(int i=0; i<tsa.current_values.size(); i++)
    {
        const subdomain_values & v = tsa.current_values[i];

        file<<"subdomain "<<i<<"\n";

        write_values("inner", v.inner, file);

        write_values(
            "local_boundary",
            v.local_boundary,
            file
        );

        write_values(
            "global_boundary",
            v.global_boundary,
            file
        );
    }

    file.flush();

    if(!file)
    throw std::runtime_error("Cannot write solution");
}



void time_scheme_io::import_solution(
    time_scheme_assembler & tsa,
    string path
) const
{
    int success = 1;

    int time_step = 0;
    double time = 0.0;
    double saved_T = 0.0;

    if(tsa.rank==0)
    {
        try
        {
            std::ifstream file(path+"/solution");

            if(!file)
            throw std::runtime_error("Cannot import solution");

            time_step = read_value<int>("time_step", file);
            time = read_value<double>("time", file);
            saved_T = read_value<double>("T", file);

            if(time_step<0 ||
               !std::isfinite(time) ||
               !std::isfinite(saved_T) ||
               saved_T<0.0 ||
               (time_step>0 && saved_T==0.0) ||
               time<tsa.tsc.start_time ||
               time>tsa.tsc.end_time)
            throw std::runtime_error("Invalid saved time");

            int amount = read_value<int>("amount", file);

            if(amount<0 || amount!=tsa.current_values.size())
            throw std::runtime_error("Invalid subdomains amount");

            vector<subdomain_values> values(amount);

            for(int i=0; i<amount; i++)
            {
                int index = read_value<int>("subdomain", file);

                if(index!=i)
                throw std::runtime_error("Invalid subdomain index");

                subdomain_values & v = values[i];

                const subdomain_values & old =
                    tsa.current_values[i];

                v.inner = read_values(
                    "inner", old.inner.size(), file
                );

                v.local_boundary = read_values(
                    "local_boundary",
                    old.local_boundary.size(),
                    file
                );

                v.global_boundary = read_values(
                    "global_boundary",
                    old.global_boundary.size(),
                    file
                );
            }

            tsa.current_values.swap(values);
        }
        catch(const std::exception & e)
        {
            std::cerr<<e.what()<<std::endl;
            success = 0;
        }
    }



    if(tsa.process_count>1)
    MPI_Bcast(
        &success, 1, MPI_INT,
        0, MPI_COMM_WORLD
    );

    if(!success)
    throw std::runtime_error("Solution import failed");



    if(tsa.process_count>1)
    {
        MPI_Bcast(
            &time_step, 1, MPI_INT,
            0, MPI_COMM_WORLD
        );

        MPI_Bcast(
            &time, 1, MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );

        MPI_Bcast(
            &saved_T, 1, MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );
    }



    tsa.current_time = time;
    tsa.current_time_step = time_step;

    tsa.distribute_initial_values();

    tsa.next_values.clear();
    tsa.next_values.resize(tsa.current_values.size());



    double remaining_time = tsa.tsc.end_time-time;

    if(remaining_time>0.0)
    {
        double new_T = tsa.compute_global_time_step();

        new_T = std::min(new_T, remaining_time);

        if(new_T!=tsa.T)
        tsa.rebuild_time_step_systems(new_T);
    }
    else if(saved_T>0.0 && saved_T!=tsa.T)
    tsa.rebuild_time_step_systems(saved_T);
    else if(saved_T==0.0)
    tsa.T = 0.0;
}



void time_scheme_io::write_solution(
    const time_scheme_assembler & tsa
) const
{
    string path = wc.WRITE_DIR+"/"+to_string(
        tsa.current_time_step
    );

    export_solution(tsa, path);
}