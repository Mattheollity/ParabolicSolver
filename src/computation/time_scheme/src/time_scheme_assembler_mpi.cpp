#include "time_scheme_assembler.h"
#include "mpi_utils.h"


namespace
{
    enum mpi_tag
    {
        MATRIX_II = 0,
        MATRIX_IB = 2,
        MATRIX_BI = 4,
        MATRIX_BB = 6,
        MATRIX_GB = 8,

        INITIAL_INNER = 10,
        INITIAL_LOCAL_BOUNDARY = 12,
        INITIAL_GLOBAL_BOUNDARY = 14,

        SCHUR = 16,
        GLOBAL_BOUNDARY = 18,
        SCHUR_RHS = 20,
        LOCAL_BOUNDARY = 22,
        INNER = 24
    };
}



void time_scheme_assembler::distribute_subdomain_systems(
    const space_matrix_assembler & sma,
    int amount
)
{
    if(process_count==1)
    {
        const auto & systems = sma.get_subdomain_systems();

        subdomains_solvers.reserve(systems.size());

        for(const subdomain_system & ss : systems)
        subdomains_solvers.emplace_back(*tsc.scheme, ss);

        return;
    }

    if(rank==0)
    {
        const auto & systems = sma.get_subdomain_systems();

        for(int i=0; i<amount; i++)
        {
            int dest = i+1;

            const subdomain_system & ss = systems[i];

            sendMatrix(ss.ii, dest, MATRIX_II, MPI_COMM_WORLD);
            sendMatrix(ss.ib, dest, MATRIX_IB, MPI_COMM_WORLD);
            sendMatrix(ss.bi, dest, MATRIX_BI, MPI_COMM_WORLD);
            sendMatrix(ss.bb, dest, MATRIX_BB, MPI_COMM_WORLD);
            sendMatrix(ss.gb_subd, dest, MATRIX_GB, MPI_COMM_WORLD);
        }
    }
    else if(rank<=amount)
    {
        subdomain_system ss;

        recvMatrix(ss.ii, 0, MATRIX_II, MPI_COMM_WORLD);
        recvMatrix(ss.ib, 0, MATRIX_IB, MPI_COMM_WORLD);
        recvMatrix(ss.bi, 0, MATRIX_BI, MPI_COMM_WORLD);
        recvMatrix(ss.bb, 0, MATRIX_BB, MPI_COMM_WORLD);
        recvMatrix(ss.gb_subd, 0, MATRIX_GB, MPI_COMM_WORLD);

        subdomains_solvers.emplace_back(
            *tsc.scheme,
            std::move(ss)
        );
    }
}



void time_scheme_assembler::distribute_initial_values()
{
    if(process_count==1) return;

    if(rank==0)
    {
        for(int i=0; i<current_values.size(); i++)
        {
            int dest = i+1;

            const subdomain_values & v = current_values[i];

            sendVector(
                v.inner, dest,
                INITIAL_INNER, MPI_COMM_WORLD
            );

            sendVector(
                v.local_boundary, dest,
                INITIAL_LOCAL_BOUNDARY, MPI_COMM_WORLD
            );

            sendVector(
                v.global_boundary, dest,
                INITIAL_GLOBAL_BOUNDARY, MPI_COMM_WORLD
            );
        }
    }
    else
    {
        current_values.resize(subdomains_solvers.size());

        if(subdomains_solvers.empty()) return;

        subdomain_values & v = current_values[0];

        recvVector(
            v.inner, 0,
            INITIAL_INNER, MPI_COMM_WORLD
        );

        recvVector(
            v.local_boundary, 0,
            INITIAL_LOCAL_BOUNDARY, MPI_COMM_WORLD
        );

        recvVector(
            v.global_boundary, 0,
            INITIAL_GLOBAL_BOUNDARY, MPI_COMM_WORLD
        );
    }
}



vector<vector<triplet>>
time_scheme_assembler::collect_schur_complements() const
{
    if(process_count==1)
    {
        vector<vector<triplet>> schur(subdomains_solvers.size());

        for(int i=0; i<subdomains_solvers.size(); i++)
        subdomains_solvers[i].build_schur_complement(schur[i]);

        return schur;
    }

    if(rank==0)
    {
        vector<vector<triplet>> schur(current_values.size());

        for(int i=0; i<schur.size(); i++)
        recvTriplets(
            schur[i], i+1,
            SCHUR, MPI_COMM_WORLD
        );

        return schur;
    }

    if(!subdomains_solvers.empty())
    {
        vector<triplet> schur;

        subdomains_solvers[0].build_schur_complement(schur);

        sendTriplets(schur, 0, SCHUR, MPI_COMM_WORLD);
    }

    return {};
}



void time_scheme_assembler::distribute_global_boundary_values()
{
    if(process_count==1) return;

    if(rank==0)
    {
        for(int i=0; i<next_values.size(); i++)
        sendVector(
            next_values[i].global_boundary,
            i+1,
            GLOBAL_BOUNDARY,
            MPI_COMM_WORLD
        );
    }
    else if(!subdomains_solvers.empty())
    {
        recvVector(
            next_values[0].global_boundary,
            0,
            GLOBAL_BOUNDARY,
            MPI_COMM_WORLD
        );
    }
}



vector<vect> time_scheme_assembler::collect_schur_rhs() const
{
    if(process_count==1)
    {
        vector<vect> rhs(subdomains_solvers.size());

        for(int i=0; i<subdomains_solvers.size(); i++)
        rhs[i] = subdomains_solvers[i].build_schur_rhs(
            current_values[i],
            next_values[i]
        );

        return rhs;
    }

    if(rank==0)
    {
        vector<vect> rhs(current_values.size());

        for(int i=0; i<rhs.size(); i++)
        recvVector(
            rhs[i], i+1,
            SCHUR_RHS, MPI_COMM_WORLD
        );

        return rhs;
    }

    if(!subdomains_solvers.empty())
    {
        vect rhs = subdomains_solvers[0].build_schur_rhs(
            current_values[0],
            next_values[0]
        );

        sendVector(rhs, 0, SCHUR_RHS, MPI_COMM_WORLD);
    }

    return {};
}



void time_scheme_assembler::distribute_local_boundary_values()
{
    if(process_count==1) return;

    if(rank==0)
    {
        for(int i=0; i<next_values.size(); i++)
        sendVector(
            next_values[i].local_boundary,
            i+1,
            LOCAL_BOUNDARY,
            MPI_COMM_WORLD
        );
    }
    else if(!subdomains_solvers.empty())
    {
        recvVector(
            next_values[0].local_boundary,
            0,
            LOCAL_BOUNDARY,
            MPI_COMM_WORLD
        );
    }
}



void time_scheme_assembler::collect_inner_values()
{
    if(process_count==1) return;

    if(rank==0)
    {
        for(int i=0; i<next_values.size(); i++)
        recvVector(
            next_values[i].inner,
            i+1,
            INNER,
            MPI_COMM_WORLD
        );
    }
    else if(!subdomains_solvers.empty())
    {
        sendVector(
            next_values[0].inner,
            0,
            INNER,
            MPI_COMM_WORLD
        );
    }
}