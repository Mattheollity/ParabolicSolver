#include "time_scheme_assembler.h"

time_scheme_assembler::time_scheme_assembler(const time_scheme_config & tsc_) : tsc(tsc_), MPI_is_active(0)
{
    
}

void time_scheme_assembler::setup_solver(const vector<subdomain_system> & systems)
{
    thread_num = 0;
    MPI_Initialized(&MPI_is_active);

    if(MPI_is_active) reset_MPI(systems);
    else
    {
        int amount = systems.size();
        subdomains_solvers.resize(amount);

        for(int i=0;i<amount;i++)
        subdomains_solvers[i] = time_step_subdomain(tsc);

        for(int i=0;i<amount;i++)
        subdomains_solvers[i].setup_time_step_solver(systems[i]);
    }
}

void time_scheme_assembler::reset_MPI(const vector<subdomain_system> & systems)
{
    subdomains_per_MPI_threads.clear();

    if(!MPI_is_active) return;

    MPI_Comm comm = MPI_COMM_WORLD;
    int N_available_procs;
    MPI_Comm_size(comm, &N_available_procs);

    if(N_available_procs==1)
    {
        MPI_is_active=0;
        return;
    }

    N_available_procs--;
    subdomains_per_MPI_threads.resize(N_available_procs);
    
    if(N_available_procs>=systems.size())
        for(int i=0;i<systems.size();i++)
        subdomains_per_MPI_threads[i].push_back(i);
    else
        for(int i=0,j=0; i<systems.size();i++,j++)
        {
            throw simple_error("NOT ENOUGHT MPI PROCS");
            if(j==subdomains_per_MPI_threads.size()) j=0;
            subdomains_per_MPI_threads[j].push_back(i);
        }
}

void time_scheme_assembler::send_subdomain_systems(const vector<subdomain_system> & systems)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    for(int thread=1; thread<=subdomains_per_MPI_threads.size(); thread++)
    {
        int amount = subdomains_per_MPI_threads[thread-1].size();
        MPI_Send(&amount, 1, MPI_INT, thread, 0, comm);
    }

    for(int thread=1; thread<=subdomains_per_MPI_threads.size(); thread++)
    for(int i=0;i<subdomains_per_MPI_threads[thread-1].size(); i++)
    {
        int index = subdomains_per_MPI_threads[thread-1][i];
        const subdomain_system & ss = systems[index];

        sendMatrix(ss.ii, thread, 1+4*i, comm);
        sendMatrix(ss.ib, thread, 2+4*i, comm);
        sendMatrix(ss.bb, thread, 3+4*i, comm);
        sendMatrix(ss.gb_subd, thread, 4+4*i, comm);   
    }
}

void time_scheme_assembler::recv_subdomain_systems()
{
    subdomains_solvers.clear();

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &thread_num);

    int amount;
    MPI_Recv(&amount, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);

    for(int i=0;i<amount;i++)
    {
        subdomain_system ss;
        recvMatrix(ss.ii, 0, 1+4*i, comm);
        recvMatrix(ss.ib, 0, 1+4*i, comm);
        recvMatrix(ss.bb, 0, 1+4*i, comm);
        recvMatrix(ss.gb_subd, 0, 1+4*i, comm);
        subdomains_solvers[i] = time_step_subdomain(tsc);
        subdomains_solvers[i].setup_time_step_solver(ss);
    }
}

void time_scheme_assembler::restart_solver(const vector<subdomain_values> & start_values)
{
    if(MPI_is_active) send_subdomain_values(start_values);
    else
    {
        int amount = start_values.size();
        for(int i=0;i<amount;i++) subdomains_solvers[i].setup_first_solutions({start_values[i]});
    }
}

void time_scheme_assembler::send_subdomain_values(const vector<subdomain_values> & values)
{
    MPI_Comm comm = MPI_COMM_WORLD;

    for(int thread=1; thread<=subdomains_per_MPI_threads.size(); thread++)
    for(int i=0;i<subdomains_per_MPI_threads[thread-1].size(); i++)
    {
        int padding = 4*subdomains_per_MPI_threads[thread-1].size();
        int index = subdomains_per_MPI_threads[thread-1][i];
        const subdomain_system & ss = systems[index];

        sendMatrix(ss.ii, thread, 1+4*i, comm);
        sendMatrix(ss.ib, thread, 2+4*i, comm);
        sendMatrix(ss.bb, thread, 3+4*i, comm);
        sendMatrix(ss.gb_subd, thread, 4+4*i, comm);   
    }
}

void time_step(vector<subdomain_values> &);