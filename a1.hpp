// File: a1.hpp
// Aashish
// Jain

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <mpi.h>

// IMPLEMENT ME!

//Added
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define TAG 0

template <typename T, typename Pred>
void mpi_extract_if(MPI_Comm comm, const std::vector<T> &in, std::vector<T> &out, Pred pred)
{

    //To store the number of ranks and number of processors
    int n_ranks, rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &n_ranks);

    MPI_Status stat;
    int size_of_T = sizeof(T);

    std::vector<int> out_vec_sizes(n_ranks);
    int out_size, goal, elements_to_retain;

    //Get out vector
    for (int i = 0; i < in.size(); i++)
        if (pred(in[i]))
            out.push_back(in[i]);

    //Storing the rounded out_size and elements to retain
    out_size = out.size();
    elements_to_retain = out_size % n_ranks;
    out_size -= elements_to_retain;

    //Gather the rounded number of elements on process
    MPI_Allgather(&out_size, 1, MPI_INT, out_vec_sizes.data(), 1, MPI_INT, comm);

    //Getting total number of elements
    //Long for total to prevent overflows
    long total_elements = 0;
    for (int i = 0; i < n_ranks; i++)
        total_elements += out_vec_sizes[i];
    goal = total_elements / n_ranks;

    //Resize out vector to goal if size is less (to store incoming data)
    if (out_size < goal)
    {
        //Some elements are retained in the process itself
        out.resize(goal + elements_to_retain);
        out_size = out_size + elements_to_retain;
    }

    //To distribute is negative if a rank needs elements
    std::vector<int> to_distribute(n_ranks);
    for (int i = 0; i < n_ranks; i++)
        to_distribute[i] = out_vec_sizes[i] - goal;

    //Position of element in out that has to be sent
    int send_offset = goal + elements_to_retain;
    //Sender loop
    for (int r = 0; r < n_ranks; ++r)
    {
        int &excess_elements = to_distribute[r];
        //If a process needs elements or has reached goal(converged)
        if (excess_elements <= 0)
            continue;

        //First fit algorithm
        //Emulate transfer to keep track of of indices to send/recieve
        //Reciever loop
        for (int j = 0, n_send, needed_elements; j < n_ranks; j++)
        {

            //To distribute is negative if a rank needs elements
            needed_elements = to_distribute[j] * -1;
            //n send will be minimum of needed and excess but needed can be negative - checking
            n_send = MIN(excess_elements, MAX(needed_elements, 0));

            //If there are elements to send then send/recieve
            if (n_send)
            {
                //If rank is the sender
                if (rank == r)
                {
                    MPI_Send(out.data() + send_offset, n_send * size_of_T, MPI_BYTE, j, TAG, comm);
                    //Update send_offset
                    send_offset += n_send;
                }
                //If rank is the reciever
                else if (rank == j)
                {
                    MPI_Recv(out.data() + out_size, n_send * size_of_T, MPI_BYTE, r, TAG, comm, &stat);
                    //Update out_size
                    out_size += n_send;
                }
            }
            //Assume that elements are sent to j
            excess_elements -= n_send;
            to_distribute[j] += n_send;
        }
    }
    //Resizing vectors with exess elements
    if (out_size > goal)
        //Some elements are retained in the process itself
        out.resize(goal + elements_to_retain);
} // mpi_extract_if

//Computationally more expensive function
template <typename T, typename Pred>
void mpi_extract_if_scatter(MPI_Comm comm, const std::vector<T> &in, std::vector<T> &out, Pred pred)
{

    //To store the rank and number of processors
    int n_proc, rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &n_proc);

    int size_of_T = sizeof(T);

    std::vector<T> out_prime;

    int out_size, elements_to_keep, elements_to_transfer;
    int *out_vec_sizes = new int[n_proc];
    int *write_indices = new int[n_proc + 1];

    //Get out_prime vector
    for (int i = 0; i < in.size(); i++)
        if (pred(in[i]))
            out_prime.push_back(in[i]);

    //Round of the out to make it a multiple of size
    out_size = out_prime.size();
    elements_to_keep = out_size % n_proc;
    out_size -= elements_to_keep;
    elements_to_transfer = out_size / n_proc;

    //Gather the rounded number of element to recieve from each process
    //gather computer and broadcast or broadcast??
    MPI_Allgather(&out_size, 1, MPI_INT, &out_vec_sizes[0], 1, MPI_INT, comm);

    //Prefix sum operation to get indices to read from the window
    write_indices[0] = 0;
    for (int i = 1; i <= n_proc; i++)
        write_indices[i] = write_indices[i - 1] + out_vec_sizes[i - 1] / n_proc;
    out.resize(write_indices[n_proc] + elements_to_keep);

    //Given POD so memcpy is used

    //Copy its own share of elements to the out
    memcpy(&out[write_indices[rank]], &out_prime[elements_to_transfer * rank], elements_to_transfer * size_of_T);
    //Copy remaining elements (elements to keep)
    memcpy(&out[write_indices[n_proc]], &out_prime[out_size], elements_to_keep * size_of_T);

    for (int r = 0; r < n_proc; r++)
    {
        long expected_bytes = out_vec_sizes[r] / n_proc * size_of_T;
        MPI_Scatter(&out_prime.front(), expected_bytes, MPI_BYTE,
                    &out.front() + write_indices[r], expected_bytes, MPI_BYTE, r, comm);
    }
} // mpi_extract_if_scatter


#endif // A1_HPP
