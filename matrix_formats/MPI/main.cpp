#include <mpi.h>
#include <memory>
#include <stdexcept>
#include <string>

#include "AbstractFormatDistributed.h"
#include "DDDNaiveDistributed.cpp"

#include <liblsb.h>
#include <stdlib.h>
#include <time.h>

#define RUNS 10


using namespace format;

// Usage: mpirun -np P ./smvm FORMAT MTX_IN VEC_IN [OTHER]

std::unique_ptr<AbstractFormatDistributed> createFormat(
        int argc, char **argv, int mpi_size, int mpi_rank, int mpi_root)
{
    std::string format(argv[1]);
    std::string mtx_in(argv[2]);
    std::string vec_in(argv[3]);

    if (format == "DDDNaive") {
        return std::unique_ptr<AbstractFormatDistributed>(
            new DDDNaiveDistributed({mtx_in, vec_in}, mpi_size,
                                    mpi_rank, mpi_root)
        );
    } // Add yours

    else {
        throw std::invalid_argument("Format doesn't match any known ones");
    }
}


int main(int argc, char **argv)
{
    if (argc < 4) {
        exit(-1);
    }

    MPI_Init(&argc, &argv);
    LSB_Init(/* proj_name */ argv[1], 0);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    LSB_Set_Rparam_int("rank", mpi_rank);
    LSB_Set_Rparam_int("runs", RUNS);

    for (int run = 0; run < RUNS; run++) {
        auto format = createFormat(argc, argv, mpi_size, mpi_rank, 0);
        format->InitData();
        format->ScatterData();
        
        // Only measure performance for the Multiplication.
        LSB_Res();
        format->PerformMultiplicationDistributed();
        LSB_Rec(run);
    }

    LSB_Finalize();
    MPI_Finalize();
    return 0;
}