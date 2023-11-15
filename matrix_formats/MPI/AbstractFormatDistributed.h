#ifndef ABSTRACT_FORMAT_DISTRIBUTED_H
#define ABSTRACT_FORMAT_DISTRIBUTED_H

#include <vector>

#include "../AbstractMatrixFormat.h"

namespace format {

class AbstractFormatDistributed : public AbstractMatrixFormat {
public:
    AbstractFormatDistributed(
            const MatrixFiles &matrix_files,
            int mpi_size,
            int mpi_rank,
            int mpi_root)
        : AbstractMatrixFormat(matrix_files)
        , mpi_size(mpi_size)
        , mpi_rank(mpi_rank)
        , mpi_root(mpi_root)
    {}

    // Delete default constructor.
    AbstractFormatDistributed() = delete;

    virtual ~AbstractFormatDistributed() = default;

    // ONLY ROOT PROCESS SHOULD DO THIS!!
    // Initialize the data, i.e. read it from file.
    virtual void InitData() = 0;

    // Send data to each process.
    virtual void ScatterData() = 0;

    // Perform multiplication.
    // Assume InitData() and ScatterData() were run before.
    // Root process should store the answer in ans.
    virtual void PerformMultiplicationDistributed() = 0;

    // Return answer in a vector so we can use our old testing method.
    std::vector<double> GetAnswer() {
        if (mpi_rank == mpi_root) {
            return std::vector<double>(ans, ans + GetDimensions().mat_rows);
        } else {
            return {};
        }
    }

    // UNUSED METHODS:
    std::vector<double> PerformMultiplicationDouble() { return {}; }
    std::vector<float> PerformMultiplicationFloat() { return {}; }
    std::vector<int32_t> PerformMultiplicationInt32() { return {}; }

protected:
    int mpi_size;   // the number of processes
    int mpi_rank;   // the rank of the current process
    int mpi_root;   // the rank of the "Master" process (in case it's not 0)

    double *ans = NULL;    // store answer in root process ONLY
};

}

#endif