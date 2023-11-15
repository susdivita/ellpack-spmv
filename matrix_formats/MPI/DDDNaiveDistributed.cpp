#include <mpi.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <assert.h>
#include <string>
#include <vector>

#include "AbstractFormatDistributed.h"


namespace format {

class DDDNaiveDistributed : public AbstractFormatDistributed {
public:
    DDDNaiveDistributed(
            const MatrixFiles &matrix_files,
            int mpi_size,
            int mpi_rank,
            int mpi_root = 0)
        : AbstractFormatDistributed(matrix_files, mpi_size, mpi_rank, mpi_root)
    {}

    ~DDDNaiveDistributed() {
        free(vec_val);
        free(l_diag_val);
        free(l_diag_num);
        free(l_diag_ptr);
        if (mpi_rank == mpi_root && ans != NULL) {
            free(ans);
        }
    }

    //////////////////////////
    // INITIALIZATION
    //////////////////////////

    void InitData() {
        // Only root process does this!
        if (mpi_rank != mpi_root) {
            return;
        }

        const MatrixFiles &matrix_files = GetMatrixFiles();
        FILE *fp_mtx, *fp_vec;
        Dimensions dim;

        fp_mtx = fopen(matrix_files.matrix_file.c_str(), "rb+");
        fp_vec = fopen(matrix_files.vector_file.c_str(), "rb+");

        // Read vector
        fscanf(fp_vec, "%d", &dim.vec_dim);
        vec_val = (double *) malloc(sizeof(double) * dim.vec_dim);
        for (int i = 0; i < dim.vec_dim; i++) {
            fscanf(fp_vec, "%lf", &vec_val[i]);
        }
        fclose(fp_vec);

        // Read matrix
        fscanf(fp_mtx, "%d %d %d", &dim.mat_rows, &dim.mat_columns, &nz_num);
        assert(dim.mat_rows == dim.mat_columns);

        max_diag = 2 * dim.mat_rows - 1;
        diag_val_capacity = nz_num; // might need to store more than these
        diag_val = (double *) malloc(sizeof(double) * diag_val_capacity);
        diag_num = (int *) malloc(sizeof(int) * max_diag);
        diag_ptr = (int *) malloc(sizeof(int) * max_diag);
        nz_diag = 0;

        for (int i = 0; i < nz_num; i++) {
            int row, col, diag;
            double val;
            fscanf(fp_mtx, "%d %d %lf", &col, &row, &val);
            row--;
            col--;
            diag = col - row;

            // Find diag in diag_num;
            int idx = 0;
            while (idx < nz_diag && diag_num[idx] != diag) {
                idx++;
            }

            if (idx == nz_diag) {
                // new diagonal discovered, reserve space for it
                diag_num[idx] = diag;
                diag_ptr[idx] = 0;
                if (idx > 0) {
                    diag_ptr[idx] = diag_ptr[idx-1] + dim.mat_rows - fabs(diag_num[idx-1]);
                }
                if (diag_ptr[idx] + dim.mat_rows - fabs(diag_num[idx]) > diag_val_capacity) {
                    diag_val = (double *) realloc(diag_val, 2 * diag_val_capacity * sizeof(double));
                    diag_val_capacity *= 2;
                }
                nz_diag++;
            }

            // already reserved space for this diagonal, just find the offset
            int offset = diag_ptr[idx] + ((diag > 0) ? row : col);
            diag_val[offset] = val;
        }

        fclose(fp_mtx);
        SetDimensions(dim);

        ans = (double *) malloc(dim.mat_rows * sizeof(double));
        memset(ans, 0, dim.mat_rows * sizeof(double));
    }

    void ScatterData() {
        int tmp[4];

        int *sendcounts = NULL;
        int *displs = NULL;

        if (mpi_rank == mpi_root) {
            N = GetDimensions().mat_rows;

            // Pack all dimensions, to do 1 single broadcast for all.
            tmp[0] = N;
            tmp[1] = nz_diag;

            MPI_Bcast(tmp, 4, MPI_INT, mpi_root, MPI_COMM_WORLD);

            // Send entire input vector.
            MPI_Bcast(vec_val, N, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);

            // Scatter the diagonals across processes.
            sendcounts = (int *) malloc(mpi_size * sizeof(int));
            displs = (int *) malloc(mpi_size * sizeof(int));
            diag_window = nz_diag / mpi_size;

            for (int i = 0; i < mpi_size - 1; i++) {
                sendcounts[i] = diag_window;
                displs[i] = i * diag_window;
            }
            sendcounts[mpi_size - 1] = diag_window + (nz_diag % mpi_size);
            displs[mpi_size - 1] = (mpi_size - 1) * diag_window;
        } else {
            MPI_Bcast(tmp, 4, MPI_INT, mpi_root, MPI_COMM_WORLD);
            N = tmp[0];
            nz_diag = tmp[1];

            // Get input vector.
            vec_val = (double *) malloc(N * sizeof(double));
            MPI_Bcast(vec_val, N, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);
        }

        // Compute diag window size for each process.
        diag_window = nz_diag / mpi_size;
        if (mpi_rank == mpi_size - 1) {
            diag_window += (nz_diag % mpi_size);
        }

        // Send diag IDs to each process.
        l_diag_num = (int *) malloc(diag_window * sizeof(int));
        MPI_Scatterv(diag_num, sendcounts, displs, MPI_INT,
                    l_diag_num, diag_window, MPI_INT, mpi_root, MPI_COMM_WORLD);

        // Compute diag_ptr.
        l_diag_ptr = (int *) malloc(diag_window * sizeof(int));
        l_diag_ptr[0] = 0;
        for (int i = 1; i < diag_window; i++) {
            l_diag_ptr[i] = l_diag_ptr[i-1] + (N - fabs(l_diag_num[i-1]));
        }
        val_window = l_diag_ptr[diag_window - 1] + (N - fabs(l_diag_num[diag_window - 1]));

        if (mpi_rank == mpi_root) {
            for (int i = 0; i < mpi_size - 1; i++) {
                sendcounts[i] = diag_ptr[(i+1) * diag_window] - diag_ptr[i*diag_window];
                displs[i] = diag_ptr[i*diag_window];
            }
            displs[mpi_size - 1] = diag_ptr[(mpi_size - 1) * diag_window];
            sendcounts[mpi_size - 1] =
                    (diag_ptr[nz_diag - 1] + N - fabs(diag_num[nz_diag - 1]))
                    - displs[mpi_size - 1];
        }

        // Send values.
        l_diag_val = (double *) malloc(val_window * sizeof(double));
        MPI_Scatterv(diag_val, sendcounts, displs, MPI_DOUBLE,
                    l_diag_val, val_window, MPI_DOUBLE, mpi_root, MPI_COMM_WORLD);


        if (mpi_rank == mpi_root) {
            free(sendcounts);
            free(displs);
            free(diag_num);
            free(diag_ptr);
            free(diag_val);
        }
    }

    /////////////////////////////
    // DISTRIBUTED COMPUTATION
    /////////////////////////////
    void PerformMultiplicationDistributed() {
        // Process-local answer vector:
        double *l_ans = (double *) malloc(N * sizeof(double));
        memset(l_ans, 0, N * sizeof(double));

        #pragma omp parallel
        {
            // Thread-private answer vector.
            double *p_ans = (double *) malloc(N * sizeof(double));
            memset(p_ans, 0, N * sizeof(double));

            #pragma omp for
            for (int d = 0; d < diag_window; d++) {
                int row=0, col=0;
                int steps;
                if (l_diag_num[d] > 0) {
                    col = l_diag_num[d];
                    steps = N - l_diag_num[d];
                } else {
                    row = -l_diag_num[d];
                    steps = N + l_diag_num[d];
                }

                for (int j = 0; j < steps; j++) {
                    p_ans[row + j] += l_diag_val[l_diag_ptr[d] + j] * vec_val[col + j];
                }
            }

            for (int i = 0; i < N; i++) {
                #pragma omp atomic
                l_ans[i] += p_ans[i];
            }
        }

        // Add all the answers.
        MPI_Reduce(l_ans, ans, N, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);

        free(l_ans);
    }

private:
    int nz_num; // non-zero values
    int nz_diag; // number of diagonals with non-zero values

    // Only in root:
    double *diag_val = NULL;
    int diag_val_capacity;
    int *diag_num = NULL, *diag_ptr = NULL;
    int max_diag;

    // Local values in each process:
    double *vec_val = NULL;
    double *l_diag_val = NULL;
    int *l_diag_num = NULL, *l_diag_ptr = NULL;
    int diag_window, val_window, N;
};

} // namespace format

