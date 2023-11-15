#include <cassert>
#include <string.h>
#include <vector>
#include <omp.h>

#include <iostream>

#include "AbstractMatrixFormat.h"

namespace format {

class DDDNaive : public AbstractMatrixFormat {
public:
    DDDNaive(const MatrixFiles &matrix_files) : AbstractMatrixFormat(matrix_files) {
        Init(matrix_files);
    }

    ~DDDNaive() {
        free(vec_val);
        free(diag_val);
        free(diag_num);
        free(diag_ptr);
    }

    std::vector<float> PerformMultiplicationFloat() {
        return {};
    }

    std::vector<int> PerformMultiplicationInt32 () {
        return {};
    }

    std::vector<double> PerformMultiplicationDouble() {
        auto dim = GetDimensions();
        const int N = dim.mat_rows;

        double *ans = (double *) malloc(N * sizeof(double));
        memset(ans, 0, N * sizeof(double));

        #pragma omp parallel
        {
            double *private_ans = (double *) malloc(N * sizeof(double));
            memset(private_ans, 0, N * sizeof(double));

            #pragma omp for
            for (int d = 0; d < nz_diag; d++) {
                int row = 0, col = 0;
                int steps;
                if (diag_num[d] > 0) {
                    col = diag_num[d];
                    steps = N - diag_num[d];
                } else {
                    row = -diag_num[d];
                    steps = N + diag_num[d];
                }

                for (int j = 0; j < steps; j++) {
                    private_ans[row + j] += diag_val[diag_ptr[d] + j] * vec_val[col + j];
                }
            }

            for (int i = 0; i < N; i++) {
                #pragma omp atomic
                ans[i] += private_ans[i];
            }

            free(private_ans);
        }

        return std::vector<double>(ans, ans + dim.mat_rows);
    }

private:
    void Init(const MatrixFiles &matrix_files) {
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
                    diag_ptr[idx] = diag_ptr[idx-1] + dim.mat_rows - abs(diag_num[idx-1]);
                }
                if (diag_ptr[idx] + dim.mat_rows - abs(diag_num[idx]) > diag_val_capacity) {
                    diag_val = (double *) realloc(diag_val, 2 * diag_val_capacity * sizeof(double));
                    diag_val_capacity *= 2;
                }
                nz_diag++;
            }

            // already reserved space for this diagonal, just find the offset
            int offset = diag_ptr[idx] + ((diag > 0) ? row : col);
            diag_val[offset] = val;
        }

        //////
        // for (int d = 0; d < nz_diag; d++) {
        //     std::cout << "DIAG: " << diag_num[d] << ", PTR: " << diag_ptr[d] << std::endl;
        //     for (int i = 0; i < dim.mat_rows - abs(diag_num[d]); i++) {
        //         std::cout << diag_val[diag_ptr[d] + i] << std::endl;
        //     }
        //     std::cout << "-----------" << std::endl;
        // }
        //////

        fclose(fp_mtx);
        SetDimensions(dim);
    }

    int nz_num; // non-zero values
    double *vec_val;
    double *diag_val;
    int diag_val_capacity;
    int *diag_num, *diag_ptr;
    int max_diag; // maximum number of diagonals
    int nz_diag; // number of diagonals with non-zero values
};

} // namespace format