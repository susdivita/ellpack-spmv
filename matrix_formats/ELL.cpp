#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>

struct Triplet {
    int row, col;
    double value;
    Triplet(int r, int c, double v) : row(r), col(c), value(v) {}
};
using index_t = std::ptrdiff_t;

namespace format {
class ELL {
public:
    ELL(const AbstractMatrixFormat::MatrixFiles &matrixFiles);
    std::vector<double> PerformMultiplicationDouble();

private:
    void buildELLrix(const std::vector<Triplet> &triplets, index_t m, index_t n);

    std::vector<double> val;  // vector containing values
    std::vector<index_t> column;  // vector containing column
    index_t maxcols;           // n of non-empty columns
    index_t m, n;              // n of rows, n of columns
    std::vector<double> x;                  // read vector x
    std::vector<double> res;
    std::vector<Triplet> triplets;  // entries
    int row, col;
    double value;
};

ELL::ELL(const AbstractMatrixFormat::MatrixFiles &matrixFiles) {
    // Read vector
    std::ifstream xFile(matrixFiles.vector_file);
    int xSize;
    xFile >> xSize;
    x.resize(xSize);
    for (int i = 0; i < xSize; ++i) {
        xFile >> x[i];
    }
    // Read matrix
    std::ifstream file(matrixFiles.matrix_file);
    file >> m >> n >> maxcols;
    column.resize(m * maxcols, -1);
    val.resize(m * maxcols, 0);
    while (file >> col >> row >> value) {
        triplets.push_back(Triplet(row - 1, col - 1, value));
    }
    buildELLrix(triplets, m, n); // Start building ELLPACK representation
    res.assign(m, 0);
}

std::vector<double> ELL::PerformMultiplicationDouble() {
    #pragma omp parallel for
    for (index_t i = 0; i < m; ++i) {
        double sum = 0.0;
        for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
            if (column[l] == -1)  break;
            sum += x[column[l]] * val[l];
        }
        #pragma omp critical
        res[i] = sum;
    }
    return res;
}

void ELL::buildELLrix(const std::vector<Triplet> &triplets, index_t m, index_t n) {
    std::vector<unsigned int> counters(m);
    std::fill(counters.begin(), counters.end(), 0);
    std::vector<std::vector<unsigned int>> counters_col(m);
    maxcols = 0;
    #pragma omp parallel
    {
        #pragma omp for reduction(max : maxcols)
        for (index_t i = 0; i < m; ++i) {
            std::vector<unsigned int> &col_idx = counters_col[i];
            for (const Triplet &tr : triplets) {
                #pragma omp cancellation point for
                if (tr.row == i) {
                    auto it = std::find(col_idx.begin(), col_idx.end(), tr.col);
                    if (it == col_idx.end()) {
                        col_idx.push_back(tr.col);
                        #pragma omp atomic update
                        ++counters[i];
                        if (counters[i] > maxcols) {
                            maxcols = counters[i];
                        }
                    }
                }
            }
        }
        #pragma omp for
        for (const Triplet &tr : triplets) {
            index_t l;
            bool found = false;
            for (l = tr.row * maxcols; l < (tr.row + 1) * maxcols; ++l) {
                if (column[l] == -1) {
                    #pragma omp atomic write
                    column[l] = tr.col;
                    #pragma omp atomic write
                    val[l] = tr.value;
                    found = true;
                    break;
                }
                if (column[l] == tr.col) {
                    #pragma omp atomic update
                    val[l] += tr.value;
                    found = true;
                    break;
                }
            }
        }
    }
}
}  // namespace format