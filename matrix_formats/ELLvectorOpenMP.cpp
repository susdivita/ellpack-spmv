#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>

// Add OpenMP header
#include <omp.h>

struct Triplet {
    int row, col;
    double value;
    Triplet(int r, int c, double v) : row(r), col(c), value(v) {}
};
using Triplets = std::vector<Triplet>;
using Vector = std::vector<double>;
using index_t = std::ptrdiff_t;

namespace format {

class ELL {
public:
    ELL(const AbstractMatrixFormat::MatrixFiles &matrixFiles);
    std::vector<double> PerformMultiplicationDouble();

private:
    void buildELLrix(const Triplets &triplets, index_t m, index_t n);

    std::vector<double> val;  // Vector containing values
    std::vector<index_t> col;  // Vector containing column
    index_t maxcols;           // Number of non-empty columns
    index_t m, n;              // Number of rows, number of columns
    Vector x;                  // Vector x
    Vector res;
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
    Triplets triplets;
    std::ifstream file(matrixFiles.matrix_file);
    file >> m >> n >> maxcols;
    col.resize(m * maxcols, -1);
    val.resize(m * maxcols, 0);
    int row, col;
    double value;
    while (file >> col >> row >> value) {
        triplets.push_back(Triplet(row - 1, col - 1, value));
    }
    buildELLrix(triplets, m, n); // Start building ELLPACK representation
    res.assign(m, 0);
}

std::vector<double> ELL::PerformMultiplicationDouble() {
    assert(!x.empty() && "Vector x is empty!");

    #pragma omp parallel for
    for (index_t i = 0; i < m; ++i) {
        double sum = 0.0;
        for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
            if (col[l] == -1)
                break;
            sum += x[col[l]] * val[l];
        }
        #pragma omp critical
        res[i] = sum;
    }

    return res;
}



void ELL::buildELLrix(const Triplets &triplets, index_t m, index_t n) {
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
            assert(0 <= tr.row && tr.row < m && 0 <= tr.col && tr.col < n && "Index out of bounds!");

            index_t l;
            bool found = false;
            for (l = tr.row * maxcols; l < (tr.row + 1) * maxcols; ++l) {
                if (col[l] == -1) {
                    #pragma omp atomic write
                    col[l] = tr.col;
                    #pragma omp atomic write
                    val[l] = tr.value;
                    found = true;
                    break;
                }

                if (col[l] == tr.col) {
                    #pragma omp atomic update
                    val[l] += tr.value;
                    found = true;
                    break;
                }
            }

            assert(found && "You did not reserve enough columns!");
        }
    }
}


}  // namespace format
