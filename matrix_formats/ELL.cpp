#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>

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
    ELL(const std::string &matrixFile, const std::string &vectorFile);
    std::vector<double> PerformMultiplicationDouble();

private:
    void buildELLrix(const Triplets &triplets, index_t m, index_t n);

    std::vector<double> val;  //< Vector containing values
    std::vector<index_t> col;  //< Vector containing column
    index_t maxcols;           //< Number of non-empty columns
    index_t m, n;              //< Number of rows, number of columns
    Vector x;                  //< Vector x
    Vector res;
};

ELL::ELL(const std::string &matrixFile, const std::string &vectorFile) {
    // Read matrix
    Triplets triplets;
    std::ifstream file(matrixFile);
    if (!file) {
        std::cerr << "Error opening matrix file" << std::endl;
        return;
    }

    file >> m >> n >> maxcols;
    col.resize(m * maxcols, -1);
    val.resize(m * maxcols, 0);

    int row, col;
    double value;
    while (file >> col >> row >> value) {
        triplets.push_back(Triplet(row - 1, col - 1, value));
    }

    buildELLrix(triplets, m, n);

    // Read vector
    std::ifstream xFile(vectorFile);
    if (!xFile) {
        std::cerr << "Error opening vector file" << std::endl;
        return;
    }

    int xSize;
    xFile >> xSize;

    x.resize(xSize);
    for (int i = 0; i < xSize; ++i) {
        xFile >> x[i];
    }

    res.assign(m, 0);
}

std::vector<double> ELL::PerformMultiplicationDouble() {
    assert(!x.empty() && "Vector x is empty!");

    for (index_t i = 0; i < m; ++i) {
        res.at(i) = 0;
        for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
            if (col[l] == -1)
                break;
            res.at(i) += x.at(col[l]) * val[l];
        }
    }
    return res;
}

void ELL::buildELLrix(const Triplets &triplets, index_t m, index_t n) {
    std::vector<unsigned int> counters(m);
    std::fill(counters.begin(), counters.end(), 0);
    std::vector<std::vector<unsigned int>> counters_col(m);
    maxcols = 0;
    for (const Triplet &tr : triplets) {
        std::vector<unsigned int> &col_idx = counters_col[tr.row];
        if (std::find(col_idx.begin(), col_idx.end(), tr.col) == col_idx.end()) {
            col_idx.push_back(tr.col);
            if (++counters[tr.row] > maxcols) {
                maxcols = counters[tr.row];
            }
        }
    }

    col.resize(m * maxcols, -1);
    val.resize(m * maxcols, 0);

    for (const Triplet &tr : triplets) {
        assert(0 <= tr.row && tr.row < m && 0 <= tr.col && tr.col < n && "Index out of bounds!");

        index_t l;
        for (l = tr.row * maxcols; l < (tr.row + 1) * maxcols; ++l) {
            if (col[l] == -1) {
                col[l] = tr.col;
                val[l] = tr.value;
                break;
            }

            if (col[l] == tr.col) {
                val[l] += tr.value;
                break;
            }
        }
        assert(l < (tr.row + 1) * maxcols && "You did not reserve enough columns!");
    }
}
}