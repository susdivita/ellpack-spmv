#include <vector>
#include <iostream>
#include <fstream>
#include "AbstractMatrixFormat.h"

namespace format {

struct Triplet {
    int row, col;
    double value;
    Triplet(int r, int c, double v) : row(r), col(c), value(v) {}
};

using Triplets = std::vector<Triplet>;
using Vector = std::vector<double>;
using index_t = std::ptrdiff_t;

class ELL : public AbstractMatrixFormat {
public:
    ELL(const MatrixFiles &matrix_files) : AbstractMatrixFormat(matrix_files) {
        Init(matrix_files);
    }

    ~ELL() {}

    std::vector<double> PerformMultiplicationDouble() {
        assert(!x.empty() && "Vector x is empty!");

        for (index_t i = 0; i < m; ++i) {
            res.at(i) = 0;
            for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
                if (columns[l] == -1) break;
                res.at(i) += x.at(columns[l]) * val[l];
                std::cout << res.at(i);
            }
        }
        return res;
    }

    std::vector<float> PerformMultiplicationFloat() { return {}; }

    std::vector<int32_t> PerformMultiplicationInt32() { return {}; }

private:
    void Init(const MatrixFiles &matrix_files) {
        // Read vector
        std::ifstream xFile("/workspaces/ellpack-spmv/inputs/test_basic_vec_in.txt");
        if (!xFile) {
            std::cerr << "Error opening x file" << std::endl;
        }

        int xSize;
        xFile >> xSize;

        x.resize(xSize);
        for (int i = 0; i < xSize; ++i) {
            xFile >> x[i];
        }
        std::cout << "vector size: " << x.size();

        // Read matrix
        Triplets triplets;
        std::ifstream file("/workspaces/ellpack-spmv/inputs/test_basic_mat_in.txt");
        if (!file) {
            std::cerr << "Error opening file" << std::endl;
        }
        
        file >> m >> n >> ntriplets;

        // Reserve space for triplets
        triplets.reserve(ntriplets);

        int row, col;
        double value;
        while (file >> col >> row >> value) {
            triplets.push_back(Triplet(row - 1, col - 1, value));
        }

        // Build Ellpack matrix
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

        std::cout << "Maxcols: " << maxcols << std::endl;
        
        columns.resize(m * maxcols, -1);
        val.resize(m * maxcols, 0);
        res.resize(m, 0);
        
        for (const Triplet &tr : triplets) {
            assert(0 <= tr.row && tr.row < m && 0 <= tr.col && tr.col < n &&
                   "Index out of bounds!");

            index_t l;
            for (l = tr.row * maxcols; l < (tr.row + 1) * maxcols; ++l) {
                if (columns[l] == -1) {
                    columns[l] = tr.col;
                    val[l] = tr.value;
                    break;
                }

                if (columns[l] == tr.col) {
                    val[l] += tr.value;
                    break;
                }
            }
            assert(l < (tr.row + 1) * maxcols &&
                   "You did not reserve enough columns!");
        }

        std::cout << "vector size: " << x.size();
        std::cout << "n size: " << n;
        std::cout << "res size: " << res.size();
        std::cout << "m size: " << m;

        std::cout << " ------------- Test of y = E*x --------------- " << std::endl;
        // E.mvmult(x, res);
        assert((x.size()) == n && "Incompatible vector x size!");
        assert(res.size() == m && "Incompatible vector y size!");
    }

    std::vector<double> val;  //< Vector containing values
    std::vector<index_t> columns;  //< Vector containing column
    index_t maxcols;  //< Number of non-empty columns
    index_t m, n;
    Vector x;  //< Vector x
    Vector res;
};

}  // namespace format
