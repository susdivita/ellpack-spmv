#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
struct Triplet
{
  int  row, col;
  double value;
      Triplet(int r, int c, double v) : row(r), col(c), value(v) {}

};
using Triplets = std::vector<Triplet>;
using Vector = std::vector<double>;
using index_t = std::ptrdiff_t;



class EllpackMat {
 public:
  EllpackMat(const Triplets &triplets, index_t m, index_t n);
  double operator()(index_t i, index_t j) const;
  void mvmult(const Vector &x, Vector &y) const;

 private:
  std::vector<double> val;  //< Vector containing values
  std::vector<index_t> col;  //< Vector containing column
  index_t maxcols;  //< Number of non-empty columns
  index_t m, n;     //< Number of rows, number of columns
};

double EllpackMat::operator()(index_t i, index_t j) const {
  assert(0 <= i && i < m && 0 <= j && j < n && "Index out of bounds!");

  for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
    if (col.at(l) == j) return val.at(l);
  }
  return 0;
}

EllpackMat::EllpackMat(const Triplets &triplets, index_t m, index_t n)
    : m(m), n(n) {

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

  col.resize(m * maxcols, -1);
  val.resize(m * maxcols, 0);

  for (const Triplet &tr : triplets) {
    assert(0 <= tr.row && tr.row < m && 0 <= tr.col && tr.col < n &&
           "Index out of bounds!");

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
    assert(l < (tr.row + 1) * maxcols &&
           "You did not reserve enough columns!");
  }
}
void EllpackMat::mvmult(const Vector &x, Vector &y) const {
  assert(x.size() == n && "Incompatible vector x size!");
  assert(y.size() == m && "Incompatible vector y size!");
  for (index_t i = 0; i < m; ++i) {
    y.at(i) = 0;
    for (index_t l = i * maxcols; l < (i + 1) * maxcols; ++l) {
      if (col[l] == -1) break;
      y.at(i) += x.at(col[l]) * val[l];
    }
  }
}


/* SAM_LISTING_BEGIN_5 */
int main() {
  Triplets triplets;
  unsigned int m, n, ntriplets;
  std::ifstream file("/workspaces/SpMMV-DPHPC/inputs/test_basic_mat_in.txt");
  if (!file) {
    std::cerr << "Error opening file" << std::endl;
    return 1;
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
  EllpackMat E(triplets, m, n);
  
  std::cout << " ------------- Test: print  E ---------------- " << std::endl;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      std::cout << std::setw(4) << E(i, j) << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << " ------------- Test of y = E*x --------------- " << std::endl;
  std::ifstream xFile("/workspaces/SpMMV-DPHPC/inputs/test_basic_vec_in.txt");
  if (!xFile) {
    std::cerr << "Error opening x file" << std::endl;
    return 1;
  }

  int xSize;
  xFile >> xSize;

  Vector x(xSize);
  for (int i = 0; i < xSize; ++i) {
    xFile >> x[i];
  }

  Vector res(m,0);
  E.mvmult(x, res);
std::cout << "Ellpack E*x =" << std::endl;
for (const auto& element : res) {
    std::cout << element << " ";
}
std::cout << std::endl;  
}
/* SAM_LISTING_END_5 */