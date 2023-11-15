#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "../matrix_formats/AbstractMatrixFormat.h"
#include "../matrix_formats/BasicMatrixFormat.cpp"
#include "../matrix_formats/CSR2.cpp"
#include "../matrix_formats/CSR.cpp"
#include "../matrix_formats/BSR.cpp"
#include "../matrix_formats/DDDNaive.cpp"
#include "../matrix_formats/ELL.cpp"
const std::size_t ULP_LIM = 10;

namespace format {

using MatrixFiles = AbstractMatrixFormat::MatrixFiles;

// Taken from docs:
// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool>
equal_within_ulps(T x, T y, std::size_t n) {
    // Since `epsilon()` is the gap size (ULP, unit in the last place)
    // of floating-point numbers in interval [1, 2), we can scale it to
    // the gap size in interval [2^e, 2^{e+1}), where `e` is the exponent
    // of `x` and `y.

    // If `x` and `y` have different gap sizes (which means they have
    // different exponents), we take the smaller one. Taking the bigger
    // one is also reasonable, I guess.
    T m = std::min(std::fabs(x), std::fabs(y));

    // Subnormal numbers have fixed exponent, which is `min_exponent - 1`.
    int exp = m < std::numeric_limits<T>::min()
                  ? std::numeric_limits<T>::min_exponent - 1
                  : std::ilogb(m);

    // We consider `x` and `y` equal if the difference between them is
    // within `n` ULPs.
    return std::fabs(x - y) <=
           n * std::ldexp(std::numeric_limits<T>::epsilon(), exp);
}

class AllTest : public testing::Test {
   protected:
    void SetUp() override {}

    void CompareWithBase(const AbstractMatrixFormat::MatrixFiles& matrix_files,
                         const std::vector<double>& other_result) {
        BasicMatrixFormat base_matrix(matrix_files);
        // Compute base
        std::vector<double> expected_result =
            base_matrix.PerformMultiplicationDouble();
        // Compare results with precision.
        for (int i = 0; i < expected_result.size(); ++i) {
            EXPECT_TRUE(equal_within_ulps(expected_result[i], other_result[i],
                                          ULP_LIM));
            if (!equal_within_ulps(expected_result[i], other_result[i], ULP_LIM)) {
                std::cout << "expected=" << expected_result[i]
                          << ", other=" << other_result[i]
                          << std::endl;
            }
        }
    }
};

TEST_F(AllTest, TestMatrixFormatBasicMatrixFormat) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/test_basic_mat_in.txt", "../inputs/test_basic_vec_in.txt"};
    // Create your matrix format here.
    // Compute matrix format
    std::vector<double> expected_result{22913.214,  10.363022,  21.811832,
                                        64087.9866, 2296.33526, 247.802478,
                                        32199.832};
    CompareWithBase(files, expected_result);
}

// ---------------------- CSR2 -------------------------

TEST_F(AllTest, TestMatrixFormatCSR2Ones) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_ones_in.txt"};
    // Create your matrix format here.
    CSR2 csr2_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr2_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatCSR2Zeros) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_zeros_in.txt"};
    // Create your matrix format here.
    CSR2 csr2_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr2_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatCSR2Random) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_random_in.txt"};
    // Create your matrix format here.
    CSR2 csr2_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr2_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

// ---------------------- CSR -------------------------

TEST_F(AllTest, TestMatrixFormatCSROnes) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_ones_in.txt"};
    // Create your matrix format here.
    CSR csr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatCSRZeros) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_zeros_in.txt"};
    // Create your matrix format here.
    CSR csr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatCSRRandom) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_random_in.txt"};
    // Create your matrix format here.
    CSR csr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        csr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}


// ---------------------- BSR -------------------------

TEST_F(AllTest, TestMatrixFormatBSRBasic) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/test_basic_mat_in.txt",
        "../inputs/test_basic_vec_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSROnes) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_ones_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSRZeros) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_zeros_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSRRandom) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_random_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSRBasicCustomBlockSize) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/test_basic_mat_in.txt",
        "../inputs/test_basic_vec_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files, 2, 2);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSROnesCustomBlockSize) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_ones_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files, 2, 2);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    for (double c : expected_result) {
        std::cout << c << " ";
    }
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSRZerosCustomBlockSize) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_zeros_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatBSRRandomCustomBlockSize) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_random_in.txt"};
    // Create your matrix format here.
    BSR bsr_matrix(files, 2, 2);
    // Compute matrix format
    std::vector<double> expected_result =
        bsr_matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

// ---------------------- DDD Naive -------------------------

TEST_F(AllTest, TestMatrixFormatDDDNaiveOnes) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_ones_in.txt"};;
    // Create your matrix format here.
    DDDNaive matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatDDDNaiveZeros) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_zeros_in.txt"};;
    // Create your matrix format here.
    DDDNaive matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

TEST_F(AllTest, TestMatrixFormatDDDNaiveRandom) {
    // Create base implementation from matrix files.
    AbstractMatrixFormat::MatrixFiles files = {
        "../inputs/Materials Problem/arc130.mtx",
        "../inputs/Materials Problem/test_vec_130_random_in.txt"};
    // Create your matrix format here.
    DDDNaive matrix(files);
    // Compute matrix format
    std::vector<double> expected_result =
        matrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}
////////////////elllpack

TEST_F(AllTest, TestMatrixFormatELL) {
    // Create matrix files
    AbstractMatrixFormat::MatrixFiles files = {
        "/workspaces/ellpack-spmv/inputs/test_basic_mat_in.txt",
        "/workspaces/ellpack-spmv/inputs/test_basic_vec_in.txt"};
    std::string matrixFile = "/workspaces/ellpack-spmv/inputs/test_basic_mat_in.txt";
    std::string vectorFile = "/workspaces/ellpack-spmv/inputs/test_basic_vec_in.txt";

    // Create your matrix format here.
    ELL ellMatrix(matrixFile, vectorFile);
    std::vector<double> expected_result = ellMatrix.PerformMultiplicationDouble();
    CompareWithBase(files, expected_result);
}

// COPY THIS FOR YOUR VERSION
// TEST_F(AllTest, TestMatrixFormatBasicMatrixFormat) {
//     // Create base implementation from matrix files.
//     AbstractMatrixFormat::MatrixFiles files = {
//         "../inputs/test_basic_mat_in.txt",
//         "../inputs/test_basic_vec_in.txt"};
//     // Create your matrix format here.
//     // Compute matrix format
//     std::vector<double> expected_result =
//     MATRIX.PerformMultiplicationDouble() CompareWithBase(files,
//     expected_result);
// }

}  // namespace format
