#include <stdio.h>
#include <stdlib.h>

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "AbstractMatrixFormat.h"

namespace format {
class BasicMatrixFormat : public AbstractMatrixFormat {
   public:
    BasicMatrixFormat(const MatrixFiles& matrix_files)
        : AbstractMatrixFormat(matrix_files) {
        ReadMatrix();
    }

    std::vector<double> PerformMultiplicationDouble() {
        std::vector<double> result(rows, 0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                result[i] += matrix[i * columns + j] * vector[j];
            }
        }
        return result;
    }

    std::vector<float> PerformMultiplicationFloat() { return {}; };

    std::vector<int32_t> PerformMultiplicationInt32() { return {}; };

    Dimensions GetDimensions() const { return {rows, columns, vector_rows}; }

   private:
    void ReadMatrix() {
        MatrixFiles matrix_files = GetMatrixFiles();
        std::ifstream matrix_file(matrix_files.matrix_file);
        std::ifstream vector_file(matrix_files.vector_file);
        if (!matrix_file || !vector_file) {
            std::cerr << "Error opening files!" + matrix_files.matrix_file +
                             ", " + matrix_files.vector_file + "\n";
            std::cerr << "Error code: " << strerror(errno);
            exit(1);
        }
        // Expect first two elements to be number of rows and columns.
        matrix_file >> rows;
        matrix_file >> columns;
        matrix_file >> total_vals;
        matrix.assign(rows * columns, 0);
        int row, col;
        double value;
        for (int i = 0; i < total_vals; ++i) {
            matrix_file >> col >> row >> value;
            // 1-indexed.
            --col;
            --row;
            matrix[row * columns + col] = value;
        }
        // Expect first line to be number of rows.
        vector_file >> vector_rows;
        // Same number of columns as vector rows is expected.
        assert(columns == vector_rows);
        vector.assign(vector_rows, 0);
        for (int i = 0; i < vector_rows; ++i) {
            vector_file >> vector[i];
        }
        matrix_file.close();
        vector_file.close();
    }

    int32_t rows;
    int32_t columns;
    int32_t total_vals;
    int32_t vector_rows;
    std::vector<double> matrix;
    std::vector<double> vector;
};
}  // namespace format