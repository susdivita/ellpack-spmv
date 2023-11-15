#include <stdio.h>
#include <stdlib.h>

#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include "AbstractMatrixFormat.h"

 struct Dimensions {
        int32_t mat_rows;
        int32_t mat_columns;
        int32_t vec_dim;
    };

    struct MatrixFiles {
        std::string matrix_file;
        std::string vector_file;
    };

    int32_t rows;
    int32_t columns;
    int32_t total_vals;
    int32_t vector_rows;
    std::vector<double> matrix;
    std::vector<double> vector;

int main() {
    std::ifstream matrix_file("/workspaces/ellpack-spmv/inputs/test_basic_mat_in.txt");
        std::ifstream vector_file("/workspaces/ellpack-spmv/inputs/test_basic_vec_in.txt");
        if (!matrix_file || !vector_file) {
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
        vector.assign(vector_rows, 1);
        for (int i = 0; i < vector_rows; ++i) {
            vector_file >> vector[i];
        }
        matrix_file.close();
        //vector_file.close();
         std::vector<double> result(rows, 0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                result[i] += matrix[i * columns + j] * vector[j];
            }
        }
    for (double c : result) {
        std::cout << c << " ";
    }
}