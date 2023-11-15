#ifndef ABSTRACT_MATRIX_FORMAT_H
#define ABSTRACT_MATRIX_FORMAT_H
#include <stdio.h>
#include <stdlib.h>

#include <string>

namespace format {
class AbstractMatrixFormat {
   public:
    struct Dimensions {
        int32_t mat_rows;
        int32_t mat_columns;
        int32_t vec_dim;
    };

    struct MatrixFiles {
        std::string matrix_file;
        std::string vector_file;
    };

    // Construction method. Each type of encoding should implement this method
    // to read the matrix and vector, and set the state with conversions as
    // needed.
    AbstractMatrixFormat(const MatrixFiles& matrix_files)
        : matrix_files_(matrix_files){};

    // Delete default constructor to not allow for mistakes.
    AbstractMatrixFormat() = delete;

    // Needed for virtualization purposes.
    virtual ~AbstractMatrixFormat() = default;

    // Returns a `Dimensions` reference that can be used to get information
    // about the size of the matrix and vector.
    Dimensions GetDimensions() const { return dimensions_; }

    // Returns a `MatrixFiles` reference that can be used to acces files.
    const MatrixFiles& GetMatrixFiles() const { return matrix_files_; }

    // Virtual method used to set the dimensions field(s). Subclasses may
    // override it.
    virtual void SetDimensions(const Dimensions& dimensions) {
        dimensions_ = dimensions;
    }

    // Virtual method used to set the matrix files. Subclasses may override it.
    virtual void SetMatrixFiles(const MatrixFiles& matrix_files) {
        matrix_files_ = matrix_files;
    }

    // Add more if needed. Types of computations are expected to be arrays.
    // TODO: Mention if this format doesn't work for specific types.
    virtual std::vector<double> PerformMultiplicationDouble() = 0;

    virtual std::vector<float> PerformMultiplicationFloat() = 0;

    virtual std::vector<int32_t> PerformMultiplicationInt32() = 0;

   private:
    MatrixFiles matrix_files_;
    Dimensions dimensions_;
};
}  // namespace format

#endif