#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <vector>
#define INT int
#define DOU double
#define DEFAULT_W 3
#define DEFAULT_H 3

#include "AbstractMatrixFormat.h"


namespace format {
class BSR : public AbstractMatrixFormat {
   public:
    BSR(const MatrixFiles &matrix_files) : AbstractMatrixFormat(matrix_files) {
        block_width = DEFAULT_W;
        block_height = DEFAULT_H;
        Init(matrix_files);
    }

    BSR(const MatrixFiles &matrix_files, int w, int h) : AbstractMatrixFormat(matrix_files) {
        block_width = w;
        block_height = h;
        Init(matrix_files);
    } 

    ~BSR() {
        free(vec_val);
        free(mtx_val);
    }

    std::vector<double> PerformMultiplicationDouble() {
        mtx_ans.assign(row_num, 0);
        #pragma omp parallel
        {
            INT row_idx = 0;
            std::vector<DOU>private_mtx_ans(row_num, 0);
            INT t = omp_get_num_threads();
            INT id = omp_get_thread_num();
            INT window_size = col_idx.size() / t;
            INT leftover = col_idx.size() % t;
            if (id < leftover) {
                window_size++;
            }
            int start = window_size * id + std::min(leftover, id);
            for (int i = start; i < start + window_size; i++) {
                // TODO: BINARY SEARCH AND TO REMOVE THE FIRST CLAUSE IN THE OR
                while ((row_ptr[row_idx] < i && row_ptr[row_idx + 1] < i) || row_ptr[row_idx + 1] == i) {
                    row_idx++;
                }
                // origin = top left corner of the block
                INT origin_row = row_idx * block_height;
                INT origin_col = col_idx[i];
                INT val_idx = i * block_height * block_width;
                for (int r = 0; r < block_height; r++) {
                    for (int c = 0; c < block_width; c++) {
                        INT row = origin_row + r;
                        INT col = origin_col + c;
                        // Handles zero padding
                        if (row < row_num && col < col_num) {
                            private_mtx_ans[row] += bsr_val[val_idx + r * block_width + c] * vec_val[col];
                        }
                    }
                }
            }
            for (int i = 0; i < row_num; i++) {
                #pragma omp atomic
                mtx_ans[i] += private_mtx_ans[i];
            }
        }

        return mtx_ans;
    }
    
    std::vector<float> PerformMultiplicationFloat() { return {}; };

    std::vector<int32_t> PerformMultiplicationInt32() { return {}; };

    private:
     void ReadVector(const MatrixFiles &matrix_files) {
        fp_vec = fopen(matrix_files.vector_file.c_str(), "rb+");
        INT vec_row;
        fscanf(fp_vec, "%d", &vec_row);
        vec_val = (DOU *) malloc(sizeof(DOU) * vec_row);
        memset(vec_val, 0, sizeof(DOU) * vec_row);
        for (INT i = 0; i < vec_row; i++) {
            fscanf(fp_vec, "%lf", &vec_val[i]);
        }
        fclose(fp_vec);
     }
     
     void ReadMatrix(const MatrixFiles &matrix_files) {
        //------------------------Read Matrix-------------------- //
        fp_mtx = fopen(matrix_files.matrix_file.c_str(), "rb+");
        fscanf(fp_mtx, "%d %d %d", &row_num, &col_num, &nzz_num);
        mtx_val = (DOU *) malloc(sizeof(DOU) * row_num * col_num);
        memset(mtx_val, 0, sizeof(DOU) * row_num * col_num);
        INT last_block = -1;
        INT row, col;
        DOU nzz;
        for (int i = 0; i < nzz_num; i++) {
            fscanf(fp_mtx, "%d %d %lf", &col, &row, &nzz);
            row--;
            col--;
            mtx_val[getIndex(row, col)] = nzz;
        }
        fclose(fp_mtx);
     }

     int getIndex(int r, int c) {
        return (col_num * r) + c;
     }


     void Init(const MatrixFiles &matrix_files) {
        ReadVector(matrix_files);
        ReadMatrix(matrix_files);
        //-------------Convert Matrix to BSR------------------
        //Number of blocks per row and col
        INT block_rows_num = row_num / block_height;
        INT block_cols_num = col_num / block_width;
        // Zero padding for when the mtx dims are not divisible by the block size
        if ((block_rows_num % block_height)) block_rows_num++;
        if ((block_cols_num % block_width)) block_cols_num++;

        INT prev_row = -block_height;
        //Iterate over the blocks
        for (int i = 0; i < block_rows_num; i++) {
            for (int j = 0; j < block_cols_num; j++) {
                bool save_block = false;
                // Top left corner of each block
                INT origin_row = i * block_height;
                INT origin_col = j * block_width;
                for (int r = 0; r < block_height; r++) {
                    for (int c = 0; c < block_width; c++) {
                        if (save_block) {
                            if (origin_row + r >= row_num || origin_col + c >= col_num) {
                                bsr_val.push_back(0);
                            } else {
                                bsr_val.push_back(mtx_val[getIndex(origin_row + r, origin_col + c)]);
                            }
                        } else if (origin_row + r < row_num && origin_col + c < col_num 
                            && mtx_val[getIndex(origin_row + r, origin_col + c)]) {

                            save_block = true;
                            // Fill the missing zeros
                            bsr_val.insert(bsr_val.end(), block_width * r + c, 0);
                            bsr_val.push_back(mtx_val[getIndex(origin_row + r, origin_col + c)]);
                            while (prev_row != origin_row) {
                                prev_row += block_height;
                                row_ptr.push_back(col_idx.size());
                            }
                            col_idx.push_back(origin_col);
                        }
                    }
                }
            }
        }
        row_ptr.push_back(col_idx.size());
    }
    
    FILE *fp_mtx;
    FILE *fp_vec;
    FILE *fp_ans;
    DOU *vec_val;
    DOU *mtx_val;
    std::vector<INT> row_ptr;
    std::vector<INT> col_idx;
    std::vector<DOU> bsr_val;
    std::vector<DOU> mtx_ans;
    INT row_num;
    INT col_num;
    INT nzz_num;
    INT block_width;
    INT block_height;
};
} //namespace format
