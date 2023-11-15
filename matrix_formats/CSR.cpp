#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#define INT int
#define DOU double

#include "AbstractMatrixFormat.h"

namespace format {
class CSR : public AbstractMatrixFormat {
   public:
    CSR(const MatrixFiles &matrix_files) : AbstractMatrixFormat(matrix_files) {
        Init(matrix_files);
    }

    ~CSR() {
        free(vec_val);
        free(row_ptr);
        free(col_idx);
        free(mtx_val);
        free(par_set);
    }

    std::vector<double> PerformMultiplicationDouble() {
        mtx_ans.assign(row_num, 0);
        for (INT i = 0; i < row_num; i++) {
            for (INT j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
                mtx_ans[i] += mtx_val[j] * vec_val[col_idx[j]];
            }
        }
        return mtx_ans;
    }

    std::vector<float> PerformMultiplicationFloat() { return {}; };

    std::vector<int32_t> PerformMultiplicationInt32() { return {}; };

   private:
    void Init(const MatrixFiles &matrix_files) {
        fp_mtx = fopen(matrix_files.matrix_file.c_str(), "rb+");
        fp_vec = fopen(matrix_files.vector_file.c_str(), "rb+");
        INT vec_row;
        fscanf(fp_vec, "%d", &vec_row);
        vec_val = (DOU *)malloc(sizeof(DOU) * vec_row);
        memset(vec_val, 0, sizeof(DOU) * vec_row);
        for (INT i = 0; i < vec_row; i++) {
            fscanf(fp_vec, "%lf", &vec_val[i]);
        }
        //-----------------------Read Matrix-------------------------//
        par_set = (INT *)malloc(sizeof(INT) * 10);
        fscanf(fp_mtx, "%d %d %d", &row_num, &col_num, &nzz_num);
        INT str_max = nzz_num;
        row_ptr = (INT *)malloc(sizeof(INT) * (row_num + 1));
        col_idx = (INT *)malloc(sizeof(INT) * str_max);
        mtx_val = (DOU *)malloc(sizeof(DOU) * str_max);
        memset(mtx_val, 0, sizeof(DOU) * str_max);
        memset(col_idx, 0, sizeof(INT) * str_max);
        INT v = -1, cnt = 0, cnt1;
        INT row, col;
        DOU nzz;
        for (INT i = 0; i < nzz_num; i++) {
            fscanf(fp_mtx, "%d %d %lf", &col, &row, &nzz);
            row--;
            col--;
            while (v < row) {
                v++;
                row_ptr[v] = cnt;
            }
            col_idx[cnt] = col;
            mtx_val[cnt] = nzz;
            cnt++;
        }
        while (v < row_num) {
            v++;
            row_ptr[v] = cnt;
        }
        row_ptr[row_num] = cnt;
    }

    FILE *fp_mtx;
    FILE *fp_vec;
    FILE *fp_ans;
    DOU *vec_val;
    INT *row_ptr;
    INT *col_idx;
    DOU *mtx_val;
    INT *par_set;
    std::vector<DOU> mtx_ans;
    INT row_num;
    INT col_num;
    INT nzz_num;
};
}  // namespace format
