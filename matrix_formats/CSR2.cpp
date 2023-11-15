#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "AbstractMatrixFormat.h"
#define INT int
#define DOU double
#define AVX __m256d
#define SSE __m128i
namespace format {
class CSR2 : public AbstractMatrixFormat {
   public:
    CSR2(const MatrixFiles &matrix_files) : AbstractMatrixFormat(matrix_files) {
        Init(matrix_files);
    }

    ~CSR2() {
        free(row_ptr);
        free(vec_val);
        free(col_idx);
        free(mtx_val);
        free(par_set);
        free(Csr2_row_ptr);
        _mm_free(Csr2_mtx_val);
    }

    std::vector<float> PerformMultiplicationFloat() { return {}; };

    std::vector<int32_t> PerformMultiplicationInt32() { return {}; };

    std::vector<double> PerformMultiplicationDouble() {
        Csr2_mid_val = (AVX *)_mm_malloc(sizeof(AVX) * par_set[6], 64);
        mtx_ans.assign(par_set[3] - 1, 0);
        //---------------------------CSR2 SpMV---------------------------//
        Avx2_Madd(Csr2_mid_val, Csr2_mtx_val, col_idx, vec_val, par_set);
        Avx2_SegSum(Csr2_row_ptr, Csr2_mid_val, mtx_ans, par_set);
        _mm_free(Csr2_mid_val);
        return mtx_ans;
    }

   private:
    void Init(const MatrixFiles &matrix_files) {
        FILE *fp_mtx;
        FILE *fp_vec;

        fp_mtx =
            fopen(matrix_files.matrix_file.c_str(), "rb+");  // Read matrix file
        fp_vec =
            fopen(matrix_files.vector_file.c_str(), "rb+");  // Read vector file

        // Construct matrix and vector
        ReadFileConvertCsr(fp_mtx, fp_vec, row_ptr, col_idx, mtx_val, vec_val,
                           par_set);
        Avx2_CsrConvertCsr2(row_ptr, mtx_val, Csr2_row_ptr, Csr2_mtx_val,
                            par_set);
    }

    void ReadFileConvertCsr(FILE *fp_mtx, FILE *fp_vec, INT *&row_ptr,
                            INT *&col_idx, DOU *&mtx_val, DOU *&vec_val,
                            INT *&par_set) {
        INT row_num, col_num, nzz_num;
        fscanf(fp_mtx, "%d %d %d", &row_num, &col_num, &nzz_num);
        //-----------------------Read Vector-------------------------//
        INT vec_row;
        vec_row = col_num;
        // srand(time(NULL));
        fscanf(fp_vec, "%d", &vec_row);
        vec_val = (DOU *)aligned_alloc(64, sizeof(DOU) * vec_row);
        // vec_val = (DOU *)malloc(sizeof(DOU)*vec_row);
        memset(vec_val, 0, sizeof(DOU) * vec_row);
        for (INT i = 1; i <= vec_row; i++) {
            fscanf(fp_vec, "%lf", &vec_val[i]);
            // vec_val[i] = ((rand()%10)+1)*1.0;
        }
        //-----------------------Read Matrix-------------------------//
        par_set = (INT *)malloc(sizeof(INT) * 10);
        // INT row_num,col_num,nzz_num;
        // fscanf(fp_mtx,"%d %d %d",&row_num,&col_num,&nzz_num);
        // printf("(%d , %d) , nzz_num = %d\n", row_num, col_num, nzz_num);
        INT mtx_width;
        INT x, cnt_ite;
        x = INT(nzz_num * 1.0 / row_num + 0.5);
        if (x <= 16) {
            mtx_width = x;
        } else {
            if (x > 128) {
                cnt_ite = 0;
                while (x > 128) {
                    if (x & 1) {
                        x++;
                    }
                    x /= 2;
                    if (cnt_ite == 4) {
                        break;
                    }
                    cnt_ite++;
                }
                mtx_width = x;
            } else {
                if (x >= 36) {
                    cnt_ite = 0;
                    while (x >= 36) {
                        if (x & 1) {
                            x++;
                        }
                        x /= 2;
                        if (cnt_ite == 2) {
                            break;
                        }
                        cnt_ite++;
                    }
                    mtx_width = x;
                } else {
                    // printf("%d\n",x);
                    if (x & 1) {
                        if (x % 10 >= 5)
                            x = (x / 10 + 1) * 10;
                        else
                            x = x / 10 * 10;
                    }
                    x /= 2;
                    mtx_width = x;
                }
            }
        }
        // mtx_width = 8;
        // printf("mtx_width = %d \n", mtx_width);
        INT mtx_high = 4;
        INT str_max;
        str_max = nzz_num + row_num * (mtx_width - 1);
        row_ptr = (INT *)aligned_alloc(64, sizeof(INT) * (row_num + 1));
        // row_ptr = (INT *)malloc(sizeof(INT)*(row_num+1));
        col_idx = (INT *)aligned_alloc(64, sizeof(INT) * str_max);
        // col_idx = (INT *)malloc(sizeof(INT)*str_max);
        mtx_val = (DOU *)aligned_alloc(64, sizeof(DOU) * str_max);
        // mtx_val = (DOU *)malloc(sizeof(DOU)*str_max);
        memset(mtx_val, 0, sizeof(DOU) * str_max);
        memset(col_idx, 0, sizeof(INT) * str_max);
        INT v = -1, cnt = 0, cnt1;
        INT row, col;
        DOU nzz;
        for (INT i = 0; i < nzz_num; i++) {
            fscanf(fp_mtx, "%d %d %lf", &col, &row, &nzz);
            row--;
            while (v < row) {
                v++;
                if (cnt % mtx_width != 0) {
                    cnt = (cnt + mtx_width - 1) / mtx_width * mtx_width;
                }
                row_ptr[v] = cnt;
            }
            col_idx[cnt] = col;
            mtx_val[cnt] = nzz;
            cnt++;
        }
        while (v < row_num) {
            v++;
            if (cnt % mtx_width != 0) {
                cnt = (cnt + mtx_width - 1) / mtx_width * mtx_width;
            }
            row_ptr[v] = cnt;
        }
        if (cnt % mtx_width != 0) {
            cnt = (cnt + mtx_width - 1) / mtx_width * mtx_width;
        }
        v++;
        row_ptr[v] = cnt;
        if (cnt % (mtx_width * mtx_high) != 0) {
            cnt = (cnt + mtx_width * mtx_high - 1) / (mtx_width * mtx_high) *
                  (mtx_width * mtx_high);
        }
        par_set[0] = cnt;
        par_set[1] = mtx_width;
        par_set[2] = mtx_high;
        par_set[3] = row_num + 1;
        par_set[5] = par_set[0] / par_set[2];
        par_set[6] = par_set[5] / par_set[1];
        par_set[7] = nzz_num;
    }

    void Avx2_CsrConvertCsr2(INT *&row_ptr, DOU *&mtx_val, INT *&Csr2_row_ptr,
                             AVX *&Csr2_mtx_val, INT *par_set) {
        INT i, x, y, z, i1, t;

        const INT i_end = par_set[3];
        const INT i_end1 = par_set[0];
        const INT mtx_width = par_set[1];
        // const INT tile_size = mtx_width<<2;

        Csr2_row_ptr = (INT *)aligned_alloc(64, sizeof(INT) * i_end);
        Csr2_mtx_val = (AVX *)_mm_malloc(sizeof(AVX) * (i_end1 >> 2), 64);

#pragma omp parallel private(i)
        {
#pragma omp for schedule(static) nowait
            for (i = 0; i < i_end; i++) {
                Csr2_row_ptr[i] = row_ptr[i] / mtx_width;
            }
        }
// i_end1 = i_end1/(mtx_high*mtx_width);
#pragma omp parallel private(i, x, y, z, t)
        {
#pragma omp for schedule(static) nowait
            for (i = 0; i < i_end1; i++) {
                /*
                for(int j=0;j<mtx_high;j++)
                {
                        int x=
                        for(int z=0;z<mtx_width;z++)
                        {
                                x
                        }
                }*/
                t = i / mtx_width;
                z = t >> 2;
                z = z * mtx_width;
                x = z + i - t * mtx_width;
                y = (i - (z << 2)) / mtx_width;
                // y = (i%tile_size)/mtx_width;
                Csr2_mtx_val[x][y] = mtx_val[i];
            }
        }
    }

    void Avx2_Madd(AVX *&Csr2_mid_val, AVX *&Csr2_mtx_val, INT *col_idx,
                   DOU *&vec_val, INT *par_set) {
        const INT i_end = par_set[6];
        const INT mtx_width = par_set[1];
        const INT mtx_high = par_set[2];
        INT i, j, xx, yy, zz;
        const INT m1 = mtx_width << 1;
        const INT m2 = m1 + mtx_width;
        AVX Csr2_mid_val1;
        AVX Csr2_col_val;
#pragma omp parallel private(i, j, xx, yy, zz, Csr2_mid_val1, Csr2_col_val)
        {
#pragma omp for schedule(static) nowait
            for (i = 0; i < i_end; i++) {
                // Csr2_mid_val[i]=_mm256_setzero_pd();
                // Csr2_mid_val1=_mm256_setzero_pd();
                xx = i * mtx_width;
                yy = xx << 2;
                Csr2_col_val = _mm256_set_pd(
                    vec_val[col_idx[yy + m2]], vec_val[col_idx[yy + m1]],
                    vec_val[col_idx[yy + mtx_width]], vec_val[col_idx[yy]]);
                Csr2_mid_val1 = _mm256_mul_pd(Csr2_mtx_val[xx], Csr2_col_val);
                xx = (xx << 1) + xx;
                zz = yy + mtx_width;
                for (j = yy + 1; j < zz; j++) {
                    // INT JJ=yy+j;
                    Csr2_col_val = _mm256_set_pd(
                        vec_val[col_idx[j + m2]], vec_val[col_idx[j + m1]],
                        vec_val[col_idx[j + mtx_width]], vec_val[col_idx[j]]);
                    // Csr2_col_val=_mm256_set_pd(vec_val[col_idx[JJ+m2]],vec_val[col_idx[JJ+m1]],vec_val[col_idx[JJ+mtx_width]],vec_val[col_idx[JJ]]);
                    // Csr2_mid_val[i] =
                    // _mm256_fmadd_pd(Csr2_mtx_val[j-xx],Csr2_col_val,Csr2_mid_val[i]);
                    Csr2_mid_val1 = _mm256_fmadd_pd(
                        Csr2_mtx_val[j - xx], Csr2_col_val, Csr2_mid_val1);
                }
                Csr2_mid_val[i] = Csr2_mid_val1;
            }
        }
    }

    void Avx2_SegSum(INT *&Csr2_row_ptr, AVX *&Csr2_mid_val,
                     std::vector<double> &mtx_ans, INT *&par_set) {
        INT i, j;
        const INT num = par_set[3] - 1;
#pragma omp parallel private(i, j)
        {
#pragma omp for schedule(static) nowait
            for (i = 0; i < num; i++) {
                DOU s = 0;
                int t1, t2;
                t1 = Csr2_row_ptr[i];
                t2 = Csr2_row_ptr[i + 1];
                for (j = t1; j < t2; j++) {
                    int t = j >> 2;
                    s += Csr2_mid_val[t][j - (t << 2)];
                }
                mtx_ans[i] = s;
            }
        }
    }

    DOU *vec_val;
    INT *row_ptr;
    INT *col_idx;
    DOU *mtx_val;
    INT *par_set;
    std::vector<double> mtx_ans;
    INT *Csr2_row_ptr;
    AVX *Csr2_mtx_val;
    AVX *Csr2_mid_val;
};
}  // namespace format
