#include "stencil/solve.h"

#include <omp.h>
#include <assert.h>
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b))

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);
    
    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;
    usz bloci = 8 ;
    usz blocj = 8 ;
    usz block = 8 ;

    f64 denom1= 1/17.0;
    f64 denom2= 1/(17.0*17.0);
    f64 denom3= 1/pow(17.0,3);
    f64 denom4= 1/pow(17.0,4);
    f64 denom5= 1/pow(17.0,5);
    f64 denom6= 1/pow(17.0,6);
    f64 denom7= 1/pow(17.0,7);
    f64 denom8= 1/pow(17.0,8);
    
#pragma omp parallel for
    for (usz i = STENCIL_ORDER; i < dim_x - STENCIL_ORDER; ++i) {
        for (usz j = STENCIL_ORDER; j < dim_y - STENCIL_ORDER; ++j) {
            for (usz k = STENCIL_ORDER; k < dim_z - STENCIL_ORDER; ++k) {
                    
                    C->cells_value[i][j][k]  = A->cells_value[i][j][k] * B->cells_value[i][j][k] +
                    
                    (A->cells_value[i + 1][j][k] * B->cells_value[i + 1][j][k] )* denom1 +
                      (A->cells_value[i - 1][j][k] * B->cells_value[i - 1][j][k] )* denom1+
                      (A->cells_value[i][j + 1][k]  * B->cells_value[i][j + 1][k] )* denom1 +
                      (A->cells_value[i][j - 1][k]  * B->cells_value[i][j - 1][k] )* denom1 +
                      (A->cells_value[i][j][k + 1]  * B->cells_value[i][j][k + 1] )* denom1 +
                      (A->cells_value[i][j][k - 1]  * B->cells_value[i][j][k - 1] )* denom1
                    +(A->cells_value[i + 2][j][k] * B->cells_value[i + 2][j][k] )* denom2 +
                      (A->cells_value[i - 2][j][k] * B->cells_value[i - 2][j][k] )* denom2 +
                      (A->cells_value[i][j + 2][k]  * B->cells_value[i][j + 2][k] )* denom2 +
                      (A->cells_value[i][j - 2][k]  * B->cells_value[i][j - 2][k] )* denom2 +
                      (A->cells_value[i][j][k + 2]  * B->cells_value[i][j][k + 2] )* denom2 +
                      (A->cells_value[i][j][k - 2]  * B->cells_value[i][j][k - 2] )* denom2
                    +(A->cells_value[i + 3][j][k] * B->cells_value[i + 3][j][k] )* denom3 +
                      (A->cells_value[i - 3][j][k] * B->cells_value[i - 3][j][k] )* denom3 +
                      (A->cells_value[i][j + 3][k]  * B->cells_value[i][j + 3][k] )* denom3 +
                      (A->cells_value[i][j - 3][k]  * B->cells_value[i][j - 3][k] )* denom3 +
                      (A->cells_value[i][j][k + 3]  * B->cells_value[i][j][k + 3] )* denom3 +
                      (A->cells_value[i][j][k - 3]  * B->cells_value[i][j][k - 3] )* denom3
                    +(A->cells_value[i + 4][j][k] * B->cells_value[i + 4][j][k] )* denom4 +
                      (A->cells_value[i - 4][j][k] * B->cells_value[i - 4][j][k] )* denom4 +
                      (A->cells_value[i][j + 4][k]  * B->cells_value[i][j + 4][k] )* denom4 +
                      (A->cells_value[i][j - 4][k]  * B->cells_value[i][j - 4][k] )* denom4 +
                      (A->cells_value[i][j][k + 4]  * B->cells_value[i][j][k + 4] )* denom4 +
                      (A->cells_value[i][j][k - 4]  * B->cells_value[i][j][k - 4] )* denom4
                    +(A->cells_value[i + 5][j][k] * B->cells_value[i + 5][j][k] )* denom5 +
                      (A->cells_value[i - 5][j][k] * B->cells_value[i - 5][j][k] )* denom5 +
                      (A->cells_value[i][j + 5][k]  * B->cells_value[i][j + 5][k] )* denom5 +
                      (A->cells_value[i][j - 5][k]  * B->cells_value[i][j - 5][k] )* denom5 +
                      (A->cells_value[i][j][k + 5]  * B->cells_value[i][j][k + 5] )* denom5 +
                      (A->cells_value[i][j][k - 5]  * B->cells_value[i][j][k - 5] )* denom5
                    +(A->cells_value[i + 6][j][k] * B->cells_value[i + 6][j][k] )* denom6 +
                      (A->cells_value[i - 6][j][k] * B->cells_value[i - 6][j][k] )* denom6 +
                      (A->cells_value[i][j + 6][k]  * B->cells_value[i][j + 6][k] )* denom6 +
                      (A->cells_value[i][j - 6][k]  * B->cells_value[i][j - 6][k] )* denom6 +
                      (A->cells_value[i][j][k + 6]  * B->cells_value[i][j][k + 6] )* denom6 +
                      (A->cells_value[i][j][k - 6]  * B->cells_value[i][j][k - 6] )
                    +(A->cells_value[i + 7][j][k] * B->cells_value[i + 7][j][k] )* denom7 +
                      (A->cells_value[i - 7][j][k] * B->cells_value[i - 7][j][k] )* denom7 +
                      (A->cells_value[i][j + 7][k]  * B->cells_value[i][j + 7][k] )* denom7 +
                      (A->cells_value[i][j - 7][k]  * B->cells_value[i][j - 7][k] )* denom7 +
                      (A->cells_value[i][j][k + 7]  * B->cells_value[i][j][k + 7] )* denom7 +
                      (A->cells_value[i][j][k - 7]  * B->cells_value[i][j][k - 7] )* denom7
                    +(A->cells_value[i + 8][j][k] * B->cells_value[i + 8][j][k] )* denom8 +
                      (A->cells_value[i - 8][j][k] * B->cells_value[i - 8][j][k] )* denom8 +
                      (A->cells_value[i][j + 8][k]  * B->cells_value[i][j + 8][k] )* denom8 +
                      (A->cells_value[i][j - 8][k]  * B->cells_value[i][j - 8][k] )* denom8 +
                      (A->cells_value[i][j][k + 8]  * B->cells_value[i][j][k + 8] )* denom8 +
                      (A->cells_value[i][j][k - 8]  * B->cells_value[i][j][k - 8] )* denom8
                    
                    ;
                    
            }
        }
    }
    
    mesh_copy_core(A, C);
}
