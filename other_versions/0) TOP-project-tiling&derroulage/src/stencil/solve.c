#include "stencil/solve.h"

#include <assert.h>
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b))

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);
    
    f64 denom1 = 1/pow(17.0,1);
    f64 denom2 = 1/pow(17.0,2);
    f64 denom3 = 1/pow(17.0,3);
    f64 denom4 = 1/pow(17.0,4);
    f64 denom5 = 1/pow(17.0,5);
    f64 denom6 = 1/pow(17.0,6);
    f64 denom7 = 1/pow(17.0,7);
    f64 denom8 = 1/pow(17.0,8);
    
    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;
    usz bloci = 8 ;
    usz blocj = 8 ;
    usz block = 8 ;
    for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += block) {
        for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += blocj) {
            for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += bloci) {
                for (usz k = kk; k < min(kk + block, dim_z - STENCIL_ORDER); k++) {
                    for (usz j = jj; j < min(jj + blocj, dim_y - STENCIL_ORDER); j++) {
                        for (usz i = ii; i < min(ii + bloci, dim_x - STENCIL_ORDER); i++) {
                                
                            C->cells_value[i][j][k]  =
                            
                                A->cells_value[i][j][k] * B->cells_value[i][j][k]
                                
                                +((A->cells_value[i + 1][j][k] * B->cells_value[i + 1][j][k] ) +
                                  (A->cells_value[i - 1][j][k] * B->cells_value[i - 1][j][k] ) +
                                  (A->cells_value[i][j + 1][k]  * B->cells_value[i][j + 1][k] ) +
                                  (A->cells_value[i][j - 1][k]  * B->cells_value[i][j - 1][k] ) +
                                  (A->cells_value[i][j][k + 1]  * B->cells_value[i][j][k + 1] ) +
                                  (A->cells_value[i][j][k - 1]  * B->cells_value[i][j][k - 1] )) * denom1
                                
                                +((A->cells_value[i + 2][j][k] * B->cells_value[i + 2][j][k] ) +
                                  (A->cells_value[i - 2][j][k] * B->cells_value[i - 2][j][k] ) +
                                  (A->cells_value[i][j + 2][k]  * B->cells_value[i][j + 2][k] ) +
                                  (A->cells_value[i][j - 2][k]  * B->cells_value[i][j - 2][k] ) +
                                  (A->cells_value[i][j][k + 2]  * B->cells_value[i][j][k + 2] ) +
                                  (A->cells_value[i][j][k - 2]  * B->cells_value[i][j][k - 2] )) * denom2
                                
                                +((A->cells_value[i + 3][j][k] * B->cells_value[i + 3][j][k] ) +
                                  (A->cells_value[i - 3][j][k] * B->cells_value[i - 3][j][k] ) +
                                  (A->cells_value[i][j + 3][k]  * B->cells_value[i][j + 3][k] ) +
                                  (A->cells_value[i][j - 3][k]  * B->cells_value[i][j - 3][k] ) +
                                  (A->cells_value[i][j][k + 3]  * B->cells_value[i][j][k + 3] ) +
                                  (A->cells_value[i][j][k - 3]  * B->cells_value[i][j][k - 3] )) * denom3
                            
                                +((A->cells_value[i + 4][j][k] * B->cells_value[i + 4][j][k] ) +
                                  (A->cells_value[i - 4][j][k] * B->cells_value[i - 4][j][k] ) +
                                  (A->cells_value[i][j + 4][k]  * B->cells_value[i][j + 4][k] ) +
                                  (A->cells_value[i][j - 4][k]  * B->cells_value[i][j - 4][k] ) +
                                  (A->cells_value[i][j][k + 4]  * B->cells_value[i][j][k + 4] ) +
                                  (A->cells_value[i][j][k - 4]  * B->cells_value[i][j][k - 4] )) * denom4
                            
                                +((A->cells_value[i + 5][j][k] * B->cells_value[i + 5][j][k] ) +
                                  (A->cells_value[i - 5][j][k] * B->cells_value[i - 5][j][k] ) +
                                  (A->cells_value[i][j + 5][k]  * B->cells_value[i][j + 5][k] ) +
                                  (A->cells_value[i][j - 5][k]  * B->cells_value[i][j - 5][k] ) +
                                  (A->cells_value[i][j][k + 5]  * B->cells_value[i][j][k + 5] ) +
                                  (A->cells_value[i][j][k - 5]  * B->cells_value[i][j][k - 5] )) * denom5
                            
                                +((A->cells_value[i + 6][j][k] * B->cells_value[i + 6][j][k] ) +
                                  (A->cells_value[i - 6][j][k] * B->cells_value[i - 6][j][k] ) +
                                  (A->cells_value[i][j + 6][k]  * B->cells_value[i][j + 6][k] ) +
                                  (A->cells_value[i][j - 6][k]  * B->cells_value[i][j - 6][k] ) +
                                  (A->cells_value[i][j][k + 6]  * B->cells_value[i][j][k + 6] ) +
                                  (A->cells_value[i][j][k - 6]  * B->cells_value[i][j][k - 6] )) * denom6
                            
                                +((A->cells_value[i + 7][j][k] * B->cells_value[i + 7][j][k] ) +
                                  (A->cells_value[i - 7][j][k] * B->cells_value[i - 7][j][k] ) +
                                  (A->cells_value[i][j + 7][k]  * B->cells_value[i][j + 7][k] ) +
                                  (A->cells_value[i][j - 7][k]  * B->cells_value[i][j - 7][k] ) +
                                  (A->cells_value[i][j][k + 7]  * B->cells_value[i][j][k + 7] ) +
                                  (A->cells_value[i][j][k - 7]  * B->cells_value[i][j][k - 7] )) * denom7
                            
                                +((A->cells_value[i + 8][j][k] * B->cells_value[i + 8][j][k] ) +
                                  (A->cells_value[i - 8][j][k] * B->cells_value[i - 8][j][k] ) +
                                  (A->cells_value[i][j + 8][k]  * B->cells_value[i][j + 8][k] ) +
                                  (A->cells_value[i][j - 8][k]  * B->cells_value[i][j - 8][k] ) +
                                  (A->cells_value[i][j][k + 8]  * B->cells_value[i][j][k + 8] ) +
                                  (A->cells_value[i][j][k - 8]  * B->cells_value[i][j][k - 8] )) * denom8;
                            
                        }
                    }
                }
	    }
	}
}

                mesh_copy_core(A, C);
            }
