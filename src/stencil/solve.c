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
#pragma omp parallel for
    for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += block) {
        for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += blocj) {
            for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += bloci) {
                for (usz k = kk; k < min(kk + block, dim_z - STENCIL_ORDER); k++) {
                    for (usz j = jj; j < min(jj + blocj, dim_y - STENCIL_ORDER); j++) {
                        for (usz i = ii; i < min(ii + bloci, dim_x - STENCIL_ORDER); i++) {
                            C->cells_value[i][j][k] = A->cells_value[i][j][k] * B->cells_value[i][j][k];
                            
                            for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                                f64 denom = 1/pow(17.0,(f64)o);
                                C->cells_value[i][j][k]  += ((A->cells_value[i + o][j][k] * B->cells_value[i + o][j][k] ) +
                                                             (A->cells_value[i - o][j][k] * B->cells_value[i - o][j][k] ) +
                                                             (A->cells_value[i][j + o][k]  * B->cells_value[i][j + o][k] ) +
                                                             (A->cells_value[i][j - o][k]  * B->cells_value[i][j - o][k] ) +
                                                             (A->cells_value[i][j][k + o]  * B->cells_value[i][j][k + o] ) +
                                                             (A->cells_value[i][j][k - o]  * B->cells_value[i][j][k - o] )) * denom;
                            }
                        }
                    }
                }
            }
        }
    }
    
    mesh_copy_core(A, C);
}
