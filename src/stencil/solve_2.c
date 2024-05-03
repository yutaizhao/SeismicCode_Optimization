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
    
#pragma omp parallel
    {
        f64 denom[STENCIL_ORDER];
        for (usz o = 0; o < STENCIL_ORDER; ++o) {
            denom[o] = 1/pow(17.0,(f64)o+1);
        }
#pragma omp for
            for (usz i = STENCIL_ORDER; i < dim_x - STENCIL_ORDER; ++i) {
                for (usz j = STENCIL_ORDER; j < dim_y - STENCIL_ORDER; ++j) {
                    for (usz k = STENCIL_ORDER; k < dim_z - STENCIL_ORDER; ++k) {
                        
                        C->cells_value[i][j][k] = A->cells_value[i][j][k] * B->cells_value[i][j][k];
                        
                        for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                            
                            C->cells_value[i][j][k]  += ((A->cells_value[i + o][j][k] * B->cells_value[i + o][j][k] ) +
                                                         (A->cells_value[i - o][j][k] * B->cells_value[i - o][j][k] ) +
                                                         (A->cells_value[i][j + o][k]  * B->cells_value[i][j + o][k] ) +
                                                         (A->cells_value[i][j - o][k]  * B->cells_value[i][j - o][k] ) +
                                                         (A->cells_value[i][j][k + o]  * B->cells_value[i][j][k + o] ) +
                                                         (A->cells_value[i][j][k - o]  * B->cells_value[i][j][k - o] )) * denom[o-1];
                            
                        }
                    }
                }
            }
    }
    
    mesh_copy_core(A, C);
}
