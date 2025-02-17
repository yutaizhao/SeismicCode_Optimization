#include "stencil/solve.h"

#include <assert.h>
#include <math.h>
#include <omp.h>

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);
    
    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;
#pragma omp parallel for
    for (usz k = STENCIL_ORDER; k < dim_z - STENCIL_ORDER; ++k) {
        for (usz j = STENCIL_ORDER; j < dim_y - STENCIL_ORDER; ++j) {
            for (usz i = STENCIL_ORDER; i < dim_x - STENCIL_ORDER; ++i) {
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
    
    mesh_copy_core(A, C);
}
