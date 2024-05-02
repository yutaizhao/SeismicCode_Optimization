#include "stencil/solve.h"
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
    usz bloci = 16 ;
    usz blocj = 2 ;
    usz block = 2 ;
    for (usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER; ii += bloci) {
        for (usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER; jj += blocj) {
            for (usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER; kk += block) {
                for (usz i = ii; i < min(ii + bloci, dim_x - STENCIL_ORDER); i++) {
                    for (usz j = jj; j < min(jj + blocj, dim_y - STENCIL_ORDER); j++) {
                        for (usz k = kk; k < min(kk + block, dim_z - STENCIL_ORDER); k++) {
			 C->cells[i][j][k].value = A->cells[i][j][k].value * B->cells[i][j][k].value;       
				for (usz o = 1; o <= STENCIL_ORDER; ++o) {
		    			f64 denom = 1/pow(17.0,(f64)o); 
                   			 C->cells[i][j][k].value += ((A->cells[i + o][j][k].value * B->cells[i + o][j][k].value) +
                                                (A->cells[i - o][j][k].value * B->cells[i - o][j][k].value) +
                                                (A->cells[i][j + o][k].value * B->cells[i][j + o][k].value) +
                                                (A->cells[i][j - o][k].value * B->cells[i][j - o][k].value) +
                                                (A->cells[i][j][k + o].value * B->cells[i][j][k + o].value) +
                                                (A->cells[i][j][k - o].value * B->cells[i][j][k - o].value)) * denom;
		
				 }
			}
		    }
		}
            }
        }
    }

    mesh_copy_core(A, C);
}
