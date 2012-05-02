#ifndef LLL
#define LLL

#include "fmpz_mat.h"
#include "fmpq_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

void fmpq_mat_swap_rows(fmpq_mat_t mat, long i, long j);
void inner_product(fmpz_t res, fmpz_t* v1, fmpz_t* v2, long n);
void compute_gso(fmpq_mat_t gso, const fmpz_mat_t mat);
void size_reduce(fmpz_mat_t mat, fmpq_mat_t gso, long i);
void fmpz_mat_lll(fmpz_mat_t B, const fmpz_mat_t A);

#ifdef __cplusplus
}
#endif

#endif
