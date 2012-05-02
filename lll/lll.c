#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"

/* swap row i and j */
void fmpq_mat_swap_rows(fmpq_mat_t mat, long i, long j)
{	fmpq_t* temp;
	temp = mat->rows[i];
	mat->rows[i] = mat->rows[j];
	mat->rows[j] = temp;
}

/* compute inner product of length-n vectors v1 and v2 */
void inner_product(fmpz_t res, fmpz_t* v1, fmpz_t* v2, long n)
{	long i;
	fmpz_zero(res);
	for(i=0; i<n; i++)
		fmpz_addmul(res, v1[i], v2[i]);
}

/* compute GSO of mat */
void compute_gso(fmpq_mat_t gso, const fmpz_mat_t mat)
{	fmpz_t inner;
	fmpz_t one;
	fmpq_mat_t gram;
	long i, j, k, n;
	n = gso->r;
	fmpz_init(inner);
	fmpz_init_set_ui(one, 1);
	fmpq_mat_init(gram, n, n);

	/* Cholesky factorization */
	for(i=0; i<n; i++)
	{	for(j=0; j<=i; j++)
		{	inner_product(inner, mat->rows[i], mat->rows[j], n);
			fmpq_set_fmpz_frac(fmpz_mat_entry(gram, i, j), inner, one);
			for(k=0; k<j; k++)
				fmpq_submul(fmpq_mat_entry(gram, i, j), fmpq_mat_entry(gso, j, k), fmpq_mat_entry(gram, i, k));
			if(i==j)
				/* squared GS norms stored on diagonal */
				fmpq_set(fmpq_mat_entry(gso, i, j), fmpq_mat_entry(gram, i, j));
			else
				fmpq_div(fmpq_mat_entry(gso, i, j), fmpq_mat_entry(gram, i, j), fmpq_mat_entry(gram, j, j));
		}
	}

	fmpz_clear(inner);
	fmpz_clear(one);
	fmpq_mat_clear(gram);
}

/* size reduce the ith row against the previous rows */
void size_reduce(fmpz_mat_t mat, fmpq_mat_t gso, long i)
{	long j, k;
	fmpz_t round, one;
	fmpq_t rat, half;
	fmpz_init(round);
	fmpz_init_set_ui(one, 1);
	fmpq_init(rat);
	fmpq_init(half);
	fmpq_set_si(half, 1, 2);

	/* size reduce row i against row j */
	for(j=i-1; j>=0; j--)
	{	/* determine reduction coefficient */
		fmpq_set(rat, fmpq_mat_entry(gso, i, j));
		fmpq_sub(rat, rat, half);
		fmpz_cdiv_q(round, fmpq_numref(rat), fmpq_denref(rat));
		/* perform reduction */
		for(k=0; k<mat->c; k++)
			fmpz_submul(fmpz_mat_entry(mat, i, k), round, fmpz_mat_entry(mat, j, k));
		/* update GSO */
		fmpq_set_fmpz_frac(rat, round, one);
		for(k=0; k<j; k++)
			fmpq_submul(fmpq_mat_entry(gso, i, k), rat, fmpq_mat_entry(gso, j, k));
		fmpq_sub(fmpq_mat_entry(gso, i, j), fmpq_mat_entry(gso, i, j), rat);
	}

	fmpz_clear(round);
	fmpz_clear(one);
	fmpq_clear(rat);
	fmpq_clear(half);
}

/* run a naive implementation of LLL on matrix A */
void fmpz_mat_lll(fmpz_mat_t B, const fmpz_mat_t A)
{	long i, j, n;
	fmpq_mat_t gso;
	fmpq_t lovasz, threequarters, temp;
	fmpq_init(lovasz);
	fmpq_init(threequarters);
	fmpq_init(temp);
	fmpq_set_si(threequarters, 3, 4);
	n = A->r;
	fmpz_mat_set(B, A);
	fmpq_mat_init(gso, n, n);
	compute_gso(gso, B);

	/* index i counts the current number of LLL-reduced rows */
	for(i=1; i<n; i++)
	{	size_reduce(B, gso, i);
		/* check Lovasz condition */
		fmpq_mul(lovasz, fmpq_mat_entry(gso, i, i-1), fmpq_mat_entry(gso, i, i-1));
		fmpq_mul(lovasz, lovasz, fmpq_mat_entry(gso, i-1, i-1));
		fmpq_add(lovasz, lovasz, fmpq_mat_entry(gso, i, i));
		fmpq_mul(temp, threequarters, fmpq_mat_entry(gso, i-1, i-1));
		fmpq_sub(temp, lovasz, temp);
		if(fmpq_sgn(temp)<0)
		{	/* swap row i-1 and row i */
			fmpz_mat_swap_rows(B, 0, i-1, i);
			/* update GSO */
			fmpq_mat_swap_rows(gso, i-1, i);
			for(j=i-1; j<n; j++)
			{	fmpq_set(temp, fmpq_mat_entry(gso, j, i-1));
				fmpq_set(fmpq_mat_entry(gso, j, i-1), fmpq_mat_entry(gso, j, i));
				fmpq_set(fmpq_mat_entry(gso, j, i), temp);
			}
			fmpq_set(temp, fmpq_mat_entry(gso, i-1, i));
			for(j=i+1; j<n; j++)
				fmpq_submul(fmpq_mat_entry(gso, j, i), temp, fmpq_mat_entry(gso, j, i-1));
			fmpq_div(temp, fmpq_mat_entry(gso, i, i), lovasz);
			fmpq_mul(fmpq_mat_entry(gso, i, i), temp, fmpq_mat_entry(gso, i-1, i-1));
			fmpq_set(fmpq_mat_entry(gso, i-1, i-1), lovasz);
			fmpq_mul(temp, temp, fmpq_mat_entry(gso, i-1, i));
			fmpq_set_si(fmpq_mat_entry(gso, i-1, i), 0, 1);
			fmpq_set(fmpq_mat_entry(gso, i, i-1), temp);
			for(j=i+1; j<n; j++)
				fmpq_addmul(fmpq_mat_entry(gso, j, i-1), temp, fmpq_mat_entry(gso, j, i));
			/* decrement LLL-reduced index */
			i = FLINT_MAX(0, i-2);
		}
	}

	fmpq_clear(lovasz);
	fmpq_clear(threequarters);
	fmpq_clear(temp);
	fmpq_mat_clear(gso);
}
