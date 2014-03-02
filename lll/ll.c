#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/* Compute Gram Matix G of vectors in B */
void GramMatrix(int d, mpz_t** G, mpz_t** B)
{	int i, j, k;
	for(i=0; i<d; i++)
		for(j=0; j<=i; j++)
			for(k=0; k<d; k++)
				mpz_addmul(G[i][j], B[i][k], B[j][k]);
}


/* Cholesky factorization: Calculate r, mu, s for GSO row i */
void CFA(int d, int i, mpz_t** G, mpf_t** r, mpf_t** mu, mpf_t* s)
{	int j, k;
	mpf_t temp;

	mpf_init(temp);

	/* Compute r[i][j] and mu[i][j] for j<=i */
	for(j=0; j<=i; j++)
	{	mpf_set_z(r[i][j], G[i][j]);
		for(k=0; k<j; k++)
		{	mpf_mul(temp, mu[j][k], r[i][k]);
			mpf_sub(r[i][j], r[i][j], temp);
		}
		if(j==i)
			mpf_set_d(mu[i][j], 1);
		else
			mpf_div(mu[i][j], r[i][j], r[j][j]);
	}
	/* Compute s[j] for j<=i */
	mpf_set_z(s[0], G[i][i]);
	for(j=1; j<=i; j++)
	{	mpf_mul(temp, mu[i][j-1], r[i][j-1]);
		mpf_sub(s[j], s[j-1], temp);
	}
	mpf_set(r[i][i], s[i]);

	mpf_clear(temp);
}

/* Size-reduce vector b_k of basis B */
void reduce(int d, int k, mpz_t** B, mpz_t** G, mpz_t** U, mpf_t** r, mpf_t** mu, mpf_t* s)
{	int i, j, loop;
	mpf_t eta, temp;
	mpz_t X[k], tempz;

	/* Initialization */
	mpf_init(temp);
	mpf_init_set_d(eta, (double)513/1024);
	mpz_init(tempz);
	for(i=0; i<k; i++)
		mpz_init(X[i]);

	do
	{	/* Compute GSO */
		CFA(d, k, G, r, mu, s);

		/* Pre-reduction computations */
		for(i=k-1; i>=0; i--)
		{	/* Compute X[i] */
			mpf_abs(temp, mu[k][i]);
			if(mpf_cmp(temp, eta)>=0)
			{	mpf_set_d(temp, (double)1/2);
				mpf_add(temp, temp, mu[k][i]);
				mpf_floor(temp, temp);
				mpz_set_f(X[i], temp);
			}
			else
				mpz_set_ui(X[i], 0);

			/* Update GSO */
			for(j=0; j<i; j++)
			{	mpf_set_z(temp, X[i]);
				mpf_mul(temp, temp, mu[i][j]);
				mpf_sub(mu[k][j], mu[k][j], temp);
			}
		}

		/* Size reduction of b_k for i<k */
		for(i=0; i<k; i++)
			for(j=0; j<d; j++)
				mpz_submul(B[k][j], X[i], B[i][j]);

		for(i=0; i<k; i++)
			for(j=0; j<d; j++)
				mpz_submul(U[k][j], X[i], U[i][j]);

		/* Update Gram matrix entry G[k][k] */
		for(j=0; j<k; j++)
		{	mpz_pow_ui(tempz, X[j], 2);
			mpz_addmul(G[k][k], tempz, G[j][j]);

			mpz_mul(tempz, X[j], G[k][j]);
			mpz_mul_2exp(tempz, tempz, 1);
			mpz_sub(G[k][k], G[k][k], tempz);

			for(i=0; i<j; i++)
			{	mpz_mul(tempz, X[i], X[j]);
				mpz_mul(tempz, tempz, G[j][i]);
				mpz_mul_2exp(tempz, tempz, 1);
				mpz_add(G[k][k], G[k][k], tempz);
			}
		}
		/* Update Gram matrix entries G[i][k] for i!=k */
		for(i=0; i<d; i++)
		{	if(i<k)
			{	for(j=0; j<=i; j++)
					mpz_submul(G[k][i], X[j], G[i][j]);
				for(j=i+1; j<k; j++)
					mpz_submul(G[k][i], X[j], G[j][i]);
			}
			else if(i>k)
			{	for(j=0; j<k; j++)
					mpz_submul(G[i][k], X[j], G[i][j]);
			}
		}

		/* Loop until all X[i] are zero */
		loop = 0;
		for(i=0; i<k; i++)
			if(mpz_cmp_ui(X[i], 0)!=0)
				loop = 1;
	} while(loop);

	/* Clean-up */
	mpf_clear(eta);
	mpf_clear(temp);
	mpz_clear(tempz);
	for(i=0; i<k-1; i++)
		mpz_clear(X[i]);
}

/* Insert vector b_kp before b_k and update Gram matrix */
void insert(int d, int k, int kp, mpz_t** B, mpz_t** G, mpz_t** U)
{	int i, j;
	mpz_t* temp;

	/* Swaps for vector insertion */
	temp = B[kp];
	for(j=kp; j>k; j--)
		B[j] = B[j-1];
	B[k] = temp;

	temp = U[kp];
	for(j=kp; j>k; j--)
		U[j] = U[j-1];
	U[k] = temp;

	/* Swaps to update Gram matrix */
	for(j=kp; j>k; j--)
	{	for(i=kp; i<d; i++)
			mpz_swap(G[i][j], G[i][j-1]);
		for(i=0; i<k; i++)
			mpz_swap(G[j][i], G[j-1][i]);
	}
	for(j=kp; j>k; j--)
		for(i=j; i>k; i--)
			mpz_swap(G[j][i], G[j-1][i-1]);
	for(j=0; 2*j<kp-k; j++)
		mpz_swap(G[k+j][k], G[kp-j][k]);
}

/* Run the LL Algorithm with factor pair (3/4, 1/2+2^(-9)) on the matrix A */
void fmpz_mat_ll(fmpz_mat_t out, const fmpz_mat_t A)
{	int d, i, j, k, kp;
	mpz_t** B, ** G, ** U;
	mpf_t** r, ** mu, * s;
	mpf_t delta, temp;

	/* Get dimension */
	d = A->r;

	/* Memory allocation */
	B = calloc(d, sizeof(mpz_t*));
	G = calloc(d, sizeof(mpz_t*));
	U = calloc(d, sizeof(mpz_t*));
	r = calloc(d, sizeof(mpf_t*));
	mu = calloc(d, sizeof(mpf_t*));
	for(i=0; i<d; i++)
	{	B[i] = calloc(d, sizeof(mpz_t));
		U[i] = calloc(d, sizeof(mpz_t));
		G[i] = calloc(i+1, sizeof(mpz_t));
		r[i] = calloc(i+1, sizeof(mpf_t));
		mu[i] = calloc(i+1, sizeof(mpf_t));
	}
	s = calloc(d, sizeof(mpf_t));

	/* Initialization */
	mpf_init(temp);
	mpf_init_set_d(delta, (double)7/8);
	for(i=0; i<d; i++)
	{	for(j=0; j<d; j++)
		{	mpz_init(B[i][j]);
			mpz_init(U[i][j]);
			if(j<=i)
			{	mpz_init(G[i][j]);
				mpf_init(r[i][j]);
				mpf_init(mu[i][j]);
			}
		}
		mpf_init(s[i]);
	}

	/* Set basis vectors B */
	for(i=0; i<d; i++)
		for(j=0; j<d; j++)
			fmpz_get_mpz(B[i][j], fmpz_mat_entry(A, i, j));

	/* Store the row operations in U */
	for(i=0; i<d; i++)
		mpz_set_ui(U[i][i], 1);

	/* Compute full Gram Matrix G */
	GramMatrix(d, G, B);

	mpf_set_z(r[0][0], G[0][0]);
	/* The first k-1 vectors in B are LLL-reduced as loop starts */
	for(k=1; k<d; k++)
	{	/* Size-reduce vector b_k */
		reduce(d, k, B, G, U, r, mu, s);

		/* Check consecutive Lovasz conditions */
		for(kp=k; k>0; k--)
		{	mpf_mul(temp, delta, r[k-1][k-1]);
			if(mpf_cmp(temp, s[k-1])<0)
				break;
		}

		/* Update GSO */
		for(i=0; i<k; i++)
		{	mpf_set(mu[k][i], mu[kp][i]);
			mpf_set(r[k][i], r[kp][i]);
		}
		mpf_set(r[k][k], s[k]);

		/* Insert vector b_kp before b_k */
		insert(d, k, kp, B, G, U);
	}

	/* Set output basis vectors */
	for(i=0; i<d; i++)
		for(j=0; j<d; j++)
			fmpz_set_mpz(fmpz_mat_entry(out, i, j), B[i][j]);

	/* Clean-up */
	for(i=0; i<d; i++)
	{	free(B[i]);
		free(U[i]);
		free(G[i]);
		free(r[i]);
		free(mu[i]);
	}
	free(B);
	free(G);
	free(U);
	free(r);
	free(mu);
	free(s);
}
