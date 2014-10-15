/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright 2009 William Hart
    Copyright 2010,2011 Fredrik Johansson
    Copyright 2014 Abhinav Baid

******************************************************************************/

#include <fplll.h>
#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz_lll.h"
#include "fmpz.h"

using namespace std;
using namespace fplll;

int main(void)
{
    for (slong dim = 10; dim <= 500; dim += 10)
    {
		double delta = 0.75;
		double eta = 0.51;
		fmpz_lll_t fl;
		flint_rand_t rnd;
		fmpz_mat_t A;
		fmpz_mat_t B;
		FLINT_TEST_INIT(state);
		fmpz_mat_init(A, dim, dim);
		fmpz_mat_init(B, dim, dim);
		fmpz_lll_context_init(fl, delta, eta, 1, 0);
		fmpz_mat_randajtai(A, state, 1.25);
		fmpz_mat_set(B, A);

		char buf[L_tmpnam];
		tmpnam(buf);
		FILE* f = fopen(buf, "w");
		fmpz_mat_fprint_pretty(f, A);
		fclose(f);

		ZZ_mat<mpz_t> M(dim, dim);
		ifstream is(buf);
		is >> M;

		timeit_t t0;
		timeit_start(t0);
		fmpz_lll_wrapper(A, NULL, fl);
		timeit_stop(t0);

		timeit_t t1;
		timeit_start(t1);
		fmpz_lll(B, NULL, fl);
		timeit_stop(t1);

		timeit_t t2;
		timeit_start(t2);
		lllReduction(M, (delta+1)/2, (2*eta+1)/4, LM_WRAPPER);
		timeit_stop(t2);

		flint_printf("dim = %wd flint/ulll/fplll %wd %wd %wd (ms) fplll-ratio %0.2f %0.2f\n", dim, t0->wall, t1->wall, t2->wall, (double)t0->wall/t2->wall, (double)t1->wall/t2->wall);

		fmpz_mat_clear(A);
		fmpz_mat_clear(B);
		flint_randclear(state);
	}

	return 0;
}
