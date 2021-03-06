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

    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2014 William Hart

******************************************************************************/

#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "fmpz_vec.h"

#define __set(B, lenB, A, lenA)      \
do {                                 \
    _fmpz_vec_set((B), (A), (lenA)); \
    (lenB) = (lenA);                 \
} while (0)

#define __rem(R, lenR, A, lenA, B, lenB)                              \
do {                                                                  \
    if ((lenA) >= (lenB))                                             \
    {                                                                 \
        fmpz_invmod(invB, B + lenB - 1, mod);                         \
        _fmpz_mod_poly_rem((R), (A), (lenA), (B), (lenB), invB, mod); \
        (lenR) = (lenB) - 1;                                          \
        FMPZ_VEC_NORM((R), (lenR));                                   \
    }                                                                 \
    else                                                              \
    {                                                                 \
        _fmpz_vec_set((R), (A), (lenA));                              \
        (lenR) = (lenA);                                              \
    }                                                                 \
} while (0)

/*
    XXX: Incidentally, this implementation currently supports aliasing.  
    But since this may change in the future, no function other than 
    fmpz_mod_poly_gcd_hgcd() should rely on this.
 */

slong _fmpz_mod_poly_gcd_hgcd(fmpz *G, const fmpz *A, slong lenA, 
                                   const fmpz *B, slong lenB, const fmpz_t mod)
{
    fmpz *J = _fmpz_vec_init(2 * lenB);
    fmpz *R = J + lenB;

    fmpz_t invB;

    slong lenG, lenJ, lenR;

    fmpz_init(invB);
    
    __rem(R, lenR, A, lenA, B, lenB);

    if (lenR == 0)
    {
        __set(G, lenG, B, lenB);
    }
    else
    {
        _fmpz_mod_poly_hgcd(NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB, R, lenR, mod);

        while (lenJ != 0)
        {
            __rem(R, lenR, G, lenG, J, lenJ);

            if (lenR == 0)
            {
                __set(G, lenG, J, lenJ);
                break;
            }
            if (lenJ < FMPZ_MOD_POLY_GCD_CUTOFF)
            {
                fmpz_invmod(invB, R + lenR - 1, mod);
                lenG = _fmpz_mod_poly_gcd_euclidean(G, J, lenJ, R, lenR, invB, mod);
                break;
            }

            _fmpz_mod_poly_hgcd(NULL, NULL, G, &(lenG), J, &(lenJ), J, lenJ, R, lenR, mod);
        }
    }

    _fmpz_vec_clear(J, 2 * lenB);
    fmpz_clear(invB);

    return lenG;
}

void fmpz_mod_poly_gcd_hgcd(fmpz_mod_poly_t G, 
                             const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_gcd_hgcd(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        fmpz_mod_poly_t tG;
        fmpz *g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            fmpz_mod_poly_zero(G);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            fmpz_mod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                fmpz_mod_poly_init2(tG, &A->p, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _fmpz_mod_poly_gcd_hgcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, &A->p);

            if (G == A || G == B)
            {
                fmpz_mod_poly_swap(tG, G);
                fmpz_mod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                fmpz_one(G->coeffs + 0);
            else
                fmpz_mod_poly_make_monic(G, G);
        }
    }
}

#undef __set
#undef __rem

