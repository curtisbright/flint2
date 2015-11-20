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

    Copyright (C) 2015 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d)
{
   slong n = M->r, i, j;
   ulong ** A = M->rows;

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         NMOD_ADDMUL(A[i][j], A[i][r], d, M->mod);
      
      for (j = r + 1; j < n; j++)
         NMOD_ADDMUL(A[i][j], A[i][r], d, M->mod); 
   }

   d = n_negmod(d, M->mod.n);

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         NMOD_ADDMUL(A[r][i], A[j][i], d, M->mod);

      for (j = r + 1; j < n; j++)
         NMOD_ADDMUL(A[r][i], A[j][i], d, M->mod);      
   }
}
