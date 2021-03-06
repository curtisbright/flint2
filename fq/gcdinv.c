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
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq.h"

void
fq_gcdinv(fq_t rop, fq_t inv, const fq_t op, const fq_ctx_t ctx)
{
    fmpz *g, *s;
    slong lenG;
    
    if (fq_is_zero(op, ctx))
    {
        fq_zero(rop, ctx);
        return;
    }
    
    if (rop == op)
    {
        g = _fmpz_vec_init(op->length);
    }
    else
    {
        fmpz_poly_fit_length(rop, op->length);
        g = rop->coeffs;
    }

    if (inv == op)
    {
        s = _fmpz_vec_init(ctx->modulus->length - 1);
    }
    else
    {
        fmpz_poly_fit_length(inv, ctx->modulus->length - 1);
        s = inv->coeffs;
    }

    lenG = _fmpz_mod_poly_gcdinv(g, s, op->coeffs, op->length,
                                 ctx->modulus->coeffs, ctx->modulus->length,
                                 fq_ctx_prime(ctx));

    if (rop == op)
    {
        _fmpz_vec_clear(rop->coeffs, rop->alloc);
        rop->coeffs = g;
        rop->alloc = op->length;
    }
    if (inv == op)
    {
        _fmpz_vec_clear(inv->coeffs, inv->alloc);
        inv->coeffs = s;
        inv->alloc = ctx->modulus->length - 1;
    }

    _fmpz_poly_set_length(rop, lenG);
    _fmpz_poly_set_length(inv, ctx->modulus->length - lenG);
    _fmpz_poly_normalise(inv);
    
    if (!fmpz_is_one(fmpz_poly_lead(rop)))
    {
        fmpz_t invG;
        fmpz_init(invG);
        fmpz_invmod(invG, fmpz_poly_lead(rop), fq_ctx_prime(ctx));
        _fmpz_mod_poly_scalar_mul_fmpz(rop->coeffs, rop->coeffs, rop->length,
                                       invG, fq_ctx_prime(ctx));
        _fmpz_mod_poly_scalar_mul_fmpz(inv->coeffs, inv->coeffs, inv->length,
                                       invG, fq_ctx_prime(ctx));
        fmpz_clear(invG);
    }

}
