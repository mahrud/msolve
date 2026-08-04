/* This file is part of msolve.
 *
 * msolve is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * msolve is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with msolve.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Mahrud Sayrafi */

/* Syzygies and free resolutions of graded modules.
 *
 * This header declares the grading-group layer and the module-order
 * vocabulary shared by the resolution engine.  See res_grading.c for the
 * arithmetic and res_frame.c for the Schreyer frame. */

#ifndef GB_RES_H
#define GB_RES_H

#include "data.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The module monomial orders (res_mord_t) live in data.h next to ht_t,
 * since they are a property of the hash table just like ht->mo. */

/* --------------------------------------------------------------------- *
 *  Degrees
 *
 *  Modules are graded by a finitely generated abelian group
 *
 *      A  =  Z^r (+) Z/t_1 (+) ... (+) Z/t_nt
 *
 *  which covers the standard grading (r = 1, nt = 0), multigradings by
 *  Z^r, the class group Cl(X) of a toric variety, and -- once an external
 *  group-arithmetic backend is attached to the vtable below -- character
 *  groups for equivariant resolutions.
 *
 *  A degree is a *struct*, never a bare array: its layout is private to
 *  the owning res_dgrp_t so that it can later grow a cached hash, a
 *  torsion handle, or a pointer into an external library without breaking
 *  any consumer.
 * --------------------------------------------------------------------- */

typedef struct res_deg_t res_deg_t;
struct res_deg_t
{
    int32_t *e;  /* [0,r) free part, [r,r+nt) torsion residues;
                  * interpretation belongs to the owning res_dgrp_t */
};

typedef struct res_dgrp_t res_dgrp_t;
struct res_dgrp_t
{
    len_t    r;      /* rank of the free part                            */
    len_t    nt;     /* number of torsion factors                        */
    len_t    len;    /* int32 slots per degree, = r + nt                 */
    len_t    nv;     /* number of ring variables                         */
    int32_t *tord;   /* orders of the torsion factors, length nt         */
    int32_t *heft;   /* length r; heft . deg(x_j) > 0 for every j        */
    int32_t *dmat;   /* len x nv, column major: column j is deg(x_j)     */
    deg_t   *vhdeg;  /* length nv, precomputed heft . deg(x_j)           */
    int      simple; /* 1 iff r == 1 && nt == 0 (the standard grading)   */

    /* Arithmetic vtable.  This is the seam an external group-arithmetic
     * library plugs into; the defaults in res_grading.c implement
     * Z^r (+) torsion directly. */
    void  (*add)(const res_dgrp_t *g, res_deg_t d,
                 const res_deg_t a, const res_deg_t b);
    void  (*sub)(const res_dgrp_t *g, res_deg_t d,
                 const res_deg_t a, const res_deg_t b);
    int   (*cmp)(const res_dgrp_t *g, const res_deg_t a, const res_deg_t b);
    hl_t  (*hash)(const res_dgrp_t *g, const res_deg_t a);
    deg_t (*heft_of)(const res_dgrp_t *g, const res_deg_t a);
};

/* Contiguous storage for many degrees of one group.  Degrees are handed
 * out as views into the pool, so a res_deg_t stays one pointer wide. */
typedef struct res_dpool_t res_dpool_t;
struct res_dpool_t
{
    const res_dgrp_t *grp;
    int32_t *data;
    hl_t     ld;   /* number of degrees in use   */
    hl_t     sz;   /* number of degrees allocated */
};

/* --- construction ---------------------------------------------------- */

/* degs is len x nv in column-major order: entries [j*len, (j+1)*len) give
 * deg(x_j), free part first then torsion residues.  heft has length r and
 * must satisfy heft . deg(x_j) > 0 for every variable j; pass NULL to use
 * the all-ones vector.  Returns NULL and reports on stderr if the heft is
 * not strictly positive on some variable. */
res_dgrp_t *res_dgrp_new(
        const len_t r,
        const len_t nt,
        const int32_t *tord,
        const len_t nv,
        const int32_t *degs,
        const int32_t *heft
        );

/* The standard grading: A = Z, every variable in degree 1. */
res_dgrp_t *res_dgrp_new_standard(
        const len_t nv
        );

void res_dgrp_free(
        res_dgrp_t **gp
        );

res_dpool_t *res_dpool_new(
        const res_dgrp_t *grp,
        const hl_t sz
        );

void res_dpool_free(
        res_dpool_t **pp
        );

/* View of degree i; the pool must already hold at least i+1 degrees. */
res_deg_t res_dpool_at(
        const res_dpool_t *p,
        const hl_t i
        );

/* Appends a zero degree and returns a view of it, growing if needed.
 * The index of the new degree is written to *idx when idx != NULL.
 * NOTE: any previously handed out res_deg_t is invalidated by a growth,
 * so hold indices rather than views across calls. */
res_deg_t res_dpool_push(
        res_dpool_t *p,
        hl_t *idx
        );

/* --- arithmetic ------------------------------------------------------ */

void res_deg_zero(
        const res_dgrp_t *g,
        res_deg_t d
        );

void res_deg_set(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a
        );

/* d <- multidegree of the monomial with exponent vector exps[0..nv), i.e.
 * the column sum of dmat weighted by the exponents.  exps points at the
 * *variable* slots, not at a raw hash-table exponent vector; callers pass
 * ev + 1 for a non-block order. */
void res_deg_of_exponents(
        const res_dgrp_t *g,
        res_deg_t d,
        const exp_t *exps
        );

/* heft . (dmat . exps), the single integer that schedules and orders the
 * computation and that gets stored in ev[DEG]. */
deg_t res_heft_of_exponents(
        const res_dgrp_t *g,
        const exp_t *exps
        );

/* --------------------------------------------------------------------- *
 *  Gröbner bases of submodules of a free module
 *
 *  export_module_f4 is the module counterpart of export_f4 in f4.h: it
 *  takes a presentation matrix instead of a list of polynomials.  The
 *  columns of that matrix generate a submodule of R^nr_rows, and the flat
 *  input arrays extend msolve's usual convention with one component per
 *  term:
 *
 *    lens   nr_gens entries, the number of terms of each column
 *    exps   nr_vars entries per term, concatenated column by column
 *    comps  one entry per term, the 1-based row that term sits in
 *    cfs    one int32_t coefficient per term
 *
 *  and, describing the ambient free module,
 *
 *    row_degs  nr_rows entries, the degree of each generator of R^nr_rows,
 *              i.e. R^nr_rows = (+)_i R(-row_degs[i]); may be NULL, which
 *              is read as all zero.  Only differences matter, so the
 *              values are normalized to start at zero internally.
 *
 *  The output uses the same layout, with bcomp giving the component of
 *  each term of the basis.  All four output arrays are allocated with the
 *  caller supplied mallocp and are released by
 *  free_module_f4_result_data.  The return value is the total number of
 *  terms, or 0 on failure, in which case nothing is allocated.
 *
 *  Restrictions, all reported on stderr rather than assumed:
 *    - prime field of characteristic 0 < p < 2^31,
 *    - degree reverse lexicographic order (mon_order 0), no elimination
 *      block, since block orders and module orders are not combined yet,
 *    - module_order is RES_MORD_POT or RES_MORD_TOP; RES_MORD_SCHREYER
 *      needs per component base monomials that only the resolution engine
 *      can supply.
 * --------------------------------------------------------------------- */

int64_t export_module_f4(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,      /* number of basis elements */
        int32_t **blen,    /* number of terms of each basis element */
        int32_t **bexp,    /* nr_vars exponents per term */
        int32_t **bcomp,   /* 1-based component of each term */
        void **bcf,        /* one int32_t coefficient per term */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t info_level
        );

void free_module_f4_result_data(
        void (*freep) (void *),
        int32_t **blen,
        int32_t **bexp,
        int32_t **bcomp,
        void **bcf
        );

#ifdef __cplusplus
}
#endif

#endif /* GB_RES_H */
