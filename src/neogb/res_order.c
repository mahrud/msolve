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

/* Monomial orders on a free module.
 *
 * These mirror the ring orders in order.c and are selected by the
 * dispatchers in io.c whenever ht->cpos != 0.  Sign conventions follow
 * order.c exactly: a positive return value means the first argument is
 * the larger monomial.
 *
 * Components are ordered so that a *smaller* component index gives the
 * larger module monomial, matching the existing signature comparison
 * matrix_row_cmp_by_increasing_signature in order.c.
 *
 * Degrees need no special treatment: ev[DEG] of a module monomial already
 * includes the degree shift of its component (see data.h), so comparing
 * ev[DEG] compares module degrees. */

#include "res.h"

/* Index range of the variable slots in a module exponent vector.  The
 * component sits at ht->cpos, which is the last slot, so variables run
 * over [1, cpos).  Block elimination orders are rejected at hash table
 * setup time, so there is no second degree slot to skip here. */

static inline int res_cmp_component(
        const exp_t * const ea,
        const exp_t * const eb,
        const len_t cpos
        )
{
    if (ea[cpos] == eb[cpos]) {
        return 0;
    }
    return ea[cpos] < eb[cpos] ? 1 : -1;
}

/* degree reverse lexicographical on the monomial parts */
static inline int res_cmp_terms_drl(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    len_t i;

    if (ea[DEG] != eb[DEG]) {
        return ea[DEG] > eb[DEG] ? 1 : -1;
    }
    for (i = ht->cpos - 1; i > 0; --i) {
        if (ea[i] != eb[i]) {
            return (int)eb[i] - (int)ea[i];
        }
    }
    return 0;
}

/* lexicographical on the monomial parts */
static inline int res_cmp_terms_lex(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    len_t i;
    const len_t cpos = ht->cpos;

    for (i = 1; i < cpos; ++i) {
        if (ea[i] != eb[i]) {
            return (int)ea[i] - (int)eb[i];
        }
    }
    return 0;
}

/* Schreyer: compare the monomials lifted along the per component base
 * monomials ht->cbase, in the base ring order.  The lift is computed on
 * the fly rather than stored, which costs two extra array reads per slot
 * and keeps the hash table free of derived data. */
static inline int res_cmp_terms_schreyer(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    len_t i;

    const exp_t * const ba = ht->ev[ht->cbase[ea[ht->cpos]]];
    const exp_t * const bb = ht->ev[ht->cbase[eb[ht->cpos]]];

    if (ea[DEG] != eb[DEG]) {
        return ea[DEG] > eb[DEG] ? 1 : -1;
    }
    if (ht->mo == 1) { /* lifted lexicographical */
        for (i = 1; i < ht->cpos; ++i) {
            const int32_t va = (int32_t)ea[i] + (int32_t)ba[i];
            const int32_t vb = (int32_t)eb[i] + (int32_t)bb[i];
            if (va != vb) {
                return va - vb;
            }
        }
        return 0;
    }
    for (i = ht->cpos - 1; i > 0; --i) { /* lifted reverse lexicographical */
        const int32_t va = (int32_t)ea[i] + (int32_t)ba[i];
        const int32_t vb = (int32_t)eb[i] + (int32_t)bb[i];
        if (va != vb) {
            return vb - va;
        }
    }
    return 0;
}

/* The module order proper. */
static inline int res_monomial_cmp_ev(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    int c;
    const len_t cpos = ht->cpos;

    switch (ht->mord) {
        case RES_MORD_POT:
            c = res_cmp_component(ea, eb, cpos);
            if (c != 0) {
                return c;
            }
            return ht->mo == 1 ? res_cmp_terms_lex(ea, eb, ht)
                               : res_cmp_terms_drl(ea, eb, ht);
        case RES_MORD_TOP:
            c = ht->mo == 1 ? res_cmp_terms_lex(ea, eb, ht)
                            : res_cmp_terms_drl(ea, eb, ht);
            if (c != 0) {
                return c;
            }
            return res_cmp_component(ea, eb, cpos);
        default: /* RES_MORD_SCHREYER */
            c = res_cmp_terms_schreyer(ea, eb, ht);
            if (c != 0) {
                return c;
            }
            return res_cmp_component(ea, eb, cpos);
    }
}

static inline int monomial_cmp_mod(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    if (a == b) {
        return 0;
    }
    return res_monomial_cmp_ev(ht->ev[a], ht->ev[b], ht);
}

static int monomial_cmp_pivots_mod(
        const hi_t a,
        const hi_t b,
        const ht_t * const ht
        )
{
#if ORDER_COLUMNS
    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
    /* first known pivots vs. tail terms */
    if (ha.idx != hb.idx) {
        return ha.idx < hb.idx ? 1 : -1;
    }
#endif
    /* columns are sorted decreasingly, hence the sign flip w.r.t.
     * monomial_cmp_mod, matching monomial_cmp_pivots_drl in order.c */
    return -res_monomial_cmp_ev(ht->ev[a], ht->ev[b], ht);
}

/* --------------------------------------------------------------------- *
 *  Entry points used by the dispatchers in io.c
 * --------------------------------------------------------------------- */

static int initial_input_cmp_mod(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;

    const hm_t ha = ((hm_t **)a)[0][OFFSET];
    const hm_t hb = ((hm_t **)b)[0][OFFSET];

    return res_monomial_cmp_ev(ht->ev[ha], ht->ev[hb], ht);
}

static int initial_gens_cmp_mod(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;

    const hm_t ha = **(hm_t **)a;
    const hm_t hb = **(hm_t **)b;

    /* reversed w.r.t. initial_input_cmp_mod, as in order.c */
    return -res_monomial_cmp_ev(ht->ev[ha], ht->ev[hb], ht);
}

static int spair_cmp_mod(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hi_t la  = ((spair_t *)a)->lcm;
    const hi_t lb  = ((spair_t *)b)->lcm;
    const ht_t *ht = (ht_t *)htp;

    const int mc = monomial_cmp_mod(la, lb, ht);
    if (mc != 0) {
        return mc < 0 ? -1 : 1;
    }
    return 0;
}

static int hcm_cmp_pivots_mod(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht = (ht_t *)htp;
    const hi_t ma  = ((hi_t *)a)[0];
    const hi_t mb  = ((hi_t *)b)[0];

    return monomial_cmp_pivots_mod(ma, mb, ht);
}
