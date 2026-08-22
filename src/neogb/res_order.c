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

/* Layout of a module exponent vector: the component sits at ht->cpos,
 * which is the last slot, and the ring part below it is laid out in
 * blocks as data.h describes.  The lexicographic comparison below is the
 * only one that walks the slots by hand, and it runs only for the single
 * block order, where [1, cpos) is exactly the variables; the block
 * orders go through cmp_blocks in order.h, which knows about the per
 * block degree slots and steps over them. */

/* --------------------------------------------------------------------- *
 *  Strategies
 *
 *  The three axes of res_strat_t, bundled at the C entry points and
 *  unpacked into ht->mord / ht->mpos and res_frame_t::strat, which is
 *  where the comparisons read them.  Keeping them as plain int fields
 *  rather than as a struct on the hash table keeps the F4 hot path the
 *  shape it has always had.
 * --------------------------------------------------------------------- */

res_strat_t res_strat_default(
        void
        )
{
    res_strat_t s;

    s.base = RES_MORD_POT;
    s.pos  = RES_POS_DOWN;
    s.lift = RES_LIFT_SCHREYER;

    return s;
}

res_stop_t res_stop_none(
        void
        )
{
    res_stop_t s;

    s.max_degree = NULL;

    return s;
}

res_strat_t res_strat_of_order(
        const int32_t module_order
        )
{
    res_strat_t s = res_strat_default();

    s.base = module_order;

    return s;
}

int res_strat_check(
        const res_strat_t * const s,
        const int for_resolution
        )
{
    if (s == NULL) {
        return 0; /* NULL means the default, which is always usable */
    }
    if (s->pos != RES_POS_DOWN && s->pos != RES_POS_UP) {
        fprintf(ERRSTREAM, "Unknown component direction %d; it is "
                "RES_POS_DOWN or RES_POS_UP.\n", s->pos);
        return 1;
    }
    if (s->lift != RES_LIFT_SCHREYER) {
        fprintf(ERRSTREAM, "Unknown lift %d; only the Schreyer induced "
                "order is implemented for the levels above zero.\n",
                s->lift);
        return 1;
    }
    if (s->base == RES_MORD_SCHREYER) {
        fprintf(ERRSTREAM, "The Schreyer order as a *base* order needs per "
                "component base monomials, which only the resolution engine "
                "can supply, and by then the Gr\u00f6bner basis is already "
                "in hand.  Note this is a different thing from the Schreyer "
                "order the levels above zero are lifted into, which is what "
                "RES_LIFT_SCHREYER means and is always in force.\n");
        return 1;
    }
    if (s->base != RES_MORD_POT && s->base != RES_MORD_TOP) {
        fprintf(ERRSTREAM, "Unknown base module order %d.\n", s->base);
        return 1;
    }
    (void)for_resolution; /* both bases carry a frame and a differential */

    return 0;
}

const char *res_strat_name(
        const res_strat_t * const s
        )
{
    const res_strat_t d = res_strat_default();
    const res_strat_t * const t = s != NULL ? s : &d;

    if (t->lift != RES_LIFT_SCHREYER) {
        return "unknown";
    }
    switch (t->base) {
        case RES_MORD_POT:
            return t->pos == RES_POS_UP ? "pot-up-schreyer"
                                        : "pot-down-schreyer";
        case RES_MORD_TOP:
            return t->pos == RES_POS_UP ? "top-up-schreyer"
                                        : "top-down-schreyer";
        default:
            return "unknown";
    }
}

/* The component key, in the direction the strategy asks for.  pos is
 * ht->mpos, constant for a whole computation, so the branch predicts
 * perfectly and the comparison still inlines into the sort. */
static inline int res_cmp_component(
        const exp_t * const ea,
        const exp_t * const eb,
        const len_t cpos,
        const int32_t pos
        )
{
    if (ea[cpos] == eb[cpos]) {
        return 0;
    }
    if (pos == RES_POS_UP) {
        return ea[cpos] > eb[cpos] ? 1 : -1;
    }
    return ea[cpos] < eb[cpos] ? 1 : -1;
}

/* Block grevlex on the monomial parts, which for a single block is
 * ordinary degree reverse lexicographical.  cmp_blocks reads only the
 * slots below ht->cpos, so the component is untouched here and stays the
 * business of res_cmp_component above. */
static inline int res_cmp_terms_ring(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    return cmp_blocks(ea, eb, ht);
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
            c = res_cmp_component(ea, eb, cpos, ht->mpos);
            if (c != 0) {
                return c;
            }
            return ht->mo == 1 ? res_cmp_terms_lex(ea, eb, ht)
                               : res_cmp_terms_ring(ea, eb, ht);
        case RES_MORD_TOP:
            c = ht->mo == 1 ? res_cmp_terms_lex(ea, eb, ht)
                            : res_cmp_terms_ring(ea, eb, ht);
            if (c != 0) {
                return c;
            }
            return res_cmp_component(ea, eb, cpos, ht->mpos);
        default: /* RES_MORD_SCHREYER */
            c = res_cmp_terms_schreyer(ea, eb, ht);
            if (c != 0) {
                return c;
            }
            return res_cmp_component(ea, eb, cpos, ht->mpos);
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
