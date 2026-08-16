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
 * Jérémy Berthomieu
 * Christian Eder
 * Mohab Safey El Din */

#ifndef GB_ORDER_H
#define GB_ORDER_H

#include "data.h"

/* The largest value an exponent, and hence a block degree slot, can
 * hold.  A weighted block degree is a sum of weight*exponent products
 * stored in one such slot, so weights are capped to leave at least
 * MIN_WEIGHTED_DEG_RANGE of usable degree; see check_block_order. */
#define EXP_T_MAX ((int32_t)((1L << (CHAR_BIT * sizeof(exp_t))) - 1))
#define MIN_WEIGHTED_DEG_RANGE 64

/* The block structure of a monomial order, as a caller hands it in.
 *
 * The variables are partitioned into nbl consecutive blocks of sizes
 * bsz[0], ..., bsz[nbl-1]; each block is ordered by the degree reverse
 * lexicographical order and the blocks are compared left to right.  So
 * nbl == 1 is plain DRL and nbl == 2 is the classical elimination order
 * that the elim_block_len entry points expose.
 *
 * bwt gives a positive weight to every variable, in variable order and
 * independent of the blocks, so a block's degree is the weighted sum of
 * its exponents.  NULL is the standard grading, every weight one.
 *
 * The struct borrows bsz and bwt; check_and_set_meta_data copies them
 * into md_t, which owns them from then on. */
typedef struct mo_block_t mo_block_t;
struct mo_block_t
{
    int32_t nbl;        /* number of blocks, >= 1 */
    const int32_t *bsz; /* block sizes, length nbl, summing to nr_vars */
    const int32_t *bwt; /* variable weights, length nr_vars, or NULL */
};

/* The total degree of a monomial: the sum of the per block degree slots.
 * This is what ht->hd[].deg caches, and it is what every degree driven
 * decision in the F4 (pair selection, the f4sat schedule, redundancy)
 * reads.  Note that ev[DEG] alone is only the first block's degree. */
static inline deg_t ht_total_degree(
        const exp_t * const e,
        const ht_t * const ht
        )
{
    const len_t nbl = ht->nbl;

    /* This runs once per monomial inserted into a hash table, which is
     * where msolve spends most of its time, so the two block counts that
     * actually occur are spelled out rather than left to a loop the
     * compiler cannot bound. */
    if (nbl == 1) {
        return (deg_t)e[DEG];
    }
    const len_t * const bst = ht->bst;
    if (nbl == 2) {
        return (deg_t)e[DEG] + (deg_t)e[bst[1]];
    }

    deg_t deg = (deg_t)e[DEG];
    for (len_t b = 1; b < nbl; ++b) {
        deg += (deg_t)e[bst[b]];
    }
    return deg;
}

/* The tail of a block comparison: reverse lexicographically inside block
 * `first`, whose degree slot the caller has already compared, then the
 * remaining blocks in full.  Splitting it out this way is what lets
 * cmp_blocks and cmp_blocks_shifted below differ only in how they
 * compare that one degree, with no second copy of the tie breaks. */
static inline int cmp_blocks_tail(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht,
        const len_t first
        )
{
    const len_t nbl         = ht->nbl;
    const len_t * const bst = ht->bst;

    for (len_t b = first; b < nbl; ++b) {
        const len_t d = bst[b];
        /* every block after the first one still needs its degree */
        if (b > first && ea[d] != eb[d]) {
            return ea[d] > eb[d] ? 1 : -1;
        }
        /* reverse lexicographically inside the block: the last differing
         * variable decides, and the smaller exponent wins */
        len_t i = bst[b+1] - 1;
        while (i > d && ea[i] == eb[i]) {
            --i;
        }
        if (i > d) {
            return (int)eb[i] - (int)ea[i];
        }
        /* this block is identical, move on to the next one */
    }
    return 0;
}

/* Compare the ring parts of ea and eb in the block order, returning a
 * value > 0 exactly when a > b.  Every other comparator in this section
 * is a sign convention away from this one.
 *
 * Weights need no code here: a block's degree slot already holds the
 * weighted sum of that block's exponents (see set_exponent_vector and
 * get_lcm), and the reverse lexicographical tie break is on raw
 * exponents either way.
 *
 * Only slots in [bst[0], bst[nbl]) are read, so this is correct verbatim
 * on a module table: the component slot lives at bst[nbl] and is left to
 * the module order in res_order.c. */
static inline int cmp_blocks(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    const len_t d = ht->bst[0];

    /* the first block's (weighted) degree */
    if (ea[d] != eb[d]) {
        return ea[d] > eb[d] ? 1 : -1;
    }
    return cmp_blocks_tail(ea, eb, ht, 0);
}

/* cmp_blocks with a degree added to the first block before it is
 * compared.  A module monomial carries its component's degree shift in
 * that slot already (see get_lcm and set_module_exponent_vector), but
 * the Schreyer frame keeps its monomials in a plain ring table where the
 * shift has to be supplied by the caller instead.  Apart from that this
 * is cmp_blocks exactly, tie breaks included. */
static inline int cmp_blocks_shifted(
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht,
        const deg_t sa,
        const deg_t sb
        )
{
    const len_t d  = ht->bst[0];
    const deg_t da = (deg_t)ea[d] + sa;
    const deg_t db = (deg_t)eb[d] + sb;

    if (da != db) {
        return da > db ? 1 : -1;
    }
    return cmp_blocks_tail(ea, eb, ht, 0);
}

/* The exponent vector slot of a single variable, together with the
 * degree slot of the block that variable lives in.  Either output may be
 * NULL.  Use ht_variable_slots below when every variable is wanted. */
static inline void ht_variable_slot(
        const ht_t * const ht,
        const len_t var_idx,
        len_t *slot,
        len_t *deg_slot
        )
{
    len_t b = 0;
    /* block b covers variables bst[b]-b .. bst[b+1]-b-2 in variable
     * numbering, so advance while var_idx is past this block's last one */
    while (b+1 < ht->nbl && var_idx >= ht->bst[b+1] - (b+1)) {
        ++b;
    }
    if (slot != NULL) {
        *slot = var_idx + b + 1;
    }
    if (deg_slot != NULL) {
        *deg_slot = ht->bst[b];
    }
}

/* Fill evi[0 .. nv-1] with the exponent vector slot of each variable, in
 * variable order, skipping the per block degree slots.  This is the one
 * place that knows how to walk around them; callers that need to read or
 * write plain nv-long exponent vectors go through here. */
static inline void ht_variable_slots(
        const ht_t * const ht,
        int32_t *evi
        )
{
    len_t ctr = 0;
    for (len_t b = 0; b < ht->nbl; ++b) {
        for (len_t i = ht->bst[b]+1; i < ht->bst[b+1]; ++i) {
            evi[ctr++] = (int32_t)i;
        }
    }
}

#endif
