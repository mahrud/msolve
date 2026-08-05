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

/* Schreyer frame construction.
 *
 * Pure combinatorics: the frame is built from monomials alone, with no
 * coefficients and no field arithmetic anywhere in this file.  See res.h
 * for the shape of the recursion and for the storage order, which the
 * construction below relies on throughout. */

#include "res.h"
#include "../msolve/streams.h"

/* --------------------------------------------------------------------- *
 *  The frame hash table
 *
 *  A plain ring table: cpos == 0, ebl == 0, evl == nv+1.  The component
 *  of a frame element is its parent index, which lives in the element
 *  rather than in the monomial, so there are far fewer components' worth
 *  of bookkeeping than a module table at every level would need.
 * --------------------------------------------------------------------- */

/* insert_in_hash_table does not grow the table, so callers reserve first;
 * this is the same contract import_input_data works under. */
static inline void res_frame_ht_reserve(
        ht_t *ht,
        const hl_t n
        )
{
    while (ht->eld + n + 1 >= ht->esz) {
        enlarge_hash_table(ht);
    }
}

/* e <- a / gcd(a,b), inserted into the table.  Degrees are the plain
 * total degree, matching the module table the Gröbner basis was computed
 * in; the graded degrees of a frame element are kept separately in hdeg
 * and mdeg. */
static inline hi_t res_frame_colon(
        ht_t *ht,
        exp_t *e,
        const hm_t a,
        const hm_t b
        )
{
    len_t i;
    deg_t d = 0;

    const len_t evl       = ht->evl;
    const exp_t * const ea = ht->ev[a];
    const exp_t * const eb = ht->ev[b];

    for (i = 1; i < evl; ++i) {
        e[i] = ea[i] > eb[i] ? (exp_t)(ea[i] - eb[i]) : 0;
        d    = d + (deg_t)e[i];
    }
    e[0] = (exp_t)d;

    return insert_in_hash_table(e, ht);
}

/* e <- a * b, inserted into the table.
 *
 * Deliberately not insert_multiplied_signature_in_hash_table: that one
 * takes hash table indices but computes the hash of the product as the
 * *sum of the two indices* (hash.c:725) rather than as the linear hash
 * sum(rn[j]*a[j]) every other path uses.  Entries it inserts are
 * therefore invisible to insert_in_hash_table, which silently duplicates
 * them under a different index -- so the same monomial ends up with two
 * hashes and nothing dedupes.  Its only other caller is the dead SBA
 * path, which uses it for every insert and lookup alike and so never
 * notices. */
static inline int res_frame_product(
        const ht_t *ht,
        exp_t *e,
        const hm_t a,
        const hm_t b
        )
{
    len_t i;

    const len_t evl        = ht->evl;
    const exp_t * const ea = ht->ev[a];
    const exp_t * const eb = ht->ev[b];

    for (i = 0; i < evl; ++i) {
        const uint32_t s = (uint32_t)ea[i] + (uint32_t)eb[i];
        if (s > UINT16_MAX) {
            return 1;
        }
        e[i] = (exp_t)s;
    }

    return 0;
}

static inline hi_t res_frame_mul(
        ht_t *ht,
        exp_t *e,
        const hm_t a,
        const hm_t b
        )
{
    if (res_frame_product(ht, e, a, b)) {
        return 0;
    }

    return insert_in_hash_table(e, ht);
}

/* Lookup-only counterpart used by verification.  A verifier must not
 * mutate the object it checks, especially when the object is malformed. */
static inline hi_t res_frame_find_product(
        const ht_t *ht,
        exp_t *e,
        const hm_t a,
        const hm_t b
        )
{
    hl_t i;
    len_t j;
    val_t h = 0;

    if (res_frame_product(ht, e, a, b)) {
        return 0;
    }
    for (j = 0; j < ht->evl; ++j) {
        h += ht->rn[j] * e[j];
    }

    hi_t k = h;
    const hi_t mod = (hi_t)(ht->hsz - 1);
    for (i = 0; i < ht->hsz; ++i) {
        k = (hi_t)((k + i) & mod);
        const hi_t hm = ht->hmap[k];
        if (hm == 0) {
            return 0;
        }
        if (ht->hd[hm].val == h
                && memcmp(e, ht->ev[hm],
                    (unsigned long)ht->evl * sizeof(exp_t)) == 0) {
            return hm;
        }
    }
    return 0;
}

/* --------------------------------------------------------------------- *
 *  Comparators
 * --------------------------------------------------------------------- */

/* Ascending total degree, then hash index.  Sorting the colon quotients
 * this way is what lets the minimalization below be a single greedy pass:
 * a proper divisor always has strictly smaller degree, so it is seen
 * first, and two monomials of equal degree divide each other only if they
 * are equal, in which case they share a hash index. */
static int res_frame_cmp_deg(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht = (ht_t *)htp;
    const hi_t x   = ((const hi_t *)a)[0];
    const hi_t y   = ((const hi_t *)b)[0];

    if (ht->hd[x].deg != ht->hd[y].deg) {
        return ht->hd[x].deg < ht->hd[y].deg ? -1 : 1;
    }
    return x < y ? -1 : (x > y ? 1 : 0);
}

/* Degree ascending, then the ring order *descending*; the hash index
 * breaks ties only for equal monomials, which cannot both survive
 * minimalization, so it is there purely to make the sort deterministic.
 *
 * The direction is the one free choice in the whole construction, every
 * choice gives a correct frame, and it decides how big the frame gets and
 * whether it stops at all.  Schreyer's theorem attaches the lead term of
 * the S-pair of l and k to whichever index is larger, so the elements
 * early in a block are the ones being divided *into*, and the colon ideal
 * that decides element k's children is generated by the quotients
 * m_l / gcd(m_l, m_k) for l < k.
 *
 * What makes this particular direction the right one is Schreyer's proof
 * of Hilbert's syzygy theorem.  Among monomials of the same degree the
 * reverse lexicographic order is descending in the exponent of the last
 * variable that distinguishes two of them, so sorting the block by
 * descending order puts the elements with the *smallest* exponent in the
 * highest variable first.  Every colon quotient m_l / gcd(m_l, m_k) then
 * avoids the highest variable, one more variable drops out at each level,
 * and the frame is guaranteed to stop by level nv.  The ascending
 * direction has no such bound: on four generic quadrics in four variables
 * it runs to level five, which the nv ceiling then silently truncates.
 *
 * This is also exactly what Macaulay2 does -- res-f4-computation.cpp
 * sorts level one with sort(1 ascending degree, -1 descending monomial
 * order), and the PreElementSorter of res-schreyer-frame.cpp sorts each
 * colon quotient by degree and then by varpower comparison, which lists
 * variables from the highest down and so agrees with descending reverse
 * lex -- and it is what the corpus in test/neogb/res compares against.
 *
 * Note the degree key is explicit rather than left to the ring order.
 * For degree reverse lexicographic the two agree, but the rule is
 * "degree ascending, order descending", not "order descending". */
static int res_frame_cmp_mono(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht = (ht_t *)htp;
    const hi_t x   = ((const hi_t *)a)[0];
    const hi_t y   = ((const hi_t *)b)[0];

    if (ht->hd[x].deg != ht->hd[y].deg) {
        return ht->hd[x].deg < ht->hd[y].deg ? -1 : 1;
    }
    const int c = monomial_cmp(x, y, ht);
    if (c != 0) {
        return c > 0 ? -1 : 1;
    }
    return x < y ? -1 : (x > y ? 1 : 0);
}

/* Level 1 in storage order: component ascending, lead monomial ascending,
 * the same shape as every level above it.  Only the second key is visible
 * in the result, since the Gröbner basis elements of two different
 * components never pair up. */
typedef struct res_frame_gen_t res_frame_gen_t;
struct res_frame_gen_t
{
    hi_t    mono;
    int32_t comp;
    int32_t src;   /* index in the Gröbner basis this lead term came from;
                    * carried along so that res_diff.c can recover the
                    * coefficients after this sort has happened */
};

static int res_frame_cmp_gen(
        const void *a,
        const void *b,
        void *htp
        )
{
    const res_frame_gen_t * const x = (const res_frame_gen_t *)a;
    const res_frame_gen_t * const y = (const res_frame_gen_t *)b;

    if (x->comp != y->comp) {
        return x->comp < y->comp ? -1 : 1;
    }
    return res_frame_cmp_mono(&x->mono, &y->mono, htp);
}

/* --------------------------------------------------------------------- *
 *  Levels
 * --------------------------------------------------------------------- */

static int res_level_enlarge(
        res_level_t *lv,
        const len_t need
        )
{
    if (lv->ld + need <= lv->sz) {
        return 0;
    }
    len_t nsz = lv->sz > 0 ? lv->sz : 16;
    while (lv->ld + need > nsz) {
        nsz = 2 * nsz;
    }
    res_felt_t *ne = (res_felt_t *)realloc(
            lv->elts, (unsigned long)nsz * sizeof(res_felt_t));
    if (ne == NULL) {
        fprintf(ERRSTREAM, "Could not enlarge a frame level to %u elements.\n",
                (unsigned)nsz);
        return 1;
    }
    lv->elts = ne;
    lv->sz   = nsz;

    return 0;
}

/* mdeg is read from raw degree slots rather than from a res_deg_t, since
 * pushing into the pool may move every view previously handed out. */
static int res_level_push(
        res_frame_t *f,
        const len_t lev,
        const hm_t mono,
        const hm_t total,
        const int32_t up,
        const int32_t root,
        const deg_t hdeg,
        const int32_t * const mdeg
        )
{
    res_level_t *lv = f->lv + lev;

    if (res_level_enlarge(lv, 1)) {
        return 1;
    }

    hl_t idx  = 0;
    res_deg_t d = res_dpool_push(lv->degs, &idx);
    if (d.e == NULL) {
        return 1;
    }
    if (mdeg != NULL) {
        memcpy(d.e, mdeg, (unsigned long)f->grp->len * sizeof(int32_t));
    }

    res_felt_t *el = lv->elts + lv->ld;
    el->mono  = mono;
    el->total = total;
    el->up    = up;
    el->root  = root;
    el->hdeg  = hdeg;
    el->mdeg  = idx;

    lv->ld++;

    return 0;
}

/* --------------------------------------------------------------------- *
 *  Construction and teardown
 * --------------------------------------------------------------------- */

res_frame_t *res_frame_new(
        const res_dgrp_t *grp,
        const md_t *md,
        const int32_t maxlevel
        )
{
    len_t i;

    if (grp == NULL || md == NULL) {
        return NULL;
    }

    res_frame_t *f = (res_frame_t *)calloc(1, sizeof(res_frame_t));
    if (f == NULL) {
        return NULL;
    }

    f->grp   = grp;
    f->nv    = md->nvars;
    f->nlv   = 0;
    f->bad   = 0;
    /* Zero is no ceiling at all.  There is no cheap a priori bound worth
     * using here: the frame is a *nonminimal* resolution, so Hilbert's
     * syzygy theorem does not bound its length, and it really does run
     * past level nv -- (z, y^2, x^2 y, x^3) in three variables reaches
     * level four, and Macaulay2 agrees element for element once its own
     * LengthLimit is raised.  What does bound it is that the frame is a
     * subcomplex of the Taylor complex on the Gröbner basis, so it stops
     * by level lv[1].ld; that is checked in res_frame_complete rather
     * than preallocated, since it is usually an enormous overestimate. */
    f->maxlv = maxlevel > 0 ? (len_t)maxlevel : 0;

    /* a private ring table: same variables and same ring order as the
     * module table the Gröbner basis lives in, but no component slot */
    md_t rmd  = *md;
    rmd.ncomp = 0;
    rmd.mord  = 0;
    rmd.nev   = 0;
    f->ht     = initialize_basis_hash_table(&rmd);
    if (f->ht == NULL) {
        free(f);
        return NULL;
    }

    f->lvsz = f->maxlv > 0 ? f->maxlv + 1 : 8;
    f->lv   = (res_level_t *)calloc(
            (unsigned long)f->lvsz, sizeof(res_level_t));
    if (f->lv == NULL) {
        full_free_hash_table(&f->ht);
        free(f);
        return NULL;
    }
    for (i = 0; i < f->lvsz; ++i) {
        f->lv[i].degs = res_dpool_new(grp, 16);
        if (f->lv[i].degs == NULL) {
            res_frame_free(&f);
            return NULL;
        }
    }

    return f;
}

/* Makes level lev addressable.  Only reached when there is no ceiling,
 * since a ceiling is preallocated in full. */
static int res_frame_reserve_levels(
        res_frame_t *f,
        const len_t lev
        )
{
    len_t i;

    if (lev < f->lvsz) {
        return 0;
    }
    len_t nsz = f->lvsz > 0 ? f->lvsz : 8;
    while (lev >= nsz) {
        nsz = 2 * nsz;
    }
    res_level_t *nl = (res_level_t *)realloc(
            f->lv, (unsigned long)nsz * sizeof(res_level_t));
    if (nl == NULL) {
        return 1;
    }
    memset(nl + f->lvsz, 0,
            (unsigned long)(nsz - f->lvsz) * sizeof(res_level_t));
    for (i = f->lvsz; i < nsz; ++i) {
        nl[i].degs = res_dpool_new(f->grp, 16);
        if (nl[i].degs == NULL) {
            /* the slots already made are owned by the frame and are freed
             * with it; only the growth failed */
            f->lv   = nl;
            f->lvsz = i;
            return 1;
        }
    }
    f->lv   = nl;
    f->lvsz = nsz;

    return 0;
}

void res_frame_free(
        res_frame_t **fp
        )
{
    len_t i;
    res_frame_t *f = *fp;

    if (f == NULL) {
        return;
    }
    if (f->lv != NULL) {
        for (i = 0; i < f->lvsz; ++i) {
            free(f->lv[i].elts);
            res_dpool_free(&f->lv[i].degs);
        }
        free(f->lv);
    }
    free(f->gbmap);
    if (f->ht != NULL) {
        full_free_hash_table(&f->ht);
    }
    free(f);
    *fp = NULL;
}

int res_frame_init(
        res_frame_t *f,
        const bs_t * const gb,
        const ht_t * const bht,
        const int32_t *row_mdegs
        )
{
    len_t i, j;
    int ret = 1;

    if (f == NULL || gb == NULL || bht == NULL || bht->cpos == 0) {
        fprintf(ERRSTREAM, "A Schreyer frame needs a module Gröbner basis.\n");
        return 1;
    }
    if (f->nlv != 0) {
        fprintf(ERRSTREAM, "The Schreyer frame is already initialized.\n");
        return 1;
    }

    ht_t *ht          = f->ht;
    const len_t evl   = ht->evl;
    const len_t nv    = f->nv;
    const len_t glen  = f->grp->len;

    f->ncomp = bht->ncomp;

    exp_t *e = (exp_t *)calloc((unsigned long)evl, sizeof(exp_t));
    res_frame_gen_t *gen = (res_frame_gen_t *)malloc(
            (unsigned long)(gb->lml > 0 ? gb->lml : 1) * sizeof(res_frame_gen_t));
    int32_t *md = (int32_t *)calloc((unsigned long)glen, sizeof(int32_t));
    if (e == NULL || gen == NULL || md == NULL) {
        goto cleanup;
    }

    /* --- level 0: the generators of the ambient free module ---------- */

    res_frame_ht_reserve(ht, (hl_t)f->ncomp + (hl_t)gb->lml + 1);

    const hi_t one = insert_in_hash_table(e, ht);

    f->nlv = 1;
    for (i = 0; i < f->ncomp; ++i) {
        const int32_t *rd = row_mdegs != NULL
            ? row_mdegs + (unsigned long)i * glen : NULL;
        if (res_level_push(f, 0, one, one, -1, (int32_t)i + 1,
                    bht->cshift[i+1], rd)) {
            goto cleanup;
        }
    }

    /* --- level 1: the lead terms of the Gröbner basis ---------------- */

    len_t ngen = 0;
    for (i = 0; i < gb->lml; ++i) {
        const hm_t *hm = gb->hm[gb->lmps[i]];
        if (hm == NULL) {
            continue; /* the zero element, carrying no lead term */
        }
        const exp_t * const be = bht->ev[hm[OFFSET]];
        deg_t d = 0;
        for (j = 1; j <= nv; ++j) {
            e[j] = be[j];
            d    = d + (deg_t)be[j];
        }
        e[0] = (exp_t)d;
        gen[ngen].mono = insert_in_hash_table(e, ht);
        gen[ngen].comp = (int32_t)be[bht->cpos];
        gen[ngen].src  = (int32_t)gb->lmps[i];
        ngen++;
    }

    sort_r(gen, (unsigned long)ngen, sizeof(res_frame_gen_t),
            res_frame_cmp_gen, ht);

    f->gbmap = (int32_t *)malloc(
            (unsigned long)(ngen > 0 ? ngen : 1) * sizeof(int32_t));
    if (f->gbmap == NULL) {
        goto cleanup;
    }
    for (i = 0; i < ngen; ++i) {
        f->gbmap[i] = gen[i].src;
    }

    /* Divisor masks are derived from the exponents seen so far, and every
     * later frame monomial is a colon quotient of these, hence bounded by
     * them.  Computing them here, once level 1 is in, is what makes the
     * divisibility tests in the minimalization below cheap. */
    calculate_divmask(ht);

    f->nlv = 2;
    for (i = 0; i < ngen; ++i) {
        const int32_t c = gen[i].comp;
        if (c < 1 || c > (int32_t)f->ncomp) {
            fprintf(ERRSTREAM, "Gröbner basis element %u has component %d, "
                    "outside the ambient free module.\n", (unsigned)i, c);
            goto cleanup;
        }
        const res_felt_t * const par = f->lv[0].elts + (c - 1);
        const res_deg_t pd = res_dpool_at(f->lv[0].degs, par->mdeg);
        res_deg_t qd = {md};
        res_deg_of_exponents(f->grp, qd, ht->ev[gen[i].mono] + 1);
        f->grp->add(f->grp, qd, qd, pd);
        const deg_t hd = par->hdeg
            + res_heft_of_exponents(f->grp, ht->ev[gen[i].mono] + 1);
        if (res_level_push(f, 1, gen[i].mono, gen[i].mono,
                    c - 1, c, hd, md)) {
            goto cleanup;
        }
    }

    ret = 0;

cleanup:
    free(e);
    free(gen);
    free(md);

    f->bad = ret;

    return ret;
}

len_t res_frame_next_level(
        res_frame_t *f
        )
{
    len_t i, k, l;
    len_t added = 0;

    if (f == NULL || f->nlv < 2
            || (f->maxlv > 0 && f->nlv > f->maxlv)) {
        return 0;
    }
    if (res_frame_reserve_levels(f, f->nlv)) {
        fprintf(ERRSTREAM, "Out of memory making room for frame level %u.\n",
                (unsigned)f->nlv);
        f->bad = 1;
        return 0;
    }

    const len_t src = f->nlv - 1;
    const len_t dst = f->nlv;

    res_level_t *sl = f->lv + src;
    if (sl->ld == 0) {
        return 0;
    }

    ht_t *ht        = f->ht;
    const len_t glen = f->grp->len;

    /* block starts: elements are stored with up ascending, so the
     * elements sharing a parent form one contiguous run */
    len_t *bstart = (len_t *)malloc((unsigned long)sl->ld * sizeof(len_t));
    hi_t *qs      = (hi_t *)malloc((unsigned long)sl->ld * sizeof(hi_t));
    hi_t *qm      = (hi_t *)malloc((unsigned long)sl->ld * sizeof(hi_t));
    exp_t *e      = (exp_t *)calloc((unsigned long)ht->evl, sizeof(exp_t));
    int32_t *md   = (int32_t *)calloc((unsigned long)glen, sizeof(int32_t));

    if (bstart == NULL || qs == NULL || qm == NULL || e == NULL || md == NULL) {
        fprintf(ERRSTREAM, "Out of memory building frame level %u.\n",
                (unsigned)dst);
        f->bad = 1;
        goto cleanup;
    }

    bstart[0] = 0;
    for (k = 1; k < sl->ld; ++k) {
        bstart[k] = sl->elts[k].up == sl->elts[k-1].up ? bstart[k-1] : k;
    }

    f->nlv = dst + 1;

    for (k = 0; k < sl->ld; ++k) {
        const len_t b  = bstart[k];
        const len_t nq = k - b;
        if (nq == 0) {
            continue;
        }

        /* one colon quotient per earlier element of the block, and one
         * lifted monomial per survivor of the minimalization */
        res_frame_ht_reserve(ht, 2 * (hl_t)nq + 2);

        const hm_t mk = sl->elts[k].mono;
        for (l = 0; l < nq; ++l) {
            qs[l] = res_frame_colon(ht, e, sl->elts[b+l].mono, mk);
        }

        /* minimal generators of the ideal quotient: greedy over ascending
         * degree, keeping what no earlier survivor divides */
        sort_r(qs, (unsigned long)nq, sizeof(hi_t), res_frame_cmp_deg, ht);

        len_t nmin = 0;
        for (l = 0; l < nq; ++l) {
            int keep = 1;
            for (i = 0; i < nmin; ++i) {
                if (check_monomial_division(qs[l], qm[i], ht)) {
                    keep = 0;
                    break;
                }
            }
            if (keep) {
                qm[nmin++] = qs[l];
            }
        }

        /* storage order inside the new block is the ring order on the own
         * monomials, which for a fixed parent is the Schreyer order */
        sort_r(qm, (unsigned long)nmin, sizeof(hi_t),
                res_frame_cmp_mono, ht);

        res_deg_t pd = res_dpool_at(sl->degs, sl->elts[k].mdeg);
        for (l = 0; l < nmin; ++l) {
            const hm_t q  = qm[l];
            const hm_t tt = res_frame_mul(ht, e, q, sl->elts[k].total);
            if (tt == 0) {
                fprintf(ERRSTREAM, "A lifted frame monomial exceeds the "
                        "16-bit exponent limit.\n");
                f->bad = 1;
                goto cleanup;
            }
            res_deg_t qd = {md};
            res_deg_of_exponents(f->grp, qd, ht->ev[q] + 1);
            f->grp->add(f->grp, qd, qd, pd);
            const deg_t hd = sl->elts[k].hdeg
                + res_heft_of_exponents(f->grp, ht->ev[q] + 1);
            if (res_level_push(f, dst, q, tt, (int32_t)k,
                        sl->elts[k].root, hd, md)) {
                f->bad = 1;
                goto cleanup;
            }
            added++;
            /* pd is a view into the source pool, which is never pushed to
             * here, so it stays valid across the loop */
        }
    }

    if (added == 0) {
        f->nlv = dst; /* an empty level is not a level */
    }

cleanup:
    free(bstart);
    free(qs);
    free(qm);
    free(e);
    free(md);

    return added;
}

int64_t res_frame_complete(
        res_frame_t *f
        )
{
    len_t i;
    int64_t total = 0;

    if (f == NULL) {
        return 0;
    }
    while (f->nlv >= 2 && (f->maxlv == 0 || f->nlv <= f->maxlv)
            && f->bad == 0) {
        if (res_frame_next_level(f) == 0) {
            break;
        }
        /* The frame is a subcomplex of the Taylor complex on level 1, so
         * it cannot outlive that many levels.  Reaching this means the
         * construction is wrong, not that the input was hard. */
        if (f->nlv > f->lv[1].ld + 2) {
            fprintf(ERRSTREAM, "The frame passed level %u, more than the "
                    "Taylor complex on %u generators allows.\n",
                    (unsigned)f->nlv, (unsigned)f->lv[1].ld);
            f->bad = 1;
        }
    }
    if (f->bad) {
        return -1;
    }
    for (i = 0; i < f->nlv; ++i) {
        total += (int64_t)f->lv[i].ld;
    }

    return total;
}

int res_frame_is_complete(
        const res_frame_t * const f
        )
{
    /* res_frame_complete stops for one of two reasons: a level came out
     * empty, which is the frame ending on its own, or the ceiling was
     * reached, which leaves nlv one past maxlv.  Only the first is a whole
     * resolution.  There is deliberately no "nv is always enough" escape
     * here: the frame is nonminimal and can be longer than nv, so a
     * ceiling is a truncation whatever it was set to. */
    if (f == NULL || f->bad != 0) {
        return 0;
    }
    return f->maxlv == 0 || f->nlv <= f->maxlv;
}

/* --------------------------------------------------------------------- *
 *  Structural check
 * --------------------------------------------------------------------- */

int res_frame_verify(
        const res_frame_t * const f
        )
{
    len_t i, k;
    int bad = 0;

    if (f == NULL) {
        return 1;
    }

    const ht_t *ht = f->ht;

    exp_t *e = (exp_t *)calloc((unsigned long)ht->evl, sizeof(exp_t));
    if (e == NULL) {
        return 1;
    }

    for (i = 0; i < f->nlv; ++i) {
        const res_level_t * const lv = f->lv + i;
        for (k = 0; k < lv->ld; ++k) {
            const res_felt_t * const el = lv->elts + k;

            if (i == 0) {
                if (el->up != -1 || el->mono != el->total
                        || el->root != (int32_t)k + 1) {
                    fprintf(ERRSTREAM, "Frame element (0,%u) is not a "
                            "generator of the ambient free module.\n",
                            (unsigned)k);
                    bad++;
                }
                continue;
            }

            if (el->up < 0 || el->up >= (int32_t)f->lv[i-1].ld) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) points at %d, "
                        "outside level %u.\n",
                        (unsigned)i, (unsigned)k, el->up, (unsigned)i - 1);
                bad++;
                continue;
            }
            const res_felt_t * const par = f->lv[i-1].elts + el->up;

            /* the lift is a lookup, not an insertion: every product below
             * was put into the table when the element was created */
            const hm_t tt = res_frame_find_product(
                    ht, e, el->mono, par->total);
            if (tt != el->total) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) has a lifted "
                        "monomial that is not its own times its parent's.\n",
                        (unsigned)i, (unsigned)k);
                bad++;
            }
            if (el->root != par->root) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) does not lift to "
                        "the same component as its parent.\n",
                        (unsigned)i, (unsigned)k);
                bad++;
            }
            if (el->hdeg != par->hdeg
                    + res_heft_of_exponents(f->grp, ht->ev[el->mono] + 1)) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) has an "
                        "inconsistent heft degree.\n",
                        (unsigned)i, (unsigned)k);
                bad++;
            }
            /* blocks are contiguous and ascending; the whole construction
             * reads the previous elements of a block off this */
            if (k > 0 && el->up < lv->elts[k-1].up) {
                fprintf(ERRSTREAM, "Level %u is not sorted by parent at "
                        "element %u.\n", (unsigned)i, (unsigned)k);
                bad++;
            }
            if (bad > 8) {
                fprintf(ERRSTREAM, "Too many frame inconsistencies.\n");
                free(e);
                return bad;
            }
        }
    }
    free(e);

    return bad;
}

/* --------------------------------------------------------------------- *
 *  Frame ranks
 * --------------------------------------------------------------------- */

deg_t res_frame_max_hdeg(
        const res_frame_t * const f
        )
{
    len_t i, k;
    deg_t m = 0;

    if (f == NULL) {
        return 0;
    }
    for (i = 0; i < f->nlv; ++i) {
        for (k = 0; k < f->lv[i].ld; ++k) {
            if (f->lv[i].elts[k].hdeg > m) {
                m = f->lv[i].elts[k].hdeg;
            }
        }
    }
    return m;
}

int64_t res_frame_betti(
        const res_frame_t * const f,
        int32_t *tab,
        const deg_t maxdeg
        )
{
    len_t i, k;
    int64_t total = 0;

    if (f == NULL || tab == NULL) {
        return 0;
    }
    memset(tab, 0, (unsigned long)f->nlv
            * (unsigned long)(maxdeg + 1) * sizeof(int32_t));

    for (i = 0; i < f->nlv; ++i) {
        for (k = 0; k < f->lv[i].ld; ++k) {
            const deg_t d = f->lv[i].elts[k].hdeg;
            if (d < 0 || d > maxdeg) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) has degree %d, "
                        "outside the requested range.\n",
                        (unsigned)i, (unsigned)k, d);
                return 0;
            }
            tab[(unsigned long)i * (maxdeg + 1) + d]++;
            total++;
        }
    }

    return total;
}
