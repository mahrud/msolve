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

/* Minimal Betti numbers by rank extraction, and Hilbert numerators.
 *
 * See res.h for the formula.  The point of the whole construction is that
 * the minimal resolution is never built: the only thing extracted from the
 * nonminimal differential is the rank of its *scalar* part, one small
 * independent matrix over F_p per (level, degree), and those ranks are
 * enough to recover every minimal Betti number.  Nothing here ever
 * back substitutes, so the elimination below stops at row echelon form.
 *
 * The Hilbert numerator is cheaper still: it is the alternating sum of the
 * frame ranks, so it needs no field arithmetic at all beyond the Gröbner
 * basis that produced the frame.  That it also equals the alternating sum
 * of the *minimal* Betti numbers is the telescoping identity
 *
 *     sum_i (-1)^i ( rank(d_i)_d + rank(d_{i+1})_d )  =  0,
 *
 * which res_selftest.c checks both ways round. */

#include "res.h"
#include "../msolve/streams.h"

/* --------------------------------------------------------------------- *
 *  Row echelon form of a scalar block
 *
 *  Rows are the frame elements at (lev, d), columns the frame elements at
 *  (lev-1, d), and the entry is the coefficient of the term of D_lev(row)
 *  whose ring monomial is 1.  Only the rank is wanted, so pivot rows are
 *  kept but never back substituted, and the column order is irrelevant.
 *
 *  Pivot rows are stored sparsely in one arena that is reset per block:
 *  they start as sparse as the differential and only fill in as far as the
 *  elimination actually makes them, which is never worse than the dense
 *  representation and is usually far better.
 * --------------------------------------------------------------------- */

/* Accumulator bound.  Every summand mult*val is below fc^2 <= 2^62, and a
 * value is reduced as soon as it reaches 2^62, so the addition itself
 * always happens below 2^63 and cannot wrap. */
#define RES_ACC_LIMIT (((uint64_t)1) << 62)

typedef struct res_ech_t res_ech_t;
struct res_ech_t
{
    uint64_t *dns;   /* dense accumulator, one slot per column       */
    int32_t  *piv;   /* piv[j] = pivot row leading in column j, or -1 */
    len_t     ncol;  /* slots allocated in dns and piv                */

    /* the pivot rows, appended into one arena */
    int32_t  *beg;   /* beg[r], end[r]: the arena range of pivot row r */
    int32_t  *end;
    len_t     nr;    /* pivot rows in use          */
    len_t     rsz;   /* slots allocated in beg/end */
    int32_t  *idx;   /* column indices, arena      */
    uint32_t *val;   /* coefficients, arena        */
    int64_t   ld;    /* arena slots in use         */
    int64_t   sz;    /* arena slots allocated      */
};

static void res_ech_clear(
        res_ech_t *e
        )
{
    free(e->dns);
    free(e->piv);
    free(e->beg);
    free(e->end);
    free(e->idx);
    free(e->val);
    memset(e, 0, sizeof(res_ech_t));
}

/* Makes room for a block of ncol columns and nrow rows.  The workspace is
 * reused across blocks, so this only ever grows. */
static int res_ech_reserve(
        res_ech_t *e,
        const len_t ncol,
        const len_t nrow
        )
{
    len_t j;

    if (ncol > e->ncol) {
        uint64_t *nd = (uint64_t *)realloc(
                e->dns, (unsigned long)ncol * sizeof(uint64_t));
        int32_t *np  = (int32_t *)realloc(
                e->piv, (unsigned long)ncol * sizeof(int32_t));
        if (nd == NULL || np == NULL) {
            free(nd);
            free(np);
            return 1;
        }
        e->dns  = nd;
        e->piv  = np;
        e->ncol = ncol;
        memset(e->dns, 0, (unsigned long)ncol * sizeof(uint64_t));
    }
    if (nrow > e->rsz) {
        int32_t *nb = (int32_t *)realloc(
                e->beg, (unsigned long)nrow * sizeof(int32_t));
        int32_t *ne = (int32_t *)realloc(
                e->end, (unsigned long)nrow * sizeof(int32_t));
        if (nb == NULL || ne == NULL) {
            free(nb);
            free(ne);
            return 1;
        }
        e->beg = nb;
        e->end = ne;
        e->rsz = nrow;
    }
    for (j = 0; j < ncol; ++j) {
        e->piv[j] = -1;
    }
    e->nr = 0;
    e->ld = 0;

    return 0;
}

static int res_ech_arena_reserve(
        res_ech_t *e,
        const int64_t want
        )
{
    if (e->ld + want <= e->sz) {
        return 0;
    }
    /* beg[] and end[] index the arena with an int32_t, so that is the
     * ceiling; a single block needing two billion nonzeros is a different
     * problem from the one this file solves. */
    if (e->ld + want > (int64_t)INT32_MAX) {
        fprintf(ERRSTREAM, "A scalar block of the differential is too "
                "large to eliminate.\n");
        return 1;
    }
    int64_t nsz = e->sz > 0 ? e->sz : 1024;
    while (nsz < e->ld + want) {
        nsz *= 2;
    }
    if (nsz > (int64_t)INT32_MAX) {
        nsz = (int64_t)INT32_MAX;
    }
    int32_t *ni  = (int32_t *)realloc(
            e->idx, (unsigned long)nsz * sizeof(int32_t));
    uint32_t *nv = (uint32_t *)realloc(
            e->val, (unsigned long)nsz * sizeof(uint32_t));
    if (ni == NULL || nv == NULL) {
        free(ni);
        free(nv);
        return 1;
    }
    e->idx = ni;
    e->val = nv;
    e->sz  = nsz;

    return 0;
}

/* The rank of the scalar part of d_lev in degree d.  rows lists the frame
 * elements at (lev, d), dpos maps a frame element at level lev-1 to its
 * position inside its own degree class -- which is its column, since every
 * scalar term of a row of degree d lands on a generator of degree d.
 * Returns -1 on failure. */
static int64_t res_ech_rank(
        res_ech_t *e,
        const res_diff_t * const rd,
        const len_t lev,
        const deg_t d,
        const int32_t * const rows,
        const len_t nrows,
        const int32_t * const dpos,
        const len_t ncols
        )
{
    len_t r;
    int64_t t;
    int32_t j, jj;

    if (ncols == 0 || nrows == 0) {
        return 0;
    }
    if (res_ech_reserve(e, ncols, nrows)) {
        return -1;
    }

    const res_frame_t * const f = rd->f;
    const ht_t * const ht       = f->ht;
    const res_felt_t * const be = f->lv[lev-1].elts;
    const uint64_t fc           = (uint64_t)rd->fc;
    int64_t rank                = 0;

    for (r = 0; r < nrows; ++r) {
        const res_dpoly_t * const p = rd->d[lev] + rows[r];
        int32_t start = (int32_t)ncols;

        /* scatter the scalar terms; every other term of the row has a ring
         * monomial of positive degree and contributes nothing here */
        for (t = 0; t < p->len; ++t) {
            if (ht->hd[p->mon[t]].deg != 0) {
                continue;
            }
            const int32_t q = p->pos[t];
            if (be[q].hdeg != d) {
                fprintf(ERRSTREAM, "A constant term of frame element "
                        "(%u,%d) sits in degree %d, not %d; the input is "
                        "not homogeneous.\n", (unsigned)lev, rows[r],
                        be[q].hdeg, d);
                return -1;
            }
            j        = dpos[q];
            e->dns[j] = p->cf[t] % fc;
            if (j < start) {
                start = j;
            }
        }

        /* reduce against the pivots found so far */
        int32_t lead = -1;
        for (j = start; j < (int32_t)ncols; ++j) {
            const uint64_t v = e->dns[j] % fc;
            if (v == 0) {
                e->dns[j] = 0;
                continue;
            }
            if (e->piv[j] < 0) {
                lead = j;
                break;
            }
            /* the pivot row is monic and leads in column j, so this also
             * clears dns[j] -- explicitly, to keep the invariant that
             * everything left behind is exactly zero */
            const uint64_t mult = fc - v;
            const int32_t pr    = e->piv[j];
            for (jj = e->beg[pr]; jj < e->end[pr]; ++jj) {
                uint64_t *slot = e->dns + e->idx[jj];
                *slot += mult * (uint64_t)e->val[jj];
                if (*slot >= RES_ACC_LIMIT) {
                    *slot %= fc;
                }
            }
            e->dns[j] = 0;
        }
        if (lead < 0) {
            continue;  /* the row reduced to zero */
        }

        /* a new pivot: normalize and append it to the arena */
        if (res_ech_arena_reserve(e, (int64_t)ncols - lead)) {
            return -1;
        }
        const uint32_t iv = mod_p_inverse_32(
                (int64_t)(e->dns[lead] % fc), (int64_t)fc);
        e->beg[e->nr] = (int32_t)e->ld;
        for (j = lead; j < (int32_t)ncols; ++j) {
            const uint64_t v = e->dns[j] % fc;
            e->dns[j] = 0;
            if (v == 0) {
                continue;
            }
            e->idx[e->ld] = j;
            e->val[e->ld] = (uint32_t)((v * (uint64_t)iv) % fc);
            e->ld++;
        }
        e->end[e->nr] = (int32_t)e->ld;
        e->piv[lead]  = (int32_t)e->nr;
        e->nr++;
        rank++;
    }

    return rank;
}

/* --------------------------------------------------------------------- *
 *  The Betti table
 * --------------------------------------------------------------------- */

res_betti_t *res_betti_new(
        const res_frame_t * const f
        )
{
    len_t i;

    if (f == NULL || f->bad) {
        return NULL;
    }

    res_betti_t *b = (res_betti_t *)calloc(1, sizeof(res_betti_t));
    if (b == NULL) {
        return NULL;
    }
    b->f      = f;
    b->nlv    = f->nlv;
    b->maxdeg = res_frame_max_hdeg(f);
    if (b->maxdeg < 0) {
        free(b);
        return NULL;
    }

    const size_t nd  = (size_t)b->maxdeg + 1;
    const size_t tab = (size_t)b->nlv * nd;

    b->frame = (int32_t *)calloc(tab, sizeof(int32_t));
    b->rank  = (int32_t *)calloc(tab, sizeof(int32_t));
    b->betti = (int32_t *)calloc(tab, sizeof(int32_t));
    b->hilb  = (int32_t *)calloc(nd, sizeof(int32_t));
    if (b->frame == NULL || b->rank == NULL
            || b->betti == NULL || b->hilb == NULL) {
        res_betti_free(&b);
        return NULL;
    }

    if (res_frame_betti(f, b->frame, b->maxdeg) <= 0) {
        res_betti_free(&b);
        return NULL;
    }

    /* Until rank extraction runs, the frame ranks are the best table there
     * is; they are the minimal Betti numbers of a resolution that happens
     * to be minimal, and an upper bound otherwise. */
    memcpy(b->betti, b->frame, tab * sizeof(int32_t));

    /* The Hilbert numerator is the alternating sum of either table, and it
     * needs no field arithmetic -- see the file header. */
    for (i = 0; i < b->nlv; ++i) {
        deg_t d;
        for (d = 0; d <= b->maxdeg; ++d) {
            const int32_t c = b->frame[(size_t)i * nd + d];
            b->hilb[d] = (i & 1) ? b->hilb[d] - c : b->hilb[d] + c;
        }
    }

    return b;
}

void res_betti_free(
        res_betti_t **bp
        )
{
    res_betti_t *b = *bp;

    if (b == NULL) {
        return;
    }
    free(b->frame);
    free(b->rank);
    free(b->betti);
    free(b->hilb);
    free(b);
    *bp = NULL;
}

int res_betti_minimalize(
        res_betti_t *b,
        const res_diff_t * const rd
        )
{
    len_t i, k, lev;
    deg_t d;
    int ret = 1;

    if (b == NULL || b->bad || rd == NULL || rd->bad) {
        return 1;
    }
    if (rd->f != b->f) {
        fprintf(ERRSTREAM, "The differential belongs to a different frame "
                "than the Betti table.\n");
        return 1;
    }

    const res_frame_t * const f = b->f;
    const len_t nlv             = b->nlv;
    const deg_t maxd            = b->maxdeg;
    const size_t nd             = (size_t)maxd + 1;

    /* Per level: the elements listed by degree, the offsets of each degree
     * inside that list, and the position of each element inside its own
     * degree class.  One counting sort per level serves every block. */
    int32_t **bkt = (int32_t **)calloc((unsigned long)nlv, sizeof(int32_t *));
    int32_t **bof = (int32_t **)calloc((unsigned long)nlv, sizeof(int32_t *));
    int32_t **dps = (int32_t **)calloc((unsigned long)nlv, sizeof(int32_t *));
    res_ech_t ech;

    memset(&ech, 0, sizeof(res_ech_t));

    if (bkt == NULL || bof == NULL || dps == NULL) {
        goto cleanup;
    }
    for (lev = 0; lev < nlv; ++lev) {
        const res_level_t * const lv = f->lv + lev;

        bkt[lev] = (int32_t *)malloc(
                (unsigned long)(lv->ld > 0 ? lv->ld : 1) * sizeof(int32_t));
        dps[lev] = (int32_t *)malloc(
                (unsigned long)(lv->ld > 0 ? lv->ld : 1) * sizeof(int32_t));
        bof[lev] = (int32_t *)calloc((unsigned long)maxd + 2, sizeof(int32_t));
        if (bkt[lev] == NULL || bof[lev] == NULL || dps[lev] == NULL) {
            goto cleanup;
        }
        for (k = 0; k < lv->ld; ++k) {
            const deg_t dd = lv->elts[k].hdeg;
            if (dd < 0 || dd > maxd) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) has degree %d, "
                        "outside the frame's own range.\n",
                        (unsigned)lev, (unsigned)k, dd);
                goto cleanup;
            }
            bof[lev][dd + 1]++;
        }
        for (d = 0; d <= maxd; ++d) {
            bof[lev][d+1] += bof[lev][d];
        }
        {
            int32_t *ctr = (int32_t *)calloc(
                    (unsigned long)maxd + 2, sizeof(int32_t));
            if (ctr == NULL) {
                goto cleanup;
            }
            for (k = 0; k < lv->ld; ++k) {
                const deg_t dd = lv->elts[k].hdeg;
                bkt[lev][bof[lev][dd] + ctr[dd]] = (int32_t)k;
                dps[lev][k] = ctr[dd];
                ctr[dd]++;
            }
            free(ctr);
        }
    }

    /* The blocks are mutually independent: (lev, d) reads d_lev and
     * nothing else.  This is the loop a device backend replaces. */
    for (lev = 1; lev < nlv; ++lev) {
        for (d = 0; d <= maxd; ++d) {
            const int32_t rb = bof[lev][d];
            const int32_t nr = bof[lev][d+1] - rb;
            const int32_t nc = bof[lev-1][d+1] - bof[lev-1][d];
            if (nr <= 0 || nc <= 0) {
                continue;
            }
            const int64_t rk = res_ech_rank(&ech, rd, lev, d,
                    bkt[lev] + rb, (len_t)nr, dps[lev-1], (len_t)nc);
            if (rk < 0) {
                goto cleanup;
            }
            b->rank[(size_t)lev * nd + d] = (int32_t)rk;
        }
    }

    /* beta_{i,d} = frame_{i,d} - rank(d_i)_d - rank(d_{i+1})_d */
    for (i = 0; i < nlv; ++i) {
        for (d = 0; d <= maxd; ++d) {
            const int32_t up = (i + 1 < nlv) ? b->rank[(size_t)(i+1)*nd + d] : 0;
            const int32_t v  = b->frame[(size_t)i*nd + d]
                - b->rank[(size_t)i*nd + d] - up;
            if (v < 0) {
                fprintf(ERRSTREAM, "A minimal Betti number is negative at "
                        "level %u degree %d; the rank corrections do not "
                        "fit inside the frame.\n", (unsigned)i, d);
                goto cleanup;
            }
            b->betti[(size_t)i*nd + d] = v;
        }
    }

    b->minimal = 1;
    ret        = 0;

cleanup:
    if (bkt != NULL) {
        for (i = 0; i < nlv; ++i) {
            free(bkt[i]);
        }
        free(bkt);
    }
    if (bof != NULL) {
        for (i = 0; i < nlv; ++i) {
            free(bof[i]);
        }
        free(bof);
    }
    if (dps != NULL) {
        for (i = 0; i < nlv; ++i) {
            free(dps[i]);
        }
        free(dps);
    }
    res_ech_clear(&ech);
    b->bad = ret;

    return ret;
}

/* --------------------------------------------------------------------- *
 *  Invariants read straight off the table
 * --------------------------------------------------------------------- */

int32_t res_betti_pdim(
        const res_betti_t * const b
        )
{
    len_t i;
    deg_t d;
    int32_t p = -1;

    if (b == NULL || b->bad) {
        return -1;
    }
    for (i = 0; i < b->nlv; ++i) {
        for (d = 0; d <= b->maxdeg; ++d) {
            if (b->betti[(size_t)i * (b->maxdeg + 1) + d] != 0) {
                p = (int32_t)i;
                break;
            }
        }
    }
    return p;
}

int32_t res_betti_reg(
        const res_betti_t * const b
        )
{
    len_t i;
    deg_t d;
    int32_t r = -1;
    int seen  = 0;

    if (b == NULL || b->bad) {
        return -1;
    }
    for (i = 0; i < b->nlv; ++i) {
        for (d = 0; d <= b->maxdeg; ++d) {
            if (b->betti[(size_t)i * (b->maxdeg + 1) + d] != 0) {
                const int32_t s = d - (int32_t)i;
                if (!seen || s > r) {
                    r    = s;
                    seen = 1;
                }
            }
        }
    }
    return seen ? r : -1;
}

int res_hilbert_invariants(
        const int32_t * const num,
        const len_t len,
        const len_t nv,
        int32_t *dim,
        int64_t *degree
        )
{
    len_t i;
    len_t c   = 0;
    int64_t s = 0;
    int allzero;

    if (num == NULL || len == 0) {
        return 1;
    }

    int64_t *k = (int64_t *)malloc((unsigned long)len * sizeof(int64_t));
    if (k == NULL) {
        return 1;
    }
    allzero = 1;
    for (i = 0; i < len; ++i) {
        k[i] = num[i];
        s   += k[i];
        if (k[i] != 0) {
            allzero = 0;
        }
    }

    /* The zero module: the numerator is identically zero and there is no
     * order of vanishing to read off. */
    if (allzero) {
        free(k);
        if (dim != NULL) {
            *dim = -1;
        }
        if (degree != NULL) {
            *degree = 0;
        }
        return 0;
    }

    /* K(t) = (1-t)^c * G(t) with G(1) != 0, so dim M = nv - c and the
     * degree of M is G(1).  Dividing by (1-t) is the running sum, and it
     * is exact exactly when K(1) = 0.  The order of vanishing cannot
     * exceed the number of variables: the Hilbert series K/(1-t)^nv of a
     * nonzero module has a pole of order dim M >= 0. */
    while (s == 0) {
        if (c >= nv) {
            fprintf(ERRSTREAM, "The Hilbert numerator vanishes to order "
                    "past the number of variables; the resolution it came "
                    "from cannot have been complete.\n");
            free(k);
            return 1;
        }
        for (i = 1; i < len; ++i) {
            k[i] += k[i-1];
        }
        k[len-1] = 0;  /* this slot held K(1), which is what let us divide */
        s = 0;
        for (i = 0; i < len; ++i) {
            s += k[i];
        }
        c++;
    }

    free(k);
    if (dim != NULL) {
        *dim = (int32_t)nv - (int32_t)c;
    }
    if (degree != NULL) {
        *degree = s;
    }

    return 0;
}
