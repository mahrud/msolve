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

/* The nonminimal differential of a Schreyer resolution.
 *
 * See res.h for what is being computed and why the reduction always
 * reaches zero.  This file is where the field arithmetic finally happens:
 * one Macaulay matrix per (level, degree), with no S-pairs, no pair
 * selection and no criteria -- the frame already decided every row, and
 * every matrix shape is known before a single coefficient is touched.
 *
 * The multiplier bookkeeping is the point.  A plain reduction throws away
 * which reducer was used with which multiple; here that record *is* the
 * answer, so every row carries a second, parallel row over the reducers,
 * exactly as exact_sparse_reduced_echelon_form_sat_ff_32 in la_ff_32.c
 * carries drm alongside drl.  Unlike that one, the columns of the two
 * sides are unrelated: monomials of F_{i-2} on one, generators of F_{i-1}
 * on the other. */

#include "res.h"
#include "../msolve/streams.h"

/* --------------------------------------------------------------------- *
 *  The module order the reduction runs in
 *
 *  A monomial of F_j is a pair (w, q) with w a ring monomial and q an
 *  index into level j.  It stands for the module monomial
 *
 *      w * total(q)  e_{root(q)}   of R^ncomp,
 *
 *  and the order on F_j is the one the strategy's base order on R^ncomp
 *  induces on those, with the index q breaking ties -- larger index
 *  larger, which is Schreyer's tie break.  That is the whole content of
 *  RES_LIFT_SCHREYER: level 0 is literally the base order (total(q) is 1
 *  and root(q) is the component) and every level above carries the order
 *  the one below induces, so the tower is consistent by construction.
 *  Schreyer's theorem needs exactly that consistency.
 *
 *  Both bases are here rather than one, because which is better is an
 *  empirical question and the frames they produce differ by a factor.
 *  The dispatch is on f->strat, constant for a whole computation, so the
 *  branches predict perfectly and the comparison still inlines into the
 *  sorts that call it -- which is why this is a switch and not the
 *  function-pointer vtable io.c uses for the ring orders.
 * --------------------------------------------------------------------- */

/* the component of R^ncomp, in the direction the strategy asks for */
static inline int res_diff_cmp_root(
        const int32_t ra,
        const int32_t rb,
        const int32_t pos
        )
{
    if (ra == rb) {
        return 0;
    }
    if (pos == RES_POS_UP) {
        return ra > rb ? 1 : -1;
    }
    return ra < rb ? 1 : -1;
}

/* The lifted monomials as *module* monomials, ignoring their components:
 * degree first, with the component's degree shift added back, then the
 * ring order's tie break.
 *
 * The shift is what the frame's own table cannot supply.  ev[DEG] of a
 * module monomial includes the shift of its component (res_module.c), and
 * that is the degree the Gröbner basis was compared by; the frame keeps
 * its monomials in a plain ring table, whose ev[DEG] is the ring degree
 * alone.  Level 0 stored the shift as the heft degree of each generator,
 * so it is read back from there. */
static inline int res_diff_cmp_lift(
        const hm_t la,
        const int32_t ra,
        const hm_t lb,
        const int32_t rb,
        const res_frame_t * const f
        )
{
    const ht_t * const ht  = f->ht;
    const exp_t * const ea = ht->ev[la];
    const exp_t * const eb = ht->ev[lb];
    len_t i;

    if (ht->mo == 1) { /* lexicographic, which looks at no degree at all */
        for (i = 1; i <= f->nv; ++i) {
            if (ea[i] != eb[i]) {
                return ea[i] > eb[i] ? 1 : -1;
            }
        }
        return 0;
    }

    /* degree reverse lexicographic */
    const deg_t da = (deg_t)ea[DEG] + f->lv[0].elts[ra-1].hdeg;
    const deg_t db = (deg_t)eb[DEG] + f->lv[0].elts[rb-1].hdeg;
    if (da != db) {
        return da > db ? 1 : -1;
    }
    for (i = f->nv; i > 0; --i) {
        if (ea[i] != eb[i]) {
            return ea[i] < eb[i] ? 1 : -1;
        }
    }

    return 0;
}

static inline int res_diff_cmp_mon(
        const hm_t la,
        const int32_t ra,
        const int32_t pa,
        const hm_t lb,
        const int32_t rb,
        const int32_t pb,
        const res_frame_t * const f
        )
{
    int c;

    if (f->strat.base == RES_MORD_TOP) {
        c = res_diff_cmp_lift(la, ra, lb, rb, f);
        if (c != 0) {
            return c;
        }
        c = res_diff_cmp_root(ra, rb, f->strat.pos);
        if (c != 0) {
            return c;
        }
    } else { /* RES_MORD_POT */
        c = res_diff_cmp_root(ra, rb, f->strat.pos);
        if (c != 0) {
            return c;
        }
        /* Same component, so the two shifts are equal and cancel: the
         * plain ring comparison on the frame's own table already is the
         * module one, and it is the cheaper of the two. */
        c = monomial_cmp(la, lb, f->ht);
        if (c != 0) {
            return c > 0 ? 1 : -1;
        }
    }

    if (pa != pb) {
        return pa > pb ? 1 : -1; /* Schreyer's tie break */
    }

    return 0;
}

/* --------------------------------------------------------------------- *
 *  Columns
 *
 *  One entry per monomial of F_{i-2} seen while building a block, found
 *  again by open addressing on (mon, pos).  lift and root are cached
 *  because the sort below touches them O(log n) times each.
 * --------------------------------------------------------------------- */

typedef struct res_col_t res_col_t;
struct res_col_t
{
    hm_t    mon;
    hm_t    lift;
    int32_t pos;
    int32_t root;
    int32_t piv;   /* reducer whose lead term this is, -1 if none */
};

typedef struct res_cmap_t res_cmap_t;
struct res_cmap_t
{
    res_col_t *col;
    int32_t    ld;
    int32_t    sz;
    int32_t   *map;   /* column index + 1, 0 marks an empty slot */
    int32_t    msz;   /* a power of two                          */
    uint32_t   msk;
    int        lift;  /* compute lift and root; off for d o d = 0,
                       * which only ever needs to tell monomials apart */
};

static inline uint32_t res_cmap_hash(
        const hm_t mon,
        const int32_t pos
        )
{
    return (uint32_t)mon * 2654435761u + (uint32_t)pos * 2246822519u;
}

static void res_cmap_clear(
        res_cmap_t *m
        )
{
    free(m->col);
    free(m->map);
    memset(m, 0, sizeof(res_cmap_t));
}

static int res_cmap_init(
        res_cmap_t *m,
        const int lift
        )
{
    memset(m, 0, sizeof(res_cmap_t));
    m->lift = lift;
    m->sz   = 256;
    m->msz  = 1024;
    m->msk  = (uint32_t)m->msz - 1;
    m->col  = (res_col_t *)malloc((unsigned long)m->sz * sizeof(res_col_t));
    m->map  = (int32_t *)calloc((unsigned long)m->msz, sizeof(int32_t));
    if (m->col == NULL || m->map == NULL) {
        res_cmap_clear(m);
        return 1;
    }
    return 0;
}

static int res_cmap_rehash(
        res_cmap_t *m
        )
{
    int32_t i;

    const int32_t nsz = 2 * m->msz;
    int32_t *nm = (int32_t *)calloc((unsigned long)nsz, sizeof(int32_t));
    if (nm == NULL) {
        return 1;
    }
    const uint32_t msk = (uint32_t)nsz - 1;
    for (i = 0; i < m->ld; ++i) {
        uint32_t k = res_cmap_hash(m->col[i].mon, m->col[i].pos) & msk;
        while (nm[k] != 0) {
            k = (k + 1) & msk;
        }
        nm[k] = i + 1;
    }
    free(m->map);
    m->map = nm;
    m->msz = nsz;
    m->msk = msk;

    return 0;
}

/* Returns the column index of (mon, pos), creating it if it is new, or -1
 * on failure.  *isnew, when not NULL, reports whether it was created. */
static int32_t res_cmap_get(
        res_cmap_t *m,
        const res_frame_t * const f,
        const len_t lev,
        const hm_t mon,
        const int32_t pos,
        exp_t *e,
        int *isnew
        )
{
    if (isnew != NULL) {
        *isnew = 0;
    }

    uint32_t k = res_cmap_hash(mon, pos) & m->msk;
    while (m->map[k] != 0) {
        const int32_t c = m->map[k] - 1;
        if (m->col[c].mon == mon && m->col[c].pos == pos) {
            return c;
        }
        k = (k + 1) & m->msk;
    }

    if (m->ld == m->sz) {
        const int32_t nsz = 2 * m->sz;
        res_col_t *nc = (res_col_t *)realloc(
                m->col, (unsigned long)nsz * sizeof(res_col_t));
        if (nc == NULL) {
            return -1;
        }
        m->col = nc;
        m->sz  = nsz;
    }

    res_col_t *c = m->col + m->ld;
    c->mon  = mon;
    c->pos  = pos;
    c->piv  = -1;
    c->lift = 0;
    c->root = 0;
    if (m->lift) {
        ht_t *ht = f->ht;
        res_frame_ht_reserve(ht, 2);
        c->lift = res_frame_mul(ht, e, mon, f->lv[lev].elts[pos].total);
        if (c->lift == 0) {
            fprintf(ERRSTREAM, "A lifted monomial of the differential "
                    "exceeds the 16-bit exponent limit.\n");
            return -1;
        }
        c->root = f->lv[lev].elts[pos].root;
    }

    m->map[k] = m->ld + 1;
    m->ld++;

    if (isnew != NULL) {
        *isnew = 1;
    }
    if (2 * m->ld >= m->msz && res_cmap_rehash(m)) {
        return -1;
    }

    return m->ld - 1;
}

/* --------------------------------------------------------------------- *
 *  Growable flat row storage
 *
 *  Rows are appended whole and never revisited, so a pair of flat arrays
 *  with one offset per row is all the structure needed.
 * --------------------------------------------------------------------- */

typedef struct res_rows_t res_rows_t;
struct res_rows_t
{
    int32_t  *col;
    uint32_t *cf;
    int64_t  *off;   /* off[r] .. off[r+1] is row r; off[0] == 0 */
    int64_t   nt;
    int64_t   tsz;
    len_t     ld;
    len_t     sz;
};

static void res_rows_clear(
        res_rows_t *r
        )
{
    free(r->col);
    free(r->cf);
    free(r->off);
    memset(r, 0, sizeof(res_rows_t));
}

static int res_rows_init(
        res_rows_t *r
        )
{
    memset(r, 0, sizeof(res_rows_t));
    r->sz  = 64;
    r->tsz = 1024;
    r->col = (int32_t *)malloc((unsigned long)r->tsz * sizeof(int32_t));
    r->cf  = (uint32_t *)malloc((unsigned long)r->tsz * sizeof(uint32_t));
    r->off = (int64_t *)malloc(((unsigned long)r->sz + 1) * sizeof(int64_t));
    if (r->col == NULL || r->cf == NULL || r->off == NULL) {
        res_rows_clear(r);
        return 1;
    }
    r->off[0] = 0;

    return 0;
}

static int res_rows_reserve(
        res_rows_t *r,
        const int64_t n
        )
{
    if (r->nt + n <= r->tsz) {
        return 0;
    }
    int64_t nsz = r->tsz;
    while (r->nt + n > nsz) {
        nsz = 2 * nsz;
    }
    int32_t *nc = (int32_t *)realloc(
            r->col, (unsigned long)nsz * sizeof(int32_t));
    if (nc == NULL) {
        return 1;
    }
    r->col = nc;
    uint32_t *nf = (uint32_t *)realloc(
            r->cf, (unsigned long)nsz * sizeof(uint32_t));
    if (nf == NULL) {
        return 1;
    }
    r->cf  = nf;
    r->tsz = nsz;

    return 0;
}

/* Closes the row started at off[ld]; every term must already be in. */
static int res_rows_close(
        res_rows_t *r
        )
{
    if (r->ld == r->sz) {
        const len_t nsz = 2 * r->sz;
        int64_t *no = (int64_t *)realloc(
                r->off, ((unsigned long)nsz + 1) * sizeof(int64_t));
        if (no == NULL) {
            return 1;
        }
        r->off = no;
        r->sz  = nsz;
    }
    r->ld++;
    r->off[r->ld] = r->nt;

    return 0;
}

/* --------------------------------------------------------------------- *
 *  Construction and teardown
 * --------------------------------------------------------------------- */

static void res_dpoly_clear(
        res_dpoly_t *p
        )
{
    free(p->mon);
    free(p->pos);
    free(p->cf);
    p->mon = NULL;
    p->pos = NULL;
    p->cf  = NULL;
    p->len = 0;
}

static int res_dpoly_alloc(
        res_dpoly_t *p,
        const len_t len
        )
{
    /* Normally this runs on a calloc'd slot and clearing is a no-op.  It
     * is not one when a level is computed a second time, which is what
     * happens after res_diff_compute_thru is abandoned part way through
     * -- an interrupt in the caller -- and then asked for again: blocks
     * below rd->thru are untouched and correct, and everything above it
     * is recomputed over the top of a partial answer. */
    res_dpoly_clear(p);

    p->len = len;
    p->mon = (hm_t *)malloc((unsigned long)(len > 0 ? len : 1) * sizeof(hm_t));
    p->pos = (int32_t *)malloc(
            (unsigned long)(len > 0 ? len : 1) * sizeof(int32_t));
    p->cf  = (uint32_t *)malloc(
            (unsigned long)(len > 0 ? len : 1) * sizeof(uint32_t));
    if (p->mon == NULL || p->pos == NULL || p->cf == NULL) {
        res_dpoly_clear(p);
        return 1;
    }
    return 0;
}

res_diff_t *res_diff_new(
        res_frame_t *f,
        const uint32_t fc
        )
{
    len_t i;

    if (f == NULL || f->bad || f->nlv < 2 || fc == 0) {
        return NULL;
    }
    if (f->gbmap == NULL) {
        fprintf(ERRSTREAM, "The frame carries no Gröbner basis map.\n");
        return NULL;
    }

    res_diff_t *rd = (res_diff_t *)calloc(1, sizeof(res_diff_t));
    if (rd == NULL) {
        return NULL;
    }
    rd->f   = f;
    rd->fc  = fc;
    rd->nlv = f->nlv;
    rd->d   = (res_dpoly_t **)calloc(
            (unsigned long)f->nlv, sizeof(res_dpoly_t *));
    if (rd->d == NULL) {
        free(rd);
        return NULL;
    }
    for (i = 1; i < f->nlv; ++i) {
        rd->d[i] = (res_dpoly_t *)calloc(
                (unsigned long)(f->lv[i].ld > 0 ? f->lv[i].ld : 1),
                sizeof(res_dpoly_t));
        if (rd->d[i] == NULL) {
            res_diff_free(&rd);
            return NULL;
        }
    }

    return rd;
}

void res_diff_free(
        res_diff_t **dp
        )
{
    len_t i, k;
    res_diff_t *rd = *dp;

    if (rd == NULL) {
        return;
    }
    if (rd->d != NULL) {
        for (i = 1; i < rd->nlv; ++i) {
            if (rd->d[i] == NULL) {
                continue;
            }
            for (k = 0; k < rd->f->lv[i].ld; ++k) {
                res_dpoly_clear(rd->d[i] + k);
            }
            free(rd->d[i]);
        }
        free(rd->d);
    }
    free(rd);
    *dp = NULL;
}

/* --------------------------------------------------------------------- *
 *  Level 1: the Gröbner basis itself
 *
 *  The basis lives in the module hash table, whose monomials carry a
 *  component slot and a degree that already includes the component's
 *  shift.  The frame's table is a plain ring table, so each term is split
 *  into a ring monomial there and a level 0 index, and re-degreed.
 * --------------------------------------------------------------------- */

static int res_diff_cmp_term(
        const void *a,
        const void *b,
        void *ctxp
        );

typedef struct res_term_t res_term_t;
struct res_term_t
{
    hm_t     mon;
    hm_t     lift;
    int32_t  pos;
    int32_t  root;
    uint32_t cf;
};

typedef struct res_tctx_t res_tctx_t;
struct res_tctx_t
{
    const res_frame_t *f;
};

/* descending, so the lead term comes first */
static int res_diff_cmp_term(
        const void *a,
        const void *b,
        void *ctxp
        )
{
    const res_tctx_t *ctx  = (const res_tctx_t *)ctxp;
    const res_term_t * const x = (const res_term_t *)a;
    const res_term_t * const y = (const res_term_t *)b;

    return -res_diff_cmp_mon(x->lift, x->root, x->pos,
            y->lift, y->root, y->pos, ctx->f);
}

int res_diff_init(
        res_diff_t *rd,
        const bs_t * const gb,
        const ht_t * const bht,
        const md_t * const md
        )
{
    len_t i, j, k;
    int ret = 1;

    if (rd == NULL || gb == NULL || bht == NULL || md == NULL) {
        return 1;
    }

    res_frame_t *f = rd->f;
    ht_t *ht       = f->ht;
    const len_t nv = f->nv;
    const uint32_t fc = rd->fc;

    if (bht->mord != RES_MORD_POT && bht->mord != RES_MORD_TOP) {
        fprintf(ERRSTREAM, "The differential runs in the Schreyer order "
                "induced by the module order the Gr\u00f6bner basis was "
                "computed in, so that order has to be one res_diff_cmp_mon "
                "implements: position over term or term over position.\n");
        return 1;
    }

    exp_t *e = (exp_t *)calloc((unsigned long)ht->evl, sizeof(exp_t));
    res_term_t *tm = NULL;
    len_t tsz = 0;
    res_tctx_t ctx = {rd->f};

    if (e == NULL) {
        goto cleanup;
    }

    for (k = 0; k < f->lv[1].ld; ++k) {
        const bl_t bi   = (bl_t)f->gbmap[k];
        const hm_t *hm  = gb->hm[bi];
        if (hm == NULL) {
            fprintf(ERRSTREAM, "A frame generator points at the zero "
                    "element of the Gröbner basis.\n");
            goto cleanup;
        }
        const len_t len = hm[LENGTH];

        if (len > tsz) {
            res_term_t *nt = (res_term_t *)realloc(
                    tm, (unsigned long)len * sizeof(res_term_t));
            if (nt == NULL) {
                goto cleanup;
            }
            tm  = nt;
            tsz = len;
        }

        /* the coefficients msolve stored, widened to 32 bits */
        for (j = 0; j < len; ++j) {
            switch (md->ff_bits) {
                case 8:
                    tm[j].cf = (uint32_t)gb->cf_8[hm[COEFFS]][j];
                    break;
                case 16:
                    tm[j].cf = (uint32_t)gb->cf_16[hm[COEFFS]][j];
                    break;
                case 32:
                    tm[j].cf = (uint32_t)gb->cf_32[hm[COEFFS]][j];
                    break;
                default:
                    fprintf(ERRSTREAM, "Unsupported coefficient width %d in "
                            "a module basis.\n", md->ff_bits);
                    goto cleanup;
            }
        }

        res_frame_ht_reserve(ht, (hl_t)len + 1);
        for (j = 0; j < len; ++j) {
            const exp_t * const be = bht->ev[hm[OFFSET+j]];
            for (i = 1; i <= nv; ++i) {
                e[i] = be[i];
            }
            e[0] = (exp_t)res_frame_hdeg(ht, e);
            const int32_t c = (int32_t)be[bht->cpos];
            if (c < 1 || c > (int32_t)f->ncomp) {
                fprintf(ERRSTREAM, "A basis term has component %d, outside "
                        "the ambient free module.\n", c);
                goto cleanup;
            }
            tm[j].mon  = insert_in_hash_table(e, ht);
            tm[j].pos  = c - 1;
            tm[j].lift = tm[j].mon;  /* total() is 1 at level 0 */
            tm[j].root = c;
        }

        sort_r(tm, (unsigned long)len, sizeof(res_term_t),
                res_diff_cmp_term, &ctx);

        if (tm[0].mon != f->lv[1].elts[k].mono
                || tm[0].pos != f->lv[1].elts[k].up) {
            fprintf(ERRSTREAM, "Gröbner basis element %u does not lead with "
                    "the frame's monomial.\n", (unsigned)k);
            goto cleanup;
        }

        /* make it monic: the reduction below subtracts multiples of these
         * rows assuming a unit lead coefficient, and the frame assumes the
         * lead coefficient of every column of the differential is one */
        if (tm[0].cf != 1) {
            const uint32_t inv =
                (uint32_t)mod_p_inverse_32((int64_t)tm[0].cf, (int64_t)fc);
            for (j = 0; j < len; ++j) {
                tm[j].cf = (uint32_t)(
                        ((uint64_t)tm[j].cf * (uint64_t)inv) % (uint64_t)fc);
            }
        }

        res_dpoly_t *p = rd->d[1] + k;
        if (res_dpoly_alloc(p, len)) {
            goto cleanup;
        }
        for (j = 0; j < len; ++j) {
            p->mon[j] = tm[j].mon;
            p->pos[j] = tm[j].pos;
            p->cf[j]  = tm[j].cf;
        }
    }

    ret = 0;
    rd->thru = 1;

cleanup:
    free(e);
    free(tm);
    rd->bad = ret;

    return ret;
}

/* --------------------------------------------------------------------- *
 *  One block: level lev, degree deg
 * --------------------------------------------------------------------- */

typedef struct res_cctx_t res_cctx_t;
struct res_cctx_t
{
    const res_col_t   *col;
    const res_frame_t *f;
};

/* descending, so column 0 carries the largest monomial */
static int res_diff_cmp_col(
        const void *a,
        const void *b,
        void *ctxp
        )
{
    const res_cctx_t *ctx = (const res_cctx_t *)ctxp;
    const res_col_t * const x = ctx->col + ((const int32_t *)a)[0];
    const res_col_t * const y = ctx->col + ((const int32_t *)b)[0];

    return -res_diff_cmp_mon(x->lift, x->root, x->pos,
            y->lift, y->root, y->pos, ctx->f);
}

/* Smallest indexed level mid element sitting over parent and whose own
 * monomial divides w; -1 if there is none.
 *
 * Smallest is not an optimization, it is the correctness condition: the
 * lead term of a column of the differential is only m_k e_{up(k)} if the
 * leading monomial is reduced by something of index below up(k), and the
 * frame guarantees such a reducer exists.  Level mid is stored with up
 * ascending, so the candidates are one contiguous run and index order is
 * scan order. */
static inline int32_t res_diff_reducer(
        const res_frame_t * const f,
        const len_t mid,
        const len_t * const pfirst,
        const len_t * const pcount,
        const int32_t parent,
        const hm_t w
        )
{
    len_t i;

    const res_felt_t * const el = f->lv[mid].elts;
    const len_t b = pfirst[parent];
    const len_t n = pcount[parent];

    for (i = 0; i < n; ++i) {
        if (check_monomial_division(w, el[b+i].mono, f->ht)) {
            return (int32_t)(b + i);
        }
    }

    return -1;
}

static int res_diff_block(
        res_diff_t *rd,
        const len_t lev,
        const int32_t * const rows,
        const len_t nrows,
        const len_t * const pfirst,
        const len_t * const pcount
        )
{
    len_t i, k;
    int32_t c;
    int64_t t;
    int ret = 1;

    res_frame_t *f    = rd->f;
    ht_t *ht          = f->ht;
    const uint32_t fc = rd->fc;
    const len_t mid   = lev - 1;
    const len_t tgt   = lev - 2;

    res_cmap_t cm;
    res_rows_t src, red;
    memset(&cm, 0, sizeof(cm));
    memset(&src, 0, sizeof(src));
    memset(&red, 0, sizeof(red));

    exp_t *e       = NULL;
    int32_t *rq    = NULL;   /* level mid index of each reducer      */
    hm_t *ru       = NULL;   /* multiplier monomial of each reducer  */
    len_t rsz      = 0;
    int32_t *ord   = NULL;
    int32_t *rank  = NULL;
    int32_t *piv   = NULL;
    uint64_t *dns  = NULL;
    uint32_t *mult = NULL;
    res_term_t *tm = NULL;

    e = (exp_t *)calloc((unsigned long)ht->evl, sizeof(exp_t));
    if (e == NULL || res_cmap_init(&cm, 1)
            || res_rows_init(&src) || res_rows_init(&red)) {
        goto cleanup;
    }

    /* --- the rows: v_k = m_k * D_{lev-1}(up(k)) --------------------- */

    for (i = 0; i < nrows; ++i) {
        const res_felt_t * const el = f->lv[lev].elts + rows[i];
        const res_dpoly_t * const p = rd->d[mid] + el->up;

        if (res_rows_reserve(&src, p->len)) {
            goto cleanup;
        }
        for (t = 0; t < p->len; ++t) {
            res_frame_ht_reserve(ht, 2);
            const hm_t w = res_frame_mul(ht, e, el->mono, p->mon[t]);
            if (w == 0) {
                fprintf(ERRSTREAM, "A monomial of the differential exceeds "
                        "the 16-bit exponent limit.\n");
                goto cleanup;
            }
            c = res_cmap_get(&cm, f, tgt, w, p->pos[t], e, NULL);
            if (c < 0) {
                goto cleanup;
            }
            src.col[src.nt] = c;
            src.cf[src.nt]  = p->cf[t];
            src.nt++;
        }
        if (res_rows_close(&src)) {
            goto cleanup;
        }
    }

    /* --- symbolic preprocessing -------------------------------------
     *
     * New columns are appended, so walking the column array forwards is
     * the work queue; nothing is ever revisited. */

    for (c = 0; c < cm.ld; ++c) {
        const hm_t w      = cm.col[c].mon;
        const int32_t par = cm.col[c].pos;
        const int32_t q   =
            res_diff_reducer(f, mid, pfirst, pcount, par, w);
        if (q < 0) {
            continue; /* no reducer; the entry has to cancel by itself */
        }

        if (red.ld == rsz) {
            const len_t nsz = rsz > 0 ? 2 * rsz : 64;
            int32_t *nq = (int32_t *)realloc(
                    rq, (unsigned long)nsz * sizeof(int32_t));
            if (nq == NULL) {
                goto cleanup;
            }
            rq = nq;
            hm_t *nu = (hm_t *)realloc(
                    ru, (unsigned long)nsz * sizeof(hm_t));
            if (nu == NULL) {
                goto cleanup;
            }
            ru  = nu;
            rsz = nsz;
        }

        res_frame_ht_reserve(ht, 2);
        const hm_t u = res_frame_colon(ht, e, w, f->lv[mid].elts[q].mono);
        rq[red.ld] = q;
        ru[red.ld] = u;
        cm.col[c].piv = (int32_t)red.ld;

        const res_dpoly_t * const p = rd->d[mid] + q;
        if (res_rows_reserve(&red, p->len)) {
            goto cleanup;
        }
        for (t = 0; t < p->len; ++t) {
            res_frame_ht_reserve(ht, 2);
            const hm_t v = res_frame_mul(ht, e, u, p->mon[t]);
            if (v == 0) {
                fprintf(ERRSTREAM, "A monomial of the differential exceeds "
                        "the 16-bit exponent limit.\n");
                goto cleanup;
            }
            const int32_t cc =
                res_cmap_get(&cm, f, tgt, v, p->pos[t], e, NULL);
            if (cc < 0) {
                goto cleanup;
            }
            red.col[red.nt] = cc;
            red.cf[red.nt]  = p->cf[t];
            red.nt++;
        }
        if (res_rows_close(&red)) {
            goto cleanup;
        }
    }

    /* --- order the columns ------------------------------------------ */

    const int32_t ncols = cm.ld;

    ord  = (int32_t *)malloc((unsigned long)ncols * sizeof(int32_t));
    rank = (int32_t *)malloc((unsigned long)ncols * sizeof(int32_t));
    piv  = (int32_t *)malloc((unsigned long)ncols * sizeof(int32_t));
    dns  = (uint64_t *)malloc((unsigned long)ncols * sizeof(uint64_t));
    mult = (uint32_t *)malloc(
            (unsigned long)(red.ld > 0 ? red.ld : 1) * sizeof(uint32_t));
    if (ord == NULL || rank == NULL || piv == NULL
            || dns == NULL || mult == NULL) {
        goto cleanup;
    }

    for (c = 0; c < ncols; ++c) {
        ord[c] = c;
    }
    res_cctx_t cctx = {cm.col, f};
    sort_r(ord, (unsigned long)ncols, sizeof(int32_t),
            res_diff_cmp_col, &cctx);
    for (c = 0; c < ncols; ++c) {
        rank[ord[c]] = c;
        piv[c]       = -1;
    }
    for (c = 0; c < ncols; ++c) {
        if (cm.col[c].piv >= 0) {
            piv[rank[c]] = cm.col[c].piv;
        }
    }

    /* --- reduce ------------------------------------------------------ */

    tm = (res_term_t *)malloc(
            (unsigned long)(red.ld + 1) * sizeof(res_term_t));
    if (tm == NULL) {
        goto cleanup;
    }
    res_tctx_t tctx = {f};

    for (i = 0; i < nrows; ++i) {
        const res_felt_t * const el = f->lv[lev].elts + rows[i];

        memset(dns, 0, (unsigned long)ncols * sizeof(uint64_t));
        memset(mult, 0, (unsigned long)red.ld * sizeof(uint32_t));

        for (t = src.off[i]; t < src.off[i+1]; ++t) {
            dns[rank[src.col[t]]] = src.cf[t];
        }

        for (c = 0; c < ncols; ++c) {
            if (dns[c] == 0) {
                continue;
            }
            const int32_t r = piv[c];
            if (r < 0) {
                fprintf(ERRSTREAM, "The differential of frame element "
                        "(%u,%u) does not reduce to zero: a monomial has no "
                        "reducer.\n", (unsigned)lev, (unsigned)rows[i]);
                goto cleanup;
            }
            const uint64_t mul = dns[c];
            mult[r] = (uint32_t)(((uint64_t)mult[r] + mul) % (uint64_t)fc);
            for (t = red.off[r]; t < red.off[r+1]; ++t) {
                const int32_t cc = rank[red.col[t]];
                dns[cc] = (dns[cc]
                        + mul * (uint64_t)(fc - red.cf[t])) % (uint64_t)fc;
            }
            if (dns[c] != 0) {
                fprintf(ERRSTREAM, "A reducer of the differential is not "
                        "monic.\n");
                goto cleanup;
            }
        }

        /* --- read the multipliers off as the column ------------------ */

        len_t nt = 0;
        tm[nt].mon  = el->mono;
        tm[nt].pos  = el->up;
        tm[nt].cf   = 1;
        tm[nt].root = f->lv[mid].elts[el->up].root;
        res_frame_ht_reserve(ht, 2);
        tm[nt].lift = res_frame_mul(
                ht, e, el->mono, f->lv[mid].elts[el->up].total);
        if (tm[nt].lift == 0) {
            goto cleanup;
        }
        nt++;

        for (k = 0; k < red.ld; ++k) {
            if (mult[k] == 0) {
                continue;
            }
            tm[nt].mon  = ru[k];
            tm[nt].pos  = rq[k];
            tm[nt].cf   = fc - mult[k];
            tm[nt].root = f->lv[mid].elts[rq[k]].root;
            res_frame_ht_reserve(ht, 2);
            tm[nt].lift = res_frame_mul(
                    ht, e, ru[k], f->lv[mid].elts[rq[k]].total);
            if (tm[nt].lift == 0) {
                goto cleanup;
            }
            nt++;
        }

        sort_r(tm, (unsigned long)nt, sizeof(res_term_t),
                res_diff_cmp_term, &tctx);

        if (tm[0].mon != el->mono || tm[0].pos != el->up) {
            fprintf(ERRSTREAM, "The differential of frame element (%u,%u) "
                    "does not lead with the frame's monomial; the reducer "
                    "selection is wrong.\n",
                    (unsigned)lev, (unsigned)rows[i]);
            goto cleanup;
        }

        res_dpoly_t *p = rd->d[lev] + rows[i];
        if (res_dpoly_alloc(p, nt)) {
            goto cleanup;
        }
        for (k = 0; k < nt; ++k) {
            p->mon[k] = tm[k].mon;
            p->pos[k] = tm[k].pos;
            p->cf[k]  = tm[k].cf;
        }
    }

    ret = 0;

cleanup:
    free(e);
    free(rq);
    free(ru);
    free(ord);
    free(rank);
    free(piv);
    free(dns);
    free(mult);
    free(tm);
    res_cmap_clear(&cm);
    res_rows_clear(&src);
    res_rows_clear(&red);

    return ret;
}

/* --------------------------------------------------------------------- *
 *  The driver
 * --------------------------------------------------------------------- */

int res_diff_compute(
        res_diff_t *rd
        )
{
    if (rd == NULL || rd->bad) {
        return 1;
    }

    return res_diff_compute_thru(rd, rd->nlv - 1);
}

int res_diff_compute_thru(
        res_diff_t *rd,
        const len_t maxlev
        )
{
    len_t i, k, lev, lo, hi;
    deg_t d;
    int ret = 1;

    if (rd == NULL || rd->bad) {
        return 1;
    }

    res_frame_t *f   = rd->f;
    const len_t nlv  = f->nlv;
    const deg_t maxd = res_frame_max_hdeg(f);

    /* Level 1 is res_diff_init's, so thru is at least 1 by the time this
     * can run; anything at or below what is already filled in is a no-op,
     * which is what makes repeated calls cheap. */
    lo = rd->thru + 1 > 2 ? rd->thru + 1 : 2;
    hi = maxlev < nlv - 1 ? maxlev : nlv - 1;
    if (hi < lo) {
        return 0;
    }

    /* per level: the elements bucketed by degree, and the level below it
     * bucketed by parent.  Both are counting sorts over data the frame
     * already stores in the right order, done once rather than once per
     * block. */
    int32_t **bkt = (int32_t **)calloc((unsigned long)nlv, sizeof(int32_t *));
    int32_t **bof = (int32_t **)calloc((unsigned long)nlv, sizeof(int32_t *));
    len_t **pfst  = (len_t **)calloc((unsigned long)nlv, sizeof(len_t *));
    len_t **pcnt  = (len_t **)calloc((unsigned long)nlv, sizeof(len_t *));

    if (bkt == NULL || bof == NULL || pfst == NULL || pcnt == NULL) {
        goto cleanup;
    }

    for (lev = lo; lev <= hi; ++lev) {
        const res_level_t * const lv = f->lv + lev;

        bkt[lev] = (int32_t *)malloc(
                (unsigned long)(lv->ld > 0 ? lv->ld : 1) * sizeof(int32_t));
        bof[lev] = (int32_t *)calloc(
                (unsigned long)maxd + 2, sizeof(int32_t));
        if (bkt[lev] == NULL || bof[lev] == NULL) {
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
                ctr[dd]++;
            }
            free(ctr);
        }

        /* level lev-1 grouped by its parent at level lev-2 */
        const res_level_t * const ml = f->lv + lev - 1;
        const len_t np = f->lv[lev-2].ld;
        pfst[lev] = (len_t *)calloc(
                (unsigned long)(np > 0 ? np : 1), sizeof(len_t));
        pcnt[lev] = (len_t *)calloc(
                (unsigned long)(np > 0 ? np : 1), sizeof(len_t));
        if (pfst[lev] == NULL || pcnt[lev] == NULL) {
            goto cleanup;
        }
        for (k = 0; k < ml->ld; ++k) {
            const int32_t up = ml->elts[k].up;
            if (up < 0 || up >= (int32_t)np) {
                fprintf(ERRSTREAM, "Frame element (%u,%u) has an out of "
                        "range parent.\n", (unsigned)(lev-1), (unsigned)k);
                goto cleanup;
            }
            if (pcnt[lev][up] == 0) {
                pfst[lev][up] = k;
            }
            pcnt[lev][up]++;
        }
    }

    /* Degree ascending, then level ascending.  Level lev in degree d
     * reduces against level lev-1 in degrees up to and including d, and
     * (d, lev-1) is exactly the step before (d, lev). */
    for (d = 0; d <= maxd; ++d) {
        for (lev = lo; lev <= hi; ++lev) {
            const int32_t b = bof[lev][d];
            const int32_t n = bof[lev][d+1] - b;
            if (n <= 0) {
                continue;
            }
            if (res_diff_block(rd, lev, bkt[lev] + b, (len_t)n,
                        pfst[lev], pcnt[lev])) {
                goto cleanup;
            }
        }
    }

    ret = 0;
    rd->thru = hi;

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
    if (pfst != NULL) {
        for (i = 0; i < nlv; ++i) {
            free(pfst[i]);
        }
        free(pfst);
    }
    if (pcnt != NULL) {
        for (i = 0; i < nlv; ++i) {
            free(pcnt[i]);
        }
        free(pcnt);
    }
    rd->bad = ret;

    return ret;
}

/* --------------------------------------------------------------------- *
 *  d o d = 0
 *
 *  The one check that needs nothing but the answer itself: expand
 *  D_{i-1}(D_i(k)) term by term into a hash map keyed on the monomials of
 *  F_{i-2} and require every coefficient to vanish.  Everything else --
 *  the frame, Schreyer's theorem, the reducer selection -- only enters
 *  through the result being wrong if any of it is.
 *
 *  It is also a polynomial multiplication per column, which on a real
 *  input costs several times the resolution itself -- measured at six
 *  times on six generic cubics in six variables -- so unlike
 *  res_frame_verify it does not run on every computation.  The cheap
 *  half, the lead term of every column against the frame, always does.
 * --------------------------------------------------------------------- */

int res_diff_verify(
        const res_diff_t * const rd,
        const int deep
        )
{
    len_t lev, k;
    int64_t s, t;
    int bad = 0;

    if (rd == NULL || rd->bad) {
        return 1;
    }

    res_frame_t *f    = rd->f;
    ht_t *ht          = f->ht;
    const uint32_t fc = rd->fc;

    res_cmap_t cm;
    exp_t *e     = (exp_t *)calloc((unsigned long)ht->evl, sizeof(exp_t));
    uint64_t *ac = NULL;
    int32_t acsz = 0;

    if (e == NULL || res_cmap_init(&cm, 0)) {
        free(e);
        return 1;
    }

    /* only what has actually been filled in: a level past thru is empty by
     * construction, not wrong, and every one of its columns would be
     * reported as a violation */
    for (lev = 1; lev <= rd->thru && bad <= 8; ++lev) {
        for (k = 0; k < f->lv[lev].ld && bad <= 8; ++k) {
            const res_dpoly_t * const p = rd->d[lev] + k;

            if (p->len == 0) {
                fprintf(ERRSTREAM, "The differential of frame element "
                        "(%u,%u) is empty.\n", (unsigned)lev, (unsigned)k);
                bad++;
                continue;
            }
            if (p->cf[0] != 1 || p->mon[0] != f->lv[lev].elts[k].mono
                    || p->pos[0] != f->lv[lev].elts[k].up) {
                fprintf(ERRSTREAM, "The differential of frame element "
                        "(%u,%u) does not lead with the frame's monomial.\n",
                        (unsigned)lev, (unsigned)k);
                bad++;
            }
            if (lev == 1 || !deep) {
                continue;
            }

            /* expand D_{lev-1}(D_lev(k)) */
            cm.ld = 0;
            memset(cm.map, 0, (unsigned long)cm.msz * sizeof(int32_t));

            int fail = 0;
            for (t = 0; t < p->len; ++t) {
                const res_dpoly_t * const q = rd->d[lev-1] + p->pos[t];
                for (s = 0; s < q->len; ++s) {
                    res_frame_ht_reserve(ht, 2);
                    const hm_t w = res_frame_mul(ht, e, p->mon[t], q->mon[s]);
                    if (w == 0) {
                        fail = 1;
                        break;
                    }
                    const int32_t c = res_cmap_get(
                            &cm, f, lev-2, w, q->pos[s], e, NULL);
                    if (c < 0) {
                        fail = 1;
                        break;
                    }
                    if (c >= acsz) {
                        const int32_t nsz = c + 64;
                        uint64_t *na = (uint64_t *)realloc(
                                ac, (unsigned long)nsz * sizeof(uint64_t));
                        if (na == NULL) {
                            fail = 1;
                            break;
                        }
                        memset(na + acsz, 0,
                                (unsigned long)(nsz - acsz) * sizeof(uint64_t));
                        ac   = na;
                        acsz = nsz;
                    }
                    ac[c] = (ac[c] + (uint64_t)p->cf[t] * (uint64_t)q->cf[s])
                        % (uint64_t)fc;
                }
                if (fail) {
                    break;
                }
            }
            int32_t c;
            if (fail) {
                bad++;
            } else {
                for (c = 0; c < cm.ld; ++c) {
                    if (ac[c] != 0) {
                        fprintf(ERRSTREAM, "d o d is not zero on frame "
                                "element (%u,%u).\n",
                                (unsigned)lev, (unsigned)k);
                        bad++;
                        break;
                    }
                }
            }
            /* the accumulator is reused, so clear exactly what was used */
            if (ac != NULL && cm.ld > 0) {
                memset(ac, 0, (unsigned long)cm.ld * sizeof(uint64_t));
            }
        }
    }

    free(e);
    free(ac);
    res_cmap_clear(&cm);

    return bad;
}
