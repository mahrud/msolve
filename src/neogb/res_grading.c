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

#include "res.h"
#include "../msolve/streams.h"

/* --------------------------------------------------------------------- *
 *  Default arithmetic for A = Z^r (+) Z/t_1 (+) ... (+) Z/t_nt
 *
 *  Two implementations of each operation: a generic one, and one
 *  specialized to the standard grading r == 1, nt == 0 so that the
 *  singly graded path pays nothing for the generality.  res_dgrp_new
 *  installs whichever applies.
 * --------------------------------------------------------------------- */

static inline int32_t res_mod_torsion(
        const int32_t v,
        const int32_t t
        )
{
    int32_t w = v % t;
    /* C99 % follows the sign of the dividend; we want [0,t) */
    return w < 0 ? w + t : w;
}

static void res_dgrp_add_gen(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a,
        const res_deg_t b
        )
{
    len_t i;
    const len_t r = g->r;

    for (i = 0; i < r; ++i) {
        d.e[i] = a.e[i] + b.e[i];
    }
    for (i = 0; i < g->nt; ++i) {
        d.e[r+i] = res_mod_torsion(a.e[r+i] + b.e[r+i], g->tord[i]);
    }
}

static void res_dgrp_sub_gen(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a,
        const res_deg_t b
        )
{
    len_t i;
    const len_t r = g->r;

    for (i = 0; i < r; ++i) {
        d.e[i] = a.e[i] - b.e[i];
    }
    for (i = 0; i < g->nt; ++i) {
        d.e[r+i] = res_mod_torsion(a.e[r+i] - b.e[r+i], g->tord[i]);
    }
}

/* Total order on A used only for bucketing and sorting: lexicographic on
 * the free part, then on the torsion residues.  It has no algebraic
 * meaning and need only be consistent. */
static int res_dgrp_cmp_gen(
        const res_dgrp_t *g,
        const res_deg_t a,
        const res_deg_t b
        )
{
    len_t i;
    const len_t len = g->len;

    for (i = 0; i < len; ++i) {
        if (a.e[i] != b.e[i]) {
            return a.e[i] < b.e[i] ? -1 : 1;
        }
    }
    return 0;
}

static hl_t res_dgrp_hash_gen(
        const res_dgrp_t *g,
        const res_deg_t a
        )
{
    len_t i;
    /* FNV-1a over the degree slots */
    hl_t h = 14695981039346656037UL;

    for (i = 0; i < g->len; ++i) {
        h ^= (hl_t)(uint32_t)a.e[i];
        h *= 1099511628211UL;
    }
    return h;
}

static deg_t res_dgrp_heft_gen(
        const res_dgrp_t *g,
        const res_deg_t a
        )
{
    len_t i;
    int64_t h = 0;

    /* torsion carries no heft */
    for (i = 0; i < g->r; ++i) {
        h += (int64_t)g->heft[i] * (int64_t)a.e[i];
    }
    return (deg_t)h;
}

/* --- specializations for the standard grading ------------------------ */

static void res_dgrp_add_z1(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a,
        const res_deg_t b
        )
{
    (void)g;
    d.e[0] = a.e[0] + b.e[0];
}

static void res_dgrp_sub_z1(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a,
        const res_deg_t b
        )
{
    (void)g;
    d.e[0] = a.e[0] - b.e[0];
}

static int res_dgrp_cmp_z1(
        const res_dgrp_t *g,
        const res_deg_t a,
        const res_deg_t b
        )
{
    (void)g;
    if (a.e[0] == b.e[0]) {
        return 0;
    }
    return a.e[0] < b.e[0] ? -1 : 1;
}

static hl_t res_dgrp_hash_z1(
        const res_dgrp_t *g,
        const res_deg_t a
        )
{
    (void)g;
    return (hl_t)(uint32_t)a.e[0];
}

static deg_t res_dgrp_heft_z1(
        const res_dgrp_t *g,
        const res_deg_t a
        )
{
    return (deg_t)(g->heft[0] * a.e[0]);
}

/* --------------------------------------------------------------------- *
 *  Construction
 * --------------------------------------------------------------------- */

res_dgrp_t *res_dgrp_new(
        const len_t r,
        const len_t nt,
        const int32_t *tord,
        const len_t nv,
        const int32_t *degs,
        const int32_t *heft
        )
{
    len_t i, j;

    if (r == 0 && nt == 0) {
        fprintf(ERRSTREAM, "Grading group must be nontrivial.\n");
        return NULL;
    }
    if (nv == 0) {
        fprintf(ERRSTREAM, "Grading needs at least one variable.\n");
        return NULL;
    }

    res_dgrp_t *g = (res_dgrp_t *)calloc(1, sizeof(res_dgrp_t));
    if (g == NULL) {
        return NULL;
    }

    g->r    = r;
    g->nt   = nt;
    g->len  = r + nt;
    g->nv   = nv;

    g->tord  = (int32_t *)calloc((unsigned long)(nt > 0 ? nt : 1), sizeof(int32_t));
    g->heft  = (int32_t *)calloc((unsigned long)(r > 0 ? r : 1), sizeof(int32_t));
    g->dmat  = (int32_t *)calloc((unsigned long)g->len * nv, sizeof(int32_t));
    g->vhdeg = (deg_t *)calloc((unsigned long)nv, sizeof(deg_t));

    if (g->tord == NULL || g->heft == NULL
            || g->dmat == NULL || g->vhdeg == NULL) {
        res_dgrp_free(&g);
        return NULL;
    }

    for (i = 0; i < nt; ++i) {
        if (tord[i] < 2) {
            fprintf(ERRSTREAM,
                    "Torsion factor %u has order %d, must be at least 2.\n",
                    (unsigned)i, tord[i]);
            res_dgrp_free(&g);
            return NULL;
        }
        g->tord[i] = tord[i];
    }

    for (i = 0; i < r; ++i) {
        g->heft[i] = heft != NULL ? heft[i] : 1;
    }

    memcpy(g->dmat, degs, (unsigned long)g->len * nv * sizeof(int32_t));
    /* normalize the torsion part of every variable degree once, so that
     * later additions only ever have to reduce a single carry */
    for (j = 0; j < nv; ++j) {
        int32_t *col = g->dmat + (unsigned long)j * g->len;
        for (i = 0; i < nt; ++i) {
            col[r+i] = res_mod_torsion(col[r+i], g->tord[i]);
        }
    }

    /* The heft must be strictly positive on every variable: this is what
     * makes ev[DEG] a legitimate scheduling degree and guarantees that the
     * degree-by-degree drivers terminate. */
    for (j = 0; j < nv; ++j) {
        const int32_t *col = g->dmat + (unsigned long)j * g->len;
        int64_t h = 0;
        for (i = 0; i < r; ++i) {
            h += (int64_t)g->heft[i] * (int64_t)col[i];
        }
        if (h <= 0) {
            fprintf(ERRSTREAM,
                    "Heft vector is not strictly positive on variable %u "
                    "(heft . deg = %ld); a graded computation is not "
                    "possible with this grading.\n", (unsigned)j, (long)h);
            res_dgrp_free(&g);
            return NULL;
        }
        /* The heft degree lands in ev[DEG], which is 16 bits wide, so a
         * variable heavier than that can never occur in a monomial the hash
         * table can hold.  Bounding it here keeps every later heft
         * computation inside deg_t rather than relying on the callers to
         * notice a truncation. */
        if (h > UINT16_MAX) {
            fprintf(ERRSTREAM,
                    "Variable %u has heft degree %ld, past the 16-bit degree "
                    "the exponent table holds.\n", (unsigned)j, (long)h);
            res_dgrp_free(&g);
            return NULL;
        }
        g->vhdeg[j] = (deg_t)h;
    }

    g->simple = (r == 1 && nt == 0);

    if (g->simple) {
        g->add     = res_dgrp_add_z1;
        g->sub     = res_dgrp_sub_z1;
        g->cmp     = res_dgrp_cmp_z1;
        g->hash    = res_dgrp_hash_z1;
        g->heft_of = res_dgrp_heft_z1;
    } else {
        g->add     = res_dgrp_add_gen;
        g->sub     = res_dgrp_sub_gen;
        g->cmp     = res_dgrp_cmp_gen;
        g->hash    = res_dgrp_hash_gen;
        g->heft_of = res_dgrp_heft_gen;
    }

    return g;
}

int32_t res_grading_len(
        const res_grading_t * const g
        )
{
    if (g == NULL) {
        return 1;
    }
    return g->r + g->nt;
}

res_dgrp_t *res_dgrp_of_grading(
        const res_grading_t * const grading,
        const len_t nv
        )
{
    if (grading == NULL) {
        return res_dgrp_new_standard(nv);
    }
    if (grading->r < 0 || grading->nt < 0) {
        fprintf(ERRSTREAM, "A grading group cannot have a negative number "
                "of free or torsion factors.\n");
        return NULL;
    }
    if (grading->r + grading->nt > RES_MTAB_MAXLEN) {
        fprintf(ERRSTREAM, "A grading group of %d factors is past the %d "
                "this engine reports.\n",
                grading->r + grading->nt, RES_MTAB_MAXLEN);
        return NULL;
    }
    if (grading->degs == NULL) {
        fprintf(ERRSTREAM, "A grading needs a degree for every variable.\n");
        return NULL;
    }
    if (grading->nt > 0 && grading->tord == NULL) {
        fprintf(ERRSTREAM, "A grading with %d torsion factors needs their "
                "orders.\n", grading->nt);
        return NULL;
    }

    return res_dgrp_new(grading->r, grading->nt, grading->tord, nv,
            grading->degs, grading->heft);
}

res_dgrp_t *res_dgrp_new_standard(
        const len_t nv
        )
{
    len_t j;
    res_dgrp_t *g = NULL;

    int32_t *degs = (int32_t *)malloc((unsigned long)nv * sizeof(int32_t));
    if (degs == NULL) {
        return NULL;
    }
    for (j = 0; j < nv; ++j) {
        degs[j] = 1;
    }
    g = res_dgrp_new(1, 0, NULL, nv, degs, NULL);
    free(degs);

    return g;
}

void res_dgrp_free(
        res_dgrp_t **gp
        )
{
    res_dgrp_t *g = *gp;

    if (g == NULL) {
        return;
    }
    free(g->tord);
    free(g->heft);
    free(g->dmat);
    free(g->vhdeg);
    free(g);
    *gp = NULL;
}

/* --------------------------------------------------------------------- *
 *  Degree pools
 * --------------------------------------------------------------------- */

res_dpool_t *res_dpool_new(
        const res_dgrp_t *grp,
        const hl_t sz
        )
{
    res_dpool_t *p = (res_dpool_t *)calloc(1, sizeof(res_dpool_t));
    if (p == NULL) {
        return NULL;
    }

    p->grp  = grp;
    p->ld   = 0;
    p->sz   = sz > 0 ? sz : 16;
    p->data = (int32_t *)calloc(
            (unsigned long)p->sz * grp->len, sizeof(int32_t));
    if (p->data == NULL) {
        free(p);
        return NULL;
    }
    return p;
}

void res_dpool_free(
        res_dpool_t **pp
        )
{
    res_dpool_t *p = *pp;

    if (p == NULL) {
        return;
    }
    free(p->data);
    free(p);
    *pp = NULL;
}

res_deg_t res_dpool_at(
        const res_dpool_t *p,
        const hl_t i
        )
{
    res_deg_t d;
    d.e = p->data + i * p->grp->len;
    return d;
}

res_deg_t res_dpool_push(
        res_dpool_t *p,
        hl_t *idx
        )
{
    res_deg_t d;
    const len_t len = p->grp->len;

    if (p->ld >= p->sz) {
        const hl_t nsz = 2 * p->sz;
        int32_t *nd = (int32_t *)realloc(
                p->data, (unsigned long)nsz * len * sizeof(int32_t));
        if (nd == NULL) {
            fprintf(ERRSTREAM, "Could not enlarge degree pool, "
                    "segmentation fault will follow.\n");
            d.e = NULL;
            return d;
        }
        memset(nd + p->sz * len, 0,
                (unsigned long)(nsz - p->sz) * len * sizeof(int32_t));
        p->data = nd;
        p->sz   = nsz;
    }

    d.e = p->data + p->ld * len;
    memset(d.e, 0, (unsigned long)len * sizeof(int32_t));
    if (idx != NULL) {
        *idx = p->ld;
    }
    p->ld++;

    return d;
}

/* --------------------------------------------------------------------- *
 *  Multidegree buckets
 *
 *  An open addressed hash set of degrees, storing them in a pool and the
 *  map only in indices, so a rehash moves no degree data.  Both the hash
 *  and the equality test go through the group's vtable, which is what lets
 *  an external group backend be bucketed without this file knowing
 *  anything about its representation.
 * --------------------------------------------------------------------- */

res_dbkt_t *res_dbkt_new(
        const res_dgrp_t *grp,
        const hl_t sz
        )
{
    hl_t m = 16;

    if (grp == NULL) {
        return NULL;
    }
    while (m < 2 * (sz > 0 ? sz : 1)) {
        m *= 2;
    }

    res_dbkt_t *b = (res_dbkt_t *)calloc(1, sizeof(res_dbkt_t));
    if (b == NULL) {
        return NULL;
    }
    b->grp  = grp;
    b->msz  = m;
    b->ld   = 0;
    b->pool = res_dpool_new(grp, sz > 0 ? sz : 16);
    b->map  = (hl_t *)calloc((unsigned long)m, sizeof(hl_t));
    if (b->pool == NULL || b->map == NULL) {
        res_dbkt_free(&b);
        return NULL;
    }

    return b;
}

void res_dbkt_free(
        res_dbkt_t **bp
        )
{
    res_dbkt_t *b = *bp;

    if (b == NULL) {
        return;
    }
    res_dpool_free(&b->pool);
    free(b->map);
    free(b);
    *bp = NULL;
}

res_deg_t res_dbkt_at(
        const res_dbkt_t * const b,
        const hl_t i
        )
{
    return res_dpool_at(b->pool, i);
}

/* The probe sequence, shared by find and insert: returns the map slot
 * holding a, or the first empty slot if a is absent. */
static inline hl_t res_dbkt_slot(
        const res_dbkt_t * const b,
        const res_deg_t a
        )
{
    const res_dgrp_t * const g = b->grp;
    const hl_t mask = b->msz - 1;
    hl_t k = g->hash(g, a) & mask;
    hl_t i;

    for (i = 0; i < b->msz; ++i) {
        const hl_t e = b->map[k];
        if (e == 0) {
            return k;
        }
        if (g->cmp(g, res_dpool_at(b->pool, e - 1), a) == 0) {
            return k;
        }
        k = (k + i + 1) & mask;
    }
    return k; /* unreachable: the map is never allowed to fill up */
}

hl_t res_dbkt_find(
        const res_dbkt_t * const b,
        const res_deg_t a
        )
{
    const hl_t e = b->map[res_dbkt_slot(b, a)];

    return e == 0 ? (hl_t)-1 : e - 1;
}

static int res_dbkt_rehash(
        res_dbkt_t *b
        )
{
    hl_t i;
    const hl_t nsz = 2 * b->msz;

    hl_t *nm = (hl_t *)calloc((unsigned long)nsz, sizeof(hl_t));
    if (nm == NULL) {
        fprintf(ERRSTREAM, "Could not enlarge the multidegree buckets.\n");
        return 1;
    }
    free(b->map);
    b->map = nm;
    b->msz = nsz;

    for (i = 0; i < b->ld; ++i) {
        b->map[res_dbkt_slot(b, res_dpool_at(b->pool, i))] = i + 1;
    }

    return 0;
}

hl_t res_dbkt_insert(
        res_dbkt_t *b,
        const res_deg_t a
        )
{
    hl_t k = res_dbkt_slot(b, a);

    if (b->map[k] != 0) {
        return b->map[k] - 1;
    }
    /* keep the map at most half full, which is what bounds the probe */
    if (2 * (b->ld + 1) > b->msz) {
        if (res_dbkt_rehash(b)) {
            return (hl_t)-1;
        }
        k = res_dbkt_slot(b, a);
    }

    hl_t idx = 0;
    res_deg_t d = res_dpool_push(b->pool, &idx);
    if (d.e == NULL) {
        return (hl_t)-1;
    }
    res_deg_set(b->grp, d, a);
    b->map[k] = idx + 1;
    b->ld++;

    return idx;
}

/* The sort is over indices rather than degrees, so the pool is permuted
 * once at the end instead of being swapped through. */
typedef struct res_dbkt_sctx_t res_dbkt_sctx_t;
struct res_dbkt_sctx_t
{
    const res_dbkt_t *b;
};

static int res_dbkt_cmp_idx(
        const void *a,
        const void *b,
        void *ctx
        )
{
    const res_dbkt_sctx_t * const c = (const res_dbkt_sctx_t *)ctx;
    const hl_t ia = *(const hl_t *)a;
    const hl_t ib = *(const hl_t *)b;
    const res_dgrp_t * const g = c->b->grp;

    return g->cmp(g, res_dpool_at(c->b->pool, ia), res_dpool_at(c->b->pool, ib));
}

int res_dbkt_sort(
        res_dbkt_t *b,
        hl_t *perm
        )
{
    hl_t i;

    if (b == NULL || b->ld == 0) {
        return 0;
    }

    const len_t len = b->grp->len;

    hl_t *ord = (hl_t *)malloc((unsigned long)b->ld * sizeof(hl_t));
    int32_t *nd = (int32_t *)malloc(
            (unsigned long)b->ld * (unsigned long)len * sizeof(int32_t));
    if (ord == NULL || nd == NULL) {
        free(ord);
        free(nd);
        fprintf(ERRSTREAM, "Could not sort the multidegree buckets.\n");
        return 1;
    }
    for (i = 0; i < b->ld; ++i) {
        ord[i] = i;
    }

    res_dbkt_sctx_t ctx = {b};
    sort_r(ord, (unsigned long)b->ld, sizeof(hl_t), res_dbkt_cmp_idx, &ctx);

    /* ord[u] is the old bucket that becomes the new bucket u, so perm --
     * which the caller applies to a table it already built -- is its
     * inverse. */
    for (i = 0; i < b->ld; ++i) {
        memcpy(nd + (unsigned long)i * len,
                b->pool->data + (unsigned long)ord[i] * len,
                (unsigned long)len * sizeof(int32_t));
        perm[ord[i]] = i;
    }
    memcpy(b->pool->data, nd, (unsigned long)b->ld * len * sizeof(int32_t));
    free(nd);
    free(ord);

    /* the map still points at the old buckets */
    memset(b->map, 0, (unsigned long)b->msz * sizeof(hl_t));
    for (i = 0; i < b->ld; ++i) {
        b->map[res_dbkt_slot(b, res_dpool_at(b->pool, i))] = i + 1;
    }

    return 0;
}

/* --------------------------------------------------------------------- *
 *  Degree arithmetic used by the frame and the differential
 * --------------------------------------------------------------------- */

void res_deg_zero(
        const res_dgrp_t *g,
        res_deg_t d
        )
{
    memset(d.e, 0, (unsigned long)g->len * sizeof(int32_t));
}

void res_deg_set(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a
        )
{
    memcpy(d.e, a.e, (unsigned long)g->len * sizeof(int32_t));
}

void res_deg_of_exponents(
        const res_dgrp_t *g,
        res_deg_t d,
        const exp_t *exps
        )
{
    len_t i, j;
    const len_t r   = g->r;
    const len_t nt  = g->nt;
    const len_t len = g->len;
    const len_t nv  = g->nv;

    if (g->simple) {
        int64_t acc = 0;
        const int32_t * const dm = g->dmat;
        for (j = 0; j < nv; ++j) {
            acc += (int64_t)exps[j] * (int64_t)dm[j];
        }
        d.e[0] = (int32_t)acc;
        return;
    }

    /* accumulate the free part in 64 bits, reduce torsion at the end */
    int64_t *acc = (int64_t *)calloc((unsigned long)r, sizeof(int64_t));
    if (acc == NULL) {
        fprintf(ERRSTREAM, "Out of memory computing a multidegree.\n");
        return;
    }
    memset(d.e, 0, (unsigned long)len * sizeof(int32_t));

    for (j = 0; j < nv; ++j) {
        const int32_t e = (int32_t)exps[j];
        if (e == 0) {
            continue;
        }
        const int32_t * const col = g->dmat + (unsigned long)j * len;
        for (i = 0; i < r; ++i) {
            acc[i] += (int64_t)e * (int64_t)col[i];
        }
        for (i = 0; i < nt; ++i) {
            d.e[r+i] = res_mod_torsion(
                    d.e[r+i] + e * col[r+i], g->tord[i]);
        }
    }
    for (i = 0; i < r; ++i) {
        d.e[i] = (int32_t)acc[i];
    }
    free(acc);
}

deg_t res_heft_of_exponents(
        const res_dgrp_t *g,
        const exp_t *exps
        )
{
    len_t j;
    int64_t h = 0;
    const len_t nv          = g->nv;
    const deg_t * const vd  = g->vhdeg;

    for (j = 0; j < nv; ++j) {
        h += (int64_t)exps[j] * (int64_t)vd[j];
    }
    return (deg_t)h;
}
