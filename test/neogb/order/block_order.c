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

/* Standalone tests for block grevlex monomial orders: the general block
 * partition, per-block weights, and the C entry points that carry them.
 *
 * Every expected value below either comes from export_f4, which is the
 * oracle for the two orders that predate blocks, or from Macaulay2 via
 * test/neogb/order/block_reference.m2.  Rerun that script after touching
 * the test systems. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../src/neogb/f4.h"
#include "../../../src/neogb/io.h"
#include "../../../src/neogb/meta_data.h"
#include "../../../src/neogb/res.h"
#include "../../../src/msolve/streams.h"

static int blk_st_fail;
static int blk_st_run;
static int blk_st_verbose;

#define BLK_CHECK(cond, msg)                                            \
    do {                                                                \
        blk_st_run++;                                                   \
        if (!(cond)) {                                                  \
            blk_st_fail++;                                              \
            fprintf(ERRSTREAM, "block_order FAIL %s:%d: %s\n",          \
                    __FILE__, __LINE__, (msg));                         \
        } else if (blk_st_verbose > 1) {                                \
            fprintf(VERBSTREAM, "block_order ok   %s\n", (msg));        \
        }                                                               \
    } while (0)

#define BLK_FC 32003

/* --------------------------------------------------------------------- *
 *  Running a basis
 * --------------------------------------------------------------------- */

typedef struct
{
    int32_t bld;    /* number of basis elements */
    int32_t *blen;  /* terms per element */
    int32_t *bexp;  /* exponents, nvars per term */
    int32_t *bcf;   /* coefficients */
    int64_t nterms; /* total terms, or <= 0 if nothing was computed */
} blk_gb_t;

/* msolve reduces the coefficient array in place and hands back memory it
 * allocated with mallocp, so each run needs a private copy of the input
 * and a matching free below. */
static void blk_run(
        blk_gb_t *g,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *cfs_src,
        const int32_t nvars,
        const int32_t ngens,
        const int32_t elim,
        const mo_block_t *blk
        )
{
    int32_t nr_terms = 0;
    for (int32_t i = 0; i < ngens; ++i) {
        nr_terms += lens[i];
    }

    int32_t *cfs = (int32_t *)malloc((size_t)nr_terms * sizeof(int32_t));
    memcpy(cfs, cfs_src, (size_t)nr_terms * sizeof(int32_t));

    void *bcf = NULL;
    g->bld    = 0;
    g->blen   = NULL;
    g->bexp   = NULL;
    g->nterms = export_f4_blocks(malloc, &g->bld, &g->blen, &g->bexp, &bcf,
            lens, exps, cfs, BLK_FC, 0 /* drl */, elim, blk,
            nvars, ngens, 12 /* ht size */, 1 /* threads */,
            0 /* max pairs */, 0 /* reset ht */, 2 /* la */, 1 /* reduce */,
            0 /* pbm */, 0 /* info */);
    g->bcf = (int32_t *)bcf;

    free(cfs);
}

static void blk_free(
        blk_gb_t *g
        )
{
    free(g->blen);
    free(g->bexp);
    free(g->bcf);
    g->blen = NULL;
    g->bexp = NULL;
    g->bcf  = NULL;
}

/* Term for term equality of two exported bases. */
static int blk_identical(
        const blk_gb_t *a,
        const blk_gb_t *b,
        const int32_t nvars
        )
{
    if (a->nterms != b->nterms || a->bld != b->bld
            || a->nterms <= 0 || a->blen == NULL || b->blen == NULL) {
        return 0;
    }
    if (memcmp(a->blen, b->blen, (size_t)a->bld * sizeof(int32_t)) != 0) {
        return 0;
    }
    const size_t ne = (size_t)a->nterms * (size_t)nvars;
    return memcmp(a->bexp, b->bexp, ne * sizeof(int32_t)) == 0
        && memcmp(a->bcf, b->bcf, (size_t)a->nterms * sizeof(int32_t)) == 0;
}

/* --------------------------------------------------------------------- *
 *  Comparing against a Macaulay2 reference
 * --------------------------------------------------------------------- */

static int32_t blk_cmp_nvars;

/* lexicographic on the exponent tuple, matching how Macaulay2 sorts the
 * lists of lists that block_reference.m2 prints */
static int blk_cmp_rows(
        const void *va,
        const void *vb
        )
{
    const int32_t *a = (const int32_t *)va;
    const int32_t *b = (const int32_t *)vb;

    for (int32_t i = 0; i < blk_cmp_nvars; ++i) {
        if (a[i] != b[i]) {
            return a[i] < b[i] ? -1 : 1;
        }
    }
    return 0;
}

/* Check the leading ideal against a reference table.  Both sides are
 * sorted, so this does not depend on the order the basis comes back in. */
static void blk_check_leads(
        const blk_gb_t *g,
        const int32_t nvars,
        const int32_t *ref,
        const int32_t nref,
        const char *what
        )
{
    if (g->nterms <= 0 || g->blen == NULL) {
        BLK_CHECK(0, what);
        return;
    }
    if (g->bld != nref) {
        BLK_CHECK(0, what);
        fprintf(ERRSTREAM, "  expected %d leading terms, got %d\n",
                nref, g->bld);
        return;
    }

    int32_t *lead = (int32_t *)malloc(
            (size_t)nref * (size_t)nvars * sizeof(int32_t));

    int32_t off = 0;
    for (int32_t i = 0; i < g->bld; ++i) {
        memcpy(lead + (size_t)i * nvars, g->bexp + (size_t)off * nvars,
                (size_t)nvars * sizeof(int32_t));
        off += g->blen[i];
    }

    blk_cmp_nvars = nvars;
    qsort(lead, (size_t)nref, (size_t)nvars * sizeof(int32_t), blk_cmp_rows);

    const int ok = memcmp(lead, ref,
            (size_t)nref * (size_t)nvars * sizeof(int32_t)) == 0;
    BLK_CHECK(ok, what);
    if (!ok) {
        for (int32_t i = 0; i < nref; ++i) {
            fprintf(ERRSTREAM, "  got {");
            for (int32_t j = 0; j < nvars; ++j) {
                fprintf(ERRSTREAM, "%d%s", lead[i * nvars + j],
                        j + 1 < nvars ? "," : "");
            }
            fprintf(ERRSTREAM, "} want {");
            for (int32_t j = 0; j < nvars; ++j) {
                fprintf(ERRSTREAM, "%d%s", ref[i * nvars + j],
                        j + 1 < nvars ? "," : "");
            }
            fprintf(ERRSTREAM, "}\n");
        }
    }

    free(lead);
}

/* Check the full term sequence of a one element basis.  For a principal
 * ideal the basis is the generator itself, so its terms come back in
 * monomial order -- this compares the whole permutation, which is a much
 * stronger statement about the comparison function than any leading
 * ideal can be. */
static void blk_check_terms(
        const blk_gb_t *g,
        const int32_t nvars,
        const int32_t *ref,
        const int32_t nref,
        const char *what
        )
{
    const int ok = g->nterms == nref && g->bld == 1
        && g->blen != NULL && g->blen[0] == nref
        && memcmp(g->bexp, ref,
                (size_t)nref * (size_t)nvars * sizeof(int32_t)) == 0;

    BLK_CHECK(ok, what);
    if (!ok && g->nterms == nref && g->bld == 1) {
        for (int32_t i = 0; i < nref; ++i) {
            fprintf(ERRSTREAM, "  term %d: got {", i);
            for (int32_t j = 0; j < nvars; ++j) {
                fprintf(ERRSTREAM, "%d%s", g->bexp[i * nvars + j],
                        j + 1 < nvars ? "," : "");
            }
            fprintf(ERRSTREAM, "} want {");
            for (int32_t j = 0; j < nvars; ++j) {
                fprintf(ERRSTREAM, "%d%s", ref[i * nvars + j],
                        j + 1 < nvars ? "," : "");
            }
            fprintf(ERRSTREAM, "}\n");
        }
    }
}

/* --------------------------------------------------------------------- *
 *  The test systems
 * --------------------------------------------------------------------- */

/* cyclic 4 over F_32003 */
static const int32_t cyc_lens[4] = {4, 4, 4, 2};
static const int32_t cyc_exps[14 * 4] = {
    1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1,
    1,1,0,0,  0,1,1,0,  0,0,1,1,  1,0,0,1,
    1,1,1,0,  0,1,1,1,  1,0,1,1,  1,1,0,1,
    1,1,1,1,  0,0,0,0
};
static const int32_t cyc_cfs[14] = {
    1, 1, 1, 1,
    1, 1, 1, 1,
    1, 1, 1, 1,
    1, -1
};

/* A single dense cubic; see the comment on blk_check_terms.
 * x0x1x2 + x0^2x3 + x1^3 + x2x3^2 + x0x3^2 + x1x2x3 + x2^2x3 + x0x1x3 */
static const int32_t den_lens[1] = {8};
static const int32_t den_exps[8 * 4] = {
    1,1,1,0,  2,0,0,1,  0,3,0,0,  0,0,1,2,
    1,0,0,2,  0,1,1,1,  0,0,2,1,  1,1,0,1
};
static const int32_t den_cfs[8] = {1, 1, 1, 1, 1, 1, 1, 1};

/* --------------------------------------------------------------------- *
 *  The legacy orders, reproduced through the block entry point
 * --------------------------------------------------------------------- */

/* A one block order is plain grevlex, so it has to agree with export_f4
 * term for term.  This is the compatibility contract that lets the DRL
 * path stay untouched. */
static void blk_test_one_block_is_drl(
        void
        )
{
    const int32_t bsz[1] = {4};
    const mo_block_t blk = {1, bsz, NULL};

    blk_gb_t legacy, blocks;
    blk_run(&legacy, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, NULL);
    blk_run(&blocks, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &blk);

    BLK_CHECK(legacy.nterms > 0, "the grevlex oracle produced a basis");
    BLK_CHECK(blk_identical(&legacy, &blocks, 4),
            "one block reproduces plain grevlex term for term");

    blk_free(&legacy);
    blk_free(&blocks);
}

/* Two blocks of sizes {k, n-k} are exactly the elimination order that
 * elim_block_len describes today, for every split point. */
static void blk_test_two_blocks_are_elimination(
        void
        )
{
    for (int32_t k = 1; k < 4; ++k) {
        const int32_t bsz[2] = {k, 4 - k};
        const mo_block_t blk = {2, bsz, NULL};

        blk_gb_t legacy, blocks;
        blk_run(&legacy, cyc_lens, cyc_exps, cyc_cfs, 4, 4, k, NULL);
        blk_run(&blocks, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &blk);

        char msg[96];
        snprintf(msg, sizeof(msg),
                "two blocks {%d,%d} reproduce the elimination order -e %d",
                k, 4 - k, k);
        BLK_CHECK(blk_identical(&legacy, &blocks, 4), msg);

        blk_free(&legacy);
        blk_free(&blocks);
    }
}

/* An explicit block description overrides elim_block_len rather than
 * combining with it; nothing else may read the old split point. */
static void blk_test_blocks_override_elim(
        void
        )
{
    const int32_t bsz[1] = {4};
    const mo_block_t blk = {1, bsz, NULL};

    blk_gb_t plain, overridden;
    blk_run(&plain, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, NULL);
    blk_run(&overridden, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 2, &blk);

    BLK_CHECK(blk_identical(&plain, &overridden, 4),
            "an explicit block order overrides elim_block_len");

    blk_free(&plain);
    blk_free(&overridden);
}

/* All weights equal to one is the unweighted order, and must not merely
 * agree on the leading ideal but on every exported term. */
static void blk_test_unit_weights_are_unweighted(
        void
        )
{
    const int32_t bsz[2] = {2, 2};
    const int32_t ones[4] = {1, 1, 1, 1};
    const mo_block_t plain = {2, bsz, NULL};
    const mo_block_t unit  = {2, bsz, ones};

    blk_gb_t a, b;
    blk_run(&a, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &plain);
    blk_run(&b, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &unit);

    BLK_CHECK(blk_identical(&a, &b, 4),
            "unit weights reproduce the unweighted block order");

    blk_free(&a);
    blk_free(&b);
}

/* --------------------------------------------------------------------- *
 *  Term order, against Macaulay2
 * --------------------------------------------------------------------- */

static void blk_test_term_order(
        void
        )
{
    /* MonomialOrder => {GRevLex => 4} */
    static const int32_t r_1[8 * 4] = {
        0,3,0,0,  1,1,1,0,  2,0,0,1,  1,1,0,1,
        0,1,1,1,  0,0,2,1,  1,0,0,2,  0,0,1,2
    };
    /* {GRevLex => 1, GRevLex => 3}, and equally {1,1,2} and {1,1,1,1} */
    static const int32_t r_13[8 * 4] = {
        2,0,0,1,  1,1,1,0,  1,1,0,1,  1,0,0,2,
        0,3,0,0,  0,1,1,1,  0,0,2,1,  0,0,1,2
    };
    /* {GRevLex => 2, GRevLex => 2}, and equally {2,1,1} */
    static const int32_t r_22[8 * 4] = {
        0,3,0,0,  2,0,0,1,  1,1,1,0,  1,1,0,1,
        1,0,0,2,  0,1,1,1,  0,0,2,1,  0,0,1,2
    };
    /* {GRevLex => 3, GRevLex => 1} */
    static const int32_t r_31[8 * 4] = {
        0,3,0,0,  1,1,1,0,  2,0,0,1,  1,1,0,1,
        0,1,1,1,  0,0,2,1,  1,0,0,2,  0,0,1,2
    };
    /* {GRevLex => {1,2,3,4}} */
    static const int32_t r_w[8 * 4] = {
        0,0,1,2,  0,0,2,1,  0,1,1,1,  1,0,0,2,
        1,1,0,1,  0,3,0,0,  1,1,1,0,  2,0,0,1
    };

    static const int32_t b1[1]  = {4};
    static const int32_t b13[2] = {1, 3};
    static const int32_t b22[2] = {2, 2};
    static const int32_t b31[2] = {3, 1};
    static const int32_t b112[3] = {1, 1, 2};
    static const int32_t b211[3] = {2, 1, 1};
    static const int32_t b1111[4] = {1, 1, 1, 1};
    static const int32_t w1234[4] = {1, 2, 3, 4};
    static const int32_t w4321[4] = {4, 3, 2, 1};

    const struct {
        mo_block_t blk;
        const int32_t *ref;
        const char *what;
    } cases[] = {
        {{1, b1,    NULL},  r_1,  "term order of one block"},
        {{2, b13,   NULL},  r_13, "term order of two blocks {1,3}"},
        {{2, b22,   NULL},  r_22, "term order of two blocks {2,2}"},
        {{2, b31,   NULL},  r_31, "term order of two blocks {3,1}"},
        {{3, b112,  NULL},  r_13, "term order of three blocks {1,1,2}"},
        {{3, b211,  NULL},  r_22, "term order of three blocks {2,1,1}"},
        {{4, b1111, NULL},  r_13, "term order of four blocks {1,1,1,1}"},
        {{1, b1,    w1234}, r_w,  "term order of one block weighted 1,2,3,4"},
        {{2, b22,   w4321}, r_22, "term order of two blocks {2,2} weighted 4,3,2,1"}
    };
    const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

    for (int i = 0; i < ncases; ++i) {
        blk_gb_t g;
        blk_run(&g, den_lens, den_exps, den_cfs, 4, 1, 0, &cases[i].blk);
        blk_check_terms(&g, 4, cases[i].ref, 8, cases[i].what);
        blk_free(&g);
    }
}

/* --------------------------------------------------------------------- *
 *  Leading ideals, against Macaulay2
 * --------------------------------------------------------------------- */

static void blk_test_leading_ideals(
        void
        )
{
    /* {GRevLex => 4}, and equally {GRevLex => 1, GRevLex => 3} */
    static const int32_t l_1[7 * 4] = {
        0,0,2,4,  0,0,3,2,  0,1,0,4,  0,1,1,2,
        0,1,2,0,  0,2,0,0,  1,0,0,0
    };
    /* {GRevLex => 2, GRevLex => 2}, and equally {1,1,2} and {1,1,1,1} */
    static const int32_t l_22[6 * 4] = {
        0,0,2,6,  0,0,3,2,  0,1,0,4,  0,1,1,0,
        0,2,0,0,  1,0,0,0
    };
    /* {GRevLex => {1,2,3,4}} */
    static const int32_t l_w[7 * 4] = {
        0,0,0,1,  0,0,2,0,  0,2,1,0,  2,1,1,0,
        2,3,0,0,  4,0,1,0,  4,2,0,0
    };
    /* {GRevLex => {1,2}, GRevLex => {3,4}} */
    static const int32_t l_22w[6 * 4] = {
        0,0,2,3,  0,0,6,2,  0,1,0,0,  1,0,0,1,
        1,0,4,0,  2,0,0,0
    };

    static const int32_t b1[1]  = {4};
    static const int32_t b13[2] = {1, 3};
    static const int32_t b22[2] = {2, 2};
    static const int32_t b112[3] = {1, 1, 2};
    static const int32_t b1111[4] = {1, 1, 1, 1};
    static const int32_t w1234[4] = {1, 2, 3, 4};

    const struct {
        mo_block_t blk;
        const int32_t *ref;
        int32_t nref;
        const char *what;
    } cases[] = {
        {{1, b1,    NULL},  l_1,   7, "cyclic4 leading ideal, one block"},
        {{2, b13,   NULL},  l_1,   7, "cyclic4 leading ideal, blocks {1,3}"},
        {{2, b22,   NULL},  l_22,  6, "cyclic4 leading ideal, blocks {2,2}"},
        {{3, b112,  NULL},  l_22,  6, "cyclic4 leading ideal, blocks {1,1,2}"},
        {{4, b1111, NULL},  l_22,  6, "cyclic4 leading ideal, blocks {1,1,1,1}"},
        {{1, b1,    w1234}, l_w,   7, "cyclic4 leading ideal, weighted 1,2,3,4"},
        {{2, b22,   w1234}, l_22w, 6, "cyclic4 leading ideal, weighted blocks {2,2}"}
    };
    const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

    for (int i = 0; i < ncases; ++i) {
        blk_gb_t g;
        blk_run(&g, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &cases[i].blk);
        blk_check_leads(&g, 4, cases[i].ref, cases[i].nref, cases[i].what);
        blk_free(&g);
    }
}

/* A block order refines grevlex within each block, so it computes a
 * Groebner basis of the same ideal.  For a zero dimensional ideal the
 * quotient dimension is an order independent invariant, so counting
 * standard monomials under each leading ideal has to give the same
 * answer every time.  That catches a comparison function consistent
 * enough to terminate but not actually a monomial order -- which no
 * comparison against a fixed reference table can, since a wrong order
 * would simply be wrong in the same way on both sides.
 *
 * cyclic4 is one dimensional, so it cannot serve here; this uses a
 * zero dimensional system of degree 16 instead.
 *
 *   x_i^2 + (sum of the other three variables) - 1,  i = 0..3 */
static const int32_t zd_lens[4] = {5, 5, 5, 5};
static const int32_t zd_exps[20 * 4] = {
    2,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1,  0,0,0,0,
    0,2,0,0,  1,0,0,0,  0,0,1,0,  0,0,0,1,  0,0,0,0,
    0,0,2,0,  1,0,0,0,  0,1,0,0,  0,0,0,1,  0,0,0,0,
    0,0,0,2,  1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,0
};
static const int32_t zd_cfs[20] = {
    1, 1, 1, 1, -1,
    1, 1, 1, 1, -1,
    1, 1, 1, 1, -1,
    1, 1, 1, 1, -1
};

static void blk_test_same_ideal(
        void
        )
{
    static const int32_t b1[1]    = {4};
    static const int32_t b13[2]   = {1, 3};
    static const int32_t b22[2]   = {2, 2};
    static const int32_t b31[2]   = {3, 1};
    static const int32_t b112[3]  = {1, 1, 2};
    static const int32_t b1111[4] = {1, 1, 1, 1};
    static const int32_t w1234[4] = {1, 2, 3, 4};

    const mo_block_t blks[] = {
        {1, b1,    NULL},
        {2, b13,   NULL},
        {2, b22,   NULL},
        {2, b31,   NULL},
        {3, b112,  NULL},
        {4, b1111, NULL},
        {1, b1,    w1234},
        {2, b22,   w1234}
    };
    const int nblks = (int)(sizeof(blks) / sizeof(blks[0]));

    /* the system has 16 solutions counted with multiplicity, so R/in(I)
     * is 16 dimensional whatever the order (Macaulay2: degree I) */
    const int32_t expected_dim = 16;
    /* the staircase fits well inside this box; escaping it is checked
     * for below rather than assumed */
    const int32_t bound = 20;

    for (int i = 0; i < nblks; ++i) {
        blk_gb_t g;
        blk_run(&g, zd_lens, zd_exps, zd_cfs, 4, 4, 0, &blks[i]);

        /* count the monomials divisible by no leading term */
        int32_t dim = 0, on_boundary = 0;
        if (g.nterms > 0 && g.blen != NULL) {
            for (int32_t e0 = 0; e0 <= bound; ++e0) {
            for (int32_t e1 = 0; e1 <= bound; ++e1) {
            for (int32_t e2 = 0; e2 <= bound; ++e2) {
            for (int32_t e3 = 0; e3 <= bound; ++e3) {
                const int32_t m[4] = {e0, e1, e2, e3};
                int standard = 1;
                int32_t off = 0;
                for (int32_t k = 0; k < g.bld && standard; ++k) {
                    const int32_t *l = g.bexp + (size_t)off * 4;
                    int divides = 1;
                    for (int32_t v = 0; v < 4; ++v) {
                        if (l[v] > m[v]) {
                            divides = 0;
                            break;
                        }
                    }
                    if (divides) {
                        standard = 0;
                    }
                    off += g.blen[k];
                }
                dim += standard;
                if (standard && (e0 == bound || e1 == bound
                            || e2 == bound || e3 == bound)) {
                    on_boundary = 1;
                }
            }}}}
        }

        char msg[96];
        snprintf(msg, sizeof(msg),
                "block order %d gives a %d dimensional quotient",
                i, expected_dim);
        BLK_CHECK(dim == expected_dim && !on_boundary, msg);
        if (dim != expected_dim || on_boundary) {
            fprintf(ERRSTREAM, "  got %d standard monomials%s\n", dim,
                    on_boundary ? " and the staircase left the box" : "");
        }

        blk_free(&g);
    }
}

/* --------------------------------------------------------------------- *
 *  Block orders on a free module
 * --------------------------------------------------------------------- */

/* A rank one free module is the ring itself, so a module Groebner basis
 * over R^1 has to agree term for term with the ideal Groebner basis of
 * the same generators -- in whatever block order both are computed in.
 * That makes export_f4_blocks an exact oracle for the whole module path
 * under blocks: the import with components, cmp_blocks reading a module
 * exponent vector, get_lcm, and the export. */
static void blk_test_module_rank_one(
        const int32_t module_order,
        const char *order_name
        )
{
    static const int32_t b1[1]    = {4};
    static const int32_t b13[2]   = {1, 3};
    static const int32_t b22[2]   = {2, 2};
    static const int32_t b112[3]  = {1, 1, 2};
    static const int32_t b1111[4] = {1, 1, 1, 1};
    static const int32_t w1234[4] = {1, 2, 3, 4};

    const mo_block_t blks[] = {
        {1, b1,    NULL},
        {2, b13,   NULL},
        {2, b22,   NULL},
        {3, b112,  NULL},
        {4, b1111, NULL},
        {1, b1,    w1234},
        {2, b22,   w1234}
    };
    const int nblks = (int)(sizeof(blks) / sizeof(blks[0]));

    /* every term of every generator sits in the single component */
    int32_t comps[14];
    for (int i = 0; i < 14; ++i) {
        comps[i] = 1;
    }

    const res_strat_t strat = res_strat_of_order(module_order);

    for (int i = 0; i < nblks; ++i) {
        blk_gb_t ideal;
        blk_run(&ideal, cyc_lens, cyc_exps, cyc_cfs, 4, 4, 0, &blks[i]);

        int32_t cfs[14];
        memcpy(cfs, cyc_cfs, sizeof(cyc_cfs));

        int32_t mbld = 0, *mblen = NULL, *mbexp = NULL, *mbcomp = NULL;
        void *mbcf = NULL;
        const int64_t mterms = export_module_f4_blocks(malloc,
                &mbld, &mblen, &mbexp, &mbcomp, &mbcf,
                cyc_lens, cyc_exps, comps, cfs, NULL /* row degrees */,
                BLK_FC, 0 /* drl */, &blks[i], &strat, NULL /* grading */,
                4 /* nvars */, 1 /* nrows */, 4 /* ngens */,
                12 /* ht size */, 1 /* threads */, 0 /* max pairs */,
                2 /* la */, 1 /* reduce */, 0 /* info */);

        blk_gb_t module;
        module.bld    = mbld;
        module.blen   = mblen;
        module.bexp   = mbexp;
        module.bcf    = (int32_t *)mbcf;
        module.nterms = mterms;

        char msg[128];
        snprintf(msg, sizeof(msg),
                "rank one %s module basis equals the ideal basis, block order %d",
                order_name, i);
        BLK_CHECK(blk_identical(&ideal, &module, 4), msg);

        /* a rank one module has only the one component */
        int ok_comp = mterms > 0;
        for (int64_t t = 0; t < mterms; ++t) {
            if (mbcomp[t] != 1) {
                ok_comp = 0;
                break;
            }
        }
        BLK_CHECK(ok_comp, "every term of a rank one basis is in component 1");

        free(mbcomp);
        blk_free(&module);
        blk_free(&ideal);
    }
}

/* The block description has to actually reach the module order, not just
 * be accepted and dropped.  Under an elimination block the leading terms
 * of a rank two module basis differ from the plain grevlex ones, so a
 * block description that went nowhere shows up as equality here. */
static void blk_test_module_blocks_take_effect(
        void
        )
{
    /* over R^2:  (x0^2 + x2*x3) e1 + x1 e2,  (x1^2) e1 + (x0*x3) e2,
     *            (x3^2) e1 + (x0*x1) e2 */
    const int32_t lens[3] = {3, 2, 2};
    const int32_t exps[7 * 4] = {
        2,0,0,0,  0,0,1,1,  0,1,0,0,
        0,2,0,0,  1,0,0,1,
        0,0,0,2,  1,1,0,0
    };
    const int32_t comps[7] = {1, 1, 2, 1, 2, 1, 2};
    const int32_t cfs_src[7] = {1, 1, 1, 1, 1, 1, 1};

    static const int32_t b1[1]  = {4};
    static const int32_t b22[2] = {2, 2};
    const mo_block_t one   = {1, b1,  NULL};
    const mo_block_t two   = {2, b22, NULL};

    const res_strat_t strat = res_strat_of_order(RES_MORD_TOP);

    int32_t bld[2] = {0, 0}, *blen[2] = {NULL, NULL};
    int32_t *bexp[2] = {NULL, NULL}, *bcomp[2] = {NULL, NULL};
    void *bcf[2] = {NULL, NULL};
    int64_t nterms[2] = {0, 0};

    const mo_block_t *which[2] = {&one, &two};
    for (int i = 0; i < 2; ++i) {
        int32_t cfs[7];
        memcpy(cfs, cfs_src, sizeof(cfs_src));
        nterms[i] = export_module_f4_blocks(malloc,
                &bld[i], &blen[i], &bexp[i], &bcomp[i], &bcf[i],
                lens, exps, comps, cfs, NULL,
                BLK_FC, 0, which[i], &strat, NULL,
                4, 2 /* nrows */, 3, 12, 1, 0, 2, 1, 0);
    }

    BLK_CHECK(nterms[0] > 0 && nterms[1] > 0,
            "both module bases were computed");

    const int same = nterms[0] == nterms[1] && bld[0] == bld[1]
        && bld[0] > 0
        && memcmp(blen[0], blen[1], (size_t)bld[0] * sizeof(int32_t)) == 0
        && memcmp(bexp[0], bexp[1],
                (size_t)nterms[0] * 4 * sizeof(int32_t)) == 0;
    BLK_CHECK(!same,
            "an elimination block changes the module basis it is passed to");

    for (int i = 0; i < 2; ++i) {
        free(blen[i]);
        free(bexp[i]);
        free(bcomp[i]);
        free(bcf[i]);
    }
}

/* ht->vwt is one array with one meaning, so order weights and grading
 * weights are the same statement.  Passing both has to be an error when
 * they disagree and fine when they agree. */
static void blk_test_module_weight_conflict(
        void
        )
{
    const int32_t lens[2] = {2, 2};
    const int32_t exps[4 * 4] = {
        2,0,0,0,  0,0,1,1,
        0,2,0,0,  1,0,0,1
    };
    const int32_t comps[4] = {1, 1, 1, 1};
    const int32_t cfs_src[4] = {1, 1, 1, 1};

    static const int32_t b1[1]     = {4};
    static const int32_t w1234[4]  = {1, 2, 3, 4};
    static const int32_t w1235[4]  = {1, 2, 3, 5};
    const mo_block_t blk = {1, b1, w1234};

    /* the same weights spelled as a grading, and a different one */
    const res_grading_t agree    = {1, 0, NULL, w1234, NULL};
    const res_grading_t disagree = {1, 0, NULL, w1235, NULL};

    const struct {
        const res_grading_t *grading;
        int want_ok;
        const char *what;
    } cases[] = {
        {&agree,    1, "matching order and grading weights are accepted"},
        {&disagree, 0, "conflicting order and grading weights are rejected"}
    };

    for (int i = 0; i < 2; ++i) {
        int32_t cfs[4];
        memcpy(cfs, cfs_src, sizeof(cfs_src));

        int32_t bld = 0, *blen = NULL, *bexp = NULL, *bcomp = NULL;
        void *bcf = NULL;
        const int64_t nterms = export_module_f4_blocks(malloc,
                &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, comps, cfs, NULL,
                BLK_FC, 0, &blk, NULL /* default strategy */,
                cases[i].grading, 4, 1, 2, 12, 1, 0, 2, 1, 0);

        BLK_CHECK((nterms > 0) == cases[i].want_ok, cases[i].what);

        free(blen);
        free(bexp);
        free(bcomp);
        free(bcf);
    }
}

/* --------------------------------------------------------------------- *
 *  Validation
 * --------------------------------------------------------------------- */

/* Bad block descriptions have to be reported through the return value.
 * These go straight at set_monomial_block_order rather than through
 * export_f4_blocks, which follows msolve's older convention of exiting
 * on corrupt input. */
static void blk_test_rejects_bad_input(
        void
        )
{
    static const int32_t ok4[2]    = {2, 2};
    static const int32_t short3[2] = {2, 1};
    static const int32_t empty[2]  = {4, 0};
    static const int32_t neg[2]    = {5, -1};
    static const int32_t w_ok[4]   = {1, 2, 3, 4};
    static const int32_t w_zero[4] = {1, 0, 3, 4};
    static const int32_t w_neg[4]  = {1, 2, -3, 4};
    static const int32_t w_huge[4] = {1, 2, 3, 1 << 20};

    const struct {
        mo_block_t blk;
        int want_ok;
        const char *what;
    } cases[] = {
        {{1, ok4,    NULL},   0, "a one block order of the wrong size is rejected"},
        {{2, ok4,    NULL},   1, "a well formed two block order is accepted"},
        {{2, ok4,    w_ok},   1, "positive weights are accepted"},
        {{0, ok4,    NULL},   0, "a zero block order is rejected"},
        {{-1, ok4,   NULL},   0, "a negative block count is rejected"},
        {{2, NULL,   NULL},   0, "missing block sizes are rejected"},
        {{2, short3, NULL},   0, "block sizes that do not sum to nvars are rejected"},
        {{2, empty,  NULL},   0, "an empty block is rejected"},
        {{2, neg,    NULL},   0, "a negative block size is rejected"},
        {{2, ok4,    w_zero}, 0, "a zero weight is rejected"},
        {{2, ok4,    w_neg},  0, "a negative weight is rejected"},
        {{2, ok4,    w_huge}, 0, "a weight that would overflow exp_t is rejected"}
    };
    const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

    /* the rejections print a diagnostic each; say so, so that a passing
     * run does not look like a failing one */
    if (blk_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "block_order: the following diagnostics are expected\n");
    }

    for (int i = 0; i < ncases; ++i) {
        md_t *st = allocate_meta_data();
        st->nvars = 4;

        const int32_t err = set_monomial_block_order(st, &cases[i].blk);
        BLK_CHECK((err == 0) == cases[i].want_ok, cases[i].what);

        free_block_order(st);
        free(st);
    }

    /* a NULL description keeps whatever the caller already had */
    md_t *st = allocate_meta_data();
    st->nvars = 4;
    st->nbl   = 7;
    BLK_CHECK(set_monomial_block_order(st, NULL) == 0
            && st->nbl == 7,
            "a NULL block description leaves the order alone");
    free_block_order(st);
    free(st);
}

/* --------------------------------------------------------------------- *
 *  Entry point
 * --------------------------------------------------------------------- */

int main(void)
{
    const int verbose = 1;

    blk_st_fail    = 0;
    blk_st_run     = 0;
    blk_st_verbose = verbose;

    blk_test_one_block_is_drl();
    blk_test_two_blocks_are_elimination();
    blk_test_blocks_override_elim();
    blk_test_unit_weights_are_unweighted();
    blk_test_term_order();
    blk_test_leading_ideals();
    blk_test_same_ideal();
    blk_test_module_rank_one(RES_MORD_POT, "POT");
    blk_test_module_rank_one(RES_MORD_TOP, "TOP");
    blk_test_module_blocks_take_effect();
    blk_test_module_weight_conflict();
    blk_test_rejects_bad_input();

    if (verbose > 0) {
        fprintf(VERBSTREAM, "block_order: %d checks, %d failures\n",
                blk_st_run, blk_st_fail);
    }
    return blk_st_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
