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

/* Standalone tests for the public grading and module Groebner-basis APIs
 * used by the resolution engine. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../src/neogb/f4.h"
#include "../../../src/neogb/res.h"
#include "../../../src/msolve/streams.h"

/* The entry points take a res_strat_t now rather than a bare res_mord_t.
 * Most of these tests care only about the base order, so this wraps one
 * up with the default position and lift.  The small rotating pool lets
 * two of them appear in one expression. */
static const res_strat_t *res_strat_p(
        const int32_t mord
        )
{
    static res_strat_t pool[4];
    static int next = 0;

    res_strat_t *s = pool + next;
    next = (next + 1) & 3;
    *s = res_strat_of_order(mord);

    return s;
}


static int res_st_fail;
static int res_st_run;
static int res_st_verbose;

#define RES_CHECK(cond, msg)                                            \
    do {                                                                \
        res_st_run++;                                                   \
        if (!(cond)) {                                                  \
            res_st_fail++;                                              \
            fprintf(ERRSTREAM, "res_selftest FAIL %s:%d: %s\n",         \
                    __FILE__, __LINE__, (msg));                         \
        } else if (res_st_verbose > 1) {                                \
            fprintf(VERBSTREAM, "res_selftest ok   %s\n", (msg));       \
        }                                                               \
    } while (0)

/* --------------------------------------------------------------------- *
 *  Grading groups
 * --------------------------------------------------------------------- */

static void res_test_standard_grading(
        void
        )
{
    exp_t e[3];
    res_dgrp_t *g = res_dgrp_new_standard(3);

    RES_CHECK(g != NULL, "standard grading is constructed");
    if (g == NULL) {
        return;
    }
    RES_CHECK(g->simple == 1, "standard grading takes the r==1 fast path");
    RES_CHECK(g->len == 1, "standard grading stores one slot per degree");

    /* x1^2 * x2 */
    e[0] = 2; e[1] = 1; e[2] = 0;
    RES_CHECK(res_heft_of_exponents(g, e) == 3,
            "standard heft degree of x1^2*x2 is 3");

    res_dpool_t *p = res_dpool_new(g, 2);
    res_deg_t a = res_dpool_push(p, NULL);
    res_deg_t b = res_dpool_push(p, NULL);
    res_deg_t c = res_dpool_push(p, NULL); /* forces a growth */

    /* A growth invalidates degree views.  Reacquire them by index before
     * using them; this also keeps the test honest about the pool API. */
    a = res_dpool_at(p, 0);
    b = res_dpool_at(p, 1);
    c = res_dpool_at(p, 2);

    res_deg_of_exponents(g, a, e);
    RES_CHECK(a.e[0] == 3, "standard multidegree of x1^2*x2 is 3");

    e[0] = 1; e[1] = 0; e[2] = 0;
    res_deg_of_exponents(g, b, e);
    g->add(g, c, a, b);
    RES_CHECK(c.e[0] == 4, "degree addition in the standard grading");
    g->sub(g, c, a, b);
    RES_CHECK(c.e[0] == 2, "degree subtraction in the standard grading");
    RES_CHECK(g->cmp(g, a, b) > 0, "degree 3 compares above degree 1");
    RES_CHECK(g->cmp(g, a, a) == 0, "degree comparison is reflexive");
    RES_CHECK(g->heft_of(g, a) == 3, "heft of the standard degree 3");

    res_dpool_free(&p);
    RES_CHECK(p == NULL, "degree pool is nulled on free");
    res_dgrp_free(&g);
    RES_CHECK(g == NULL, "grading group is nulled on free");
}

static void res_test_multigrading(
        void
        )
{
    /* P^1 x P^1 style bigrading on 4 variables */
    const int32_t degs[8] = {
        1, 0,   /* deg(x1) */
        1, 0,   /* deg(x2) */
        0, 1,   /* deg(x3) */
        0, 1    /* deg(x4) */
    };
    const int32_t heft[2] = {1, 1};
    exp_t e[4];

    res_dgrp_t *g = res_dgrp_new(2, 0, NULL, 4, degs, heft);
    RES_CHECK(g != NULL, "bigrading is constructed");
    if (g == NULL) {
        return;
    }
    RES_CHECK(g->simple == 0, "bigrading uses the generic path");
    RES_CHECK(g->len == 2, "bigrading stores two slots per degree");

    res_dpool_t *p = res_dpool_new(g, 4);
    res_deg_t a = res_dpool_push(p, NULL);

    /* x1 * x3^2 has bidegree (1,2) and heft degree 3 */
    e[0] = 1; e[1] = 0; e[2] = 2; e[3] = 0;
    res_deg_of_exponents(g, a, e);
    RES_CHECK(a.e[0] == 1 && a.e[1] == 2, "bidegree of x1*x3^2 is (1,2)");
    RES_CHECK(res_heft_of_exponents(g, e) == 3,
            "heft degree of x1*x3^2 is 3");
    RES_CHECK(g->heft_of(g, a) == 3, "heft_of agrees with heft_of_exponents");

    /* two monomials of the same heft degree but different bidegree must
     * not collide: this is exactly what rank bucketing relies on */
    res_deg_t b = res_dpool_push(p, NULL);
    e[0] = 2; e[1] = 1; e[2] = 0; e[3] = 0;
    res_deg_of_exponents(g, b, e);
    RES_CHECK(g->heft_of(g, b) == 3, "x1^2*x2 also has heft degree 3");
    RES_CHECK(g->cmp(g, a, b) != 0,
            "distinct bidegrees of equal heft degree stay distinct");
    RES_CHECK(g->hash(g, a) != g->hash(g, b),
            "distinct bidegrees hash apart");

    res_dpool_free(&p);
    res_dgrp_free(&g);
}

static void res_test_torsion_grading(
        void
        )
{
    /* A = Z (+) Z/3, as for the class group of a fake weighted space */
    const int32_t degs[4] = {
        1, 1,   /* deg(x1) = (1, 1 mod 3) */
        1, 2    /* deg(x2) = (1, 2 mod 3) */
    };
    const int32_t tord[1] = {3};
    const int32_t heft[1] = {1};
    exp_t e[2];

    res_dgrp_t *g = res_dgrp_new(1, 1, tord, 2, degs, heft);
    RES_CHECK(g != NULL, "grading with torsion is constructed");
    if (g == NULL) {
        return;
    }
    RES_CHECK(g->simple == 0, "torsion forces the generic path");
    RES_CHECK(g->len == 2, "Z (+) Z/3 stores two slots per degree");

    res_dpool_t *p = res_dpool_new(g, 4);
    res_deg_t a = res_dpool_push(p, NULL);

    /* x1*x2 has degree (2, (1+2) mod 3) = (2, 0) */
    e[0] = 1; e[1] = 1;
    res_deg_of_exponents(g, a, e);
    RES_CHECK(a.e[0] == 2, "free part of deg(x1*x2) is 2");
    RES_CHECK(a.e[1] == 0, "torsion part of deg(x1*x2) is 0 mod 3");

    /* x1^2*x2^2 has degree (4, (2+4) mod 3) = (4, 0) */
    e[0] = 2; e[1] = 2;
    res_deg_of_exponents(g, a, e);
    RES_CHECK(a.e[0] == 4 && a.e[1] == 0,
            "torsion reduces correctly on higher powers");

    /* torsion must wrap under addition, never grow */
    res_deg_t b = res_dpool_push(p, NULL);
    res_deg_t c = res_dpool_push(p, NULL);
    e[0] = 1; e[1] = 0;
    res_deg_of_exponents(g, b, e);          /* (1, 1) */
    RES_CHECK(b.e[1] == 1, "deg(x1) has torsion part 1");
    g->add(g, c, b, b);
    RES_CHECK(c.e[0] == 2 && c.e[1] == 2, "deg(x1^2) = (2, 2)");
    g->add(g, c, c, b);
    RES_CHECK(c.e[0] == 3 && c.e[1] == 0, "deg(x1^3) = (3, 0), wrapped");
    g->sub(g, c, b, c);
    RES_CHECK(c.e[0] == -2 && c.e[1] == 1,
            "subtraction keeps torsion in [0,3)");
    RES_CHECK(g->heft_of(g, b) == 1, "torsion carries no heft");

    res_dpool_free(&p);
    res_dgrp_free(&g);
}

static void res_test_grading_rejects_bad_heft(
        void
        )
{
    /* deg(x2) = 0, so no heft vector can be positive on it */
    const int32_t degs[2] = {1, 0};
    const int32_t heft[1] = {1};

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "res_selftest: one grading error message is expected next\n");
    }
    res_dgrp_t *g = res_dgrp_new(1, 0, NULL, 2, degs, heft);
    RES_CHECK(g == NULL, "a non positive heft is rejected");
    if (g != NULL) {
        res_dgrp_free(&g);
    }
}

/* --------------------------------------------------------------------- *
 *  Module Groebner bases, end to end through the C entry point
 * --------------------------------------------------------------------- */

/* A rank one free module is the ring itself, so a module Groebner basis
 * over R^1 has to agree term for term with the ideal Groebner basis of
 * the same generators.  That makes export_f4 an exact oracle for
 * export_module_f4 and exercises the whole module path -- import with
 * components, the module order, get_lcm, the export -- against code that
 * is already known to be right. */
static void res_test_module_gb_rank_one(
        const int32_t module_order,
        const char *what
        )
{
    /* x1^2 + x2*x3, x2^2 + x1*x3, x3^2 + x1*x2 over F_32003 */
    const int32_t lens[3] = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,
        0,2,0,  1,0,1,
        0,0,2,  1,1,0
    };
    const int32_t cfs_src[6] = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6]   = {1, 1, 1, 1, 1, 1};
    const uint32_t fc = 32003;

    /* msolve reduces the coefficient array in place, so each run needs a
     * private copy */
    int32_t cfs_a[6], cfs_b[6];
    memcpy(cfs_a, cfs_src, sizeof(cfs_src));
    memcpy(cfs_b, cfs_src, sizeof(cfs_src));

    int32_t ibld = 0, *iblen = NULL, *ibexp = NULL;
    void *ibcf = NULL;
    const int64_t interms = export_f4(malloc, &ibld, &iblen, &ibexp, &ibcf,
            lens, exps, cfs_a, fc, 0 /* drl */, 0 /* no elim block */,
            3 /* nvars */, 3 /* ngens */, 12 /* ht size */, 1 /* threads */,
            0 /* max pairs */, 0 /* reset ht */, 2 /* la */, 1 /* reduce */,
            0 /* pbm */, 0 /* info */);

    int32_t mbld = 0, *mblen = NULL, *mbexp = NULL, *mbcomp = NULL;
    void *mbcf = NULL;
    const int64_t mterms = export_module_f4(malloc,
            &mbld, &mblen, &mbexp, &mbcomp, &mbcf,
            lens, exps, comps, cfs_b, NULL /* row degrees */,
            fc, 0 /* drl */, res_strat_p(module_order), NULL, NULL /* stop */,
            3 /* nvars */, 1 /* nrows */, 3 /* ngens */,
            12 /* ht size */, 1 /* threads */, 0 /* max pairs */,
            2 /* la */, 1 /* reduce */, 0 /* info */);

    RES_CHECK(interms > 0, "the ideal oracle produced a basis");
    RES_CHECK(mterms == interms, what);
    RES_CHECK(mbld == ibld, "rank one module and ideal bases agree in size");

    if (mterms == interms && mbld == ibld
            && iblen != NULL && mblen != NULL) {
        int ok_len = 1, ok_exp = 1, ok_cf = 1, ok_comp = 1;
        int64_t t;
        for (t = 0; t < mbld; ++t) {
            if (mblen[t] != iblen[t]) {
                ok_len = 0;
            }
        }
        for (t = 0; t < interms * 3; ++t) {
            if (mbexp[t] != ibexp[t]) {
                ok_exp = 0;
            }
        }
        for (t = 0; t < interms; ++t) {
            if (((int32_t *)mbcf)[t] != ((int32_t *)ibcf)[t]) {
                ok_cf = 0;
            }
            if (mbcomp[t] != 1) {
                ok_comp = 0;
            }
        }
        RES_CHECK(ok_len, "rank one basis element lengths agree");
        RES_CHECK(ok_exp, "rank one basis exponents agree");
        RES_CHECK(ok_cf, "rank one basis coefficients agree");
        RES_CHECK(ok_comp, "every term of a rank one basis is in component 1");
    }

    free_f4_julia_result_data(free, &iblen, &ibexp, &ibcf,
            (int64_t)ibld, (int32_t)fc);
    free_module_f4_result_data(free, &mblen, &mbexp, &mbcomp, &mbcf);
}

/* Generators living in disjoint components generate no S-pairs across
 * them, so the basis is the union of the component-wise ideal bases.
 * Here <x1^2, x1*x2> sits in component 1 and <x3> in component 2. */
static void res_test_module_gb_split_components(
        void
        )
{
    const int32_t lens[3] = {1, 1, 1};
    const int32_t exps[9] = {
        2,0,0,   /* x1^2 */
        1,1,0,   /* x1*x2 */
        0,0,1    /* x3    */
    };
    int32_t cfs[3]  = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 2};
    const uint32_t fc = 32003;

    int32_t bld = 0, *blen = NULL, *bexp = NULL, *bcomp = NULL;
    void *bcf = NULL;
    const int64_t nterms = export_module_f4(malloc,
            &bld, &blen, &bexp, &bcomp, &bcf,
            lens, exps, comps, cfs, NULL,
            fc, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */,
            3 /* nvars */, 2 /* nrows */, 3 /* ngens */,
            12, 1, 0, 2, 1, 0);

    RES_CHECK(nterms == 3, "the split module basis has three terms");
    RES_CHECK(bld == 3, "the split module basis has three elements");
    if (bld == 3 && nterms == 3) {
        /* x1^2 and x1*x2 are already a Groebner basis of their ideal and
         * x3 of its own, so nothing is added or removed */
        int n_c1 = 0, n_c2 = 0;
        int64_t t;
        for (t = 0; t < nterms; ++t) {
            if (bcomp[t] == 1) n_c1++;
            if (bcomp[t] == 2) n_c2++;
        }
        RES_CHECK(n_c1 == 2, "two basis terms stayed in component 1");
        RES_CHECK(n_c2 == 1, "one basis term stayed in component 2");
    }
    free_module_f4_result_data(free, &blen, &bexp, &bcomp, &bcf);
}

/* A genuine rank two computation whose answer is known independently:
 * the columns (x1, x2) and (x2, x1) of a 2x2 matrix over k[x1,x2,x3].
 * Their S-pair is x2*(x1,x2) - x1*(x2,x1) = (0, x2^2 - x1^2) under any
 * order, so the basis must acquire an element supported in one component
 * only. */
static void res_test_module_gb_rank_two(
        void
        )
{
    const int32_t lens[2] = {2, 2};
    const int32_t exps[12] = {
        1,0,0,  0,1,0,   /* first column:  x1*e1 + x2*e2 */
        0,1,0,  1,0,0    /* second column: x2*e1 + x1*e2 */
    };
    int32_t cfs[4] = {1, 1, 1, 1};
    const int32_t comps[4] = {1, 2, 1, 2};
    const uint32_t fc = 32003;

    int32_t bld = 0, *blen = NULL, *bexp = NULL, *bcomp = NULL;
    void *bcf = NULL;
    const int64_t nterms = export_module_f4(malloc,
            &bld, &blen, &bexp, &bcomp, &bcf,
            lens, exps, comps, cfs, NULL,
            fc, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */,
            3 /* nvars */, 2 /* nrows */, 2 /* ngens */,
            12, 1, 0, 2, 1, 0);

    RES_CHECK(nterms > 0, "the rank two module basis is nonempty");
    RES_CHECK(bld >= 3,
            "the rank two basis grew beyond its two generators, so a "
            "cross-generator S-pair really was formed");
    if (bld > 0 && nterms > 0) {
        /* every element must be homogeneous of degree one or two and live
         * in components 1 and 2 only */
        int comps_ok = 1;
        int64_t t;
        for (t = 0; t < nterms; ++t) {
            if (bcomp[t] < 1 || bcomp[t] > 2) {
                comps_ok = 0;
            }
        }
        RES_CHECK(comps_ok, "every basis term has a valid component");

        /* look for an element all of whose terms sit in one component:
         * that is the reduced S-pair (0, x2^2 - x1^2) */
        int found_pure = 0;
        int64_t off = 0;
        for (t = 0; t < bld; ++t) {
            int pure = 1;
            int64_t u;
            for (u = 1; u < blen[t]; ++u) {
                if (bcomp[off+u] != bcomp[off]) {
                    pure = 0;
                }
            }
            if (pure && blen[t] > 1) {
                found_pure = 1;
            }
            off += blen[t];
        }
        RES_CHECK(found_pure,
                "the basis contains the reduced S-pair, supported in a "
                "single component");
    }
    free_module_f4_result_data(free, &blen, &bexp, &bcomp, &bcf);
}

/* Two classical resolutions exercise the module entry point with realistic
 * maps rather than synthetic generators.
 *
 * For the twisted cubic I = (xz-y^2, xw-yz, yw-z^2), Hilbert--Burch gives
 *
 *   0 <- R/I <- R <- R(-2)^3 <- R(-3)^2 <- 0,
 *
 * with columns (z,-y,x) and (-w,z,-y) in the second differential.  The
 * catalecticant A = {{x,y,z},{y,z,w}} has resolution
 *
 *   0 <- coker A <- R^2 <- R(-1)^3 <- R(-3) <- 0,
 *
 * whose last map is the signed vector of the same three quadrics.  Under
 * POT, reducing the columns of A exposes those maximal minors as three
 * pure component-2 basis elements, giving a useful cross-check between the
 * two computations. */
static void res_test_classical_resolutions(
        void
        )
{
    const uint32_t fc = 32003;

    /* --- first differential of the twisted cubic -------------------- */
    const int32_t tc_lens[3] = {2, 2, 2};
    const int32_t tc_exps[24] = {
        1,0,1,0,  0,2,0,0,   /* xz - y^2 */
        1,0,0,1,  0,1,1,0,   /* xw - yz  */
        0,1,0,1,  0,0,2,0    /* yw - z^2 */
    };
    const int32_t tc_comps[6] = {1, 1, 1, 1, 1, 1};
    int32_t tc_cfs[6] = {1, -1, 1, -1, 1, -1};
    const int32_t tc_row_degs[1] = {0};

    const int32_t tc_gb_lens[3] = {2, 2, 2};
    const int32_t tc_gb_exps[24] = {
        0,0,2,0,  0,1,0,1,   /* z^2 - yw */
        0,1,1,0,  1,0,0,1,   /* yz  - xw */
        0,2,0,0,  1,0,1,0    /* y^2 - xz */
    };
    const int32_t tc_gb_comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t tc_gb_cfs[6] = {1, 32002, 1, 32002, 1, 32002};

    int32_t tc_bld = 0, *tc_blen = NULL, *tc_bexp = NULL;
    int32_t *tc_bcomp = NULL;
    void *tc_bcf = NULL;
    const int64_t tc_nterms = export_module_f4(malloc,
            &tc_bld, &tc_blen, &tc_bexp, &tc_bcomp, &tc_bcf,
            tc_lens, tc_exps, tc_comps, tc_cfs, tc_row_degs,
            fc, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */, 4, 1, 3, 12, 1, 0, 2, 1, 0);

    const int tc_shape_ok = tc_bld == 3 && tc_nterms == 6
        && tc_blen != NULL && tc_bexp != NULL
        && tc_bcomp != NULL && tc_bcf != NULL;
    RES_CHECK(tc_shape_ok,
            "the twisted cubic first differential has three quadratic "
            "Groebner generators");
    RES_CHECK(tc_shape_ok
            && memcmp(tc_blen, tc_gb_lens, sizeof(tc_gb_lens)) == 0
            && memcmp(tc_bexp, tc_gb_exps, sizeof(tc_gb_exps)) == 0
            && memcmp(tc_bcomp, tc_gb_comps, sizeof(tc_gb_comps)) == 0
            && memcmp(tc_bcf, tc_gb_cfs, sizeof(tc_gb_cfs)) == 0,
            "the twisted cubic first differential matches its reference "
            "reduced basis");

    /* --- Hilbert--Burch differential of the twisted cubic ------------ */
    const int32_t hb_lens[2] = {3, 3};
    const int32_t hb_exps[24] = {
        0,0,1,0,  0,1,0,0,  1,0,0,0,   /* (z,-y,x)  */
        0,0,0,1,  0,0,1,0,  0,1,0,0    /* (-w,z,-y) */
    };
    const int32_t hb_comps[6] = {1, 2, 3, 1, 2, 3};
    int32_t hb_cfs[6] = {1, -1, 1, -1, 1, -1};
    const int32_t hb_row_degs[3] = {2, 2, 2};

    const int32_t hb_gb_lens[3] = {4, 3, 3};
    const int32_t hb_gb_exps[40] = {
        0,0,2,0,  0,1,0,1,  0,1,1,0,  1,0,0,1,
        0,0,0,1,  0,0,1,0,  0,1,0,0,
        0,0,1,0,  0,1,0,0,  1,0,0,0
    };
    const int32_t hb_gb_comps[10] = {2,2,3,3, 1,2,3, 1,2,3};
    const int32_t hb_gb_cfs[10] = {
        1,32002,32002,1, 1,32002,1, 1,32002,1
    };

    int32_t hb_bld = 0, *hb_blen = NULL, *hb_bexp = NULL;
    int32_t *hb_bcomp = NULL;
    void *hb_bcf = NULL;
    const int64_t hb_nterms = export_module_f4(malloc,
            &hb_bld, &hb_blen, &hb_bexp, &hb_bcomp, &hb_bcf,
            hb_lens, hb_exps, hb_comps, hb_cfs, hb_row_degs,
            fc, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */, 4, 3, 2, 12, 1, 0, 2, 1, 0);

    const int hb_shape_ok = hb_bld == 3 && hb_nterms == 10
        && hb_blen != NULL && hb_bexp != NULL
        && hb_bcomp != NULL && hb_bcf != NULL;
    RES_CHECK(hb_shape_ok,
            "the twisted cubic Hilbert-Burch columns acquire one module "
            "Groebner generator");
    RES_CHECK(hb_shape_ok
            && memcmp(hb_blen, hb_gb_lens, sizeof(hb_gb_lens)) == 0
            && memcmp(hb_bexp, hb_gb_exps, sizeof(hb_gb_exps)) == 0
            && memcmp(hb_bcomp, hb_gb_comps, sizeof(hb_gb_comps)) == 0
            && memcmp(hb_bcf, hb_gb_cfs, sizeof(hb_gb_cfs)) == 0,
            "the twisted cubic Hilbert-Burch differential matches its "
            "reference module basis");

    /* --- presentation of coker {{x,y,z},{y,z,w}} -------------------- */
    const int32_t cat_lens[3] = {2, 2, 2};
    const int32_t cat_exps[24] = {
        1,0,0,0,  0,1,0,0,   /* (x,y) */
        0,1,0,0,  0,0,1,0,   /* (y,z) */
        0,0,1,0,  0,0,0,1    /* (z,w) */
    };
    const int32_t cat_comps[6] = {1, 2, 1, 2, 1, 2};
    int32_t cat_cfs[6] = {1, 1, 1, 1, 1, 1};
    const int32_t cat_row_degs[2] = {0, 0};

    const int32_t cat_gb_lens[6] = {2, 2, 2, 2, 2, 2};
    const int32_t cat_gb_exps[48] = {
        0,0,2,0,  0,1,0,1,
        0,1,1,0,  1,0,0,1,
        0,2,0,0,  1,0,1,0,
        0,0,1,0,  0,0,0,1,
        0,1,0,0,  0,0,1,0,
        1,0,0,0,  0,1,0,0
    };
    const int32_t cat_gb_comps[12] = {2,2, 2,2, 2,2, 1,2, 1,2, 1,2};
    const int32_t cat_gb_cfs[12] = {
        1,32002, 1,32002, 1,32002, 1,1, 1,1, 1,1
    };

    int32_t cat_bld = 0, *cat_blen = NULL, *cat_bexp = NULL;
    int32_t *cat_bcomp = NULL;
    void *cat_bcf = NULL;
    const int64_t cat_nterms = export_module_f4(malloc,
            &cat_bld, &cat_blen, &cat_bexp, &cat_bcomp, &cat_bcf,
            cat_lens, cat_exps, cat_comps, cat_cfs, cat_row_degs,
            fc, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */, 4, 2, 3, 12, 1, 0, 2, 1, 0);

    const int cat_shape_ok = cat_bld == 6 && cat_nterms == 12
        && cat_blen != NULL && cat_bexp != NULL
        && cat_bcomp != NULL && cat_bcf != NULL;
    RES_CHECK(cat_shape_ok,
            "the catalecticant cokernel presentation has six module "
            "Groebner generators");
    RES_CHECK(cat_shape_ok
            && memcmp(cat_blen, cat_gb_lens, sizeof(cat_gb_lens)) == 0
            && memcmp(cat_bexp, cat_gb_exps, sizeof(cat_gb_exps)) == 0
            && memcmp(cat_bcomp, cat_gb_comps, sizeof(cat_gb_comps)) == 0
            && memcmp(cat_bcf, cat_gb_cfs, sizeof(cat_gb_cfs)) == 0,
            "the catalecticant cokernel presentation matches its reference "
            "module basis");

    /* Its first three POT basis elements are the maximal minors, supported
     * in component 2.  They must agree term-for-term with the twisted
     * cubic basis above; these are the entries of the last differential
     * R(-3) -> R(-1)^3 in the cokernel resolution. */
    RES_CHECK(tc_shape_ok && cat_shape_ok
            && memcmp(cat_bexp, tc_bexp, sizeof(tc_gb_exps)) == 0
            && memcmp(cat_bcf, tc_bcf, sizeof(tc_gb_cfs)) == 0,
            "the cokernel resolution recovers the twisted cubic maximal "
            "minors");

    free_module_f4_result_data(
            free, &tc_blen, &tc_bexp, &tc_bcomp, &tc_bcf);
    free_module_f4_result_data(
            free, &hb_blen, &hb_bexp, &hb_bcomp, &hb_bcf);
    free_module_f4_result_data(
            free, &cat_blen, &cat_bexp, &cat_bcomp, &cat_bcf);
}

/* Bad input must be refused rather than trusted. */
static void res_test_module_gb_rejects_bad_input(
        void
        )
{
    const int32_t lens[1] = {1};
    const int32_t exps[3] = {1, 0, 0};
    int32_t cfs[1] = {1};
    const int32_t bad_comps[1] = {3};   /* only 2 rows exist */
    const int32_t ok_comps[1]  = {1};

    int32_t bld = 0, *blen = NULL, *bexp = NULL, *bcomp = NULL;
    void *bcf = NULL;

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "res_selftest: three module input errors are expected next\n");
    }
    RES_CHECK(export_module_f4(malloc, &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, bad_comps, cfs, NULL, 32003, 0, res_strat_p(RES_MORD_POT), NULL, NULL /* stop */,
                3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "an out of range component is rejected");
    RES_CHECK(blen == NULL && bexp == NULL && bcomp == NULL && bcf == NULL,
            "a rejected call allocates nothing");

    RES_CHECK(export_module_f4(malloc, &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, ok_comps, cfs, NULL, 32003, 0,
                res_strat_p(RES_MORD_SCHREYER), NULL, NULL /* stop */,
                3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "the Schreyer order is refused as a *base* order");

    RES_CHECK(export_module_f4(malloc, &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, ok_comps, cfs, NULL, 0 /* char 0 */, 0,
                res_strat_p(RES_MORD_POT), NULL, NULL /* stop */, 3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "characteristic zero is refused");
}

/* --------------------------------------------------------------------- *
 *  Schreyer frames
 *
 *  Every reference table below was read off Macaulay2's
 *
 *      betti res(M, Strategy => Nonminimal)
 *
 *  over the same prime, which is the resolution whose ranks the frame
 *  predicts.  They are written as (level, degree, rank) triples, one per
 *  nonzero entry, terminated by a negative level; Macaulay2 prints the
 *  same numbers with the rows indexed by the slanted degree, degree minus
 *  level. */

static void res_check_frame(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t max_level,
        const int32_t *ref
        )
{
    int32_t i, l, d;

    int32_t nterms = 0;
    for (i = 0; i < ngens; ++i) {
        nterms += lens[i];
    }
    /* msolve reduces the coefficient array in place */
    int32_t *cfs = (int32_t *)malloc((unsigned long)nterms * sizeof(int32_t));
    memcpy(cfs, cfs_src, (unsigned long)nterms * sizeof(int32_t));

    int32_t nlv = 0, maxdeg = 0, *betti = NULL;
    const int64_t nelts = export_module_frame(malloc, &nlv, &maxdeg, &betti,
            lens, exps, comps, cfs, row_degs, 32003, 0 /* drl */,
            res_strat_p(RES_MORD_POT), NULL, nv, nrows, ngens, max_level,
            12 /* ht size */, 1 /* threads */, 0 /* max pairs */,
            2 /* la */, 0 /* info */);

    int64_t rtotal = 0;
    int32_t rlevels = 0;
    for (i = 0; ref[i] >= 0; i += 3) {
        rtotal += ref[i+2];
        if (ref[i] + 1 > rlevels) {
            rlevels = ref[i] + 1;
        }
    }

    RES_CHECK(nelts == rtotal && betti != NULL, what);
    if (betti == NULL || nelts != rtotal) {
        free(cfs);
        free_module_frame_result_data(free, &betti);
        return;
    }
    RES_CHECK(nlv == rlevels, "the frame has the expected number of levels");

    /* every reference entry is present, and nothing else is */
    int ok = nlv == rlevels;
    for (i = 0; ref[i] >= 0 && ok; i += 3) {
        if (ref[i+1] > maxdeg
                || betti[ref[i]*(maxdeg+1) + ref[i+1]] != ref[i+2]) {
            ok = 0;
        }
    }
    int64_t seen = 0;
    for (l = 0; l < nlv; ++l) {
        for (d = 0; d <= maxdeg; ++d) {
            seen += betti[l*(maxdeg+1) + d];
        }
    }
    RES_CHECK(ok && seen == rtotal, "the frame ranks match their reference "
            "Macaulay2 nonminimal Betti table");

    /* The frame is a subcomplex of the Taylor complex on level 1, so it
     * stops by that many levels.  Note this is *not* nv: the frame is a
     * nonminimal resolution and res_test_frame_past_nv is a three variable
     * example that reaches level four. */
    int32_t ngb = 0;
    if (nlv > 1) {
        for (d = 0; d <= maxdeg; ++d) {
            ngb += betti[1*(maxdeg+1) + d];
        }
    }
    RES_CHECK(nlv <= ngb + 1, "the frame stops by the length of the Taylor "
            "complex on the Gröbner basis");

    free(cfs);
    free_module_frame_result_data(free, &betti);
}

/* The Koszul complex, the one case where the frame is forced: the lead
 * terms are the variables themselves. */
static void res_test_frame_koszul(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,1,3,  2,2,3,  3,3,1,  -1
    };

    res_check_frame("the Koszul frame of (x,y,z) has 1+3+3+1 elements",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, ref);
}

/* The Hilbert-Burch case: the Groebner basis is already a minimal
 * generating set, so the frame is the minimal resolution. */
static void res_test_frame_twisted_cubic(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0,  0,2,0,0,   /* xz - y^2 */
        1,0,0,1,  0,1,1,0,   /* xw - yz  */
        0,1,0,1,  0,0,2,0    /* yw - z^2 */
    };
    const int32_t cfs[6]   = {1, -1, 1, -1, 1, -1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,2,3,  2,3,2,  -1
    };

    res_check_frame("the twisted cubic frame is its Hilbert-Burch resolution",
            lens, exps, comps, cfs, NULL, 4, 1, 3, 0, ref);
}

/* A genuinely nonminimal frame: three quadrics that form a complete
 * intersection, whose minimal resolution is 1,3,3,1 in degrees 0,2,4,6
 * while the frame is 1,6,8,3.  This is the case that would silently pass
 * if the frame merely reproduced minimal Betti numbers. */
static void res_test_frame_nonminimal(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,   /* x^2 + yz */
        0,2,0,  1,0,1,   /* y^2 + xz */
        0,0,2,  1,1,0    /* z^2 + xy */
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,
        1,2,3,  1,3,2,  1,4,1,
        2,3,2,  2,4,4,  2,5,2,
        3,5,2,  3,6,1,
        -1
    };

    res_check_frame("a complete intersection of three quadrics has the "
            "nonminimal frame 1,6,8,3",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, ref);
}

/* Discriminating test for the within-block storage order, which is degree
 * ascending and then ring order *descending*.  For (z^2, y^2 z, y^3) that
 * gives the frame 1,3,3,1, which is what Macaulay2 reports; flipping the
 * direction in res_frame_cmp_mono gives 1,3,2 instead -- a valid frame,
 * a smaller one, and the wrong one.  Direction is not about size: on
 * (z, y^3, x y^2, x^2) the sizes come out the other way round, 1,4,5,2
 * for descending against 1,4,6,4,1 for ascending.
 *
 * The classical (x^2, xy, y^3) does *not* discriminate -- both directions
 * give 1,3,2 -- so it is kept only as a second small case. */
static void res_test_frame_block_order(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {2,0,0,  1,1,0,  0,3,0};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,2,2,  1,3,1,  2,3,1,  2,4,1,  -1
    };

    const int32_t dexps[9] = {0,0,2,  0,2,1,  0,3,0};
    const int32_t dref[]   = {
        0,0,1,  1,2,1,  1,3,2,  2,4,2,  2,5,1,  3,5,1,  -1
    };

    res_check_frame("(z^2, y^2 z, y^3) has the frame 1,3,3,1, so blocks "
            "are ordered by descending monomial",
            lens, dexps, comps, cfs, NULL, 3, 1, 3, 0, dref);

    res_check_frame("(x^2, xy, y^3) has the frame 1,3,2",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, ref);
}

/* The frame is a *nonminimal* resolution, so Hilbert's syzygy theorem
 * does not bound its length: (z, y^2, x^2 y, x^3) resolves in three
 * variables but its frame runs to level four, 1,4,6,4,1.  Macaulay2 shows
 * the same table once LengthLimit is raised past its own default of nv,
 * which is exactly the truncation this used to hide behind. */
static void res_test_frame_past_nv(
        void
        )
{
    const int32_t lens[4]   = {1, 1, 1, 1};
    const int32_t exps[12]  = {0,0,1,  0,2,0,  2,1,0,  3,0,0};
    const int32_t cfs[4]    = {1, 1, 1, 1};
    const int32_t comps[4]  = {1, 1, 1, 1};
    const int32_t ref[]     = {
        0,0,1,
        1,1,1,  1,2,1,  1,3,2,
        2,3,1,  2,4,4,  2,5,1,
        3,5,3,  3,6,1,
        4,6,1,
        -1
    };

    res_check_frame("the frame of (z, y^2, x^2 y, x^3) reaches level four "
            "in three variables",
            lens, exps, comps, cfs, NULL, 3, 1, 4, 0, ref);
}

/* A rank two module rather than an ideal: the cokernel of the
 * catalecticant, whose level 0 already has two elements. */
static void res_test_frame_module(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,0,0,  0,1,0,0,   /* (x,y) */
        0,1,0,0,  0,0,1,0,   /* (y,z) */
        0,0,1,0,  0,0,0,1    /* (z,w) */
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
    const int32_t rd[2]    = {0, 0};
    const int32_t ref[]    = {
        0,0,2,
        1,1,3,  1,2,3,
        2,2,3,  2,3,2,
        3,3,1,
        -1
    };

    res_check_frame("the catalecticant cokernel has the frame 2,6,5,1",
            lens, exps, comps, cfs, rd, 4, 2, 3, 0, ref);
}

/* Degree shifts of the ambient free module have to reach the frame, not
 * just the Groebner basis: the columns of {{x^2,y^2},{z,w}} are
 * homogeneous only because the second row sits in degree one.  Getting
 * the shift wrong moves half the table without changing any rank. */
static void res_test_frame_row_degrees(
        void
        )
{
    const int32_t lens[2]  = {2, 2};
    const int32_t exps[16] = {
        2,0,0,0,  0,0,1,0,   /* (x^2, z) */
        0,2,0,0,  0,0,0,1    /* (y^2, w) */
    };
    const int32_t cfs[4]   = {1, 1, 1, 1};
    const int32_t comps[4] = {1, 2, 1, 2};
    const int32_t rd[2]    = {0, 1};
    const int32_t ref[]    = {
        0,0,1,  0,1,1,
        1,2,2,  1,4,1,
        2,4,1,
        -1
    };

    res_check_frame("row degrees shift the frame of coker {{x2,y2},{z,w}}",
            lens, exps, comps, cfs, rd, 4, 2, 2, 0, ref);
}

/* Truncating at a level must give exactly the head of the full frame. */
static void res_test_frame_truncation(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,1,3,  2,2,3,  -1
    };

    res_check_frame("truncating the Koszul frame at level 2 keeps 1,3,3",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 2, ref);
}

static void res_test_frame_rejects_bad_input(
        void
        )
{
    const int32_t lens[1]  = {1};
    const int32_t exps[3]  = {1, 0, 0};
    int32_t cfs[1]         = {1};
    const int32_t comps[1] = {1};

    int32_t nlv = 0, maxdeg = 0, *betti = NULL;

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "res_selftest: two frame input errors are expected next\n");
    }
    RES_CHECK(export_module_frame(malloc, &nlv, &maxdeg, &betti,
                lens, exps, comps, cfs, NULL, 32003, 0,
                res_strat_p(RES_MORD_SCHREYER), NULL,
                3, 1, 1, 0, 12, 1, 0, 2, 0) == 0,
            "the Schreyer order is refused as a base order by the frame");
    RES_CHECK(betti == NULL && nlv == 0,
            "a rejected frame call allocates nothing");

    RES_CHECK(export_module_frame(malloc, &nlv, &maxdeg, &betti,
                lens, exps, comps, cfs, NULL, 0 /* char 0 */, 0,
                res_strat_p(RES_MORD_POT), NULL, 3, 1, 1, 0, 12, 1, 0, 2, 0) == 0,
            "characteristic zero is refused by the frame");
}

static void res_test_frame_rejects_exponent_overflow(
        void
        )
{
    const int32_t lens[2]  = {1, 1};
    const int32_t exps[4]  = {40000, 0,  0, 40000};
    int32_t cfs[2]         = {1, 1};
    const int32_t comps[2] = {1, 1};
    int32_t nlv = 0, maxdeg = 0, *betti = NULL;

    RES_CHECK(export_module_frame(malloc, &nlv, &maxdeg, &betti,
                lens, exps, comps, cfs, NULL, 32003, 0, res_strat_p(RES_MORD_POT), NULL,
                2, 1, 2, 0, 12, 1, 0, 2, 0) == 0,
            "a lifted monomial exceeding the exponent representation is refused");
    RES_CHECK(betti == NULL && nlv == 0,
            "exponent overflow publishes no partial frame");
}

/* --------------------------------------------------------------------- *
 *  Nonminimal differentials and syzygies
 *
 *  Every reference below was read off Macaulay2 over the same prime; the
 *  script that produced them is test/neogb/res/res_reference.m2, which
 *  also checks the things this file cannot check cheaply in C -- that
 *  each complex is exact and resolves the module it claims to, and that
 *  the syzygies generate the same module as Macaulay2's syz.
 * --------------------------------------------------------------------- */

#define RES_FC 32003

/* --------------------------------------------------------------------- *
 *  M6: gradings by a finitely generated abelian group
 *
 *  Three things have to be true for a grading to be more than a label.
 *  The Gröbner basis has to be computed under the *heft* order it induces,
 *  which is what makes the degree by degree schedules terminate; the Betti
 *  table has to come out bucketed by multidegree, which is strictly finer
 *  than the heft degree and is the thing Macaulay2 has no minimalBetti
 *  for; and the heft indexed table has to stay exactly what it was, since
 *  a heft class is a disjoint union of multidegree classes.  The tests
 *  below check all three, and every one of them is discriminating: run
 *  them with the grading replaced by the standard one and they fail.
 * --------------------------------------------------------------------- */

static void res_test_degree_buckets(
        void
        )
{
    /* Z^2, so that two distinct degrees can share a heft */
    const int32_t degs[4] = {1, 0,  0, 1};
    const int32_t heft[2] = {1, 1};
    int32_t buf[2];
    hl_t perm[4];

    res_dgrp_t *g = res_dgrp_new(2, 0, NULL, 2, degs, heft);
    res_dbkt_t *b = g != NULL ? res_dbkt_new(g, 2) : NULL;

    RES_CHECK(g != NULL && b != NULL, "a multidegree bucket set is built");
    if (g == NULL || b == NULL) {
        res_dbkt_free(&b);
        res_dgrp_free(&g);
        return;
    }

    res_deg_t d = {buf};

    buf[0] = 2; buf[1] = 1;
    const hl_t u0 = res_dbkt_insert(b, d);
    buf[0] = 1; buf[1] = 2;
    const hl_t u1 = res_dbkt_insert(b, d);
    buf[0] = 2; buf[1] = 1;
    const hl_t u2 = res_dbkt_insert(b, d);

    RES_CHECK(u0 != u1, "two degrees of the same heft get distinct buckets");
    RES_CHECK(u0 == u2, "re-inserting a degree returns its bucket");
    RES_CHECK(b->ld == 2, "the bucket set holds only the distinct degrees");

    buf[0] = 5; buf[1] = 5;
    RES_CHECK(res_dbkt_find(b, d) == (hl_t)-1,
            "a degree never inserted is not found");

    /* force a rehash, then check nothing was lost */
    int32_t k;
    for (k = 0; k < 40; ++k) {
        buf[0] = k; buf[1] = -k;
        RES_CHECK(res_dbkt_insert(b, d) != (hl_t)-1,
                "bucket insertion survives a rehash");
        res_st_run--;  /* one check for the whole loop, not forty */
    }
    res_st_run++;
    buf[0] = 2; buf[1] = 1;
    RES_CHECK(res_dbkt_find(b, d) == u0,
            "a bucket assigned before a rehash still resolves after it");

    res_dbkt_free(&b);

    /* sorting relabels the buckets into the group's own order */
    b = res_dbkt_new(g, 2);
    if (b != NULL) {
        buf[0] = 3; buf[1] = 0;
        const hl_t a0 = res_dbkt_insert(b, d);
        buf[0] = 1; buf[1] = 0;
        const hl_t a1 = res_dbkt_insert(b, d);
        buf[0] = 2; buf[1] = 0;
        const hl_t a2 = res_dbkt_insert(b, d);

        RES_CHECK(res_dbkt_sort(b, perm) == 0, "the bucket set sorts");
        RES_CHECK(perm[a1] == 0 && perm[a2] == 1 && perm[a0] == 2,
                "sorting puts the buckets in the group's own order");
        RES_CHECK(res_dbkt_at(b, 0).e[0] == 1 && res_dbkt_at(b, 2).e[0] == 3,
                "the sorted buckets hold the sorted degrees");
        buf[0] = 2; buf[1] = 0;
        RES_CHECK(res_dbkt_find(b, d) == 1,
                "lookup still works after a sort");
        res_dbkt_free(&b);
    }

    res_dgrp_free(&g);
}

/* Runs export_module_betti under a grading and checks the multigraded
 * table entry for entry, the heft table entry for entry, and that the one
 * really is the fibrewise sum of the other.  ref is triples-with-a-degree:
 * level, then dlen degree slots, then the count, terminated by a level of
 * -1.  href is the heft table in the usual (level, degree, count) triples. */
static void res_check_graded_betti(
        const char *what,
        const res_grading_t * const grading,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t *ref,
        const int32_t *href,
        const int32_t rpdim,
        const int32_t rreg,
        const int32_t rdim,
        const int64_t rdeg
        )
{
    int32_t i, l;
    int64_t nterms = 0;

    for (i = 0; i < ngens; ++i) {
        nterms += lens[i];
    }
    int32_t *cfs = (int32_t *)malloc((unsigned long)nterms * sizeof(int32_t));
    memcpy(cfs, cfs_src, (unsigned long)nterms * sizeof(int32_t));

    const int32_t dlen = res_grading_len(grading);

    int32_t nlv = 0, maxdeg = 0, shift = 0;
    int32_t pdim = -2, reg = -2, dim = -2;
    int64_t deg  = -1;
    int32_t *tab = NULL, *num = NULL;
    res_mtable_t mt;

    memset(&mt, 0, sizeof(res_mtable_t));

    const int64_t nelts = export_module_betti(malloc, &nlv, &maxdeg, &shift,
            &tab, &num, &pdim, &reg, &dim, &deg, &mt,
            lens, exps, comps, cfs, row_degs, RES_FC, 0 /* drl */,
            NULL, grading, nv, nrows, ngens, 0 /* full */, 1 /* minimal */,
            1 /* verify */, 12, 1, 0, 2, 0);
    free(cfs);

    RES_CHECK(nelts > 0 && tab != NULL && mt.betti != NULL, what);
    if (nelts <= 0 || tab == NULL || mt.betti == NULL) {
        free_module_betti_result_data(free, &tab, &num);
        free_module_mtable_data(free, &mt);
        return;
    }

    RES_CHECK(mt.dlen == dlen,
            "the multigraded table reports the grading's own degree length");
    RES_CHECK(mt.nlevels == nlv,
            "both tables report the same levels");

    /* every reference entry is present in the multigraded table, and the
     * table holds nothing besides them */
    int ok = 1;
    int64_t want = 0, seen = 0;
    for (i = 0; ref[i] >= 0; i += dlen + 2) {
        int32_t u, found = -1;
        want += ref[i + dlen + 1];
        for (u = 0; u < mt.ndegs; ++u) {
            int32_t j, eq = 1;
            for (j = 0; j < dlen; ++j) {
                if (mt.degs[u*dlen + j] != ref[i + 1 + j]) {
                    eq = 0;
                }
            }
            if (eq) {
                found = u;
            }
        }
        if (found < 0 || ref[i] >= mt.nlevels
                || mt.betti[ref[i]*mt.ndegs + found] != ref[i + dlen + 1]) {
            ok = 0;
        }
    }
    for (l = 0; l < mt.nlevels; ++l) {
        for (i = 0; i < mt.ndegs; ++i) {
            seen += mt.betti[l*mt.ndegs + i];
        }
    }
    RES_CHECK(ok && seen == want,
            "the multigraded Betti table matches its Macaulay2 reference");

    /* the heft table, which has to stay exactly what it always was */
    int hok = 1;
    int64_t hwant = 0, hseen = 0;
    for (i = 0; href[i] >= 0; i += 3) {
        hwant += href[i+2];
        if (href[i] >= nlv || href[i+1] > maxdeg
                || tab[href[i]*(maxdeg+1) + href[i+1]] != href[i+2]) {
            hok = 0;
        }
    }
    for (l = 0; l < nlv; ++l) {
        int32_t d;
        for (d = 0; d <= maxdeg; ++d) {
            hseen += tab[l*(maxdeg+1) + d];
        }
    }
    RES_CHECK(hok && hseen == hwant,
            "the heft indexed Betti table matches its Macaulay2 reference");

    /* and the two are the same numbers: summing the multigraded table over
     * a heft fibre has to reproduce the heft table.  This is the check
     * that fails if the rank extraction blocks by the wrong thing. */
    int fib = 1;
    int32_t *acc = (int32_t *)calloc(
            (unsigned long)nlv * (unsigned long)(maxdeg+1), sizeof(int32_t));
    if (acc != NULL) {
        for (l = 0; l < mt.nlevels; ++l) {
            for (i = 0; i < mt.ndegs; ++i) {
                const int32_t h = mt.heft[i];
                if (h < 0 || h > maxdeg) {
                    fib = 0;
                    continue;
                }
                acc[l*(maxdeg+1) + h] += mt.betti[l*mt.ndegs + i];
            }
        }
        for (l = 0; l < nlv; ++l) {
            int32_t d;
            for (d = 0; d <= maxdeg; ++d) {
                if (acc[l*(maxdeg+1) + d] != tab[l*(maxdeg+1) + d]) {
                    fib = 0;
                }
            }
        }
        free(acc);
    }
    RES_CHECK(fib, "the heft table is the multigraded one summed over each "
            "heft fibre");

    /* the same identity for the Hilbert numerator */
    int nok = 1;
    if (num != NULL) {
        int32_t *hn = (int32_t *)calloc(
                (unsigned long)maxdeg + 1, sizeof(int32_t));
        if (hn != NULL) {
            for (i = 0; i < mt.ndegs; ++i) {
                if (mt.heft[i] >= 0 && mt.heft[i] <= maxdeg) {
                    hn[mt.heft[i]] += mt.hilbnum[i];
                }
            }
            for (i = 0; i <= maxdeg; ++i) {
                if (hn[i] != num[i]) {
                    nok = 0;
                }
            }
            free(hn);
        }
    }
    RES_CHECK(nok, "the heft Hilbert numerator is the multigraded one "
            "specialized along the heft");

    RES_CHECK(pdim == rpdim, "the projective dimension matches Macaulay2");
    RES_CHECK(reg == rreg, "the regularity matches Macaulay2");
    RES_CHECK(dim == rdim, "the Krull dimension matches Macaulay2");
    RES_CHECK(deg == rdeg, "the degree matches Macaulay2");

    free_module_betti_result_data(free, &tab, &num);
    free_module_mtable_data(free, &mt);
}

/* R = k[x,y,z] with deg x = 1, deg y = 2, deg z = 3, and I = (x^2 y, yz, xz).
 * The heft degree is the weighted degree, so the Gröbner basis itself is
 * computed in weighted DRL; with the weights dropped the generators would
 * sit in degrees 3, 2, 2 instead of 4, 5, 4 and every number below would
 * change.  Macaulay2 reference in test/neogb/res/res_reference.m2. */
static void res_test_weighted_betti(
        void
        )
{
    const int32_t vdegs[3] = {1, 2, 3};
    const res_grading_t grading = {1, 0, NULL, vdegs, NULL};

    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {2,1,0,  0,1,1,  1,0,1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t cfs[3]   = {1, 1, 1};

    /* level, degree, count -- r = 1, so the two tables coincide */
    const int32_t ref[]  = {0,0,1,  1,4,2,  1,5,1,  2,6,1,  2,7,1,  -1};
    const int32_t href[] = {0,0,1,  1,4,2,  1,5,1,  2,6,1,  2,7,1,  -1};

    res_check_graded_betti(
            "a weighted grading resolves through the C entry point",
            &grading, lens, exps, comps, cfs, NULL, 3, 1, 3,
            ref, href, 2 /* pdim */, 5 /* reg */, 1 /* dim */, 14 /* deg */);
}

/* The same weighted ring, but an ideal whose *order* depends on the
 * weights: y^4 - x^5 z and z^2 - x^2 y^2 are weighted homogeneous of
 * degrees 8 and 6 and are not homogeneous at all for the standard grading,
 * so the weights are not decoration here -- they decide which term of each
 * generator leads, hence the Gröbner basis, hence the frame.  Under plain
 * degree reverse lexicographic the lead of the first would be x^5 z rather
 * than y^4.
 *
 * This is the test that fails if the heft degree is not carried in ev[DEG]:
 * the input then looks inhomogeneous and is refused outright.  Macaulay2
 * reference in test/neogb/res/res_reference.m2. */
static void res_test_weighted_order(
        void
        )
{
    const int32_t vdegs[3] = {1, 2, 3};
    const res_grading_t grading = {1, 0, NULL, vdegs, NULL};

    const int32_t lens[2]  = {2, 2};
    const int32_t exps[12] = {
        0,4,0, 5,0,1,   /* y^4 - x^5 z,   weighted degree 8 */
        0,0,2, 2,2,0    /* z^2 - x^2 y^2, weighted degree 6 */
    };
    const int32_t comps[4] = {1, 1, 1, 1};
    const int32_t cfs[4]   = {1, RES_FC-1, 1, RES_FC-1};

    const int32_t ref[]  = {0,0,1,  1,6,1,  1,8,1,  2,14,1,  -1};
    const int32_t href[] = {0,0,1,  1,6,1,  1,8,1,  2,14,1,  -1};

    res_check_graded_betti(
            "a weighted grading decides the order, not just the degrees",
            &grading, lens, exps, comps, cfs, NULL, 3, 1, 2,
            ref, href, 2 /* pdim */, 12 /* reg */, 1 /* dim */, 48 /* deg */);
}

/* T = k[a,b,c,d] over P^1 x P^1 and J = (ac, bd, ad).  This is the case the
 * whole milestone is about: level 2 carries two generators of the same heft
 * degree 3 sitting in different multidegrees, (1,2) and (2,1), so the heft
 * table reports a single 2 where the multigraded table reports two 1s.  A
 * rank extraction that blocked by heft degree would still get the heft
 * table right and could not produce this one at all. */
static void res_test_multigraded_betti(
        void
        )
{
    const int32_t vdegs[8] = {1,0,  1,0,  0,1,  0,1};
    const int32_t heft[2]  = {1, 1};
    const res_grading_t grading = {2, 0, NULL, vdegs, heft};

    const int32_t lens[3]   = {1, 1, 1};
    const int32_t exps[12]  = {1,0,1,0,  0,1,0,1,  1,0,0,1};
    const int32_t comps[3]  = {1, 1, 1};
    const int32_t cfs[3]    = {1, 1, 1};

    /* level, deg0, deg1, count */
    const int32_t ref[]  = {0,0,0,1,  1,1,1,3,  2,1,2,1,  2,2,1,1,  -1};
    const int32_t href[] = {0,0,1,  1,2,3,  2,3,2,  -1};

    res_check_graded_betti(
            "a Z^2 grading resolves through the C entry point",
            &grading, lens, exps, comps, cfs, NULL, 4, 1, 3,
            ref, href, 2 /* pdim */, 1 /* reg */, 2 /* dim */, 3 /* deg */);
}

/* The same ring and module presented with a degree shift: N = coker[ac bd]
 * over R^1, whose resolution is 1, 2, 1 in multidegrees (0,0), (1,1),
 * (2,2). */
static void res_test_multigraded_module(
        void
        )
{
    const int32_t vdegs[8] = {1,0,  1,0,  0,1,  0,1};
    const int32_t heft[2]  = {1, 1};
    const res_grading_t grading = {2, 0, NULL, vdegs, heft};

    const int32_t lens[2]  = {1, 1};
    const int32_t exps[8]  = {1,0,1,0,  0,1,0,1};
    const int32_t comps[2] = {1, 1};
    const int32_t cfs[2]   = {1, 1};

    const int32_t ref[]  = {0,0,0,1,  1,1,1,2,  2,2,2,1,  -1};
    const int32_t href[] = {0,0,1,  1,2,2,  2,4,1,  -1};

    res_check_graded_betti(
            "a Z^2 graded cokernel resolves through the C entry point",
            &grading, lens, exps, comps, cfs, NULL, 4, 1, 2,
            ref, href, 2 /* pdim */, 2 /* reg */, 2 /* dim */, 4 /* deg */);
}

/* Torsion.  R = k[x,y] graded by Z (+) Z/2 with deg x = (1,0) and
 * deg y = (1,1); the heft ignores the torsion, so x and y both have heft
 * degree 1 and the ring order is msolve's ordinary DRL, but the two
 * generators of I = (x^2, xy) land in different multidegrees, (2,0) and
 * (2,1).  Macaulay2 has no grading group with torsion, so the reference
 * here is by hand: the resolution is the Koszul-like 1, 2, 1 with the
 * syzygy y*x^2 - x*(xy) in multidegree (3,1).
 *
 * This is the test that the torsion arithmetic is real rather than carried
 * along: with nt = 0 the two level 1 generators share a bucket and the
 * multigraded table collapses onto the heft one. */
static void res_test_torsion_betti(
        void
        )
{
    const int32_t tord[1]  = {2};
    const int32_t vdegs[4] = {1,0,  1,1};
    const int32_t heft[1]  = {1};
    const res_grading_t grading = {1, 1, tord, vdegs, heft};

    const int32_t lens[2]  = {1, 1};
    const int32_t exps[4]  = {2,0,  1,1};
    const int32_t comps[2] = {1, 1};
    const int32_t cfs[2]   = {1, 1};

    /* level, free part, torsion residue, count */
    const int32_t ref[]  = {0,0,0,1,  1,2,0,1,  1,2,1,1,  2,3,1,1,  -1};
    const int32_t href[] = {0,0,1,  1,2,2,  2,3,1,  -1};

    res_check_graded_betti(
            "a grading with torsion resolves through the C entry point",
            &grading, lens, exps, comps, cfs, NULL, 2, 1, 2,
            ref, href, 2 /* pdim */, 1 /* reg */, 1 /* dim */, 1 /* deg */);
}

/* The standard grading spelled out as a res_grading_t has to give exactly
 * what NULL gives -- not merely the same Betti numbers, but the same
 * arrays.  This is what makes every test above a test of the *grading*
 * rather than of a second code path. */
static void res_test_explicit_standard_grading(
        void
        )
{
    const int32_t vdegs[4] = {1, 1, 1, 1};
    const res_grading_t grading = {1, 0, NULL, vdegs, NULL};

    /* the twisted cubic, whose resolution is 1, 3, 2 */
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0, 0,2,0,0,   /* xz - y^2 */
        1,0,0,1, 0,1,1,0,   /* xw - yz  */
        0,1,0,1, 0,0,2,0    /* yw - z^2 */
    };
    const int32_t comps[6] = {1,1, 1,1, 1,1};
    const int32_t cfs[6]   = {1, RES_FC-1, 1, RES_FC-1, 1, RES_FC-1};

    int32_t a_nlv = 0, a_max = 0, b_nlv = 0, b_max = 0;
    int32_t *a_tab = NULL, *a_num = NULL, *b_tab = NULL, *b_num = NULL;
    int32_t cfsa[6], cfsb[6];

    memcpy(cfsa, cfs, sizeof(cfs));
    memcpy(cfsb, cfs, sizeof(cfs));

    const int64_t na = export_module_betti(malloc, &a_nlv, &a_max, NULL,
            &a_tab, &a_num, NULL, NULL, NULL, NULL, NULL,
            lens, exps, comps, cfsa, NULL, RES_FC, 0, NULL, NULL,
            4, 1, 3, 0, 1, 0, 12, 1, 0, 2, 0);
    const int64_t nb = export_module_betti(malloc, &b_nlv, &b_max, NULL,
            &b_tab, &b_num, NULL, NULL, NULL, NULL, NULL,
            lens, exps, comps, cfsb, NULL, RES_FC, 0, NULL, &grading,
            4, 1, 3, 0, 1, 0, 12, 1, 0, 2, 0);

    RES_CHECK(na == nb && a_nlv == b_nlv && a_max == b_max,
            "the explicit standard grading agrees with NULL on the shape");
    if (a_tab != NULL && b_tab != NULL && a_nlv == b_nlv && a_max == b_max) {
        RES_CHECK(memcmp(a_tab, b_tab,
                    (unsigned long)a_nlv * (a_max+1) * sizeof(int32_t)) == 0,
                "the explicit standard grading gives an identical table");
        RES_CHECK(memcmp(a_num, b_num,
                    (unsigned long)(a_max+1) * sizeof(int32_t)) == 0,
                "the explicit standard grading gives an identical numerator");
    }

    free_module_betti_result_data(free, &a_tab, &a_num);
    free_module_betti_result_data(free, &b_tab, &b_num);
}

/* A grading the engine cannot run, and input the grading cannot grade. */
static void res_test_grading_rejects_bad_input(
        void
        )
{
    const int32_t lens[2]  = {1, 1};
    const int32_t exps[4]  = {2,0,  1,1};
    const int32_t comps[2] = {1, 1};
    int32_t cfs[2]         = {1, 1};

    const int32_t vdegs[4] = {1,0,  1,1};
    const int32_t heft[1]  = {1};

    int32_t nlv = 0, maxdeg = 0;

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "res_selftest: five grading errors are expected next\n");
    }

    /* a torsion factor of order one is not a torsion factor */
    const int32_t bad_tord[1] = {1};
    const res_grading_t g1 = {1, 1, bad_tord, vdegs, heft};
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, RES_FC, 0, NULL, &g1,
                2, 1, 2, 0, 0, 0, 12, 1, 0, 2, 0) == 0,
            "a torsion factor of order below two is refused");

    /* torsion without the orders to reduce by */
    const res_grading_t g2 = {1, 1, NULL, vdegs, heft};
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, RES_FC, 0, NULL, &g2,
                2, 1, 2, 0, 0, 0, 12, 1, 0, 2, 0) == 0,
            "torsion factors without their orders are refused");

    /* no degrees at all */
    const res_grading_t g3 = {1, 0, NULL, NULL, heft};
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, RES_FC, 0, NULL, &g3,
                2, 1, 2, 0, 0, 0, 12, 1, 0, 2, 0) == 0,
            "a grading without variable degrees is refused");

    /* a variable of heft degree zero: the degree by degree schedule would
     * never terminate, so this is refused rather than attempted */
    const int32_t zdegs[2] = {1, 0};
    const res_grading_t g4 = {1, 0, NULL, zdegs, heft};
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, RES_FC, 0, NULL, &g4,
                2, 1, 2, 0, 0, 0, 12, 1, 0, 2, 0) == 0,
            "a variable of heft degree zero is refused");

    /* Input that is homogeneous for the heft but not for the grading:
     * x^2 + xy has heft degree 2 throughout, and multidegrees (2,0) and
     * (1,1), so it is graded for Z and not for Z^2.  Catching this needs
     * the multidegree of every term, which is exactly what the heft
     * homogeneity test the engine had before cannot see. */
    const int32_t mdegs[4] = {1,0,  0,1};
    const int32_t mheft[2] = {1, 1};
    const res_grading_t g5 = {2, 0, NULL, mdegs, mheft};
    const int32_t ilens[1]  = {2};
    const int32_t iexps[4]  = {2,0,  1,1};
    const int32_t icomps[2] = {1, 1};
    int32_t icfs[2]         = {1, 1};
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                ilens, iexps, icomps, icfs, NULL, RES_FC, 0, NULL, &g5,
                2, 1, 1, 0, 0, 0, 12, 1, 0, 2, 0) == 0,
            "input homogeneous only for the heft is refused as multigraded");

    /* ... and the very same input *is* accepted for the heft grading it is
     * homogeneous for, which is what makes the refusal above about the
     * grading rather than about the input */
    icfs[0] = 1; icfs[1] = 1;
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL,
                ilens, iexps, icomps, icfs, NULL, RES_FC, 0, NULL, NULL,
                2, 1, 1, 0, 0, 0, 12, 1, 0, 2, 0) > 0,
            "the same input is accepted under the standard grading");
}

/* The resolution kept alive reports multidegrees too, and they have to be
 * the ones the one shot entry point tabulates. */
static void res_test_comp_multidegrees(
        void
        )
{
    const int32_t vdegs[8] = {1,0,  1,0,  0,1,  0,1};
    const int32_t heft[2]  = {1, 1};
    const res_grading_t grading = {2, 0, NULL, vdegs, heft};

    const int32_t lens[3]   = {1, 1, 1};
    const int32_t exps[12]  = {1,0,1,0,  0,1,0,1,  1,0,0,1};
    const int32_t comps[3]  = {1, 1, 1};
    int32_t cfs[3]          = {1, 1, 1};

    res_comp_t *c = res_comp_new(lens, exps, comps, cfs, NULL, RES_FC, 0,
            NULL, &grading, 4, 1, 3, 0, 12, 1, 0, 2, 0);

    RES_CHECK(c != NULL, "a multigraded resolution handle is built");
    if (c == NULL) {
        return;
    }

    RES_CHECK(res_comp_glen(c) == 2,
            "the handle reports the grading's degree length");

    int32_t shift[2] = {-1, -1};
    RES_CHECK(res_comp_multidegshift(c, shift) == 0
            && shift[0] == 0 && shift[1] == 0,
            "an unshifted ambient module has a zero multidegree shift");

    /* level 1 is the Gröbner basis of (ac, bd, ad), every element of which
     * sits in multidegree (1,1) */
    const int32_t rk = res_comp_rank(c, 1);
    RES_CHECK(rk == 3, "the multigraded frame has three level 1 elements");
    if (rk > 0) {
        int32_t *md = (int32_t *)malloc(
                (unsigned long)rk * 2 * sizeof(int32_t));
        int32_t *hd = (int32_t *)malloc(
                (unsigned long)rk * sizeof(int32_t));
        int ok = md != NULL && hd != NULL
            && res_comp_multidegrees(c, 1, md) == 0
            && res_comp_degrees(c, 1, hd) == 0;
        int32_t k;
        for (k = 0; ok && k < rk; ++k) {
            if (md[2*k] != 1 || md[2*k+1] != 1 || hd[k] != 2) {
                ok = 0;
            }
        }
        RES_CHECK(ok, "every level 1 generator sits in multidegree (1,1), "
                "of heft degree 2");
        free(md);
        free(hd);
    }

    res_comp_free(&c);
    RES_CHECK(c == NULL, "the multigraded handle is nulled on free");
}


typedef struct res_res_t res_res_t;
struct res_res_t
{
    int32_t  nlv;
    int32_t *ranks;
    int32_t *degs;
    int32_t *dlen;
    int32_t *dexp;
    int32_t *dcomp;
    int32_t *dcf;
    int64_t  nterms;
    int32_t *cfs;   /* the working copy msolve reduced in place */
};

static int64_t res_run_resolution_stop(
        res_res_t *r,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t max_level,
        const int32_t syz_of,
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop
        )
{
    int32_t i, nt = 0;
    void *cf = NULL;

    memset(r, 0, sizeof(res_res_t));
    for (i = 0; i < ngens; ++i) {
        nt += lens[i];
    }
    r->cfs = (int32_t *)malloc((unsigned long)nt * sizeof(int32_t));
    memcpy(r->cfs, cfs_src, (unsigned long)nt * sizeof(int32_t));

    r->nterms = export_module_resolution(malloc, &r->nlv, &r->ranks,
            &r->degs, &r->dlen, &r->dexp, &r->dcomp, &cf,
            lens, exps, comps, r->cfs, row_degs, RES_FC, 0 /* drl */,
            strat, grading, stop, nv, nrows, ngens, max_level, syz_of,
            1 /* verify d o d = 0 exactly */,
            12 /* ht size */, 1 /* threads */, 0 /* max pairs */,
            2 /* la */, 0 /* info */);
    r->dcf = (int32_t *)cf;

    return r->nterms;
}

static int64_t res_run_resolution(
        res_res_t *r,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t max_level,
        const int32_t syz_of,
        const res_strat_t * const strat
        )
{
    return res_run_resolution_stop(r, lens, exps, comps, cfs_src, row_degs,
            nv, nrows, ngens, max_level, syz_of, strat, NULL, NULL);
}

static void res_free_resolution(
        res_res_t *r
        )
{
    void *cf = r->dcf;
    free_module_resolution_result_data(free, &r->ranks, &r->degs,
            &r->dlen, &r->dexp, &r->dcomp, &cf);
    r->dcf = NULL;
    free(r->cfs);
    r->cfs = NULL;
}

/* d_{i-1} o d_i = 0, checked by substituting a point of F_p^n into every
 * entry and multiplying the resulting scalar matrices.  A nonzero
 * composite is a nonzero polynomial matrix of degree well under a
 * hundred, so one point of F_32003 already misses it with probability
 * below 1/300; two independent points make that negligible, and this way
 * the test needs no polynomial arithmetic of its own -- it shares no code
 * at all with the engine it is checking. */
static int res_composite_is_zero(
        const res_res_t * const r,
        const int32_t nv,
        const uint32_t * const pt
        )
{
    int32_t i, k, j, c;
    int64_t cl = 0, ce = 0, cc = 0;
    int ok = 1;

    if (r->nlv < 3) {
        return 1;
    }

    /* mat[i] is the ranks[i-1] x ranks[i] matrix of d_i at the point */
    uint32_t **mat = (uint32_t **)calloc(
            (unsigned long)r->nlv, sizeof(uint32_t *));

    for (i = 1; i < r->nlv; ++i) {
        const int32_t nr = r->ranks[i-1];
        const int32_t nc = r->ranks[i];
        mat[i] = (uint32_t *)calloc(
                (unsigned long)(nr > 0 && nc > 0 ? nr * nc : 1),
                sizeof(uint32_t));
        for (k = 0; k < nc; ++k) {
            const int32_t len = r->dlen[cl++];
            for (j = 0; j < len; ++j) {
                uint64_t v = (uint64_t)(uint32_t)r->dcf[cc];
                for (c = 0; c < nv; ++c) {
                    int32_t e = r->dexp[ce++];
                    while (e-- > 0) {
                        v = (v * (uint64_t)pt[c]) % RES_FC;
                    }
                }
                const int32_t row = r->dcomp[cc] - 1;
                uint32_t *slot = mat[i] + (int64_t)row * nc + k;
                *slot = (uint32_t)(((uint64_t)*slot + v) % RES_FC);
                cc++;
            }
        }
    }

    for (i = 2; i < r->nlv; ++i) {
        const int32_t nr = r->ranks[i-2];
        const int32_t nm = r->ranks[i-1];
        const int32_t nc = r->ranks[i];
        for (j = 0; j < nr; ++j) {
            for (k = 0; k < nc; ++k) {
                uint64_t s = 0;
                for (c = 0; c < nm; ++c) {
                    s += (uint64_t)mat[i-1][(int64_t)j * nm + c]
                        * (uint64_t)mat[i][(int64_t)c * nc + k];
                    s %= RES_FC;
                }
                if (s != 0) {
                    ok = 0;
                }
            }
        }
    }

    for (i = 1; i < r->nlv; ++i) {
        free(mat[i]);
    }
    free(mat);

    return ok;
}

/* ranks, the multiset of degrees at every level, and d o d = 0 */
static void res_check_resolution(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t max_level,
        const int32_t syz_of,
        const int32_t *ranks,   /* terminated by a negative entry */
        const int32_t *degs     /* sum(ranks) entries, ascending per level */
        )
{
    int32_t i, k;
    res_res_t r;

    const int64_t n = res_run_resolution(&r, lens, exps, comps, cfs,
            row_degs, nv, nrows, ngens, max_level, syz_of, res_strat_p(RES_MORD_POT));

    int32_t nlv = 0;
    while (ranks[nlv] >= 0) {
        nlv++;
    }

    RES_CHECK(n > 0 && r.nlv == nlv, what);
    if (n <= 0 || r.nlv != nlv) {
        res_free_resolution(&r);
        return;
    }

    int ok = 1;
    for (i = 0; i < nlv; ++i) {
        if (r.ranks[i] != ranks[i]) {
            ok = 0;
        }
    }
    RES_CHECK(ok, "the free modules of the resolution have the expected "
            "ranks");

    /* degrees, compared as a sorted list per level: the order inside a
     * level is the frame's storage order, which is not part of the
     * contract */
    int okd = 1;
    int32_t off = 0;
    for (i = 0; i < nlv && okd; ++i) {
        for (k = 0; k < ranks[i]; ++k) {
            int32_t seen = 0, j;
            for (j = 0; j < ranks[i]; ++j) {
                if (r.degs[off+j] == degs[off+k]) {
                    seen++;
                }
            }
            int32_t want = 0;
            for (j = 0; j < ranks[i]; ++j) {
                if (degs[off+j] == degs[off+k]) {
                    want++;
                }
            }
            if (seen != want) {
                okd = 0;
            }
        }
        off += ranks[i];
    }
    RES_CHECK(okd, "every free module is generated in the expected degrees");

    const uint32_t pt1[8] = {2, 3, 5, 7, 11, 13, 17, 19};
    const uint32_t pt2[8] = {9161, 421, 30011, 7, 12289, 251, 4099, 65};
    RES_CHECK(res_composite_is_zero(&r, nv, pt1)
            && res_composite_is_zero(&r, nv, pt2),
            "the differential composes to zero");

    res_free_resolution(&r);
}

/* Term for term against a reference, for the cases small enough that the
 * whole differential can be written down. */
static void res_check_differential(
        const char *what,
        const res_res_t * const r,
        const int32_t nv,
        const int32_t *dlen,
        const int32_t *dexp,
        const int32_t *dcomp,
        const int32_t *dcf
        )
{
    int32_t i;
    int64_t ncols = 0, nterms = 0;

    for (i = 1; i < r->nlv; ++i) {
        ncols += r->ranks[i];
    }
    for (i = 0; i < ncols; ++i) {
        nterms += dlen[i];
    }

    int ok = nterms == r->nterms;
    for (i = 0; i < ncols && ok; ++i) {
        if (r->dlen[i] != dlen[i]) {
            ok = 0;
        }
    }
    for (i = 0; i < nterms && ok; ++i) {
        if (r->dcomp[i] != dcomp[i] || r->dcf[i] != dcf[i]) {
            ok = 0;
        }
    }
    for (i = 0; i < nterms * nv && ok; ++i) {
        if (r->dexp[i] != dexp[i]) {
            ok = 0;
        }
    }
    RES_CHECK(ok, what);
}

/* The Koszul complex on (x,y,z), the one resolution that is forced all
 * the way down to the sign of every entry. */
static void res_test_resolution_koszul(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ranks[]  = {1, 3, 3, 1, -1};
    const int32_t degs[]   = {0,  1,1,1,  2,2,2,  3};

    /* d_1 = (x, y, z), the Gröbner basis in the frame's storage order,
     * which sorts a block by descending monomial */
    const int32_t rlen[7]  = {1, 1, 1,  2, 2, 2,  3};
    const int32_t rexp[]   = {
        1,0,0,   0,1,0,   0,0,1,
        1,0,0, 0,1,0,   1,0,0, 0,0,1,   0,1,0, 0,0,1,
        1,0,0, 0,1,0, 0,0,1
    };
    const int32_t rcomp[]  = {1, 1, 1,  2,1,  3,1,  3,2,  3,2,1};
    const int32_t p        = RES_FC;
    const int32_t rcf[]    = {
        1, 1, 1,
        1, p-1,   1, p-1,   1, p-1,
        1, p-1, 1
    };

    res_check_resolution("the Koszul resolution of (x,y,z) is 1,3,3,1",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, RES_SYZ_OF_GB,
            ranks, degs);

    res_res_t r;
    res_run_resolution(&r, lens, exps, comps, cfs, NULL,
            3, 1, 3, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT));
    if (r.nterms == 12 && r.nlv == 4) {
        res_check_differential("the Koszul differential matches its "
                "reference entry for entry", &r, 3, rlen, rexp, rcomp, rcf);
    } else {
        RES_CHECK(0, "the Koszul differential has twelve terms");
    }
    res_free_resolution(&r);
}

/* Hilbert--Burch: the frame is already the minimal resolution, and the
 * second differential is the classical 3 x 2 matrix. */
static void res_test_resolution_twisted_cubic(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0,  0,2,0,0,   /* xz - y^2 */
        1,0,0,1,  0,1,1,0,   /* xw - yz  */
        0,1,0,1,  0,0,2,0    /* yw - z^2 */
    };
    const int32_t cfs[6]   = {1, -1, 1, -1, 1, -1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ranks[]  = {1, 3, 2, -1};
    const int32_t degs[]   = {0,  2,2,2,  3,3};

    res_check_resolution("the twisted cubic resolves as 1,3,2",
            lens, exps, comps, cfs, NULL, 4, 1, 3, 0, RES_SYZ_OF_GB,
            ranks, degs);

    /* d_1 is y^2-xz, yz-xw, z^2-yw and d_2 the Hilbert--Burch columns
     * y e_2 - z e_1 - x e_3 and y e_3 - z e_2 + w e_1 */
    const int32_t p        = RES_FC;
    const int32_t rlen[5]  = {2, 2, 2,  3, 3};
    const int32_t rexp[]   = {
        0,2,0,0, 1,0,1,0,
        0,1,1,0, 1,0,0,1,
        0,0,2,0, 0,1,0,1,
        0,1,0,0, 0,0,1,0, 1,0,0,0,
        0,1,0,0, 0,0,1,0, 0,0,0,1
    };
    const int32_t rcomp[]  = {1,1, 1,1, 1,1,  2,1,3,  3,2,1};
    const int32_t rcf[]    = {
        1,p-1, 1,p-1, 1,p-1,
        1,p-1,p-1,
        1,p-1,1
    };

    res_res_t r;
    res_run_resolution(&r, lens, exps, comps, cfs, NULL,
            4, 1, 3, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT));
    if (r.nterms == 12 && r.nlv == 3) {
        res_check_differential("the twisted cubic differential is its "
                "Hilbert-Burch matrix", &r, 4, rlen, rexp, rcomp, rcf);
    } else {
        RES_CHECK(0, "the twisted cubic differential has twelve terms");
    }
    res_free_resolution(&r);
}

/* The case where the resolution is genuinely nonminimal: the frame 1,6,8,3
 * of a complete intersection of three quadrics whose minimal resolution is
 * 1,3,3,1.  This is also the example that forces the schedule -- level 2
 * in degree 3 reduces against a level 1 generator of degree 3, with
 * multiplier 1, so a slanted degree driver would need it before it had
 * been computed. */
static void res_test_resolution_nonminimal(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,   /* x^2 + yz */
        0,2,0,  1,0,1,   /* y^2 + xz */
        0,0,2,  1,1,0    /* z^2 + xy */
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ranks[]  = {1, 6, 8, 3, -1};
    const int32_t degs[]   = {
        0,
        2,2,2,3,3,4,
        3,3,4,4,4,4,5,5,
        5,5,6
    };

    res_check_resolution("a complete intersection of three quadrics has "
            "the nonminimal resolution 1,6,8,3",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, RES_SYZ_OF_GB,
            ranks, degs);
}

/* A rank two module rather than an ideal. */
static void res_test_resolution_module(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,0,0,  0,1,0,0,   /* (x,y) */
        0,1,0,0,  0,0,1,0,   /* (y,z) */
        0,0,1,0,  0,0,0,1    /* (z,w) */
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
    const int32_t rd[2]    = {0, 0};
    const int32_t ranks[]  = {2, 6, 5, 1, -1};
    const int32_t degs[]   = {0,0,  1,1,1,2,2,2,  2,2,2,3,3,  3};

    res_check_resolution("the catalecticant cokernel resolves as 2,6,5,1",
            lens, exps, comps, cfs, rd, 4, 2, 3, 0, RES_SYZ_OF_GB,
            ranks, degs);
}

/* Degree shifts of the ambient free module have to reach the coefficients
 * too, not just the frame: this module is homogeneous only because the
 * second row sits in degree one. */
static void res_test_resolution_row_degrees(
        void
        )
{
    const int32_t lens[2]  = {2, 2};
    const int32_t exps[16] = {
        2,0,0,0,  0,0,1,0,   /* (x^2, z) */
        0,2,0,0,  0,0,0,1    /* (y^2, w) */
    };
    const int32_t cfs[4]   = {1, 1, 1, 1};
    const int32_t comps[4] = {1, 2, 1, 2};
    const int32_t rd[2]    = {0, 1};
    const int32_t ranks[]  = {2, 3, 1, -1};
    const int32_t degs[]   = {0,1,  2,2,4,  4};

    res_check_resolution("row degrees shift the resolution of "
            "coker {{x2,y2},{z,w}}",
            lens, exps, comps, cfs, rd, 4, 2, 2, 0, RES_SYZ_OF_GB,
            ranks, degs);
}

/* Truncation must give exactly the head of the full resolution. */
static void res_test_resolution_truncation(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ranks1[] = {1, 3, -1};
    const int32_t degs1[]  = {0,  1,1,1};
    const int32_t ranks2[] = {1, 3, 3, -1};
    const int32_t degs2[]  = {0,  1,1,1,  2,2,2};

    res_check_resolution("truncating the Koszul resolution at level 1 "
            "leaves the Gröbner basis",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 1, RES_SYZ_OF_GB,
            ranks1, degs1);
    res_check_resolution("truncating at level 2 leaves the single syzygy "
            "matrix",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 2, RES_SYZ_OF_GB,
            ranks2, degs2);
}

/* The ranks of the resolution are the frame ranks, which M3 already
 * checked against Macaulay2; agreeing with them ties the two entry points
 * together rather than letting them drift apart. */
static void res_test_resolution_matches_frame(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,
        0,2,0,  1,0,1,
        0,0,2,  1,1,0
    };
    const int32_t cfs_src[6] = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6]   = {1, 1, 1, 1, 1, 1};

    int32_t cfs[6];
    memcpy(cfs, cfs_src, sizeof(cfs_src));

    int32_t nlv = 0, maxdeg = 0, *betti = NULL;
    export_module_frame(malloc, &nlv, &maxdeg, &betti,
            lens, exps, comps, cfs, NULL, RES_FC, 0, res_strat_p(RES_MORD_POT), NULL,
            3, 1, 3, 0, 12, 1, 0, 2, 0);

    res_res_t r;
    res_run_resolution(&r, lens, exps, comps, cfs_src, NULL,
            3, 1, 3, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT));

    int ok = betti != NULL && r.nterms > 0 && r.nlv == nlv;
    int32_t i, d, off = 0;
    for (i = 0; i < nlv && ok; ++i) {
        int32_t sum = 0;
        for (d = 0; d <= maxdeg; ++d) {
            sum += betti[i*(maxdeg+1) + d];
            /* and degree by degree, not just in total */
            int32_t seen = 0, k;
            for (k = 0; k < r.ranks[i]; ++k) {
                if (r.degs[off+k] == d) {
                    seen++;
                }
            }
            if (seen != betti[i*(maxdeg+1) + d]) {
                ok = 0;
            }
        }
        if (sum != r.ranks[i]) {
            ok = 0;
        }
        off += r.ranks[i];
    }
    RES_CHECK(ok, "the resolution has exactly the frame's ranks, degree by "
            "degree");

    free_module_frame_result_data(free, &betti);
    res_free_resolution(&r);
}

/* --- syzygies of the input generators -------------------------------- */

/* The three Koszul relations, and nothing else: here the input generators
 * are already a Gröbner basis and already minimal, so the graph module
 * trick has to reproduce Macaulay2's syz exactly. */
static void res_test_syz_of_input_koszul(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ranks[]  = {1, 3, 3, -1};
    const int32_t degs[]   = {0,  1,1,1,  2,2,2};
    const int32_t p        = RES_FC;

    /* d_1 is the caller's matrix, in the caller's order */
    const int32_t rlen[6]  = {1, 1, 1,  2, 2, 2};
    const int32_t rexp[]   = {
        1,0,0,   0,1,0,   0,0,1,
        0,0,1, 0,1,0,   0,0,1, 1,0,0,   0,1,0, 1,0,0
    };
    const int32_t rcomp[]  = {1, 1, 1,  2,3,  1,3,  1,2};
    const int32_t rcf[]    = {1, 1, 1,  1,p-1,  1,p-1,  1,p-1};

    res_check_resolution("the syzygies of (x,y,z) are the three Koszul "
            "relations",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 2, RES_SYZ_OF_INPUT,
            ranks, degs);

    res_res_t r;
    res_run_resolution(&r, lens, exps, comps, cfs, NULL,
            3, 1, 3, 2, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT));
    if (r.nterms == 9 && r.nlv == 3) {
        res_check_differential("the Koszul syzygy matrix matches its "
                "reference entry for entry", &r, 3, rlen, rexp, rcomp, rcf);
    } else {
        RES_CHECK(0, "the Koszul syzygy matrix has nine terms");
    }
    res_free_resolution(&r);
}

/* The catalecticant {{x,y,z},{y,z,w}} has a single syzygy, the vector of
 * signed maximal minors; this is Macaulay2's syz A term for term. */
static void res_test_syz_of_input_catalecticant(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,0,0,  0,1,0,0,
        0,1,0,0,  0,0,1,0,
        0,0,1,0,  0,0,0,1
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
    const int32_t rd[2]    = {0, 0};
    const int32_t ranks[]  = {2, 3, 1, -1};
    const int32_t degs[]   = {0,0,  1,1,1,  3};
    const int32_t p        = RES_FC;

    const int32_t rlen[4]  = {2, 2, 2,  6};
    const int32_t rexp[]   = {
        1,0,0,0, 0,1,0,0,
        0,1,0,0, 0,0,1,0,
        0,0,1,0, 0,0,0,1,
        /* z^2 - yw, -yz + xw, y^2 - xz */
        0,0,2,0, 0,1,0,1,
        0,1,1,0, 1,0,0,1,
        0,2,0,0, 1,0,1,0
    };
    const int32_t rcomp[]  = {1,2, 1,2, 1,2,  1,1,2,2,3,3};
    const int32_t rcf[]    = {1,1, 1,1, 1,1,  1,p-1, p-1,1, 1,p-1};

    res_check_resolution("the catalecticant has one syzygy, in degree 3",
            lens, exps, comps, cfs, rd, 4, 2, 3, 2, RES_SYZ_OF_INPUT,
            ranks, degs);

    res_res_t r;
    res_run_resolution(&r, lens, exps, comps, cfs, rd,
            4, 2, 3, 2, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT));
    if (r.nterms == 12 && r.nlv == 3) {
        res_check_differential("the catalecticant syzygy is the vector of "
                "signed maximal minors", &r, 4, rlen, rexp, rcomp, rcf);
    } else {
        RES_CHECK(0, "the catalecticant syzygy matrix has twelve terms");
    }
    res_free_resolution(&r);
}

/* Where the two flavours part company.  For the twisted cubic the input
 * generators are a Gröbner basis, so the Schreyer syzygies of the basis
 * are the two Hilbert--Burch columns, while the graph module produces a
 * Gröbner basis of the same syzygy module and that has a third, redundant
 * element in degree 4.  Both generate what Macaulay2's syz generates --
 * res_reference.m2 checks the module equality -- but only the first is
 * minimal.  Minimalization is M5. */
static void res_test_syz_of_input_is_a_basis_not_a_minimal_one(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0,  0,2,0,0,
        1,0,0,1,  0,1,1,0,
        0,1,0,1,  0,0,2,0
    };
    const int32_t cfs[6]   = {1, -1, 1, -1, 1, -1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ranks[]  = {1, 3, 3, -1};
    const int32_t degs[]   = {0,  2,2,2,  3,3,4};

    res_check_resolution("the syzygies of the twisted cubic generators are "
            "a Gröbner basis of the syzygy module, so 3 rather than 2",
            lens, exps, comps, cfs, NULL, 4, 1, 3, 2, RES_SYZ_OF_INPUT,
            ranks, degs);
}

/* --- stopping conditions --------------------------------------------- *
 *
 *  res_stop_t, one field at a time.  The references here are msolve's own
 *  untruncated answers rather than tables, so the tests stay honest if the
 *  corpus examples move; what is hard coded is only the *relation* between
 *  a truncated run and a complete one, which is the contract.
 * --------------------------------------------------------------------- */

typedef struct res_mgb_t res_mgb_t;
struct res_mgb_t
{
    int32_t  bld;
    int32_t *blen;
    int32_t *bexp;
    int32_t *bcomp;
    int32_t *bcf;
    int64_t  nterms;
    int32_t *cfs;   /* the working copy msolve reduced in place */
};

static int64_t res_run_module_gb(
        res_mgb_t *g,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const res_grading_t * const grading,
        const res_stop_t * const stop
        )
{
    int32_t i, nt = 0;
    void *cf = NULL;

    memset(g, 0, sizeof(res_mgb_t));
    for (i = 0; i < ngens; ++i) {
        nt += lens[i];
    }
    g->cfs = (int32_t *)malloc((unsigned long)nt * sizeof(int32_t));
    memcpy(g->cfs, cfs_src, (unsigned long)nt * sizeof(int32_t));

    g->nterms = export_module_f4(malloc, &g->bld, &g->blen, &g->bexp,
            &g->bcomp, &cf, lens, exps, comps, g->cfs, row_degs,
            RES_FC, 0 /* drl */, res_strat_p(RES_MORD_POT), grading, stop,
            nv, nrows, ngens, 12 /* ht size */, 1 /* threads */,
            0 /* max pairs */, 2 /* la */, 1 /* reduce */, 0 /* info */);
    g->bcf = (int32_t *)cf;

    return g->nterms;
}

static void res_free_module_gb(
        res_mgb_t *g
        )
{
    void *cf = g->bcf;
    free_module_f4_result_data(free, &g->blen, &g->bexp, &g->bcomp, &cf);
    g->bcf = NULL;
    free(g->cfs);
    g->cfs = NULL;
}

/* Element k of a basis, as an offset into the flat term arrays. */
static int64_t res_gb_offset(
        const res_mgb_t * const g,
        const int32_t k
        )
{
    int64_t o = 0;
    int32_t i;

    for (i = 0; i < k; ++i) {
        o += g->blen[i];
    }
    return o;
}

/* Whether element a of one basis and element b of another are the same
 * element of the free module, term for term. */
static int res_gb_elements_agree(
        const res_mgb_t * const x,
        const int32_t a,
        const res_mgb_t * const y,
        const int32_t b,
        const int32_t nv
        )
{
    const int32_t len = x->blen[a];
    const int64_t oa  = res_gb_offset(x, a);
    const int64_t ob  = res_gb_offset(y, b);
    int64_t t;

    if (len != y->blen[b]) {
        return 0;
    }
    for (t = 0; t < len; ++t) {
        int32_t c;
        if (x->bcomp[oa+t] != y->bcomp[ob+t]
                || x->bcf[oa+t] != y->bcf[ob+t]) {
            return 0;
        }
        for (c = 0; c < nv; ++c) {
            if (x->bexp[(oa+t)*nv + c] != y->bexp[(ob+t)*nv + c]) {
                return 0;
            }
        }
    }
    return 1;
}

/* The heft degree of element k, read off its lead term: the input is
 * homogeneous, so any term would do. */
static int32_t res_gb_degree(
        const res_mgb_t * const g,
        const int32_t k,
        const int32_t nv,
        const int32_t * const vdegs, /* NULL means the standard grading */
        const int32_t * const heft
        )
{
    const int64_t o = res_gb_offset(g, k);
    int32_t c, d = 0;

    for (c = 0; c < nv; ++c) {
        const int32_t e = g->bexp[o*nv + c];
        if (vdegs == NULL) {
            d += e;
        } else {
            /* one weight per variable, folded through the heft */
            int32_t j, w = 0;
            for (j = 0; j < 2; ++j) {
                w += heft[j] * vdegs[c*2 + j];
            }
            d += e * w;
        }
    }
    return d;
}

/* The whole point of a ceiling: what a truncated run has in hand when it
 * stops is exactly the part of the complete basis it had reached.  F4
 * selects S-pairs by degree, and for homogeneous input a term of degree d
 * can only be reduced by a lead term of degree at most d, so the final
 * interreduction cannot tell the two runs apart either. */
static void res_check_degree_limit(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs,
        const int32_t nv,
        const int32_t ngens,
        const res_grading_t * const grading,
        const int32_t *vdegs, /* NULL, or grading->degs again for the test */
        const int32_t *heft
        )
{
    res_mgb_t full;
    int32_t i, k, d, dmin = INT32_MAX, dmax = -1;

    if (res_run_module_gb(&full, lens, exps, comps, cfs, NULL,
                nv, 1, ngens, grading, NULL) <= 0) {
        RES_CHECK(0, what);
        return;
    }

    for (k = 0; k < full.bld; ++k) {
        const int32_t dk = res_gb_degree(&full, k, nv, vdegs, heft);
        if (dk < dmin) { dmin = dk; }
        if (dk > dmax) { dmax = dk; }
    }
    RES_CHECK(dmax > dmin, what);

    for (d = dmin; d <= dmax; ++d) {
        /* the ceiling is a multidegree; the one used here is d times the
         * multidegree of the first variable, whose heft is d when the
         * grading is standard and the caller's own scale otherwise */
        int32_t md[RES_MTAB_MAXLEN];
        res_stop_t stop = res_stop_none();
        res_mgb_t part;
        int32_t nexp = 0;
        int ok = 1;

        memset(md, 0, sizeof(md));
        if (grading == NULL) {
            md[0] = d;
        } else {
            for (i = 0; i < 2; ++i) {
                md[i] = d * vdegs[i]; /* deg of x_1^d, so heft d * heft(x_1) */
            }
        }
        stop.max_degree = md;

        if (res_run_module_gb(&part, lens, exps, comps, cfs, NULL,
                    nv, 1, ngens, grading, &stop) <= 0) {
            RES_CHECK(0, "a run under a degree ceiling produced a basis");
            res_free_module_gb(&full);
            return;
        }

        for (k = 0; k < full.bld; ++k) {
            if (res_gb_degree(&full, k, nv, vdegs, heft) <= d) {
                nexp++;
            }
        }
        RES_CHECK(part.bld == nexp, "a degree ceiling keeps exactly the "
                "elements of the complete basis that are of no greater "
                "degree");

        for (k = 0; k < part.bld; ++k) {
            int found = 0;
            if (res_gb_degree(&part, k, nv, vdegs, heft) > d) {
                ok = 0;
            }
            for (i = 0; i < full.bld && !found; ++i) {
                found = res_gb_elements_agree(&part, k, &full, i, nv);
            }
            if (!found) {
                ok = 0;
            }
        }
        RES_CHECK(ok, "every element of a truncated basis is an element of "
                "the complete one, of degree within the ceiling");

        res_free_module_gb(&part);
    }

    res_free_module_gb(&full);
}

static void res_test_gb_degree_limit(
        void
        )
{
    /* x1^2 + x2*x3, x2^2 + x1*x3, x3^2 + x1*x2, whose basis reaches
     * degree 3 */
    const int32_t lens[3] = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,
        0,2,0,  1,0,1,
        0,0,2,  1,1,0
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};

    res_check_degree_limit("the basis of three quadrics spans more than one "
            "degree, so a ceiling has something to cut",
            lens, exps, comps, cfs, 3, 3, NULL, NULL, NULL);
}

/* The ceiling is a multidegree and is honoured as its heft, which is the
 * coarsening res_stop_t documents.  Here the heft is (1,2) rather than
 * (1,1), so a ceiling that is small in the first coordinate and one that
 * is small in the second are *not* interchangeable -- only their hefts
 * matter, and the test below pins that down by asking for two different
 * multidegrees of the same heft and demanding the same answer. */
static void res_test_gb_degree_limit_multigraded(
        void
        )
{
    /* R = k[x,y,z,w] with deg x = deg z = (1,0) and deg y = deg w = (0,1),
     * hefted by (1,2); I = (xy, zw, xw - zy) is bihomogeneous of
     * multidegree (1,1) and heft degree 3 */
    const int32_t vdegs[8] = {1,0,  0,1,  1,0,  0,1};
    const int32_t heft[2]  = {1, 2};
    const res_grading_t grading = {2, 0, NULL, vdegs, heft};

    const int32_t lens[3]  = {1, 1, 2};
    const int32_t exps[16] = {
        1,1,0,0,
        0,0,1,1,
        1,0,0,1,  0,1,1,0
    };
    const int32_t cfs[4]   = {1, 1, 1, RES_FC-1};
    const int32_t comps[4] = {1, 1, 1, 1};

    res_check_degree_limit("a hefted Z^2 grading truncates by heft degree",
            lens, exps, comps, cfs, 4, 3, &grading, vdegs, heft);

    /* Two ceilings of the same heft must truncate at the same place, and
     * the heft picked here is one below the top degree of the complete
     * basis, so both of them really do cut something. */
    res_mgb_t full;
    int32_t k, dmax = -1;
    if (res_run_module_gb(&full, lens, exps, comps, cfs, NULL,
                4, 1, 3, &grading, NULL) <= 0) {
        RES_CHECK(0, "the multigraded example has a basis to truncate");
        return;
    }
    for (k = 0; k < full.bld; ++k) {
        const int32_t dk = res_gb_degree(&full, k, 4, vdegs, heft);
        if (dk > dmax) { dmax = dk; }
    }

    /* the heft is (1,2), so (h,0) and (h-2,1) are different multidegrees
     * of the same heft h */
    const int32_t h  = dmax - 1;
    const int32_t a[2] = {h, 0};
    const int32_t b[2] = {h - 2, 1};
    res_stop_t sa = res_stop_none(), sb = res_stop_none();
    res_mgb_t ga, gb;
    sa.max_degree = a;
    sb.max_degree = b;

    const int64_t na = res_run_module_gb(&ga, lens, exps, comps, cfs, NULL,
            4, 1, 3, &grading, &sa);
    const int64_t nb = res_run_module_gb(&gb, lens, exps, comps, cfs, NULL,
            4, 1, 3, &grading, &sb);

    RES_CHECK(h >= 2 && na > 0 && ga.bld < full.bld,
            "a ceiling one heft degree below the top of the basis cuts "
            "something");
    RES_CHECK(na > 0 && na == nb && ga.bld == gb.bld,
            "two ceilings of the same heft truncate at the same place");
    if (na > 0 && na == nb && ga.bld == gb.bld) {
        int ok = 1;
        for (k = 0; k < ga.bld; ++k) {
            if (!res_gb_elements_agree(&ga, k, &gb, k, 4)) {
                ok = 0;
            }
        }
        RES_CHECK(ok, "and produce the same basis term for term");
    }

    res_free_module_gb(&ga);
    res_free_module_gb(&gb);
    res_free_module_gb(&full);
}

/* A ceiling at or below the degree shift of the ambient free module leaves
 * the round loop nothing it could select, and the only way to ask for one
 * is to have the scale wrong.  It is refused rather than answered. */
static void res_test_gb_degree_limit_rejects_bad_scale(
        void
        )
{
    const int32_t lens[2]  = {1, 1};
    const int32_t exps[6]  = {1,0,0,  0,1,0};
    const int32_t cfs[2]   = {1, 1};
    const int32_t comps[2] = {1, 1};
    const int32_t rd[1]    = {5}; /* R(-5), so the shift is 5 */

    const int32_t at[1]  = {5};
    const int32_t under[1] = {4};
    res_stop_t sa = res_stop_none(), su = res_stop_none();
    res_mgb_t ga, gu;
    sa.max_degree = at;
    su.max_degree = under;

    /* x + 1 is not homogeneous, so its schedule is by sugar degree */
    const int32_t ilens[2]  = {2, 1};
    const int32_t iexps[9]  = {1,0,0,  0,0,0,  0,1,0};
    const int32_t icfs[3]   = {1, 1, 1};
    const int32_t icomps[3] = {1, 1, 1};
    const int32_t two[1]    = {2};
    res_stop_t si = res_stop_none();
    res_mgb_t gi;
    si.max_degree = two;

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM,
                "res_selftest: three degree limit errors are expected next\n");
    }
    RES_CHECK(res_run_module_gb(&ga, lens, exps, comps, cfs, rd,
                3, 1, 2, NULL, &sa) == 0,
            "a ceiling at the degree shift is refused");
    RES_CHECK(res_run_module_gb(&gu, lens, exps, comps, cfs, rd,
                3, 1, 2, NULL, &su) == 0,
            "a ceiling below the degree shift is refused");
    RES_CHECK(res_run_module_gb(&gi, ilens, iexps, icomps, icfs, NULL,
                3, 1, 2, NULL, &si) == 0,
            "a ceiling on inhomogeneous input is refused, since the degree "
            "it would bound is a sugar degree");
    /* and the same input without a ceiling is computed as it always was */
    RES_CHECK(res_run_module_gb(&gi, ilens, iexps, icomps, icfs, NULL,
                3, 1, 2, NULL, NULL) > 0,
            "while the same inhomogeneous input without one still computes");

    res_free_module_gb(&ga);
    res_free_module_gb(&gu);
    res_free_module_gb(&gi);
}

/* SyzygyLimit.  The three Koszul relations of (x,y,z) all turn up in the
 * same round, so the round loop cannot stop after exactly one of them --
 * the cap on the export is what makes the count come out.  Both halves
 * matter and this test pins the visible half; that the round loop really
 * does stop early is what res_test_syz_limit_stops_early checks. */
static void res_test_syz_limit(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    int32_t lim;

    for (lim = 1; lim <= 5; ++lim) {
        const int32_t want = lim < 3 ? lim : 3;
        res_stop_t stop = res_stop_none();
        res_res_t r;
        stop.syz_limit = lim;

        res_run_resolution_stop(&r, lens, exps, comps, cfs, NULL,
                3, 1, 3, 2, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT),
                NULL, &stop);

        RES_CHECK(r.nlv == 3 && r.ranks[0] == 1 && r.ranks[1] == 3
                && r.ranks[2] == want,
                "a syzygy limit reports that many Koszul relations, or all "
                "of them if it asks for more than there are");

        res_free_resolution(&r);
    }
}

/* That the limit is a genuine early stop and not only a cap on the output.
 * The twisted cubic's syzygy module has a Gröbner basis of three elements,
 * two in degree 3 and a redundant one in degree 4; asking for two stops in
 * degree 3, so the degree 4 element is never computed, and the reported
 * degrees say so. */
static void res_test_syz_limit_stops_early(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0,  0,2,0,0,
        1,0,0,1,  0,1,1,0,
        0,1,0,1,  0,0,2,0
    };
    const int32_t cfs[6]   = {1, RES_FC-1, 1, RES_FC-1, 1, RES_FC-1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};

    res_stop_t stop = res_stop_none();
    res_res_t r;
    stop.syz_limit = 2;

    res_run_resolution_stop(&r, lens, exps, comps, cfs, NULL,
            4, 1, 3, 2, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT),
            NULL, &stop);

    RES_CHECK(r.nlv == 3 && r.ranks[2] == 2,
            "asking for two of the twisted cubic's three syzygies gives two");
    if (r.nlv == 3 && r.ranks[2] == 2) {
        /* degs is laid out level by level: 1 + 3 entries come first */
        RES_CHECK(r.degs[4] == 3 && r.degs[5] == 3,
                "and they are the two of degree 3, so the degree 4 element "
                "was never reached");
    }

    res_free_resolution(&r);
}

/* Bad input must be refused rather than trusted. */
static void res_test_resolution_rejects_bad_input(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};

    /* x + 1 is not homogeneous */
    const int32_t ilens[1]  = {2};
    const int32_t iexps[6]  = {1,0,0,  0,0,0};
    const int32_t icfs[2]   = {1, 1};
    const int32_t icomps[2] = {1, 1};

    const int32_t xlens[1]  = {1};
    const int32_t xexps[1]  = {65536};
    const int32_t xcfs[1]   = {1};
    const int32_t xcomps[1] = {1};
    const int32_t xrows[2]  = {INT32_MIN, INT32_MAX};

    res_res_t r;

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM, "res_selftest: six resolution input errors are "
                "expected next\n");
    }

    RES_CHECK(res_run_resolution(&r, lens, exps, comps, cfs, NULL,
                3, 1, 3, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_SCHREYER)) == 0,
            "the Schreyer order is refused as a base by the resolution");
    RES_CHECK(r.ranks == NULL && r.dlen == NULL && r.dcf == NULL,
            "a rejected resolution allocates nothing");
    res_free_resolution(&r);

    RES_CHECK(res_run_resolution(&r, lens, exps, comps, cfs, NULL,
                3, 1, 3, 3, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT)) == 0,
            "syzygies of the input do not resolve past level 2");
    res_free_resolution(&r);

    RES_CHECK(res_run_resolution(&r, ilens, iexps, icomps, icfs, NULL,
                3, 1, 1, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT)) == 0,
            "an inhomogeneous ideal has no graded resolution");
    res_free_resolution(&r);

    RES_CHECK(res_run_resolution(&r, ilens, iexps, icomps, icfs, NULL,
                3, 1, 1, 2, RES_SYZ_OF_INPUT, res_strat_p(RES_MORD_POT)) == 0,
            "an inhomogeneous ideal has no graded syzygy matrix either");
    res_free_resolution(&r);

    RES_CHECK(res_run_resolution(&r, xlens, xexps, xcomps, xcfs, NULL,
                1, 1, 1, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT)) == 0,
            "an exponent outside the hash table range is refused");
    res_free_resolution(&r);

    RES_CHECK(res_run_resolution(&r, xlens, xcomps, xcomps, xcfs, xrows,
                1, 2, 1, 0, RES_SYZ_OF_GB, res_strat_p(RES_MORD_POT)) == 0,
            "row shifts outside the hash table range are refused");
    res_free_resolution(&r);
}

/* --------------------------------------------------------------------- *
 *  Minimal Betti numbers and Hilbert information
 *
 *  Every reference table below came from Macaulay2's
 *
 *      minimalBetti M, poincare M, pdim M, regularity M, dim M, degree M
 *
 *  over the same prime; see res_reference.m2.  The tables are written as
 *  (level, degree, value) triples exactly as the frame ones are, so the
 *  two can be read side by side -- the interesting examples are the ones
 *  where they differ.
 * --------------------------------------------------------------------- */

static void res_check_betti(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t *ref,
        const int32_t rpdim,
        const int32_t rreg,
        const int32_t rdim,
        const int64_t rdeg
        )
{
    int32_t i, l, d;

    int32_t nterms = 0;
    for (i = 0; i < ngens; ++i) {
        nterms += lens[i];
    }
    /* msolve reduces the coefficient array in place, and this input is
     * used twice */
    int32_t *cfs = (int32_t *)malloc((unsigned long)nterms * sizeof(int32_t));

    int32_t nlv = 0, maxdeg = 0, shift = 0;
    int32_t pdim = -2, reg = -2, dim = -2;
    int64_t deg  = -1;
    int32_t *tab = NULL, *num = NULL;

    memcpy(cfs, cfs_src, (unsigned long)nterms * sizeof(int32_t));
    const int64_t nelts = export_module_betti(malloc, &nlv, &maxdeg, &shift,
            &tab, &num, &pdim, &reg, &dim, &deg, NULL,
            lens, exps, comps, cfs, row_degs, 32003, 0 /* drl */,
            res_strat_p(RES_MORD_POT), NULL, nv, nrows, ngens, 0 /* full */, 1 /* minimal */,
            1 /* verify */, 12, 1, 0, 2, 0);

    RES_CHECK(nelts > 0 && tab != NULL && num != NULL, what);
    if (nelts <= 0 || tab == NULL || num == NULL) {
        free(cfs);
        free_module_betti_result_data(free, &tab, &num);
        return;
    }

    /* every reference entry is present, and nothing else is */
    int ok = 1;
    int64_t seen = 0, want = 0;
    for (i = 0; ref[i] >= 0; i += 3) {
        want += ref[i+2];
        if (ref[i] >= nlv || ref[i+1] > maxdeg
                || tab[ref[i]*(maxdeg+1) + ref[i+1]] != ref[i+2]) {
            ok = 0;
        }
    }
    for (l = 0; l < nlv; ++l) {
        for (d = 0; d <= maxdeg; ++d) {
            seen += tab[l*(maxdeg+1) + d];
        }
    }
    RES_CHECK(ok && seen == want,
            "the minimal Betti table matches its Macaulay2 reference");

    RES_CHECK(pdim == rpdim, "the projective dimension matches Macaulay2");
    RES_CHECK(reg == rreg, "the regularity matches Macaulay2");
    RES_CHECK(dim == rdim, "the Krull dimension matches Macaulay2");
    RES_CHECK(deg == rdeg, "the degree matches Macaulay2");

    /* The numerator the engine reports is the alternating sum of the
     * *frame* ranks, computed before any field arithmetic happens; this
     * recomputes it from the minimal table, which is a different set of
     * numbers whenever the frame is nonminimal.  The two agreeing is the
     * telescoping identity of res.h, and it fails the moment a rank
     * correction lands at the wrong level or degree. */
    int hok = 1;
    for (d = 0; d <= maxdeg; ++d) {
        int32_t s = 0;
        for (l = 0; l < nlv; ++l) {
            const int32_t c = tab[l*(maxdeg+1) + d];
            s = (l & 1) ? s - c : s + c;
        }
        if (s != num[d]) {
            hok = 0;
        }
    }
    RES_CHECK(hok, "the Hilbert numerator is the alternating sum of the "
            "minimal Betti numbers, not just of the frame ranks");

    /* the same run without the differential must give the frame ranks,
     * which dominate the minimal ones entry for entry and have the same
     * alternating sum */
    int32_t fnlv = 0, fmaxdeg = 0;
    int32_t *ftab = NULL, *fnum = NULL;

    memcpy(cfs, cfs_src, (unsigned long)nterms * sizeof(int32_t));
    const int64_t fnelts = export_module_betti(malloc, &fnlv, &fmaxdeg, NULL,
            &ftab, &fnum, NULL, NULL, NULL, NULL, NULL,
            lens, exps, comps, cfs, row_degs, 32003, 0, res_strat_p(RES_MORD_POT),
            NULL, nv, nrows, ngens, 0, 0 /* frame only */, 0, 12, 1, 0, 2, 0);

    RES_CHECK(fnelts == nelts && fnlv == nlv && fmaxdeg == maxdeg
            && ftab != NULL && fnum != NULL,
            "the frame only run agrees on the shape of the table");
    if (ftab != NULL && fnum != NULL && fnlv == nlv && fmaxdeg == maxdeg) {
        int dom = 1, same = 1;
        int64_t ftot = 0;
        for (l = 0; l < nlv; ++l) {
            for (d = 0; d <= maxdeg; ++d) {
                const int32_t fv = ftab[l*(maxdeg+1) + d];
                ftot += fv;
                if (fv < tab[l*(maxdeg+1) + d]) {
                    dom = 0;
                }
            }
        }
        for (d = 0; d <= maxdeg; ++d) {
            if (fnum[d] != num[d]) {
                same = 0;
            }
        }
        RES_CHECK(dom, "the frame ranks dominate the minimal Betti numbers");
        RES_CHECK(ftot == nelts,
                "the frame ranks add up to the size of the frame");
        RES_CHECK(same, "both runs report the same Hilbert numerator");
    }
    free_module_betti_result_data(free, &ftab, &fnum);

    free(cfs);
    free_module_betti_result_data(free, &tab, &num);
}

static void res_test_betti_koszul(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,1,3,  2,2,3,  3,3,1,  -1
    };

    res_check_betti("the Koszul complex is already minimal, 1,3,3,1",
            lens, exps, comps, cfs, NULL, 3, 1, 3, ref, 3, 0, 0, 1);
}

/* The case the whole milestone is for: the frame is 1,6,8,3 and the
 * minimal resolution is 1,3,3,1, so eleven generators have to be
 * cancelled by rank corrections.  Everything that could be off by a level
 * or a degree shows up here. */
static void res_test_betti_nonminimal(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[18] = {
        2,0,0,  0,1,1,   /* x^2 + yz */
        0,2,0,  1,0,1,   /* y^2 + xz */
        0,0,2,  1,1,0    /* z^2 + xy */
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,2,3,  2,4,3,  3,6,1,  -1
    };

    res_check_betti("three quadrics minimalize from the frame 1,6,8,3 to "
            "the Koszul shape 1,3,3,1",
            lens, exps, comps, cfs, NULL, 3, 1, 3, ref, 3, 3, 0, 8);
}

static void res_test_betti_twisted_cubic(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,1,0,  0,2,0,0,
        1,0,0,1,  0,1,1,0,
        0,1,0,1,  0,0,2,0
    };
    const int32_t cfs[6]   = {1, -1, 1, -1, 1, -1};
    const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,2,3,  2,3,2,  -1
    };

    res_check_betti("the twisted cubic is the Hilbert-Burch resolution, of "
            "dimension two and degree three",
            lens, exps, comps, cfs, NULL, 4, 1, 3, ref, 2, 1, 2, 3);
}

static void res_test_betti_block_order(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {2,0,0,  1,1,0,  0,3,0};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};
    const int32_t ref[]    = {
        0,0,1,  1,2,2,  1,3,1,  2,3,1,  2,4,1,  -1
    };

    res_check_betti("(x^2, xy, y^3) is minimal already, and z leaves it "
            "one dimensional",
            lens, exps, comps, cfs, NULL, 3, 1, 3, ref, 2, 2, 1, 4);
}

/* A rank two module: level 0 carries two generators and the minimal
 * resolution 2,3,1 sits inside the frame 2,6,5,1, so a whole level
 * disappears. */
static void res_test_betti_module(
        void
        )
{
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,0,0,  0,1,0,0,
        0,1,0,0,  0,0,1,0,
        0,0,1,0,  0,0,0,1
    };
    const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
    const int32_t rd[2]    = {0, 0};
    const int32_t ref[]    = {
        0,0,2,  1,1,3,  2,3,1,  -1
    };

    res_check_betti("the catalecticant cokernel minimalizes from 2,6,5,1 "
            "to 2,3,1",
            lens, exps, comps, cfs, rd, 4, 2, 3, ref, 2, 1, 2, 3);
}

/* Row degrees, and a global shift on top of them: only the differences
 * matter, so the table is the same and degshift carries the offset. */
static void res_test_betti_row_degrees(
        void
        )
{
    const int32_t lens[2]  = {2, 2};
    const int32_t exps[16] = {
        2,0,0,0,  0,0,1,0,
        0,2,0,0,  0,0,0,1
    };
    int32_t cfs[4]         = {1, 1, 1, 1};
    const int32_t comps[4] = {1, 2, 1, 2};
    const int32_t rd[2]    = {0, 1};
    const int32_t sd[2]    = {5, 6};
    const int32_t ref[]    = {
        0,0,1,  0,1,1,  1,2,2,  -1
    };

    res_check_betti("coker {{x2,y2},{z,w}} has the minimal table 2,2 with "
            "the second generator in degree one",
            lens, exps, comps, cfs, rd, 4, 2, 2, ref, 1, 1, 3, 3);

    int32_t nlv = 0, maxdeg = 0, shift = 0, pdim = -2, reg = -2, dim = -2;
    int64_t deg  = -1;
    int32_t *tab = NULL;

    int32_t c2[4] = {1, 1, 1, 1};
    const int64_t n = export_module_betti(malloc, &nlv, &maxdeg, &shift,
            &tab, NULL, &pdim, &reg, &dim, &deg, NULL,
            lens, exps, comps, c2, sd, 32003, 0, res_strat_p(RES_MORD_POT),
            NULL, 4, 2, 2, 0, 1, 1, 12, 1, 0, 2, 0);

    RES_CHECK(n > 0 && tab != NULL && shift == 5,
            "a global shift of the row degrees is reported, not baked in");
    if (tab != NULL) {
        int ok = nlv >= 2 && maxdeg >= 2
            && tab[0*(maxdeg+1) + 0] == 1
            && tab[0*(maxdeg+1) + 1] == 1
            && tab[1*(maxdeg+1) + 2] == 2;
        RES_CHECK(ok, "shifting every row degree by five leaves the table "
                "where it was");
    }
    /* The table is indexed by degree and so has to stay shifted, but the
     * scalars are reported in the caller's own degrees: only regularity
     * sees the difference. */
    RES_CHECK(reg == 1 + 5, "the regularity comes back on the caller's own "
            "degree scale, not the shifted one");
    RES_CHECK(pdim == 1 && dim == 3 && deg == 3,
            "projective dimension, dimension and degree do not move with "
            "the shift");
    free_module_betti_result_data(free, &tab, NULL);
}

/* Something with a table wide enough to exercise the elimination: the
 * cube of the maximal ideal in three variables, 1,10,15,6. */
static void res_test_betti_monomial(
        void
        )
{
    const int32_t lens[10]   = {1,1,1,1,1,1,1,1,1,1};
    const int32_t exps[30]   = {
        3,0,0,  2,1,0,  2,0,1,  1,2,0,  1,1,1,
        1,0,2,  0,3,0,  0,2,1,  0,1,2,  0,0,3
    };
    const int32_t cfs[10]    = {1,1,1,1,1,1,1,1,1,1};
    const int32_t comps[10]  = {1,1,1,1,1,1,1,1,1,1};
    const int32_t ref[]      = {
        0,0,1,  1,3,10,  2,4,15,  3,5,6,  -1
    };

    res_check_betti("the cube of the maximal ideal in three variables has "
            "the Eagon-Northcott table 1,10,15,6",
            lens, exps, comps, cfs, NULL, 3, 1, 10, ref, 3, 2, 0, 10);
}

/* Six random cubics in four variables.  Every example above has scalar
 * blocks whose pivots happen to be plus or minus one, so the elimination
 * gets the right rank even without normalizing its pivot rows; here the
 * coefficients are generic and it does not.  This is the only check in the
 * file that fails when res_ech_rank stops dividing by the pivot. */
static void res_test_betti_generic_cubics(
        void
        )
{
    const int32_t lens[6]   = {3, 3, 3, 3, 3, 3};
    const int32_t exps[72]  = {
        0,1,2,0,  1,1,0,1,  0,2,0,1,
        1,0,0,2,  1,1,1,0,  1,2,0,0,
        1,0,1,1,  0,1,1,1,  1,0,0,2,
        2,1,0,0,  1,0,2,0,  2,0,0,1,
        1,0,2,0,  1,1,1,0,  2,0,0,1,
        0,0,1,2,  1,1,1,0,  2,0,1,0
    };
    const int32_t cfs[18]   = {
         3277, 10825, 23704,
        22284, 19561, 23260,
        19176, 20404, 27057,
         3298,  4024,  8445,
         1838, 10295,  1669,
        12895, 14444, 10377
    };
    const int32_t comps[18] = {
        1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1, 1,1,1
    };
    const int32_t ref[]     = {
        0,0,1,
        1,3,6,
        2,5,12,  2,6,2,
        3,6,6,   3,7,6,
        4,8,3,
        -1
    };

    res_check_betti("six generic cubics minimalize from the frame "
            "1,18,39,30,8 to 1,6,14,12,3",
            lens, exps, comps, cfs, NULL, 4, 1, 6, ref, 4, 4, 1, 8);
}

/* Truncation keeps the head of the table, but Hilbert information needs
 * the whole alternating sum and has to be refused. */
static void res_test_betti_truncation(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    int32_t cfs[3]         = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};

    int32_t nlv = 0, maxdeg = 0;
    int32_t *tab = NULL, *num = NULL;

    int32_t c1[3] = {1, 1, 1};
    const int64_t n = export_module_betti(malloc, &nlv, &maxdeg, NULL,
            &tab, NULL, NULL, NULL, NULL, NULL, NULL,
            lens, exps, comps, c1, NULL, 32003, 0, res_strat_p(RES_MORD_POT),
            NULL, 3, 1, 3, 2 /* truncate */, 1, 1, 12, 1, 0, 2, 0);

    /* The frame is built to level three so that beta_2 knows about d_3,
     * hence eight elements, but only levels zero to two are reported. */
    RES_CHECK(n == 8 && nlv == 3 && tab != NULL,
            "truncating the Koszul table at level two keeps 1,3,3");
    if (tab != NULL) {
        RES_CHECK(tab[0*(maxdeg+1)+0] == 1 && tab[1*(maxdeg+1)+1] == 3
                && tab[2*(maxdeg+1)+2] == 3,
                "the truncated table is the head of the full one");
    }
    free_module_betti_result_data(free, &tab, &num);

    if (res_st_verbose > 0) {
        fprintf(VERBSTREAM, "res_selftest: three Betti input errors are "
                "expected next\n");
    }
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL,
                &tab, &num, NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, 32003, 0, res_strat_p(RES_MORD_POT),
                NULL, 3, 1, 3, 2, 1, 1, 12, 1, 0, 2, 0) == 0,
            "a truncated resolution has no Hilbert numerator");
    RES_CHECK(tab == NULL && num == NULL,
            "a rejected Betti call allocates nothing");
    RES_CHECK(export_module_betti(malloc, &nlv, &maxdeg, NULL,
                &tab, NULL, NULL, NULL, NULL, NULL, NULL,
                lens, exps, comps, cfs, NULL, 32003, 0,
                res_strat_p(RES_MORD_SCHREYER),
                NULL, 3, 1, 3, 0, 1, 1, 12, 1, 0, 2, 0) == 0,
            "the Schreyer order is refused as a base by the Betti table");
}

/* --------------------------------------------------------------------- *
 *  A resolution kept alive
 *
 *  The handle has to be indistinguishable from the one shot entry point:
 *  same ranks, same degrees in the same storage order, and the same
 *  differential term for term.  What makes it worth its own tests is the
 *  laziness -- asking for level i computes levels 2 to i, so asking for
 *  the levels in descending order makes the driver fill in a prefix it
 *  was never asked for, and asking twice must change nothing.
 * --------------------------------------------------------------------- */

static void res_check_comp(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens,
        const int32_t max_level,
        const int descending
        )
{
    int32_t i, k, lev;
    int64_t t;
    res_res_t r;

    const int64_t n = res_run_resolution(&r, lens, exps, comps, cfs_src,
            row_degs, nv, nrows, ngens, max_level, RES_SYZ_OF_GB,
            res_strat_p(RES_MORD_POT));
    if (n <= 0) {
        RES_CHECK(0, what);
        res_free_resolution(&r);
        return;
    }

    int32_t nt = 0;
    for (i = 0; i < ngens; ++i) {
        nt += lens[i];
    }
    int32_t *cfs = (int32_t *)malloc((unsigned long)nt * sizeof(int32_t));
    memcpy(cfs, cfs_src, (unsigned long)nt * sizeof(int32_t));

    res_comp_t *c = res_comp_new(lens, exps, comps, cfs, row_degs, RES_FC,
            0 /* drl */, res_strat_p(RES_MORD_POT), NULL, nv, nrows, ngens, max_level,
            12 /* ht size */, 1 /* threads */, 0 /* max pairs */,
            2 /* la */, 0 /* info */);

    RES_CHECK(c != NULL && res_comp_nlevels(c) == r.nlv, what);
    if (c == NULL || res_comp_nlevels(c) != r.nlv) {
        res_comp_free(&c);
        free(cfs);
        res_free_resolution(&r);
        return;
    }

    /* the shape of the resolution, with no differential computed at all */
    int ok = 1, okd = 1;
    int32_t off = 0;
    for (i = 0; i < r.nlv; ++i) {
        if (res_comp_rank(c, i) != r.ranks[i]) {
            ok = 0;
            break;
        }
        int32_t *dg = (int32_t *)malloc(
                (unsigned long)(r.ranks[i] > 0 ? r.ranks[i] : 1)
                * sizeof(int32_t));
        if (res_comp_degrees(c, i, dg)) {
            okd = 0;
        }
        for (k = 0; k < r.ranks[i]; ++k) {
            if (dg[k] != r.degs[off+k]) {
                okd = 0;
            }
        }
        off += r.ranks[i];
        free(dg);
    }
    RES_CHECK(ok, "the free modules of the handle have the frame's ranks");
    RES_CHECK(okd, "the free modules of the handle have the frame's degrees, "
            "in the frame's own storage order");
    RES_CHECK(res_comp_rank(c, -1) == -1 && res_comp_rank(c, r.nlv) == -1,
            "a level outside the resolution has no free module");

    /* offsets of each level in the one shot entry point's flat arrays */
    int64_t *coff = (int64_t *)calloc((unsigned long)r.nlv + 1,
            sizeof(int64_t));
    int64_t *toff = (int64_t *)calloc((unsigned long)r.nlv + 1,
            sizeof(int64_t));
    for (i = 1; i < r.nlv; ++i) {
        coff[i+1] = coff[i] + r.ranks[i];
        toff[i+1] = toff[i];
        for (k = 0; k < r.ranks[i]; ++k) {
            toff[i+1] += r.dlen[coff[i] + k];
        }
    }

    int okt = 1;
    for (i = 1; i < r.nlv; ++i) {
        lev = descending ? r.nlv - i : i;

        int32_t *dlen = NULL, *dexp = NULL, *dcomp = NULL;
        void *dcf = NULL;
        const int64_t nterms = res_comp_differential(malloc, c, lev,
                &dlen, &dexp, &dcomp, &dcf);
        const int32_t *cf = (const int32_t *)dcf;

        if (nterms != toff[lev+1] - toff[lev] || dlen == NULL
                || dexp == NULL || dcomp == NULL || cf == NULL) {
            okt = 0;
        } else {
            for (k = 0; k < r.ranks[lev]; ++k) {
                if (dlen[k] != r.dlen[coff[lev] + k]) {
                    okt = 0;
                }
            }
            for (t = 0; t < nterms; ++t) {
                if (dcomp[t] != r.dcomp[toff[lev] + t]
                        || cf[t] != r.dcf[toff[lev] + t]) {
                    okt = 0;
                }
            }
            for (t = 0; t < nterms * nv; ++t) {
                if (dexp[t] != r.dexp[toff[lev] * nv + t]) {
                    okt = 0;
                }
            }
        }
        free_module_differential_data(free, &dlen, &dexp, &dcomp, &dcf);
        RES_CHECK(dlen == NULL && dexp == NULL && dcomp == NULL
                && dcf == NULL, "releasing a differential clears the "
                "caller's pointers");
    }
    RES_CHECK(okt, descending
            ? "every differential matches the one shot entry point term for "
              "term, asked for from the top down"
            : "every differential matches the one shot entry point term for "
              "term");

    /* asking again for one already computed changes nothing */
    if (r.nlv > 1) {
        int32_t *dlen = NULL, *dexp = NULL, *dcomp = NULL;
        void *dcf = NULL;
        const int64_t again = res_comp_differential(malloc, c, 1,
                &dlen, &dexp, &dcomp, &dcf);
        int oka = again == toff[2] - toff[1];
        for (t = 0; t < again && oka; ++t) {
            if (dcomp[t] != r.dcomp[t]
                    || ((const int32_t *)dcf)[t] != r.dcf[t]) {
                oka = 0;
            }
        }
        RES_CHECK(oka, "asking twice for the same differential gives the "
                "same answer");
        free_module_differential_data(free, &dlen, &dexp, &dcomp, &dcf);
    }

    RES_CHECK(res_comp_differential(malloc, c, 0, NULL, NULL, NULL, NULL)
            == 0, "there is no differential at level zero");

    free(coff);
    free(toff);
    res_comp_free(&c);
    RES_CHECK(c == NULL, "releasing the handle clears the caller's pointer");
    free(cfs);
    res_free_resolution(&r);
}

static void res_test_comp_koszul(
        void
        )
{
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t cfs[3]   = {1, 1, 1};
    const int32_t comps[3] = {1, 1, 1};

    res_check_comp("the handle resolves the Koszul complex",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, 0);
    res_check_comp("the Koszul complex, from the top down",
            lens, exps, comps, cfs, NULL, 3, 1, 3, 0, 1);
}

static void res_test_comp_twisted_cubic(
        void
        )
{
    /* y^2-xz, yz-xw, z^2-yw in x,y,z,w */
    const int32_t lens[3]   = {2, 2, 2};
    const int32_t exps[24]  = {
        0,2,0,0,  1,0,1,0,
        0,1,1,0,  1,0,0,1,
        0,0,2,0,  0,1,0,1};
    const int32_t cfs[6]    = {1, RES_FC-1, 1, RES_FC-1, 1, RES_FC-1};
    const int32_t comps[6]  = {1, 1, 1, 1, 1, 1};

    res_check_comp("the handle resolves the twisted cubic",
            lens, exps, comps, cfs, NULL, 4, 1, 3, 0, 0);
    res_check_comp("the twisted cubic, from the top down",
            lens, exps, comps, cfs, NULL, 4, 1, 3, 0, 1);
}

/* A frame strictly larger than the minimal resolution, so the levels
 * really do have work to do, and one that runs past nv. */
static void res_test_comp_nonminimal(
        void
        )
{
    /* z, y^2, x^2 y, x^3 in x,y,z: the frame is 1,4,6,4,1 */
    const int32_t lens[4]   = {1, 1, 1, 1};
    const int32_t exps[12]  = {0,0,1,  0,2,0,  2,1,0,  3,0,0};
    const int32_t cfs[4]    = {1, 1, 1, 1};
    const int32_t comps[4]  = {1, 1, 1, 1};

    res_check_comp("the handle resolves a frame that runs past nv",
            lens, exps, comps, cfs, NULL, 3, 1, 4, 0, 0);
    res_check_comp("a frame that runs past nv, from the top down",
            lens, exps, comps, cfs, NULL, 3, 1, 4, 0, 1);
}

static void res_test_comp_module(
        void
        )
{
    /* coker {{x2,y2},{z,w}}: a rank two ambient free module that is
     * homogeneous only because the second row sits in degree one, so
     * components and row degrees both have to reach the handle */
    const int32_t lens[2]   = {2, 2};
    const int32_t exps[16]  = {
        2,0,0,0,  0,0,1,0,
        0,2,0,0,  0,0,0,1};
    const int32_t cfs[4]    = {1, 1, 1, 1};
    const int32_t comps[4]  = {1, 2, 1, 2};
    const int32_t rd[2]     = {0, 1};

    res_check_comp("the handle resolves a shifted rank two module",
            lens, exps, comps, cfs, rd, 4, 2, 2, 0, 0);
    res_check_comp("a shifted rank two module, from the top down",
            lens, exps, comps, cfs, rd, 4, 2, 2, 0, 1);
}

static void res_test_comp_degshift_and_truncation(
        void
        )
{
    /* the catalecticant cokernel, whose frame is 2,6,5,1, with every row
     * degree pushed up by three so that the normalization has something
     * to report */
    const int32_t lens[3]  = {2, 2, 2};
    const int32_t exps[24] = {
        1,0,0,0,  0,1,0,0,
        0,1,0,0,  0,0,1,0,
        0,0,1,0,  0,0,0,1};
    int32_t cfs[6]         = {1, 1, 1, 1, 1, 1};
    const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
    const int32_t rd[2]    = {3, 3};

    res_comp_t *c = res_comp_new(lens, exps, comps, cfs, rd, RES_FC, 0,
            res_strat_p(RES_MORD_POT), NULL, 4, 2, 3, 0, 12, 1, 0, 2, 0);
    RES_CHECK(c != NULL, "the handle takes shifted row degrees");
    if (c == NULL) {
        return;
    }
    RES_CHECK(res_comp_degshift(c) == 3,
            "the handle reports the shift it normalized the row degrees by");
    RES_CHECK(res_comp_is_complete(c),
            "a frame with no ceiling on it ends on its own");
    const int32_t full = res_comp_nlevels(c);
    RES_CHECK(full == 4 && res_comp_rank(c, 1) == 6,
            "the handle's frame is the catalecticant's 2,6,5,1");
    res_comp_free(&c);

    /* the same input, cut off */
    int32_t cfs2[6] = {1, 1, 1, 1, 1, 1};
    c = res_comp_new(lens, exps, comps, cfs2, rd, RES_FC, 0, res_strat_p(RES_MORD_POT), NULL,
            4, 2, 3, 2, 12, 1, 0, 2, 0);
    RES_CHECK(c != NULL && res_comp_nlevels(c) == 3,
            "max_level truncates the handle's frame");
    RES_CHECK(c != NULL && full > 3 && !res_comp_is_complete(c),
            "a cut off frame does not report itself as complete");
    if (c != NULL) {
        int32_t *dlen = NULL, *dexp = NULL, *dcomp = NULL;
        void *dcf = NULL;
        RES_CHECK(res_comp_differential(malloc, c, 2, &dlen, &dexp,
                    &dcomp, &dcf) > 0,
                "the top level of a truncated handle still differentiates");
        free_module_differential_data(free, &dlen, &dexp, &dcomp, &dcf);
        RES_CHECK(res_comp_differential(malloc, c, 3, &dlen, &dexp,
                    &dcomp, &dcf) == 0,
                "past the truncation there is no differential to ask for");
        RES_CHECK(dlen == NULL && dcf == NULL,
                "a refused differential allocates nothing");
    }
    res_comp_free(&c);
}

static void res_test_comp_rejects_bad_input(
        void
        )
{
    int32_t cfs[3]         = {1, 1, 1};
    const int32_t lens[3]  = {1, 1, 1};
    const int32_t exps[9]  = {1,0,0,  0,1,0,  0,0,1};
    const int32_t comps[3] = {1, 1, 1};

    RES_CHECK(res_comp_new(lens, exps, comps, cfs, NULL, RES_FC, 0,
                res_strat_p(RES_MORD_SCHREYER), NULL,
                3, 1, 3, 0, 12, 1, 0, 2, 0) == NULL,
            "the handle refuses the Schreyer order as a base");
    res_strat_t bad = res_strat_default();
    bad.pos = 42;
    RES_CHECK(res_comp_new(lens, exps, comps, cfs, NULL, RES_FC, 0,
                &bad, NULL, 3, 1, 3, 0, 12, 1, 0, 2, 0) == NULL,
            "the handle refuses an unknown component direction");
    bad = res_strat_default();
    bad.lift = 42;
    RES_CHECK(res_comp_new(lens, exps, comps, cfs, NULL, RES_FC, 0,
                &bad, NULL, 3, 1, 3, 0, 12, 1, 0, 2, 0) == NULL,
            "the handle refuses an unknown lift");
    RES_CHECK(res_comp_new(lens, exps, comps, cfs, NULL, RES_FC, 0,
                res_strat_p(RES_MORD_POT), NULL, 3, 1, 3, -1, 12, 1, 0, 2, 0) == NULL,
            "the handle refuses a negative truncation level");

    /* x + y^2 is not homogeneous, so it has no graded resolution */
    int32_t icfs[2]         = {1, 1};
    const int32_t ilens[1]  = {2};
    const int32_t iexps[6]  = {1,0,0,  0,2,0};
    const int32_t icomps[2] = {1, 1};
    RES_CHECK(res_comp_new(ilens, iexps, icomps, icfs, NULL, RES_FC, 0,
                res_strat_p(RES_MORD_POT), NULL, 3, 1, 1, 0, 12, 1, 0, 2, 0) == NULL,
            "the handle refuses inhomogeneous input");

    res_comp_t *nc = NULL;
    res_comp_free(&nc);
    RES_CHECK(res_comp_nlevels(NULL) == 0 && res_comp_rank(NULL, 0) == -1
            && res_comp_degrees(NULL, 0, NULL) != 0,
            "every query of a null handle is answered rather than crashing");
}

/* --------------------------------------------------------------------- *
 *  Strategies
 *
 *  The order the resolution runs in is a parameter, and the four bases
 *  and directions produce genuinely different Gröbner bases, frames and
 *  differentials -- the whole point of making it a parameter.  What they
 *  must *not* change is anything that is an invariant of the module, and
 *  that is what these tests pin down:
 *
 *    - the minimal Betti numbers, entry for entry;
 *    - the Hilbert numerator, which is the alternating sum of the frame
 *      ranks under every strategy even though the ranks themselves are
 *      not equal;
 *    - projective dimension, regularity, Krull dimension and degree.
 *
 *  Anything wrong with a comparator shows up here immediately, because a
 *  broken order gives a Gröbner basis that is not one and the ranks stop
 *  telescoping.  d o d = 0 is checked exactly for every strategy too, on
 *  the resolution rather than the table.
 * --------------------------------------------------------------------- */

static const res_strat_t res_strat_matrix[4] = {
    {RES_MORD_POT, RES_POS_DOWN, RES_LIFT_SCHREYER},
    {RES_MORD_POT, RES_POS_UP,   RES_LIFT_SCHREYER},
    {RES_MORD_TOP, RES_POS_DOWN, RES_LIFT_SCHREYER},
    {RES_MORD_TOP, RES_POS_UP,   RES_LIFT_SCHREYER}
};

static void res_check_strategies(
        const char *what,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const int32_t *cfs_src,
        const int32_t *row_degs,
        const int32_t nv,
        const int32_t nrows,
        const int32_t ngens
        )
{
    int32_t k, i, d, nt = 0;
    char msg[160];

    for (i = 0; i < ngens; ++i) {
        nt += lens[i];
    }

    int32_t rnlv = 0, rmax = 0, rshift = 0, rpdim = 0, rreg = 0, rdim = 0;
    int64_t rdeg = 0;
    int32_t *rbetti = NULL, *rhilb = NULL;

    for (k = 0; k < 4; ++k) {
        const res_strat_t * const sp = res_strat_matrix + k;
        const char * const nm = res_strat_name(sp);

        int32_t *cfs = (int32_t *)malloc((unsigned long)nt * sizeof(int32_t));
        memcpy(cfs, cfs_src, (unsigned long)nt * sizeof(int32_t));

        int32_t nlv = 0, maxdeg = 0, shift = 0, pdim = 0, reg = 0, dim = 0;
        int64_t deg = 0;
        int32_t *betti = NULL, *hilb = NULL;
        const int64_t n = export_module_betti(malloc, &nlv, &maxdeg, &shift,
                &betti, &hilb, &pdim, &reg, &dim, &deg, NULL,
                lens, exps, comps, cfs, row_degs, RES_FC, 0, sp,
                NULL, nv, nrows, ngens, 0 /* no ceiling */, 1 /* minimal */,
                0 /* the exact d o d = 0 check runs below instead */,
                12, 1, 0, 2, 0);
        free(cfs);

        snprintf(msg, sizeof(msg), "%s resolves under %s", what, nm);
        RES_CHECK(n > 0 && betti != NULL && hilb != NULL, msg);
        if (n <= 0 || betti == NULL || hilb == NULL) {
            free_module_betti_result_data(free, &betti, &hilb);
            continue;
        }

        if (k == 0) {
            rnlv = nlv; rmax = maxdeg; rshift = shift;
            rpdim = pdim; rreg = reg; rdim = dim; rdeg = deg;
            rbetti = betti; rhilb = hilb;
            continue;
        }

        /* the tables may run to different lengths, the frames being
         * different; every entry they do not share has to be zero */
        int okb = shift == rshift;
        const int32_t md = maxdeg > rmax ? maxdeg : rmax;
        const int32_t ml = nlv > rnlv ? nlv : rnlv;
        for (i = 0; i < ml && okb; ++i) {
            for (d = 0; d <= md; ++d) {
                const int32_t a = (i < rnlv && d <= rmax)
                    ? rbetti[(size_t)i * (rmax + 1) + d] : 0;
                const int32_t b = (i < nlv && d <= maxdeg)
                    ? betti[(size_t)i * (maxdeg + 1) + d] : 0;
                if (a != b) {
                    okb = 0;
                    break;
                }
            }
        }
        snprintf(msg, sizeof(msg),
                "%s has the same minimal Betti numbers under %s", what, nm);
        RES_CHECK(okb, msg);

        int okh = 1;
        for (d = 0; d <= md; ++d) {
            const int32_t a = d <= rmax ? rhilb[d] : 0;
            const int32_t b = d <= maxdeg ? hilb[d] : 0;
            if (a != b) {
                okh = 0;
            }
        }
        snprintf(msg, sizeof(msg),
                "%s has the same Hilbert numerator under %s", what, nm);
        RES_CHECK(okh, msg);

        snprintf(msg, sizeof(msg),
                "%s has the same pdim, regularity, dimension and degree "
                "under %s", what, nm);
        RES_CHECK(pdim == rpdim && reg == rreg && dim == rdim && deg == rdeg,
                msg);

        free_module_betti_result_data(free, &betti, &hilb);
    }
    free_module_betti_result_data(free, &rbetti, &rhilb);

    /* and the differential itself, checked exactly, under each strategy */
    for (k = 0; k < 4; ++k) {
        res_res_t r;
        const int64_t n = res_run_resolution(&r, lens, exps, comps, cfs_src,
                row_degs, nv, nrows, ngens, 0, RES_SYZ_OF_GB,
                res_strat_matrix + k);
        snprintf(msg, sizeof(msg), "%s is a complex under %s", what,
                res_strat_name(res_strat_matrix + k));
        const uint32_t pt[8] = {2, 3, 5, 7, 11, 13, 17, 19};
        RES_CHECK(n > 0 && res_composite_is_zero(&r, nv, pt), msg);
        res_free_resolution(&r);
    }
}

static void res_test_strategies(
        void
        )
{
    /* an ideal: every strategy has to agree here for a trivial reason,
     * there being one component, which makes this the control */
    {
        const int32_t lens[3]  = {2, 2, 2};
        const int32_t exps[24] = {
            0,2,0,0,  1,0,1,0,
            0,1,1,0,  1,0,0,1,
            0,0,2,0,  0,1,0,1};
        const int32_t cfs[6]   = {1, RES_FC-1, 1, RES_FC-1, 1, RES_FC-1};
        const int32_t comps[6] = {1, 1, 1, 1, 1, 1};
        res_check_strategies("the twisted cubic",
                lens, exps, comps, cfs, NULL, 4, 1, 3);
    }

    /* the catalecticant cokernel, rank two and no degree shift */
    {
        const int32_t lens[3]  = {2, 2, 2};
        const int32_t exps[24] = {
            1,0,0,0,  0,1,0,0,
            0,1,0,0,  0,0,1,0,
            0,0,1,0,  0,0,0,1};
        const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
        const int32_t comps[6] = {1, 2, 1, 2, 1, 2};
        const int32_t rd[2]    = {0, 0};
        res_check_strategies("the catalecticant cokernel",
                lens, exps, comps, cfs, rd, 4, 2, 3);
    }

    /* a shifted ambient free module, where the two bases really diverge:
     * the component's degree shift is in the degree the Gröbner basis was
     * compared by, so term over position has to add it back and position
     * over term must not */
    {
        const int32_t lens[2]  = {2, 2};
        const int32_t exps[16] = {
            2,0,0,0,  0,0,1,0,
            0,2,0,0,  0,0,0,1};
        const int32_t cfs[4]   = {1, 1, 1, 1};
        const int32_t comps[4] = {1, 2, 1, 2};
        const int32_t rd[2]    = {0, 1};
        res_check_strategies("coker {{x2,y2},{z,w}} with a shifted row",
                lens, exps, comps, cfs, rd, 4, 2, 2);
    }

    /* rank three, so the direction of the component key has three values
     * to permute rather than two */
    {
        const int32_t lens[3]  = {2, 2, 2};
        const int32_t exps[18] = {
            1,0,0,  0,1,0,
            0,1,0,  0,0,1,
            1,0,0,  0,0,1};
        const int32_t cfs[6]   = {1, 1, 1, 1, 1, 1};
        const int32_t comps[6] = {1, 2, 2, 3, 1, 3};
        const int32_t rd[3]    = {0, 0, 0};
        res_check_strategies("a rank three module",
                lens, exps, comps, cfs, rd, 3, 3, 3);
    }
}

static void res_test_strategy_names(
        void
        )
{
    res_strat_t s = res_strat_default();

    RES_CHECK(s.base == RES_MORD_POT && s.pos == RES_POS_DOWN
            && s.lift == RES_LIFT_SCHREYER,
            "the default strategy is position over term, component down, "
            "Schreyer above level zero");
    RES_CHECK(strcmp(res_strat_name(NULL), "pot-down-schreyer") == 0,
            "a null strategy names itself as the default");
    RES_CHECK(strcmp(res_strat_name(&s), "pot-down-schreyer") == 0,
            "the default strategy has the name the default has");
    s.base = RES_MORD_TOP;
    s.pos  = RES_POS_UP;
    RES_CHECK(strcmp(res_strat_name(&s), "top-up-schreyer") == 0,
            "each axis appears in the name");
    RES_CHECK(res_strat_check(&s, 1) == 0,
            "term over position with the component up is a usable strategy");
    RES_CHECK(res_strat_check(NULL, 1) == 0,
            "a null strategy is usable, being the default");
    s.base = RES_MORD_SCHREYER;
    RES_CHECK(res_strat_check(&s, 1) != 0,
            "the Schreyer order is not usable as a *base* order");
    RES_CHECK(strcmp(res_strat_name(&s), "unknown") == 0,
            "an unusable strategy has no name");

    const res_strat_t t = res_strat_of_order(RES_MORD_TOP);
    RES_CHECK(t.base == RES_MORD_TOP && t.pos == RES_POS_DOWN
            && t.lift == RES_LIFT_SCHREYER,
            "a bare module order becomes that base with the other axes "
            "left at their defaults");
}

/* --------------------------------------------------------------------- *
 *  Entry point
 * --------------------------------------------------------------------- */

int main(void)
{
    const int verbose = 1;

    res_st_fail    = 0;
    res_st_run     = 0;
    res_st_verbose = verbose;

    res_test_standard_grading();
    res_test_multigrading();
    res_test_torsion_grading();
    res_test_grading_rejects_bad_heft();
    res_test_module_gb_rank_one(RES_MORD_POT,
            "a rank one module basis equals the ideal basis under POT");
    res_test_module_gb_rank_one(RES_MORD_TOP,
            "a rank one module basis equals the ideal basis under TOP");
    res_test_module_gb_split_components();
    res_test_module_gb_rank_two();
    res_test_classical_resolutions();
    res_test_module_gb_rejects_bad_input();
    res_test_frame_koszul();
    res_test_frame_twisted_cubic();
    res_test_frame_nonminimal();
    res_test_frame_block_order();
    res_test_frame_past_nv();
    res_test_frame_module();
    res_test_frame_row_degrees();
    res_test_frame_truncation();
    res_test_frame_rejects_bad_input();
    res_test_frame_rejects_exponent_overflow();
    res_test_resolution_koszul();
    res_test_resolution_twisted_cubic();
    res_test_resolution_nonminimal();
    res_test_resolution_module();
    res_test_resolution_row_degrees();
    res_test_resolution_truncation();
    res_test_resolution_matches_frame();
    res_test_syz_of_input_koszul();
    res_test_syz_of_input_catalecticant();
    res_test_syz_of_input_is_a_basis_not_a_minimal_one();
    res_test_gb_degree_limit();
    res_test_gb_degree_limit_multigraded();
    res_test_gb_degree_limit_rejects_bad_scale();
    res_test_syz_limit();
    res_test_syz_limit_stops_early();
    res_test_resolution_rejects_bad_input();
    res_test_betti_koszul();
    res_test_betti_nonminimal();
    res_test_betti_twisted_cubic();
    res_test_betti_block_order();
    res_test_betti_module();
    res_test_betti_row_degrees();
    res_test_betti_monomial();
    res_test_betti_generic_cubics();
    res_test_betti_truncation();
    res_test_comp_koszul();
    res_test_comp_twisted_cubic();
    res_test_comp_nonminimal();
    res_test_comp_module();
    res_test_comp_degshift_and_truncation();
    res_test_comp_rejects_bad_input();
    res_test_degree_buckets();
    res_test_weighted_betti();
    res_test_weighted_order();
    res_test_multigraded_betti();
    res_test_multigraded_module();
    res_test_torsion_betti();
    res_test_explicit_standard_grading();
    res_test_grading_rejects_bad_input();
    res_test_comp_multidegrees();
    res_test_strategy_names();
    res_test_strategies();

    if (verbose > 0) {
        fprintf(VERBSTREAM, "res_selftest: %d checks, %d failures\n",
                res_st_run, res_st_fail);
    }
    return res_st_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
