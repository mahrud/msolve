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
            fc, 0 /* drl */, module_order,
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
            fc, 0, RES_MORD_POT,
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
            fc, 0, RES_MORD_POT,
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
            fc, 0, RES_MORD_POT, 4, 1, 3, 12, 1, 0, 2, 1, 0);

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
            fc, 0, RES_MORD_POT, 4, 3, 2, 12, 1, 0, 2, 1, 0);

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
            fc, 0, RES_MORD_POT, 4, 2, 3, 12, 1, 0, 2, 1, 0);

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
                lens, exps, bad_comps, cfs, NULL, 32003, 0, RES_MORD_POT,
                3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "an out of range component is rejected");
    RES_CHECK(blen == NULL && bexp == NULL && bcomp == NULL && bcf == NULL,
            "a rejected call allocates nothing");

    RES_CHECK(export_module_f4(malloc, &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, ok_comps, cfs, NULL, 32003, 0,
                RES_MORD_SCHREYER, 3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "the Schreyer order is refused at this entry point");

    RES_CHECK(export_module_f4(malloc, &bld, &blen, &bexp, &bcomp, &bcf,
                lens, exps, ok_comps, cfs, NULL, 0 /* char 0 */, 0,
                RES_MORD_POT, 3, 2, 1, 12, 1, 0, 2, 1, 0) == 0,
            "characteristic zero is refused");
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

    if (verbose > 0) {
        fprintf(VERBSTREAM, "res_selftest: %d checks, %d failures\n",
                res_st_run, res_st_fail);
    }
    return res_st_fail == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
