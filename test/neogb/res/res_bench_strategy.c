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

/* Which order strategy should a resolution run in?
 *
 * The choice is a res_strat_t parameter precisely because it is an
 * empirical question, and this is what measures it: the same random
 * corpus resolved under every strategy, reporting the size of the frame,
 * the number of terms in the nonminimal differential, and wall time.
 *
 * The frame size is the number that matters most.  It is the rank of
 * every free module in the nonminimal resolution, so it sets the shape of
 * every Macaulay matrix above level one, and the linear algebra is where
 * the time goes on anything large.
 *
 * Correctness is not assumed: the alternating sum of the ranks is the
 * Euler characteristic of the module and cannot depend on the strategy,
 * so a disagreement means a comparator is wrong rather than merely
 * different.  res_selftest.c makes the far stronger check that the whole
 * minimal Betti table agrees; this one only needs a cheap net.
 *
 * Not a test: it is in noinst_PROGRAMS rather than check_PROGRAMS, since
 * what it reports is a measurement and not a pass or a fail.
 *
 *   ./res_bench_strategy [trials] [seed]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../../../src/neogb/res.h"
#include "../../../src/msolve/streams.h"

#define BENCH_FC 32003
#define BENCH_NSTRAT 4

static const res_strat_t bench_strats[BENCH_NSTRAT] = {
    {RES_MORD_POT, RES_POS_DOWN, RES_LIFT_SCHREYER},
    {RES_MORD_POT, RES_POS_UP,   RES_LIFT_SCHREYER},
    {RES_MORD_TOP, RES_POS_DOWN, RES_LIFT_SCHREYER},
    {RES_MORD_TOP, RES_POS_UP,   RES_LIFT_SCHREYER}
};

static double bench_now(
        void
        )
{
    struct timespec t;

    clock_gettime(CLOCK_MONOTONIC, &t);

    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static uint64_t bench_state = 88172645463325252ull;

static uint32_t bench_rnd(
        const uint32_t n
        )
{
    bench_state ^= bench_state << 13;
    bench_state ^= bench_state >> 7;
    bench_state ^= bench_state << 17;

    return (uint32_t)(bench_state % n);
}

/* One random homogeneous submodule of R^nr, with the generators of R^nr
 * in degrees sh and every column in degree deg.  A column's terms in
 * component c carry ring degree deg - sh[c], which is what makes it
 * homogeneous; msolve's flat input wants the module monomials of one
 * column distinct, so duplicates are dropped rather than merged. */
typedef struct bench_in_t bench_in_t;
struct bench_in_t
{
    int32_t *lens, *exps, *comps, *cfs, *rd;
    int32_t  nv, nr, ng, nt;
};

static void bench_free(
        bench_in_t *in
        )
{
    free(in->lens);
    free(in->exps);
    free(in->comps);
    free(in->cfs);
    free(in->rd);
    memset(in, 0, sizeof(bench_in_t));
}

static void bench_gen(
        bench_in_t *in,
        const int32_t nv,
        const int32_t nr,
        const int32_t ng,
        const int32_t maxshift,
        const int32_t extradeg
        )
{
    int32_t i, j, k, q, p, mx = 0;

    memset(in, 0, sizeof(bench_in_t));
    in->nv = nv;
    in->nr = nr;
    in->ng = ng;

    in->rd = (int32_t *)calloc((unsigned long)nr, sizeof(int32_t));
    for (i = 0; i < nr; ++i) {
        in->rd[i] = maxshift > 0 ? (int32_t)bench_rnd((uint32_t)maxshift) : 0;
        if (in->rd[i] > mx) {
            mx = in->rd[i];
        }
    }
    const int32_t deg = mx + 2 + extradeg;
    const size_t cap  = (size_t)ng * (size_t)nr * 3;

    in->lens  = (int32_t *)calloc((unsigned long)ng, sizeof(int32_t));
    in->exps  = (int32_t *)calloc(cap * (size_t)nv, sizeof(int32_t));
    in->comps = (int32_t *)calloc(cap, sizeof(int32_t));
    in->cfs   = (int32_t *)calloc(cap, sizeof(int32_t));

    for (j = 0; j < ng; ++j) {
        int32_t added = 0;
        for (i = 0; i < nr; ++i) {
            const int32_t d = deg - in->rd[i];
            if (d < 0) {
                continue;
            }
            const int32_t howmany = 1 + (int32_t)bench_rnd(2);
            for (k = 0; k < howmany; ++k) {
                int32_t e[16];
                for (q = 0; q < nv; ++q) {
                    e[q] = 0;
                }
                for (q = 0; q < d; ++q) {
                    e[bench_rnd((uint32_t)nv)]++;
                }
                int dup = 0;
                for (p = in->nt - added; p < in->nt && !dup; ++p) {
                    if (in->comps[p] != i + 1) {
                        continue;
                    }
                    dup = 1;
                    for (q = 0; q < nv; ++q) {
                        if (in->exps[(size_t)p * nv + q] != e[q]) {
                            dup = 0;
                            break;
                        }
                    }
                }
                if (dup) {
                    continue;
                }
                for (q = 0; q < nv; ++q) {
                    in->exps[(size_t)in->nt * nv + q] = e[q];
                }
                in->comps[in->nt] = i + 1;
                in->cfs[in->nt]   = 1 + (int32_t)bench_rnd(BENCH_FC - 1);
                in->nt++;
                added++;
            }
        }
        in->lens[j] = added;
    }
}

typedef struct bench_out_t bench_out_t;
struct bench_out_t
{
    int32_t nlev;
    int64_t frame;   /* generators of every free module, level 0 included */
    int64_t terms;   /* terms of the nonminimal differential              */
    int64_t euler;   /* sum_i (-1)^i rank F_i, an invariant of the module  */
    double  secs;
    int     ok;
};

static void bench_run(
        bench_out_t *o,
        const bench_in_t * const in,
        const res_strat_t * const s
        )
{
    int32_t i;

    memset(o, 0, sizeof(bench_out_t));

    int32_t *cfs = (int32_t *)malloc((unsigned long)in->nt * sizeof(int32_t));
    memcpy(cfs, in->cfs, (unsigned long)in->nt * sizeof(int32_t));

    int32_t nlv = 0, *ranks = NULL, *degs = NULL, *dlen = NULL, *dexp = NULL;
    int32_t *dcomp = NULL;
    void *dcf = NULL;

    const double t0 = bench_now();
    const int64_t nterms = export_module_resolution(malloc, &nlv, &ranks,
            &degs, &dlen, &dexp, &dcomp, &dcf,
            in->lens, in->exps, in->comps, cfs, in->rd, BENCH_FC,
            0 /* drl */, s, in->nv, in->nr, in->ng, 0 /* no ceiling */,
            RES_SYZ_OF_GB, 0 /* the cheap structural check always runs */,
            17, 1, 0, 2, 0);
    o->secs = bench_now() - t0;
    free(cfs);

    if (nterms > 0 && ranks != NULL) {
        o->ok    = 1;
        o->nlev  = nlv;
        o->terms = nterms;
        for (i = 0; i < nlv; ++i) {
            o->frame += ranks[i];
            o->euler += (i % 2 ? -1 : 1) * (int64_t)ranks[i];
        }
    }
    free_module_resolution_result_data(free, &ranks, &degs, &dlen, &dexp,
            &dcomp, &dcf);
}

int main(
        int argc,
        char **argv
        )
{
    int k, t;

    const int trials = argc > 1 ? atoi(argv[1]) : 40;
    if (argc > 2) {
        bench_state = (uint64_t)strtoull(argv[2], NULL, 10) | 1ull;
    }

    int64_t frame[BENCH_NSTRAT], terms[BENCH_NSTRAT];
    int     wins[BENCH_NSTRAT];
    double  secs[BENCH_NSTRAT];
    for (k = 0; k < BENCH_NSTRAT; ++k) {
        frame[k] = 0;
        terms[k] = 0;
        wins[k]  = 0;
        secs[k]  = 0;
    }
    int done = 0, disagree = 0;

    printf("# %d trials, resolving each under every strategy\n", trials);
    printf("# %-5s %-3s %-3s %-3s", "case", "nv", "nr", "ng");
    for (k = 0; k < BENCH_NSTRAT; ++k) {
        printf("  %-22s", res_strat_name(bench_strats + k));
    }
    printf("\n");

    for (t = 0; t < trials; ++t) {
        bench_in_t in;
        bench_gen(&in, 3 + (int32_t)bench_rnd(2), 2 + (int32_t)bench_rnd(3),
                3 + (int32_t)bench_rnd(4), 4, (int32_t)bench_rnd(3));

        bench_out_t o[BENCH_NSTRAT];
        int allok = 1;
        for (k = 0; k < BENCH_NSTRAT; ++k) {
            bench_run(o + k, &in, bench_strats + k);
            if (!o[k].ok) {
                allok = 0;
            }
        }
        if (allok) {
            for (k = 1; k < BENCH_NSTRAT; ++k) {
                if (o[k].euler != o[0].euler) {
                    disagree++;
                }
            }
            int best = 0;
            for (k = 1; k < BENCH_NSTRAT; ++k) {
                if (o[k].frame < o[best].frame) {
                    best = k;
                }
            }
            wins[best]++;
            done++;
            printf("  %-5d %-3d %-3d %-3d", t, in.nv, in.nr, in.ng);
            for (k = 0; k < BENCH_NSTRAT; ++k) {
                frame[k] += o[k].frame;
                terms[k] += o[k].terms;
                secs[k]  += o[k].secs;
                printf("  %2d/%5ld/%7ld/%.3f", o[k].nlev, (long)o[k].frame,
                        (long)o[k].terms, o[k].secs);
            }
            printf("\n");
        }
        bench_free(&in);
    }

    printf("\n# %d resolved, %d Euler characteristic disagreements "
            "(anything but zero is a bug, not a difference)\n",
            done, disagree);
    printf("# %-22s %10s %12s %10s %8s %8s\n", "strategy", "frame",
            "diff terms", "seconds", "vs best", "smallest");
    int64_t best = 0;
    for (k = 0; k < BENCH_NSTRAT; ++k) {
        if (best == 0 || frame[k] < best) {
            best = frame[k];
        }
    }
    for (k = 0; k < BENCH_NSTRAT; ++k) {
        printf("  %-22s %10ld %12ld %10.3f %8.2f %8d\n",
                res_strat_name(bench_strats + k), (long)frame[k],
                (long)terms[k], secs[k],
                best > 0 ? (double)frame[k] / (double)best : 0.0, wins[k]);
    }

    return 0;
}
