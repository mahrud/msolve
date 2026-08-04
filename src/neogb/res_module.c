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

/* Gröbner bases of submodules of a free module.
 *
 * The entry point here mirrors export_f4 in f4.c, but takes a presentation
 * matrix rather than a list of polynomials: the columns of the matrix
 * generate a submodule of the free module R^nr_rows, and every term
 * carries the row it sits in.  See res.h for the flat array convention. */

#include "res.h"
#include "../msolve/streams.h"

/* Like set_exponent_vector in io.c, but for a module monomial: writes the
 * component into the trailing slot and folds that component's degree
 * shift into ev[DEG], which is the invariant the rest of the hash table
 * relies on (see data.h).  Block elimination orders are not supported
 * together with modules, so there is only one degree slot to fill. */
static inline void set_module_exponent_vector(
        exp_t *ev,
        const int32_t *iev,
        const int32_t *icomp,
        const int32_t idx,
        const ht_t *ht
        )
{
    len_t i;

    const len_t nv    = ht->nv;
    const exp_t comp  = (exp_t)icomp[idx];

    ev[DEG] = 0;
    for (i = 0; i < nv; ++i) {
        ev[i+1]  = (exp_t)(iev+(nv*idx))[i];
        ev[DEG]  = (exp_t)(ev[DEG] + ev[i+1]);
    }
    ev[DEG]      = (exp_t)(ev[DEG] + ht->cshift[comp]);
    ev[ht->cpos] = comp;
}

/* The module counterpart of import_input_data in io.c.  Only the monomial
 * import differs; coefficients go through the very same helper. */
static void import_module_input_data(
        bs_t *bs,
        md_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *vcfs,
        const int *invalid_gens
        )
{
    int32_t i, j;
    len_t k;
    hm_t *hm;
    len_t ctr = 0;

    ht_t *ht = bs->ht;

    int32_t off       = 0;
    int32_t init_off  = 0;

    len_t ngens = stop - start;

    for (i = 0; i < start; ++i) {
        init_off += lens[i];
    }

    check_enlarge_basis(bs, ngens, st);

    /* import monomials */
    exp_t *e = ht->ev[0]; /* scratch, as in import_input_data */
    off = init_off;
    for (i = start; i < stop; ++i) {
        if (invalid_gens == NULL || invalid_gens[i] == 0) {
            while (lens[i] >= ht->esz-ht->eld) {
                enlarge_hash_table(ht);
                e = ht->ev[0];
            }
            hm = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
            bs->hm[ctr] = hm;

            hm[COEFFS]  = ctr;
            hm[PRELOOP] = (lens[i] % UNROLL);
            hm[LENGTH]  = lens[i];

            bs->red[ctr] = 0;

            for (j = off; j < off+lens[i]; ++j) {
                set_module_exponent_vector(e, exps, comps, j, ht);
                hm[j-off+OFFSET] = insert_in_hash_table(e, ht);
            }
            ctr++;
        }
        off += lens[i];
    }

    /* import coefficients and sort each generator into the module order */
    ctr = import_input_coefficients(
            bs, st, start, stop, lens, vcfs, invalid_gens, init_off);

    ngens = ctr;

    /* The degree of a module element is that of its lead term, which
     * already includes the component shift.  There is no block order to
     * worry about here, so this is the nev == 0 branch of
     * import_input_data. */
    for (i = 0; i < ngens; ++i) {
        hm = bs->hm[i];
        bs->hm[i][DEG] = ht->hd[hm[OFFSET]].deg;
    }

    /* homogeneous means every term of every generator has the same
     * (shifted) degree; this is what makes the degree by degree strategy
     * in F4 exact rather than heuristic */
    st->homogeneous = 1;
    for (i = 0; i < ngens; ++i) {
        hm = bs->hm[i];
        const deg_t deg = ht->hd[hm[OFFSET]].deg;
        k = hm[LENGTH] + OFFSET;
        for (j = OFFSET+1; j < (int32_t)k; ++j) {
            if (deg != ht->hd[hm[j]].deg) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
done:

    st->ngens      = ngens;
    st->ngens_input = ngens;
    bs->ld         = ngens;
}

/* Exports the basis in the same flat layout as export_data in io.c, with
 * one extra array giving the component of each term. */
static int64_t export_module_data(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        int32_t **bcomp,
        void **bcf,
        void *(*mallocp) (size_t),
        const bs_t * const bs,
        const ht_t * const ht,
        const md_t * const md
        )
{
    len_t i, j, k;
    hm_t *dt;

    const len_t nv   = ht->nv;
    const len_t lml  = bs->lml;

    int64_t nterms = 0;
    int64_t nelts  = 0;

    for (i = 0; i < lml; ++i) {
        if (bs->hm[bs->lmps[i]] != NULL) {
            nterms += (int64_t)bs->hm[bs->lmps[i]][LENGTH];
        } else {
            nterms++; /* one term for the zero element */
        }
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        fprintf(ERRSTREAM,
                "Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)(*mallocp)(
            (unsigned long)nelts * sizeof(int32_t));
    int32_t *exp  = (int32_t *)(*mallocp)(
            (unsigned long)nterms * (unsigned long)nv * sizeof(int32_t));
    int32_t *comp = (int32_t *)(*mallocp)(
            (unsigned long)nterms * sizeof(int32_t));
    int32_t *cf   = (int32_t *)(*mallocp)(
            (unsigned long)nterms * sizeof(int32_t));

    int64_t cl = 0, ce = 0, cc = 0;

    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        if (bs->hm[bi] == NULL) {
            cf[cc]   = 0;
            comp[cc] = 0;
            for (k = 0; k < nv; ++k) {
                exp[ce++] = 0;
            }
            cc += 1;
            len[cl++] = 1;
            continue;
        }
        len[cl] = bs->hm[bi][LENGTH];
        /* msolve picks the coefficient width from the characteristic in
         * set_ff_bits, so a small prime lands in cf_8 or cf_16 and only a
         * large one in cf_32; the exported array is int32_t either way */
        switch (md->ff_bits) {
            case 8:
                for (j = 0; j < len[cl]; ++j) {
                    cf[cc+j] = (int32_t)bs->cf_8[bs->hm[bi][COEFFS]][j];
                }
                break;
            case 16:
                for (j = 0; j < len[cl]; ++j) {
                    cf[cc+j] = (int32_t)bs->cf_16[bs->hm[bi][COEFFS]][j];
                }
                break;
            case 32:
                for (j = 0; j < len[cl]; ++j) {
                    cf[cc+j] = (int32_t)bs->cf_32[bs->hm[bi][COEFFS]][j];
                }
                break;
            default:
                fprintf(ERRSTREAM, "Unsupported coefficient width %d in a "
                        "module basis.\n", md->ff_bits);
                return 0;
        }
        dt = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 1; k <= nv; ++k) {
                exp[ce++] = (int32_t)ht->ev[dt[j]][k];
            }
            comp[cc+j] = (int32_t)ht->ev[dt[j]][ht->cpos];
        }
        cc += len[cl];
        cl++;
    }

    *bload = (int32_t)nelts;
    *blen  = len;
    *bexp  = exp;
    *bcomp = comp;
    *bcf   = (void *)cf;

    return nterms;
}

/* Everything the module entry points have in common: validate the input,
 * build the meta data, import the presentation matrix and run F4.  On
 * success the caller owns the returned basis, *stp, and the shared hash
 * table data of the basis's hash table, and releases them with
 *
 *     free_shared_hash_data(gb->ht); free_basis(&gb); free(st);
 *
 * On failure NULL is returned, *stp is NULL, and everything this function
 * allocated has been released.
 *
 * module_order is restricted to POT and TOP: the Schreyer order needs per
 * component base monomials, which only the frame in res_frame.c can
 * supply, and by then the Gröbner basis is already in hand. */
static bs_t *module_gb_from_input(
        md_t **stp,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t info_level
        )
{
    int32_t i, j;

    *stp = NULL;

    if (field_char == 0 || field_char >= ((uint32_t)1 << 31)) {
        fprintf(ERRSTREAM, "Module Groebner bases need a prime field of "
                "characteristic less than 2^31.\n");
        return NULL;
    }
    if (nr_rows < 1 || nr_gens < 1 || nr_vars < 1) {
        fprintf(ERRSTREAM, "Empty module input.\n");
        return NULL;
    }
    if (module_order != RES_MORD_POT && module_order != RES_MORD_TOP) {
        fprintf(ERRSTREAM, "Only position over term and term over position "
                "are available here; the Schreyer order needs per component "
                "base monomials that only the resolution engine can supply.\n");
        return NULL;
    }
    if (mon_order != 0) {
        fprintf(ERRSTREAM, "Module Groebner bases are only implemented for "
                "the degree reverse lexicographic order.\n");
        return NULL;
    }

    /* Validate all dimensions before anything is allocated.  The importer
     * below uses int32_t offsets, while the hash table stores both individual
     * exponents and shifted total degrees in exp_t.  Reject values that would
     * otherwise be silently narrowed to a different monomial. */
    int64_t nt = 0;
    for (i = 0; i < nr_gens; ++i) {
        if (lens[i] < 1) {
            fprintf(ERRSTREAM, "Generator %d has no terms.\n", i);
            return NULL;
        }
        nt += lens[i];
        if (nt > INT32_MAX) {
            fprintf(ERRSTREAM, "Module input has too many terms.\n");
            return NULL;
        }
    }
    int32_t mn = 0;
    if (row_degs != NULL) {
        mn = row_degs[0];
        for (i = 1; i < nr_rows; ++i) {
            if (row_degs[i] < mn) {
                mn = row_degs[i];
            }
        }
        for (i = 0; i < nr_rows; ++i) {
            const int64_t shift = (int64_t)row_degs[i] - (int64_t)mn;
            if (shift > UINT16_MAX) {
                fprintf(ERRSTREAM, "Row degree %d is too far from the "
                        "minimum row degree to fit in the exponent table.\n",
                        i);
                return NULL;
            }
        }
    }
    for (j = 0; j < nt; ++j) {
        if (comps[j] < 1 || comps[j] > nr_rows) {
            fprintf(ERRSTREAM, "Term %d has component %d, outside the range "
                    "1 to %d.\n", j, comps[j], nr_rows);
            return NULL;
        }
        int64_t deg = row_degs != NULL
            ? (int64_t)row_degs[comps[j]-1] - (int64_t)mn : 0;
        for (i = 0; i < nr_vars; ++i) {
            const int32_t exponent = exps[(int64_t)j * nr_vars + i];
            if (exponent < 0 || exponent > UINT16_MAX) {
                fprintf(ERRSTREAM, "Term %d has exponent %d, outside the "
                        "16-bit exponent range.\n", j, exponent);
                return NULL;
            }
            deg += exponent;
        }
        if (deg > UINT16_MAX) {
            fprintf(ERRSTREAM, "Term %d has shifted total degree %ld, outside "
                    "the 16-bit exponent range.\n", j, (long)deg);
            return NULL;
        }
    }

    md_t *st = allocate_meta_data();
    if (st == NULL) {
        return NULL;
    }

    int32_t elim_block_len = 0, nr_nf = 0, reset_ht = 0, use_signatures = 0;
    int32_t pbm_file = 0, truncate_lifting = 0;
    int32_t l_mon_order = mon_order, l_nr_vars = nr_vars, l_nr_gens = nr_gens;
    int32_t l_ht_size = ht_size, l_nr_threads = nr_threads;
    int32_t l_max_nr_pairs = max_nr_pairs, l_la_option = la_option;
    int32_t l_reduce_gb = reduce_gb, l_info_level = info_level;
    uint32_t l_field_char = field_char;

    int *invalid_gens = NULL;
    int res = validate_input_data(&invalid_gens, cfs, lens, &l_field_char,
            &l_mon_order, &elim_block_len, &l_nr_vars, &l_nr_gens, &nr_nf,
            &l_ht_size, &l_nr_threads, &l_max_nr_pairs, &reset_ht,
            &l_la_option, &use_signatures, &l_reduce_gb, &truncate_lifting,
            &l_info_level);
    if (res == -1) {
        free(invalid_gens);
        free(st);
        return NULL;
    }
    if (check_and_set_meta_data(st, lens, exps, cfs, invalid_gens,
                l_field_char, l_mon_order, elim_block_len, l_nr_vars,
                l_nr_gens, nr_nf, l_ht_size, l_nr_threads, l_max_nr_pairs,
                reset_ht, l_la_option, use_signatures, l_reduce_gb, pbm_file,
                truncate_lifting, l_info_level)) {
        free(invalid_gens);
        free(st);
        return NULL;
    }

    /* this is what makes the hash table a module one; it has to happen
     * before initialize_basis, which is where the table is built */
    st->ncomp = nr_rows;
    st->mord  = module_order;

    bs_t *bs  = initialize_basis(st, NULL);
    ht_t *bht = bs->ht;

    /* Degree shifts of the ambient free module.  msolve keeps degrees in a
     * uint16_t, so they are normalized to start at zero; a global shift of
     * all row degrees changes neither the module nor its Groebner basis,
     * only the absolute degrees, which the caller can put back. */
    if (row_degs != NULL) {
        for (i = 0; i < nr_rows; ++i) {
            bht->cshift[i+1] = (deg_t)((int64_t)row_degs[i] - (int64_t)mn);
        }
    }

    import_module_input_data(bs, st, 0, st->ngens_input,
            lens, exps, comps, cfs, invalid_gens);

    print_initial_statistics(VERBSTREAM, st);

    calculate_divmask(bht);

    sort_r(bs->hm, (unsigned long)bs->ld, sizeof(hm_t *),
            initial_input_cmp, bht);
    normalize_initial_basis(bs, st->fc);

    int32_t err = 0;
    bs_t *gb = core_gba(bs, st, &err, (len_t)field_char);

    free(invalid_gens);

    if (gb == NULL || err > 0) {
        fprintf(ERRSTREAM, "Module Groebner basis computation failed.\n");
        if (gb != NULL) {
            free_shared_hash_data(gb->ht);
            free_basis(&gb);
        }
        free(st);
        return NULL;
    }

    *stp = st;

    return gb;
}

/* Frames and resolutions are graded objects: every degree they report,
 * and the whole degree by degree schedule the differential runs on, is
 * meaningless if the input is not homogeneous.  A module Gröbner basis is
 * not, which is why this guard sits here rather than in
 * module_gb_from_input.  st->homogeneous is set by
 * import_module_input_data and already accounts for the component shifts. */
static int module_input_is_graded(
        const md_t * const st
        )
{
    if (st->homogeneous) {
        return 1;
    }
    fprintf(ERRSTREAM, "The input is not homogeneous, so it has no graded "
            "free resolution; check the degree shifts of the ambient free "
            "module.\n");

    return 0;
}

void free_module_f4_result_data(
        void (*freep) (void *),
        int32_t **blen,
        int32_t **bexp,
        int32_t **bcomp,
        void **bcf
        )
{
    if (*blen != NULL) {
        (*freep)(*blen);
        *blen = NULL;
    }
    if (*bexp != NULL) {
        (*freep)(*bexp);
        *bexp = NULL;
    }
    if (*bcomp != NULL) {
        (*freep)(*bcomp);
        *bcomp = NULL;
    }
    if (*bcf != NULL) {
        (*freep)(*bcf);
        *bcf = NULL;
    }
}

int64_t export_module_f4(
        void *(*mallocp) (size_t),
        int32_t *bld,
        int32_t **blen,
        int32_t **bexp,
        int32_t **bcomp,
        void **bcf,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t info_level
        )
{
    int64_t nterms = 0;
    md_t *st = NULL;

    *bld   = 0;
    *blen  = NULL;
    *bexp  = NULL;
    *bcomp = NULL;
    *bcf   = NULL;

    bs_t *gb = module_gb_from_input(&st, lens, exps, comps, cfs, row_degs,
            field_char, mon_order, module_order, nr_vars, nr_rows, nr_gens,
            ht_size, nr_threads, max_nr_pairs, la_option, reduce_gb,
            info_level);
    if (gb == NULL) {
        return 0;
    }

    ht_t *bht = gb->ht;

    nterms = export_module_data(
            bld, blen, bexp, bcomp, bcf, mallocp, gb, bht, st);

    free_shared_hash_data(bht);
    free_basis(&gb);
    free(st);

    return nterms;
}

/* --------------------------------------------------------------------- *
 *  Schreyer frame of a presentation matrix
 * --------------------------------------------------------------------- */

void free_module_frame_result_data(
        void (*freep) (void *),
        int32_t **betti
        )
{
    if (*betti != NULL) {
        (*freep)(*betti);
        *betti = NULL;
    }
}

int64_t export_module_frame(
        void *(*mallocp) (size_t),
        int32_t *nlevels,
        int32_t *maxdeg,
        int32_t **betti,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t max_level,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level
        )
{
    int32_t i;
    int64_t nelts = 0;
    md_t *st      = NULL;

    *nlevels = 0;
    *maxdeg  = 0;
    *betti   = NULL;

    /* The frame is read off the lead terms of a *reduced* basis, which are
     * the minimal generators of the module of lead terms; a non-reduced
     * basis would carry redundant elements into level 1 and inflate every
     * level above it. */
    bs_t *gb = module_gb_from_input(&st, lens, exps, comps, cfs, row_degs,
            field_char, mon_order, module_order, nr_vars, nr_rows, nr_gens,
            ht_size, nr_threads, max_nr_pairs, la_option, 1 /* reduce */,
            info_level);
    if (gb == NULL) {
        return 0;
    }

    ht_t *bht = gb->ht;

    if (!module_input_is_graded(st)) {
        free_shared_hash_data(bht);
        free_basis(&gb);
        free(st);
        return 0;
    }

    /* Multigraded frames are M6; here the grading is the standard one and
     * the row degrees supply the only shifts, matching bht->cshift. */
    res_dgrp_t *grp = res_dgrp_new_standard(nr_vars);
    int32_t *rowmd  = (int32_t *)calloc(
            (unsigned long)nr_rows, sizeof(int32_t));
    res_frame_t *f  = grp != NULL ? res_frame_new(grp, st, max_level) : NULL;

    if (grp == NULL || rowmd == NULL || f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }
    for (i = 0; i < nr_rows; ++i) {
        rowmd[i] = (int32_t)bht->cshift[i+1];
    }

    if (res_frame_init(f, gb, bht, rowmd)) {
        goto cleanup;
    }
    nelts = res_frame_complete(f);
    if (nelts < 0 || res_frame_verify(f)) {
        nelts = 0;
        goto cleanup;
    }

    /* Tabulate into plain memory first: mallocp comes without a matching
     * free, so nothing allocated with it may be abandoned on a failure
     * path.  The caller only ever sees a table that is already known
     * good. */
    const deg_t mxd = res_frame_max_hdeg(f);
    if (mxd < 0) {
        fprintf(ERRSTREAM, "A frame has a negative maximum degree.\n");
        nelts = 0;
        goto cleanup;
    }
    const size_t nrows = (size_t)f->nlv;
    const size_t ncols = (size_t)mxd + 1;
    if (ncols > SIZE_MAX / sizeof(int32_t)
            || nrows > SIZE_MAX / (ncols * sizeof(int32_t))) {
        fprintf(ERRSTREAM, "The frame rank table is too large to allocate.\n");
        nelts = 0;
        goto cleanup;
    }
    const size_t tsz = nrows * ncols * sizeof(int32_t);
    int32_t *tmp = (int32_t *)malloc(tsz);
    if (tmp == NULL || res_frame_betti(f, tmp, mxd) != nelts) {
        fprintf(ERRSTREAM, "Frame ranks do not add up to the frame size.\n");
        free(tmp);
        nelts = 0;
        goto cleanup;
    }

    int32_t *tab = (int32_t *)(*mallocp)(tsz);
    if (tab == NULL) {
        free(tmp);
        nelts = 0;
        goto cleanup;
    }
    memcpy(tab, tmp, tsz);
    free(tmp);

    *nlevels = (int32_t)f->nlv;
    *maxdeg  = (int32_t)mxd;
    *betti   = tab;

cleanup:
    res_frame_free(&f);
    res_dgrp_free(&grp);
    free(rowmd);
    free_shared_hash_data(bht);
    free_basis(&gb);
    free(st);

    return nelts;
}

/* --------------------------------------------------------------------- *
 *  Nonminimal free resolutions and syzygies
 * --------------------------------------------------------------------- */

void free_module_resolution_result_data(
        void (*freep) (void *),
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        )
{
    if (*ranks != NULL) {
        (*freep)(*ranks);
        *ranks = NULL;
    }
    if (*degs != NULL) {
        (*freep)(*degs);
        *degs = NULL;
    }
    if (*dlen != NULL) {
        (*freep)(*dlen);
        *dlen = NULL;
    }
    if (*dexp != NULL) {
        (*freep)(*dexp);
        *dexp = NULL;
    }
    if (*dcomp != NULL) {
        (*freep)(*dcomp);
        *dcomp = NULL;
    }
    if (*dcf != NULL) {
        (*freep)(*dcf);
        *dcf = NULL;
    }
}

/* Everything is already computed and checked by the time this runs, so
 * the sizes are exact and nothing can fail after the first allocation --
 * which matters, since mallocp comes without a matching free. */
static int64_t export_resolution_data(
        void *(*mallocp) (size_t),
        int32_t *nlevels,
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf,
        const res_diff_t * const rd,
        const len_t nv
        )
{
    len_t i, k;
    int64_t t;

    const res_frame_t * const f = rd->f;
    const len_t nlv = f->nlv;

    int64_t ngens = 0, ncols = 0, nterms = 0;
    for (i = 0; i < nlv; ++i) {
        ngens += (int64_t)f->lv[i].ld;
        if (i == 0) {
            continue;
        }
        ncols += (int64_t)f->lv[i].ld;
        for (k = 0; k < f->lv[i].ld; ++k) {
            nterms += (int64_t)rd->d[i][k].len;
        }
    }

    if (nterms > (int64_t)INT32_MAX / (nv > 0 ? (int64_t)nv : 1)) {
        fprintf(ERRSTREAM,
                "The resolution is too large to store in flat arrays.\n");
        return 0;
    }

    int32_t *rk = (int32_t *)(*mallocp)((unsigned long)nlv * sizeof(int32_t));
    int32_t *dg = (int32_t *)(*mallocp)(
            (unsigned long)ngens * sizeof(int32_t));
    int32_t *dl = (int32_t *)(*mallocp)(
            (unsigned long)(ncols > 0 ? ncols : 1) * sizeof(int32_t));
    int32_t *de = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1)
            * (unsigned long)nv * sizeof(int32_t));
    int32_t *dc = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1) * sizeof(int32_t));
    int32_t *cf = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1) * sizeof(int32_t));

    int64_t cg = 0, cl = 0, ce = 0, cc = 0;

    for (i = 0; i < nlv; ++i) {
        rk[i] = (int32_t)f->lv[i].ld;
        for (k = 0; k < f->lv[i].ld; ++k) {
            dg[cg++] = (int32_t)f->lv[i].elts[k].hdeg;
        }
    }
    for (i = 1; i < nlv; ++i) {
        for (k = 0; k < f->lv[i].ld; ++k) {
            const res_dpoly_t * const p = rd->d[i] + k;
            dl[cl++] = (int32_t)p->len;
            for (t = 0; t < p->len; ++t) {
                const exp_t * const ev = f->ht->ev[p->mon[t]];
                len_t j;
                for (j = 1; j <= nv; ++j) {
                    de[ce++] = (int32_t)ev[j];
                }
                dc[cc] = p->pos[t] + 1;
                cf[cc] = (int32_t)p->cf[t];
                cc++;
            }
        }
    }

    *nlevels = (int32_t)nlv;
    *ranks   = rk;
    *degs    = dg;
    *dlen    = dl;
    *dexp    = de;
    *dcomp   = dc;
    *dcf     = (void *)cf;

    return nterms;
}

/* --- syzygies of the input generators -------------------------------- *
 *
 * The graph module trick.  Adjoin one new component per input generator
 * and one extra term e_{nr_rows+j} to generator j; under position over
 * term with the original components first, an element of the Gröbner
 * basis whose original components all vanish is exactly a relation among
 * the generators, read off its new components.  No frame and no
 * differential are involved -- this is a module Gröbner basis and nothing
 * else, which is the whole reason msolve can do it at all despite
 * discarding the change of basis from the input to the Gröbner basis. */
static int64_t module_syz_of_input(
        void *(*mallocp) (size_t),
        int32_t *nlevels,
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level,
        const int32_t max_level
        )
{
    int32_t i, j;
    int64_t t, nterms = 0;
    md_t *st = NULL;
    bs_t *gb = NULL;

    int32_t *lens2 = NULL, *exps2 = NULL, *comps2 = NULL, *cfs2 = NULL;
    int32_t *rdeg  = NULL, *sylen = NULL;
    int64_t *syidx = NULL;

    if (nr_rows < 1 || nr_gens < 1 || nr_vars < 1) {
        fprintf(ERRSTREAM, "Empty module input.\n");
        return 0;
    }

    int64_t nt = 0;
    for (i = 0; i < nr_gens; ++i) {
        if (lens[i] < 1) {
            fprintf(ERRSTREAM, "Generator %d has no terms.\n", i);
            return 0;
        }
        nt += lens[i];
    }
    for (t = 0; t < nt; ++t) {
        if (comps[t] < 1 || comps[t] > nr_rows) {
            fprintf(ERRSTREAM, "Term %ld has component %d, outside the range "
                    "1 to %d.\n", (long)t, comps[t], nr_rows);
            return 0;
        }
    }

    if (nr_gens > INT32_MAX - nr_rows || nt > INT32_MAX - nr_gens) {
        fprintf(ERRSTREAM, "The graph module is too large to index.\n");
        goto cleanup;
    }
    const int32_t nr2 = nr_rows + nr_gens;
    const int64_t nt2 = nt + nr_gens;
    if ((uint64_t)nt2 > SIZE_MAX / sizeof(int32_t)
            || (uint64_t)nt2 > SIZE_MAX / (uint64_t)nr_vars
                / sizeof(int32_t)) {
        fprintf(ERRSTREAM, "The graph module is too large to allocate.\n");
        goto cleanup;
    }

    lens2  = (int32_t *)malloc((unsigned long)nr_gens * sizeof(int32_t));
    exps2  = (int32_t *)calloc(
            (size_t)nt2 * (size_t)nr_vars,
            sizeof(int32_t));
    comps2 = (int32_t *)malloc(
            (size_t)nt2 * sizeof(int32_t));
    cfs2   = (int32_t *)malloc(
            (size_t)nt2 * sizeof(int32_t));
    rdeg   = (int32_t *)malloc((unsigned long)nr2 * sizeof(int32_t));
    if (lens2 == NULL || exps2 == NULL || comps2 == NULL
            || cfs2 == NULL || rdeg == NULL) {
        goto cleanup;
    }

    /* Degrees, normalized to start at zero exactly as export_module_f4
     * does, so that both entry points report the same table. */
    int32_t mn = 0;
    if (row_degs != NULL) {
        mn = row_degs[0];
        for (i = 1; i < nr_rows; ++i) {
            if (row_degs[i] < mn) {
                mn = row_degs[i];
            }
        }
    }
    for (i = 0; i < nr_rows; ++i) {
        rdeg[i] = row_degs != NULL ? row_degs[i] - mn : 0;
    }

    const int32_t *icf = (const int32_t *)cfs;
    int64_t off = 0, off2 = 0;
    for (i = 0; i < nr_gens; ++i) {
        lens2[i] = lens[i] + 1;
        /* the generator's degree, read off its first term; for the graded
         * input this machinery is about, every term agrees */
        int32_t d = rdeg[comps[off] - 1];
        for (j = 0; j < nr_vars; ++j) {
            d += exps[off * nr_vars + j];
        }
        rdeg[nr_rows + i] = d;

        for (t = 0; t < lens[i]; ++t) {
            for (j = 0; j < nr_vars; ++j) {
                exps2[(off2 + t) * nr_vars + j] = exps[(off + t) * nr_vars + j];
            }
            comps2[off2 + t] = comps[off + t];
            cfs2[off2 + t]   = icf[off + t];
        }
        /* ... and the tautological term e_{nr_rows + i} */
        comps2[off2 + lens[i]] = nr_rows + i + 1;
        cfs2[off2 + lens[i]]   = 1;

        off  += lens[i];
        off2 += lens2[i];
    }

    gb = module_gb_from_input(&st, lens2, exps2, comps2, cfs2, rdeg,
            field_char, mon_order, module_order, nr_vars, nr2, nr_gens,
            ht_size, nr_threads, max_nr_pairs, la_option, 1 /* reduce */,
            info_level);
    if (gb == NULL) {
        goto cleanup;
    }
    if (!module_input_is_graded(st)) {
        goto cleanup;
    }

    ht_t *bht = gb->ht;

    /* pick out the elements supported in the adjoined components */
    sylen = (int32_t *)malloc(
            (unsigned long)(gb->lml > 0 ? gb->lml : 1) * sizeof(int32_t));
    syidx = (int64_t *)malloc(
            (unsigned long)(gb->lml > 0 ? gb->lml : 1) * sizeof(int64_t));
    if (sylen == NULL || syidx == NULL) {
        goto cleanup;
    }

    int32_t nsyz = 0;
    int64_t syterms = 0;
    for (i = 0; i < (int32_t)gb->lml; ++i) {
        const bl_t bi   = gb->lmps[i];
        const hm_t *hm  = gb->hm[bi];
        if (hm == NULL) {
            continue;
        }
        /* Position over term puts the original components first, so the
         * lead term alone would decide.  Every term is looked at anyway:
         * that turns "the original components vanish" from something
         * inferred from the order into something seen, and it is exactly
         * the statement that this column is a relation among the
         * generators, which is why the syzygies of the input need no
         * separate d_1 o d_2 = 0 check. */
        int pure = 1;
        for (t = 0; t < (int64_t)hm[LENGTH]; ++t) {
            if ((int32_t)bht->ev[hm[OFFSET+t]][bht->cpos] <= nr_rows) {
                pure = 0;
                break;
            }
        }
        if (!pure) {
            continue;
        }
        syidx[nsyz] = bi;
        sylen[nsyz] = (int32_t)hm[LENGTH];
        syterms    += (int64_t)hm[LENGTH];
        nsyz++;
    }

    const int32_t nlv = (max_level == 1 || nsyz == 0) ? 2 : 3;
    if (nlv == 2) {
        syterms = 0;
        nsyz    = 0;
    }

    nterms = nt + syterms;
    if (nterms > (int64_t)INT32_MAX / (int64_t)nr_vars) {
        fprintf(ERRSTREAM,
                "The resolution is too large to store in flat arrays.\n");
        nterms = 0;
        goto cleanup;
    }

    const int64_t ngens = (int64_t)nr_rows + nr_gens + nsyz;
    const int64_t ncols = (int64_t)nr_gens + nsyz;

    int32_t *rk = (int32_t *)(*mallocp)((unsigned long)nlv * sizeof(int32_t));
    int32_t *dg = (int32_t *)(*mallocp)(
            (unsigned long)ngens * sizeof(int32_t));
    int32_t *dl = (int32_t *)(*mallocp)(
            (unsigned long)ncols * sizeof(int32_t));
    int32_t *de = (int32_t *)(*mallocp)(
            (unsigned long)nterms * (unsigned long)nr_vars * sizeof(int32_t));
    int32_t *dc = (int32_t *)(*mallocp)(
            (unsigned long)nterms * sizeof(int32_t));
    int32_t *cf = (int32_t *)(*mallocp)(
            (unsigned long)nterms * sizeof(int32_t));

    rk[0] = nr_rows;
    rk[1] = nr_gens;
    if (nlv == 3) {
        rk[2] = nsyz;
    }

    int64_t cg = 0, cl = 0, ce = 0, cc = 0;
    for (i = 0; i < nr2; ++i) {
        dg[cg++] = rdeg[i];
    }

    /* d_1 is the caller's matrix, only reduced into [0, p) */
    const uint32_t fc = st->fc;
    off = 0;
    for (i = 0; i < nr_gens; ++i) {
        dl[cl++] = lens[i];
        for (t = 0; t < lens[i]; ++t) {
            for (j = 0; j < nr_vars; ++j) {
                de[ce++] = exps[(off + t) * nr_vars + j];
            }
            dc[cc] = comps[off + t];
            int64_t v = (int64_t)icf[off + t] % (int64_t)fc;
            cf[cc]    = (int32_t)(v < 0 ? v + (int64_t)fc : v);
            cc++;
        }
        off += lens[i];
    }

    /* d_2 is the syzygies, with the adjoined components folded back onto
     * the generators they index */
    for (i = 0; i < nsyz; ++i) {
        const hm_t *hm = gb->hm[syidx[i]];
        dg[cg++] = (int32_t)bht->hd[hm[OFFSET]].deg;
        dl[cl++] = sylen[i];
        for (t = 0; t < sylen[i]; ++t) {
            const exp_t * const ev = bht->ev[hm[OFFSET+t]];
            for (j = 1; j <= nr_vars; ++j) {
                de[ce++] = (int32_t)ev[j];
            }
            dc[cc] = (int32_t)ev[bht->cpos] - nr_rows;
            switch (st->ff_bits) {
                case 8:
                    cf[cc] = (int32_t)gb->cf_8[hm[COEFFS]][t];
                    break;
                case 16:
                    cf[cc] = (int32_t)gb->cf_16[hm[COEFFS]][t];
                    break;
                default:
                    cf[cc] = (int32_t)gb->cf_32[hm[COEFFS]][t];
                    break;
            }
            cc++;
        }
    }

    *nlevels = nlv;
    *ranks   = rk;
    *degs    = dg;
    *dlen    = dl;
    *dexp    = de;
    *dcomp   = dc;
    *dcf     = (void *)cf;

cleanup:
    free(lens2);
    free(exps2);
    free(comps2);
    free(cfs2);
    free(rdeg);
    free(sylen);
    free(syidx);
    if (gb != NULL) {
        free_shared_hash_data(gb->ht);
        free_basis(&gb);
    }
    free(st);

    return nterms;
}

int64_t export_module_resolution(
        void *(*mallocp) (size_t),
        int32_t *nlevels,
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t module_order,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t max_level,
        const int32_t syz_of,
        const int32_t verify,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level
        )
{
    int32_t i;
    int64_t nterms = 0;
    md_t *st       = NULL;

    *nlevels = 0;
    *ranks   = NULL;
    *degs    = NULL;
    *dlen    = NULL;
    *dexp    = NULL;
    *dcomp   = NULL;
    *dcf     = NULL;

    if (module_order != RES_MORD_POT) {
        fprintf(ERRSTREAM, "Resolutions are only implemented for the "
                "position over term module order: the Schreyer order the "
                "differential runs in is the one that induces, and term "
                "over position would need the component degree shifts "
                "which the frame's ring hash table does not carry.\n");
        return 0;
    }
    if (max_level < 0) {
        fprintf(ERRSTREAM, "A negative truncation level makes no sense.\n");
        return 0;
    }

    if (syz_of == RES_SYZ_OF_INPUT) {
        if (max_level > 2) {
            fprintf(ERRSTREAM, "Syzygies of the input generators stop at "
                    "level 2; resolving further is the Gröbner basis story "
                    "again, so ask for RES_SYZ_OF_GB.\n");
            return 0;
        }
        return module_syz_of_input(mallocp, nlevels, ranks, degs,
                dlen, dexp, dcomp, dcf, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, module_order, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                info_level, max_level == 0 ? 2 : max_level);
    }
    if (syz_of != RES_SYZ_OF_GB) {
        fprintf(ERRSTREAM, "Unknown syzygy flavour %d.\n", syz_of);
        return 0;
    }

    /* As for the frame: the lead terms have to be the minimal generators
     * of the module of lead terms, which is what reducing gives. */
    bs_t *gb = module_gb_from_input(&st, lens, exps, comps, cfs, row_degs,
            field_char, mon_order, module_order, nr_vars, nr_rows, nr_gens,
            ht_size, nr_threads, max_nr_pairs, la_option, 1 /* reduce */,
            info_level);
    if (gb == NULL) {
        return 0;
    }

    ht_t *bht = gb->ht;

    if (!module_input_is_graded(st)) {
        free_shared_hash_data(bht);
        free_basis(&gb);
        free(st);
        return 0;
    }

    res_dgrp_t *grp  = res_dgrp_new_standard(nr_vars);
    int32_t *rowmd   = (int32_t *)calloc(
            (unsigned long)nr_rows, sizeof(int32_t));
    res_frame_t *f   = grp != NULL ? res_frame_new(grp, st, max_level) : NULL;
    res_diff_t *rd   = NULL;

    if (grp == NULL || rowmd == NULL || f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }
    for (i = 0; i < nr_rows; ++i) {
        rowmd[i] = (int32_t)bht->cshift[i+1];
    }

    if (res_frame_init(f, gb, bht, rowmd)) {
        goto cleanup;
    }
    if (res_frame_complete(f) < 0 || res_frame_verify(f)) {
        goto cleanup;
    }

    rd = res_diff_new(f, st->fc);
    if (rd == NULL) {
        fprintf(ERRSTREAM, "Could not set up the differential.\n");
        goto cleanup;
    }
    if (res_diff_init(rd, gb, bht, st) || res_diff_compute(rd)) {
        goto cleanup;
    }
    if (res_diff_verify(rd, verify != 0)) {
        fprintf(ERRSTREAM, "The computed differential is not a complex.\n");
        goto cleanup;
    }

    nterms = export_resolution_data(mallocp, nlevels, ranks, degs,
            dlen, dexp, dcomp, dcf, rd, (len_t)nr_vars);

cleanup:
    res_diff_free(&rd);
    res_frame_free(&f);
    res_dgrp_free(&grp);
    free(rowmd);
    free_shared_hash_data(bht);
    free_basis(&gb);
    free(st);

    return nterms;
}
