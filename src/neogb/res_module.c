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

    /* validate the components before anything is allocated */
    int64_t nt = 0;
    for (i = 0; i < nr_gens; ++i) {
        if (lens[i] < 1) {
            fprintf(ERRSTREAM, "Generator %d has no terms.\n", i);
            return NULL;
        }
        nt += lens[i];
    }
    for (j = 0; j < nt; ++j) {
        if (comps[j] < 1 || comps[j] > nr_rows) {
            fprintf(ERRSTREAM, "Term %d has component %d, outside the range "
                    "1 to %d.\n", j, comps[j], nr_rows);
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
        int32_t mn = row_degs[0];
        for (i = 1; i < nr_rows; ++i) {
            if (row_degs[i] < mn) {
                mn = row_degs[i];
            }
        }
        for (i = 0; i < nr_rows; ++i) {
            bht->cshift[i+1] = (deg_t)(row_degs[i] - mn);
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
