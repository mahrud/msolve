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
 * relies on (see data.h).
 *
 * With a block order there is one degree slot per block; the component
 * shift goes into the first of them, matching get_lcm in hash.c, so that
 * the module degree is the shifted one in exactly one place and the
 * remaining blocks carry their plain ring degrees. */
static inline void set_module_exponent_vector(
        exp_t *ev,
        const int32_t *iev,
        const int32_t *icomp,
        const int32_t idx,
        const ht_t *ht
        )
{
    len_t i;

    const len_t nv          = ht->nv;
    const len_t * const bst = ht->bst;
    const exp_t comp        = (exp_t)icomp[idx];
    const deg_t * const vwt = ht->vwt;
    const int32_t * const e = iev + (nv * idx);

    len_t ctr = 0;
    for (len_t b = 0; b < ht->nbl; ++b) {
        const len_t d   = bst[b];
        const len_t end = bst[b+1];
        /* a block's degree slot is the *heft* degree of that block:
         * under the standard grading every weight is one and this is
         * msolve's usual total degree, and vwt is NULL so the
         * multiplication is not even compiled into the path input takes */
        deg_t deg = 0;
        for (i = d+1; i < end; ++i) {
            ev[i] = (exp_t)e[ctr++];
            deg += vwt == NULL ? (deg_t)ev[i] : vwt[i] * (deg_t)ev[i];
        }
        ev[d] = (exp_t)deg;
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

    /* the exponent vector slot of each variable, which under a block
     * order is not simply k+1: the per block degree slots sit in between */
    int32_t *evi = (int32_t *)malloc((unsigned long)nv * sizeof(int32_t));
    if (evi == NULL) {
        return 0;
    }
    ht_variable_slots(ht, evi);

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
        free(evi);
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
                free(evi);
                return 0;
        }
        dt = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 0; k < nv; ++k) {
                exp[ce++] = (int32_t)ht->ev[dt[j]][evi[k]];
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

    free(evi);

    return nterms;
}

/* Everything the module entry points have in common: the Gröbner basis,
 * the grading it was computed under, and the normalized row degrees.  All
 * of it is owned by the caller and released together by
 * module_input_clear, which is safe on a partially filled struct. */
typedef struct module_input_t module_input_t;
struct module_input_t
{
    md_t       *st;
    bs_t       *gb;
    res_dgrp_t *grp;
    int32_t    *rowmd;    /* nr_rows * grp->len, normalized as below   */
    int32_t     degshift[RES_MTAB_MAXLEN];  /* what was subtracted     */
    int         graded;   /* every generator is multihomogeneous       */
};

static void module_input_clear(
        module_input_t *mi
        )
{
    if (mi->gb != NULL) {
        free_shared_hash_data(mi->gb->ht);
        free_basis(&mi->gb);
    }
    res_dgrp_free(&mi->grp);
    free(mi->rowmd);
    free(mi->st);
    memset(mi, 0, sizeof(module_input_t));
}

/* The heft degree of a raw multidegree.  heft_of never writes through the
 * view, so handing it a pointer into a const array is sound. */
static inline deg_t module_row_heft(
        const res_dgrp_t * const g,
        const int32_t * const d
        )
{
    res_deg_t a;

    a.e = (int32_t *)d;

    return g->heft_of(g, a);
}

/* The same in 64 bits, for validating a caller's row degrees *before* any
 * of them is subtracted.  deg_t is an int32_t and the caller's degrees are
 * arbitrary int32_t, so both the heft and the difference of two row degrees
 * can overflow -- and a row degree of INT32_MIN against one of INT32_MAX is
 * exactly the input a caller would hand over to see it refused.  Refusing
 * it must not itself be undefined. */
static inline int64_t module_row_heft64(
        const res_dgrp_t * const g,
        const int32_t * const d
        )
{
    len_t i;
    int64_t h = 0;

    for (i = 0; i < g->r; ++i) {
        h += (int64_t)g->heft[i] * (int64_t)d[i];
    }

    return h;
}

/* The multidegree of one input term: the degree of its monomial plus the
 * degree of the row it sits in.  res_deg_of_exponents works on the hash
 * table's 16-bit exponents, and this is the raw int32 input, so the column
 * sum is spelled out here rather than shared.  Torsion is reduced by the
 * group's own addition, which is the only place that knows how. */
static void module_term_multidegree(
        const res_dgrp_t * const g,
        int32_t *out,
        const int32_t * const exps,
        const int32_t * const rowmd
        )
{
    len_t i, k;
    const len_t r    = g->r;
    const len_t glen = g->len;
    const len_t nv   = g->nv;

    int64_t acc[RES_MTAB_MAXLEN];

    memset(acc, 0, (unsigned long)glen * sizeof(int64_t));
    memset(out, 0, (unsigned long)glen * sizeof(int32_t));
    for (k = 0; k < nv; ++k) {
        const int32_t e = exps[k];
        if (e == 0) {
            continue;
        }
        const int32_t * const col = g->dmat + (unsigned long)k * glen;
        for (i = 0; i < r; ++i) {
            acc[i] += (int64_t)e * (int64_t)col[i];
        }
        for (i = r; i < glen; ++i) {
            /* reduced every step, and the product formed in 64 bits, so
             * neither the residue nor the multiplication can run away */
            const int32_t t = g->tord[i-r];
            out[i] = res_mod_torsion(
                    out[i] + (int32_t)(((int64_t)e * col[i]) % t), t);
        }
    }
    for (i = 0; i < r; ++i) {
        out[i] = (int32_t)acc[i];
    }

    res_deg_t a = {out};
    res_deg_t b = {(int32_t *)rowmd};
    g->add(g, a, a, b);
}

/* Validate the input, build the grading, normalize the row degrees, build
 * the meta data, import the presentation matrix and run F4.  Returns 0 on
 * success, filling *mi; on failure everything allocated is released and
 * *mi is zeroed, so module_input_clear on it is harmless either way.
 *
 * The strategy is validated here by res_strat_check, so every entry point
 * that reaches this one is checked before any work happens; NULL means
 * res_strat_default().  Its base is restricted to POT and TOP, the
 * Schreyer order as a base needing per component base monomials which
 * only the frame in res_frame.c can supply, and by then the Gröbner basis
 * is already in hand.
 *
 * The grading is validated the same way, and NULL means the standard one.
 * Row degrees are normalized by subtracting the multidegree of the row of
 * least heft, so that every component shift is a nonnegative heft degree
 * that fits the 16-bit slot the hash table keeps it in -- which is the
 * same normalization as before, since under the standard grading the heft
 * of a row degree is the row degree. */
static int module_gb_from_input(
        module_input_t *mi,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const mo_block_t * const blk,
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop,
        const int32_t syz_comp_lo,
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

    memset(mi, 0, sizeof(module_input_t));

    if (res_strat_check(strat, 1)) {
        return 1;
    }
    const res_strat_t sdef = res_strat_default();
    const res_strat_t * const sp = strat != NULL ? strat : &sdef;

    if (field_char == 0 || field_char >= ((uint32_t)1 << 31)) {
        fprintf(ERRSTREAM, "Module Groebner bases need a prime field of "
                "characteristic less than 2^31.\n");
        return 1;
    }
    if (nr_rows < 1 || nr_gens < 1 || nr_vars < 1) {
        fprintf(ERRSTREAM, "Empty module input.\n");
        return 1;
    }
    if (mon_order != 0) {
        fprintf(ERRSTREAM, "Module Groebner bases are only implemented for "
                "the degree reverse lexicographic order.\n");
        return 1;
    }

    /* --- the grading ------------------------------------------------- */

    res_dgrp_t *grp = res_dgrp_of_grading(grading, nr_vars);
    if (grp == NULL) {
        return 1;
    }
    mi->grp = grp;

    const len_t glen = grp->len;

    /* Validate all dimensions before anything else is allocated.  The
     * importer below uses int32_t offsets, while the hash table stores both
     * individual exponents and shifted heft degrees in exp_t.  Reject values
     * that would otherwise be silently narrowed to a different monomial. */
    int64_t nt = 0;
    for (i = 0; i < nr_gens; ++i) {
        if (lens[i] < 1) {
            fprintf(ERRSTREAM, "Generator %d has no terms.\n", i);
            goto fail;
        }
        nt += lens[i];
        if (nt > INT32_MAX) {
            fprintf(ERRSTREAM, "Module input has too many terms.\n");
            goto fail;
        }
    }

    /* --- row degrees, normalized to the row of least heft ------------- */

    int32_t *rowmd = (int32_t *)calloc(
            (unsigned long)nr_rows * (unsigned long)glen, sizeof(int32_t));
    if (rowmd == NULL) {
        goto fail;
    }
    mi->rowmd = rowmd;

    if (row_degs != NULL) {
        int32_t i0  = 0;
        int64_t hmn = module_row_heft64(grp, row_degs);
        for (i = 1; i < nr_rows; ++i) {
            const int64_t h = module_row_heft64(
                    grp, row_degs + (int64_t)i * glen);
            if (h < hmn) {
                hmn = h;
                i0  = i;
            }
        }

        /* Everything is checked in 64 bits first, against the row that was
         * picked, and only then subtracted: both the heft difference and
         * the free part of the difference are int32_t quantities that a
         * caller's degrees can overflow, and the whole point of the check
         * is to be handed degrees that do. */
        const int32_t * const base = row_degs + (int64_t)i0 * glen;
        for (i = 0; i < nr_rows; ++i) {
            const int32_t * const rd = row_degs + (int64_t)i * glen;
            const int64_t shift = module_row_heft64(grp, rd) - hmn;
            if (shift < 0 || shift > UINT16_MAX) {
                fprintf(ERRSTREAM, "Row %d has heft degree %ld relative to "
                        "the lightest row, which does not fit in the "
                        "exponent table.\n", i, (long)shift);
                goto fail;
            }
            for (j = 0; j < grp->r; ++j) {
                const int64_t v = (int64_t)rd[j] - (int64_t)base[j];
                if (v < INT32_MIN || v > INT32_MAX) {
                    fprintf(ERRSTREAM, "Row %d has a degree too far from the "
                            "lightest row to normalize against it.\n", i);
                    goto fail;
                }
            }
        }

        memcpy(mi->degshift, base, (unsigned long)glen * sizeof(int32_t));

        res_deg_t vbase = {mi->degshift};
        for (i = 0; i < nr_rows; ++i) {
            res_deg_t d = {rowmd + (int64_t)i * glen};
            res_deg_t s = {(int32_t *)(row_degs + (int64_t)i * glen)};
            grp->sub(grp, d, s, vbase);
        }
    }

    /* --- every term: bounds, and the multidegree it sits in ----------- */

    const deg_t * const vhdeg = grp->vhdeg;

    for (j = 0; j < nt; ++j) {
        if (comps[j] < 1 || comps[j] > nr_rows) {
            fprintf(ERRSTREAM, "Term %d has component %d, outside the range "
                    "1 to %d.\n", j, comps[j], nr_rows);
            goto fail;
        }
        int64_t deg = module_row_heft(grp, rowmd + (int64_t)(comps[j]-1) * glen);
        for (i = 0; i < nr_vars; ++i) {
            const int32_t exponent = exps[(int64_t)j * nr_vars + i];
            if (exponent < 0 || exponent > UINT16_MAX) {
                fprintf(ERRSTREAM, "Term %d has exponent %d, outside the "
                        "16-bit exponent range.\n", j, exponent);
                goto fail;
            }
            deg += (int64_t)exponent * (int64_t)vhdeg[i];
        }
        if (deg > UINT16_MAX) {
            fprintf(ERRSTREAM, "Term %d has shifted heft degree %ld, outside "
                    "the 16-bit exponent range.\n", j, (long)deg);
            goto fail;
        }
    }

    /* Multihomogeneity.  With the standard grading this is exactly the
     * heft homogeneity import_module_input_data goes on to compute, so
     * nothing changes there; with any other it is strictly stronger, and it
     * is what every graded output here means.  Checking it on the input
     * rather than on the basis is the cheap place: the Gröbner basis of a
     * multihomogeneous module is multihomogeneous. */
    {
        int32_t *da = (int32_t *)calloc((unsigned long)glen, sizeof(int32_t));
        int32_t *db = (int32_t *)calloc((unsigned long)glen, sizeof(int32_t));
        if (da == NULL || db == NULL) {
            free(da);
            free(db);
            goto fail;
        }
        mi->graded  = 1;
        int64_t off = 0;
        for (i = 0; i < nr_gens; ++i) {
            for (j = 0; j < lens[i]; ++j) {
                const int64_t t = off + j;
                int32_t * const cur = j == 0 ? da : db;
                module_term_multidegree(grp, cur,
                        exps + t * nr_vars,
                        rowmd + (int64_t)(comps[t]-1) * glen);
                if (j > 0 && memcmp(da, db,
                            (unsigned long)glen * sizeof(int32_t)) != 0) {
                    mi->graded = 0;
                    break;
                }
            }
            if (!mi->graded) {
                break;
            }
            off += lens[i];
        }
        free(da);
        free(db);
    }

    md_t *st = allocate_meta_data();
    if (st == NULL) {
        goto fail;
    }
    mi->st = st;

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
        goto fail;
    }
    if (check_and_set_meta_data(st, lens, exps, cfs, invalid_gens,
                l_field_char, l_mon_order, elim_block_len, l_nr_vars,
                l_nr_gens, nr_nf, l_ht_size, l_nr_threads, l_max_nr_pairs,
                reset_ht, l_la_option, use_signatures, l_reduce_gb, pbm_file,
                truncate_lifting, l_info_level)) {
        free(invalid_gens);
        goto fail;
    }

    /* A general block order replaces the single block default.  cmp_blocks
     * reads only the slots below ht->cpos, so the module order keeps the
     * component key it always had and only its ring part changes; see
     * res_cmp_terms_ring. */
    if (set_monomial_block_order(st, blk)) {
        free(invalid_gens);
        goto fail;
    }

    /* this is what makes the hash table a module one; it has to happen
     * before initialize_basis, which is where the table is built */
    st->ncomp = nr_rows;
    st->mord  = sp->base;
    st->mpos  = sp->pos;
    st->mlift = sp->lift;

    bs_t *bs  = initialize_basis(st, NULL);
    ht_t *bht = bs->ht;

    /* Variable weights, which a weighted block order may already have
     * put here; see res_install_weights. */
    if (res_install_weights(bht, grp)) {
        free(invalid_gens);
        free_shared_hash_data(bht);
        free_basis(&bs);
        goto fail;
    }

    /* Degree shifts of the ambient free module, as heft degrees, which is
     * what ev[DEG] carries.  They are normalized to start at zero because
     * msolve keeps a degree in a uint16_t; a global shift of all row
     * degrees changes neither the module nor its Gröbner basis, only the
     * absolute degrees, which the caller puts back from degshift. */
    for (i = 0; i < nr_rows; ++i) {
        bht->cshift[i+1] = module_row_heft(grp, rowmd + (int64_t)i * glen);
    }

    import_module_input_data(bs, st, 0, st->ngens_input,
            lens, exps, comps, cfs, invalid_gens);

    print_initial_statistics(VERBSTREAM, st);

    calculate_divmask(bht);

    sort_r(bs->hm, (unsigned long)bs->ld, sizeof(hm_t *),
            initial_input_cmp, bht);
    normalize_initial_basis(bs, st->fc);

    /* A degree ceiling.  The caller's multidegree is on the caller's own
     * scale and the computation runs on the normalized one, so subtract
     * what the rows were shifted by; and the schedule is by heft, so the
     * ceiling that reaches the round loop is a heft degree.  See
     * res_stop_t for why the coarsening is the safe direction.
     *
     * A ceiling at or below the shift leaves the round loop nothing it
     * could ever select.  That is refused rather than answered with the
     * input generators, because the only way to ask for it is to have the
     * scale wrong -- forgetting the degree shift, or handing over a heft
     * where a multidegree was wanted.
     *
     * And a ceiling only says what it appears to say about a *homogeneous*
     * computation.  For inhomogeneous input the degree the round loop
     * schedules by is a sugar degree that can fall -- which is what
     * md->min_deg_in_first_deg_fall is there to notice -- so "complete
     * through degree d" is not a true statement about the result.  Refused,
     * rather than silently meaning something else. */
    if (stop != NULL && stop->max_degree != NULL) {
        const int64_t h = module_row_heft64(grp, stop->max_degree)
            - module_row_heft64(grp, mi->degshift);
        if (!st->homogeneous) {
            fprintf(ERRSTREAM, "A degree limit needs homogeneous input: the "
                    "degree an inhomogeneous computation schedules by is a "
                    "sugar degree and can fall, so a ceiling on it does not "
                    "mean the basis is complete through that degree.\n");
            free(invalid_gens);
            free_shared_hash_data(bht);
            free_basis(&bs);
            goto fail;
        }
        if (h <= 0) {
            fprintf(ERRSTREAM, "The degree limit is at or below the degree "
                    "of the lightest generator of the ambient free module, "
                    "so no S-pair could be selected.\n");
            free(invalid_gens);
            free_shared_hash_data(bht);
            free_basis(&bs);
            goto fail;
        }
        st->max_gb_degree = (deg_t)(h > (int64_t)INT32_MAX ? INT32_MAX : h);
    }

    /* A syzygy limit only means anything where there is an adjoined block
     * to read relations off, which the caller says by passing its lower
     * boundary; everywhere else the round loop has nothing to count. */
    if (stop != NULL && stop->syz_limit > 0 && syz_comp_lo > 0) {
        st->syz_limit   = stop->syz_limit;
        st->syz_comp_lo = syz_comp_lo;
    }

    int32_t err = 0;
    bs_t *gb = core_gba(bs, st, &err, (len_t)field_char);

    free(invalid_gens);

    if (gb == NULL || err > 0) {
        fprintf(ERRSTREAM, "Module Groebner basis computation failed.\n");
        mi->gb = gb;
        goto fail;
    }

    mi->gb = gb;

    return 0;

fail:
    module_input_clear(mi);

    return 1;
}

/* Frames and resolutions are graded objects: every degree they report,
 * and the whole degree by degree schedule the differential runs on, is
 * meaningless if the input is not homogeneous.  A module Gröbner basis is
 * not, which is why this guard sits here rather than in
 * module_gb_from_input.  st->homogeneous is set by
 * import_module_input_data and already accounts for the component shifts. */
static int module_input_is_graded(
        const module_input_t * const mi
        )
{
    if (mi->st->homogeneous && mi->graded) {
        return 1;
    }
    fprintf(ERRSTREAM, "The input is not homogeneous, so it has no graded "
            "free resolution; check the degree shifts of the ambient free "
            "module%s.\n", mi->grp->simple
            ? "" : " and the degrees of the variables");

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

int64_t export_module_f4_blocks(
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
        const mo_block_t * const blk,
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop,
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
    module_input_t mi;

    *bld   = 0;
    *blen  = NULL;
    *bexp  = NULL;
    *bcomp = NULL;
    *bcf   = NULL;

    if (module_gb_from_input(&mi, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, blk, strat, grading, stop, 0, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                reduce_gb, info_level)) {
        return 0;
    }

    nterms = export_module_data(
            bld, blen, bexp, bcomp, bcf, mallocp, mi.gb, mi.gb->ht, mi.st);

    module_input_clear(&mi);

    return nterms;
}

/* The single block order is the general one with no block description. */
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
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop,
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
    return export_module_f4_blocks(mallocp, bld, blen, bexp, bcomp, bcf,
            lens, exps, comps, cfs, row_degs, field_char, mon_order, NULL,
            strat, grading, stop, nr_vars, nr_rows, nr_gens, ht_size,
            nr_threads, max_nr_pairs, la_option, reduce_gb, info_level);
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
        const res_strat_t * const strat,
        const res_grading_t * const grading,
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
    int64_t nelts = 0;
    module_input_t mi;
    res_frame_t *f = NULL;

    *nlevels = 0;
    *maxdeg  = 0;
    *betti   = NULL;

    /* The frame is read off the lead terms of a *reduced* basis, which are
     * the minimal generators of the module of lead terms; a non-reduced
     * basis would carry redundant elements into level 1 and inflate every
     * level above it. */
    if (module_gb_from_input(&mi, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, NULL, strat, grading, NULL, 0, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                1 /* reduce */, info_level)) {
        return 0;
    }

    ht_t *bht = mi.gb->ht;

    if (!module_input_is_graded(&mi)) {
        module_input_clear(&mi);
        return 0;
    }

    f = res_frame_new(mi.grp, mi.st, max_level);
    if (f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }

    if (res_frame_init(f, mi.gb, bht, mi.rowmd)) {
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
    module_input_clear(&mi);

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
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop,
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
    module_input_t mi;

    int32_t *lens2 = NULL, *exps2 = NULL, *comps2 = NULL, *cfs2 = NULL;
    int32_t *rdeg  = NULL, *sylen = NULL;
    int64_t *syidx = NULL;
    res_dgrp_t *grp = NULL;

    memset(&mi, 0, sizeof(module_input_t));

    if (nr_rows < 1 || nr_gens < 1 || nr_vars < 1) {
        fprintf(ERRSTREAM, "Empty module input.\n");
        return 0;
    }

    /* The grading is needed here, before the Gröbner basis, because the
     * adjoined component of generator j carries that generator's own
     * multidegree and there is nowhere else to compute it. */
    grp = res_dgrp_of_grading(grading, nr_vars);
    if (grp == NULL) {
        return 0;
    }
    const len_t glen = grp->len;

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
    rdeg   = (int32_t *)calloc(
            (unsigned long)nr2 * (unsigned long)glen, sizeof(int32_t));
    if (lens2 == NULL || exps2 == NULL || comps2 == NULL
            || cfs2 == NULL || rdeg == NULL) {
        goto cleanup;
    }

    /* Row degrees are passed through unnormalized; module_gb_from_input
     * applies the same normalization every other entry point gets, so both
     * report the same table. */
    if (row_degs != NULL) {
        memcpy(rdeg, row_degs,
                (unsigned long)nr_rows * (unsigned long)glen * sizeof(int32_t));
    }

    const int32_t *icf = (const int32_t *)cfs;
    int64_t off = 0, off2 = 0;
    for (i = 0; i < nr_gens; ++i) {
        lens2[i] = lens[i] + 1;
        /* the generator's multidegree, read off its first term; for the
         * graded input this machinery is about, every term agrees */
        module_term_multidegree(grp, rdeg + (int64_t)(nr_rows + i) * glen,
                exps + off * nr_vars,
                rdeg + (int64_t)(comps[off] - 1) * glen);

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

    if (module_gb_from_input(&mi, lens2, exps2, comps2, cfs2, rdeg,
                field_char, mon_order, NULL, strat, grading, stop, nr_rows, nr_vars, nr2,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                1 /* reduce */, info_level)) {
        goto cleanup;
    }
    if (!module_input_is_graded(&mi)) {
        goto cleanup;
    }

    bs_t * const gb = mi.gb;
    md_t * const st = mi.st;
    ht_t *bht       = gb->ht;

    /* pick out the elements supported in the adjoined components */
    sylen = (int32_t *)malloc(
            (unsigned long)(gb->lml > 0 ? gb->lml : 1) * sizeof(int32_t));
    syidx = (int64_t *)malloc(
            (unsigned long)(gb->lml > 0 ? gb->lml : 1) * sizeof(int64_t));
    if (sylen == NULL || syidx == NULL) {
        goto cleanup;
    }

    /* Keep only the first syz_rows rows of the syzygy matrix, i.e. only
     * the terms sitting in the first syz_rows adjoined components.  This
     * is a projection of each syzygy onto those coordinates, so a column
     * whose every term is dropped is the zero syzygy and goes with them.
     * 0 keeps everything, and so does any bound at or above nr_gens. */
    const int32_t syrows =
        (stop != NULL && stop->syz_rows > 0 && stop->syz_rows < nr_gens)
        ? stop->syz_rows : nr_gens;
    /* The round loop has already stopped on this, so what is left here is
     * to report no more than was asked for: the basis it stopped with can
     * hold one or two syzygies past the limit, the last round having
     * produced several at once. */
    const int32_t sylimit =
        (stop != NULL && stop->syz_limit > 0) ? stop->syz_limit : INT32_MAX;

    int32_t nsyz = 0;
    int64_t syterms = 0;
    for (i = 0; i < (int32_t)gb->lml && nsyz < sylimit; ++i) {
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
        int32_t klen = 0;
        for (t = 0; t < (int64_t)hm[LENGTH]; ++t) {
            const int32_t c = (int32_t)bht->ev[hm[OFFSET+t]][bht->cpos];
            if (c <= nr_rows) {
                pure = 0;
                break;
            }
            if (c - nr_rows <= syrows) {
                klen++;
            }
        }
        if (!pure || klen == 0) {
            continue;
        }
        syidx[nsyz] = bi;
        sylen[nsyz] = klen;
        syterms    += (int64_t)klen;
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
    /* the *normalized* heft degrees, which is what every other entry point
     * reports and what bht->cshift was built from */
    for (i = 0; i < nr2; ++i) {
        dg[cg++] = module_row_heft(grp, mi.rowmd + (int64_t)i * glen);
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
     * the generators they index.  The generator degree is read off the
     * lead term even when syrows drops it: the input is graded, so every
     * term of a syzygy carries the same heft degree, component shifts
     * included, and that is what the projection preserves. */
    for (i = 0; i < nsyz; ++i) {
        const hm_t *hm = gb->hm[syidx[i]];
        dg[cg++] = (int32_t)bht->hd[hm[OFFSET]].deg;
        dl[cl++] = sylen[i];
        for (t = 0; t < (int64_t)hm[LENGTH]; ++t) {
            const exp_t * const ev = bht->ev[hm[OFFSET+t]];
            if ((int32_t)ev[bht->cpos] - nr_rows > syrows) {
                continue;
            }
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
    res_dgrp_free(&grp);
    module_input_clear(&mi);

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
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const res_stop_t * const stop,
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
    int64_t nterms = 0;
    module_input_t mi;
    res_frame_t *f = NULL;
    res_diff_t *rd = NULL;

    memset(&mi, 0, sizeof(module_input_t));

    *nlevels = 0;
    *ranks   = NULL;
    *degs    = NULL;
    *dlen    = NULL;
    *dexp    = NULL;
    *dcomp   = NULL;
    *dcf     = NULL;

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
                field_char, mon_order, strat, grading, stop, nr_vars,
                nr_rows, nr_gens, ht_size, nr_threads, max_nr_pairs,
                la_option, info_level, max_level == 0 ? 2 : max_level);
    }
    if (syz_of != RES_SYZ_OF_GB) {
        fprintf(ERRSTREAM, "Unknown syzygy flavour %d.\n", syz_of);
        return 0;
    }

    /* As for the frame: the lead terms have to be the minimal generators
     * of the module of lead terms, which is what reducing gives. */
    if (module_gb_from_input(&mi, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, NULL, strat, grading, stop, 0, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                1 /* reduce */, info_level)) {
        return 0;
    }

    ht_t *bht = mi.gb->ht;

    if (!module_input_is_graded(&mi)) {
        module_input_clear(&mi);
        return 0;
    }

    f = res_frame_new(mi.grp, mi.st, max_level);
    if (f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }

    if (res_frame_init(f, mi.gb, bht, mi.rowmd)) {
        goto cleanup;
    }
    if (res_frame_complete(f) < 0 || res_frame_verify(f)) {
        goto cleanup;
    }

    rd = res_diff_new(f, mi.st->fc);
    if (rd == NULL) {
        fprintf(ERRSTREAM, "Could not set up the differential.\n");
        goto cleanup;
    }
    if (res_diff_init(rd, mi.gb, bht, mi.st) || res_diff_compute(rd)) {
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
    module_input_clear(&mi);

    return nterms;
}

/* --------------------------------------------------------------------- *
 *  Minimal Betti numbers and Hilbert information
 * --------------------------------------------------------------------- */

void free_module_betti_result_data(
        void (*freep) (void *),
        int32_t **betti,
        int32_t **hilbnum
        )
{
    if (betti != NULL && *betti != NULL) {
        (*freep)(*betti);
        *betti = NULL;
    }
    if (hilbnum != NULL && *hilbnum != NULL) {
        (*freep)(*hilbnum);
        *hilbnum = NULL;
    }
}

void free_module_mtable_data(
        void (*freep) (void *),
        res_mtable_t *mtab
        )
{
    if (mtab == NULL) {
        return;
    }
    if (mtab->degs != NULL) {
        (*freep)(mtab->degs);
        mtab->degs = NULL;
    }
    if (mtab->heft != NULL) {
        (*freep)(mtab->heft);
        mtab->heft = NULL;
    }
    if (mtab->betti != NULL) {
        (*freep)(mtab->betti);
        mtab->betti = NULL;
    }
    if (mtab->hilbnum != NULL) {
        (*freep)(mtab->hilbnum);
        mtab->hilbnum = NULL;
    }
    mtab->ndegs = 0;
}

int64_t export_module_betti(
        void *(*mallocp) (size_t),
        int32_t *nlevels,
        int32_t *maxdeg,
        int32_t *degshift,
        int32_t **betti,
        int32_t **hilbnum,
        int32_t *pdim,
        int32_t *reg,
        int32_t *dimension,
        int64_t *degree,
        res_mtable_t *mtab,
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t * const strat,
        const res_grading_t * const grading,
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t max_level,
        const int32_t minimal,
        const int32_t verify,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level
        )
{
    int64_t nelts   = 0;
    module_input_t mi;
    res_betti_t *bt = NULL;
    res_frame_t *f  = NULL;
    res_diff_t *rd  = NULL;
    int32_t *tab    = NULL;
    int32_t *num    = NULL;

    memset(&mi, 0, sizeof(module_input_t));

    *nlevels = 0;
    *maxdeg  = 0;
    if (degshift != NULL) {
        *degshift = 0;
    }
    if (betti != NULL) {
        *betti = NULL;
    }
    if (hilbnum != NULL) {
        *hilbnum = NULL;
    }
    if (mtab != NULL) {
        memset(mtab, 0, sizeof(res_mtable_t));
    }

    if (max_level < 0) {
        fprintf(ERRSTREAM, "A negative truncation level makes no sense.\n");
        return 0;
    }

    if (module_gb_from_input(&mi, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, NULL, strat, grading, NULL, 0, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                1 /* reduce */, info_level)) {
        return 0;
    }

    ht_t *bht = mi.gb->ht;

    if (!module_input_is_graded(&mi)) {
        module_input_clear(&mi);
        return 0;
    }

    /* module_gb_from_input normalizes the row degrees by subtracting the
     * multidegree of the lightest row.  The tables are arrays indexed by
     * degree and so have to stay on that scale, with degshift reporting the
     * offset; the scalars below do not, and are in the caller's own
     * degrees. */
    const len_t glen     = mi.grp->len;
    const int32_t dshift = module_row_heft(mi.grp, mi.degshift);
    if (degshift != NULL) {
        *degshift = dshift;
    }

    /* One level past what is reported: beta_{i,d} reads the rank of
     * d_{i+1} as well as of d_i, so the top level of a truncated table is
     * only right if the level above it was built.  Macaulay2's
     * minimalBetti does the same. */
    f = res_frame_new(mi.grp, mi.st, max_level > 0 ? max_level + 1 : 0);
    if (f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }

    if (res_frame_init(f, mi.gb, bht, mi.rowmd)) {
        goto cleanup;
    }
    nelts = res_frame_complete(f);
    if (nelts < 0 || res_frame_verify(f)) {
        nelts = 0;
        goto cleanup;
    }

    const int complete = res_frame_is_complete(f);
    if (!complete && (hilbnum != NULL || dimension != NULL
                || degree != NULL || pdim != NULL || reg != NULL)) {
        fprintf(ERRSTREAM, "Invariants of the whole module need the whole "
                "resolution -- the Hilbert numerator is an alternating sum "
                "over every level -- and the frame was truncated at level "
                "%d.\n", max_level);
        nelts = 0;
        goto cleanup;
    }

    bt = res_betti_new(f);
    if (bt == NULL) {
        fprintf(ERRSTREAM, "Could not tabulate the frame ranks.\n");
        nelts = 0;
        goto cleanup;
    }

    if (minimal) {
        rd = res_diff_new(f, mi.st->fc);
        if (rd == NULL) {
            fprintf(ERRSTREAM, "Could not set up the differential.\n");
            nelts = 0;
            goto cleanup;
        }
        if (res_diff_init(rd, mi.gb, bht, mi.st) || res_diff_compute(rd)) {
            nelts = 0;
            goto cleanup;
        }
        if (res_diff_verify(rd, verify != 0)) {
            fprintf(ERRSTREAM, "The computed differential is not a "
                    "complex.\n");
            nelts = 0;
            goto cleanup;
        }
        if (res_betti_minimalize(bt, rd)) {
            nelts = 0;
            goto cleanup;
        }
    }

    /* Every scalar is derived before a single byte is handed out, so that
     * mallocp -- which comes without a matching free -- is only ever
     * called once nothing downstream of it can fail. */
    /* The extra level exists only to make the top of the table right; it
     * is not part of what was asked for. */
    const len_t rlv  = (max_level > 0 && (len_t)max_level + 1 < bt->nlv)
        ? (len_t)max_level + 1 : bt->nlv;
    const size_t nd  = (size_t)bt->maxdeg + 1;
    const size_t tsz = (size_t)rlv * nd * sizeof(int32_t);
    int32_t xdim     = -1;
    int64_t xdeg     = 0;

    if (dimension != NULL || degree != NULL) {
        if (res_hilbert_invariants(bt->hilb, (len_t)nd, (len_t)nr_vars,
                    &xdim, &xdeg)) {
            nelts = 0;
            goto cleanup;
        }
    }

    if (betti != NULL) {
        tab = (int32_t *)(*mallocp)(tsz);
    }
    if (hilbnum != NULL) {
        num = (int32_t *)(*mallocp)(nd * sizeof(int32_t));
    }
    if ((betti != NULL && tab == NULL) || (hilbnum != NULL && num == NULL)) {
        fprintf(ERRSTREAM, "Could not allocate the Betti table.\n");
        nelts = 0;
        goto cleanup;
    }

    if (tab != NULL) {
        memcpy(tab, bt->betti, tsz);
        *betti = tab;
    }
    if (num != NULL) {
        memcpy(num, bt->hilb, nd * sizeof(int32_t));
        *hilbnum = num;
    }

    /* The multigraded table, whose rows are the same levels the heft table
     * reports and whose columns are the multidegrees that occur.  Allocated
     * last and all at once: mallocp has no matching free, so nothing here
     * may be abandoned, and by this point nothing downstream can fail. */
    if (mtab != NULL) {
        const size_t mtsz = (size_t)rlv * (size_t)bt->ndeg * sizeof(int32_t);

        mtab->nlevels = (int32_t)rlv;
        mtab->ndegs   = (int32_t)bt->ndeg;
        mtab->dlen    = (int32_t)glen;
        memcpy(mtab->degshift, mi.degshift,
                (unsigned long)glen * sizeof(int32_t));

        if (bt->ndeg > 0) {
            mtab->degs = (int32_t *)(*mallocp)(
                    (size_t)bt->ndeg * (size_t)glen * sizeof(int32_t));
            mtab->heft = (int32_t *)(*mallocp)(
                    (size_t)bt->ndeg * sizeof(int32_t));
            mtab->betti = (int32_t *)(*mallocp)(mtsz);
            /* The multigraded numerator is the same alternating sum over
             * every level the heft indexed one is, so a truncated frame
             * does not know it either; it is left NULL rather than filled
             * in wrongly.  Asking for hilbnum outright is refused above,
             * but mtab carries the two together. */
            mtab->hilbnum = complete ? (int32_t *)(*mallocp)(
                    (size_t)bt->ndeg * sizeof(int32_t)) : NULL;
            if (mtab->degs == NULL || mtab->heft == NULL
                    || mtab->betti == NULL
                    || (complete && mtab->hilbnum == NULL)) {
                fprintf(ERRSTREAM,
                        "Could not allocate the multigraded Betti table.\n");
                nelts = 0;
                goto cleanup;
            }
            memcpy(mtab->degs, bt->mdegs,
                    (size_t)bt->ndeg * (size_t)glen * sizeof(int32_t));
            memcpy(mtab->betti, bt->mbetti, mtsz);
            if (complete) {
                memcpy(mtab->hilbnum, bt->mhilb,
                        (size_t)bt->ndeg * sizeof(int32_t));
            }
            for (hl_t u = 0; u < bt->ndeg; ++u) {
                mtab->heft[u] = (int32_t)bt->mheft[u];
            }
        }
    }
    const int32_t xpd = res_betti_pdim(bt);
    if (pdim != NULL) {
        *pdim = xpd;
    }
    if (reg != NULL) {
        /* -1 is the zero module here, not a degree, so it is not shifted */
        *reg = xpd < 0 ? -1 : res_betti_reg(bt) + dshift;
    }
    if (dimension != NULL) {
        *dimension = xdim;
    }
    if (degree != NULL) {
        *degree = xdeg;
    }

    *nlevels = (int32_t)rlv;
    *maxdeg  = (int32_t)bt->maxdeg;

cleanup:
    res_betti_free(&bt);
    res_diff_free(&rd);
    res_frame_free(&f);
    module_input_clear(&mi);

    return nelts;
}

/* --------------------------------------------------------------------- *
 *  A resolution kept alive
 *
 *  Everything the one shot entry points do in a single call, split at the
 *  one place the work naturally divides: the frame, which is
 *  combinatorial and answers every question about the shape of the
 *  resolution, and the differential, which is where the field arithmetic
 *  is and which nobody should pay for until they ask for a matrix.
 *
 *  The Gröbner basis does not survive res_comp_new.  Level 1 of the
 *  differential is the only thing that ever reads it -- res_frame_init
 *  takes the lead terms and res_diff_init the coefficients -- so both run
 *  eagerly and the basis and its hash table are released before the handle
 *  is returned.  res_diff_init is O(the basis) and buys the caller the
 *  right to hold a large resolution without also holding the Gröbner basis
 *  it came from, which for the inputs this interface exists for is the
 *  bigger of the two.
 * --------------------------------------------------------------------- */

struct res_comp_t
{
    md_t        *st;   /* the meta data the whole computation ran under;
                        * nothing below reads it once the frame and level
                        * 1 are built, but it is what a later res_comp_*
                        * would resume from, and it is a few hundred
                        * bytes against a resolution                    */
    res_dgrp_t  *grp;
    res_frame_t *f;
    res_diff_t  *rd;   /* NULL when there is no level to differentiate */
    int32_t      degshift;  /* the heft of mdegshift                   */
    int32_t      mdegshift[RES_MTAB_MAXLEN];
    int32_t      nv;
};

void res_comp_free(
        res_comp_t **cp
        )
{
    res_comp_t *c = *cp;

    if (c == NULL) {
        return;
    }
    res_diff_free(&c->rd);
    res_frame_free(&c->f);
    res_dgrp_free(&c->grp);
    free(c->st);
    free(c);
    *cp = NULL;
}

res_comp_t *res_comp_new(
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs,
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t * const strat,
        const res_grading_t * const grading,
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
    module_input_t mi;
    res_comp_t *c = NULL;

    memset(&mi, 0, sizeof(module_input_t));

    if (max_level < 0) {
        fprintf(ERRSTREAM, "A negative truncation level makes no sense.\n");
        return NULL;
    }

    /* reduced, as for the frame: the lead terms have to be the minimal
     * generators of the module of lead terms */
    if (module_gb_from_input(&mi, lens, exps, comps, cfs, row_degs,
                field_char, mon_order, NULL, strat, grading, NULL, 0, nr_vars, nr_rows,
                nr_gens, ht_size, nr_threads, max_nr_pairs, la_option,
                1 /* reduce */, info_level)) {
        return NULL;
    }

    ht_t * const bht = mi.gb->ht;

    if (!module_input_is_graded(&mi)) {
        goto cleanup;
    }

    c = (res_comp_t *)calloc(1, sizeof(res_comp_t));
    if (c == NULL) {
        goto cleanup;
    }
    c->nv       = nr_vars;
    c->degshift = module_row_heft(mi.grp, mi.degshift);
    memcpy(c->mdegshift, mi.degshift,
            (unsigned long)mi.grp->len * sizeof(int32_t));

    /* The grading group outlives the Gröbner basis: the frame holds a
     * pointer to it and every multidegree the handle reports is a view
     * into a pool it owns, so it moves from mi into the handle here
     * rather than being freed with the rest of the input. */
    c->grp   = mi.grp;
    mi.grp   = NULL;

    c->f = res_frame_new(c->grp, mi.st, max_level);
    if (c->f == NULL) {
        fprintf(ERRSTREAM, "Could not set up the Schreyer frame.\n");
        goto cleanup;
    }

    if (res_frame_init(c->f, mi.gb, bht, mi.rowmd)) {
        goto cleanup;
    }
    if (res_frame_complete(c->f) < 0 || res_frame_verify(c->f)) {
        goto cleanup;
    }

    /* A frame with only level 0 has nothing to differentiate; that is the
     * zero submodule, and a legitimate answer rather than a failure. */
    if (c->f->nlv >= 2) {
        c->rd = res_diff_new(c->f, mi.st->fc);
        if (c->rd == NULL) {
            fprintf(ERRSTREAM, "Could not set up the differential.\n");
            goto cleanup;
        }
        if (res_diff_init(c->rd, mi.gb, bht, mi.st)) {
            goto cleanup;
        }
    }

    c->st  = mi.st;
    mi.st  = NULL;

    module_input_clear(&mi);

    return c;

cleanup:
    res_comp_free(&c);
    module_input_clear(&mi);

    return NULL;
}

int32_t res_comp_nlevels(
        const res_comp_t * const c
        )
{
    if (c == NULL || c->f == NULL || c->f->bad) {
        return 0;
    }

    return (int32_t)c->f->nlv;
}

int32_t res_comp_degshift(
        const res_comp_t * const c
        )
{
    return c == NULL ? 0 : c->degshift;
}

int res_comp_is_complete(
        const res_comp_t * const c
        )
{
    if (c == NULL || c->f == NULL) {
        return 0;
    }

    return res_frame_is_complete(c->f);
}

int32_t res_comp_rank(
        const res_comp_t * const c,
        const int32_t level
        )
{
    const int32_t nlv = res_comp_nlevels(c);

    if (level < 0 || level >= nlv) {
        return -1;
    }

    return (int32_t)c->f->lv[level].ld;
}

int res_comp_degrees(
        const res_comp_t * const c,
        const int32_t level,
        int32_t *degs
        )
{
    len_t k;
    const int32_t rk = res_comp_rank(c, level);

    if (rk < 0 || degs == NULL) {
        return 1;
    }
    for (k = 0; k < (len_t)rk; ++k) {
        degs[k] = (int32_t)c->f->lv[level].elts[k].hdeg;
    }

    return 0;
}

int32_t res_comp_glen(
        const res_comp_t * const c
        )
{
    if (c == NULL || c->grp == NULL) {
        return 0;
    }

    return (int32_t)c->grp->len;
}

int res_comp_multidegrees(
        const res_comp_t * const c,
        const int32_t level,
        int32_t *mdegs
        )
{
    len_t k;
    const int32_t rk = res_comp_rank(c, level);

    if (rk < 0 || mdegs == NULL) {
        return 1;
    }

    const res_level_t * const lv = c->f->lv + level;
    const len_t glen = c->grp->len;

    for (k = 0; k < (len_t)rk; ++k) {
        const res_deg_t d = res_dpool_at(lv->degs, lv->elts[k].mdeg);
        memcpy(mdegs + (int64_t)k * glen, d.e,
                (unsigned long)glen * sizeof(int32_t));
    }

    return 0;
}

int res_comp_multidegshift(
        const res_comp_t * const c,
        int32_t *shift
        )
{
    if (c == NULL || c->grp == NULL || shift == NULL) {
        return 1;
    }
    memcpy(shift, c->mdegshift,
            (unsigned long)c->grp->len * sizeof(int32_t));

    return 0;
}

void free_module_differential_data(
        void (*freep) (void *),
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        )
{
    if (dlen != NULL && *dlen != NULL) {
        (*freep)(*dlen);
        *dlen = NULL;
    }
    if (dexp != NULL && *dexp != NULL) {
        (*freep)(*dexp);
        *dexp = NULL;
    }
    if (dcomp != NULL && *dcomp != NULL) {
        (*freep)(*dcomp);
        *dcomp = NULL;
    }
    if (dcf != NULL && *dcf != NULL) {
        (*freep)(*dcf);
        *dcf = NULL;
    }
}

int64_t res_comp_differential(
        void *(*mallocp) (size_t),
        res_comp_t *c,
        const int32_t level,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        )
{
    len_t k;
    int64_t t, nterms = 0;

    if (dlen == NULL || dexp == NULL || dcomp == NULL || dcf == NULL) {
        return 0;
    }
    *dlen  = NULL;
    *dexp  = NULL;
    *dcomp = NULL;
    *dcf   = NULL;

    const int32_t rk = res_comp_rank(c, level);
    if (rk < 0 || level < 1) {
        fprintf(ERRSTREAM, "There is no differential at level %d; the "
                "resolution has levels 0 to %d.\n",
                level, res_comp_nlevels(c) - 1);
        return 0;
    }
    if (c->rd == NULL || c->rd->bad) {
        fprintf(ERRSTREAM, "The differential is unusable.\n");
        return 0;
    }

    /* the prefix below this level, once; already computed levels cost
     * nothing */
    if (res_diff_compute_thru(c->rd, (len_t)level)) {
        fprintf(ERRSTREAM, "Could not compute the differential up to "
                "level %d.\n", level);
        return 0;
    }
    if (res_diff_verify(c->rd, 0)) {
        fprintf(ERRSTREAM, "The computed differential is not a complex.\n");
        return 0;
    }

    const res_frame_t * const f = c->f;
    const len_t nv = (len_t)c->nv;

    for (k = 0; k < (len_t)rk; ++k) {
        nterms += (int64_t)c->rd->d[level][k].len;
    }
    if (nterms > (int64_t)INT32_MAX / (nv > 0 ? (int64_t)nv : 1)) {
        fprintf(ERRSTREAM,
                "The differential is too large to store in flat arrays.\n");
        return 0;
    }

    /* nothing can fail after this point, mallocp coming without a
     * matching free */
    int32_t *dl = (int32_t *)(*mallocp)(
            (unsigned long)(rk > 0 ? rk : 1) * sizeof(int32_t));
    int32_t *de = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1)
            * (unsigned long)nv * sizeof(int32_t));
    int32_t *dc = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1) * sizeof(int32_t));
    int32_t *cf = (int32_t *)(*mallocp)(
            (unsigned long)(nterms > 0 ? nterms : 1) * sizeof(int32_t));

    int64_t ce = 0, cc = 0;
    for (k = 0; k < (len_t)rk; ++k) {
        const res_dpoly_t * const p = c->rd->d[level] + k;
        dl[k] = (int32_t)p->len;
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

    *dlen  = dl;
    *dexp  = de;
    *dcomp = dc;
    *dcf   = (void *)cf;

    return nterms;
}
