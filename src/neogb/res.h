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

/* Syzygies and free resolutions of graded modules.
 *
 * This header declares the grading-group layer and the module-order
 * vocabulary shared by the resolution engine.  See res_grading.c for the
 * arithmetic and res_frame.c for the Schreyer frame. */

#ifndef GB_RES_H
#define GB_RES_H

#include "data.h"
#include "order.h"

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------- *
 *  Order strategies
 *
 *  There is more than one sensible order on a free module, they are not
 *  equivalent -- the Gröbner basis, the frame ranks and every matrix above
 *  them change with the choice -- and which one is best is an empirical
 *  question this engine exists to answer.  So the choice is a *parameter*,
 *  carried as one struct from the C entry points down to the innermost
 *  comparison, rather than a constant compiled into res_diff.c.
 *
 *  Three independent axes:
 *
 *    base  how two monomials of the ambient free module R^ncomp compare.
 *          res_mord_t, in data.h next to ht_t, since the Gröbner basis in
 *          levels 0 and 1 is computed under it by msolve's own F4 and the
 *          hash table is what carries it there.
 *
 *    pos   which component index counts as the larger one.  msolve has
 *          always taken the smaller index to be larger (RES_POS_DOWN);
 *          Macaulay2's Position => Up is the other one.  Note Macaulay2's
 *          interface currently obtains Up by numbering components in
 *          reverse as it marshals, so setting RES_POS_UP there as well
 *          would compose to Down.
 *
 *    lift  how a level above 0 is ordered.  Only RES_LIFT_SCHREYER exists
 *          and it is not really optional: Schreyer's theorem determining
 *          the lead terms of the syzygies from the level below is exactly
 *          what makes the frame combinatorial, and it determines the order
 *          too.  The axis is here because it is a genuine choice
 *          mathematically -- one could resolve level by level in any
 *          module order -- and naming it keeps that visible.
 *
 *  Extending this is meant to be adding a field, not changing a signature:
 *  the entry points take a res_strat_t pointer rather than a bare int for
 *  that reason.  A user supplied weight matrix is the next field.
 * --------------------------------------------------------------------- */

typedef enum {
    RES_POS_DOWN = 0, /* smaller component index is the larger monomial */
    RES_POS_UP   = 1  /* larger component index is the larger monomial  */
} res_mpos_t;

typedef enum {
    RES_LIFT_SCHREYER = 0 /* levels above 0 carry the order level 0 induces */
} res_mlift_t;

typedef struct res_strat_t res_strat_t;
struct res_strat_t
{
    int32_t base; /* res_mord_t  */
    int32_t pos;  /* res_mpos_t  */
    int32_t lift; /* res_mlift_t */
};

/* Position over term, smaller component larger, Schreyer above level 0.
 * This is what every entry point does when handed a NULL strategy, and it
 * is what the engine did before strategies existed. */
res_strat_t res_strat_default(
        void
        );

/* The strategy a bare res_mord_t used to mean, for callers that have one
 * and want the defaults for the rest. */
res_strat_t res_strat_of_order(
        const int32_t module_order
        );

/* Whether the strategy is one the engine can run.  for_resolution asks
 * the stricter question, since the frame and the differential support
 * fewer bases than a plain module Gröbner basis does.  Returns 0 when it
 * is usable and reports the reason on stderr otherwise. */
int res_strat_check(
        const res_strat_t * const s,
        const int for_resolution
        );

/* A short stable name, "pot-down-schreyer" and so on, for test output and
 * benchmark tables.  Points into static storage. */
const char *res_strat_name(
        const res_strat_t * const s
        );

/* --------------------------------------------------------------------- *
 *  Stopping conditions
 *
 *  Macaulay2's gb, syz and res ask for less than a complete answer in
 *  several ways; these are the three the engine can honour today.  A NULL
 *  res_stop_t asks for everything, which is what every entry point did
 *  before this existed, and each field switches itself off with a value
 *  that cannot be meant literally.
 *
 *    max_degree  A ceiling on the degree of the Gröbner basis, as one
 *                multidegree -- res_grading_len(grading) entries, so a
 *                single integer under the standard grading -- on the
 *                caller's own scale, not the internally normalized one.
 *                NULL means no ceiling.
 *
 *                A multidegree is *not* what the computation schedules
 *                by.  msolve selects S-pairs by the single integer
 *                heft . deg, which is what makes the degree by degree
 *                schedule terminate, so a multidegree ceiling is honoured
 *                as its heft: everything of multidegree at most
 *                max_degree is computed, and so is everything else of no
 *                greater heft.  That is a coarsening in the safe
 *                direction, it is what Macaulay2 does to its own
 *                DegreeLimit before the engine ever sees it
 *                (degreeToHeft in m2/gb.m2), and under the standard
 *                grading the two coincide.
 *
 *                A truncated basis is a Gröbner basis of nothing in
 *                particular: it generates the same module only in degrees
 *                up to the ceiling.  Callers that care must say so.
 *
 *                A ceiling needs homogeneous input, and inhomogeneous
 *                input with one is refused.  What the round loop
 *                schedules by there is a sugar degree, which can fall --
 *                md->min_deg_in_first_deg_fall exists to notice it -- so
 *                a ceiling on it would not mean the basis is complete
 *                through that degree, which is the only thing a caller
 *                could want it to mean.
 *
 *    syz_limit   Stop once this many syzygies of the input generators are
 *                known, 0 or less meaning no limit.  Only the syzygy
 *                entry points look at it -- a plain Gröbner basis
 *                computes no syzygies to count.
 *
 *    syz_rows    Keep only the first syz_rows rows of the syzygy matrix,
 *                0 or less meaning all of them.  A syzygy of nr_gens
 *                generators is a vector in R^nr_gens; this projects it
 *                onto its first syz_rows coordinates and drops the
 *                columns that become zero, which is what Macaulay2 wants
 *                when it computes a Gröbner basis of generators together
 *                with relations and cares only about the coefficients on
 *                the generators.  It is an output filter and saves no
 *                work: pruning the computation itself needs the discarded
 *                components ordered last, which is the elimination block
 *                question the module orders do not yet answer.
 *
 *                A projected syzygy matrix is a submatrix, so it is no
 *                longer a complex: d_1 o d_2 is the composite of the
 *                presentation with only some of the relations and has no
 *                reason to vanish.  The d o d = 0 verification is
 *                therefore not run on this path, and callers must not
 *                expect it to hold.
 *
 *  None of these makes a computation resumable.  Asking again for a
 *  larger ceiling recomputes from the input.
 * --------------------------------------------------------------------- */

typedef struct res_stop_t res_stop_t;
struct res_stop_t
{
    const int32_t *max_degree; /* res_grading_len entries, or NULL */
    int32_t        syz_limit;  /* <= 0 means no limit  */
    int32_t        syz_rows;   /* <= 0 means all rows  */
};

/* Everything, i.e. what a NULL res_stop_t means. */
res_stop_t res_stop_none(
        void
        );

/* --------------------------------------------------------------------- *
 *  Degrees
 *
 *  Modules are graded by a finitely generated abelian group
 *
 *      A  =  Z^r (+) Z/t_1 (+) ... (+) Z/t_nt
 *
 *  which covers the standard grading (r = 1, nt = 0), multigradings by
 *  Z^r, the class group Cl(X) of a toric variety, and -- once an external
 *  group-arithmetic backend is attached to the vtable below -- character
 *  groups for equivariant resolutions.
 *
 *  A degree is a *struct*, never a bare array: its layout is private to
 *  the owning res_dgrp_t so that it can later grow a cached hash, a
 *  torsion handle, or a pointer into an external library without breaking
 *  any consumer.
 * --------------------------------------------------------------------- */

typedef struct res_deg_t res_deg_t;
struct res_deg_t
{
    int32_t *e;  /* [0,r) free part, [r,r+nt) torsion residues;
                  * interpretation belongs to the owning res_dgrp_t */
};

typedef struct res_dgrp_t res_dgrp_t;
struct res_dgrp_t
{
    len_t    r;      /* rank of the free part                            */
    len_t    nt;     /* number of torsion factors                        */
    len_t    len;    /* int32 slots per degree, = r + nt                 */
    len_t    nv;     /* number of ring variables                         */
    int32_t *tord;   /* orders of the torsion factors, length nt         */
    int32_t *heft;   /* length r; heft . deg(x_j) > 0 for every j        */
    int32_t *dmat;   /* len x nv, column major: column j is deg(x_j)     */
    deg_t   *vhdeg;  /* length nv, precomputed heft . deg(x_j)           */
    int      simple; /* 1 iff r == 1 && nt == 0 (the standard grading)   */

    /* Arithmetic vtable.  This is the seam an external group-arithmetic
     * library plugs into; the defaults in res_grading.c implement
     * Z^r (+) torsion directly. */
    void  (*add)(const res_dgrp_t *g, res_deg_t d,
                 const res_deg_t a, const res_deg_t b);
    void  (*sub)(const res_dgrp_t *g, res_deg_t d,
                 const res_deg_t a, const res_deg_t b);
    int   (*cmp)(const res_dgrp_t *g, const res_deg_t a, const res_deg_t b);
    hl_t  (*hash)(const res_dgrp_t *g, const res_deg_t a);
    deg_t (*heft_of)(const res_dgrp_t *g, const res_deg_t a);
};

/* Contiguous storage for many degrees of one group.  Degrees are handed
 * out as views into the pool, so a res_deg_t stays one pointer wide. */
typedef struct res_dpool_t res_dpool_t;
struct res_dpool_t
{
    const res_dgrp_t *grp;
    int32_t *data;
    hl_t     ld;   /* number of degrees in use   */
    hl_t     sz;   /* number of degrees allocated */
};

/* --- construction ---------------------------------------------------- */

/* degs is len x nv in column-major order: entries [j*len, (j+1)*len) give
 * deg(x_j), free part first then torsion residues.  heft has length r and
 * must satisfy heft . deg(x_j) > 0 for every variable j; pass NULL to use
 * the all-ones vector.  Returns NULL and reports on stderr if the heft is
 * not strictly positive on some variable. */
res_dgrp_t *res_dgrp_new(
        const len_t r,
        const len_t nt,
        const int32_t *tord,
        const len_t nv,
        const int32_t *degs,
        const int32_t *heft
        );

/* The standard grading: A = Z, every variable in degree 1. */
res_dgrp_t *res_dgrp_new_standard(
        const len_t nv
        );

void res_dgrp_free(
        res_dgrp_t **gp
        );

/* --- how a caller describes a grading ------------------------------- *
 *
 *  res_dgrp_t is the engine's internal object, owning storage and a
 *  vtable.  This is the flat description a caller hands in: no
 *  allocation, no ownership, valid only for the duration of the call.
 *  Every C entry point takes one, and NULL means the standard grading --
 *  A = Z with every variable in degree 1 -- which is what the engine did
 *  before gradings were a parameter.
 *
 *  Everything a grading needs is here, so extending it is adding a field:
 *  a name for an external group backend, a presentation of a non-abelian
 *  group, a cached heft.  The entry points take a pointer for that reason,
 *  exactly as they do for res_strat_t.
 * --------------------------------------------------------------------- */

typedef struct res_grading_t res_grading_t;
struct res_grading_t
{
    int32_t r;           /* rank of the free part of A                   */
    int32_t nt;          /* number of torsion factors                    */
    const int32_t *tord; /* nt orders, each >= 2; NULL iff nt == 0       */
    const int32_t *degs; /* (r+nt) x nv, column major: entries
                          * [j*(r+nt), (j+1)*(r+nt)) are deg(x_j), free
                          * part first then torsion residues             */
    const int32_t *heft; /* r entries with heft . deg(x_j) > 0 for every
                          * j; NULL means the all ones vector            */
};

/* The number of int32 slots one degree of this grading occupies, which is
 * how many entries per row the row_degs argument of every entry point
 * carries and how many per generator res_comp_multidegrees writes.  A NULL
 * grading is the standard one and reports 1. */
int32_t res_grading_len(
        const res_grading_t * const g
        );

/* Builds the engine's degree group from the caller's description.
 * grading == NULL asks for res_dgrp_new_standard(nv).  Returns NULL and
 * reports on stderr if the description is inconsistent -- a torsion order
 * below two, or a heft that is not strictly positive on some variable,
 * which is what would make the degree by degree schedule non-terminating. */
res_dgrp_t *res_dgrp_of_grading(
        const res_grading_t * const grading,
        const len_t nv
        );

res_dpool_t *res_dpool_new(
        const res_dgrp_t *grp,
        const hl_t sz
        );

void res_dpool_free(
        res_dpool_t **pp
        );

/* View of degree i; the pool must already hold at least i+1 degrees. */
res_deg_t res_dpool_at(
        const res_dpool_t *p,
        const hl_t i
        );

/* Appends a zero degree and returns a view of it, growing if needed.
 * The index of the new degree is written to *idx when idx != NULL.
 * NOTE: any previously handed out res_deg_t is invalidated by a growth,
 * so hold indices rather than views across calls. */
res_deg_t res_dpool_push(
        res_dpool_t *p,
        hl_t *idx
        );

/* --- multidegree buckets --------------------------------------------- *
 *
 *  Everything graded that this engine reports -- frame ranks, minimal
 *  Betti numbers, the Hilbert numerator -- is a table indexed by degree.
 *  Under the standard grading a degree is a small nonnegative integer and
 *  the table is a dense array; under a grading by Z^r (+) torsion it is
 *  not, and the degrees that actually occur are a sparse, unpredictable
 *  subset of a lattice that may well be negative in places.
 *
 *  So the general table is indexed by *bucket*: an index into the set of
 *  distinct degrees that occur, discovered as the frame is walked.  This
 *  is a hash set over the group's own hash and comparison, which is why
 *  those are in the vtable -- an external group backend gets bucketing for
 *  free by supplying them.
 *
 *  res_dbkt_sort puts the buckets in the group's own order once, after
 *  which bucket indices are stable and every table built on them can be
 *  reported directly.  Bucket 0 is not special.
 * --------------------------------------------------------------------- */

typedef struct res_dbkt_t res_dbkt_t;
struct res_dbkt_t
{
    const res_dgrp_t *grp;
    res_dpool_t      *pool; /* the distinct degrees, in insertion order  */
    hl_t             *map;  /* open addressed, entries are 1 + bucket    */
    hl_t              msz;  /* a power of two                            */
    hl_t              ld;   /* number of distinct degrees                */
};

res_dbkt_t *res_dbkt_new(
        const res_dgrp_t *grp,
        const hl_t sz
        );

void res_dbkt_free(
        res_dbkt_t **bp
        );

/* The bucket of a, or -1 when a has not been inserted. */
hl_t res_dbkt_find(
        const res_dbkt_t * const b,
        const res_deg_t a
        );

/* The bucket of a, inserting it if it is new; -1 on allocation failure. */
hl_t res_dbkt_insert(
        res_dbkt_t *b,
        const res_deg_t a
        );

/* A view of the degree in bucket i.  Invalidated by a later insert, so
 * hold buckets rather than views. */
res_deg_t res_dbkt_at(
        const res_dbkt_t * const b,
        const hl_t i
        );

/* Reorders the buckets into the group's own order, writing the new bucket
 * of every old bucket into perm, which must have b->ld entries and which
 * the caller uses to relabel any table it has already built.  Returns 0 on
 * success.  No further insert may follow a sort. */
int res_dbkt_sort(
        res_dbkt_t *b,
        hl_t *perm
        );

/* --- arithmetic ------------------------------------------------------ */

void res_deg_zero(
        const res_dgrp_t *g,
        res_deg_t d
        );

void res_deg_set(
        const res_dgrp_t *g,
        res_deg_t d,
        const res_deg_t a
        );

/* d <- multidegree of the monomial with exponent vector exps[0..nv), i.e.
 * the column sum of dmat weighted by the exponents.  exps points at the
 * *variable* slots, not at a raw hash-table exponent vector; callers pass
 * ev + 1 for a non-block order. */
void res_deg_of_exponents(
        const res_dgrp_t *g,
        res_deg_t d,
        const exp_t *exps
        );

/* heft . (dmat . exps), the single integer that schedules and orders the
 * computation and that gets stored in ev[DEG]. */
deg_t res_heft_of_exponents(
        const res_dgrp_t *g,
        const exp_t *exps
        );

/* Install the grading's variable weights on a hash table, or verify that
 * a weighted block order already put the same ones there.  Call before
 * the secondary hash tables are built, since they share vwt by pointer.
 * Returns 0 on success. */
int res_install_weights(
        ht_t *ht,
        const res_dgrp_t * const grp
        );

/* --------------------------------------------------------------------- *
 *  Schreyer frames
 *
 *  The frame is the skeleton of a free resolution: the lead terms of
 *  every differential, and nothing else.  Schreyer's theorem determines
 *  the lead terms at level i+1 from those at level i, so the whole thing
 *  is combinatorics -- not one field operation is performed here.
 *
 *  Level 0 is the ambient free module R^ncomp, level 1 is the Gröbner
 *  basis of the input submodule, and level i+1 is obtained from level i
 *  as follows.  Write the elements of level i as m_k e_{up(k)}, with m_k a
 *  ring monomial.  Schreyer's theorem says the lead terms of the syzygies
 *  are the monomials (lcm(m_k,m_l)/m_k) e_k for l != k with up(l) == up(k),
 *  the larger of the two indices carrying the lead term.  For fixed k the
 *  ring monomials so obtained generate the ideal quotient
 *
 *      J_k = ( <m_l : up(l) == up(k), l < k> : m_k ),
 *
 *  and the level i+1 elements sitting over k are exactly its minimal
 *  generators.  "l < k" is the storage index order, which is also the tie
 *  break of the Schreyer order, so the two are consistent by construction.
 *
 *  Storage order within a level is (up ascending, own monomial ascending).
 *  The elements over a given parent are therefore contiguous, which is what
 *  makes the ideal quotient above a scan of one block, and within a block
 *  the order is the Schreyer order, since all of its elements share the
 *  same lift below.
 *
 *  Only the within-block order is visible in the result: elements in
 *  different blocks never pair up, so the number and the degrees of the
 *  frame elements do not depend on how the blocks themselves are arranged.
 *  The direction of the within-block order, on the other hand, changes the
 *  frame -- every direction gives a valid one, but ascending gives the
 *  small one, and is what Macaulay2 uses; see res_frame.c.
 * --------------------------------------------------------------------- */

typedef struct res_felt_t res_felt_t;
struct res_felt_t
{
    hm_t    mono;   /* own ring monomial, in the frame hash table       */
    hm_t    total;  /* mono lifted all the way to level 0, same table   */
    int32_t up;     /* index of the parent at the previous level,
                     * -1 at level 0                                    */
    int32_t root;   /* 1-based component of R^ncomp this element lifts
                     * to; constant along the chain of parents          */
    deg_t   hdeg;   /* heft degree, component shifts included           */
    hl_t    mdeg;   /* index of the multidegree in the level's pool     */
};

/* total is what makes a Schreyer comparison at any level cost one ring
 * monomial comparison: comparing u*e_a with v*e_b is comparing the ring
 * monomials u*total(a) and v*total(b), and, if those agree and the two
 * elements lift to the same component of R^ncomp, the indices a and b,
 * larger index being the larger monomial.  Without it every comparison
 * would walk the chain of parents. */

typedef struct res_level_t res_level_t;
struct res_level_t
{
    res_felt_t  *elts;
    len_t        ld;    /* number of elements            */
    len_t        sz;    /* number of elements allocated  */
    res_dpool_t *degs;  /* one multidegree per element   */
};

typedef struct res_frame_t res_frame_t;
struct res_frame_t
{
    res_level_t      *lv;     /* lv[0 .. lvsz)                         */
    len_t             nlv;    /* number of levels built, so lv[0,nlv)  */
    len_t             lvsz;   /* levels allocated                      */
    len_t             maxlv;  /* highest level index that may be built,
                               * or 0 for no ceiling at all            */
    ht_t             *ht;     /* ring monomials of the frame; a plain
                               * table, the component lives in up      */
    const res_dgrp_t *grp;
    len_t             nv;
    len_t             ncomp;  /* rank of the ambient free module       */
    res_strat_t       strat;  /* the order strategy the whole tower runs
                               * in: strat.base and strat.pos are what the
                               * Gröbner basis in levels 0 and 1 was
                               * computed under, cross checked against the
                               * module hash table in res_frame_init, and
                               * strat.lift is how every level above 0 is
                               * derived from it                        */
    int32_t          *gbmap;  /* for each level 1 element, the index in
                               * the Gröbner basis it was read off, so
                               * that res_diff.c can find its
                               * coefficients again after the frame has
                               * sorted the lead terms; length lv[1].ld */
    int               bad;    /* set once a level could not be built;
                               * a truncated frame must never be read
                               * as a complete one                     */
};

/* md supplies the number of variables, the ring order and the initial
 * hash table size; a private ring hash table is built from it.  Levels 0
 * to maxlevel may be built; pass maxlevel <= 0 for no ceiling, which is
 * the only way to get a complete frame.
 *
 * There is no a priori level that is "always enough".  The frame is a
 * nonminimal resolution, so Hilbert's syzygy theorem does not apply to
 * it: (z, y^2, x^2 y, x^3) in three variables has the frame 1,4,6,4,1 and
 * so reaches level four.  Any ceiling is therefore a truncation, and
 * res_frame_is_complete is what tells the two apart. */
res_frame_t *res_frame_new(
        const res_dgrp_t *grp,
        const md_t *md,
        const int32_t maxlevel
        );

void res_frame_free(
        res_frame_t **fp
        );

/* Fills levels 0 and 1 from a module Gröbner basis and the module hash
 * table it lives in.  row_mdegs holds the multidegrees of the generators
 * of R^ncomp, ncomp * grp->len int32 slots with component i at offset
 * i*grp->len; NULL is read as all zero.  Heft degrees come from
 * bht->cshift, which is what the Gröbner basis itself was graded by.
 * Returns 0 on success. */
int res_frame_init(
        res_frame_t *f,
        const bs_t * const gb,
        const ht_t * const bht,
        const int32_t *row_mdegs
        );

/* Builds one more level from the last one; returns the number of elements
 * added, which is 0 exactly when the frame is complete. */
len_t res_frame_next_level(
        res_frame_t *f
        );

/* Repeats res_frame_next_level until a level comes out empty or maxlv is
 * reached.  Returns the total number of frame elements, or -1 if a level
 * could not be built, in which case the frame is truncated at an
 * arbitrary point and its ranks mean nothing. */
int64_t res_frame_complete(
        res_frame_t *f
        );

/* Whether the frame ended on its own -- a level came out empty -- rather
 * than being cut off at maxlv.  Alternating sums over a cut off frame are
 * missing their tail, so anything Hilbert refuses to run on one. */
int res_frame_is_complete(
        const res_frame_t * const f
        );

/* Checks the invariants the rest of the engine reads the frame through:
 * every element's lifted monomial really is its own times its parent's,
 * it lifts to the same component of R^ncomp as its parent, its heft
 * degree is its parent's plus its own, and each level is grouped by
 * parent.  Returns the number of violations, so 0 means the frame is
 * sound.  Cheap -- the lift is a hash lookup of a monomial the frame
 * already inserted -- so it runs on every frame rather than under a
 * debug flag. */
int res_frame_verify(
        const res_frame_t * const f
        );

deg_t res_frame_max_hdeg(
        const res_frame_t * const f
        );

/* Writes the frame ranks into tab, which must have f->nlv * (maxdeg+1)
 * int32 entries: tab[i*(maxdeg+1) + d] is the number of level i elements
 * of heft degree d.  Returns the total number of frame elements. */
int64_t res_frame_betti(
        const res_frame_t * const f,
        int32_t *tab,
        const deg_t maxdeg
        );

/* --------------------------------------------------------------------- *
 *  The nonminimal differential
 *
 *  The frame gives the lead term of every entry of the resolution; this
 *  layer fills in the rest.  Write D_i(k) for the image in F_{i-1} of the
 *  k-th basis element of F_i.  Level 1 is the module Gröbner basis, and
 *  for i >= 2 the defining property D_{i-1}(D_i(k)) = 0 reads
 *
 *      m_k * D_{i-1}(up(k))  =  sum_j c_j u_j D_{i-1}(q_j),
 *
 *  so D_i(k) = m_k e_{up(k)} - sum_j c_j u_j e_{q_j}.  The left hand side
 *  is a multiple of a generator, hence lies in the module the D_{i-1}(q)
 *  generate, and Schreyer's theorem says those are a Gröbner basis of it;
 *  the reduction therefore always reaches zero, and the coefficients it
 *  used *are* the differential.  Not one S-pair is formed: the frame is
 *  the entire schedule.
 *
 *  The one rule that makes the lead term come out right is to reduce the
 *  leading monomial of the left hand side by the *smallest indexed*
 *  admissible reducer.  The frame guarantees one with index < up(k)
 *  exists -- that is exactly what minimality of m_k in the colon ideal
 *  says -- and every tail term is then automatically below m_k e_{up(k)}
 *  in the Schreyer order.
 *
 *  Scheduling is by degree ascending and, within a degree, by level
 *  ascending.  Level i in degree d reduces against level i-1 in degrees
 *  *up to and including* d: a reducer used with multiplier 1 has the same
 *  degree as the row it reduces.  That is why the slanted degree d - i is
 *  not by itself a valid schedule -- it would put level i-1 in degree d,
 *  whose slanted degree is one larger, after the row that needs it.
 * --------------------------------------------------------------------- */

/* Accumulator bound.  Every summand mult*val is below fc^2 <= 2^62, and a
 * value is reduced as soon as it reaches 2^62, so the addition itself
 * always happens below 2^63 and cannot wrap.
 *
 * Used by both res_diff.c (delayed reduction of a differential row) and
 * res_betti.c (delayed reduction of a rank-extraction row).  gb.c's unity
 * build includes res_diff.c before res_betti.c, so this has to live here
 * rather than in either .c file. */
#define RES_ACC_LIMIT (((uint64_t)1) << 62)

/* One column of a differential: an element of the free module one level
 * down, as a sorted list of terms cf * mon * e_pos with mon a ring
 * monomial in the frame's hash table and pos an index into that level. */
typedef struct res_dpoly_t res_dpoly_t;
struct res_dpoly_t
{
    len_t     len;
    hm_t     *mon;
    int32_t  *pos;
    uint32_t *cf;
};

typedef struct res_diff_t res_diff_t;
struct res_diff_t
{
    res_frame_t  *f;    /* not owned                                    */
    uint32_t      fc;   /* field characteristic                         */
    res_dpoly_t **d;    /* d[i][k] = D_i(k) for 1 <= i < f->nlv         */
    len_t         nlv;
    len_t         thru;  /* levels 1 to thru are filled in; 0 before
                          * res_diff_init has run                       */
    len_t         nthrds; /* set from md->nthrds in res_diff_init, clamped
                           * to at least 1; rows within one res_diff_block
                           * are independent once the pivot set is fixed,
                           * which is the only parallelism available -- the
                           * degree-then-level schedule that calls it is
                           * genuinely sequential                        */
    int           bad;
};

/* The differential of a completed frame, in the order the frame's own
 * strategy induces.  Both bases are supported: under position over term
 * components are compared first, so the ring comparison happens inside a
 * single component where the degree shifts cancel and the frame's plain
 * ring table is enough; under term over position degrees come first and
 * the shift is added back from level 0, where res_frame_init stored it as
 * the heft degree of each generator of R^ncomp. */
res_diff_t *res_diff_new(
        res_frame_t *f,
        const uint32_t fc
        );

void res_diff_free(
        res_diff_t **dp
        );

/* Level 1: the Gröbner basis itself, re-expressed over the frame's ring
 * hash table and made monic.  gb and bht are the same ones res_frame_init
 * was given. */
int res_diff_init(
        res_diff_t *rd,
        const bs_t * const gb,
        const ht_t * const bht,
        const md_t * const md
        );

/* Levels 2 upwards.  Returns 0 on success. */
int res_diff_compute(
        res_diff_t *rd
        );

/* Levels rd->thru+1 up to and including maxlev, leaving the rest empty.
 * Returns 0 on success, and 0 immediately if there is nothing to do.
 *
 * A prefix is the only truncation that makes sense: the block at level i
 * in degree d reduces against level i-1 in degrees up to and including d,
 * so D_i cannot be had without all of D_2 ... D_{i-1}.  Within that
 * constraint the work is exactly the same as computing the levels in one
 * go -- the schedule is degree ascending and then level ascending, and
 * restricting it to a range of levels visits the same blocks in the same
 * relative order, each one seeing the same already reduced data below it.
 * So res_diff_compute is res_diff_compute_thru(rd, rd->nlv - 1), and
 * calling it in several steps costs no more than calling it once. */
int res_diff_compute_thru(
        res_diff_t *rd,
        const len_t maxlev
        );

/* Checks the differential against the frame: every column leads with the
 * frame's monomial, with coefficient one, over the frame's parent.  That
 * is O(number of columns) and is what the rest of the engine reads the
 * differential through.
 *
 * With deep != 0 it also expands D_{i-1}(D_i(k)) term by term and
 * requires every coefficient to vanish.  That is the real check -- it
 * needs nothing but the answer itself, and everything the frame and
 * Schreyer's theorem and the reducer selection contribute shows up in it
 * -- but it is a polynomial multiplication per column and costs several
 * times the resolution it checks, so it is off by default.
 *
 * Returns the number of violations, so 0 means the differential is a
 * complex. */
int res_diff_verify(
        const res_diff_t * const rd,
        const int deep
        );

/* --------------------------------------------------------------------- *
 *  Minimal Betti numbers, by rank extraction
 *
 *  The frame gives the ranks of a *nonminimal* resolution.  The minimal
 *  ones come from it without ever building the minimal complex: for each
 *  level i and degree d let (d_i)_d be the part of the differential whose
 *  ring monomial is 1, a plain scalar matrix over F_p from the generators
 *  of F_i in degree d to those of F_{i-1} in degree d.  Then
 *
 *      beta_{i,d} = frame_{i,d} - rank (d_i)_d - rank (d_{i+1})_d,
 *
 *  because a constant entry of the differential cancels one generator
 *  against one generator one level down, and the ranks count how many such
 *  cancellations are independent.  Only ranks are ever needed -- no back
 *  substitution, no minimalized complex -- which is what makes this
 *  dramatically cheaper than resolving and then pruning.
 *
 *  The blocks are mutually independent across (i, d) and dense-ish, which
 *  is what makes them the natural device offload target.
 *
 *  The Hilbert numerator is cheaper again.  Writing the Hilbert series of
 *  the resolved module as K(t) / (1-t)^nv, the numerator K is
 *
 *      K_d = sum_i (-1)^i beta_{i,d} = sum_i (-1)^i frame_{i,d},
 *
 *  the two agreeing because the rank corrections telescope away.  So K
 *  needs no field arithmetic beyond the Gröbner basis: the frame alone
 *  determines it, and res_betti_new fills it in without ever looking at a
 *  differential.  This is Macaulay2's poincare, up to the degree shift
 *  applied to row_degs.  It is only correct for a *complete* frame, since
 *  a truncated one is missing the tail of the alternating sum.
 * --------------------------------------------------------------------- */

/*  Under a grading by Z^r (+) torsion every table here exists twice.
 *
 *  The *heft* tables -- frame, rank, betti, hilb -- are indexed by the
 *  single integer heft . a that schedules the computation, and are what a
 *  caller that does not care about the finer grading reads.  They remain
 *  correct under any grading: a heft class is a disjoint union of
 *  multidegree classes, the differential is multihomogeneous so its scalar
 *  part is block diagonal with respect to them, and ranks therefore add.
 *
 *  The *multigraded* tables -- mframe, mrank, mbetti, mhilb -- are indexed
 *  by bucket, mdegs listing the distinct multidegrees in the group's own
 *  order.  Under the standard grading the two indexings coincide and the
 *  multigraded tables are a permutation of the heft ones.
 *
 *  Rank extraction always runs on the *finer* blocks and sums up, which is
 *  not a compromise but a strict improvement: the blocks are smaller, so
 *  the elimination is cheaper, and the multigraded answer -- which
 *  Macaulay2 has no minimalBetti for -- comes out at no extra cost. */

typedef struct res_betti_t res_betti_t;
struct res_betti_t
{
    const res_frame_t *f;  /* not owned                                  */
    const res_dgrp_t *grp; /* not owned; f->grp                          */
    len_t    nlv;          /* levels, so 0 <= i < nlv                    */
    deg_t    maxdeg;       /* largest heft degree, so 0 <= d <= maxdeg   */
    int32_t *frame;        /* nlv*(maxdeg+1), the nonminimal ranks       */
    int32_t *rank;         /* nlv*(maxdeg+1), rank of the scalar part of
                            * d_i in degree d; row 0 is zero, there
                            * being no d_0                               */
    int32_t *betti;        /* nlv*(maxdeg+1) minimal Betti numbers, but
                            * a copy of frame until res_betti_minimalize
                            * has run                                    */
    int32_t *hilb;         /* maxdeg+1 Hilbert numerator coefficients    */

    hl_t     ndeg;         /* distinct multidegrees, so 0 <= u < ndeg    */
    int32_t *mdegs;        /* ndeg*grp->len, in the group's own order    */
    deg_t   *mheft;        /* ndeg, the heft degree of each bucket       */
    int32_t *mframe;       /* nlv*ndeg                                   */
    int32_t *mrank;        /* nlv*ndeg                                   */
    int32_t *mbetti;       /* nlv*ndeg                                   */
    int32_t *mhilb;        /* ndeg                                       */
    hl_t   **bkt;          /* bkt[i][k] is the bucket of frame element
                            * (i,k); this is the blocking rank extraction
                            * runs on, kept rather than rebuilt          */

    int      minimal;      /* 1 once res_betti_minimalize has run        */
    int      bad;
};

/* Tabulates the frame ranks and the Hilbert numerator.  The frame must be
 * complete for the numerator to mean anything. */
res_betti_t *res_betti_new(
        const res_frame_t * const f
        );

void res_betti_free(
        res_betti_t **bp
        );

/* Extracts the ranks of the scalar parts of rd and turns the frame ranks
 * into minimal Betti numbers.  rd must be the differential of the very
 * frame the table was built from.  Returns 0 on success. */
int res_betti_minimalize(
        res_betti_t *b,
        const res_diff_t * const rd
        );

/* Projective dimension and Castelnuovo-Mumford regularity, the largest i
 * and the largest d - i carrying a nonzero entry.  Both are -1 for the
 * zero module, and both are upper bounds rather than the real thing until
 * res_betti_minimalize has run. */
int32_t res_betti_pdim(
        const res_betti_t * const b
        );

int32_t res_betti_reg(
        const res_betti_t * const b
        );

/* Krull dimension and degree (multiplicity) from a Hilbert numerator of
 * len coefficients over nv variables: writing num(t) = (1-t)^c * G(t) with
 * G(1) != 0, the dimension is nv - c and the degree is G(1).  Both are
 * unaffected by the degree shift, which only multiplies num by a power of
 * t.  The zero module reports dimension -1 and degree 0.  Either output
 * may be NULL.  Returns 0 on success. */
int res_hilbert_invariants(
        const int32_t * const num,
        const len_t len,
        const len_t nv,
        int32_t *dim,
        int64_t *degree
        );

/* --------------------------------------------------------------------- *
 *  Gröbner bases of submodules of a free module
 *
 *  export_module_f4 is the module counterpart of export_f4 in f4.h: it
 *  takes a presentation matrix instead of a list of polynomials.  The
 *  columns of that matrix generate a submodule of R^nr_rows, and the flat
 *  input arrays extend msolve's usual convention with one component per
 *  term:
 *
 *    lens   nr_gens entries, the number of terms of each column
 *    exps   nr_vars entries per term, concatenated column by column
 *    comps  one entry per term, the 1-based row that term sits in
 *    cfs    one int32_t coefficient per term
 *
 *  and, describing the ambient free module,
 *
 *    row_degs  the degree of each generator of R^nr_rows, i.e.
 *              R^nr_rows = (+)_i R(-row_degs[i]); res_grading_len(grading)
 *              entries per row, concatenated row by row, so nr_rows
 *              entries under the standard grading.  May be NULL, which is
 *              read as all zero.  Only differences matter, so the degrees
 *              are normalized internally by subtracting the degree of the
 *              row of least heft; export_module_betti reports that shift.
 *
 *  stop asks for less than the whole basis; see res_stop_t, of which only
 *  max_degree means anything here, a plain Gröbner basis computing no
 *  syzygies to count or rows to project.  The ceiling is stated on the
 *  caller's scale, i.e. against the row_degs above and not against their
 *  normalization.  NULL, and a NULL max_degree, ask for everything.
 *
 *  The output uses the same layout, with bcomp giving the component of
 *  each term of the basis.  All four output arrays are allocated with the
 *  caller supplied mallocp and are released by
 *  free_module_f4_result_data.  The return value is the total number of
 *  terms, or 0 on failure, in which case nothing is allocated.
 *
 *  Under a grading other than the standard one the ring order is the
 *  *heft* degree reverse lexicographic order: monomials are compared by
 *  heft . deg first and ties broken reverse lexicographically.  That is a
 *  term order precisely because res_dgrp_of_grading insists the heft be
 *  strictly positive on every variable, which is also what makes the
 *  degree by degree schedules terminate.  With the standard grading the
 *  heft degree is the total degree and this is msolve's usual DRL.
 *
 *  export_module_f4_blocks takes a block description on top of that (see
 *  order.h), so the ring part of the module order can be a block grevlex
 *  order with any number of blocks and per variable weights; NULL asks
 *  for the single block order, which is what export_module_f4 passes.
 *  The component key is unaffected -- the blocks refine only the ring
 *  part -- so POT and TOP mean exactly what they meant before.
 *
 *  Weights are one notion, not two: if the block description carries
 *  weights and the grading is not the standard one, the two have to
 *  agree, and disagreeing is an error rather than a silent choice.
 *
 *  Restrictions, all reported on stderr rather than assumed:
 *    - prime field of characteristic 0 < p < 2^31,
 *    - degree reverse lexicographic order (mon_order 0); blocks refine
 *      DRL, so a block description does not combine with mon_order 1,
 *    - strat is checked by res_strat_check: its base is RES_MORD_POT or
 *      RES_MORD_TOP, and NULL asks for the default.  RES_MORD_SCHREYER
 *      needs per component base monomials that only the resolution engine
 *      can supply.
 *
 *  The frame, resolution and Betti entry points below keep the single
 *  block order: their Schreyer lift is only implemented for it.
 * --------------------------------------------------------------------- */

int64_t export_module_f4_blocks(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,      /* number of basis elements */
        int32_t **blen,    /* number of terms of each basis element */
        int32_t **bexp,    /* nr_vars exponents per term */
        int32_t **bcomp,   /* 1-based component of each term */
        void **bcf,        /* one int32_t coefficient per term */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const mo_block_t * const blk, /* NULL means one block */
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
        const res_stop_t *stop,       /* NULL means no limits    */
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t info_level
        );

int64_t export_module_f4(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,      /* number of basis elements */
        int32_t **blen,    /* number of terms of each basis element */
        int32_t **bexp,    /* nr_vars exponents per term */
        int32_t **bcomp,   /* 1-based component of each term */
        void **bcf,        /* one int32_t coefficient per term */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
        const res_stop_t *stop,       /* NULL means no limits    */
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t info_level
        );

void free_module_f4_result_data(
        void (*freep) (void *),
        int32_t **blen,
        int32_t **bexp,
        int32_t **bcomp,
        void **bcf
        );

/* --------------------------------------------------------------------- *
 *  Schreyer frame of a presentation matrix
 *
 *  Takes exactly the input of export_module_f4 and returns the ranks of
 *  the frame instead of the Gröbner basis: betti[i*(*maxdeg+1) + d] is the
 *  number of level i frame elements of degree d, for 0 <= i < *nlevels and
 *  0 <= d <= *maxdeg.  Level 0 is the ambient free module, level 1 the
 *  Gröbner basis of the submodule, and level i+1 the Schreyer syzygies of
 *  level i, so these are the ranks of a nonminimal free resolution.
 *
 *  Degrees are shifted so that the smallest of row_degs is zero, exactly
 *  as in export_module_f4; the caller can shift them back.
 *
 *  max_level truncates the computation after that many levels; pass 0 to
 *  run to completion.  Note the frame is a *nonminimal* resolution and can
 *  be longer than nr_vars -- (z, y^2, x^2 y, x^3) in three variables has
 *  the frame 1,4,6,4,1 -- so nr_vars is not a safe ceiling.  betti is
 *  allocated with mallocp and
 *  released by free_module_frame_result_data.  The return value is the
 *  total number of frame elements, or 0 on failure, in which case nothing
 *  is allocated.
 * --------------------------------------------------------------------- */

int64_t export_module_frame(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *nlevels,
        int32_t *maxdeg,
        int32_t **betti,
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t max_level,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level
        );

void free_module_frame_result_data(
        void (*freep) (void *),
        int32_t **betti
        );

/* --------------------------------------------------------------------- *
 *  Nonminimal free resolutions and syzygies
 *
 *  Two things are called a syzygy and they are not the same object:
 *
 *    RES_SYZ_OF_GB     syzygies of the *Gröbner basis* of the submodule.
 *                      This is what the Schreyer frame produces and what
 *                      the resolution is built out of; d_1 is the Gröbner
 *                      basis and d_2 its first syzygies.
 *    RES_SYZ_OF_INPUT  syzygies of the *generators the caller gave*, which
 *                      is what Macaulay2's syz returns.  msolve discards
 *                      the change of basis from the input to the Gröbner
 *                      basis, so these are computed instead by the graph
 *                      module trick: take the Gröbner basis of the columns
 *                      (f_j, e_j) of R^nr_rows (+) R^nr_gens under a
 *                      position over term order with the R^nr_rows block
 *                      first, and keep the elements whose R^nr_rows part
 *                      vanishes.  Here d_1 is the input matrix itself.
 *
 *  In both cases the result is a complex of free modules
 *
 *      F_0 <-- F_1 <-- ... <-- F_{nlevels-1},
 *
 *  reported as
 *
 *    ranks   nlevels entries, the rank of each F_i
 *    degs    sum(ranks) entries, the degree of every generator of every
 *            F_i, concatenated level by level; shifted so that the
 *            smallest row degree is zero, exactly as in export_module_f4
 *    dlen    one entry per generator of F_1, ..., F_{nlevels-1}, again
 *            concatenated level by level: the number of terms of that
 *            column of the differential
 *    dexp    nr_vars exponents per term
 *    dcomp   one entry per term, the 1-based generator of F_{i-1} it sits
 *            in
 *    dcf     one int32_t coefficient per term
 *
 *  max_level truncates at that level, so max_level = 2 is the single
 *  syzygy matrix and 0 runs to completion.  There is no level that is
 *  always enough: the complex reported here is the nonminimal one and can
 *  run past nr_vars, unlike the minimal resolution Hilbert's syzygy
 *  theorem bounds.  RES_SYZ_OF_INPUT accepts no max_level above 2, since
 *  resolving beyond the first syzygies is the Gröbner basis story again.
 *
 *  verify != 0 additionally runs the exact d o d = 0 check of
 *  res_diff_verify over the whole complex, refusing to report a result
 *  that is not a complex.  It costs several times the resolution itself,
 *  so it is a request, not the default.  It has no effect on
 *  RES_SYZ_OF_INPUT, where d_1 o d_2 = 0 is structural: a column is
 *  reported only once every one of its terms has been seen to sit in the
 *  adjoined components, and that is exactly the statement that its
 *  original components cancel.  A row bound breaks that -- see stop below.
 *
 *  The differential is the *nonminimal* one: its ranks are the frame
 *  ranks of export_module_frame, not the minimal Betti numbers.
 *
 *  stop asks for less than the whole thing; see res_stop_t.  All three
 *  fields reach RES_SYZ_OF_INPUT, and what they mean here is:
 *
 *    max_degree  a ceiling on the degree of the Gröbner basis the
 *                syzygies are read off, so on the graph module for
 *                RES_SYZ_OF_INPUT and on the presentation itself for
 *                RES_SYZ_OF_GB.  What comes back is then a truncation of
 *                the complex rather than the complex.
 *
 *    syz_limit   stops as soon as that many syzygies are known and
 *                reports that many.  It is a genuine early stop, not a
 *                cap on the output: the round loop counts the relations
 *                the basis has picked up and quits.  A round is always
 *                finished, so more than syz_limit of them can turn up at
 *                once, and the extras are dropped rather than reported.
 *                Only RES_SYZ_OF_INPUT counts anything.
 *
 *    syz_rows    projects each column of d_2 onto the first syz_rows
 *                generators of F_1 and drops the columns that become
 *                zero, so ranks[2] can come back smaller.  The result is
 *                a submatrix of the syzygy matrix and *not* a complex;
 *                verify is ignored, as it would have to be.  Only
 *                RES_SYZ_OF_INPUT projects anything.
 *
 *  All six arrays are allocated with mallocp and released by
 *  free_module_resolution_result_data.  The return value is the total
 *  number of terms, or 0 on failure, in which case nothing is allocated.
 *  A resolution with no terms at all cannot occur, since every column of
 *  d_1 is nonzero.
 * --------------------------------------------------------------------- */

typedef enum {
    RES_SYZ_OF_GB    = 0,
    RES_SYZ_OF_INPUT = 1
} res_syz_t;

int64_t export_module_resolution(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *nlevels,
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf,
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
        const res_stop_t *stop,       /* NULL means no limits    */
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
        );

void free_module_resolution_result_data(
        void (*freep) (void *),
        int32_t **ranks,
        int32_t **degs,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        );

/* --------------------------------------------------------------------- *
 *  Minimal Betti numbers and Hilbert information
 *
 *  Takes exactly the input of export_module_f4 and reports invariants of
 *  the module it presents, that is of R^nr_rows modulo the submodule the
 *  columns generate.  Nothing here materializes a resolution: with
 *  minimal != 0 the nonminimal differential is computed and immediately
 *  reduced to the ranks of its scalar parts, and with minimal == 0 not one
 *  field operation happens past the Gröbner basis.
 *
 *    betti     (*nlevels)*(*maxdeg+1) entries, the number of generators of
 *              the minimal resolution at level i in degree d at index
 *              i*(*maxdeg+1) + d.  With minimal == 0 these are the frame
 *              ranks instead, which bound the minimal ones from above.
 *              Level 0 is R^nr_rows, so an ideal starts with a single 1.
 *    hilbnum   *maxdeg+1 coefficients of the numerator of the Hilbert
 *              series over (1-t)^nr_vars, low degree first; Macaulay2
 *              calls this poincare.  It is the alternating sum of either
 *              table and so costs nothing beyond the frame, but it needs
 *              a complete resolution: max_level must not truncate.
 *    pdim      projective dimension, the largest level with a generator
 *    reg       Castelnuovo-Mumford regularity, the largest d - i
 *    dimension Krull dimension of the module, nr_vars minus the order of
 *              vanishing of the numerator at t = 1
 *    degree    degree (multiplicity), the numerator divided by that many
 *              factors of (1-t) and evaluated at 1
 *
 *  betti and hilbnum are arrays indexed by degree, so they have to start
 *  at zero: their degrees are shifted so that the smallest of row_degs is
 *  zero, exactly as everywhere else here, and degshift reports the shift
 *  that was applied.  The caller's own degrees are betti's plus degshift,
 *  and its own numerator is hilbnum times t^degshift.  The four scalars
 *  are under no such constraint and are reported in the caller's own
 *  degrees: reg has degshift already added, and pdim, dimension and degree
 *  do not depend on it at all.  The zero module reports dimension -1,
 *  degree 0 and pdim = reg = -1.
 *
 *  All of that is in terms of the single *heft* degree, and it stays
 *  correct under any grading -- a heft class is a disjoint union of
 *  multidegree classes.  The finer table lives in mtab, below, which is
 *  computed only when one is supplied.  dimension and degree in particular
 *  need no generalization: they are read off the heft numerator, and
 *  because dividing by (1 - t^w) contributes exactly one zero at t = 1 per
 *  variable whatever the weight w, the same "write K = (1-t)^c G, then
 *  dim = nv - c and degree = G(1)" recovers what Macaulay2 reports for
 *  weighted and multigraded rings alike.
 *
 *  max_level truncates the reported table at that level.  The frame is
 *  still built one level further, since beta at level i reads the rank of
 *  d_{i+1}, so the top row of a truncated table is exact rather than an
 *  upper bound.  pdim, reg, dimension, degree and hilbnum are refused on a
 *  truncated table: they are invariants of the whole module and there is
 *  no level past which the frame is known to be empty -- it is a
 *  nonminimal resolution and really can run past nr_vars.
 *
 *  Every return value except nlevels and maxdeg may be NULL and is then
 *  not computed.  betti and hilbnum are allocated with mallocp and are
 *  released by free_module_betti_result_data.  The return value is the
 *  number of frame elements built, which for a truncated table counts the
 *  extra level too, or 0 on failure, in which case nothing is allocated.
 * --------------------------------------------------------------------- */

/* The multigraded half of export_module_betti: the same Betti numbers and
 * the same Hilbert numerator, indexed by the multidegrees that actually
 * occur rather than by a single heft degree.  Macaulay2 has no
 * minimalBetti for a multigraded module, so this is the novel output.
 *
 *   dlen      int32 slots per multidegree, res_grading_len(grading)
 *   ndegs     number of distinct multidegrees occurring anywhere in the
 *             resolution, so a bucket u runs over 0 <= u < ndegs
 *   degs      ndegs * dlen, the multidegrees themselves, in the grading
 *             group's own order, free part first then torsion residues
 *   heft      ndegs, the heft degree of each bucket, i.e. the column of
 *             the heft indexed table each bucket contributes to
 *   betti     nlevels * ndegs at index i*ndegs + u
 *   hilbnum   ndegs, the multigraded Hilbert numerator: the coefficient of
 *             the monomial degs[u] is hilbnum[u].  NULL on a truncated
 *             table, the numerator being an alternating sum over every
 *             level, exactly as max_level refuses the heft indexed one
 *   degshift  dlen, the multidegree subtracted from every row degree, so
 *             the caller's own degrees are degs plus this
 *
 * degs, heft, betti and hilbnum are allocated with mallocp and released by
 * free_module_mtable_data; degshift is a fixed size array inside the
 * struct, so nothing has to be freed when only it is wanted.  The struct
 * is a struct rather than another eight parameters so that a later grading
 * feature -- a torsion label, a group element handle -- is a new field
 * instead of a new signature. */

#define RES_MTAB_MAXLEN 64  /* slots the inline degshift holds */

typedef struct res_mtable_t res_mtable_t;
struct res_mtable_t
{
    int32_t  nlevels;
    int32_t  ndegs;
    int32_t  dlen;
    int32_t  degshift[RES_MTAB_MAXLEN];
    int32_t *degs;
    int32_t *heft;
    int32_t *betti;
    int32_t *hilbnum;
};

void free_module_mtable_data(
        void (*freep) (void *),
        res_mtable_t *mtab
        );

int64_t export_module_betti(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *nlevels,
        int32_t *maxdeg,
        int32_t *degshift,   /* may be NULL */
        int32_t **betti,     /* may be NULL */
        int32_t **hilbnum,   /* may be NULL */
        int32_t *pdim,       /* may be NULL */
        int32_t *reg,        /* may be NULL */
        int32_t *dimension,  /* may be NULL */
        int64_t *degree,     /* may be NULL */
        res_mtable_t *mtab,  /* may be NULL */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
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
        );

void free_module_betti_result_data(
        void (*freep) (void *),
        int32_t **betti,
        int32_t **hilbnum
        );

/* --------------------------------------------------------------------- *
 *  A resolution kept alive: ranks now, differentials on demand
 *
 *  export_module_resolution is a single shot -- it resolves everything and
 *  hands back the whole complex in flat arrays.  That is the wrong shape
 *  for a caller that wants to look at a large resolution before deciding
 *  what to materialize, which is what Macaulay2 does: rawResolutionGetFree
 *  asks only for the rank and the degrees of F_i, and rawResolutionGetMatrix
 *  asks for one differential at a time, each of them possibly much later.
 *
 *  res_comp_t is that computation, kept across calls.  res_comp_new runs
 *  the module Gröbner basis and builds the whole Schreyer frame, which is
 *  combinatorial -- no field arithmetic happens past the Gröbner basis --
 *  and so answers every question about the shape of the resolution
 *  immediately.  The differential is filled in only when asked for, one
 *  request at a time, and remembered.
 *
 *  The complex reported is the *nonminimal* one, exactly as with
 *  export_module_resolution: the ranks are the frame ranks, not the
 *  minimal Betti numbers, and there is no level past which the frame is
 *  known to be empty.  For minimal Betti numbers use export_module_betti,
 *  which never materializes a differential at all.
 *
 *  Level 0 is the ambient free module R^nr_rows, level 1 the Gröbner basis
 *  of the submodule -- not the input generators, which msolve does not
 *  keep a change of basis to -- and level i+1 the Schreyer syzygies of
 *  level i.
 *
 *  The differential is only ever computed as a prefix: asking for D_i
 *  computes D_2 ... D_i, because the block at level i in degree d reduces
 *  against level i-1 in degrees up to and including d.  Asking for them
 *  out of order is allowed and costs the same as asking in order; asking
 *  for one already computed costs nothing.
 *
 *  Ownership: res_comp_new returns a handle the caller owns and must
 *  release with res_comp_free.  Nothing else in this header allocates
 *  anything that outlives the call, so this is the one place a caller has
 *  to keep track of a lifetime.  The Gröbner basis and its hash table are
 *  released inside res_comp_new, once level 1 of the differential has been
 *  read off them; what the handle holds afterwards is the frame, its own
 *  ring hash table, and whatever differential has been asked for.
 * --------------------------------------------------------------------- */

typedef struct res_comp_t res_comp_t;

/* Takes exactly the input of export_module_f4, plus max_level, which
 * truncates the frame at that level; pass 0 for no ceiling.  Returns NULL
 * on failure, reporting the reason on stderr, exactly as the export_module
 * entry points do. */
res_comp_t *res_comp_new(
        const int32_t *lens,
        const int32_t *exps,
        const int32_t *comps,
        const void *cfs,
        const int32_t *row_degs, /* may be NULL */
        const uint32_t field_char,
        const int32_t mon_order,
        const res_strat_t *strat,     /* NULL means the default  */
        const res_grading_t *grading, /* NULL means the standard */
        const int32_t nr_vars,
        const int32_t nr_rows,
        const int32_t nr_gens,
        const int32_t max_level,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t la_option,
        const int32_t info_level
        );

void res_comp_free(
        res_comp_t **cp
        );

/* The number of levels the frame has, so the free modules are F_0 to
 * F_{nlevels-1} and the differentials D_1 to D_{nlevels-1}.  0 if the
 * handle is unusable. */
int32_t res_comp_nlevels(
        const res_comp_t * const c
        );

/* The shift that was applied to row_degs, i.e. the smallest of them: every
 * degree reported here is the caller's minus this. */
int32_t res_comp_degshift(
        const res_comp_t * const c
        );

/* Whether the frame ended on its own rather than being cut off at
 * max_level.  A truncated resolution is still exact everywhere it is
 * reported; it just does not end. */
int res_comp_is_complete(
        const res_comp_t * const c
        );

/* The rank of F_level, or -1 if level is out of range. */
int32_t res_comp_rank(
        const res_comp_t * const c,
        const int32_t level
        );

/* The degrees of the generators of F_level, res_comp_rank(c, level)
 * entries, written into degs, which the caller supplies.  Returns 0 on
 * success.  Needs no field arithmetic: this is the frame. */
int res_comp_degrees(
        const res_comp_t * const c,
        const int32_t level,
        int32_t *degs
        );

/* The number of int32 slots one multidegree occupies here, so 1 under the
 * standard grading.  0 if the handle is unusable. */
int32_t res_comp_glen(
        const res_comp_t * const c
        );

/* The multidegrees of the generators of F_level, res_comp_glen(c) entries
 * each and res_comp_rank(c, level) of them, written into mdegs.  These are
 * the shifted degrees, exactly as res_comp_degrees reports the shifted
 * heft degree; res_comp_multidegshift is the shift.  Returns 0 on success. */
int res_comp_multidegrees(
        const res_comp_t * const c,
        const int32_t level,
        int32_t *mdegs
        );

/* The multidegree that was subtracted from every row degree, res_comp_glen
 * entries written into shift.  Its heft is res_comp_degshift. */
int res_comp_multidegshift(
        const res_comp_t * const c,
        int32_t *shift
        );

/* D_level, as one column per generator of F_level in the flat layout
 * export_module_resolution uses:
 *
 *   dlen   one entry per generator of F_level, the number of terms
 *   dexp   nr_vars exponents per term, concatenated column by column
 *   dcomp  one entry per term, the 1-based generator of F_{level-1} it
 *          sits in
 *   dcf    one int32_t coefficient per term
 *
 * Computes whatever is missing below level as well, and remembers all of
 * it.  The four arrays are allocated with mallocp and are released by
 * free_module_differential_data; the handle keeps its own copy and stays
 * usable.  Returns the number of terms, or 0 on failure, in which case
 * nothing is allocated.  D_level of a zero F_level is 0 terms and is
 * reported as a failure only if level itself is out of range, so check
 * res_comp_rank first. */
int64_t res_comp_differential(
        void *(*mallocp) (size_t),
        res_comp_t *c,
        const int32_t level,
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        );

void free_module_differential_data(
        void (*freep) (void *),
        int32_t **dlen,
        int32_t **dexp,
        int32_t **dcomp,
        void **dcf
        );

#ifdef __cplusplus
}
#endif

#endif /* GB_RES_H */
