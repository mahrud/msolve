# Syzygies and minimal free resolutions on an msolve backend

## Context

**Goal.** Implement state-of-the-art computation of syzygies and minimal free resolutions
for graded modules over graded k-algebras, k = Z/pZ, using msolve's F4 engine
(`libneogb`) as the Gröbner-basis and linear-algebra substrate. Multigradings are a
high-priority extension; simple non-commutative rings (exterior algebras, quiver path
algebras) and other FLINT-backed finite fields are lower priority. The linear algebra
should be swappable across CPU/GPU/NPU.

### State of the art: confirmed, with caveats

The claimed pipeline — **Schreyer frame → nonminimal differential → minimal Betti by rank
extraction, in an F4-style linear algebra framework** — is correct for the primary case
(standard-graded k[x₁…xₙ] over 𝔽_p). It is what Macaulay2's `res(…, Strategy => Nonminimal)`
and `minimalBetti` implement (Schreyer–Stillman, `Macaulay2/e/schreyer-resolution/`),
descending from La Scala–Stillman, *Strategies for Computing Minimal Free Resolutions*
(J. Symbolic Comput. 1998).

Mechanism: Schreyer's theorem determines the lead terms of the syzygies at level i+1 from
the lead terms at level i, so the **entire resolution skeleton (the frame) is combinatorial** —
computed with zero field arithmetic. The differential is then filled in level-by-level,
degree-by-degree by Macaulay-matrix reduction over k. Minimal Betti numbers come from ranks
of the *scalar* (degree-0) part of the differentials:

    β^min_{i,j} = frame_{i,j} − rank(D_i)_{i,j} − rank(D_{i+1})_{i,j}

Only **ranks** are ever needed, never the minimalized complex. That is why `minimalBetti`
is orders of magnitude faster than `res` followed by pruning.

Three caveats that shape this plan:

1. **Not universally dominant.** Erocal–Motsak–Schreyer–Steenpaß, *Refined Algorithms to
   Compute Syzygies* (arXiv:1502.01654; Singular's `fres`) uses LiftHybrid/LiftTree —
   polynomial arithmetic with term-skipping and cached partial lifts, not batched linear
   algebra — and wins on sparse/monomial-flavored input. Frame + LA is SOTA **when linear
   algebra is the bottleneck**, which is precisely the regime that justifies GPU offload.
   This reinforces the plan rather than undercutting it, but it means benchmarks must
   include sparse inputs where we expect to lose.

2. **"F4-style" is a favorable misnomer.** Step 3 has **no S-pairs and no pair selection** —
   the frame *is* the schedule, and every matrix shape is known before any arithmetic
   happens. Real F4's bottleneck is irregular, memory-latency-bound symbolic preprocessing
   (cf. arXiv:2601.06765). The resolution engine has none of that, which makes it a far
   better offload target than F4 itself. **Do not reuse `select_spairs_by_minimal_degree`
   or `update_basis_f4`.**

3. **Every existing implementation is prime-field and singly-graded.** M2's `minimalBetti`
   explicitly refuses multigraded input. The multigraded case is unimplemented *anywhere*,
   not merely unimplemented in msolve — it is the highest-value novel contribution here.

### What msolve actually provides

Less than one might hope; worth stating plainly so the plan is honestly scoped.

- **No module/component support anywhere.** `ht->ev` is `[deg, e₁…eₙ]` (`src/neogb/data.h:101-117`).
- The **SBA path is a dead end as-is**: it stores signature monomials but discards syzygies,
  is 32-bit/homogeneous/serial, and carries two allocation bugs
  (`basis.c:209-251` never reallocs `sm`/`si`; `basis.c:404-416` memcpys into unallocated
  `sm`/`si`). `engine.c:137,154,170` `exit(1)` whenever `use_signatures != 0`.
- `core_nf` reduces polynomials mod a GB but **discards the reduction coefficients**.

The genuinely reusable asset is `exact_sparse_reduced_echelon_form_sat_ff_32`
(`la_ff_32.c:2772-2996`): it maintains a dense row `drl` over monomial columns **and** a
parallel row `drm` over multiplier columns, applies every reduction to both, and emits an
explicit coefficient vector when `drl` reaches zero. That is a working syzygy kernel; it
needs generalizing from one saturation element to a rank-r free module. Its row kernel is
`reduce_dense_row_by_known_pivots_sparse_sat_ff_31_bit` (`la_ff_32.c:1245`), with
`adjust_multiplier_sparse_matrix_row_ff_32` (`la_ff_32.c:257`) keeping the two sides
consistent under renormalization.

So msolve contributes: hash table + monomial arithmetic, the sparse 𝔽_p row-reduction
kernels (with AVX2/AVX-512/NEON), the level-0/1 Gröbner basis, and the I/O plumbing.
Everything module-theoretic is new code.

### Decisions taken

| Axis             | Decision                                                                            |
|------------------|-------------------------------------------------------------------------------------|
| Placement        | **In-tree fork**, new `src/neogb/res_*.c` added to `gb.c`'s include list            |
| Milestone 1      | **Both, staged** — one core (frame + differential), two consumers (ranks, complex)  |
| Module monomials | **Component slot in `ht_t`/`ev`**, plus a lifted `total` monomial per frame element |
| Multigrading     | **Designed in from day one**; degrees are an opaque **struct**                      |
| Input            | **Modules from day one**, not just ideals — an ideal is the rank-1 case             |
| Outputs          | **Single syzygy matrix** is a first-class output alongside full resolutions         |

---

## Architecture

### Placement and build

`libneogb` is a **unity build**: `src/neogb/gb.c:29-51` `#include`s every `.c` file and
`src/neogb/Makefile.am:2` compiles only `gb.c`. Consequently nearly everything is `static`.
New files are added to that include list, after `la_ff_32.c` and before `engine.c`:

```c
#include "res_grading.c"   /* ZZ^r degrees, heft, multidegree buckets */
#include "res_order.c"     /* Schreyer/POT/TOP module orders          */
#include "res_frame.c"     /* Schreyer frame construction (monomial only) */
#include "res_diff.c"      /* nonminimal differential, degree-by-degree driver */
#include "res_betti.c"     /* rank extraction, Betti table, Hilbert numerator */
#include "res_module.c"    /* module GB, frame, resolution and Betti entry points */
#include "res_la.c"        /* backend vtable + CPU reference kernels */
```

All but `res_la.c`, which is M8, exist.  The flat-array C API that was sketched as
`res_export.c` lives in `res_module.c` next to the entry points that need
`module_gb_from_input`, rather than in a file of its own.

Public headers go in `src/neogb/res.h`, added to `libneogb.h` and to `Makefile.am:4-6`.
Note `libneogb.h` has **no `extern "C"` guard** (unlike `msolve.h`) — add one to `res.h`.

Keep every change to existing files behind a predicate that is false in the non-module
path, so the plain F4 path stays byte-identical and upstream rebases stay cheap.

### Module monomials in the hash table

Place the component slot at the **end** of the exponent vector, not the front:

```
non-module: ev = [deg, e₁ … eₙ]              evl = nv + 1        (unchanged)
module:     ev = [deg, e₁ … eₙ, comp]        evl = nv + 2
block+mod:  ev = [deg_b1, e₁…e_k, deg_b2, e_{k+1}…eₙ, comp]      evl = nv + 3
```

Why the end, rather than shifting `DEG`: `DEG` is `OFFSET-6 == 0` (`data.h:58`) and is
**deliberately overloaded** — `ev[DEG]` is the degree slot of an exponent vector while
`hm[DEG]` is the degree field of a polynomial's `hm_t` array. Changing it breaks both at
every one of ~25 call sites. A trailing slot leaves `ev[DEG]` at index 0 and all existing
degree logic untouched. It also leaves the divmask machinery untouched, because `ht->dv[]`
holds **explicit variable-slot indices** (`hash.c:75-91`, consumed at `hash.c:410, 441-453`),
so a slot not listed in `dv` is simply invisible to divmasks.

New `ht_t` fields (`data.h:120-144`):

```c
len_t  cpos;    /* index of the component slot in ev; 0 = not a module table */
len_t  ncomp;   /* number of components */
hm_t  *cbase;   /* per-component Schreyer base monomial (hash in the ring table) */
deg_t *cshift;  /* per-component heft-degree shift */
int32_t mord;   /* module order: RES_SCHREYER | RES_POT | RES_TOP */
```

Three properties make this work:

- **Hashing is linear** — `h = Σ rn[j]·a[j]` (`hash.c:836-838`) — so multiplying a module
  monomial by a ring monomial (which has `comp == 0`) is still hash *addition*. No change
  to `insert_multiplied_signature_in_hash_table` (`hash.c:639`).
- **`ev[DEG]` holds the heft degree including the component's shift**, baked in at
  construction. Vector addition then keeps degrees correct automatically, since ring
  monomials contribute shift 0. Every existing degree computation works unchanged.
- Two module monomials are never added, so the `comp` slots never collide.

Two things genuinely must change, both guarded on `ht->cpos != 0`:

1. **Divisibility.** `check_monomial_division` (`hash.c:478-507`) and
   `check_monomial_division_in_update` (`hash.c:509`) loop over *all* slots and would treat
   `comp` as an exponent, wrongly accepting a divisor with a smaller component. The module
   variant must require `eb[cpos] == 0` (ring divisor) or `ea[cpos] == eb[cpos]`.
2. **Order.** The reverse-lex loops run `i = ht->evl-1; while (i > 1 && ea[i] == eb[i]) --i;`
   (`order.c:388-392, 421-425, 452-456`, and the `monomial_cmp_*` family from `order.c:459`).
   In module mode they must start at `evl-2` and handle the component per `ht->mord`.

Extend the five dispatchers in `io.c:789-866` with a `use_module_order(ht)` predicate
alongside the existing `use_block_order` / `use_lex_order` (`io.c:768-780`). **Follow this
convention rather than adding globals** — commit `4122be5` deliberately moved order
selection off globals and onto `ht_t` to make concurrent computations re-entrant, and
`set_function_pointers` / `reset_function_pointers` are now no-op stubs (`io.c:1301-1306`,
`io.c:1333-1340`). `md_t.use_signatures` already reserves `1=SCHREYER, 2=POT, 3=DPOT`
(`data.h:409-413`) but nothing branches on it; promote that to a real enum in `res.h`.

### Multigrading, and the road to equivariant resolutions

Degrees are graded by an arbitrary finitely generated abelian group A = Z^r ⊕ T with T
finite (Cl(X) for a toric X; character groups for equivariant resolutions). A degree is an
**opaque struct**, never a bare array, so that a dedicated group-arithmetic library can back
it later without touching a single call site:

```c
/* One element of the grading group A. Layout is defined by the owning
 * res_dgrp_t and is opaque to every consumer. Kept a struct (not a bare
 * int32_t*) so it can grow a cached hash, a torsion handle, or a pointer
 * into an external group library without an API break. */
typedef struct res_deg_t {
  int32_t *e;   /* [0,r) free part; [r,r+nt) torsion residues */
} res_deg_t;

/* The grading group itself: owns arithmetic, storage layout, and the heft. */
typedef struct res_dgrp_t {
  len_t    r;      /* rank of the free part                          */
  len_t    nt;     /* number of torsion factors                      */
  int32_t *tord;   /* orders of the torsion factors, length nt       */
  len_t    len;    /* int32 slots per degree = r + nt                */
  int32_t *heft;   /* length r; heft . deg(x_j) > 0 for every j      */
  int32_t *dmat;   /* len x n degree matrix; column j = deg(x_j)     */
  /* vtable — the seam for an external group-arithmetic backend */
  void (*add)(const struct res_dgrp_t *, res_deg_t *d, const res_deg_t *a, const res_deg_t *b);
  void (*sub)(const struct res_dgrp_t *, res_deg_t *d, const res_deg_t *a, const res_deg_t *b);
  int  (*cmp)(const struct res_dgrp_t *, const res_deg_t *a, const res_deg_t *b);
  hl_t (*hash)(const struct res_dgrp_t *, const res_deg_t *a);
  deg_t(*heft_of)(const struct res_dgrp_t *, const res_deg_t *a);
} res_dgrp_t;
```

Torsion arithmetic is reduction mod `tord[i]` in the default backend — enough for Cl(X) of
a toric variety — and the vtable is where a real group library plugs in for the equivariant
case. Provide specialized `add`/`cmp`/`hash` for the common `r==1, nt==0` case so the
singly-graded path pays nothing for the generality.

The key economy: **the multidegree of a monomial is never stored** — it is determined by the
exponent vector as `dmat·a + cdeg[comp]` and computed on demand. Only `ev[DEG]`, holding a
single int32 **heft degree**, is stored, and that is what schedules and orders the
computation. Since msolve's DRL already uses `ev[DEG]` as the degree, DRL becomes heft-DRL
for free. Only **frame elements** store a full `res_deg_t` — there are vastly fewer of them
than monomials. Rank extraction buckets by degree through `cmp`/`hash`.

### Module input and the single-syzygy-matrix output

Input is a set of elements of a free module F = ⊕ᵢ R(−dᵢ); an ideal is the rank-1 case with
no degree shift. **This forces one thing msolve cannot currently do: a module Gröbner
basis.** msolve's F4 computes ideal GBs only.

The good news is that a module F4 falls out of the component slot almost for free. With
component-aware divisibility and order in place, the only other component-aware pieces are:

- `get_lcm` (`hash.c:1306`) — must return "no pair" when components differ, since two module
  monomials in different components have no lcm.
- `insert_and_update_spairs` (`update.c:64-190`) — skip cross-component pairs entirely; the
  Gebauer–Möller criteria then apply unchanged within each component.
- `find_multiplied_reducer` (`symbol.c:493-556`) — already goes through
  `check_monomial_division`, so it becomes correct automatically.

Everything else — symbolic preprocessing, column assignment, the linear algebra, basis
update — is component-agnostic and needs no change. This is why the component slot is worth
the invasiveness: it buys the module GB, the Schreyer frame, and multigrading with one
mechanism.

**Two distinct notions of "syzygy" must be kept straight**, and both are exposed:

1. **Syzygies of the Gröbner basis** — what the Schreyer frame produces, and what the
   resolution consumes. Cheap, and a direct byproduct of the machinery.
2. **Syzygies of the user's generators** — what M2's `syz` returns. These require the
   change of basis from input generators to GB elements, which msolve discards. Compute them
   via the standard graph-module trick: take the GB of `{(gᵢ, eᵢ)} ⊂ F ⊕ Rᵐ` under a POT
   order with the F-block first, and read off the elements with zero F-component. This needs
   only the module order machinery already being built.

Expose both through one entry point, since the difference is invisible in the rank-1
homogeneous ideal case that most users will hit first:

```c
res_t *res_compute(const res_input_t *in, const res_opts_t *opts, int32_t *err);
/* opts->maxlevel : 1 = single syzygy matrix; n = resolution truncated at level n */
/* opts->minimize : RES_NONMINIMAL | RES_MINIMAL_BETTI | RES_MINIMAL           */
/* opts->syz_of   : RES_SYZ_OF_GB | RES_SYZ_OF_INPUT                          */
```

### Schreyer frame

```c
typedef struct {
  hm_t     mono;   /* own monomial, in the module hash table at this level */
  hm_t     total;  /* fully lifted monomial, in the ring hash table       */
  int32_t  up;     /* index into level i-1                                */
  deg_t    hdeg;   /* heft degree                                         */
  int32_t *mdeg;   /* ZZ^r degree                                         */
} res_felt_t;
```

`total` is the load-bearing trick: it is the image of the element's lead term lifted all the
way down to level 0, so **Schreyer comparison at any level reduces to a ring
`monomial_cmp(total_a, total_b)`, tie-broken by descending `up`**. Without it, every
comparison walks the level chain and thrashes cache. (This is effectively what Stillman's
M2 implementation does.)

Construction is pure combinatorics — no coefficients, no field arithmetic:
level 0 = the free module generators; level 1 = the GB lead terms from msolve; level i+1
from level i by the induced Schreyer computation. Reuse `check_monomial_division`
(`hash.c:478`) and the degree-by-degree enumeration pattern of `quotient_basis`
(`f4sat.c:302-372`), which already walks standard monomials by degree using
`check_lm_divisibility_and_insert_in_hash_table` (`hash.c:549`).

### Nonminimal differential

~~Drive by **slanted degree** `s = d − i` ascending (La Scala–Stillman), and within a slanted
degree by level.~~ **Slanted degree is not a valid schedule.** A reducer used with
multiplier 1 has the *same* degree as the row it reduces, so level *i* in degree *d* can
need level *i−1* in degree *d*, whose slanted degree is `s+1` — one strand that has not been
computed yet, in either direction of the inner loop. Drive by **degree ascending, then level
ascending** instead: `(i−1, d)` is exactly the step before `(i, d)`. See the M4 status note;
a mutation to the slanted-degree order fails three selftest checks.

For each (level i, degree d): rows = frame elements at (i, d), columns = module monomials at
level i−1 in degree d. Reduce against the level i−1 pivot data using a **rank-r
generalization of the sat dual-row trick**: dense row `drl` over monomial columns, parallel
row `drm` over expression columns, both updated by every reduction
(`la_ff_32.c:2896-2971`; row kernel `la_ff_32.c:1245`; renormalization
`la_ff_32.c:257`). The existing code manages two simultaneous column spaces by
pointer-swapping a single `md->hcm` field (`f4sat.c:673-677`) — replace that hack with two
explicit column maps in the new `res_block_t`.

`f4sat.c:632` already contains the codebase's only explicit degree-driver loop with state
(`next_deg`) persisted across iterations — a useful template for the slanted-degree loop.

### Rank extraction

For each (i, multidegree j) the degree-0 part of d_i is a scalar matrix over 𝔽_p — rows the
generators of F_i in degree j, columns those of F_{i-1} in degree j, entries the
coefficients of the terms whose ring monomial is 1 — and the formula above applies.
These matrices are **dense-ish and mutually independent across (i, j)** — the
natural GPU/NPU target. A rank-only kernel can **skip back-substitution entirely**, unlike
`exact_dense_linear_algebra_ff_32` (`la_ff_32.c:3390`), whose interreduction pass
(`la_ff_32.c:3354`) exists only to produce RREF. Note there is currently no function named
`rank` anywhere in `src/neogb` — the de-facto rank is `mat->np`, the count of non-NULL pivot
slots (`la_ff_32.c:3644-3665`).

### Pluggable linear algebra

Two separate seams, deliberately:

**(a) F4's own LA** — for completeness only. Add `int32_t la_backend` beside `laopt`
(`data.h:394`), extend `dispatch_linear_algebra` (`io.c:878`), and relax the `{1,2,42,44}`
whitelist in `validate_input_data` (`io.c:700-704`). Gate on `trace_level == NO_TRACER`,
since a device backend will not maintain `mat->rba` / `row[MULT]` / `row[BINDEX]`.

**(b) The resolution engine's LA** — the one that matters. Keep the vtable **narrow**; only
two operations are hot, and both have shapes known before arithmetic begins:

```c
typedef struct res_la_backend_t {
  const char *name;
  int (*available)(void);                                  /* runtime probe */
  int (*reduce_block)(res_block_t *blk, md_t *st);         /* rows vs. known pivots */
  int (*rank)(const res_dense_t *m, len_t *rank_out, md_t *st);
} res_la_backend_t;
```

Register backends in `res_la.c` and select by name/env. Implement the CPU reference first,
delegating to the existing `reduce_dense_row_by_known_pivots_sparse_*_ff_32` kernels.

One portability caveat to design around: msolve bakes `-mavx512f` / `-mavx2` at *configure*
time via an `AX_EXT` **build-host CPUID probe** (`m4/ax_ext.m4:234-243, 316-325`), so
binaries are not portable and kernel selection is a compile-time `#if`, not runtime
dispatch. The new backend registry should probe at **runtime** — otherwise "plug-n-play"
is a lie on any machine that is not the build machine.

---

## Status

| # | State |
|---|---|
| **M1** | **Done.** `src/neogb/res.h`, `res_grading.c` (degree groups over Z^r ⊕ T with an arithmetic vtable), `res_order.c` (Schreyer/POT/TOP), component slot + component-aware divisibility in `hash.c`, module dispatch in `io.c`. |
| **M2** | **Done.** Module-aware `get_lcm`, `insert_and_update_spairs` and `find_multiplied_reducer`; `res_module.c` with the `export_module_f4` C entry point taking a presentation matrix. Macaulay2 calls it through `rawMsolveModuleGB` in `Macaulay2/e/interface/msolve.{h,cpp}`. |
| **M3** | **Done.** `res_frame.c`: `res_frame_new` / `res_frame_init` / `res_frame_next_level` / `res_frame_complete`, frame ranks via `res_frame_betti`, and the `export_module_frame` C entry point. `res_module.c`'s validate-import-F4 preamble is now shared by both entry points as `module_gb_from_input`. |
| **M4** | **Done.** `res_diff.c`: `res_diff_new` / `res_diff_init` / `res_diff_compute` / `res_diff_verify`, one Macaulay matrix per (level, degree) with a parallel multiplier row. `export_module_resolution` in `res_module.c` covers both `RES_SYZ_OF_GB` and `RES_SYZ_OF_INPUT` and, with `max_level = 2`, is the single syzygy matrix. Cross-checked in Macaulay2 by `test/neogb/res/res_reference.m2`. |
| **M5** | **Done.** `res_betti.c`: `res_betti_new` / `res_betti_minimalize` / `res_betti_pdim` / `res_betti_reg` / `res_hilbert_invariants`, and the `export_module_betti` C entry point. Minimal Betti numbers by rank extraction, plus the Hilbert numerator (Macaulay2's `poincare`), projective dimension, regularity, Krull dimension and degree. Fixing the frame's block order is part of this milestone; see the notes. |
| **M7** | **Done.** `res_comp_t` in `res_module.c` is a resolution kept alive — `res_comp_new` / `res_comp_free` / `res_comp_nlevels` / `res_comp_rank` / `res_comp_degrees` / `res_comp_degshift` / `res_comp_is_complete` / `res_comp_differential`, with `res_diff_compute_thru` doing the lazy prefix. Macaulay2 drives it as an ordinary `ResolutionComputation` through `rawMsolveResolution` and the existing `rawResolutionGetFree` / `rawResolutionGetMatrix`, and the `Msolve` package wraps it as `msolveResolution`. Prune-to-minimal is `complex MsolveResolution` plus `Complexes`' own `minimize`; see the notes. |
| **M6** | **Done.** `res_grading_t` in `res.h` is the flat description a caller hands in and `res_dgrp_of_grading` turns into the engine's degree group; every entry point takes one and `NULL` is the standard grading. `ht->vwt` carries the variable weights so the Gröbner basis is computed in the *heft* order the grading induces, `res_dbkt_t` in `res_grading.c` buckets multidegrees through the group's own hash and comparison, and `res_betti.c` extracts ranks blocked by multidegree, reporting both the heft table and the multigraded one through `res_mtable_t`. See the notes. |
| **Orders** | **Done.** `res_strat_t` in `res.h` and `res_order.c`: base order (POT/TOP), component direction (up/down) and lift (Schreyer), threaded from the C entry points to `res_diff_cmp_mon`. `test/neogb/res/res_bench_strategy.c` measures them. The default is unchanged; see the notes for the measurement and why. |
| **M8, M9** | Not started. M9 is the option surface Macaulay2's `gb`/`syz`/`res` expect; see the notes. |

Verification in place: `neogb_res_selftest` (515 checks, run by `make check`), the
64 pre-existing diff tests still pass, and a cyclic-8 Gröbner basis is byte-identical
to the pristine 0.10.1 baseline with no measurable slowdown.
`test/neogb/res/res_reference.m2` is the Macaulay2 script every reference number in the
selftest came from; it also checks the two things C cannot check cheaply — that each
complex is *exact*, and that the syzygies generate the same module as `syz`.

M3's frame ranks agree entry for entry with Macaulay2's
`betti res(M, Strategy => Nonminimal)` on a nine-example corpus — monomial ideals,
the twisted cubic, a complete intersection of three quadrics where the frame
(1,6,8,3) is strictly larger than the minimal resolution (1,3,3,1), the
catalecticant cokernel, and a rank-two module with a nonzero degree shift.

The one free choice in the construction is the direction of the storage order
inside a block, and it is not cosmetic: Schreyer's theorem attaches the lead term
of an S-pair to the *larger* of the two indices, so the elements early in a block
are the ones being divided into, and their colon quotients are what the next level
is made of. ~~Ascending is what Macaulay2 does.~~ **Corrected in M5:** the order is
degree ascending and then the ring order **descending**, which is Macaulay2's
`sort(1, -1)` at level 1 (`res-f4-computation.cpp:132`) and its `PreElementSorter`
above (`res-schreyer-frame.cpp:32`, comparing varpower monomials, which list
variables from the highest down). `res_test_frame_block_order` is the discriminating
test — `(z², y²z, y³)` has the frame 1,3,3,1 descending and 1,3,2 ascending, and
Macaulay2 says 1,3,3,1. See the M5 notes for how the ascending version was caught.
The classical `(x², xy, y³)` does *not* discriminate; both directions give 1,3,2.

The subtlest bug found so far, and the reason the rank-two test earns its keep:
`find_multiplied_reducer` (`symbol.c:493`) **open codes** its divisibility test rather
than calling `check_monomial_division`, so the component guard added to the latter did
not apply there. A reducer from an earlier component was accepted, its multiplier
picked up a nonzero component, and each round produced monomials one component higher
than the last — the hash table grew to 131072 entries of the form `x·e_c` for
c = 0, 1, 2, …  Any future divisibility rule has to be added in **both** places.

Every check in `res_selftest.c` that concerns a module order is *discriminating*: it
fails when the module dispatch in `io.c` is disabled. This matters because the
component slot sits at the top of the reverse-lexicographic loop, so plain DRL breaks
ties by component all on its own and in the same direction — a naive test passes even
with the module code switched off. Verified by mutation.

### M4 notes

The reduction is where the two notions of syzygy finally separate in the code, and each
one taught something.

**The schedule is by degree, not by slanted degree**, and the correction is not
cosmetic — see the crossed-out paragraph above. The complete intersection of three
quadrics already exhibits it: the level 2 column of degree 3 reduces against a level 1
generator of degree 3 with multiplier 1, and under a slanted-degree driver that
generator's strand (`s = 2`) comes after the row that needs it (`s = 1`), in either
direction of the inner loop. Reordering the driver to slanted degree fails three
selftest checks. Degree ascending then level ascending is correct because a reducer's
degree never exceeds the row's, and `(i−1, d)` precedes `(i, d)`.

**Position over term is the only module order the differential supports.** The
Schreyer order the reduction runs in is the one that level 0's order induces, and at
level 0 that order *is* POT — `total` is 1 and the index is the component. Under term
over position the induced order would have to compare degrees that include the
component shifts, and the frame's ring hash table deliberately does not carry them.
Refused at the entry point rather than silently producing a non-complex.

**`RES_SYZ_OF_INPUT` returns a Gröbner basis of the syzygy module, not a minimal
generating set.** The graph-module trick computes in `R^nr_rows ⊕ R^nr_gens` under an
order the syzygy module did not choose, so its reduced GB can carry redundant
generators: the twisted cubic gets 3 columns where M2's `syz` gives 2, three quadrics 4
where M2 gives 3. `res_reference.m2` checks the thing that is actually true —
`image ours == image (syz f)` — for every corpus example. When the input generators are
already a reduced Gröbner basis *and* already minimal syzygy generators, the two agree
term for term; the Koszul complex and the catalecticant are checked that way.
Minimalizing is not a patch here, and M5 does not fix it either: rank extraction
minimalizes Betti *numbers* without ever building a minimal generating set,
so `RES_SYZ_OF_INPUT` still returns the Gröbner basis it always did.

**The exact `d ∘ d = 0` check does not run by default.** It is a polynomial
multiplication per column and measured at six times the cost of the resolution itself
(118 s versus 20 s on six generic cubics in six variables, 6.4 M terms), so
`res_diff_verify` takes a `deep` flag and `export_module_resolution` a `verify`
parameter. The cheap half — every column leads with the frame's monomial, with
coefficient 1, over the frame's parent — always runs, and is what caught the
reducer-selection mutation below. `make check` passes `verify = 1` throughout.

**Non-homogeneous input is now refused** by both the frame and the resolution entry
points. `st->homogeneous` is already computed by `import_module_input_data` and already
accounts for the component shifts; without the guard an inhomogeneous ideal produced a
confusing "a reducer is not monic" failure from deep inside the linear algebra.

Mutation testing: wrong sign on the multipliers fails 9 checks, taking the
*largest*-indexed admissible reducer instead of the smallest fails 6, and the
slanted-degree schedule fails 3. The reducer one is the subtle rule — the frame
guarantees an admissible reducer of index below `up(k)` exists precisely because `m_k`
is a minimal generator of the colon ideal, and picking any other one moves the lead
term off the frame.

### M5 notes

Rank extraction itself is the easy half. The formula holds as written, the blocks are
small and independent, and the elimination stops at row echelon form — no back
substitution, since only `rank` is wanted. `res_ech_rank` keeps a dense accumulator over
the columns and stores pivot rows sparsely in one arena reset per block, so fill-in is
paid for only where it happens.

**The Hilbert numerator is free.** `K_d = Σ_i (−1)^i frame_{i,d}` and
`K_d = Σ_i (−1)^i β_{i,d}` are the same number, because the rank corrections telescope,
so the numerator needs no field arithmetic at all past the Gröbner basis —
`res_betti_new` fills it in from the frame alone and `export_module_betti` with
`minimal = 0` never touches a differential. That the two agree is also the sharpest
cheap check there is on the rank extraction, since it fails the moment a correction
lands at the wrong level or degree; the selftest recomputes the sum from the minimal
table and compares it against the one the engine got from the frame.

Krull dimension and degree come out of the numerator by dividing out `(1 − t)`, which
is a running sum, until the value at 1 is nonzero: `dim = nv − c` and `degree = G(1)`.
Both are invariant under the degree shift, which only multiplies `K` by a power of `t`.

**The frame's block order was wrong, and M5 is what caught it.** Randomised
cross-checking against Macaulay2 — twenty ideals, not the hand-picked corpus — turned up
four whose `poincare` disagreed. The numerator is an alternating sum over the whole
frame, so disagreeing on it means the frame is not a resolution of the thing it claims
to resolve, which is a far louder signal than any Betti number. Two separate defects
were behind it:

1. **The within-block sort ran ascending instead of descending** (above). Every example
   in the M3 corpus is insensitive to the direction, which is why nine of them agreed
   with Macaulay2 anyway.
2. **`max_level = 0` meant `nr_vars`, and `nr_vars` is not enough.** Hilbert's syzygy
   theorem bounds the *minimal* resolution, not the frame. `(z, y², x²y, x³)` in three
   variables has the frame 1,4,6,4,1 — Macaulay2 agrees once its own `LengthLimit` is
   raised past its default of `nv`, which hides the same thing. The ceiling silently
   truncated the frame and the missing tail is exactly what made the alternating sums
   come out wrong. `res_frame_t` now grows its level array and `max_level = 0` means no
   ceiling; `res_frame_is_complete` no longer has an "nv is always enough" escape.
   The real bound is that the frame is a subcomplex of the Taylor complex on the
   Gröbner basis, so it stops by `lv[1].ld` levels; that is asserted rather than
   preallocated, since it is an enormous overestimate.

With both fixed, twenty random ideals in three to six variables agree with Macaulay2 on
frame ranks, `minimalBetti`, `poincare`, `pdim`, `regularity`, `dim` and `degree`, and
so do seven hand-picked ones including two modules and a shifted presentation.

**`max_level` builds one level more than it reports.** `β_{i,d}` reads the rank of
`d_{i+1}` as well as of `d_i`, so a table truncated at level *L* is only exact if level
*L+1* was built. Macaulay2's `minimalBetti` does the same. Everything that is an
invariant of the whole module — the numerator, `pdim`, `reg`, `dim`, `degree` — is
refused outright on a truncated frame rather than reported as a bound.

**Scalars are reported in the caller's own degrees, tables are not.** `betti` and
`hilbnum` are arrays indexed by degree and have to start at zero, so they keep the
shift that normalises `row_degs`, with `degshift` reporting it. `reg` has the shift
added back; `pdim`, `dim` and `degree` do not depend on it. Only regularity can tell
the difference, and it did — the shifted-presentation cross-check is what found it.

Mutation testing, all against `neogb_res_selftest`: dropping the `rank(d_{i+1})`
correction fails 10 checks, dropping `rank(d_i)` 13, counting every term instead of
only the constant ones 8, flipping the sign of the alternating sum 17, flipping the
block order tie break 4, dropping the degree key from the block order 3, and leaving
pivot rows un-normalised 1. That last one is the reason
`res_test_betti_generic_cubics` exists: every smaller example in the file has scalar
blocks whose pivots happen to be ±1, so the elimination gets the right rank without
ever dividing, and six random cubics in four variables is the smallest case found that
does not.

### M7 notes: a resolution kept alive

Everything up to M5 is shaped as a single call — marshal in, run, marshal out,
free — and that is the wrong shape for the question Macaulay2 actually asks. It
asks for the *shape* of a resolution first: `rawResolutionGetFree` wants the
rank and the degrees of F_i and nothing else. Only later, and only for the
levels it turns out to care about, does it call `rawResolutionGetMatrix`.
Answering the first question by materializing the whole complex defeats the
point on exactly the inputs a resolution engine exists for.

**The split is already in the mathematics.** The frame is combinatorial — no
field arithmetic past the Gröbner basis — and it determines every free module in
the resolution. The differential is where the linear algebra is. So
`res_comp_new` builds the Gröbner basis and the whole frame, and from that
moment `res_comp_rank` and `res_comp_degrees` are lookups; `res_comp_differential`
is the only thing that costs anything.

**A prefix is the only truncation that makes sense.** The block at level i in
degree d reduces against level i−1 in degrees up to and including d, so D_i
cannot be had without all of D_2 … D_{i−1}. Within that constraint laziness is
free: the schedule is degree ascending and then level ascending, and restricting
it to a range of levels visits the same blocks in the same relative order, each
seeing the same already-reduced data below it. `res_diff_compute_thru(rd, L)`
therefore costs exactly what computing levels 2…L in one go would, and
`res_diff_compute` is just `res_diff_compute_thru(rd, nlv-1)`. `rd->thru` records
how far it got, and asking again for something below it returns immediately.
The selftest checks this the hard way: it asks for the levels from the top down,
which makes the driver fill in a prefix it was never asked for, and requires the
result to match `export_module_resolution` term for term.

**The Gröbner basis does not survive construction.** Level 1 of the differential
is the only thing that ever reads it — `res_frame_init` takes the lead terms and
`res_diff_init` the coefficients — so both run eagerly and the basis and its hash
table are released before the handle is returned. `res_diff_init` is O(the
basis), and it buys the caller the right to hold a large resolution without also
holding the Gröbner basis it came from, which on the inputs this exists for is
the bigger of the two.

**Lifetime across the garbage collector.** This is the one place in the interface
where msolve state outlives a call, so it is the one place a lifetime has to be
tracked. Macaulay2's `ResolutionComputation` derives from `our_gc_cleanup`, so
`intern_res` installs a finalizer that deletes the object, and the destructor is
what calls `res_comp_free`. Nothing msolve allocated is ever handed to the
collector: a differential's flat arrays are copied into Macaulay2 objects and
released inside `get_matrix`.

**Interrupts are recoverable here**, unlike everywhere else in this interface. A
`SIGINT` inside `res_comp_differential` longjmps out mid-`res_diff_compute_thru`,
leaving levels at or below `rd->thru` complete and untouched and levels above it
partly filled. Asking again recomputes from `rd->thru + 1` straight over the top
of that, which is correct because a block only ever reads levels below it.
`res_dpoly_alloc` clears its slot first so the second pass does not leak the
first one's columns. The handle itself survives, so a Ctrl-C costs the work since
the last completed level and nothing else.

**Prune to minimal is Macaulay2's job, and that is the point.** The milestone
asks for the minimal complex, and the cheapest way to get it is not to write a
minimalizer in C: `complex MsolveResolution` in the `Msolve` package
materializes every differential and hands back an ordinary `Complex`, and
`Complexes`' `minimize` does the rest. That method is deliberately the one
place the laziness is given up — `C_i` and `length C` stay free, and nothing
else calls it implicitly.

What this buys is not a fast path but a *check*. Rank extraction (M5) reports
minimal Betti numbers without ever constructing a minimal generating set, which
makes it the one output in this engine with no independent witness on the C
side. Minimalizing the materialized complex is that witness, and it shares no
code with rank extraction at all. `res_reference.m2` now runs the three-way
comparison on every M7 example — `betti minimize complex C`, Macaulay2's
`minimalBetti`, and `unpackMsolveBetti rawMsolveMinimalBetti` — and all three
agree, including on the rank-two module with a shifted row, where msolve's
nonminimal resolution (2,2) and Macaulay2's (2,3,1) genuinely differ.

One caveat worth stating where someone will hit it: minimizing a complex
truncated by `LengthLimit` is right below the cut and wrong at it. The minimal
rank at level n depends on what cancels against level n+1, and a truncated
complex has not seen it.

**One thing to expect and not mistake for a bug**: the ranks reported here are
the frame's, and the frame is built from the Gröbner basis, which depends on the
module order. msolve resolves under position over term; Macaulay2's own
`res(…, Strategy => Nonminimal)` does not. So the two nonminimal resolutions of
the same module can genuinely differ — `coker {{x²,y²},{z,0}}` over
`k[x,y,z]` with rows in degrees 0 and 1 comes back 2,2 from msolve and 2,3,1
from Macaulay2, and both are right, the module being free of rank two. Only the
minimal Betti numbers are an invariant, and those are M5's, extracted from ranks
without any of this being materialized.

### Order strategies, and why the default did not change

The order a resolution runs in is now a **parameter**, `res_strat_t` in `res.h`,
carried from the C entry points down to the innermost comparison. Three axes:

| axis | values | where it acts |
|---|---|---|
| `base` | `RES_MORD_POT`, `RES_MORD_TOP` | the order on R^ncomp, i.e. what msolve's F4 computes the level 0/1 Gröbner basis under; `ht->mord` |
| `pos` | `RES_POS_DOWN`, `RES_POS_UP` | which component index counts as larger; `ht->mpos` |
| `lift` | `RES_LIFT_SCHREYER` | how a level above 0 is ordered — the induced Schreyer order, the only one implemented |

`NULL` means `res_strat_default()`, which is `pot-down-schreyer` — exactly what
the engine did before strategies existed, so the refactor is behaviour
preserving by construction. Extending this is meant to be *adding a field*, not
changing a signature; a user supplied weight matrix is the next one, which is
why the entry points take a pointer to the struct rather than a bare int.

**What had to change for TOP.** Only the comparison of two monomials of F_j,
`res_diff_cmp_mon`. Everything else in the frame — the colon quotients, the
block sort, the minimalization — only ever compares *within* a block, and
elements of one block share a parent and so a component, where POT and TOP agree
and the component's degree shift is a common constant. Across components they do
not agree, and there the frame's own hash table is not enough: `ev[DEG]` of a
module monomial includes the shift of its component (`res_module.c:52`), and
that is the degree the Gröbner basis was compared by, whereas the frame keeps
its monomials in a plain ring table whose `ev[DEG]` is the ring degree alone.
The shift is read back from level 0, where `res_frame_init` already stored it as
the heft degree of each generator of R^ncomp. That was the whole obstruction —
an unfinished feature, not a theorem.

The dispatch is a `switch` on a field that is constant for a whole computation,
not the function-pointer vtable `io.c` uses for the ring orders: it predicts
perfectly and still inlines into the sorts that call it.

**What the strategy must not change**, and what `res_test_strategies` checks
across the whole matrix on four inputs: the minimal Betti numbers entry for
entry, the Hilbert numerator, projective dimension, regularity, Krull dimension
and degree — and `d ∘ d = 0` exactly, under each. These are invariants of the
module; the frame ranks are not, and they differ by a factor. A broken
comparator gives a "Gröbner basis" that is not one and the ranks stop
telescoping, so this catches it immediately.

**The measurement, and the surprise.** `test/neogb/res/res_bench_strategy.c`
resolves a random corpus under every strategy. Three seeds, 50–60 random
homogeneous modules each (3–4 variables, rank 2–4, 3–6 generators, random degree
shifts):

| | frame (vs best) | differential terms | seconds | smallest frame |
|---|---|---|---|---|
| `pot-down-schreyer` | 2.07 / 2.62 / 2.33 | 1.40M / 1.50M / 2.45M | 1.88 / 3.61 / 4.83 | 9 / 5 / 4 |
| `pot-up-schreyer` | 2.14 / 2.36 / 2.18 | 1.39M / 1.38M / 2.20M | 1.86 / 3.32 / 3.61 | 5 / 3 / 5 |
| `top-down-schreyer` | **1.00** | 2.45M / 1.10M / 3.30M | 4.28 / 3.63 / 10.2 | **42 / 36 / 35** |
| `top-up-schreyer` | **1.00** | 2.43M / 1.11M / 3.31M | 4.20 / 3.62 / 9.90 | 4 / 6 / 6 |

Term over position gives a frame **about twice as small**, consistently, and
wins on frame size in roughly 80% of cases. It is also, on this corpus,
**slower** — because the columns of its differential are much denser: about 350
terms per generator against POT's 97. A smaller frame is not the same as less
linear algebra, which is precisely the thing that could not have been settled by
argument.

So the default stays `pot-down-schreyer`. The point of the milestone is that
changing it is now one line — `msolveResolutionStrategy()` in
`Macaulay2/e/interface/msolve.cpp`, or the `strat` argument in C — and that
there is a harness to justify the change with. Two caveats before anyone reads
the table as final: the corpus is small *dense* random modules, the regime where
frame size matters least and where the LiftTree family wins anyway (see the
caveats at the top), and none of it is threaded yet. The frame ratio is the
number likely to survive to real input; the wall clock is not.

For ideals the question is empty: with one component POT and TOP coincide, and a
rank-one run of the benchmark gives frame ratios of exactly 1.000 across all 30
cases.

### M6 notes: what a grading actually has to touch

The degree-group layer landed with M1 — `res_dgrp_t`, the arithmetic vtable, the
degree pools, and a full `res_deg_t` on every frame element. What M6 adds is
everything *between* that and the answer: nothing was reaching it, since every
entry point called `res_dgrp_new_standard`.

**How a caller says it.** `res_grading_t` — rank of the free part, torsion
orders, the `(r+nt) × nv` degree matrix, a heft — passed by pointer to every
entry point, `NULL` meaning the standard grading. Deliberately the same shape as
`res_strat_t`: extending it is adding a field, not changing a signature, which
is what the eventual non-abelian backend needs. `row_degs` generalizes with it,
carrying `res_grading_len(grading)` entries per row.

**Three places a grading has to reach, and only three.**

1. **`ev[DEG]` must be the heft degree.** This is what schedules the whole
   computation and what DRL compares first, so a grading that does not reach it
   is decoration. `ht->vwt` holds the per-variable weights and is `NULL` for the
   standard grading, so the unweighted path is provably the one it always was.
   Only two sites had to change: `set_module_exponent_vector`, and the module
   branch of `get_lcm` (`hash.c:1412`) — the one place in the engine that
   rebuilds a degree from exponents rather than adding two that were already
   right. Everything else is additive on slot 0 and so weight-agnostic for free.
   Weighted DRL is a term order exactly because `res_dgrp_new` insists the heft
   be strictly positive on every variable, which is the same condition that
   makes the degree-by-degree drivers terminate.

2. **The frame's own ring table needs the same weights.** This one is easy to
   miss and was: `res_frame_new` builds a *private* table for the frame's
   monomials, and `res_frame_colon` / `res_frame_init` / `res_diff_init` filled
   its degree slot with the plain total degree. Under a weighted grading that
   makes the Schreyer order the differential reduces under a *different* order
   from the one the level-1 Gröbner basis is a Gröbner basis for. `res_frame_hdeg`
   is now the single place that slot is computed.

3. **Rank extraction must block by multidegree, not by heft degree.** The
   differential is multihomogeneous, so its scalar part is block diagonal with
   respect to multidegree; the finer blocks are the diagonal blocks of the
   coarser ones. Blocking finer is therefore a *strict improvement* — smaller
   eliminations, and the multigraded table falls out — never a compromise, and
   the heft table stays exactly what it was because ranks add over a fibre.
   `res_dbkt_t` is the hash set that discovers the multidegrees, going through
   `grp->hash` and `grp->cmp`, which is precisely why those are in the vtable:
   an external group backend gets bucketing by supplying them.

**What did not have to change: `res_hilbert_invariants`.** Krull dimension and
degree are read off the *heft* numerator, and dividing by `1 - t^w` contributes
exactly one zero at `t = 1` whatever the weight `w`, so "write K = (1-t)^c·G,
then dim = nv − c and degree = G(1)" is already right. Checked against Macaulay2
on `k[x,y,z]` with degrees (1,2,3) — `dim 1`, `degree 18` — and on P¹×P¹, where
Macaulay2 collapses along the heft too.

**Homogeneity got stronger.** `st->homogeneous` compares heft degrees, which
under a multigrading is strictly weaker than being multihomogeneous: `x² + xy`
is homogeneous for Z and not for Z². `module_gb_from_input` now computes the
multidegree of every input term and the graded entry points refuse input that is
only heft-homogeneous. With the standard grading the two coincide, so nothing
changes there.

**Verification.** `res_test_weighted_*`, `res_test_multigraded_*`,
`res_test_torsion_betti`, `res_test_degree_buckets`,
`res_test_explicit_standard_grading` and `res_test_comp_multidegrees` — 81 new
checks, references in `res_reference.m2`. Each is *discriminating*, confirmed by
mutation: ignoring the grading entirely gives 18 failures, leaving `ev[DEG]`
unweighted 1, leaving the frame's table unweighted 1, never reducing torsion 4,
and bucketing by heft rather than multidegree 3.

The case that carries the milestone is P¹×P¹ with `J = (ac, bd, ad)`: level 2
holds two generators of the same heft degree 3 in multidegrees (1,2) and (2,1),
so the heft table reports a single `2` where the multigraded table reports two
`1`s. Macaulay2 has `multigraded betti` on a resolution it already built, but no
`minimalBetti` here — so on a large example this is doing something Macaulay2
cannot.

Two tests earn their keep beyond the obvious. `res_test_weighted_order` uses
`(y⁴ − x⁵z, z² − x²y²)`, weighted homogeneous of degrees 8 and 6 and *not*
homogeneous at all for the standard grading: the weights decide which term of
each generator leads, so dropping them makes the engine refuse the input rather
than answer it differently. And `res_test_explicit_standard_grading` requires the
standard grading spelled out as a `res_grading_t` to produce byte-identical
arrays to `NULL` — which is what makes every other test a test of the grading
rather than of a second code path.

**What the Macaulay2 binding now does with it.** The one-shot Betti and Hilbert
entry points — `rawMsolveMinimalBetti` and `rawMsolvePoincare` — hand a
`res_grading_t` over rather than inferring the grading from the order, so
`collectGrading` in `msolve.cpp` marshals the ring's own degrees and heft across
and the exponents go through untouched. Three things follow, all of which were
open when M6 landed:

- A multigraded ring is resolved as readily as a singly graded one, and there is
  no separate multigraded entry point, because there is no separate multigraded
  table: `rawMsolveMinimalBetti` reports the multidegrees themselves, and
  `unpackMsolveBetti` in `m2/betti.m2` reads them into keys `(i, d, h)` whose `d`
  is always a degree vector — of length one under a singly graded ring, which is
  what the slanted layout the engine's own resolutions use had been spelling as
  a single number all along. The heft-indexed table M2's own `minimalBetti`
  reports is one `applyKeys(B, (i,d,h) -> (i,{h},h), plus)` away, and `poincare`
  comes back in all `r` variables of the degrees ring.

  It is a *second* reader rather than a second case of `unpackEngineBetti`
  because the two layouts have nothing in common past the keys they build. The
  slanted one is dense and indexed by degree minus level, which only a degree
  that is a single number can be; a multidegree minus a level is not a
  multidegree, and the entries are sparse in the lattice — under the seven-fold
  grading of the `3*P` example, 10520 nonzero entries spread over 12712 buckets
  and nine levels, with every bucket that carries anything carrying it at
  exactly one level. So `rawMsolveMinimalBetti` sends the nonzero entries and
  nothing else, grouped by level:

      [r, len, n_0, ..., n_(len-1), (d_1, ..., d_r, heft, beta) x n_i per level i]

  which is 94691 ints for that example against 216108 for a record per bucket,
  and packs as a straight scan of `mbetti`, which is stored level-major already.
  What comes back is a plain `BettiTally`, printing folded onto the heft degree
  until the caller asks for `multigraded`, as for any other table.
- The exponent substitution is retired on those paths, so a weighted ring no
  longer approaches the 16-bit exponent ceiling through its exponents; only the
  degree accumulated from them is weighted, and that is checked.
- The monomial order stops mattering to them. msolve works in the heft order its
  grading induces, and minimal Betti numbers and the Hilbert numerator are
  invariants of the module rather than of the order, so a ring whose grevlex
  weights have nothing to do with its degrees is answered instead of refused.

What still infers the grading from the order is `rawMsolveResolution`, whose
result is a live `ResolutionComputation` handing back Macaulay2 `FreeModule`s one
at a time; `res_comp_multidegrees` reports what those would need, but building
them over a general degree monoid is work on the M2 side, so
`checkSinglyGradedByWeights` keeps that entry point singly graded for now.

### M9 notes: what the Macaulay2 option tables ask for

M5 finished the *outputs*. What is still entirely missing is *control over the
computation*, and the shape of the gap is easiest to read off Macaulay2's own option
tables for `gb`, `syz` and `res`. Between them the four module entry points accept
`max_level`, `syz_of`, `minimal`, `verify`, `ht_size`, `nr_threads`, `max_nr_pairs`,
`la_option`, `reduce_gb` and `info_level` — and nothing else.

The single structural fact behind most of the table: **msolve's F4 has no stop
conditions.** The round loop in `f4.c` runs a degree at a time until the pair set is
empty; there is no early exit, no partial result, and no way to ask for one. So
`BasisElementLimit`, `CodimensionLimit`, `PairLimit`, `StopWithMinimalGenerators`,
`SubringLimit` and `SyzygyLimit` are all the same piece of work — a stop-condition
mechanism in the round loop plus a way to report that the result is partial — not six
independent features. Note `max_nr_pairs` is *not* one of them: it sets `st->mnsel`,
the most pairs selected in a single round (`symbol.c:268`), which is a batching knob.
Mapping `PairLimit` onto it would silently return a complete basis computed slowly
instead of a partial one.

Ranked by value over cost:

1. **`DegreeLimit` / `HardDegreeLimit`.** The cheapest, because both loops are already
   degree-driven: `symbol.c` selects the pairs of minimal degree, and `res_diff_compute`
   is scheduled one `(level, degree)` at a time. A ceiling is a `break`, and the
   truncated-frame bookkeeping `export_module_betti` already does for `max_level`
   transfers directly.
2. **`Hilbert`.** M5 made msolve *produce* the numerator; teaching the F4 to *consume*
   one closes the loop. The classical win — knowing how many elements a degree must
   contribute and stopping the reduction once they appear — is exactly the case where
   msolve currently does the most redundant work, and it is the one option here that
   makes existing computations faster rather than merely stoppable.
3. **`ChangeMatrix`.** Nearly free, and already half implemented: the graph module
   `(f_j, e_j)` carries the expression of every Gröbner basis element in terms of the
   input in its adjoined block, and `RES_SYZ_OF_INPUT` computes that Gröbner basis and
   then discards precisely the elements a change matrix would keep. This is blocked on
   the same missing capability as the syzygy ordering (below), so the two should be done
   together.
4. **`SyzygyRows`.** Free as an output filter, since the retained rows are adjoined
   components; making it prune *work* rather than output needs those components ordered
   last, which is again the elimination-block question.

And the one that keeps coming back: `res.h` refuses to combine an elimination block
with a module order. That refusal is what forces `RES_SYZ_OF_INPUT` onto POT — `TOP`
compares `res_cmp_terms_drl` first (`res_order.c:147-153`), and for a homogeneous graph
module element every term shares `ev[DEG]`, so the degree ties and the lead can land in
the adjoined block while the original part is still nonzero, destroying the invariant
the whole trick rests on. Lifting the refusal — an elimination block over the
`R^nr_rows` components with TOP on the rest — fixes three things at once: `msolveSyzygy`
returns a basis for Macaulay2's default term-over-position-up order, `ChangeMatrix`
becomes available, and `SyzygyRows` can prune work.

`Macaulay2/packages/Msolve.m2` carries the full option-by-option table next to
`msolveDefaultOptions`, including the ones that are M2-side or meaningless for F4
(`Algorithm`, `MaxReductionCount`, `ParallelizeByDegree`, `SortStrategy`,
`StopBeforeComputation`).

## Milestones

| # | Deliverable | Validation |
|---|---|---|
| **M0** | Build + test harness; corpus of examples with reference Betti tables from M2 | `make check` green; corpus committed under `test/res/` |
| **M1** | `res_dgrp_t` degree groups (Z^r ⊕ T, vtable, r=1/T=0 fast path); module monomials in `ht_t` (comp slot, divisibility, orders, dispatchers) | Existing msolve test suite **bit-identical**; unit tests for degree arithmetic, module divisibility, Schreyer/POT/TOP comparison |
| **M2** | **Module F4** — GB of a submodule of a free module (`get_lcm`, `insert_and_update_spairs`) | Module GBs match M2 `gb` on module input; ideal case still bit-identical |
| **M3** | Schreyer frame construction (monomial only) + frame Betti table | Frame ranks match M2 `res(…, Strategy => Nonminimal)` |
| **M4** | Nonminimal differential (degree-by-degree driver) + **single syzygy matrix** output, both `SYZ_OF_GB` and `SYZ_OF_INPUT` | `d ∘ d = 0`; complexes exact in M2; syzygies generate the same module as M2 `syz` |
| **M5** | Rank extraction → `minimalBetti` equivalent, plus the Hilbert numerator and the invariants that fall out of it | Betti tables match M2 `minimalBetti` across the corpus, and `poincare` / `pdim` / `regularity` / `dim` / `degree` match on randomised input too |
| **M6** | Multigraded bucketing, then torsion in the degree group | Cross-check against M2 `res` + `betti` on multigraded examples (M2 has no `minimalBetti` here — this is the novel result); the heft-indexed table must be unchanged, and the standard grading passed explicitly must be byte-identical to passing none |
| **M7** | Materialize and export the full complex (flat-array C API, `mallocp` convention), and **incrementally**: a resolution kept alive, answering for the free modules from the frame and computing a differential only when one is asked for | Round-trip into M2; the incremental handle agrees term for term with the one-shot entry point, whatever order the levels are asked for; `prune` to minimal and compare |
| **M8** | LA backend vtable + CPU reference; then CUDA | Identical results across backends; benchmark sweep incl. sparse inputs where LiftTree should win |
| **M9** | **Computation controls** — degree ceilings, the Hilbert function hint, a change matrix, and the stop conditions Macaulay2's `gb`/`syz`/`res` options ask for | Each control reproduces the corresponding `gb(…, Option => v)` / `res(…, Option => v)` result in M2; the unconstrained path stays bit-identical |
| **Later** | Equivariant grading via an external group library → graded quotient rings (truncated, infinite frame) → exterior algebras (BGG/Tate) → FLINT fields → path algebras (Anick/Green–Solberg) | |

M1 and M2 are the risky ones — they touch msolve core. Keep the non-module path provably
unchanged, guarded on `ht->cpos != 0`.

---

## Verification

**Build.** `./autogen.sh && ./configure && make -j` from the repo root. FLINT, GMP and MPFR
are already hard dependencies (`configure.ac:25-41`), though `src/neogb/` currently uses
none of them — an `nmod_mat` backend would be the first FLINT use in neogb, and it links
already (`msolve.pc.in:10`).

**Regression guard for M1.** Run the existing suite under `test/neogb/` and `test/msolve/`
before and after, and diff outputs. The non-module path must be unchanged, not merely
passing.

**Cross-validation against Macaulay2.** This is the primary correctness signal. Use the
`/M2` skill to drive a live M2 session (much faster than `M2 --script`, since rings and
complexes persist):

```
I = ideal(...)                       -- corpus example
minimalBetti I                       -- reference
C = res(I, Strategy => Nonminimal)   -- reference nonminimal frame ranks
betti C
```

then compare against the new engine's output on the same input over the same prime. Build
the corpus from cases with known behavior: monomial ideals, generic determinantal ideals,
points in P^n, and (for M5) points in P^1 × P^1 and other toric examples where multigraded
Betti tables are documented.

**Composition check.** Independent of M2: verify `d_i ∘ d_{i+1} = 0` for every materialized
differential, and that frame ranks minus the two rank corrections are non-negative.
`res_diff_verify(rd, 1)` does the first exactly; `res_selftest.c` does it a second and
independent way, by substituting two points of 𝔽_p^n and multiplying the resulting scalar
matrices, which shares no code at all with the engine it checks. `res_betti_minimalize`
does the second, refusing to report a negative entry rather than clamping it.

The Hilbert numerator is a third such check and the cheapest of the three: it is the
alternating sum of the frame ranks *and* of the minimal Betti numbers, so computing it
both ways and comparing costs nothing and catches any rank correction that landed at the
wrong level or degree. It is also the check that caught the frame's block order — see
the M5 notes.

**Benchmarks.** Include sparse/monomial-flavored inputs where Singular's `fres` (LiftTree)
is expected to win — knowing where the matrix approach loses is part of the deliverable, and
it determines whether a hybrid strategy is worth adding later.

---

## Known landmines

- **Unity build ordering.** New `.c` files must be inserted into `gb.c:29-51` in dependency
  order. Nothing outside that list is compiled.
- **`src/neogb/la.h` is entirely dead** — 1302 lines inside `#if 0`, includes a nonexistent
  `types.h`, referenced by no file. It is nonetheless the only blocked / OpenMP-task-graph LA
  prior art in the tree (`la.h:772-780, 1103-1278`) and is worth reading before M7.
- **`insert_multiplied_signature_in_hash_table` (`hash.c:707`) computes the wrong hash.**
  It takes two hash *table indices* and multiplies their exponent vectors correctly, but
  sets the product's hash to `h1 + h2` — the sum of the two indices — instead of the
  linear hash `Σ rn[j]·a[j]` that `insert_in_hash_table` and every lookup path use.
  Entries it inserts are therefore invisible to every other insert, which silently adds
  a *second* copy of the same monomial under a different index, so nothing dedupes and
  two equal monomials compare unequal by index. Its only other caller is the dead SBA
  path (`sba.c:452`), which uses it for all inserts and lookups alike and so never
  notices. **Do not use it to multiply monomials.** `res_frame.c`'s `res_frame_mul`
  is the correct version. This bit the `total` field during M3 and was caught only by
  `res_frame_verify`, which is why that check runs on every frame rather than under a
  debug flag.
- **SBA path bugs**, if it is ever touched: `basis.c:209-251` (`check_enlarge_basis` never
  reallocs `sm`/`si`), `basis.c:404-416` (`copy_basis_mod_p` memcpys into unallocated
  `sm`/`si`), `update.c:386` (literal syntax error inside `#if 0`),
  `engine.c:137,154,170` (`exit(1)` on `use_signatures != 0`).
- **`md->max_gb_degree`** (`data.h:400`) is a working truncation hook checked at
  `symbol.c:231`, but it is hardwired to `INT32_MAX` at `f4.c:392` with a `TODO`. Useful for
  degree-truncated GBs at level 1; note `core_f4` then runs its finalization and frees
  `md->ps`, so it is not resumable as written.
- **`core_f4`'s loop body is duplicated six times in `modular.c`** (lines 203, 384, 670,
  1020-1060, 1148, 1475). Any refactor of the F4 driver must account for all six.
- **Aliasing in the LA layer**: `mat->rr` rows point into `bs->cf_32` (shared, not owned)
  while `mat->tr` rows own their coefficients via `mat->cf_32`, and
  `la_ff_32.c:2623-2628` rewrites `rr[i][COEFFS]`. Respect this on free.
