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

| Axis | Decision |
|---|---|
| Placement | **In-tree fork**, new `src/neogb/res_*.c` added to `gb.c`'s include list |
| Milestone 1 | **Both, staged** — one core (frame + differential), two consumers (ranks, complex) |
| Module monomials | **Component slot in `ht_t`/`ev`**, plus a lifted `total` monomial per frame element |
| Multigrading | **Designed in from day one**; degrees are an opaque **struct** over an arbitrary f.g. abelian group A = Z^r ⊕ T, heft-degree scheduling, r=1/T=0 specialized |
| Input | **Modules from day one**, not just ideals — an ideal is the rank-1 case |
| Outputs | **Single syzygy matrix** is a first-class output alongside full resolutions |

---

## Architecture

### Placement and build

`libneogb` is a **unity build**: `src/neogb/gb.c:29-51` `#include`s every `.c` file and
`src/neogb/Makefile.am:2` compiles only `gb.c`. Consequently nearly everything is `static`.
New files are added to that include list, after `la_ff_32.c` and before `engine.c`:

```c
#include "res_grading.c"   /* ZZ^r degrees, heft, multidegree buckets */
#include "res_frame.c"     /* Schreyer frame construction (monomial only) */
#include "res_la.c"        /* backend vtable + CPU reference kernels */
#include "res_diff.c"      /* nonminimal differential, slanted-degree driver */
#include "res_betti.c"     /* rank extraction, Betti table */
#include "res_export.c"    /* materialize complex, flat-array C API */
```

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

Drive by **slanted degree** `s = d − i` ascending (La Scala–Stillman), and within a slanted
degree by level. This is what bounds memory — it makes each level's data available exactly
when needed and lets finished strands be freed.

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

For each (i, multidegree j) the degree-0 part of d_i is a scalar matrix over 𝔽_p; apply the
formula above. These matrices are **dense-ish and mutually independent across (i, j)** — the
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
| **M3–M8** | Not started. |

Verification in place: `neogb_res_selftest` (108 checks, run by `make check`), the
64 pre-existing diff tests still pass, and a cyclic-8 Gröbner basis is byte-identical
to the pristine 0.10.1 baseline with no measurable slowdown.

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

## Milestones

| # | Deliverable | Validation |
|---|---|---|
| **M0** | Build + test harness; corpus of examples with reference Betti tables from M2 | `make check` green; corpus committed under `test/res/` |
| **M1** | `res_dgrp_t` degree groups (Z^r ⊕ T, vtable, r=1/T=0 fast path); module monomials in `ht_t` (comp slot, divisibility, orders, dispatchers) | Existing msolve test suite **bit-identical**; unit tests for degree arithmetic, module divisibility, Schreyer/POT/TOP comparison |
| **M2** | **Module F4** — GB of a submodule of a free module (`get_lcm`, `insert_and_update_spairs`) | Module GBs match M2 `gb` on module input; ideal case still bit-identical |
| **M3** | Schreyer frame construction (monomial only) + frame Betti table | Frame ranks match M2 `res(…, Strategy => Nonminimal)` |
| **M4** | Nonminimal differential (slanted-degree driver) + **single syzygy matrix** output, both `SYZ_OF_GB` and `SYZ_OF_INPUT` | `d ∘ d = 0`; syzygy matrix matches M2 `syz` |
| **M5** | Rank extraction → `minimalBetti` equivalent | Betti tables match M2 `minimalBetti` across the corpus |
| **M6** | Multigraded bucketing, then torsion in the degree group | Cross-check against M2 `res` + `betti` on multigraded examples (M2 has no `minimalBetti` here — this is the novel result) |
| **M7** | Materialize and export the full complex (flat-array C API, `mallocp` convention) | Round-trip into M2; `prune` to minimal and compare |
| **M8** | LA backend vtable + CPU reference; then CUDA | Identical results across backends; benchmark sweep incl. sparse inputs where LiftTree should win |
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
