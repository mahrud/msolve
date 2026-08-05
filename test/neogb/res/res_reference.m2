-- Macaulay2 reference for test/neogb/res/res_selftest.c.
--
-- Every rank, degree and syzygy matrix hard coded in res_selftest.c came
-- from running this file, over the same prime the selftest uses.  Run it
-- with
--
--     M2 --script test/neogb/res/res_reference.m2
--
-- The selftest itself needs no Macaulay2: it checks ranks and degrees
-- against the numbers printed here and verifies d o d = 0 by substituting
-- random points.  What it cannot check in C, and what this file is for,
-- is that the complexes are *exact* and that the syzygies generate the
-- same module Macaulay2's syz does.  Reproducing those two takes the
-- resolutions msolve computed; emit them from C as matrices over the same
-- ring and feed them to the `verify` function at the bottom.

needsPackage "Complexes";

p = 32003;
R = ZZ/p[x,y,z,w];
Q = ZZ/p[x,y,z];

-- ---------------------------------------------------------------------
-- Frame ranks: what export_module_frame reports, and what the resolution
-- must reproduce level by level and degree by degree.  Macaulay2 prints
-- these with the rows indexed by the slanted degree, degree minus level.
-- ---------------------------------------------------------------------

use Q;
print "-- Koszul (x,y,z): res(..., Nonminimal)";
print betti res(ideal(x,y,z), Strategy => Nonminimal);

use R;
Itc = ideal(x*z-y^2, x*w-y*z, y*w-z^2);
print "-- twisted cubic";
print betti res(Itc, Strategy => Nonminimal);

use Q;
Jq = ideal(x^2+y*z, y^2+x*z, z^2+x*y);
print "-- three quadrics: minimal is 1,3,3,1 but the frame is 1,6,8,3";
print betti res(Jq, Strategy => Nonminimal);
print betti res Jq;

use R;
A = matrix{{x,y,z},{y,z,w}};
print "-- catalecticant cokernel";
print betti res(coker A, Strategy => Nonminimal);

B = matrix{{x^2, y^2}, {z, w}};
print "-- coker {{x2,y2},{z,w}}, with the second row in degree one";
print betti res(coker(map(R^{0,-1}, R^{-2,-2}, B)), Strategy => Nonminimal);

use Q;
print "-- m^3 in three variables";
print betti res(ideal flatten entries basis(3, Q), Strategy => Nonminimal);

-- ---------------------------------------------------------------------
-- Syzygies.  Two different objects, and the selftest keeps them apart:
--
--   RES_SYZ_OF_GB     the Schreyer syzygies of the Gröbner basis, which
--                     is the second differential of the resolution above
--   RES_SYZ_OF_INPUT  syz of the generators the caller wrote down
--
-- They agree exactly when the input generators already are a reduced
-- Gröbner basis and already are minimal generators of the syzygy module.
-- The Koszul complex and the catalecticant are such cases and are
-- checked against these matrices term for term; the twisted cubic is not,
-- and there msolve returns a Gröbner basis of the syzygy module, which
-- has one redundant element more than Macaulay2's minimal answer.
-- M5 minimalizes Betti *numbers*, by rank extraction; it does not
-- minimalize a generating set, so that difference is still there.
-- ---------------------------------------------------------------------

use Q;
print "-- syz of (x,y,z): the three Koszul relations";
print syz gens ideal(x,y,z);

use R;
print "-- syz of the catalecticant: the signed maximal minors";
print syz A;

print "-- syz of the twisted cubic generators: Macaulay2 gives 2 columns,";
print "-- msolve's graph module gives a Gröbner basis with 3";
print betti syz gens Itc;
print syz gens Itc;

use Q;
print "-- syz of the three quadrics: Macaulay2 gives 3, msolve gives 4";
print betti syz gens Jq;

-- ---------------------------------------------------------------------
-- Checking a resolution msolve actually computed.
--
-- D is the list of differentials d_1, d_2, ... as Macaulay2 matrices over
-- the same ring, in the order export_module_resolution reports them, and
-- I is the module the complex is supposed to resolve.  Degrees are
-- dropped first: msolve reports them in its own array and they are
-- checked separately, so leaving them on the matrices would only make
-- Macaulay2 refuse to compose maps whose free modules it graded
-- differently.
-- ---------------------------------------------------------------------

strip = m -> map((ring m)^(numrows m), (ring m)^(numcols m), entries m);

verify = (name, D, I) -> (
    E := apply(D, strip);
    C := complex E;
    print("== " | name);
    print("   d o d = 0                : " |
        toString all(#E-1, i -> E#i * E#(i+1) == 0));
    print("   higher homology vanishes : " |
        toString all(toList(1..#E), i -> prune HH_i C == 0));
    print("   resolves the right thing : " | toString (coker E#0 == I));
    print("   ranks                    : " |
        toString prepend(numrows E#0, apply(E, numcols)));
    );

-- and, for the syzygy flavours, that the module is the right one even
-- when the generating set is not minimal
sameSyzygies = (name, S, f) -> print("   " | name | " : " |
    toString (image strip S == image strip syz f));

-- ---------------------------------------------------------------------
-- Minimal Betti numbers and Hilbert information.
--
-- What export_module_betti reports, in the same (level, degree, value)
-- form the frame tables above use.  The interesting rows are the ones
-- where this differs from the frame: three quadrics minimalize 1,6,8,3
-- down to 1,3,3,1 and the catalecticant 2,6,5,1 down to 2,3,1.
--
-- poincare is the numerator of the Hilbert series over (1-T)^nv, which
-- the engine gets from the *frame* alone -- no field arithmetic at all --
-- and which the selftest re-derives from the minimal table as a check
-- that the rank corrections landed at the right level and degree.
-- ---------------------------------------------------------------------

report = (name, M) -> (
    print("== " | name);
    print("   minimalBetti : " | toString new HashTable from minimalBetti M);
    print("   poincare     : " | toString poincare M);
    print("   pdim         : " | toString pdim M);
    print("   regularity   : " | toString regularity M);
    print("   dim          : " | toString dim M);
    print("   degree       : " | toString degree M);
    );

use Q;
report("koszul (x,y,z)", Q^1/ideal(x,y,z));
report("(x2,xy,y3)", Q^1/ideal(x^2,x*y,y^3));
report("three quadrics", Q^1/Jq);
report("m^3 in 3 vars", Q^1/ideal flatten entries basis(3,Q));
report("degree 4 monomials in x,y", Q^1/ideal(x^4,x^3*y,x^2*y^2,x*y^3,y^4));

use R;
report("twisted cubic", R^1/Itc);
report("catalecticant cokernel", coker A);
report("shifted coker {{x2,y2},{z,w}}",
    coker(map(R^{0,-1}, R^{-2,-2}, B)));

-- Six random cubics in four variables.  Every smaller example has scalar
-- blocks whose pivots are plus or minus one, so the rank extraction gets
-- away without normalizing its pivot rows; this one does not, and is the
-- only check in res_selftest.c that catches that.
S = ZZ/p[x_1..x_4];
Icub = ideal(
     3277*x_2*x_3^2 + 10825*x_1*x_2*x_4 + 23704*x_2^2*x_4,
    22284*x_1*x_4^2 + 19561*x_1*x_2*x_3 + 23260*x_1*x_2^2,
    19176*x_1*x_3*x_4 + 20404*x_2*x_3*x_4 + 27057*x_1*x_4^2,
     3298*x_1^2*x_2 + 4024*x_1*x_3^2 + 8445*x_1^2*x_4,
     1838*x_1*x_3^2 + 10295*x_1*x_2*x_3 + 1669*x_1^2*x_4,
    12895*x_3*x_4^2 + 14444*x_1*x_2*x_3 + 10377*x_1^2*x_3);
report("six generic cubics", S^1/Icub);
print("   frame        : " | toString new HashTable from
    betti res(S^1/Icub, Strategy => Nonminimal));

-- ---------------------------------------------------------------------
-- The frame is a *nonminimal* resolution, so Hilbert's syzygy theorem
-- does not bound its length.  Macaulay2 truncates at nv by default, which
-- hides this; raise LengthLimit and both engines agree that the frame of
-- (z, y^2, x^2 y, x^3) reaches level four in three variables.
--
-- The same examples pin down the one free choice in the frame, the order
-- of a block: degree ascending and then monomial order *descending*, per
-- res-f4-computation.cpp's sort(1, -1) and the PreElementSorter of
-- res-schreyer-frame.cpp.  Flipping it gives 1,3,2 for the first ideal
-- below and 1,4,6,4,1 for the third, neither of which is what Macaulay2
-- reports.
-- ---------------------------------------------------------------------

use Q;
for I in {ideal(z^2, y^2*z, y^3),
          ideal(z, y^2, x^2*y, x^3),
          ideal(z, y^3, x*y^2, x^2)} do (
    print("-- " | toString I);
    print betti res(I, Strategy => Nonminimal, LengthLimit => 8);
    );

-- ---------------------------------------------------------------------
-- A resolution kept alive (M7): res_comp_t, driven from Macaulay2 as an
-- ordinary ResolutionComputation through rawMsolveResolution and the
-- existing rawResolutionGetFree / rawResolutionGetMatrix.  This is what
-- the Msolve package's msolveResolution wraps.
--
-- What the C selftest cannot check is what a complex is *for*: that it is
-- exact, and that it resolves the module it was handed.  That is what
-- this section is for.  It also exercises the laziness -- the levels are
-- asked for from the top down, so the driver fills in a prefix it was
-- never asked for -- and requires that to agree with asking bottom up.
--
-- Expect the ranks to be the *frame's*, hence dependent on the Gröbner
-- basis and so on the module order.  msolve resolves under position over
-- term and Macaulay2 does not, so the two nonminimal resolutions of one
-- module can genuinely differ: coker {{x2,y2},{z,0}} with rows in degrees
-- 0 and 1 comes back 2,2 from msolve and 2,3,1 from Macaulay2, and both
-- are right -- the module is free of rank two, and the minimal Betti
-- numbers, 2,2, are what actually agree.
-- ---------------------------------------------------------------------

needsPackage "Msolve";

reportResolution = (name, M) -> (
    C := msolveResolution M;
    n := length C;
    frees := apply(n+1, i -> C_i);
    -- top down, which makes the differential driver fill in a prefix
    diffs := new MutableHashTable;
    scan(reverse toList (1 .. n), j -> diffs#j = C.dd_j);
    -- and again bottom up, which must be identical
    C2 := msolveResolution M;
    sameUpDown := all(1 .. n, j -> diffs#j == C2.dd_j)
        and frees == apply(n+1, i -> C2_i);
    print("-- " | name);
    print("   ranks        : " | toString apply(frees, numgens));
    print("   degrees      : " | toString apply(frees, F -> flatten degrees F));
    print("   d o d = 0    : " | toString all(2 .. n, j -> diffs#(j-1) * diffs#j == 0));
    print("   exact        : " | toString all(2 .. n, j -> ker diffs#(j-1) == image diffs#j));
    print("   resolves M   : " | toString (image diffs#1 == image M));
    print("   order free   : " | toString sameUpDown);
    print("   M2 nonminimal: " | toString new HashTable from
        betti res(coker M, Strategy => Nonminimal, LengthLimit => 8));
    );

use Q;  -- ZZ/p[x,y,z]
reportResolution("(z, y^2, x^2 y, x^3), whose frame runs past nv",
    gens ideal(z, y^2, x^2*y, x^3));
reportResolution("a rank two module with a shifted row",
    map(Q^{0,-1}, Q^{-2,-2}, {{x^2, y^2}, {z, 0}}));

use R;  -- ZZ/p[x,y,z,w]
reportResolution("the twisted cubic",
    gens minors_2 matrix {{x,y,z}, {y,z,w}});

-- LengthLimit truncates the frame, and the differentials below the cut
-- are the same ones the untruncated computation reports
Ctr = msolveResolution(ideal(z, y^2, x^2*y, x^3), LengthLimit => 2);
print("-- truncated at level 2: ranks " | toString apply(1 + length Ctr, i -> numgens Ctr_i));

-- ---------------------------------------------------------------------
-- Gradings by a finitely generated abelian group (M6).
--
-- Every number below is hard coded in res_selftest.c's res_test_weighted_*
-- and res_test_multigraded_* .  Two things are worth reading off:
--
--  * the *heft* table, which is what a caller that does not care about the
--    finer grading sees, and which stays what it always was -- a heft class
--    is a disjoint union of multidegree classes;
--  * the *multigraded* table, which is the novel output.  Macaulay2 will
--    print it with `multigraded betti`, but it has no minimalBetti here, so
--    on a large example msolve is doing something Macaulay2 cannot.
--
-- The discriminating case is the P^1 x P^1 one: level 2 carries two
-- generators of the same heft degree 3 in different multidegrees, (1,2)
-- and (2,1), so the heft table reports a single 2 where the multigraded
-- table reports two 1s.  Rank extraction that blocked by heft degree could
-- not produce the second table at all.
--
-- dim and degree need no multigraded generalization: they are read off the
-- heft numerator, and because dividing by (1 - t^w) contributes exactly one
-- zero at t = 1 per variable whatever the weight w, the same "write
-- K = (1-t)^c G, then dim = nv - c and degree = G(1)" recovers what
-- Macaulay2 reports for weighted and multigraded rings alike.  That is why
-- res_hilbert_invariants did not have to change for this milestone.
-- ---------------------------------------------------------------------

reportGraded = (name, M) -> (
    C := res M;
    print("-- " | name);
    print("   betti        : " | toString betti C);
    print("   multigraded  : " | toString (multigraded betti C));
    print("   poincare     : " | toString poincare M);
    print("   pdim/reg     : " | toString pdim M | " " | toString regularity M);
    print("   dim/degree   : " | toString dim M | " " | toString degree M);
    );

-- weighted: deg x = 1, deg y = 2, deg z = 3
Rw = ZZ/p[x,y,z, Degrees => {1,2,3}];
reportGraded("weighted (1,2,3), I = (x^2 y, y z, x z)",
    Rw^1/ideal(x^2*y, y*z, x*z));

-- the same ring, but an ideal whose *order* depends on the weights: these
-- two are weighted homogeneous of degrees 8 and 6 and are not homogeneous
-- at all for the standard grading, so the weights decide which term leads
Iw = ideal(y^4 - x^5*z, z^2 - x^2*y^2);
print("-- weighted-homogeneous? " | toString isHomogeneous Iw
    | ", standard-homogeneous? false");
print("   gb lead terms: " | toString leadTerm gens gb Iw);
reportGraded("weighted (1,2,3), I = (y^4 - x^5 z, z^2 - x^2 y^2)", Rw^1/Iw);

-- multigraded: P^1 x P^1
Tm = ZZ/p[a,b,c,d, Degrees => {{1,0},{1,0},{0,1},{0,1}}];
reportGraded("P^1 x P^1, J = (ac, bd, ad)", Tm^1/ideal(a*c, b*d, a*d));
reportGraded("P^1 x P^1, N = coker [ac  bd]", coker matrix{{a*c, b*d}});

-- Torsion has no Macaulay2 counterpart: its degree monoid is ZZ^r with no
-- torsion, so the reference for res_test_torsion_betti is by hand.  With
-- R = k[x,y] graded by ZZ (+) ZZ/2, deg x = (1,0) and deg y = (1,1), the
-- heft ignores the torsion and I = (x^2, xy) has both generators in heft
-- degree 2 -- but in multidegrees (2,0) and (2,1).  The syzygy
-- y*x^2 - x*(xy) sits in (3,1).  The heft table below is what Macaulay2
-- does see, over the ungraded ring, and the multigraded refinement of it
-- is what msolve adds.
use Q;
print("-- torsion reference: ZZ (+) ZZ/2 on k[x,y], I = (x^2, xy)");
print("   heft table (Macaulay2, standard grading): "
    | toString betti res (Q^1/ideal(x^2, x*y)));
print("   multigraded (by hand): level 1 in (2,0) and (2,1), level 2 in (3,1)");
