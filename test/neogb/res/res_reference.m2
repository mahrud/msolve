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
-- Minimalization is M5.
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
