-- How long does the nonminimal differential take, and does it still give the
-- same answer?  A measurement rather than a test, and the one the steps in
-- doc/plan-res-perf.md are checked against.
--
-- Deliberately differential heavy: the Groebner basis and the whole Schreyer
-- frame together are under two seconds here, so essentially the entire run is
-- res_diff_compute_thru.  The frame is 1, 205, 922, 1698, 1584, 745, 141.
--
--   M2 --script test/neogb/res/res_bench_diff.m2
--
-- The Betti table it prints must be identical across every change to the
-- differential; the timing is what varies.
debug Core
needsPackage "Msolve"
R = ZZ/32003[a..h]
I = ideal random(R^1, R^{6:-3});
G = raw gens I;

print "-- Groebner basis and Schreyer frame only:";
elapsedTime rawMsolvePoincare(G, 1, 0);

print "-- differential and minimalization as well:";
elapsedTime w = rawMsolveMinimalBetti(G, 0, 1, 0);

print "-- the answer, which must not change:";
print unpackMsolveBetti w
