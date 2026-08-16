-- Reference leading terms for test/neogb/order/block_order.c.
--
-- Every expected leading ideal in that test came from running this
-- script; rerun it after changing the test systems and paste the output
-- back into the reference tables.  Macaulay2 spells msolve's block
-- grevlex orders directly:
--
--   MonomialOrder => {GRevLex => k1, GRevLex => k2, ...}   -- blocks
--   MonomialOrder => {GRevLex => {w1, ..., wn}}            -- weights
--
-- Note that msolve orders its variables the same way, first block first,
-- so the variable lists below need no permutation.

p = 32003

-- print the leading terms of a Groebner basis as msolve exponent vectors:
-- one line per generator, exponents in variable order, space separated
printLeadExponents = (I) -> (
    R := ring I;
    n := numgens R;
    L := flatten entries leadTerm gens gb I;
    -- sorted so the output does not depend on the order gb returns them in
    E := sort apply(L, m -> first exponents m);
    << "  /* " << #E << " leading terms */" << endl;
    for e in E do << "    {" << demark(",", apply(e, toString)) << "}," << endl;
    )

-- The test system: cyclic 4 over F_p.
cyclic4 = (R) -> (
    x := gens R;
    ideal(
        x_0 + x_1 + x_2 + x_3,
        x_0*x_1 + x_1*x_2 + x_2*x_3 + x_3*x_0,
        x_0*x_1*x_2 + x_1*x_2*x_3 + x_2*x_3*x_0 + x_3*x_0*x_1,
        x_0*x_1*x_2*x_3 - 1)
    )

<< "=== cyclic4, 1 block (plain grevlex)" << endl
R1 = ZZ/p[x_0..x_3, MonomialOrder => {GRevLex => 4}];
printLeadExponents cyclic4 R1

<< "=== cyclic4, 2 blocks {1,3}" << endl
R2 = ZZ/p[x_0..x_3, MonomialOrder => {GRevLex => 1, GRevLex => 3}];
printLeadExponents cyclic4 R2

<< "=== cyclic4, 2 blocks {2,2}" << endl
R3 = ZZ/p[x_0..x_3, MonomialOrder => {GRevLex => 2, GRevLex => 2}];
printLeadExponents cyclic4 R3

<< "=== cyclic4, 3 blocks {1,1,2}" << endl
R4 = ZZ/p[x_0..x_3, MonomialOrder => {GRevLex => 1, GRevLex => 1, GRevLex => 2}];
printLeadExponents cyclic4 R4

<< "=== cyclic4, 4 blocks {1,1,1,1} (this is lex)" << endl
R5 = ZZ/p[x_0..x_3,
    MonomialOrder => {GRevLex => 1, GRevLex => 1, GRevLex => 1, GRevLex => 1}];
printLeadExponents cyclic4 R5

<< "=== cyclic4, 1 weighted block, weights {1,2,3,4}" << endl
R6 = ZZ/p[x_0..x_3, MonomialOrder => {GRevLex => {1,2,3,4}}];
printLeadExponents cyclic4 R6

<< "=== cyclic4, 2 weighted blocks {2,2}, weights {1,2,3,4}" << endl
R7 = ZZ/p[x_0..x_3,
    MonomialOrder => {GRevLex => {1,2}, GRevLex => {3,4}}];
printLeadExponents cyclic4 R7

-- print every term of a polynomial as msolve exponent vectors, in the
-- order the ring sorts them.  NOT sorted afterwards: the sequence is the
-- point.
printTermExponents = (f) -> (
    E := apply(terms f, m -> first exponents m);
    << "  /* " << #E << " terms, in order */" << endl;
    for e in E do << "    {" << demark(",", apply(e, toString)) << "}," << endl;
    )

-- A single dense form: for a principal ideal the Groebner basis is the
-- generator itself, so msolve exports it with its terms in monomial order.
-- Listing all of them makes this a full permutation test of the comparison
-- function, which no ideal-level check can be -- and unlike a leading
-- ideal it cannot collapse between block splits.
dense = (R) -> (
    x := gens R;
    ideal(x_0*x_1*x_2 + x_0^2*x_3 + x_1^3 + x_2*x_3^2
        + x_0*x_3^2 + x_1*x_2*x_3 + x_2^2*x_3 + x_0*x_1*x_3)
    )

for spec in {
    {"1 block", {GRevLex => 4}},
    {"2 blocks {1,3}", {GRevLex => 1, GRevLex => 3}},
    {"2 blocks {2,2}", {GRevLex => 2, GRevLex => 2}},
    {"2 blocks {3,1}", {GRevLex => 3, GRevLex => 1}},
    {"3 blocks {1,1,2}", {GRevLex => 1, GRevLex => 1, GRevLex => 2}},
    {"3 blocks {2,1,1}", {GRevLex => 2, GRevLex => 1, GRevLex => 1}},
    {"4 blocks {1,1,1,1}", {GRevLex => 1, GRevLex => 1, GRevLex => 1, GRevLex => 1}},
    {"1 weighted block {1,2,3,4}", {GRevLex => {1,2,3,4}}},
    {"2 weighted blocks {2,2} w {4,3,2,1}", {GRevLex => {4,3}, GRevLex => {2,1}}}
    } do (
    << "=== dense form, " << spec#0 << endl;
    S := ZZ/p[x_0..x_3, MonomialOrder => spec#1];
    printTermExponents first flatten entries gens gb dense S;
    )
