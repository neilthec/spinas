# Five-point QED Feynman diagrams with fixed photon references

This derivation concerns the all-plus photon amplitude implemented by
`SM/eeAAA.cpp`. It uses the supplied personal notes, *Massless Limit & Feynman*,
pp. 3-7, for the reference assignments and Weyl-block starting point.
Christensen et al., *Challenges with Internal Photons in Constructive QED*,
arXiv:2209.15018v3, Sec. I A and Appendix A, supplies the massive-photon-limit
context. The paper's four-point results do not themselves derive these six
five-point expressions.

## Conventions and the error

SPINAS stores physical, positive-energy momenta:

    p1 + p2 = k3 + k4 + k5
    Qa = p1 - ka,  Qc = p2 - kc
    D1a = Qa^2 - m^2,  D2c = Qc^2 - m^2.

The ordering `abc` runs from particle 1 to particle 2. The exchanged momenta
along the fermion arrow are `-Qa` and `Qc`. No internal-shell condition is
used. The ordinary external Dirac equations remain applicable.

The fixed spurious references in the notes are

    q3 = |4>, q4 = |5>, q5 = |3>.

The fixed physical references have angle bras

    <q3| = [3|p2, <q4| = [4|p1, <q5| = [5|p1.

These assignments are made once for the full amplitude. In particular,
`q4` remains the reference for photon 4 even if photon 4 is first or last in
the fermion ordering. The reference for photon 5 remains tied to photon 5.
Different photons may use different reference types.

The old code obtained the other diagrams by permuting the reference choices
with the photon ordering. For example, its `354` expression corresponds to
`q3=5, q5=4, q4=3`, rather than the fixed assignment above. The physical code
likewise moved the `p2` reference onto whichever photon became first.
This combines different external photon wavefunctions in one purported
Feynman sum. It is not an allowed use of permutation symmetry.
There was also an independent transcription error in the old physical `453`
term: its first product used `s15s` where the ordering required `s14s`.

Permuting full external states is allowed. To generate diagrams for one fixed
external state, however, permute only the vertex ordering and propagator
momenta. Keep the external wavefunctions fixed. Constructive diagram-local
numerator choices are a separate operation and are not modified here.

The formulas below are evaluated in SPINAS's physical-momentum convention.
Do not insert all-incoming photon spinors into them without also crossing
momenta and polarization tensors. The all-incoming `tau` denominators in
the notes become the `D` denominators above. The reference denominators must
also be crossed; their signs cannot be inferred from their names alone.

## Direct Weyl-block derivation

Let `P(v)` denote SPINAS's lower two-by-two momentum matrix and `Pbar(v)`
its adjugate (the upper matrix). Let `ai`, `si` be the right angle and right
square columns and `<i|`, `[i|` their corresponding rows.
For an angle-reference column `qi`, its dual row is `<qi|`. Define

    Ri = <qi i>,  R = R3 R4 R5
    Ei.lower = sqrt(2) |qi>[i| / Ri
    Ei.upper = sqrt(2) |i]<qi| / Ri
    Ei = [[0, Ei.lower], [Ei.upper, 0]].

The fixed external columns and rows in SPINAS conventions are

    U1 = (|1>, |1])^T,  V2 = (<2|, [2|).

For each ordering the Feynman amplitude is

    Mabc = -e^3 V2 Ec (slash(Qc)+m) Eb (-slash(Qa)+m) Ea U1
           / (D1a D2c),

where

    slash(Q) = [[0, P(Q)], [Pbar(Q), 0]].

This follows from the two fermion propagators with their arrow-directed
momenta and the three photon vertices. The overall sign is kept fixed in
both the matrix calculation and its scalar expansion. No fitted sign, scale,
constructive numerator, or on-cut replacement enters this calculation.

Multiplying from the right, suppressing the common photon normalizations,
gives two scalar coefficients after the middle vertex:

    Xab = -[b|Qa|qa>[a1] + m[ba]<qa1>
    Yab = m<qb qa>[a1] - <qb|Qa|a]<qa1>.

The remaining propagator and left photon vertex give

    Abc = <2 qc>[c|Qc|qb> + m[2c]<qc qb>
    Bbc = m<2 qc>[cb] + [2c]<qc|Qc|b].

Consequently the complete numerator is

    Nabc = Abc Xab + Bbc Yab
    Mabc = -2 sqrt(2) e^3 Nabc / (R D1a D2c).

The first index of `Abc`/`Bbc` is the middle photon `b`; the second is the
last photon `c`. Their `Qc` always means `p2-kc`. `Xab`/`Yab` always use
`Qa=p1-ka`. These are functions of labeled external references, not references
chosen anew for each ordering.

For direct comparison with the eight terms in the notes, set

    L=<2 qc>, S=[2c], T=[a1], C=<qa1>,
    U=[c|Qc|qb>, V=[b|Qa|qa>,
    W=<qb|Qa|a], Z=<qc|Qc|b],
    cb=[cb], ba=[ba], rr=<qb qa>, rc=<qc qb>.

Then

    Nabc = -L U V T + m L U ba C + m^2 L cb rr T - m L cb W C
           -m S rc V T + m^2 S rc ba C + m S Z rr T - S Z W C.

`amp_feynman` evaluates the original Weyl blocks. `amp_feynman_r` evaluates
these eight scalar terms separately with the same cyclic references. It is
an independent expansion of the diagrams, rather than an alias for the
compact constructive answer.

## All six amplitudes

All entries have the SAME `R`; there are no extra permutation signs to attach
to these expressions. Signs from bracket orientation are already included
in the definitions above.

| Diagram | Numerator | Physical propagators |
| --- | --- | --- |
| 345 | `A45 X34 + B45 Y34` | `D13 D25` |
| 354 | `A54 X35 + B54 Y35` | `D13 D24` |
| 435 | `A35 X43 + B35 Y43` | `D14 D25` |
| 453 | `A53 X45 + B53 Y45` | `D14 D23` |
| 534 | `A34 X53 + B34 Y53` | `D15 D24` |
| 543 | `A43 X54 + B43 Y54` | `D15 D23` |

Multiply every row by `-2 sqrt(2) e^3/R` and divide by the stated physical
propagators. Their sum is the full Feynman amplitude.

### Spurious assignment, explicitly

Set `h3=4`, `h4=5`, `h5=3`, once. The coefficients above become

    Xab = -[b|Qa|ha>[a1] + m[ba]<ha1>
    Yab = m<hb ha>[a1] - <hb|Qa|a]<ha1>
    Abc = <2 hc>[c|Qc|hb> + m[2c]<hc hb>
    Bbc = m<2 hc>[cb] + [2c]<hc|Qc|b]
    R = <43><54><35> = <34><35><45>.

In particular, the `354` diagram has

    X35 = -[5|Q3|4>[31] + m[53]<41>
    Y35 = m<34>[31] - <3|Q3|3]<41>
    A54 = <25>[4|Q4|3> + m[24]<53>
    B54 = m<25>[45] + [24]<5|Q4|5],

where `Q3=p1-k3` in `X35,Y35` and `Q4=p2-k4` in `A54,B54`.
Unlike a wholesale `4 <-> 5` swap of the original expression, photon 3 still
uses reference 4, photon 4 still uses reference 5, and photon 5 still uses
reference 3. The table and the same definitions explicitly specify every
other ordering without an ambiguous swap prescription.

### Physical assignment, explicitly

Set `r3=2`, `r4=1`, `r5=1`, once. The angle rows are `[i|p_ri`.
In SPINAS component conventions the corresponding columns are

    |qi>column = -P(p_ri)|i]column
    <qi|row    = [i|row Pbar(p_ri).

The minus sign is the index-raising sign needed to make these columns and
rows dual; it is not an extra diagram sign. This yields explicit coefficients
with no auxiliary spinors:

    Xab = [b|Qa p_ra|a][a1] + m[ba][a|p_ra|1>
    Yab = -m[b|p_rb p_ra|a][a1] - [b|p_rb Qa|a][a|p_ra|1>
    Abc = [c|p_rc|2>[c|Qc p_rb|b] - m[2c][c|p_rc p_rb|b]
    Bbc = -m[c|p_rc|2>[cb] + [2c][c|p_rc Qc|b].

The reference denominator is fixed for ALL SIX rows:

    Rphys = [3|p2|3>[4|p1|4>[5|p1|5>
          = -D23 D14 D15

in the stored physical-momentum convention, since `[i|pj|i>=2 ki.pj`
and `Dji=-2 pj.ki`. Thus each physical row may equivalently be written as

    Mabc.physical = +2 sqrt(2) e^3 (Abc Xab + Bbc Yab)
                    / (D23 D14 D15 D1a D2c).

The reference factor `D23 D14 D15` never moves with the ordering.
The six complete denominator products, retaining repeated factors, are:

| Diagram | Denominator after using `Rphys=-D23 D14 D15` |
| --- | --- |
| 345 | `D23 D14 D15 D13 D25` |
| 354 | `D23 D14 D15 D13 D24` |
| 435 | `D23 D14 D15 D14 D25` |
| 453 | `D23 D14 D15 D14 D23` |
| 534 | `D23 D14 D15 D15 D24` |
| 543 | `D23 D14 D15 D15 D23` |

Repeated denominator factors here come from a reference denominator coinciding
with a propagator. They are not additional internal lines. Keep the complete
numerator; cancellations in the full amplitude remove gauge artifacts.
The implementation evaluates the reference contractions directly to avoid
confusing reference and propagator signs.

## Implementation and use

The existing public names retain their signatures and all-plus-only scope:

- `amp_feynman`: direct matrix contraction, fixed cyclic references.
- `amp_feynman_r`: expanded eight-term contraction, fixed cyclic references.
- `amp_feynman_x`: direct matrix contraction, fixed physical references.

The added `feynman_diagrams` interface returns individual diagrams in order
`345,354,435,453,534,543`. For example:

```cpp
using Ref = spinas::eeAAA::feynman_reference;
const std::array<Ref,3> physical = {{{2,true},{1,true},{1,true}}};
const std::array<Ref,3> mixed    = {{{2,true},{5,false},{1,true}}};
auto diagrams = amplitude.feynman_diagrams(ds1, ds2, physical);
auto mixed_diagrams = amplitude.feynman_diagrams(ds1, ds2, mixed);
auto expanded = amplitude.feynman_diagrams(ds1, ds2, physical, true);
auto ward4 = amplitude.feynman_diagrams(ds1, ds2, physical, false, 4);
```

A reference entry `{j,false}` means `qi=|j>`, with a different massless
external leg. `{j,true}` means `<qi|=[i|pj`. Every entry remains fixed over
all six diagrams. `ward_leg=3,4,5` replaces that polarization by its momentum
in the matrix calculation. The sum of the resulting Ward diagrams vanishes;
individual diagrams generally do not.

An optional `evaluation_scales` output stores absolute matrix-product bounds
before cancellations. This is a numerical diagnostic, not an alternate
amplitude or an analytic error certificate.

## Numerical stability and validation

The original `sqrt(E-|p|)` in `source/particle.cpp` loses relative precision
when the electron mass is tiny compared with its energy. The on-shell identity

    sqrt(E-|p|) = m / sqrt(E+|p|)

avoids that subtraction for nonnegative masses. The zero-mass case is handled
explicitly, including zero energy. This is an external-spinor numerical fix;
it makes no internal-shell substitution. A dedicated test checks a mass/energy
ratio of `1e-9`, where the original subtraction rounds to zero on this host.

Feynman diagrams can still suffer severe cancellation at `m=0.0005`,
`sqrt(s)=600`, especially for the suppressed `ds1=ds2=+1` state. Near-zero
spin components also make relative errors meaningless. The regression uses
both a relative bound and a roundoff allowance proportional to the uncanceled
matrix-product scale. It reports which comparisons are cancellation limited.
It does not claim `1e-9` relative accuracy for such amplitudes.
The displayed ratios are now true complex ratios, rather than separate ratios
of real and imaginary parts.

Validation on the native arm64 host:

- Existing `H=S2` checks: 8,640 amplitude comparisons, zero failures.
- Fixed-reference Feynman/H checks: 6,912 comparisons, covering all fermion
  spins, all photon orders, cyclic/physical/mixed/reverse references, three
  frames, and three mass/energy settings; zero failures.
- Individual matrix/expanded diagram checks: 41,472; zero failures.
- Ward checks: 20,736; zero failures.
- Well-resolved Feynman/H comparisons: 3,648; maximum relative error
  approximately `4.43e-14`. The other 3,264 are reported as cancellation limited.
- Maximum error divided by epsilon times the uncanceled evaluation scale:
  approximately `6.73` (regression allowance: `256`).
- The independent SymPy script `docs/eeAAA-feynman-check.py` uses exact rational
  Weyl spinors at a conserved, external-shell point. It compares all diagrams
  with the eight-term formulas, checks the sum against H, and checks all Ward
  replacements, for all four fermion spins and four reference assignments.
  There are 184 exact equality assertions plus a negative control that rejects
  moving the cyclic references with the ordering. The explicitly substituted
  physical coefficients above are also checked diagram by diagram.
  These finite exact fixtures supplement the general block derivation; they
  are not a generic symbolic proof over phase space.

The preexisting `amp_permutation` is a different legacy S1 formula. It is not
used as the reference answer, and this change does not repair or endorse it.
The compact `amp_nonspurious` and `amp_spurious` functions remain independent
of the Feynman implementation.

The 33 matrix/vector/particle/propagator/spinor-product Boost test cases pass,
including the new small-spinor regression. The complete core suite additionally
reports four event-parser assertions; the same four fail with the unchanged
source (two mass-square tolerance assertions and two missing-exception checks).
The configured Boost dynamic library on this host is x86_64 while the compiler
is arm64, so native core checks used Boost's header-only runner without changing
the repository's CMake configuration. `test_SM` builds normally through CMake
and exits successfully. Its current main returns after `test_eeAAA`; this is
not a claim to have run the later SM process tests that main leaves unreachable.

For the twelve printed high-energy fixtures, the corrected Feynman/constructive
complex ratios differ from one by roughly `1e-5` to `5e-4`, rather than factors
of `1e10` to `1e12`. The residual is consistent with the conditioning of a sum
whose individual diagrams can be twelve orders larger than the answer.
