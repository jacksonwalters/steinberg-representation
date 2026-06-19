# Representations of finite general linear groups

This repository contains small SageMath computations of representations of
finite general linear groups. The original focus is the Steinberg
representation of `GL_3(F_q)`, with concrete scripts for `q = 2, 3` and a
parameterized script for general prime powers `q`. It also includes
projective-line and representation-family scripts for `GL_2(F_q)`.

The original `GL_3(F_2)` computation corresponds to the expository paper:
[https://jacksonwalters.com/docs/notes/steinberg_representation_GL(3,2).pdf](https://jacksonwalters.com/docs/notes/steinberg_representation_GL(3,2).pdf)

The representations are realized geometrically from the spherical building of
`F_q^3`. In this rank-2 case, the building is the incidence graph whose
vertices are:

- the 1-dimensional subspaces of `F_q^3`;
- the 2-dimensional subspaces of `F_q^3`;
- edges given by inclusions `line <= plane`.

The Steinberg representation is computed as the top homology of this building:

```text
St = H_1(Delta; C) = ker(partial: C_1(Delta; C) -> C_0(Delta; C)).
```

For `GL_3(F_2)`, the computation finds an 8-dimensional Steinberg module. For
`GL_3(F_3)`, the same construction gives a 27-dimensional Steinberg module. In
general this rank-2 construction gives `dim St = q^3`.

## Files

- `gl_3_q/steinberg_representation_gl_3_2.py` is the original `GL_3(F_2)`
  SageMath script, now with concise default output.
- `gl_3_q/steinberg_representation_gl_3_3.py` is the `GL_3(F_3)` SageMath
  script.
- `gl_3_q/steinberg_representation_gl_3_q.py` is the general `GL_3(F_q)`
  script. It accepts `q` as a command-line argument and keeps larger
  computations opt-in or automatic only for small enough inputs.
- `gl_3_q/fano_plane_representations_gl_3_2.py` specializes to `GL_3(F_2)` as
  the automorphism group of the Fano plane. It constructs the point, line, flag,
  reduced point/line, boundary-image, and Steinberg representations as explicit
  matrices, and records the two complex 3-dimensional irreducible characters
  over `Q(zeta_7)`.
- `gl_2_q/projective_line.py` contains reusable `P^1(F_q)` geometry and
  projective-line action helpers for `GL_2(F_q)`.
- `gl_2_q/steinberg_representation_gl_2_q.py` constructs the `GL_2(F_q)`
  Steinberg representation as reduced `H_0(P^1(F_q))`.
- `gl_2_q/representations_gl_2_q.py` organizes the standard complex
  representation families of `GL_2(F_q)`, constructs determinant twists of
  Steinberg, and constructs principal series by explicit induction from the
  Borel subgroup.
- `gl_2_q/whittaker_basis_gl_2_q.py` changes the cuspidal projector model into
  a Whittaker basis indexed by `F_q^*`, giving compact `(q - 1) x (q - 1)`
  cuspidal matrices while retaining the Gelfand-Graev construction as a
  correctness oracle.
- `gl_2_q/direct_whittaker_cuspidal_gl_2_q.py` constructs the same cuspidal
  Whittaker model directly from the nonsplit-torus character formula and Bruhat
  decomposition, avoiding the large Gelfand-Graev projector space.

## What the scripts do

The scripts:

1. Construct the points and planes in `F_q^3`.
2. Build the incidence edges of the spherical building.
3. Form the boundary map `C_1 -> C_0`.
4. Compute `ker(partial)`, giving the Steinberg representation.
5. Construct apartment cycles from bases of `F_q^3`.
6. Define or test the natural `GL_3(F_q)` action on the building and cycles.
7. Count apartment cycles up to orientation.
8. Compute the span of apartment cycles.
9. Compute or verify the character of the resulting representation.

The `GL_3(F_3)` and general `GL_3(F_q)` scripts compute the Steinberg character
from fixed vertex and edge counts instead of building a `q^3 x q^3` matrix for
every group element. For a connected graph,

```text
chi_St(g) = # fixed edges - # fixed vertices + 1.
```

## Running the computation

Install SageMath (tested with SageMath 10.8), then run from this directory:

```bash
sage gl_3_q/steinberg_representation_gl_3_2.py
sage gl_3_q/steinberg_representation_gl_3_3.py
sage gl_3_q/steinberg_representation_gl_3_q.py 3
sage gl_3_q/fano_plane_representations_gl_3_2.py
sage gl_2_q/steinberg_representation_gl_2_q.py 5
sage gl_2_q/representations_gl_2_q.py 5
```

For another prime power, pass `q`:

```bash
sage gl_3_q/steinberg_representation_gl_3_q.py 5
```

The general script suppresses bulky data by default. Useful options are:

```bash
sage gl_3_q/steinberg_representation_gl_3_q.py 3 --character-table
sage gl_3_q/steinberg_representation_gl_3_q.py 7 --apartment-span yes
sage gl_3_q/steinberg_representation_gl_3_q.py 7 --character yes
sage gl_3_q/steinberg_representation_gl_3_q.py 3 --weyl-orbits
sage gl_3_q/fano_plane_representations_gl_3_2.py --character-table
sage gl_3_q/fano_plane_representations_gl_3_2.py --incidence
sage gl_3_q/fano_plane_representations_gl_3_2.py --generators
sage gl_2_q/steinberg_representation_gl_2_q.py 5 --character-table
sage gl_2_q/representations_gl_2_q.py 5 --all-principal-series
sage gl_2_q/representations_gl_2_q.py 4 --cuspidal-parameters
sage gl_2_q/representations_gl_2_q.py 5 --cuspidal 1
sage gl_2_q/representations_gl_2_q.py 3 --all-cuspidals
sage gl_2_q/whittaker_basis_gl_2_q.py 5 --cuspidal 1 --generators
sage gl_2_q/direct_whittaker_cuspidal_gl_2_q.py 7 --cuspidal 1 --generators
```

## `GL_2(F_q)` scripts

For `GL_2(F_q)`, the spherical building is the projective line
`P^1(F_q)`. The Steinberg representation is the reduced zero-th homology:

```text
St = \widetilde H_0(P^1(F_q); C) = ker(sum: C[P^1(F_q)] -> C).
```

Since `#P^1(F_q) = q + 1`, this gives `dim St = q`, and the projective-line
permutation module decomposes as

```text
C[P^1(F_q)] = 1 + St.
```

The broader `GL_2(F_q)` representation script records the standard complex
irreducible families:

- determinant characters, dimension `1`;
- determinant twists of Steinberg, dimension `q`;
- principal series `Ind_B^G(chi_1 tensor chi_2)` with `chi_1 != chi_2`,
  dimension `q + 1`;
- cuspidal series from nonsplit-torus character orbits, dimension `q - 1`.

The script checks the family count by verifying that the sum of squares of the
listed dimensions is `|GL_2(F_q)|`. It constructs the projective-line module,
Steinberg twists, and principal series explicitly.

Cuspidal representations are constructed by a Gelfand-Graev projector. For a
regular character `theta` of the nonsplit torus

```text
T = F_{q^2}^*,
```

the script forms

```text
Gamma = Ind_U^G(psi),
```

where `U = {[[1, x], [0, 1]] : x in F_q}` and `psi` is a nontrivial additive
character of `F_q`. It then applies the central idempotent

```text
e_theta = (dim pi_theta / |G|) sum_g chi_theta(g^{-1}) Gamma(g).
```

The image of `e_theta` has dimension `q - 1`; restricting `Gamma(g)` to this
image gives explicit cuspidal matrices. This exact projector method is useful
for small `q`, but it uses a Gelfand-Graev space of dimension `|GL_2(F_q)|/q`,
so larger values of `q` become slower than the principal-series construction.

The Whittaker-basis script keeps this projector model as a source of truth and
then decomposes the cuspidal image under the upper-unipotent subgroup. For each
`y in F_q^*`, it applies the Fourier idempotent

```text
P_y = (1/q) sum_x psi(-yx) pi(n(x)),
```

where `n(x) = [[1, x], [0, 1]]`. Each image is one-dimensional, and these lines
form a basis in which

```text
pi(n(x)) e_y = psi(yx) e_y.
```

This gives a genuine Whittaker/Kirillov-style basis indexed by `F_q^*`, while
still checking the resulting traces against the nonsplit-torus cuspidal
character formula.

The direct Whittaker script avoids constructing `Gamma`. It uses the explicit
Borel action on basis vectors `e_y`, together with Bruhat decomposition

```text
g = n(a/c) diag(-det(g)/c, -c) w n(d/c)        if c != 0,
```

for `g = [[a, b], [c, d]]`. The Weyl matrix for
`w = [[0, 1], [-1, 0]]` is recovered by Fourier inversion from the identity

```text
tr(pi(diag(a, 1) w n(s))) = sum_y W[a*y, y] psi(y*s).
```

Those traces are supplied by the nonsplit-torus cuspidal character formula, so
the construction stays in dimension `q - 1` throughout.

## Fano-plane representations for `GL_3(F_2)`

Since `F_2^*` is trivial, `GL_3(F_2) = PGL_3(F_2)` acts faithfully on the Fano
plane. The Fano-plane script builds this action on:

- the 7 projective points;
- the 7 projective lines;
- the 21 incident point-line flags.

It then forms the same building boundary map used in the Steinberg scripts,

```text
partial: C[flags] -> C[points] + C[lines],
```

and verifies the equivariant decomposition

```text
C[flags] = im(partial) + St,
dim C[flags] = 21 = 13 + 8.
```

The point and line permutation representations each split as

```text
C[points] = 1 + chi_6,
C[lines]  = 1 + chi_6,
```

where `chi_6` is the irreducible 6-dimensional augmentation representation. The
flag representation has character decomposition

```text
C[flags] = 1 + 2*chi_6 + St.
```

The script checks these decompositions using character inner products over the
six conjugacy classes of `GL_3(F_2)`. It also adds the two conjugate
3-dimensional characters with order-7 values

```text
alpha = zeta_7 + zeta_7^2 + zeta_7^4,
alpha_bar = zeta_7^3 + zeta_7^5 + zeta_7^6,
```

and verifies the character identities

```text
Sym^2(chi_3) = chi_6,
chi_3 * chi_3bar = 1 + St.
```

## Expected output

Default output is concise: it reports counts, dimensions, apartment checks, and
character irreducibility checks without printing large basis vectors or full
character tables.

For `GL_3(F_2)`:

- number of points: `7`;
- number of planes: `7`;
- number of incidence edges: `21`;
- dimension of `ker(partial)`: `8`;
- number of unique apartments up to orientation: `28`;
- Weyl group size: `6`;
- number of Weyl-group orbits on apartments: `7`;
- orbit sizes: `[3, 6, 1, 6, 6, 3, 3]`;
- dimension of the `GL_3(F_2)`-orbit span of an apartment cycle: `8`;
- Steinberg character values by class order/size:
  `8, 0, -1, 0, 1, 1`;
- inner product of the Steinberg character with itself: `1`.

For the Fano-plane representation script:

- group order: `168`;
- Fano-plane counts: `7` points, `7` lines, `21` flags;
- building boundary rank: `13`;
- `dim St = 8`;
- `dim im(partial) = 13`;
- point and line permutation modules split as `1 + chi_6`;
- flag permutation character splits as `1 + 2*chi_6 + St`;
- character inner products satisfy
  `<chi_6, chi_6> = 1`, `<St, St> = 1`, and `<chi_6, St> = 0`;
- the two 3-dimensional characters satisfy
  `<chi_3, chi_3> = 1`, `<chi_3, chi_3bar> = 0`,
  `Sym^2(chi_3) = chi_6`, and `chi_3*chi_3bar = 1 + St`.

For `GL_2(F_5)`:

- group order: `480`;
- projective-line size: `6`;
- `dim St = 5`;
- `C[P^1] = 1 + St`;
- character inner product `<St, St> = 1`;
- representation-family counts:
  `4` determinant characters, `4` Steinberg twists, `6` principal series,
  and `10` cuspidal representations;
- sum of squares of listed dimensions: `480`;
- each constructed principal series has character inner product `1`;
- `--cuspidal 1` builds a 96-dimensional Gelfand-Graev representation,
  projects to rank `4`, and verifies the resulting cuspidal character has
  inner product `1`.

For `GL_3(F_3)`:

- group order: `11232`;
- number of points: `13`;
- number of planes: `13`;
- number of incidence edges: `52`;
- boundary rank: `25`;
- dimension of `ker(partial)`: `27`;
- number of unique apartments up to orientation: `234`;
- rank of the apartment-cycle span: `27`;
- Weyl group size: `6`;
- number of Weyl-group orbits on apartments: `49`;
- number of conjugacy classes: `24`;
- inner product of the Steinberg character with itself: `1`.

The `GL_3(F_3)` script ran end-to-end in under 20 seconds on the machine used
to create this repository.

For general `q`, the script checks the formulas:

```text
# points = # planes = q^2 + q + 1
# edges = (q + 1)(q^2 + q + 1)
dim St = #edges - #vertices + 1 = q^3
# apartments up to orientation = q^3(q + 1)(q^2 + q + 1) / 6
```

## Mathematical context

For a finite group of Lie type, the Steinberg representation can be constructed
from the top reduced homology of the associated spherical building. Here the
groups are `GL_3(F_q)`, so the building is one-dimensional, and its chambers
are flags

```text
line <= plane <= F_q^3.
```

Apartments come from choices of ordered bases of `F_q^3`. Each basis determines
a hexagonal cycle in the incidence graph, and the alternating sum of its
oriented chambers gives an apartment cycle. The scripts check that the span of
these apartment cycles has the expected Steinberg dimension.

Scalar matrices act trivially on the building, so this geometric action factors
through `PGL_3(F_q)`.

## Notes

This is exploratory research code rather than a packaged library. The
calculation is intentionally explicit so the construction of the building,
boundary map, group action, apartment cycles, and character can all be inspected
directly.

The repository is distributed under the MIT license; see the repository-level
`LICENSE` file.
