# Steinberg Representations of `GL_3(F_2)` and `GL_3(F_3)`

This folder contains small SageMath computations of the Steinberg
representations of the finite general linear groups `GL_3(F_2)` and
`GL_3(F_3)`.

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
`GL_3(F_3)`, the same construction gives a 27-dimensional Steinberg module.

## Files

- `steinberg_representation_gl_3_2.py` is the original `GL_3(F_2)` SageMath
  script.
- `steinberg_representation_gl_3_3.py` is the `GL_3(F_3)` SageMath script.

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

The `GL_3(F_3)` script computes the Steinberg character from fixed vertex and
edge counts instead of building a `27 x 27` matrix for every group element. For
a connected graph,

```text
chi_St(g) = # fixed edges - # fixed vertices + 1.
```

## Running the computation

Install SageMath (tested with SageMath 10.8), then run from this directory:

```bash
sage steinberg_representation_gl_3_2.py
sage steinberg_representation_gl_3_3.py
```

## Expected output

Some printed vectors depend on the random group element selected during the
`GL_3(F_2)` run, but the structural checks should be stable.

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
- inner product of the Steinberg character with itself: `1.0`.

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
- inner product of the Steinberg character with itself: `1`.

The `GL_3(F_3)` script ran end-to-end in under 20 seconds on the machine used
to create this repository.

## Mathematical context

For a finite group of Lie type, the Steinberg representation can be constructed
from the top reduced homology of the associated spherical building. Here the
groups are `GL_3(F_q)` for `q = 2, 3`, so the building is one-dimensional, and
its chambers are flags

```text
line <= plane <= F_q^3.
```

Apartments come from choices of ordered bases of `F_q^3`. Each basis determines
a hexagonal cycle in the incidence graph, and the alternating sum of its
oriented chambers gives an apartment cycle. The scripts check that the span of
these apartment cycles has the expected Steinberg dimension.

## Notes

This is exploratory research code rather than a packaged library. The
calculation is intentionally explicit so the construction of the building,
boundary map, group action, apartment cycles, and character can all be inspected
directly.

The repository is distributed under the MIT license; see the repository-level
`LICENSE` file.
