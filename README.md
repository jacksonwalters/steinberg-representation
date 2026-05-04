# Steinberg Representation of `GL_3(F_2)`

This folder contains a small SageMath computation of the Steinberg
representation of the finite general linear group `GL_3(F_2)`.

This corresponds to the expository paper: [https://jacksonwalters.com/docs/notes/steinberg_representation_GL(3,2).pdf](https://jacksonwalters.com/docs/notes/steinberg_representation_GL(3,2).pdf)

The representation is realized geometrically from the spherical building of
`F_2^3`. In this rank-2 case, the building is the incidence graph whose
vertices are:

- the 1-dimensional subspaces of `F_2^3`;
- the 2-dimensional subspaces of `F_2^3`;
- edges given by inclusions `line <= plane`.

The Steinberg representation is computed as the top homology of this building:

```text
St = H_1(Delta; C) = ker(partial: C_1(Delta; C) -> C_0(Delta; C)).
```

For `GL_3(F_2)`, the computation finds an 8-dimensional Steinberg module.

## Files

- `steinberg_representation.py` is the runnable SageMath script.
- `steinberg_representation.ipynb` is the notebook version of the same
  computation.

The script is the best entry point if you want to reproduce the calculation
from a terminal. The notebook is useful for inspecting the computation cell by
cell.

## What the script does

The script:

1. Constructs the 7 points and 7 planes in `F_2^3`.
2. Builds the 21 incidence edges of the spherical building.
3. Forms the boundary map `C_1 -> C_0`.
4. Computes `ker(partial)`, giving the Steinberg representation.
5. Constructs apartment cycles from bases of `F_2^3`.
6. Defines the natural `GL_3(F_2)` action on the building and on cycles.
7. Verifies the group action relation `g(h(v)) = (gh)(v)`.
8. Counts apartment cycles up to orientation.
9. Computes the span of a `GL_3(F_2)` orbit of an apartment cycle.
10. Computes the character of the resulting representation.

## Running the computation

Install SageMath (tested with SageMath 10.8), then run from this directory:

```bash
sage steinberg_representation.py
```

Or run it from the repository root:

```bash
sage representation_theory/finite_general_linear_group/steinberg_representation/steinberg_representation.py
```

To open the notebook with a SageMath kernel:

```bash
sage -n jupyter steinberg_representation.ipynb
```

## Expected output

Some printed vectors depend on the random group element selected during the run,
but the structural checks should be stable:

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

The character values by conjugacy-class order and class size are:

```text
order  1, size   1 : chi = 8
order  2, size  21 : chi = 0
order  3, size  56 : chi = 1
order  4, size  42 : chi = 0
order  7, size  24 : chi = 1
order  7, size  24 : chi = 1
```

The final character inner product `1.0` verifies irreducibility in this
example.

## Mathematical context

For a finite group of Lie type, the Steinberg representation can be constructed
from the top reduced homology of the associated spherical building. Here the
group is `GL_3(F_2)`, so the building is one-dimensional, and its chambers are
flags

```text
line <= plane <= F_2^3.
```

Apartments come from choices of ordered bases of `F_2^3`. Each basis determines
a hexagonal cycle in the incidence graph, and the alternating sum of its
oriented chambers gives an apartment cycle. The script checks that the
`GL_3(F_2)`-span of such a cycle has dimension 8, matching the computed
Steinberg module.

## Notes

This is exploratory research code rather than a packaged library. The
calculation is intentionally explicit so the construction of the building,
boundary map, group action, apartment cycles, and character can all be inspected
directly.

The repository is distributed under the MIT license; see the repository-level
`LICENSE` file.
