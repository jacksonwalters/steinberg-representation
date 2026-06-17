from itertools import permutations

from sage.all import *


q = 3
VERBOSE = False

F = GF(q)
V = VectorSpace(F, 3)
G = GL(3, F)


def canonical_cycle_tuple(cycle):
    cycle_tuple = tuple(cycle)
    neg_cycle_tuple = tuple(-x for x in cycle_tuple)
    return min(cycle_tuple, neg_cycle_tuple)


def normalized_vector_tuple(v):
    coords = list(v)
    for x in coords:
        if x != 0:
            inv = x ** (-1)
            return tuple(int(inv * y) for y in coords)
    raise ValueError("zero vector has no projective normalization")


def line_repr(L):
    return normalized_vector_tuple(L.basis()[0])


def linear_form_str(a):
    terms = []
    for i, ai in enumerate(a):
        if ai == 1:
            terms.append(f"x{i + 1}")
        elif ai != 0:
            terms.append(f"{ai}*x{i + 1}")
    return " + ".join(terms) if terms else "0"


def plane_normal(H):
    for a in V:
        if a != 0 and all(a * v == 0 for v in H.basis()):
            return normalized_vector_tuple(a)
    raise ValueError("could not find a normal vector")


def plane_repr(H):
    return "{" + linear_form_str(plane_normal(H)) + " = 0}"


def edge_repr(edge):
    p, H = edge
    return f"{line_repr(p)} <= {plane_repr(H)}"


def pretty_print_cycle(cycle, edges):
    for j, coeff in enumerate(cycle):
        if coeff != 0:
            print(f"{coeff:+} * [{edge_repr(edges[j])}]")


print(f"Computing the Steinberg representation of GL_3(F_{q})")
print(f"|GL(3,{q})|:", G.order())

# The vertices of the rank-2 building are the proper nonzero subspaces
# of F_q^3, hence the projective points and projective lines.
points = list(V.subspaces(1))
planes = list(V.subspaces(2))

print("number of points:", len(points))
print("number of planes:", len(planes))

edges = [(p, H) for p in points for H in planes if p.is_subspace(H)]
print("number of incidence edges:", len(edges))

point_index = {points[i]: i for i in range(len(points))}
plane_index = {planes[i]: i for i in range(len(planes))}
edge_index = {edges[i]: i for i in range(len(edges))}

num_vertices = len(points) + len(planes)
num_edges = len(edges)

boundary = Matrix(QQ, num_vertices, num_edges)

for (p, H), j in edge_index.items():
    boundary[point_index[p], j] = -1
    boundary[plane_index[H] + len(points), j] = 1

ker = boundary.right_kernel()
steinberg_basis = ker.basis()

print("boundary rank:", boundary.rank())
print("dimension of ker(partial):", ker.dimension())
print("expected Steinberg dimension q^3:", q**3)
assert ker.dimension() == q**3

if VERBOSE:
    print("first Steinberg basis vector:")
    print(steinberg_basis[0])


def apartment_cycle_from_columns(columns):
    lines = [V.subspace([V(v)]) for v in columns]
    cycle = vector(QQ, len(edges))

    for w in permutations([0, 1, 2]):
        sign = Permutation([i + 1 for i in w]).signature()
        p = lines[w[0]]
        H = lines[w[0]] + lines[w[1]]
        cycle[edge_index[(p, H)]] += sign

    return cycle


def apartment_cycle_from_g(g):
    M = g.matrix()
    return apartment_cycle_from_columns([M.column(i) for i in range(3)])


v_A = apartment_cycle_from_g(G.one())
print("standard apartment cycle is in ker(partial):", (boundary * v_A).is_zero())
print("standard apartment cycle support size:", len([x for x in v_A if x != 0]))
assert (boundary * v_A).is_zero()


def act_on_subspace(M, S):
    W = S.ambient_vector_space()
    basis = S.basis_matrix().transpose()
    vecs = [W(M * basis.column(i)) for i in range(basis.ncols())]
    return W.subspace(vecs)


def edge_permutation(M):
    return [
        edge_index[(act_on_subspace(M, p), act_on_subspace(M, H))]
        for p, H in edges
    ]


def action_on_vector(M, v):
    perm = edge_permutation(M)
    v_new = vector(QQ, len(v))

    for i, val in enumerate(v):
        if val != 0:
            v_new[perm[i]] += val

    return v_new


g = G.random_element().matrix()
h = G.random_element().matrix()
v0 = steinberg_basis[0]

lhs = action_on_vector(g, action_on_vector(h, v0))
rhs = action_on_vector(g * h, v0)
assert lhs == rhs, "G-action property failed"
print("G-action property verified on a random pair")


unique_apartments = set()

for g in G:
    unique_apartments.add(canonical_cycle_tuple(apartment_cycle_from_g(g)))

expected_apartments = G.order() // ((q - 1) ** 3 * factorial(3))

print("number of unique apartments up to orientation:", len(unique_apartments))
print("expected number of apartments:", expected_apartments)
assert len(unique_apartments) == expected_apartments

apartment_span_rank = matrix(QQ, [list(v) for v in unique_apartments]).rank()
print("rank of apartment-cycle span:", apartment_span_rank)
assert apartment_span_rank == ker.dimension()

Weyl_group = []

for perm in permutations([0, 1, 2]):
    M = matrix(F, 3, 3, 0)
    for i, p_i in enumerate(perm):
        M[p_i, i] = 1
    Weyl_group.append(M)

print("Weyl group size:", len(Weyl_group))


def apartment_orbit(apartment):
    orbit_set = set()
    v = vector(QQ, apartment)

    for w in Weyl_group:
        orbit_set.add(canonical_cycle_tuple(action_on_vector(w, v)))

    return orbit_set


orbits = []
visited = set()

for apartment in unique_apartments:
    if apartment in visited:
        continue

    orbit = apartment_orbit(apartment)
    orbits.append(orbit)
    visited.update(orbit)

orbit_sizes = sorted(len(orbit) for orbit in orbits)
print("number of Weyl-group orbits on apartments:", len(orbits))
if VERBOSE:
    print("Weyl-group orbit sizes:", orbit_sizes)


def fixed_counts(M):
    point_images = {p: act_on_subspace(M, p) for p in points}
    plane_images = {H: act_on_subspace(M, H) for H in planes}

    fixed_vertices = sum(1 for p in points if point_images[p] == p)
    fixed_vertices += sum(1 for H in planes if plane_images[H] == H)

    fixed_edges = sum(
        1
        for p, H in edges
        if point_images[p] == p and plane_images[H] == H
    )

    return fixed_vertices, fixed_edges


def steinberg_trace(M):
    fixed_vertices, fixed_edges = fixed_counts(M)
    # For a connected graph, H_1 - H_0 = C_1 - C_0 and H_0 is trivial.
    return fixed_edges - fixed_vertices + 1


print("Computing Steinberg character from fixed simplex counts...")

rows = []

for C in G.conjugacy_classes():
    rep = C.representative()
    fixed_vertices, fixed_edges = fixed_counts(rep.matrix())
    rows.append(
        {
            "order": rep.order(),
            "size": len(C),
            "fixed_vertices": fixed_vertices,
            "fixed_edges": fixed_edges,
            "trace": steinberg_trace(rep.matrix()),
        }
    )

rows.sort(key=lambda row: (row["order"], row["size"], row["trace"]))

if VERBOSE:
    for row in rows:
        print(
            f"order {row['order']:2}, size {row['size']:4} : "
            f"fixed vertices = {row['fixed_vertices']:2}, "
            f"fixed edges = {row['fixed_edges']:2}, "
            f"chi = {row['trace']:3}"
        )

inner = sum(ZZ(row["size"]) * ZZ(row["trace"]) ** 2 for row in rows) / ZZ(G.order())
print("Number of conjugacy classes:", len(rows))
print("Inner product of Steinberg character with itself:", inner)
assert inner == 1
