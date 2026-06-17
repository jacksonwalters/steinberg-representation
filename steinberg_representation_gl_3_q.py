from argparse import ArgumentParser
from itertools import combinations, permutations

from sage.all import *


APARTMENT_SPAN_AUTO_LIMIT = 5000
CHARACTER_AUTO_Q_LIMIT = 5


def parse_args():
    parser = ArgumentParser(
        description="Construct the Steinberg representation of GL_3(F_q)."
    )
    parser.add_argument(
        "q",
        nargs="?",
        type=int,
        default=3,
        help="prime power q for the finite field F_q",
    )
    parser.add_argument(
        "--apartment-span",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute the rank of the span of apartment cycles",
    )
    parser.add_argument(
        "--character",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute the Steinberg character by fixed simplex counts",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print one row for each conjugacy class when computing the character",
    )
    parser.add_argument(
        "--weyl-orbits",
        action="store_true",
        help="compute Weyl-group orbits on apartments when apartments are generated",
    )
    return parser.parse_args()


def canonical_cycle_tuple(cycle):
    cycle_tuple = tuple(cycle)
    neg_cycle_tuple = tuple(-x for x in cycle_tuple)
    return min(cycle_tuple, neg_cycle_tuple)


def build_data(q):
    if q < 2:
        raise ValueError("q must be a prime power greater than 1")

    F = GF(q, name="a")
    V = VectorSpace(F, 3)
    G = GL(3, F)

    points = list(V.subspaces(1))
    planes = list(V.subspaces(2))
    edges = [(p, H) for p in points for H in planes if p.is_subspace(H)]

    point_index = {p: i for i, p in enumerate(points)}
    plane_index = {H: i for i, H in enumerate(planes)}
    edge_index = {edge: i for i, edge in enumerate(edges)}

    boundary = Matrix(QQ, len(points) + len(planes), len(edges))
    for (p, H), j in edge_index.items():
        boundary[point_index[p], j] = -1
        boundary[plane_index[H] + len(points), j] = 1

    return {
        "q": q,
        "F": F,
        "V": V,
        "G": G,
        "points": points,
        "planes": planes,
        "edges": edges,
        "edge_index": edge_index,
        "boundary": boundary,
        "kernel": boundary.right_kernel(),
    }


def expected_apartment_count(q):
    return q**3 * (q + 1) * (q**2 + q + 1) // factorial(3)


def apartment_cycle_from_lines(lines, edge_index, num_edges):
    cycle = vector(QQ, num_edges)

    for w in permutations([0, 1, 2]):
        sign = Permutation([i + 1 for i in w]).signature()
        p = lines[w[0]]
        H = lines[w[0]] + lines[w[1]]
        cycle[edge_index[(p, H)]] += sign

    return cycle


def standard_apartment_cycle(data):
    V = data["V"]
    lines = [V.subspace([v]) for v in V.basis()]
    return apartment_cycle_from_lines(lines, data["edge_index"], len(data["edges"]))


def noncollinear_line_triples(points):
    for lines in combinations(points, 3):
        if (lines[0] + lines[1] + lines[2]).dimension() == 3:
            yield lines


def generate_apartment_cycles(data):
    return {
        canonical_cycle_tuple(
            apartment_cycle_from_lines(
                lines,
                data["edge_index"],
                len(data["edges"]),
            )
        )
        for lines in noncollinear_line_triples(data["points"])
    }


def act_on_subspace(M, S):
    W = S.ambient_vector_space()
    basis = S.basis_matrix().transpose()
    vecs = [W(M * basis.column(i)) for i in range(basis.ncols())]
    return W.subspace(vecs)


def edge_permutation(M, data):
    return [
        data["edge_index"][(act_on_subspace(M, p), act_on_subspace(M, H))]
        for p, H in data["edges"]
    ]


def action_on_vector(M, v, data):
    perm = edge_permutation(M, data)
    v_new = vector(QQ, len(v))

    for i, val in enumerate(v):
        if val != 0:
            v_new[perm[i]] += val

    return v_new


def check_random_action(data):
    G = data["G"]
    basis = data["kernel"].basis()
    g = G.random_element().matrix()
    h = G.random_element().matrix()
    v = basis[0]

    lhs = action_on_vector(g, action_on_vector(h, v, data), data)
    rhs = action_on_vector(g * h, v, data)
    return lhs == rhs


def weyl_group(data):
    F = data["F"]
    W = []

    for perm in permutations([0, 1, 2]):
        M = matrix(F, 3, 3, 0)
        for i, p_i in enumerate(perm):
            M[p_i, i] = 1
        W.append(M)

    return W


def apartment_orbit(apartment, W, data):
    v = vector(QQ, apartment)
    return {canonical_cycle_tuple(action_on_vector(w, v, data)) for w in W}


def weyl_orbit_count(apartments, data):
    W = weyl_group(data)
    orbits = []
    visited = set()

    for apartment in apartments:
        if apartment in visited:
            continue

        orbit = apartment_orbit(apartment, W, data)
        orbits.append(orbit)
        visited.update(orbit)

    return len(W), len(orbits), sorted(len(orbit) for orbit in orbits)


def fixed_counts(M, data):
    points = data["points"]
    planes = data["planes"]
    edges = data["edges"]

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


def steinberg_trace(M, data):
    fixed_vertices, fixed_edges = fixed_counts(M, data)
    return fixed_edges - fixed_vertices + 1


def character_rows(data):
    rows = []

    for C in data["G"].conjugacy_classes():
        rep = C.representative()
        fixed_vertices, fixed_edges = fixed_counts(rep.matrix(), data)
        rows.append(
            {
                "order": rep.order(),
                "size": len(C),
                "fixed_vertices": fixed_vertices,
                "fixed_edges": fixed_edges,
                "trace": steinberg_trace(rep.matrix(), data),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], row["trace"]))
    return rows


def should_run_apartment_span(mode, expected_apartments):
    if mode == "yes":
        return True
    if mode == "no":
        return False
    return expected_apartments <= APARTMENT_SPAN_AUTO_LIMIT


def should_run_character(mode, q):
    if mode == "yes":
        return True
    if mode == "no":
        return False
    return q <= CHARACTER_AUTO_Q_LIMIT


def main():
    args = parse_args()
    data = build_data(args.q)

    q = data["q"]
    G = data["G"]
    boundary = data["boundary"]
    kernel = data["kernel"]
    expected_vertices_per_type = q**2 + q + 1
    expected_edges = (q + 1) * expected_vertices_per_type
    expected_apartments = expected_apartment_count(q)

    print(f"Computing St(GL_3(F_{q}))")
    print(f"group order: {G.order()}")
    print(
        "building: "
        f"points={len(data['points'])}, "
        f"planes={len(data['planes'])}, "
        f"edges={len(data['edges'])}"
    )
    print(
        "expected: "
        f"points=planes={expected_vertices_per_type}, "
        f"edges={expected_edges}"
    )

    assert len(data["points"]) == expected_vertices_per_type
    assert len(data["planes"]) == expected_vertices_per_type
    assert len(data["edges"]) == expected_edges

    print(f"boundary rank: {boundary.rank()}")
    print(f"dim ker(partial): {kernel.dimension()} (expected {q**3})")
    assert kernel.dimension() == q**3

    v_A = standard_apartment_cycle(data)
    print(f"standard apartment cycle in ker(partial): {(boundary * v_A).is_zero()}")
    assert (boundary * v_A).is_zero()
    print(f"standard apartment support size: {len([x for x in v_A if x != 0])}")

    action_ok = check_random_action(data)
    print(f"G-action check on a random pair: {action_ok}")
    assert action_ok

    print(f"expected apartments up to orientation: {expected_apartments}")

    apartments = None
    if should_run_apartment_span(args.apartment_span, expected_apartments):
        apartments = generate_apartment_cycles(data)
        print(f"generated apartments up to orientation: {len(apartments)}")
        assert len(apartments) == expected_apartments

        apartment_span_rank = matrix(QQ, [list(v) for v in apartments]).rank()
        print(f"apartment-cycle span rank: {apartment_span_rank}")
        assert apartment_span_rank == kernel.dimension()
    else:
        print(
            "apartment-cycle span rank: skipped "
            f"(use --apartment-span yes to compute {expected_apartments} apartments)"
        )

    if args.weyl_orbits:
        if apartments is None:
            apartments = generate_apartment_cycles(data)
            print(f"generated apartments up to orientation: {len(apartments)}")
            assert len(apartments) == expected_apartments

        weyl_size, orbit_count, orbit_sizes = weyl_orbit_count(apartments, data)
        print(f"Weyl group size: {weyl_size}")
        print(f"Weyl-group orbit count on apartments: {orbit_count}")
        if len(orbit_sizes) <= 20:
            print(f"Weyl-group orbit sizes: {orbit_sizes}")

    if args.character_table or should_run_character(args.character, q):
        rows = character_rows(data)

        if args.character_table:
            for row in rows:
                print(
                    f"order {row['order']:3}, size {row['size']:7}: "
                    f"fixed vertices = {row['fixed_vertices']:3}, "
                    f"fixed edges = {row['fixed_edges']:3}, "
                    f"chi = {row['trace']:4}"
                )

        inner = (
            sum(ZZ(row["size"]) * ZZ(row["trace"]) ** 2 for row in rows)
            / ZZ(G.order())
        )
        print(f"character conjugacy classes: {len(rows)}")
        print(f"character inner product <St, St>: {inner}")
        assert inner == 1
    else:
        print("character: skipped (use --character yes to compute)")


if __name__ == "__main__":
    main()
