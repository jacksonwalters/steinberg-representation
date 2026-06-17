from argparse import ArgumentParser

from sage.all import *


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct Fano-plane representations of GL_3(F_2), "
            "equivalently the full collineation group of the Fano plane."
        )
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print character values on the conjugacy classes of GL_3(F_2)",
    )
    parser.add_argument(
        "--incidence",
        action="store_true",
        help="print the seven projective lines as triples of Fano-plane points",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's group generators in the main modules",
    )
    return parser.parse_args()


def normalized_vector_tuple(v):
    coords = list(v)
    for x in coords:
        if x != 0:
            inv = x**(-1)
            return tuple(int(inv * y) for y in coords)
    raise ValueError("zero vector has no projective normalization")


def point_repr(point):
    return normalized_vector_tuple(point.basis()[0])


def line_repr(line, points):
    return tuple(
        sorted(point_repr(point) for point in points if point.is_subspace(line))
    )


def build_fano_data():
    F = GF(2)
    V = VectorSpace(F, 3)
    G = GL(3, F)

    points = list(V.subspaces(1))
    lines = list(V.subspaces(2))
    flags = [(point, line) for point in points for line in lines if point.is_subspace(line)]

    point_index = {point: i for i, point in enumerate(points)}
    line_index = {line: i for i, line in enumerate(lines)}
    flag_index = {flag: i for i, flag in enumerate(flags)}

    # This is the building boundary C_1 -> C_0, with each flag oriented
    # from its point vertex to its line vertex.
    boundary = Matrix(QQ, len(points) + len(lines), len(flags))
    for (point, line), j in flag_index.items():
        boundary[point_index[point], j] = -1
        boundary[len(points) + line_index[line], j] = 1

    return {
        "F": F,
        "V": V,
        "G": G,
        "points": points,
        "lines": lines,
        "flags": flags,
        "point_index": point_index,
        "line_index": line_index,
        "flag_index": flag_index,
        "boundary": boundary,
        "steinberg": boundary.right_kernel(),
        "boundary_image": boundary.column_space(),
    }


def matrix_of(g):
    return g.matrix() if hasattr(g, "matrix") else g


def act_on_subspace(M, S):
    W = S.ambient_vector_space()
    basis = S.basis_matrix().transpose()
    vecs = [W(M * basis.column(i)) for i in range(basis.ncols())]
    return W.subspace(vecs)


def permutation_matrix_from_images(images):
    n = len(images)
    P = Matrix(QQ, n, n, 0)
    for i, j in enumerate(images):
        P[j, i] = 1
    return P


def point_images(M, data):
    M = matrix_of(M)
    return [
        data["point_index"][act_on_subspace(M, point)]
        for point in data["points"]
    ]


def line_images(M, data):
    M = matrix_of(M)
    return [
        data["line_index"][act_on_subspace(M, line)]
        for line in data["lines"]
    ]


def flag_images(M, data):
    M = matrix_of(M)
    return [
        data["flag_index"][
            (
                act_on_subspace(M, point),
                act_on_subspace(M, line),
            )
        ]
        for point, line in data["flags"]
    ]


def vertex_images(M, data):
    point_part = point_images(M, data)
    line_part = [len(data["points"]) + j for j in line_images(M, data)]
    return point_part + line_part


def point_permutation_matrix(M, data):
    return permutation_matrix_from_images(point_images(M, data))


def line_permutation_matrix(M, data):
    return permutation_matrix_from_images(line_images(M, data))


def flag_permutation_matrix(M, data):
    return permutation_matrix_from_images(flag_images(M, data))


def vertex_permutation_matrix(M, data):
    return permutation_matrix_from_images(vertex_images(M, data))


def augmentation_matrix(permutation_matrix):
    n = permutation_matrix.nrows()
    columns = []

    for i in range(n - 1):
        v = vector(QQ, n)
        v[i] = 1
        v[n - 1] = -1

        image = permutation_matrix * v
        assert sum(image) == 0
        columns.append(list(image[: n - 1]))

    return matrix(QQ, columns).transpose()


def restricted_matrix(linear_map, invariant_subspace):
    columns = []

    for basis_vector in invariant_subspace.basis():
        image = linear_map * basis_vector
        columns.append(invariant_subspace.coordinates(image))

    return matrix(QQ, columns).transpose()


def representation_matrices(M, data):
    point_perm = point_permutation_matrix(M, data)
    line_perm = line_permutation_matrix(M, data)
    flag_perm = flag_permutation_matrix(M, data)
    vertex_perm = vertex_permutation_matrix(M, data)

    return {
        "points": point_perm,
        "lines": line_perm,
        "flags": flag_perm,
        "vertices": vertex_perm,
        "point_augmentation": augmentation_matrix(point_perm),
        "line_augmentation": augmentation_matrix(line_perm),
        "boundary_image": restricted_matrix(vertex_perm, data["boundary_image"]),
        "steinberg": restricted_matrix(flag_perm, data["steinberg"]),
    }


def fixed_count(images):
    return sum(1 for i, j in enumerate(images) if i == j)


def character_rows(data):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        M = representative.matrix()

        fixed_points = fixed_count(point_images(M, data))
        fixed_lines = fixed_count(line_images(M, data))
        fixed_flags = fixed_count(flag_images(M, data))
        fixed_vertices = fixed_points + fixed_lines

        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "fixed_points": fixed_points,
                "fixed_lines": fixed_lines,
                "fixed_flags": fixed_flags,
                "point_permutation": fixed_points,
                "line_permutation": fixed_lines,
                "point_augmentation": fixed_points - 1,
                "line_augmentation": fixed_lines - 1,
                "vertices": fixed_vertices,
                "boundary_image": fixed_vertices - 1,
                "flags": fixed_flags,
                "steinberg": fixed_flags - fixed_vertices + 1,
            }
        )

    rows.sort(
        key=lambda row: (
            row["order"],
            row["size"],
            row["fixed_points"],
            row["fixed_flags"],
        )
    )
    return rows


def inner_product(rows, left_key, right_key):
    group_order = sum(ZZ(row["size"]) for row in rows)
    return (
        sum(
            ZZ(row["size"]) * ZZ(row[left_key]) * ZZ(row[right_key])
            for row in rows
        )
        / group_order
    )


def verify_action_property(data, representation):
    G = data["G"]
    g = G.random_element().matrix()
    h = G.random_element().matrix()
    return representation(g * h, data) == representation(g, data) * representation(h, data)


def verify_representations(data, rows):
    G = data["G"]
    boundary = data["boundary"]
    steinberg = data["steinberg"]
    boundary_image = data["boundary_image"]

    assert G.order() == 168
    assert len(data["points"]) == 7
    assert len(data["lines"]) == 7
    assert len(data["flags"]) == 21
    assert boundary.rank() == 13
    assert steinberg.dimension() == 8
    assert boundary_image.dimension() == 13

    assert verify_action_property(data, point_permutation_matrix)
    assert verify_action_property(data, line_permutation_matrix)
    assert verify_action_property(data, flag_permutation_matrix)

    g = G.random_element().matrix()
    assert boundary * flag_permutation_matrix(g, data) == vertex_permutation_matrix(g, data) * boundary

    for conjugacy_class in G.conjugacy_classes():
        M = conjugacy_class.representative().matrix()
        matrices = representation_matrices(M, data)

        row = next(
            row
            for row in rows
            if row["order"] == conjugacy_class.representative().order()
            and row["size"] == len(conjugacy_class)
            and row["fixed_flags"] == fixed_count(flag_images(M, data))
        )

        assert matrices["point_augmentation"].trace() == row["point_augmentation"]
        assert matrices["line_augmentation"].trace() == row["line_augmentation"]
        assert matrices["boundary_image"].trace() == row["boundary_image"]
        assert matrices["steinberg"].trace() == row["steinberg"]
        assert row["flags"] == row["boundary_image"] + row["steinberg"]
        assert row["flags"] == 1 + 2 * row["point_augmentation"] + row["steinberg"]

    assert inner_product(rows, "point_augmentation", "point_augmentation") == 1
    assert inner_product(rows, "line_augmentation", "line_augmentation") == 1
    assert inner_product(rows, "steinberg", "steinberg") == 1
    assert inner_product(rows, "point_augmentation", "steinberg") == 0
    assert inner_product(rows, "flags", "flags") == 6


def print_incidence(data):
    print("Fano-plane lines as triples of points:")
    for i, line in enumerate(data["lines"], start=1):
        print(f"  L{i}: {line_repr(line, data['points'])}")


def print_character_table(rows):
    print("Characters by conjugacy class:")
    print(
        "  "
        "ord size | fix_pts fix_lines fix_flags | "
        "chi_6 chi_St chi_flags"
    )
    for row in rows:
        print(
            f"  {row['order']:3} {row['size']:4} | "
            f"{row['fixed_points']:7} {row['fixed_lines']:9} {row['fixed_flags']:9} | "
            f"{row['point_augmentation']:5} {row['steinberg']:6} {row['flags']:9}"
        )


def print_generator_matrices(data):
    print("Matrices for Sage's GL(3,2) generators:")

    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        matrices = representation_matrices(M, data)
        print(f"generator {i} in GL_3(F_2):")
        print(M)
        print("  point permutation module C[points]:")
        print(matrices["points"])
        print("  6-dimensional point augmentation module:")
        print(matrices["point_augmentation"])
        print("  8-dimensional Steinberg module:")
        print(matrices["steinberg"])


def main():
    args = parse_args()
    data = build_fano_data()
    rows = character_rows(data)
    verify_representations(data, rows)

    boundary = data["boundary"]

    print("Constructing Fano-plane representations of GL_3(F_2)")
    print(f"group order: {data['G'].order()}")
    print(
        "Fano plane: "
        f"points={len(data['points'])}, "
        f"lines={len(data['lines'])}, "
        f"flags={len(data['flags'])}"
    )
    print(f"building boundary rank: {boundary.rank()}")
    print(f"dim ker(partial) = dim Steinberg module: {data['steinberg'].dimension()}")
    print(f"dim im(partial) in C_0: {data['boundary_image'].dimension()}")
    print("point permutation module: dim 7 = 1 + irreducible dim 6")
    print("line permutation module:  dim 7 = 1 + irreducible dim 6")
    print("flag permutation module:  dim 21 = im(partial) dim 13 + Steinberg dim 8")
    print("character decomposition check: C[flags] = 1 + 2*chi_6 + St")
    print(
        "inner products: "
        f"<chi_6, chi_6>={inner_product(rows, 'point_augmentation', 'point_augmentation')}, "
        f"<St, St>={inner_product(rows, 'steinberg', 'steinberg')}, "
        f"<chi_6, St>={inner_product(rows, 'point_augmentation', 'steinberg')}"
    )

    if args.incidence:
        print_incidence(data)

    if args.character_table:
        print_character_table(rows)

    if args.generators:
        print_generator_matrices(data)


if __name__ == "__main__":
    main()
