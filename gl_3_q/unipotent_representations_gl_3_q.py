from argparse import ArgumentParser

from sage.all import *

from representation_utils_gl_3_q import (
    determinant_character_value,
    hermitian_inner_product,
    matrix_of,
    multiplicative_character_data,
    print_character_table,
)
from steinberg_representation_gl_3_q import (
    act_on_subspace,
    build_data,
    edge_permutation,
    steinberg_trace,
)


CHARACTER_AUTO_Q_LIMIT = 5


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct the unipotent families for GL_3(F_q): determinant "
            "characters, the projective-plane augmentation module, and "
            "Steinberg twists."
        )
    )
    parser.add_argument(
        "q",
        nargs="?",
        type=int,
        default=3,
        help="prime power q for the finite field F_q",
    )
    parser.add_argument(
        "--twist",
        type=int,
        default=0,
        metavar="A",
        help=(
            "determinant twist exponent: chi_A sends a multiplicative "
            "generator of F_q^* to zeta_(q-1)^A"
        ),
    )
    parser.add_argument(
        "--character",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute character inner products for the three unipotent twists",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print the three unipotent-twist characters by conjugacy class",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's GL(3,q) generators",
    )
    return parser.parse_args()


def should_run_character(mode, q):
    if mode == "yes":
        return True
    if mode == "no":
        return False
    return q <= CHARACTER_AUTO_Q_LIMIT


def projective_point_images(M, data):
    M = matrix_of(M)
    if "point_index" not in data:
        data["point_index"] = {point: i for i, point in enumerate(data["points"])}

    return [
        data["point_index"][act_on_subspace(M, point)]
        for point in data["points"]
    ]


def permutation_matrix_from_images(images, base_ring=QQ):
    n = len(images)
    P = Matrix(base_ring, n, n, 0)

    for i, j in enumerate(images):
        P[j, i] = base_ring(1)

    return P


def projective_point_permutation_matrix(M, data, base_ring=QQ):
    return permutation_matrix_from_images(
        projective_point_images(M, data),
        base_ring=base_ring,
    )


def augmentation_matrix_from_permutation(permutation_matrix):
    n = permutation_matrix.nrows()
    base_ring = permutation_matrix.base_ring()
    columns = []

    for i in range(n - 1):
        v = vector(base_ring, n)
        v[i] = 1
        v[n - 1] = -1

        image = permutation_matrix * v
        assert sum(image) == 0
        columns.append(list(image[: n - 1]))

    return matrix(base_ring, columns).transpose()


def middle_unipotent_matrix(M, data, base_ring=QQ):
    return augmentation_matrix_from_permutation(
        projective_point_permutation_matrix(M, data, base_ring=base_ring)
    )


def restricted_matrix(linear_map, invariant_subspace):
    columns = []

    for basis_vector in invariant_subspace.basis():
        image = linear_map * basis_vector
        columns.append(invariant_subspace.coordinates(image))

    return matrix(linear_map.base_ring(), columns).transpose()


def edge_permutation_matrix(M, data, base_ring=QQ):
    return permutation_matrix_from_images(
        edge_permutation(matrix_of(M), data),
        base_ring=base_ring,
    )


def steinberg_matrix(M, data):
    return restricted_matrix(edge_permutation_matrix(M, data), data["kernel"])


def determinant_twist_scalar(M, exponent, character_data):
    return determinant_character_value(M, exponent, character_data)


def one_dimensional_matrix(M, exponent, data, character_data=None):
    character_data = character_data or multiplicative_character_data(data["F"])
    K = character_data["base_ring"]
    return Matrix(K, 1, 1, [determinant_twist_scalar(M, exponent, character_data)])


def middle_unipotent_twist_matrix(M, exponent, data, character_data=None):
    character_data = character_data or multiplicative_character_data(data["F"])
    K = character_data["base_ring"]
    scalar = determinant_twist_scalar(M, exponent, character_data)
    return scalar * middle_unipotent_matrix(M, data, base_ring=K)


def steinberg_twist_matrix(M, exponent, data, character_data=None):
    character_data = character_data or multiplicative_character_data(data["F"])
    K = character_data["base_ring"]
    scalar = determinant_twist_scalar(M, exponent, character_data)
    return scalar * Matrix(K, steinberg_matrix(M, data))


def fixed_count(images):
    return sum(1 for i, j in enumerate(images) if i == j)


def character_rows(data, exponent, character_data=None):
    character_data = character_data or multiplicative_character_data(data["F"])
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        M = representative.matrix()
        scalar = determinant_twist_scalar(M, exponent, character_data)
        fixed_points = fixed_count(projective_point_images(M, data))
        st_trace = steinberg_trace(M, data)

        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "linear": scalar,
                "middle_unipotent": scalar * (fixed_points - 1),
                "steinberg": scalar * st_trace,
                "trace": scalar * (fixed_points - 1),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], str(row["middle_unipotent"])))
    return rows


def verify_unipotent_models(data, exponent, rows=None):
    G = data["G"]
    character_data = multiplicative_character_data(data["F"])
    g = G.random_element().matrix()
    h = G.random_element().matrix()

    assert one_dimensional_matrix(g * h, exponent, data, character_data) == (
        one_dimensional_matrix(g, exponent, data, character_data)
        * one_dimensional_matrix(h, exponent, data, character_data)
    )
    assert middle_unipotent_twist_matrix(g * h, exponent, data, character_data) == (
        middle_unipotent_twist_matrix(g, exponent, data, character_data)
        * middle_unipotent_twist_matrix(h, exponent, data, character_data)
    )
    assert steinberg_twist_matrix(g * h, exponent, data, character_data) == (
        steinberg_twist_matrix(g, exponent, data, character_data)
        * steinberg_twist_matrix(h, exponent, data, character_data)
    )

    if rows is not None:
        assert hermitian_inner_product(rows, "linear", "linear") == 1
        assert hermitian_inner_product(rows, "middle_unipotent", "middle_unipotent") == 1
        assert hermitian_inner_product(rows, "steinberg", "steinberg") == 1
        assert hermitian_inner_product(rows, "linear", "middle_unipotent") == 0
        assert hermitian_inner_product(rows, "linear", "steinberg") == 0
        assert hermitian_inner_product(rows, "middle_unipotent", "steinberg") == 0


def print_unipotent_character_table(rows):
    print("Unipotent-twist characters by conjugacy class:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"linear = {row['linear']}, "
            f"middle = {row['middle_unipotent']}, "
            f"St = {row['steinberg']}"
        )


def print_generator_matrices(data, exponent, character_data):
    print("Matrices for Sage's GL(3,q) generators in the unipotent models:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_3(F_{data['q']}):")
        print(M)
        print("  determinant character:")
        print(one_dimensional_matrix(M, exponent, data, character_data))
        print("  middle unipotent twist from C[P^2] / constants:")
        print(middle_unipotent_twist_matrix(M, exponent, data, character_data))
        print("  Steinberg twist:")
        print(steinberg_twist_matrix(M, exponent, data, character_data))


def main():
    args = parse_args()
    data = build_data(args.q)
    character_data = multiplicative_character_data(data["F"])
    exponent = args.twist % character_data["order"]

    rows = None
    if should_run_character(args.character, data["q"]) or args.character_table:
        rows = character_rows(data, exponent, character_data)

    verify_unipotent_models(data, exponent, rows)

    q = data["q"]
    print(f"Constructing unipotent representations of GL_3(F_{q})")
    print(f"determinant twist exponent: {exponent}")
    print(f"projective plane points: {len(data['points'])} = q^2 + q + 1")
    print(f"linear dimension: 1")
    print(f"middle unipotent dimension: {len(data['points']) - 1} = q(q + 1)")
    print(f"Steinberg dimension: {data['kernel'].dimension()} = q^3")

    if rows is not None:
        print(f"character conjugacy classes: {len(rows)}")
        print(f"<linear, linear>: {hermitian_inner_product(rows, 'linear', 'linear')}")
        print(
            "<middle, middle>: "
            f"{hermitian_inner_product(rows, 'middle_unipotent', 'middle_unipotent')}"
        )
        print(f"<St, St>: {hermitian_inner_product(rows, 'steinberg', 'steinberg')}")
    else:
        print("characters: skipped (use --character yes to compute)")

    if args.character_table:
        print_unipotent_character_table(rows)

    if args.generators:
        print_generator_matrices(data, exponent, character_data)


if __name__ == "__main__":
    main()
