from argparse import ArgumentParser

from sage.all import *

from projective_line import build_projective_line_data, matrix_of
from representations_gl_2_q import (
    additive_character_data,
    additive_character_value,
    build_cuspidal_representation,
    cuspidal_character_rows,
    cuspidal_matrix,
    cuspidal_parameter_orbits,
    hermitian_inner_product,
    nonsplit_torus_character_value,
)


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct a cuspidal GL_2(F_q) representation in a Whittaker "
            "basis indexed by F_q^*. The script first builds the existing "
            "Gelfand-Graev/projector model, then changes basis."
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
        "--cuspidal",
        type=int,
        metavar="M",
        help=(
            "regular exponent M for a character theta_M of F_{q^2}^*; "
            "defaults to the first cuspidal parameter"
        ),
    )
    parser.add_argument(
        "--basis",
        action="store_true",
        help="print Whittaker basis vectors in the original cuspidal image basis",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print cuspidal character values computed in the Whittaker basis",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's GL(2,q) generators in the Whittaker basis",
    )
    return parser.parse_args()


def unipotent_matrix(F, x):
    return matrix(F, [[1, x], [0, 1]])


def scalar_matrix(F, z):
    return matrix(F, [[z, 0], [0, z]])


def nonzero_field_elements(F):
    return [x for x in F if x != 0]


def normalized_vector(v):
    for entry in v:
        if entry != 0:
            return v / entry
    raise ValueError("cannot normalize the zero vector")


def unipotent_character_projector(y, cuspidal):
    """Project to the U-character n(x) |-> psi(y*x) inside a cuspidal module."""
    data = cuspidal["data"]
    F = data["F"]
    q = data["q"]
    K = cuspidal["base_ring"]
    additive_character = additive_character_data(F, K)
    dimension = cuspidal["image"].dimension()
    projector = Matrix(K, dimension, dimension, 0)

    for x in F:
        coefficient = additive_character_value(-y * x, additive_character)
        projector += coefficient * cuspidal_matrix(unipotent_matrix(F, x), cuspidal)

    return projector / K(q)


def build_whittaker_model(cuspidal):
    """Build a Whittaker basis indexed by F_q^* from the projector model."""
    data = cuspidal["data"]
    F = data["F"]
    K = cuspidal["base_ring"]
    labels = nonzero_field_elements(F)
    vectors = []
    projectors = {}

    for y in labels:
        projector = unipotent_character_projector(y, cuspidal)
        image = projector.column_space()

        assert projector * projector == projector
        assert image.dimension() == 1

        vector = normalized_vector(image.basis()[0])
        vectors.append(vector)
        projectors[y] = projector

    basis_matrix = Matrix(K, vectors).transpose()
    assert basis_matrix.rank() == len(labels)

    return {
        "data": data,
        "cuspidal": cuspidal,
        "base_ring": K,
        "labels": labels,
        "basis_matrix": basis_matrix,
        "basis_inverse": basis_matrix.inverse(),
        "projectors": projectors,
    }


def whittaker_matrix(M, whittaker):
    old_matrix = cuspidal_matrix(M, whittaker["cuspidal"])
    return whittaker["basis_inverse"] * old_matrix * whittaker["basis_matrix"]


def whittaker_character_rows(data, whittaker):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "trace": whittaker_matrix(
                    representative.matrix(),
                    whittaker,
                ).trace(),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], str(row["trace"])))
    return rows


def diagonal_matrix_from_entries(K, entries):
    return diagonal_matrix(K, list(entries))


def verify_unipotent_action(whittaker):
    data = whittaker["data"]
    F = data["F"]
    K = whittaker["base_ring"]
    additive_character = additive_character_data(F, K)

    for x in F:
        actual = whittaker_matrix(unipotent_matrix(F, x), whittaker)
        expected = diagonal_matrix_from_entries(
            K,
            [
                additive_character_value(y * x, additive_character)
                for y in whittaker["labels"]
            ],
        )
        assert actual == expected


def verify_central_character(whittaker):
    data = whittaker["data"]
    F = data["F"]
    K = whittaker["base_ring"]
    torus_character = whittaker["cuspidal"]["torus_character"]
    identity = identity_matrix(K, len(whittaker["labels"]))

    for z in nonzero_field_elements(F):
        actual = whittaker_matrix(scalar_matrix(F, z), whittaker)
        expected = nonsplit_torus_character_value(z, torus_character) * identity
        assert actual == expected


def verify_representation_action(whittaker):
    G = whittaker["data"]["G"]
    g = G.random_element().matrix()
    h = G.random_element().matrix()

    assert whittaker_matrix(g * h, whittaker) == (
        whittaker_matrix(g, whittaker) * whittaker_matrix(h, whittaker)
    )


def verify_character_rows(data, whittaker):
    formula_rows = cuspidal_character_rows(
        data,
        whittaker["cuspidal"]["torus_character"],
    )
    whittaker_rows = whittaker_character_rows(data, whittaker)

    assert hermitian_inner_product(formula_rows, "trace", "trace") == 1

    for formula_row, whittaker_row in zip(formula_rows, whittaker_rows):
        assert formula_row["order"] == whittaker_row["order"]
        assert formula_row["size"] == whittaker_row["size"]
        assert formula_row["trace"] == whittaker_row["trace"]

    return whittaker_rows


def verify_whittaker_model(data, whittaker):
    verify_unipotent_action(whittaker)
    verify_central_character(whittaker)
    verify_representation_action(whittaker)
    return verify_character_rows(data, whittaker)


def print_basis_vectors(whittaker):
    print("Whittaker basis vectors in the original cuspidal image basis:")
    for y, vector in zip(whittaker["labels"], whittaker["basis_matrix"].columns()):
        print(f"  y = {y}: {vector}")


def print_character_table(rows):
    print("Cuspidal character by conjugacy class, computed in the Whittaker basis:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"trace = {row['trace']}"
        )


def print_generator_matrices(whittaker):
    data = whittaker["data"]

    print("Matrices for Sage's GL(2,q) generators in the Whittaker basis:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = matrix_of(generator)
        print(f"generator {i} in GL_2(F_{data['q']}):")
        print(M)
        print(whittaker_matrix(M, whittaker))


def default_cuspidal_exponent(q):
    return cuspidal_parameter_orbits(q)[0][0]


def main():
    args = parse_args()
    data = build_projective_line_data(args.q)
    q = data["q"]
    exponent = args.cuspidal if args.cuspidal is not None else default_cuspidal_exponent(q)

    cuspidal = build_cuspidal_representation(data, exponent)
    whittaker = build_whittaker_model(cuspidal)
    rows = verify_whittaker_model(data, whittaker)

    print(f"Constructing Whittaker basis for a cuspidal representation of GL_2(F_{q})")
    print(f"cuspidal exponent: {cuspidal['exponent']}")
    print(f"Whittaker labels: {[str(y) for y in whittaker['labels']]}")
    print(f"dimension: {len(whittaker['labels'])}")
    print("verified: U acts diagonally by n(x) |-> psi(y*x)")
    print("verified: central scalars act by theta restricted to F_q^*")
    print("verified: traces match the nonsplit-torus cuspidal character formula")
    print(
        f"<pi_{cuspidal['exponent']}, pi_{cuspidal['exponent']}>: "
        f"{hermitian_inner_product(rows, 'trace', 'trace')}"
    )

    if args.basis:
        print_basis_vectors(whittaker)

    if args.character_table:
        print_character_table(rows)

    if args.generators:
        print_generator_matrices(whittaker)


if __name__ == "__main__":
    main()
