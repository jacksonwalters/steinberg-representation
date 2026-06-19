from argparse import ArgumentParser

from sage.all import *

from projective_line import build_projective_line_data, matrix_of
from representations_gl_2_q import (
    additive_character_data,
    additive_character_value,
    cuspidal_base_ring,
    cuspidal_character_rows,
    cuspidal_character_value,
    cuspidal_parameter_orbits,
    hermitian_inner_product,
    nonsplit_torus_character_data,
    nonsplit_torus_character_value,
)


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct cuspidal GL_2(F_q) matrices directly in a Whittaker "
            "basis indexed by F_q^*, without first building the Gelfand-Graev "
            "projector model."
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
        "--character-table",
        action="store_true",
        help="print cuspidal character values computed in the direct model",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's GL(2,q) generators",
    )
    parser.add_argument(
        "--all-elements",
        action="store_true",
        help="verify the representation property on every pair of group elements",
    )
    return parser.parse_args()


def nonzero_field_elements(F):
    return [x for x in F if x != 0]


def unipotent_matrix(F, x):
    return matrix(F, [[1, x], [0, 1]])


def diagonal_matrix_gl2(F, a, d):
    return matrix(F, [[a, 0], [0, d]])


def weyl_matrix_gl2(F):
    return matrix(F, [[0, 1], [-1, 0]])


def direct_model_base_ring(data):
    return cuspidal_base_ring(data)


def default_cuspidal_exponent(q):
    return cuspidal_parameter_orbits(q)[0][0]


def label_index(labels):
    return {label: i for i, label in enumerate(labels)}


def build_direct_whittaker_model(data, exponent):
    K = direct_model_base_ring(data)
    torus_character = nonsplit_torus_character_data(data, exponent, K)
    labels = nonzero_field_elements(data["F"])

    model = {
        "data": data,
        "base_ring": K,
        "torus_character": torus_character,
        "exponent": torus_character["exponent"],
        "labels": labels,
        "label_index": label_index(labels),
        "additive_character": additive_character_data(data["F"], K),
    }
    model["weyl_matrix"] = compute_weyl_matrix(model)
    return model


def central_character_value(z, model):
    return nonsplit_torus_character_value(z, model["torus_character"])


def unipotent_action_matrix(x, model):
    K = model["base_ring"]
    additive_character = model["additive_character"]

    return diagonal_matrix(
        K,
        [
            additive_character_value(y * x, additive_character)
            for y in model["labels"]
        ],
    )


def diagonal_action_matrix(a, d, model):
    K = model["base_ring"]
    labels = model["labels"]
    index = model["label_index"]
    A = Matrix(K, len(labels), len(labels), 0)
    scalar = central_character_value(d, model)

    for j, y in enumerate(labels):
        A[index[y * d / a], j] = scalar

    return A


def upper_triangular_action_matrix(M, model):
    F = model["data"]["F"]
    M = matrix_of(M)
    a = F(M[0, 0])
    b = F(M[0, 1])
    d = F(M[1, 1])

    return (
        unipotent_action_matrix(b / d, model)
        * diagonal_action_matrix(a, d, model)
    )


def compute_weyl_matrix(model):
    """Compute the Weyl action from character values by finite Fourier inversion.

    If W is the matrix of w = [[0, 1], [-1, 0]], then

        tr(pi(diag(a, 1) w n(s))) = sum_y W[a*y, y] psi(y*s).

    Fourier inversion over the additive group of F_q therefore determines every
    matrix coefficient W[a*y, y].
    """
    data = model["data"]
    F = data["F"]
    q = data["q"]
    K = model["base_ring"]
    labels = model["labels"]
    index = model["label_index"]
    additive_character = model["additive_character"]
    W = Matrix(K, len(labels), len(labels), 0)
    w = weyl_matrix_gl2(F)

    for a in labels:
        for y in labels:
            coefficient = K(0)
            for s in F:
                g = diagonal_matrix_gl2(F, a, F(1)) * w * unipotent_matrix(F, s)
                coefficient += cuspidal_character_value(
                    g,
                    model["torus_character"],
                ) * additive_character_value(-y * s, additive_character)
            W[index[a * y], index[y]] = coefficient / K(q)

    return W


def direct_whittaker_matrix(M, model):
    F = model["data"]["F"]
    M = matrix_of(M)
    a = F(M[0, 0])
    b = F(M[0, 1])
    c = F(M[1, 0])
    d = F(M[1, 1])

    if c == 0:
        return upper_triangular_action_matrix(M, model)

    determinant = a * d - b * c
    left_unipotent = unipotent_action_matrix(a / c, model)
    diagonal = diagonal_action_matrix(-determinant / c, -c, model)
    right_unipotent = unipotent_action_matrix(d / c, model)

    return left_unipotent * diagonal * model["weyl_matrix"] * right_unipotent


def direct_character_rows(data, model):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "trace": direct_whittaker_matrix(
                    representative.matrix(),
                    model,
                ).trace(),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], str(row["trace"])))
    return rows


def verify_weyl_relations(model):
    F = model["data"]["F"]
    K = model["base_ring"]
    dimension = len(model["labels"])
    identity = identity_matrix(K, dimension)
    scalar = central_character_value(F(-1), model)

    assert model["weyl_matrix"] * model["weyl_matrix"] == scalar * identity


def verify_borel_relations(model):
    F = model["data"]["F"]

    for x in F:
        actual = direct_whittaker_matrix(unipotent_matrix(F, x), model)
        expected = unipotent_action_matrix(x, model)
        assert actual == expected

    for a in nonzero_field_elements(F):
        for d in nonzero_field_elements(F):
            actual = direct_whittaker_matrix(diagonal_matrix_gl2(F, a, d), model)
            expected = diagonal_action_matrix(a, d, model)
            assert actual == expected


def verify_representation_on_generators(data, model):
    generators = [generator.matrix() for generator in data["G"].gens()]
    generators.append(weyl_matrix_gl2(data["F"]))

    for g in generators:
        for h in generators:
            assert direct_whittaker_matrix(g * h, model) == (
                direct_whittaker_matrix(g, model)
                * direct_whittaker_matrix(h, model)
            )


def verify_representation_on_all_elements(data, model):
    matrices = [g.matrix() for g in data["G"]]
    representations = {
        tuple(g.list()): direct_whittaker_matrix(g, model)
        for g in matrices
    }

    for g in matrices:
        for h in matrices:
            assert representations[tuple((g * h).list())] == (
                representations[tuple(g.list())] * representations[tuple(h.list())]
            )


def verify_character_rows(data, model):
    formula_rows = cuspidal_character_rows(data, model["torus_character"])
    direct_rows = direct_character_rows(data, model)

    assert hermitian_inner_product(formula_rows, "trace", "trace") == 1

    for formula_row, direct_row in zip(formula_rows, direct_rows):
        assert formula_row["order"] == direct_row["order"]
        assert formula_row["size"] == direct_row["size"]
        assert formula_row["trace"] == direct_row["trace"]

    return direct_rows


def verify_direct_model(data, model, all_elements=False):
    verify_weyl_relations(model)
    verify_borel_relations(model)
    verify_representation_on_generators(data, model)
    if all_elements:
        verify_representation_on_all_elements(data, model)
    return verify_character_rows(data, model)


def print_character_table(rows):
    print("Cuspidal character by conjugacy class, computed in the direct model:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"trace = {row['trace']}"
        )


def print_generator_matrices(model):
    data = model["data"]

    print("Matrices for Sage's GL(2,q) generators in the direct Whittaker model:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = matrix_of(generator)
        print(f"generator {i} in GL_2(F_{data['q']}):")
        print(M)
        print(direct_whittaker_matrix(M, model))


def main():
    args = parse_args()
    data = build_projective_line_data(args.q)
    q = data["q"]
    exponent = args.cuspidal if args.cuspidal is not None else default_cuspidal_exponent(q)
    model = build_direct_whittaker_model(data, exponent)
    rows = verify_direct_model(data, model, all_elements=args.all_elements)

    print(f"Constructing direct Whittaker cuspidal model for GL_2(F_{q})")
    print(f"cuspidal exponent: {model['exponent']}")
    print(f"Whittaker labels: {[str(y) for y in model['labels']]}")
    print(f"dimension: {len(model['labels'])}")
    print("verified: Borel action is explicit on F_q^*")
    print("verified: Weyl matrix is recovered by finite Fourier inversion")
    print("verified: traces match the nonsplit-torus cuspidal character formula")
    if args.all_elements:
        print("verified: representation property on all element pairs")
    else:
        print("verified: representation property on generator pairs")
    print(
        f"<pi_{model['exponent']}, pi_{model['exponent']}>: "
        f"{hermitian_inner_product(rows, 'trace', 'trace')}"
    )

    if args.character_table:
        print_character_table(rows)

    if args.generators:
        print_generator_matrices(model)


if __name__ == "__main__":
    main()
