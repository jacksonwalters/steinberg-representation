from argparse import ArgumentParser

from sage.all import *

from representation_utils_gl_3_q import (
    build_gl3_data,
    build_induced_representation,
    character_rows,
    coerce_matrix_to_field,
    determinant_character_value,
    hermitian_inner_product,
    induced_dimension,
    induced_matrix,
    levi_2_1_blocks,
    multiplicative_character_data,
    multiplicative_character_value,
    parabolic_2_1_point_coset_data,
    print_character_table,
    scalar_matrix,
    verify_induced_representation,
)
from projective_line import build_projective_line_data, steinberg_matrix as gl2_steinberg_matrix


CHARACTER_AUTO_Q_LIMIT = 4


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct the two irreducible constituents attached to a split "
            "torus character triple (chi, chi, mu) with chi != mu."
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
        "--repeated",
        type=int,
        default=0,
        metavar="A",
        help="exponent A for the repeated character chi_A",
    )
    parser.add_argument(
        "--single",
        type=int,
        default=1,
        metavar="B",
        help="exponent B for the single character chi_B",
    )
    parser.add_argument(
        "--which",
        choices=["small", "large", "both"],
        default="both",
        help="which repeated-split constituent to construct",
    )
    parser.add_argument(
        "--character",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute characters and verify irreducibility",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print character values by conjugacy class",
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


def validate_exponents(data, repeated_exponent, single_exponent):
    order = data["q"] - 1
    if order == 1:
        raise ValueError("there are no distinct split characters for q = 2")

    repeated_exponent %= order
    single_exponent %= order
    if repeated_exponent == single_exponent:
        raise ValueError("the repeated and single characters must be distinct")

    return repeated_exponent, single_exponent


def build_repeated_split_small(data, repeated_exponent, single_exponent):
    repeated_exponent, single_exponent = validate_exponents(
        data,
        repeated_exponent,
        single_exponent,
    )
    character_data = multiplicative_character_data(data["F"])
    K = character_data["base_ring"]

    def coefficient_matrix(p):
        A, d = levi_2_1_blocks(p, data["F"])
        value = determinant_character_value(A, repeated_exponent, character_data)
        value *= multiplicative_character_value(d, single_exponent, character_data)
        return scalar_matrix(value, K)

    induced = build_induced_representation(
        data=data,
        subgroup_matrices=None,
        coefficient_matrix=coefficient_matrix,
        coefficient_dimension=1,
        base_ring=K,
        name="repeated_split_small",
        coset_data=parabolic_2_1_point_coset_data(data),
    )
    induced["repeated_exponent"] = repeated_exponent
    induced["single_exponent"] = single_exponent
    induced["character_data"] = character_data
    return induced


def build_repeated_split_large(data, repeated_exponent, single_exponent):
    repeated_exponent, single_exponent = validate_exponents(
        data,
        repeated_exponent,
        single_exponent,
    )
    character_data = multiplicative_character_data(data["F"])
    K = character_data["base_ring"]
    gl2_data = build_projective_line_data(data["q"])

    def coefficient_matrix(p):
        A, d = levi_2_1_blocks(p, data["F"])
        A_gl2 = coerce_matrix_to_field(A, gl2_data["F"])
        gl2_matrix = gl2_steinberg_matrix(A_gl2, gl2_data, base_ring=K)
        scalar = determinant_character_value(A, repeated_exponent, character_data)
        scalar *= multiplicative_character_value(d, single_exponent, character_data)
        return scalar * gl2_matrix

    induced = build_induced_representation(
        data=data,
        subgroup_matrices=None,
        coefficient_matrix=coefficient_matrix,
        coefficient_dimension=data["q"],
        base_ring=K,
        name="repeated_split_large",
        coset_data=parabolic_2_1_point_coset_data(data),
    )
    induced["repeated_exponent"] = repeated_exponent
    induced["single_exponent"] = single_exponent
    induced["character_data"] = character_data
    induced["gl2_data"] = gl2_data
    return induced


def repeated_split_matrix(M, model):
    return induced_matrix(M, model)


def repeated_split_character_rows(data, model):
    return character_rows(data, lambda M: repeated_split_matrix(M, model))


def construct_models(data, repeated_exponent, single_exponent, which):
    models = []
    if which in ("small", "both"):
        models.append(build_repeated_split_small(data, repeated_exponent, single_exponent))
    if which in ("large", "both"):
        models.append(build_repeated_split_large(data, repeated_exponent, single_exponent))
    return models


def print_generator_matrices(model):
    data = model["data"]
    print(f"Matrices for Sage's GL(3,q) generators in {model['name']}:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_3(F_{data['q']}):")
        print(M)
        print(repeated_split_matrix(M, model))


def main():
    args = parse_args()
    data = build_gl3_data(args.q)
    models = construct_models(data, args.repeated, args.single, args.which)

    print(f"Constructing repeated-split representations of GL_3(F_{data['q']})")
    for model in models:
        rows = None
        if should_run_character(args.character, data["q"]) or args.character_table:
            rows = repeated_split_character_rows(data, model)

        verify_induced_representation(model, rows)

        print(
            f"{model['name']}: repeated chi_{model['repeated_exponent']}, "
            f"single chi_{model['single_exponent']}"
        )
        print(f"  cosets G/P_{{2,1}}: {len(model['representatives'])}")
        print(f"  dimension: {induced_dimension(model)}")
        print("  verified: induced action respects multiplication")

        if rows is not None:
            print(f"  character conjugacy classes: {len(rows)}")
            print(f"  <pi, pi>: {hermitian_inner_product(rows, 'trace', 'trace')}")
        else:
            print("  character: skipped (use --character yes to compute)")

        if args.character_table:
            print_character_table(rows, model["name"])

        if args.generators:
            print_generator_matrices(model)


if __name__ == "__main__":
    main()
