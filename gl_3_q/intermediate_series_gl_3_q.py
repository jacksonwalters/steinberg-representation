from argparse import ArgumentParser

from sage.all import *

from representation_utils_gl_3_q import (
    build_gl3_data,
    build_induced_representation,
    character_rows,
    hermitian_inner_product,
    induced_dimension,
    induced_matrix,
    levi_2_1_blocks,
    multiplicative_character_data,
    multiplicative_character_value,
    parabolic_2_1_point_coset_data,
    print_character_table,
    verify_induced_representation,
)

from direct_whittaker_cuspidal_gl_2_q import (
    build_direct_whittaker_model,
    default_cuspidal_exponent,
    direct_whittaker_matrix,
    verify_direct_model,
)
from projective_line import build_projective_line_data


CHARACTER_AUTO_Q_LIMIT = 3


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct the GL_3(F_q) intermediate series induced from a "
            "cuspidal representation of the Levi GL_2(F_q) x GL_1(F_q)."
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
            "regular exponent M for the GL_2 cuspidal character of "
            "F_{q^2}^*; defaults to the first orbit representative"
        ),
    )
    parser.add_argument(
        "--character-twist",
        type=int,
        default=0,
        metavar="A",
        help="exponent A for the GL_1(F_q) character on the second Levi factor",
    )
    parser.add_argument(
        "--character",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute the character and verify irreducibility",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print intermediate-series values by conjugacy class",
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


def build_intermediate_series(data, cuspidal_exponent, character_exponent=0):
    gl2_data = build_projective_line_data(data["q"])
    if cuspidal_exponent is None:
        cuspidal_exponent = default_cuspidal_exponent(data["q"])

    gl2_model = build_direct_whittaker_model(gl2_data, cuspidal_exponent)
    verify_direct_model(gl2_data, gl2_model)

    K = gl2_model["base_ring"]
    character_data = multiplicative_character_data(data["F"], base_ring=K)
    character_exponent %= character_data["order"]

    def coefficient_matrix(p):
        A, d = levi_2_1_blocks(p, data["F"])
        scalar = multiplicative_character_value(d, character_exponent, character_data)
        return scalar * direct_whittaker_matrix(A, gl2_model)

    induced = build_induced_representation(
        data=data,
        subgroup_matrices=None,
        coefficient_matrix=coefficient_matrix,
        coefficient_dimension=data["q"] - 1,
        base_ring=K,
        name="gl2_cuspidal_intermediate",
        coset_data=parabolic_2_1_point_coset_data(data),
    )
    induced["gl2_data"] = gl2_data
    induced["gl2_model"] = gl2_model
    induced["cuspidal_exponent"] = gl2_model["exponent"]
    induced["character_exponent"] = character_exponent
    induced["character_data"] = character_data
    return induced


def intermediate_series_matrix(M, intermediate_series):
    return induced_matrix(M, intermediate_series)


def intermediate_series_character_rows(data, intermediate_series):
    return character_rows(data, lambda M: intermediate_series_matrix(M, intermediate_series))


def print_generator_matrices(intermediate_series):
    data = intermediate_series["data"]
    print("Matrices for Sage's GL(3,q) generators in the intermediate series:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_3(F_{data['q']}):")
        print(M)
        print(intermediate_series_matrix(M, intermediate_series))


def main():
    args = parse_args()
    data = build_gl3_data(args.q)
    intermediate_series = build_intermediate_series(
        data,
        args.cuspidal,
        character_exponent=args.character_twist,
    )

    rows = None
    if should_run_character(args.character, data["q"]) or args.character_table:
        rows = intermediate_series_character_rows(data, intermediate_series)

    verify_induced_representation(intermediate_series, rows)

    q = data["q"]
    print(f"Constructing the GL_2-cuspidal intermediate series of GL_3(F_{q})")
    print(
        "Levi representation: "
        f"pi_{intermediate_series['cuspidal_exponent']} tensor "
        f"chi_{intermediate_series['character_exponent']}"
    )
    print(f"GL_2 cuspidal dimension: {q - 1}")
    print(f"cosets G/P_{{2,1}}: {len(intermediate_series['representatives'])}")
    print(f"dimension: {induced_dimension(intermediate_series)} = q^3 - 1")
    print("verified: GL_2 cuspidal Whittaker model")
    print("verified: induced action respects multiplication")

    if rows is not None:
        print(f"character conjugacy classes: {len(rows)}")
        print(f"<pi, pi>: {hermitian_inner_product(rows, 'trace', 'trace')}")
    else:
        print("character: skipped (use --character yes to compute)")

    if args.character_table:
        print_character_table(rows, "Intermediate-series")

    if args.generators:
        print_generator_matrices(intermediate_series)


if __name__ == "__main__":
    main()
