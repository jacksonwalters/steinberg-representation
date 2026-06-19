from argparse import ArgumentParser
from itertools import combinations

from sage.all import *

from representation_utils_gl_3_q import (
    build_gl3_data,
    build_induced_representation,
    borel_flag_coset_data,
    character_rows,
    hermitian_inner_product,
    induced_dimension,
    induced_matrix,
    multiplicative_character_data,
    multiplicative_character_value,
    print_character_table,
    scalar_matrix,
    verify_induced_representation,
)


CHARACTER_AUTO_Q_LIMIT = 4


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct the irreducible principal series of GL_3(F_q) from "
            "three distinct split-torus characters."
        )
    )
    parser.add_argument(
        "q",
        nargs="?",
        type=int,
        default=4,
        help="prime power q for the finite field F_q",
    )
    parser.add_argument(
        "--characters",
        nargs=3,
        type=int,
        metavar=("A", "B", "C"),
        help=(
            "construct Ind_B^G(chi_A tensor chi_B tensor chi_C), where chi_E "
            "sends a multiplicative generator of F_q^* to zeta_(q-1)^E"
        ),
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
        help="print principal-series values by conjugacy class",
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


def principal_series_triples(q):
    return combinations(range(q - 1), 3)


def default_character_triple(q):
    try:
        return next(principal_series_triples(q))
    except StopIteration as exc:
        raise ValueError(
            "GL_3(F_q) has no three-distinct-character principal series unless q >= 4"
        ) from exc


def build_principal_series(data, exponents):
    character_data = multiplicative_character_data(data["F"])
    order = character_data["order"]
    exponents = tuple(exponent % order for exponent in exponents)

    if len(set(exponents)) != 3:
        raise ValueError("principal series in this module require three distinct characters")

    K = character_data["base_ring"]

    def coefficient_matrix(b):
        value = K(1)
        for i, exponent in enumerate(exponents):
            value *= multiplicative_character_value(b[i, i], exponent, character_data)
        return scalar_matrix(value, K)

    induced = build_induced_representation(
        data=data,
        subgroup_matrices=None,
        coefficient_matrix=coefficient_matrix,
        coefficient_dimension=1,
        base_ring=K,
        name="principal_series",
        coset_data=borel_flag_coset_data(data),
    )
    induced["exponents"] = exponents
    induced["character_data"] = character_data
    return induced


def principal_series_matrix(M, principal_series):
    return induced_matrix(M, principal_series)


def principal_series_character_rows(data, principal_series):
    return character_rows(data, lambda M: principal_series_matrix(M, principal_series))


def print_generator_matrices(principal_series):
    data = principal_series["data"]
    print("Matrices for Sage's GL(3,q) generators in the principal series:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_3(F_{data['q']}):")
        print(M)
        print(principal_series_matrix(M, principal_series))


def main():
    args = parse_args()
    data = build_gl3_data(args.q)

    if args.characters is None:
        try:
            exponents = default_character_triple(data["q"])
        except ValueError as exc:
            print(exc)
            print("Try q = 4 or pass --characters A B C with three distinct exponents.")
            return
    else:
        exponents = tuple(args.characters)

    principal_series = build_principal_series(data, exponents)

    rows = None
    if should_run_character(args.character, data["q"]) or args.character_table:
        rows = principal_series_character_rows(data, principal_series)

    verify_induced_representation(principal_series, rows)

    q = data["q"]
    print(f"Constructing a principal series of GL_3(F_{q})")
    print(
        "representation: "
        f"Ind_B^G(chi_{principal_series['exponents'][0]} tensor "
        f"chi_{principal_series['exponents'][1]} tensor "
        f"chi_{principal_series['exponents'][2]})"
    )
    print(f"cosets G/B: {len(principal_series['representatives'])}")
    print(f"dimension: {induced_dimension(principal_series)} = (q + 1)(q^2 + q + 1)")
    print("verified: induced action respects multiplication")

    if rows is not None:
        print(f"character conjugacy classes: {len(rows)}")
        print(f"<pi, pi>: {hermitian_inner_product(rows, 'trace', 'trace')}")
    else:
        print("character: skipped (use --character yes to compute)")

    if args.character_table:
        print_character_table(rows, "Principal-series")

    if args.generators:
        print_generator_matrices(principal_series)


if __name__ == "__main__":
    main()
