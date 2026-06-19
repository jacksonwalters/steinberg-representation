from argparse import ArgumentParser

from sage.all import *

from representation_utils_gl_3_q import (
    family_counts,
    family_square_sum,
    gl3_order,
    root_of_unity,
)
from representations_gl_2_q import primitive_element_by_search


def parse_args():
    parser = ArgumentParser(
        description=(
            "Organize the genuine cuspidal parameters for GL_3(F_q), coming "
            "from regular characters of the cubic nonsplit torus F_{q^3}^*."
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
        "--parameters",
        action="store_true",
        help="print all Frobenius orbits of regular character exponents",
    )
    parser.add_argument(
        "--character",
        type=int,
        metavar="M",
        help="construct the torus character theta_M of F_{q^3}^*",
    )
    return parser.parse_args()


def cubic_extension(F):
    R = PolynomialRing(F, "x")
    return F.extension(R.irreducible_element(3), "c")


def is_regular_cubic_exponent(q, exponent):
    modulus = q**3 - 1
    exponent %= modulus
    return exponent % (q**2 + q + 1) != 0


def cuspidal_parameter_orbits(q):
    modulus = q**3 - 1
    visited = set()
    orbits = []

    for exponent in range(modulus):
        if not is_regular_cubic_exponent(q, exponent):
            continue

        orbit = tuple(sorted({exponent % modulus, (q * exponent) % modulus, (q**2 * exponent) % modulus}))
        if orbit in visited:
            continue

        visited.add(orbit)
        orbits.append(orbit)

    return sorted(orbits)


def gl3_cuspidal_degree(q):
    return (q - 1) * (q**2 - 1)


def gl3_cuspidal_square_sum(q):
    return family_counts(q)["gl3_cuspidal"] * gl3_cuspidal_degree(q) ** 2


def cubic_torus_character_data(q, exponent, base_ring=None):
    F = GF(q, name="a")
    E = cubic_extension(F)
    modulus = q**3 - 1
    exponent %= modulus

    if not is_regular_cubic_exponent(q, exponent):
        raise ValueError(
            "GL_3 cuspidal parameters must be regular; exponents divisible by "
            "q^2 + q + 1 factor through the norm to F_q^*"
        )

    K = base_ring if base_ring is not None else CyclotomicField(modulus)
    primitive = primitive_element_by_search(E)

    return {
        "q": q,
        "F": F,
        "E": E,
        "modulus": modulus,
        "exponent": exponent,
        "base_ring": K,
        "primitive": primitive,
        "log": {primitive**k: k for k in range(modulus)},
        "zeta": root_of_unity(K, modulus),
    }


def cubic_torus_character_value(x, character_data):
    E = character_data["E"]
    x = E(x)
    if x == 0:
        raise ValueError("torus characters are only defined on nonzero elements")

    return character_data["zeta"] ** (
        character_data["exponent"] * character_data["log"][x]
    )


def print_parameter_orbits(q):
    orbits = cuspidal_parameter_orbits(q)
    print(f"regular character exponents for F_{q**3}^* are taken modulo {q**3 - 1}")
    print("Frobenius acts by m -> q*m, so cuspidal parameters are size-three orbits")
    print(f"excluded exponents: multiples of q^2 + q + 1 = {q**2 + q + 1}")
    for orbit in orbits:
        print(f"  {orbit}")
    print(f"number of GL_3 cuspidal parameters: {len(orbits)}")
    assert len(orbits) == family_counts(q)["gl3_cuspidal"]


def main():
    args = parse_args()
    q = args.q
    orbits = cuspidal_parameter_orbits(q)
    degree = gl3_cuspidal_degree(q)

    print(f"Organizing genuine cuspidal parameters for GL_3(F_{q})")
    print(f"cubic nonsplit torus: F_{q**3}^*")
    print(f"regular Frobenius orbits: {len(orbits)}")
    print(f"cuspidal representation degree: {degree} = (q - 1)(q^2 - 1)")
    print(f"cuspidal square contribution: {gl3_cuspidal_square_sum(q)}")
    print(f"full family square sum: {family_square_sum(q)}")
    print(f"|GL_3(F_{q})|: {gl3_order(q)}")
    assert family_square_sum(q) == gl3_order(q)

    if args.parameters:
        print_parameter_orbits(q)

    if args.character is not None:
        character_data = cubic_torus_character_data(q, args.character)
        print(f"constructed theta_{character_data['exponent']} on F_{q**3}^*")
        print(f"primitive torus generator: {character_data['primitive']}")
        print(
            "theta(generator) = "
            f"{cubic_torus_character_value(character_data['primitive'], character_data)}"
        )


if __name__ == "__main__":
    main()
