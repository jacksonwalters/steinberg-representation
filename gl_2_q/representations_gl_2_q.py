from argparse import ArgumentParser

from sage.all import *

from projective_line import (
    build_projective_line_data,
    character_rows as projective_character_rows,
    fixed_projective_points,
    inner_product,
    matrix_of,
    projective_line_permutation_matrix,
    steinberg_matrix,
    steinberg_trace,
    verify_projective_line_data,
)


CHARACTER_AUTO_Q_LIMIT = 13
PRINCIPAL_SERIES_AUTO_GROUP_ORDER_LIMIT = 50000


def parse_args():
    parser = ArgumentParser(
        description=(
            "Construct basic GL_2(F_q) representations and organize the "
            "standard complex irreducible families."
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
        "--character",
        choices=["auto", "yes", "no"],
        default="auto",
        help="compute the projective-line and Steinberg characters",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print projective-line and Steinberg values by conjugacy class",
    )
    parser.add_argument(
        "--principal-series",
        nargs=2,
        type=int,
        metavar=("A", "B"),
        help=(
            "construct Ind_B^G(chi_A tensor chi_B), where chi_E sends a "
            "multiplicative generator of F_q^* to zeta_(q-1)^E"
        ),
    )
    parser.add_argument(
        "--all-principal-series",
        action="store_true",
        help="construct every irreducible principal series and check its character norm",
    )
    parser.add_argument(
        "--cuspidal-parameters",
        action="store_true",
        help="print cuspidal parameter orbits for characters of F_{q^2}^*",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's group generators in constructed modules",
    )
    return parser.parse_args()


def should_run_character(mode, q):
    if mode == "yes":
        return True
    if mode == "no":
        return False
    return q <= CHARACTER_AUTO_Q_LIMIT


def family_counts(q):
    return {
        "linear": q - 1,
        "steinberg_twists": q - 1,
        "principal_series": (q - 1) * (q - 2) // 2,
        "cuspidal": q * (q - 1) // 2,
    }


def family_square_sum(q):
    counts = family_counts(q)
    return (
        counts["linear"]
        + counts["steinberg_twists"] * q**2
        + counts["principal_series"] * (q + 1) ** 2
        + counts["cuspidal"] * (q - 1) ** 2
    )


def multiplicative_character_data(F):
    order = F.order() - 1

    if order == 1:
        return {
            "order": order,
            "base_ring": QQ,
            "generator": F(1),
            "log": {F(1): 0},
            "zeta": QQ(1),
        }

    generator = F.multiplicative_generator()
    logs = {generator**k: k for k in range(order)}
    K = CyclotomicField(order)

    return {
        "order": order,
        "base_ring": K,
        "generator": generator,
        "log": logs,
        "zeta": K.gen(),
    }


def multiplicative_character_value(a, exponent, character_data):
    if a == 0:
        raise ValueError("multiplicative characters are only defined on F_q^*")

    order = character_data["order"]
    if order == 1:
        return QQ(1)

    exponent = exponent % order
    return character_data["zeta"] ** (exponent * character_data["log"][a])


def determinant_character_value(M, exponent, character_data):
    return multiplicative_character_value(matrix_of(M).det(), exponent, character_data)


def one_dimensional_character(M, exponent, character_data):
    return determinant_character_value(M, exponent, character_data)


def steinberg_twist_matrix(M, exponent, data, character_data=None):
    if character_data is None:
        character_data = multiplicative_character_data(data["F"])

    scalar = determinant_character_value(M, exponent, character_data)
    return scalar * steinberg_matrix(M, data, base_ring=character_data["base_ring"])


def steinberg_twist_trace(M, exponent, data, character_data=None):
    if character_data is None:
        character_data = multiplicative_character_data(data["F"])

    return determinant_character_value(M, exponent, character_data) * steinberg_trace(M, data)


def matrix_key(M):
    return tuple(matrix_of(M).list())


def upper_borel_matrices(data):
    return [
        g.matrix()
        for g in data["G"]
        if g.matrix()[1, 0] == 0
    ]


def left_coset_data(data):
    all_matrices = [g.matrix() for g in data["G"]]
    borel = upper_borel_matrices(data)
    remaining = {matrix_key(M) for M in all_matrices}
    representatives = []
    decomposition = {}

    for M in all_matrices:
        key = matrix_key(M)
        if key not in remaining:
            continue

        representatives.append(M)
        coset_index = len(representatives) - 1

        for b in borel:
            X = b * M
            X_key = matrix_key(X)
            decomposition[X_key] = (coset_index, b)
            remaining.discard(X_key)

    assert len(representatives) == data["q"] + 1
    assert not remaining

    return {
        "borel": borel,
        "representatives": representatives,
        "decomposition": decomposition,
    }


def principal_series_pairs(q):
    for a in range(q - 1):
        for b in range(a + 1, q - 1):
            yield a, b


def build_principal_series(data, exponent_a, exponent_b, cosets=None):
    character_data = multiplicative_character_data(data["F"])
    order = character_data["order"]

    if order == 1 or exponent_a % order == exponent_b % order:
        raise ValueError("principal series is irreducible here only for distinct characters")

    if cosets is None:
        cosets = left_coset_data(data)

    return {
        "data": data,
        "exponent_a": exponent_a % order,
        "exponent_b": exponent_b % order,
        "character_data": character_data,
        "base_ring": character_data["base_ring"],
        "representatives": cosets["representatives"],
        "decomposition": cosets["decomposition"],
    }


def borel_character(b, principal_series):
    character_data = principal_series["character_data"]
    return (
        multiplicative_character_value(
            b[0, 0],
            principal_series["exponent_a"],
            character_data,
        )
        * multiplicative_character_value(
            b[1, 1],
            principal_series["exponent_b"],
            character_data,
        )
    )


def principal_series_matrix(M, principal_series):
    M = matrix_of(M)
    representatives = principal_series["representatives"]
    decomposition = principal_series["decomposition"]
    K = principal_series["base_ring"]
    A = Matrix(K, len(representatives), len(representatives), 0)

    for i, representative in enumerate(representatives):
        j, b = decomposition[matrix_key(representative * M)]
        A[i, j] = borel_character(b, principal_series)

    return A


def principal_series_character_rows(data, principal_series):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "trace": principal_series_matrix(
                    representative.matrix(),
                    principal_series,
                ).trace(),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], str(row["trace"])))
    return rows


def conjugate_value(x):
    return x.conjugate() if hasattr(x, "conjugate") else x


def hermitian_inner_product(rows, left_key, right_key):
    group_order = sum(ZZ(row["size"]) for row in rows)
    return (
        sum(
            ZZ(row["size"])
            * row[left_key]
            * conjugate_value(row[right_key])
            for row in rows
        )
        / group_order
    )


def verify_principal_series(data, principal_series, rows=None):
    G = data["G"]
    g = G.random_element().matrix()
    h = G.random_element().matrix()

    assert principal_series_matrix(g * h, principal_series) == (
        principal_series_matrix(g, principal_series)
        * principal_series_matrix(h, principal_series)
    )
    assert principal_series_matrix(G.one().matrix(), principal_series).trace() == data["q"] + 1

    if rows is not None:
        assert hermitian_inner_product(rows, "trace", "trace") == 1


def cuspidal_parameter_orbits(q):
    modulus = q**2 - 1
    fixed_step = q + 1
    visited = set()
    orbits = []

    for exponent in range(modulus):
        if exponent % fixed_step == 0:
            continue

        orbit = tuple(sorted({exponent % modulus, (q * exponent) % modulus}))
        if orbit in visited:
            continue

        visited.add(orbit)
        orbits.append(orbit)

    return sorted(orbits)


def print_family_summary(data):
    q = data["q"]
    counts = family_counts(q)
    group_order = data["G"].order()
    square_sum = family_square_sum(q)

    print("irreducible-family counts over C:")
    print(f"  determinant characters: {counts['linear']} of dimension 1")
    print(f"  Steinberg twists:       {counts['steinberg_twists']} of dimension {q}")
    print(f"  principal series:       {counts['principal_series']} of dimension {q + 1}")
    print(f"  cuspidal series:        {counts['cuspidal']} of dimension {q - 1}")
    print(f"sum of squares of dimensions: {square_sum}")
    print(f"group order: {group_order}")
    assert square_sum == group_order


def print_projective_character_table(rows):
    print("Projective-line and Steinberg characters by conjugacy class:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"fixed P^1 = {row['fixed_points']:3}, "
            f"chi_perm = {row['permutation']:4}, "
            f"chi_St = {row['steinberg']:4}"
        )


def print_cuspidal_parameters(q):
    orbits = cuspidal_parameter_orbits(q)
    print(f"cuspidal parameter orbits for characters of F_{q**2}^*:")
    print(f"  characters are indexed by exponents modulo {q**2 - 1}")
    print(f"  fixed exponents modulo {q + 1} factor through the norm and are excluded")
    print("  remaining exponents are paired by m ~ q*m")
    for orbit in orbits:
        print(f"  {orbit}")
    print(f"number of cuspidal parameters: {len(orbits)}")
    assert len(orbits) == family_counts(q)["cuspidal"]


def print_generator_matrices(data, principal_series=None):
    character_data = multiplicative_character_data(data["F"])

    print("Matrices for Sage's GL(2,q) generators:")
    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_2(F_{data['q']}):")
        print(M)
        print("  projective-line permutation module C[P^1]:")
        print(projective_line_permutation_matrix(M, data))
        print("  Steinberg module St:")
        print(steinberg_matrix(M, data))
        print("  first determinant twist of St:")
        print(steinberg_twist_matrix(M, 1, data, character_data))

        if principal_series is not None:
            print(
                "  principal series "
                f"Ind_B^G(chi_{principal_series['exponent_a']} tensor "
                f"chi_{principal_series['exponent_b']}):"
            )
            print(principal_series_matrix(M, principal_series))


def construct_requested_principal_series(args, data):
    counts = family_counts(data["q"])
    if counts["principal_series"] == 0:
        return None

    if args.principal_series is not None:
        exponent_a, exponent_b = args.principal_series
        return build_principal_series(data, exponent_a, exponent_b)

    if data["G"].order() <= PRINCIPAL_SERIES_AUTO_GROUP_ORDER_LIMIT:
        exponent_a, exponent_b = next(principal_series_pairs(data["q"]))
        return build_principal_series(data, exponent_a, exponent_b)

    return None


def main():
    args = parse_args()
    data = build_projective_line_data(args.q)
    q = data["q"]

    projective_rows = None
    if should_run_character(args.character, q) or args.character_table:
        projective_rows = projective_character_rows(data)

    verify_projective_line_data(data, projective_rows)

    print(f"Organizing representations of GL_2(F_{q})")
    print(f"projective line size: {len(data['points'])}")
    print(f"dim St from reduced H_0(P^1): {data['steinberg'].dimension()}")
    print("geometric decomposition: C[P^1] = 1 + St")
    if projective_rows is not None:
        print(f"character conjugacy classes: {len(projective_rows)}")
        print(f"<St, St>: {inner_product(projective_rows, 'steinberg', 'steinberg')}")
    else:
        print("projective-line character: skipped (use --character yes to compute)")

    print_family_summary(data)

    principal_series = construct_requested_principal_series(args, data)
    if principal_series is not None:
        ps_rows = principal_series_character_rows(data, principal_series)
        verify_principal_series(data, principal_series, ps_rows)
        print(
            "constructed principal series: "
            f"Ind_B^G(chi_{principal_series['exponent_a']} tensor "
            f"chi_{principal_series['exponent_b']})"
        )
        print(f"principal-series dimension: {q + 1}")
        print(f"<principal series, itself>: {hermitian_inner_product(ps_rows, 'trace', 'trace')}")
    elif family_counts(q)["principal_series"] == 0:
        print("principal series: none for q = 2")
    else:
        print(
            "principal series: skipped by default for this group size "
            "(use --principal-series A B)"
        )

    print(
        "cuspidal representations: counted by nonsplit-torus character "
        "orbits; matrix construction is not included in this script"
    )

    if args.all_principal_series and family_counts(q)["principal_series"] > 0:
        cosets = left_coset_data(data)
        print("checking all irreducible principal series:")
        for exponent_a, exponent_b in principal_series_pairs(q):
            ps = build_principal_series(data, exponent_a, exponent_b, cosets=cosets)
            rows = principal_series_character_rows(data, ps)
            verify_principal_series(data, ps, rows)
            print(
                f"  chi_{exponent_a}, chi_{exponent_b}: "
                f"<pi, pi> = {hermitian_inner_product(rows, 'trace', 'trace')}"
            )

    if args.character_table:
        print_projective_character_table(projective_rows)

    if args.cuspidal_parameters:
        print_cuspidal_parameters(q)

    if args.generators:
        print_generator_matrices(data, principal_series)


if __name__ == "__main__":
    main()
