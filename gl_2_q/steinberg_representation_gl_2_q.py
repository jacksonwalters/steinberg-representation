from argparse import ArgumentParser

from projective_line import (
    build_projective_line_data,
    character_rows,
    inner_product,
    point_repr,
    projective_line_permutation_matrix,
    steinberg_matrix,
    verify_projective_line_data,
)


CHARACTER_AUTO_Q_LIMIT = 13


def parse_args():
    parser = ArgumentParser(
        description="Construct the Steinberg representation of GL_2(F_q)."
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
        help="compute the Steinberg character by fixed projective-point counts",
    )
    parser.add_argument(
        "--character-table",
        action="store_true",
        help="print one row for each conjugacy class when computing the character",
    )
    parser.add_argument(
        "--points",
        action="store_true",
        help="print the projective line P^1(F_q)",
    )
    parser.add_argument(
        "--generators",
        action="store_true",
        help="print matrices for Sage's group generators",
    )
    return parser.parse_args()


def should_run_character(mode, q):
    if mode == "yes":
        return True
    if mode == "no":
        return False
    return q <= CHARACTER_AUTO_Q_LIMIT


def print_points(data):
    print("projective points:")
    for i, point in enumerate(data["points"], start=1):
        print(f"  P{i}: {point_repr(point)}")


def print_character_table(rows):
    print("Steinberg character by conjugacy class:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"fixed points = {row['fixed_points']:3}, "
            f"chi_St = {row['steinberg']:4}"
        )


def print_generator_matrices(data):
    print("Matrices for Sage's GL(2,q) generators:")

    for i, generator in enumerate(data["G"].gens(), start=1):
        M = generator.matrix()
        print(f"generator {i} in GL_2(F_{data['q']}):")
        print(M)
        print("  permutation module C[P^1]:")
        print(projective_line_permutation_matrix(M, data))
        print("  Steinberg module St:")
        print(steinberg_matrix(M, data))


def main():
    args = parse_args()
    data = build_projective_line_data(args.q)
    q = data["q"]

    rows = None
    if should_run_character(args.character, q) or args.character_table:
        rows = character_rows(data)

    verify_projective_line_data(data, rows)

    print(f"Computing St(GL_2(F_{q}))")
    print(f"group order: {data['G'].order()}")
    print(f"projective line points: {len(data['points'])} (expected {q + 1})")
    print(f"rank of sum map C[P^1] -> C: {data['sum_map'].rank()}")
    print(f"dim ker(sum) = dim Steinberg module: {data['steinberg'].dimension()}")
    print("permutation module decomposition: C[P^1] = 1 + St")

    if rows is not None:
        print(f"character conjugacy classes: {len(rows)}")
        print(f"character inner product <St, St>: {inner_product(rows, 'steinberg', 'steinberg')}")
        print(
            "permutation character inner product "
            f"<C[P^1], C[P^1]>: {inner_product(rows, 'permutation', 'permutation')}"
        )
    else:
        print("character: skipped (use --character yes to compute)")

    if args.points:
        print_points(data)

    if args.character_table:
        print_character_table(rows)

    if args.generators:
        print_generator_matrices(data)


if __name__ == "__main__":
    main()
