from argparse import ArgumentParser

from representation_utils_gl_3_q import build_gl3_data, print_family_summary


def parse_args():
    parser = ArgumentParser(
        description=(
            "Organize the standard complex irreducible families of GL_3(F_q) "
            "and point to the modular construction scripts."
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
        "--scripts",
        action="store_true",
        help="print the construction scripts corresponding to the family table",
    )
    parser.add_argument(
        "--cuspidal-parameters",
        action="store_true",
        help="print the genuine GL_3 cuspidal parameter orbits",
    )
    return parser.parse_args()


def print_script_map():
    print("construction scripts:")
    print("  unipotent twists:")
    print("    sage gl_3_q/unipotent_representations_gl_3_q.py q")
    print("  three-distinct-character principal series:")
    print("    sage gl_3_q/principal_series_gl_3_q.py q --characters A B C")
    print("  repeated split characters (chi, chi, mu):")
    print("    sage gl_3_q/repeated_split_representations_gl_3_q.py q --repeated A --single B")
    print("  GL_2-cuspidal intermediate series:")
    print("    sage gl_3_q/intermediate_series_gl_3_q.py q --cuspidal M --character-twist A")
    print("  genuine GL_3 cuspidal parameters:")
    print("    sage gl_3_q/cuspidal_parameters_gl_3_q.py q --parameters")
    print("  Steinberg from the Tits building:")
    print("    sage gl_3_q/steinberg_representation_gl_3_q.py q")


def main():
    args = parse_args()
    data = build_gl3_data(args.q)

    print(f"Organizing representations of GL_3(F_{data['q']})")
    print(f"group order: {data['G'].order()}")
    print_family_summary(data["q"])

    if args.scripts:
        print_script_map()

    if args.cuspidal_parameters:
        from cuspidal_parameters_gl_3_q import print_parameter_orbits

        print_parameter_orbits(data["q"])


if __name__ == "__main__":
    main()
