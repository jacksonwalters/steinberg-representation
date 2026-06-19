from pathlib import Path
import sys

from sage.all import *


REPO_ROOT = Path(__file__).resolve().parents[1]
GL2_DIR = REPO_ROOT / "gl_2_q"
if str(GL2_DIR) not in sys.path:
    sys.path.insert(0, str(GL2_DIR))


def matrix_of(g):
    return g.matrix() if hasattr(g, "matrix") else g


def matrix_key(M):
    return tuple(matrix_of(M).list())


def build_gl3_data(q):
    if q < 2:
        raise ValueError("q must be a prime power greater than 1")

    F = GF(q, name="a")
    return {
        "q": q,
        "F": F,
        "V": VectorSpace(F, 3),
        "G": GL(3, F),
    }


def gl3_order(q):
    return q**3 * (q - 1) ** 3 * (q + 1) * (q**2 + q + 1)


def family_counts(q):
    r = q - 1
    return {
        "linear": r,
        "middle_unipotent": r,
        "steinberg": r,
        "repeated_split_small": r * (r - 1),
        "repeated_split_large": r * (r - 1),
        "principal_series": binomial(r, 3),
        "gl2_cuspidal_intermediate": q * r**2 // 2,
        "gl3_cuspidal": q * (q - 1) * (q + 1) // 3,
    }


def family_square_sum(q):
    counts = family_counts(q)
    D = (q + 1) * (q**2 + q + 1)

    return (
        counts["linear"]
        + counts["middle_unipotent"] * (q * (q + 1)) ** 2
        + counts["steinberg"] * q**6
        + counts["repeated_split_small"] * (q**2 + q + 1) ** 2
        + counts["repeated_split_large"] * (q * (q**2 + q + 1)) ** 2
        + counts["principal_series"] * D**2
        + counts["gl2_cuspidal_intermediate"] * (q**3 - 1) ** 2
        + counts["gl3_cuspidal"] * ((q - 1) * (q**2 - 1)) ** 2
    )


def print_family_summary(q):
    counts = family_counts(q)

    print(f"irreducible-family counts for GL_3(F_{q}) over C:")
    print(f"  determinant characters:             {counts['linear']} of dimension 1")
    print(
        "  middle unipotent twists:           "
        f"{counts['middle_unipotent']} of dimension {q * (q + 1)}"
    )
    print(f"  Steinberg twists:                   {counts['steinberg']} of dimension {q**3}")
    print(
        "  repeated split, small constituent: "
        f"{counts['repeated_split_small']} of dimension {q**2 + q + 1}"
    )
    print(
        "  repeated split, large constituent: "
        f"{counts['repeated_split_large']} of dimension {q * (q**2 + q + 1)}"
    )
    print(
        "  three distinct split characters:   "
        f"{counts['principal_series']} of dimension {(q + 1) * (q**2 + q + 1)}"
    )
    print(
        "  GL_2 cuspidal intermediate series: "
        f"{counts['gl2_cuspidal_intermediate']} of dimension {q**3 - 1}"
    )
    print(
        "  GL_3 cuspidal series:              "
        f"{counts['gl3_cuspidal']} of dimension {(q - 1) * (q**2 - 1)}"
    )
    print(f"sum of squares of dimensions: {family_square_sum(q)}")
    print(f"|GL_3(F_{q})|: {gl3_order(q)}")
    assert family_square_sum(q) == gl3_order(q)


def nonzero_field_elements(F):
    return [x for x in F if x != 0]


def root_of_unity(base_ring, order):
    if order == 1:
        return base_ring(1)
    if hasattr(base_ring, "zeta"):
        return base_ring.zeta(order)
    return CyclotomicField(order).gen()


def multiplicative_character_data(F, base_ring=None):
    order = F.order() - 1

    if order == 1:
        return {
            "F": F,
            "order": order,
            "base_ring": base_ring if base_ring is not None else QQ,
            "generator": F(1),
            "log": {F(1): 0},
            "zeta": QQ(1),
        }

    generator = F.multiplicative_generator()
    K = base_ring if base_ring is not None else CyclotomicField(order)

    return {
        "F": F,
        "order": order,
        "base_ring": K,
        "generator": generator,
        "log": {generator**k: k for k in range(order)},
        "zeta": root_of_unity(K, order),
    }


def multiplicative_character_value(a, exponent, character_data):
    F = character_data["F"]
    a = F(a)
    if a == 0:
        raise ValueError("multiplicative characters are only defined on F_q^*")

    order = character_data["order"]
    if order == 1:
        return character_data["base_ring"](1)

    return character_data["zeta"] ** (
        (exponent % order) * character_data["log"][a]
    )


def determinant_character_value(M, exponent, character_data):
    return multiplicative_character_value(matrix_of(M).det(), exponent, character_data)


def scalar_matrix(value, base_ring):
    return Matrix(base_ring, 1, 1, [base_ring(value)])


def coerce_matrix_to_field(M, F):
    M = matrix_of(M)
    return Matrix(F, M.nrows(), M.ncols(), [F(x) for x in M.list()])


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


def left_coset_data_from_subgroup(data, subgroup_matrices):
    all_matrices = [g.matrix() for g in data["G"]]
    remaining = {matrix_key(M) for M in all_matrices}
    representatives = []
    decomposition = {}

    for M in all_matrices:
        key = matrix_key(M)
        if key not in remaining:
            continue

        representatives.append(M)
        coset_index = len(representatives) - 1

        for h in subgroup_matrices:
            X = h * M
            X_key = matrix_key(X)
            decomposition[X_key] = (coset_index, h)
            remaining.discard(X_key)

    assert len(representatives) * len(subgroup_matrices) == data["G"].order()
    assert not remaining

    return {
        "subgroup": subgroup_matrices,
        "representatives": representatives,
        "decomposition": decomposition,
    }


def upper_borel_matrices(data):
    return [
        g.matrix()
        for g in data["G"]
        if g.matrix()[1, 0] == 0
        and g.matrix()[2, 0] == 0
        and g.matrix()[2, 1] == 0
    ]


def upper_block_parabolic_2_1_matrices(data):
    return [
        g.matrix()
        for g in data["G"]
        if g.matrix()[2, 0] == 0 and g.matrix()[2, 1] == 0
    ]


def levi_2_1_blocks(M, F):
    M = matrix_of(M)
    A = Matrix(F, 2, 2, [F(M[i, j]) for i in range(2) for j in range(2)])
    d = F(M[2, 2])
    return A, d


def build_induced_representation(
    data,
    subgroup_matrices,
    coefficient_matrix,
    coefficient_dimension,
    base_ring,
    name,
    coset_data=None,
):
    cosets = (
        coset_data
        if coset_data is not None
        else left_coset_data_from_subgroup(data, subgroup_matrices)
    )
    return {
        "data": data,
        "name": name,
        "base_ring": base_ring,
        "coefficient_matrix": coefficient_matrix,
        "coefficient_dimension": coefficient_dimension,
        "representatives": cosets["representatives"],
        "decomposition": cosets["decomposition"],
        "decompose": cosets.get("decompose"),
        "coefficient_cache": {},
    }


def induced_dimension(induced):
    return len(induced["representatives"]) * induced["coefficient_dimension"]


def cached_coefficient_matrix(h, induced):
    key = matrix_key(h)
    cache = induced["coefficient_cache"]
    if key not in cache:
        cache[key] = induced["coefficient_matrix"](h)
    return cache[key]


def induced_matrix(M, induced):
    M = matrix_of(M)
    representatives = induced["representatives"]
    decomposition = induced["decomposition"]
    d = induced["coefficient_dimension"]
    K = induced["base_ring"]
    A = Matrix(K, len(representatives) * d, len(representatives) * d, 0)

    for i, representative in enumerate(representatives):
        if induced["decompose"] is None:
            j, h = decomposition[matrix_key(representative * M)]
        else:
            j, h = induced["decompose"](representative * M)
        block = cached_coefficient_matrix(h, induced)
        for row in range(d):
            for col in range(d):
                A[i * d + row, j * d + col] = block[row, col]

    return A


def character_rows(data, matrix_function, trace_key="trace"):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                trace_key: matrix_function(representative.matrix()).trace(),
            }
        )

    rows.sort(key=lambda row: (row["order"], row["size"], str(row[trace_key])))
    return rows


def verify_induced_representation(induced, rows=None):
    data = induced["data"]
    G = data["G"]
    g = G.random_element().matrix()
    h = G.random_element().matrix()

    assert induced_matrix(g * h, induced) == induced_matrix(g, induced) * induced_matrix(h, induced)
    assert induced_matrix(G.one().matrix(), induced).trace() == induced_dimension(induced)

    if rows is not None:
        assert hermitian_inner_product(rows, "trace", "trace") == 1


def print_character_table(rows, label):
    print(f"{label} character by conjugacy class:")
    for row in rows:
        print(
            f"  order {row['order']:3}, size {row['size']:7}: "
            f"trace = {row['trace']}"
        )


def line_from_vector(V, v):
    return V.subspace([V(v)])


def plane_from_vectors(V, v, w):
    return V.subspace([V(v), V(w)])


def vector_not_in_subspace(V, subspace):
    for v in V:
        if v != 0 and not line_from_vector(V, v).is_subspace(subspace):
            return V(v)
    raise ValueError("could not find vector outside subspace")


def vector_in_plane_not_in_line(V, plane, line):
    for v in plane:
        if v != 0 and not line_from_vector(V, v).is_subspace(line):
            return V(v)
    raise ValueError("could not find vector in plane outside line")


def upper_triangular_3x3(M):
    return M[1, 0] == 0 and M[2, 0] == 0 and M[2, 1] == 0


def upper_block_2_1(M):
    return M[2, 0] == 0 and M[2, 1] == 0


def borel_flag_coset_data(data):
    """Represent B\\G by complete flags in the row space.

    Left multiplication by the upper triangular Borel preserves the flag
    span(row_3) <= span(row_2, row_3).
    """
    F = data["F"]
    V = VectorSpace(F, 3)
    representatives = []
    flag_index = {}

    for line in V.subspaces(1):
        r3 = V(line.basis()[0])
        for plane in V.subspaces(2):
            if not line.is_subspace(plane):
                continue

            r2 = vector_in_plane_not_in_line(V, plane, line)
            r1 = vector_not_in_subspace(V, plane)
            R = Matrix(F, [r1, r2, r3])
            assert R.det() != 0

            flag_index[(line, plane)] = len(representatives)
            representatives.append(R)

    def decompose(X):
        X = matrix_of(X)
        row3 = V(X.row(2))
        row2 = V(X.row(1))
        line = line_from_vector(V, row3)
        plane = plane_from_vectors(V, row2, row3)
        j = flag_index[(line, plane)]
        h = X * representatives[j].inverse()
        assert upper_triangular_3x3(h)
        return j, h

    return {
        "subgroup": None,
        "representatives": representatives,
        "decomposition": None,
        "decompose": decompose,
    }


def parabolic_2_1_point_coset_data(data):
    """Represent P_{2,1}\\G by lines in the row space.

    Left multiplication by P_{2,1} preserves span(row_3).
    """
    F = data["F"]
    V = VectorSpace(F, 3)
    representatives = []
    line_index = {}

    for line in V.subspaces(1):
        r3 = V(line.basis()[0])
        r1 = vector_not_in_subspace(V, line)
        plane = V.subspace([r1, r3])
        r2 = vector_not_in_subspace(V, plane)
        R = Matrix(F, [r1, r2, r3])
        assert R.det() != 0

        line_index[line] = len(representatives)
        representatives.append(R)

    def decompose(X):
        X = matrix_of(X)
        row3 = V(X.row(2))
        line = line_from_vector(V, row3)
        j = line_index[line]
        h = X * representatives[j].inverse()
        assert upper_block_2_1(h)
        return j, h

    return {
        "subgroup": None,
        "representatives": representatives,
        "decomposition": None,
        "decompose": decompose,
    }
