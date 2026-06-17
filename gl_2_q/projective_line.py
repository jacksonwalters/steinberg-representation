from sage.all import *


def build_projective_line_data(q):
    if q < 2:
        raise ValueError("q must be a prime power greater than 1")

    F = GF(q, name="a")
    V = VectorSpace(F, 2)
    G = GL(2, F)
    points = list(V.subspaces(1))
    point_index = {point: i for i, point in enumerate(points)}

    sum_map = Matrix(QQ, 1, len(points), [1] * len(points))
    steinberg = sum_map.right_kernel()

    return {
        "q": q,
        "F": F,
        "V": V,
        "G": G,
        "points": points,
        "point_index": point_index,
        "sum_map": sum_map,
        "steinberg": steinberg,
    }


def matrix_of(g):
    return g.matrix() if hasattr(g, "matrix") else g


def normalized_vector_tuple(v):
    coords = list(v)
    for x in coords:
        if x != 0:
            inv = x**(-1)
            return tuple(inv * y for y in coords)
    raise ValueError("zero vector has no projective normalization")


def point_repr(point):
    return normalized_vector_tuple(point.basis()[0])


def act_on_subspace(M, S):
    M = matrix_of(M)
    W = S.ambient_vector_space()
    basis = S.basis_matrix().transpose()
    vecs = [W(M * basis.column(i)) for i in range(basis.ncols())]
    return W.subspace(vecs)


def projective_line_images(M, data):
    M = matrix_of(M)
    return [
        data["point_index"][act_on_subspace(M, point)]
        for point in data["points"]
    ]


def permutation_matrix_from_images(images, base_ring=QQ):
    n = len(images)
    P = Matrix(base_ring, n, n, 0)

    for i, j in enumerate(images):
        P[j, i] = base_ring(1)

    return P


def projective_line_permutation_matrix(M, data, base_ring=QQ):
    return permutation_matrix_from_images(
        projective_line_images(M, data),
        base_ring=base_ring,
    )


def fixed_projective_points(M, data):
    return sum(1 for i, j in enumerate(projective_line_images(M, data)) if i == j)


def augmentation_matrix_from_permutation(permutation_matrix):
    n = permutation_matrix.nrows()
    base_ring = permutation_matrix.base_ring()
    columns = []

    for i in range(n - 1):
        v = vector(base_ring, n)
        v[i] = 1
        v[n - 1] = -1

        image = permutation_matrix * v
        assert sum(image) == 0
        columns.append(list(image[: n - 1]))

    return matrix(base_ring, columns).transpose()


def steinberg_matrix(M, data, base_ring=QQ):
    permutation = projective_line_permutation_matrix(M, data, base_ring=base_ring)
    return augmentation_matrix_from_permutation(permutation)


def steinberg_trace(M, data):
    return fixed_projective_points(M, data) - 1


def character_rows(data):
    rows = []

    for conjugacy_class in data["G"].conjugacy_classes():
        representative = conjugacy_class.representative()
        fixed_points = fixed_projective_points(representative.matrix(), data)

        rows.append(
            {
                "order": representative.order(),
                "size": len(conjugacy_class),
                "fixed_points": fixed_points,
                "permutation": fixed_points,
                "trivial": 1,
                "steinberg": fixed_points - 1,
            }
        )

    rows.sort(
        key=lambda row: (
            row["order"],
            row["size"],
            row["fixed_points"],
        )
    )
    return rows


def inner_product(rows, left_key, right_key):
    group_order = sum(ZZ(row["size"]) for row in rows)
    return (
        sum(
            ZZ(row["size"]) * ZZ(row[left_key]) * ZZ(row[right_key])
            for row in rows
        )
        / group_order
    )


def verify_action_property(data, representation_matrix):
    G = data["G"]
    g = G.random_element().matrix()
    h = G.random_element().matrix()
    return representation_matrix(g * h, data) == (
        representation_matrix(g, data) * representation_matrix(h, data)
    )


def verify_projective_line_data(data, rows=None):
    q = data["q"]
    points = data["points"]
    sum_map = data["sum_map"]
    steinberg = data["steinberg"]

    assert len(points) == q + 1
    assert sum_map.rank() == 1
    assert steinberg.dimension() == q
    assert verify_action_property(data, projective_line_permutation_matrix)
    assert verify_action_property(data, steinberg_matrix)

    g = data["G"].random_element().matrix()
    permutation = projective_line_permutation_matrix(g, data)
    assert sum_map * permutation == sum_map

    if rows is not None:
        assert inner_product(rows, "trivial", "steinberg") == 0
        assert inner_product(rows, "steinberg", "steinberg") == 1
        assert inner_product(rows, "permutation", "permutation") == 2
