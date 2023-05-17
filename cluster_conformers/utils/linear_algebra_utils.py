"""
Functions for performing linear algebra on matrices pertaining to protein CA distances.

Used by other scripts in the peptide_analysis library and the ClusterConformers() class
in the find_conformers() script. These functions perform the majority of the
calculations performed in the package, with the scripts they call handling specific
logic, such as deciding which chains to call these functions on.
"""

from math import dist

from numpy import (
    absolute,
    asarray,
    delete,
    intersect1d,
    isnan,
    ndarray,
    s_,
    sum,
    triu_indices,
)


def euclidean(
    coords_3D_1: "tuple[float, float, float]", coords_3D_2: "tuple[float, float, float]"
) -> float:
    """
    Finds the Euclidean distance between two points in 3D space. Now needed as math
    module in Python3.x <= 3.7 lacks the math.dist function. Alternatives from Scipy and
    Numpy are incredibly slow in comparison.

    :param coords_3D_1: Point 1 = (x, y, z)
    :type coords_3D_1: tuple[float, float, float]
    :param coords_3D_2: Point 2 = (x, y, z)
    :type coords_3D_2: tuple[float, float, float]
    :return: Euclidean distance between two Cartesian coordinates
    :rtype: float
    """
    return dist(coords_3D_1, coords_3D_2)


def generate_ca_matx(x: "list[float]", y: "list[float]", z: "list[float]") -> ndarray:
    """
    Takes three lists corresponding to x,y,z coordinates and returns a matrix of the
    pair-wise Euclidian distances between all coordinates. Each input must be
    1*N-dimensional and an N*N-dimensional Numpy matrix is returned.
    Requirement: `len(x) == len(y) == len(z)`.

    :param x: Cartesian x coords
    :type x: list[float]
    :param y: Cartesian y coords
    :type y: list[float]
    :param z: Cartesian z coords
    :type z: list[float]
    :return: List of Euclidean distances
    :rtype: list[float]
    """

    # Populated and returned
    ca_dist_matx = []

    # Row loop
    for x_a, y_a, z_a in zip(x, y, z):
        distance_row = []
        res1_coords = (x_a, y_a, z_a)

        # Column loop
        for x_b, y_b, z_b in zip(x, y, z):
            res2_coords = (x_b, y_b, z_b)

            # Euclidian distance
            dist = euclidean(res1_coords, res2_coords)

            # Add element to row
            distance_row.append(dist)

        # Add row to matrix
        ca_dist_matx.append(distance_row)

    # For faster linear algebra later
    return asarray(ca_dist_matx)


def matx_subtract(matx1: ndarray, matx2: ndarray) -> ndarray:
    """
    Returns difference between matx1 and matx2. Sign of the elements in the returned
    difference matrix is ignored -- all elements in matrix > 0.

    :param matx1: Matrix 1
    :type matx1: np.ndarray
    :param matx2: Matrix 2
    :type matx2: np.ndarray
    :return: Matrix of absolute differences
    :rtype: np.ndarray
    """

    return absolute(matx1 - matx2)


def matx_trim(matx: ndarray, end: int) -> ndarray:
    """
    Truncates the final (end) rows and columns for a given matrix.
    end = the new end of the matrix. If the number of rows/columns in the matrix exceeds
    `end`, the function will remove them, returning the truncated matrix.

    :param matx: Square matrix
    :type matx: np.ndarray
    :param end: Final index position before truncation
    :type end: int
    :raises IndexError: Negative index position parsed
    :return: Trimmed (end truncated) matrix
    :rtype: np.ndarray
    """

    # Only trim if end index is positive
    if end > 0:
        for axis in (0, 1):
            matx = delete(matx, s_[end:], axis)
        return matx
    else:
        raise IndexError("Index to trim end matrix rows and columns must be positive.")


def upper_triangle(matx: ndarray, res_mask: int = 0, numeric: bool = True) -> ndarray:
    """
    Masks the diagonal of a given matrix within a given window size. Returns a 1D list
    of the linearised elements.

    :param matx: Square matrix
    :type matx: np.ndarray
    :param res_mask: Position above diagonal to begin transformation, defaults to 0
    :type res_mask: int, optional
    :param numeric: Remove np.NaN values, defaults to True
    :type numeric: bool, optional
    :return: Upper triangle as a 1D array
    :rtype: np.array
    """

    # Apply residue mask
    out_vector = matx[triu_indices(matx.shape[0], k=res_mask)]

    # Remove NaNs
    if numeric:
        out_vector = out_vector[~isnan(out_vector)]

    return out_vector


def trim_to_smallest(mtx1: ndarray, mtx2: ndarray) -> "tuple[ndarray, ndarray]":
    """
    Given two matrices, function truncates the larger one by removing exactly the number
    of rows and columns (from its end) by which it is larger. Either mtx1 or mtx2 can be
    the smaller of the two.

    :param mtx1: Square matrix 1
    :type mtx1: np.ndarray
    :param mtx2: Square matrix 2
    :type mtx2: np.ndarray
    :return: Matrix 1 and 2, truncated to whichever is smallest
    :rtype: tuple( np.ndarray, np.ndarray )
    """

    # Identify smaller matrix
    len_mtx1 = mtx1.shape[0]
    len_mtx2 = mtx2.shape[0]
    smallest = min(len_mtx1, len_mtx2)

    # Truncate the larger matrix.
    if len_mtx1 > smallest:
        mtx1 = matx_trim(mtx1, smallest)
    else:
        mtx2 = matx_trim(mtx2, smallest)

    return mtx1, mtx2


def find_max(matx: ndarray, starting_max: float = 0.0) -> float:
    """
    Returns the maximum element from a matrix of floats, if any are larger than the
    starting_max value. Otherwise, starting_max is returned.

    :param matx: N*M matrix
    :type matx: np.ndarray
    :param starting_max: Initial maximum value to beat, defaults to 0.0
    :type starting_max: float, optional
    :return: New maximum value (if >starting_max)
    :rtype: float
    """

    # Remove NaN values
    return max(matx[~isnan(matx)].max(), starting_max)


def len_overlap(unps1: list, unps2: list) -> int:
    """
    Parse in a distance difference matrix and get the length of the overlap.
    NB: Replace with function that takes two UNP lists and returns the Nump 1D
    intersect.

    :param unps1: UniProt residue IDs for protein 1
    :type unps1: list
    :param unps2: UniProt residue IDs for protein 2
    :type unps2: list
    :return: Length of residue overlap
    :rtype: int
    """

    overlap = intersect1d(asarray(unps1), asarray(unps2)).shape[0]

    return overlap


def calc_score(dd_matx: ndarray, unps1: set, unps2: set, res_mask: int = 1) -> float:
    """
    Calcultes the distance-based score used for agglomerative clustering.

    :param dd_matx: Square matrix of absolute distance differences
    :type dd_matx: np.ndarray
    :param unps1: UniProt residue IDs for protein 1
    :type unps1: set
    :param unps2: UniProt residue IDs for protein 2
    :type unps2: set
    :param res_mask: Number of elements above diagonal to ignore, defaults to 1
    :type res_mask: int, optional
    :return: CA distance-based clustering score
    :rtype: float
    """

    # Convert upper triangle of distance difference matrix to 1D vector
    upper_tri_vector = upper_triangle(dd_matx, res_mask=res_mask)

    upper_triangle_sum = sum(upper_tri_vector)

    # Length of overlap
    pdbe1_pdbe2_ovlp = len(set(unps1).intersection(set(unps2)))

    # Calculate score
    score = (pdbe1_pdbe2_ovlp**2) / (len(unps1) * len(unps2)) * upper_triangle_sum

    return score
