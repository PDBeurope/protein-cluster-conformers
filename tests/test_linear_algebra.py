"""
TODO:
 - consider refactoring some of the functions to minimise the number of asserts per
    method
 - Add decorators
 - Simplify some of the test objects
"""

# import os
# import pathlib

# For testing the test script
# import sys
import unittest
from unittest import TestCase, mock

import numpy as np

# Import functions to test
from cluster_conformers.utils import linear_algebra_utils

# # Modify Python path for custom imports
# sys.path[0] = str(  # Override default Python import part
#     pathlib.Path(
#         os.path.realpath(__file__)  # Path to this file
#     ).parent.parent  # Path to folder this file is in
# )


# # Global variables
# PATH_HOME = pathlib.Path.home()
# PATH_HOME = PATH_HOME.joinpath(
#     "EMBL-EBI", "funclan_work", "contact_map_difference", "tests"
# )
# PATH_BASE = str(PATH_HOME)

# Mock matrices
MOCK_MATX_1 = np.array(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]  # Positive floats
)

MOCK_MATX_2 = np.array(
    [[-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0], [-7.0, -8.0, -9.0]]  # Negative floats
)

MOCK_MATX_3 = np.array(
    [  # Positive floats (for testing returned negatives)
        [9.0, 8.0, 7.0],
        [6.0, 5.0, 4.0],
        [3.0, 2.0, 1.0],
    ]
)

MOCK_MATX_4 = np.array(
    [[np.NaN, 2.0, 3.0], [4.0, np.NaN, 6.0], [7.0, 8.0, np.NaN]]  # Inclusion of np.NaN
)

MOCK_MATX_5 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])  # Positive floats

MOCK_MATX_DENARY_REPEATED = np.array(
    [  # Integer matrix
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 2, 3, 4, 5, 6, 7, 8, 9],
    ]
)

MOCK_MATX_DENARY_EVEN_REPEATED = np.array(
    [  # Even integer matrix
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
    ]
)


class TestLinearAlgebraUtils(TestCase):
    """
    Runs tests on functions included in the cluster_conformers.utils.io_utils script.
    """

    def test_euclidean(self) -> None:
        """
        Tests whether the Euclidean distance function for 3D coordinates returns
        expected results.
        """

        # Check distance from a point to same point is zero
        test_euclidean_distance = linear_algebra_utils.euclidean((1, 2, 3), (1, 2, 3))
        self.assertEqual(
            test_euclidean_distance,
            0.0,
            msg="Distance between the same two 3D points should be zero, not "
            f"{test_euclidean_distance}",
        )

        # Check floats are handled
        test_euclidean_distance = linear_algebra_utils.euclidean(
            (2.5, 3.5, 4.5), (1.1, 2.2, 3.3)
        )
        self.assertAlmostEqual(
            test_euclidean_distance,
            2.2561028345357,
            msg="Positive floats should be handled in the input tuples.",
        )

        # Check negative numbers are handled
        test_euclidean_distance = linear_algebra_utils.euclidean(
            (-2.5, -3.5, -4.5), (1.1, 2.2, 3.3)
        )
        self.assertAlmostEqual(
            test_euclidean_distance,
            10.309704166464,
            msg="Negative floats should be handled. Negative floats in the first "
            "argument causes miscalculation of distance. ",
        )

        # Check negative numbers are handled
        test_euclidean_distance = linear_algebra_utils.euclidean(
            (2.5, 3.5, 4.5), (-1.1, -2.2, -3.3)
        )
        self.assertAlmostEqual(
            test_euclidean_distance,
            10.309704166464,
            msg="Negative floats should be handled. Negative floats in the second "
            "argument causes miscalculation of distance. ",
        )

        # Check negative numbers handled
        test_euclidean_distance = linear_algebra_utils.euclidean(
            (-2.5, -3.5, -4.5), (-1.1, -2.2, -3.3)
        )
        self.assertAlmostEqual(
            test_euclidean_distance,
            2.2561028345357,
            msg="Negative floats should be handled. Negative floats in the first and "
            "second argument causes miscalculation of distance. ",
        )

        # Check NaNs are handled
        test_euclidean_distance = linear_algebra_utils.euclidean(
            (np.NaN, 3.5, 4.5), (np.NaN, 2.2, 3.3)
        )
        self.assertTrue(
            np.isnan(test_euclidean_distance),
            msg="numpy.NaN values should be handled and return numpy.NaN if one "
            "element is numpy.NaN. The function is currently returning "
            f"{test_euclidean_distance}",
        )

    @mock.patch(
        "cluster_conformers.utils.linear_algebra_utils.euclidean",
        side_effect=[
            np.NaN,
            np.NaN,
            np.NaN,
            np.NaN,
            np.NaN,
            np.NaN,
            np.NaN,
            0.0,
            7.1414,
            5018.0,
            7.4833,
            29.172,
            np.NaN,
            7.1414,
            0.0,
            5011.0,
            8.1854,
            22.181,
            np.NaN,
            5018.0,
            5011.0,
            0.0,
            5014.4,
            4989.1,
            np.NaN,
            7.4833,
            8.1854,
            5014.4,
            0.0,
            25.981,
            np.NaN,
            29.172,
            22.181,
            4989.1,
            25.981,
            0.0,
        ],
    )
    def test_generate_ca_matx(self, mock_euclidean):
        """
        Test for the function that generates a distance matrix from 3 lists of
        equal-length (length=N) atomic coordinates. The 3 lists should contain the
        Cartesian x,y,z coordinates of (CA) atoms. Returned should be a N*N-dimensional
        matrix containing the pairwise Euclidean distances between each parsed atom.
        """
        # Coordinate inputs
        x_coords_expected = [np.NaN, 1, 2.0, 300.0, -5.0, 0]
        y_coords_expected = [np.NaN, 3.0, 2, -30.0, 5.0, 0]
        z_coords_expected = [np.NaN, -9.0, -2.0, 5000.0, -5.0, 20]

        # Run the function
        test_matx = linear_algebra_utils.generate_ca_matx(
            x_coords_expected, y_coords_expected, z_coords_expected
        )

        # Test the output
        self.assertTrue(
            np.array_equal(
                test_matx,
                np.array(
                    [
                        [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN],
                        [np.NaN, 0.0, 7.1414, 5018.0, 7.4833, 29.172],
                        [np.NaN, 7.1414, 0.0, 5011.0, 8.1854, 22.181],
                        [np.NaN, 5018.0, 5011.0, 0.0, 5014.4, 4989.1],
                        [np.NaN, 7.4833, 8.1854, 5014.4, 0.0, 25.981],
                        [np.NaN, 29.172, 22.181, 4989.1, 25.981, 0.0],
                    ]
                ),
                equal_nan=True,
            ),
            msg="Function is not constructing the correct array from input lists. ",
        )

    def test_matx_subtract(self):
        """
        Tests the matx_subtract() function. The tested function should take two
        N*N-dimensional matrices.
        """

        # Check negative numbers are handled
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_subtract(MOCK_MATX_1, MOCK_MATX_2),
                np.array(
                    [  # Expected output
                        [2.0, 4.0, 6.0],
                        [8.0, 10.0, 12.0],
                        [14.0, 16.0, 18.0],
                    ]
                ),
            ),
            msg="The matx_subtract() function cannot handle negative values in given "
            "matrices properly.",
        )

        # Check subtraction of identical matrices returns empty matrix
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_subtract(MOCK_MATX_1, MOCK_MATX_1),
                np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),  # Expected output
            ),
            msg="The matx_subtract() function should return empty (filled with zeros) "
            "matrix when parsing identical matrices.",
        )

        # Check negative elements are not returned when b_ij > a_ij, for a_ij - b_ij
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_subtract(MOCK_MATX_1, MOCK_MATX_3),
                np.array(
                    [  # Expected output
                        [8.0, 6.0, 4.0],
                        [2.0, 0.0, 2.0],
                        [4.0, 6.0, 8.0],
                    ]
                ),
            ),
            msg="The matx_subtract() function should return a matrix of only positive "
            "elements, even when b_ij > a_ij, for a_ij - b_ij.",
        )

        # Check np.NaN values are handled correctly
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_subtract(MOCK_MATX_1, MOCK_MATX_4),
                np.array(
                    [  # Expectded output
                        [np.NaN, 0.0, 0.0],
                        [0.0, np.NaN, 0.0],
                        [0.0, 0.0, np.NaN],
                    ]
                ),
                equal_nan=True,
            ),
            msg="The matx_subtract() function should handle np.NaN values. Subtraction "
            "of np.NaN from float (or float from np.NaN) should return np.NaN for "
            "element position.",
        )

    def test_matx_trim(self):
        """
        Tests the function which removes a specified number of rows and columns from
        the end of the parsed matrix. Should not be allowed to parse negative trim
        ranges.
        """
        # Remove single row/column
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_trim(MOCK_MATX_DENARY_REPEATED, end=8),
                np.array(
                    [
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                        [1, 2, 3, 4, 5, 6, 7, 8],
                    ]
                ),
            ),
            msg="Array is not being trimmed from the end row-column properly. In this "
            "case, a single row-column pair from the end should be removed.",
        )

        # Remove multiple rows/columns
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.matx_trim(MOCK_MATX_DENARY_REPEATED, end=4),
                np.array([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]),
            ),
            msg="Array is not being trimmed properly from the end rows-columns "
            "properly. In this case, multiple rows/columns should have been removed.",
        )

        # Should raise IndexError. Prevents slicing from matx start
        return_msg = (
            "Parsing in a negative index position will attempt to slice "
            + "from the beginning of the matrix. This is not the desired "
            + "and so an IndexError must be returned"
        )
        with self.assertRaises(IndexError, msg=return_msg):
            linear_algebra_utils.matx_trim(MOCK_MATX_1, -1)

    def test_upper_triangle(self):
        """
        Function should convert a matrix into a vector. The `numeric` argument should be
        set to False when np.NaN values should be included.
        """
        # Linearise with no mask
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(MOCK_MATX_1, numeric=True),
                np.array([1.0, 2.0, 3.0, 5.0, 6.0, 9.0]),
            ),
            msg="Function is not correcty linearising a simple matrix. ",
        )

        # Linearise and apply mask
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(
                    MOCK_MATX_1, res_mask=1, numeric=True
                ),
                np.array([2.0, 3.0, 6.0]),
            ),
            msg="Function is not applying the correct residue mask (1 in this case) "
            "to return the upper triangle as a 1*N Numpy array.",
        )

        # Linearise and apply large mask
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(
                    MOCK_MATX_DENARY_REPEATED, res_mask=4, numeric=True
                ),
                np.array(
                    [
                        5.0,
                        6.0,
                        7.0,
                        8.0,
                        9.0,
                        6.0,
                        7.0,
                        8.0,
                        9.0,
                        7.0,
                        8.0,
                        9.0,
                        8.0,
                        9.0,
                        9.0,
                    ]
                ),
            ),
            msg="Function is not applying the correct residue mask (1 in this case) "
            "to return the upper triangle as a 1*N Numpy array.",
        )
        # Linearise and remove np.NaNs
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(
                    MOCK_MATX_4, res_mask=0, numeric=True
                ),
                np.array([2.0, 3.0, 6.0]),
            ),
            msg="Function does not remove np.NaN elements from matrix.",
        )

        # Linearise and keep np.NaNs
        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(
                    MOCK_MATX_4, res_mask=0, numeric=False
                ),
                np.array([np.NaN, 2.0, 3.0, np.NaN, 6.0, np.NaN]),
                equal_nan=True,
            ),
            msg="Function removes np.NaN values even when forced not to (using "
            "numeric=False).",
        )

        # Linearise a string-containing matrix
        string_based_matrix = np.array(
            [
                ["label1", "label2", "label3"],
                ["label4", "label5", "label6"],
                ["label7", "label8", "label9"],
            ]
        )

        self.assertTrue(
            np.array_equal(
                linear_algebra_utils.upper_triangle(string_based_matrix, numeric=False),
                np.array(["label1", "label2", "label3", "label5", "label6", "label9"]),
            ),
            msg="Function removes np.NaN values even when forced not to (using "
            "numeric=False).",
        )

    @mock.patch(
        "cluster_conformers.utils.linear_algebra_utils.matx_trim",
        side_effect=[
            np.array(
                [
                    # Output should match the dimension of the smaller matx. 9*9 here
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    [2, 4, 6, 8, 10, 12, 14, 16, 18],
                ]
            )
        ]
        * 2,
    )
    def test_trim_to_smallest(self, mock_matx_trim):
        """
        Tests the function which trims the largest of two matrices from the C-terminal
        end of the CA matrix (the end of the Numpy array in this case).
        """
        # Run trim function
        not_trimmed_matx, trimmed_matx = linear_algebra_utils.trim_to_smallest(
            MOCK_MATX_DENARY_REPEATED,  # Smaller
            MOCK_MATX_DENARY_EVEN_REPEATED,  # Larger
        )

        # Check the smaller matrix was not changed
        self.assertTrue(
            np.array_equal(not_trimmed_matx, MOCK_MATX_DENARY_REPEATED),
            msg="The smaller matrix was erroneously trimmed.",
        )

        # Check the larger matrix was trimmed from the end correctly
        self.assertTrue(
            np.array_equal(
                trimmed_matx,
                np.array(
                    [  # Output should match the dimension of the smaller matx. 9*9 here
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    ]
                ),
            ),
            msg="Larger matrix was not trimmed, or not trimmed correctly. ",
        )

        # Check whether the function can trim the larger matrix, when parsed in first
        # Run trim function
        trimmed_matx, not_trimmed_matx = linear_algebra_utils.trim_to_smallest(
            MOCK_MATX_DENARY_EVEN_REPEATED,  # Larger
            MOCK_MATX_DENARY_REPEATED,  # Smaller
        )

        # Check the smaller matrix was not changed
        self.assertTrue(
            np.array_equal(not_trimmed_matx, MOCK_MATX_DENARY_REPEATED),
            msg="The smaller matrix was erroneously trimmed.",
        )

        # Check the larger matrix was trimmed from the end correctly
        self.assertTrue(
            np.array_equal(
                trimmed_matx,
                np.array(
                    [  # Output should match the dimension of the smaller matx. 9*9 here
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                        [2, 4, 6, 8, 10, 12, 14, 16, 18],
                    ]
                ),
            ),
            msg="Larger matrix was not trimmed, or not trimmed correctly. ",
        )

    def test_find_max(self):
        """
        Function to test the find_max() function, which should return the maximum
        value inside a given matrix. An optional argument can be parsed to set an
        existing max value; if this parsed maximum is exceeded by an element in the
        given matrix, it is updated and returned from the function.
        """
        # Simple check
        self.assertEqual(
            linear_algebra_utils.find_max(MOCK_MATX_1),
            9,
            msg="Function cannot find the largest value in a simple matrix.",
        )

        # Replace smaller starting_max with 9
        self.assertEqual(
            linear_algebra_utils.find_max(MOCK_MATX_1, starting_max=2),
            9,
            msg="Function cannot find the maxmimum value and replace it with the "
            "smaller initial max value (in this case 2) parsed into the funcion.",
        )

        # starting_max is greater than any element in matrix so should not replace 100
        self.assertEqual(
            linear_algebra_utils.find_max(MOCK_MATX_1, starting_max=100),
            100,
            msg="The parsed in starting_max is larger than any element in the parsed "
            "matrix, although the function is still returning the highest value in "
            "the matrix by mistake.",
        )

        # Check function can handle np.NaN values
        self.assertEqual(
            linear_algebra_utils.find_max(MOCK_MATX_4),
            8,
            msg="Function cannot handle np.NaN values contained in parsed matrix. ",
        )

    def test_len_overlap(self):
        """
        Tests the function which calculates the length of the overlap between two parsed
        lists of UniProt residue indices. The lists do not necessarily need to be
        ordered but are often loaded in ascending order by default.
        """

        unps1_test = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        unps2_test = [3, 4, 5, 8]
        unps3_test = [10, 11, 12, 13]
        unps4_test = [4, 8, 15, 20]

        # Check unps2_test <intersect> unps1_test = len(unps2_test)
        self.assertEqual(
            linear_algebra_utils.len_overlap(unps1=unps1_test, unps2=unps2_test),
            4,
            msg="Cannot find the intersect between elements: larger-smaller",
        )

        # Check unps2_test <intersect> unps1_test = len(unps2_test). Reversed
        self.assertEqual(
            linear_algebra_utils.len_overlap(unps1=unps2_test, unps2=unps1_test),
            4,
            msg="Cannot find the intersect between elements: smaller-larger",
        )

        # Check unps1_test <intersect> unps3_test = 0
        self.assertEqual(
            linear_algebra_utils.len_overlap(unps1=unps1_test, unps2=unps3_test),
            0,
            msg="Erroneous behaviour: Elements in list A do not appear in list B but "
            "the length of the returned intersect is non-zero...",
        )

        # Check unps2_test <intersect> unps4_test = 2
        self.assertEqual(
            linear_algebra_utils.len_overlap(unps1=unps2_test, unps2=unps4_test),
            2,
            msg="These two lists contain two shared elements and the rest of their "
            "contents are different. Function is not identifying these two elements",
        )

    @mock.patch(
        "cluster_conformers.utils.linear_algebra_utils.upper_triangle",
        side_effect=[
            np.array([2.0, 3.0, 6.0]),
            np.array([2.0, 3.0, 6.0]),
            np.array(
                [2.0, 6.0]  # np.NaN in the middle gets removed by upper_triangle() func
            ),
            np.array([1.0, 2.0, 3.0, 5.0, 6.0, 9.0]),
            np.array(
                [
                    8,
                    10,
                    12,
                    14,
                    16,
                    18,
                    20,
                    22,
                    10,
                    12,
                    14,
                    16,
                    18,
                    20,
                    22,
                    12,
                    14,
                    16,
                    18,
                    20,
                    22,
                    14,
                    16,
                    18,
                    20,
                    22,
                    16,
                    18,
                    20,
                    22,
                    18,
                    20,
                    22,
                    20,
                    22,
                    22,
                ]
            ),
        ],
    )
    @mock.patch(
        "cluster_conformers.utils.linear_algebra_utils.len_overlap",
        side_effect=[3, 2, 3, 3, 3],
    )
    def test_calc_score(self, mock_upper_triangle, mock_len_overlap):
        """
        Tests the score function, eventually used to cluster protein conformations.
        """

        # Default function
        self.assertEqual(
            linear_algebra_utils.calc_score(
                dd_matx=MOCK_MATX_1, unps1={1, 2, 3}, unps2={1, 2, 3}
            ),
            11,
            msg="Returns incorrect score when UniProt overlap is perfect.",
        )

        # Partial overlap
        self.assertEqual(
            linear_algebra_utils.calc_score(
                dd_matx=MOCK_MATX_1, unps1={1, 2}, unps2={1, 2, 3}
            ),
            11 * 2 / 3,
            msg="Returns incorrect score when UniProt overlap is partial.",
        )

        # Handling np.NaN values
        self.assertEqual(
            linear_algebra_utils.calc_score(
                dd_matx=np.array([[1, 2, np.NaN], [4, 5, 6], [7, 8, 9]]),
                unps1={1, 2, 3},
                unps2={1, 2, 3},
            ),
            8,
            msg="Function cannot handle np.NaN in upper triangle of matrix.",
        )

        # Removing mask
        self.assertEqual(
            linear_algebra_utils.calc_score(
                dd_matx=MOCK_MATX_1, unps1={1, 2, 3}, unps2={1, 2, 3}, res_mask=0
            ),
            26,
            msg="Removing the residue mask returns incorrect score.",
        )

        # Increasing mask (larger test matrix used)
        self.assertEqual(
            linear_algebra_utils.calc_score(
                dd_matx=MOCK_MATX_DENARY_EVEN_REPEATED,
                unps1={1, 2, 3},
                unps2={1, 2, 3},
                res_mask=3,
            ),
            624,
            msg="Increasing the residue mask returns incorrect score.",
        )


# Run unit tests on call of script
if __name__ == "__main__":

    unittest.main()
