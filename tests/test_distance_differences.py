# Import of modified TestCase class
import os
import pathlib

# For testing the test script
# import sys
import unittest
from unittest import mock

import numpy as np

# Import functions to test
from cluster_conformers import distance_differences

from .test_case import TestCaseModified

# # Modify Python path for custom imports
# sys.path[0] = str(  # Override default Python import part
#     pathlib.Path(
#         os.path.realpath(__file__)  # Path to this file
#     ).parent.parent  # Path to folder this file is in
# )


# Paths
# PATH_HOME = pathlib.Path.home()
# PATH_BASE = PATH_HOME.joinpath(
#     "EMBL-EBI", "funclan_work", "contact_map_difference", "tests"
# )
PATH_BASE = pathlib.Path("./tests")
PATH_SERIALISED_MATXS = PATH_BASE.joinpath("mock_serial_matxs", "small_full")

PATH_BASE = str(PATH_BASE)
PATH_MOCK_MMCIFS = PATH_BASE + "/mock_data/test_mmcifs/"
PATH_SAVE_CA = pathlib.Path(f"{PATH_BASE}/test_output/ca_distances")

# Mock input data
TEST_UNP = "A12345"


class TestDistanceDifferences(TestCaseModified):
    """
    Tests for the script that handles calculations for finding distance differences
    """

    # Decorators needed
    @mock.patch(
        "cluster_conformers.distance_differences.linear_algebra_utils.trim_to_smallest"
    )
    @mock.patch(
        "cluster_conformers.distance_differences.linear_algebra_utils.matx_subtract"
    )
    def test_generate_matx_diff(self, mock_matx_subtract, mock_trim_to_smallest):
        """
        Test for function used to find the absolute difference between two matrices,
        trim the larger of the two to match the dimensions of the smaller input, flatten
        the results elements to zero below a cutoff and save the output.
        """
        # Mocks
        mock_ca_matx1 = np.full(
            (9, 9),
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
        )
        mock_ca_matx2 = np.array(
            [
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
                [10, 20, np.nan, 40, 50, 60, 70, 80],
            ]
        )
        expected_dd_matx = np.array(
            [
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
                [9.0, 18.0, np.nan, 36.0, 45.0, 54.0, 63.0, 72.0],
            ]
        )

        # Patch return objects
        mock_trim_to_smallest.return_value = (
            np.full((8, 8), [1, 2, 3, 4, 5, 6, 7, 8]),
            mock_ca_matx2,
        )

        # Path absolute difference matrix
        mock_matx_subtract.return_value = expected_dd_matx

        # Remove saved files from previous tests
        existing_files = [f"{PATH_SAVE_CA}/{TEST_UNP}_1pdb_A_to_2pdb_B_difference.npz"]

        for file in existing_files:
            os.system(f"rm {file}")

        """
        ##### CHECK 1 #####
        Check function can return distance difference matrix to variable only"""
        # Run function
        dd_matx = distance_differences.generate_matx_diff(
            mock_ca_matx1,
            mock_ca_matx2,
            # return_matx=True
        )
        # Check returned matrix is exactly as expected
        self.assertTrue(
            np.array_equal(dd_matx, expected_dd_matx, equal_nan=True),
            msg="Function does not return the correct distance difference matrix.",
        )

        """
        ##### CHECK 2 #####
        Check function can save the correct distance difference matrix only"""
        distance_differences.generate_matx_diff(
            mock_ca_matx1,
            mock_ca_matx2,
            # PATH_SAVE_CA.joinpath(f"{TEST_UNP}_1pdb_A_to_2pdb_B_difference"),
            # return_matx=False,
        )

        # Check a file has been saved with the expected name. No need to check contents
        # as the test for this function is already performed.
        # self.assertIsFile(f"{PATH_SAVE_CA}/{TEST_UNP}_1pdb_A_to_2pdb_B_difference.npz")

        # Remove saved matrix for next test
        os.system(f"rm {PATH_SAVE_CA}/{TEST_UNP}_1pdb_A_to_2pdb_B_difference.npz")

        """
        ##### CHECK 3 #####
        Check function can return and save the correct distance difference matrix"""
        # Run function
        dd_matx = distance_differences.generate_matx_diff(
            mock_ca_matx1,
            mock_ca_matx2,
            # PATH_SAVE_CA.joinpath(f"{TEST_UNP}_1pdb_A_to_2pdb_B_difference"),
            # return_matx=True,
        )
        # Check returned matrix is exactly as expected
        self.assertTrue(
            np.array_equal(dd_matx, expected_dd_matx, equal_nan=True),
            msg="Function does not return the correct distance difference matrix.",
        )
        # Check a file has been saved with the expected name
        # self.assertIsFile(f"{PATH_SAVE_CA}/{TEST_UNP}_1pdb_A_to_2pdb_B_difference.npz")

    def test_make_heatmap_kwargs(self):
        """
        Test the function responsible for making heatmap kwargs.
        """

        # Run test
        test_heatmap_kwargs = distance_differences.make_heatmap_kwargs()

        # Check output is a dictionary
        self.assertIsInstance(
            test_heatmap_kwargs,
            dict,
            msg="Heatmap kwargs should be a dictionary, not "
            f"{type(test_heatmap_kwargs)}",
        )

    # Other plotting functions are not that necessary.


# Run unit tests on call of script
if __name__ == "__main__":

    unittest.main()
