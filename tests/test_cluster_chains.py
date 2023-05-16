"""
Script for testing all functions related to handling mid-level matrix manipulation for
clustering monomeric protein conformers.
"""

# import os
import pathlib

# For testing the test script
# import sys
import unittest
from unittest import mock

import numpy as np

from cluster_conformers import cluster_chains

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

MOCK_SCORE_MATX = np.array(
    [
        [0.0, 33000.0, 50000.0, 25000.0, 32000.0, 45000.0, 23000.0],
        [33000.0, 0.0, 41000.0, 24000.0, 29000.0, 26626.0, 33219.0],
        [50000.0, 41000.0, 0.0, 22000.0, 30000.0, 44000.0, 21000.0],
        [25000.0, 24000.0, 22000.0, 0.0, 27000.0, 24000.0, 24000.0],
        [32000.0, 29000.0, 30000.0, 27000.0, 0.0, 40000.0, 23000.0],
        [45000.0, 26626.0, 44000.0, 24000.0, 40000.0, 0.0, 30000.0],
        [23000.0, 33219.0, 21000.0, 24000.0, 23000.0, 30000.0, 0.0],
    ]
)

MOCK_MODEL_DISTANCES = np.array(
    [21000.0, 23000.0, 26626.0, 26666.66666667, 32500.0, 34321.9]
)

MOCK_MODEL_CHILDREN = np.array([[2, 6], [3, 7], [1, 5], [4, 8], [0, 10], [9, 11]])

MOCK_EXPECTED_MODEL = mock.Mock(
    name="Mock AgglomerativeClustering model",
    **{
        "distances_": MOCK_MODEL_DISTANCES,
        # "labels_": np.array( [5, 4, 6, 2, 3, 1, 0] ),
        "children_": MOCK_MODEL_CHILDREN,
    },
)


class TestClusterChains(TestCaseModified):
    """
    Tests for functions in the cluster.py script

    TODO:
     - Add test to assess clustering success
        -- might need an example with more chains
    """

    # @mock.patch("cluster_conformers.cluster_chains.io_utils.serial_load")
    # @mock.patch(
    #     "cluster_conformers.cluster_chains.distance_differences.generate_matx_diff"
    # )
    # @mock.patch("cluster_conformers.cluster_chains.io_utils.load_matrix")
    # @mock.patch("cluster_conformers.cluster_chains.linear_algebra_utils.calc_score")
    # def test_build_inputs(
    #     self,
    #     mock_calc_score,
    #     mock_load_matrix,
    #     mock_distance_differences,
    #     mock_serial_load,
    # ):
    #     """
    #     Test for the function that is responsible for building the input to the
    #     clustering function. It tests if the score and label matrices are generated as
    #     expected.

    #     :param mock_calc_score: Mock output to the calc_score() function
    #     :type mock_calc_score: list[int]
    #     :param mock_load_matrix: Mock output from the load_matrix() function, which
    #         should load a compressed, serialised file. Set to None here because the
    #         score variable is patched instead.
    #     :type mock_load_matrix: None
    #     :param mock_distance_differences: Mock distance difference matrix. Set to None
    #         as it simply gets parsed into the score function.
    #     :type mock_distance_differences: None
    #     :param mock_serial_load: Mock return obj from the function that loads the list
    #         of UniProt residue IDs. Returns None as just parsed into score function.
    #     :type mock_serial_load: None
    #     """

    #     # Patched outputs
    #     mock_serial_load.return_value = None
    #     mock_distance_differences.return_value = None
    #     mock_load_matrix.return_value = None
    #     mock_calc_score.side_effect = [9, 18, 9]

    #     test_score_matx, test_label_matx = cluster_chains.build_inputs(
    #         unp="A12345",
    #         pdbe_chain_ids=["1pdb_A", "2pdb_B", "3pdb_C"],
    #         unp_res_ids={
    #             "1pdb_A": pathlib.Path(),
    #             "2pdb_B": pathlib.Path(),
    #             "3pdb_C": pathlib.Path(),
    #         },
    #         path_save=PATH_SAVE_CA,
    #         ca_matxs={
    #             "1pdb_A": pathlib.Path(),
    #             "2pdb_B": pathlib.Path(),
    #             "3pdb_C": pathlib.Path(),
    #         },
    #     )

    #     # Check score matrix is generated correctly
    #     self.assertTrue(
    #         np.array_equal(
    #             test_score_matx,
    #             np.array([[0.0, 9.0, 18.0], [9.0, 0.0, 9.0], [18.0, 9.0, 0.0]]),
    #         ),
    #         msg="Output score matrix does not match the expected result. Check whether "
    #         "this test has been patched meaningfully if refactoring performed is "
    #         "considerable. ",
    #     )

    #     # Check label matrix is generated correctly
    #     self.assertTrue(
    #         np.array_equal(
    #             test_label_matx,
    #             np.array(
    #                 [
    #                     ["1pdb_A_to_1pdb_A", "1pdb_A_to_2pdb_B", "1pdb_A_to_3pdb_C"],
    #                     ["2pdb_B_to_1pdb_A", "2pdb_B_to_2pdb_B", "2pdb_B_to_3pdb_C"],
    #                     ["3pdb_C_to_1pdb_A", "3pdb_C_to_2pdb_B", "3pdb_C_to_3pdb_C"],
    #                 ]
    #             ),
    #         ),
    #         msg="Output label matrix is not being generated correctly. Each element "
    #         "should pair with the corresponding element in the score matrix. I.e. "
    #         "an N*N score matrix should have a corresponding N*N label matrix. ",
    #     )

    def test_linkage_matx(self):
        """
        Tests the function that returns a linkage matrix from a SKLearn model.
        """

        # Adding attributes to expected input
        mock_model = mock.Mock(
            name="Mock AgglomerativeClustering model",
            **{
                "distances_": np.array(
                    [21000.0, 23000.0, 26626.0, 26666.66666667, 32500.0, 34321.9]
                ),
                "labels_": np.array([5, 4, 6, 2, 3, 1, 0]),
                "children_": np.array(
                    [[2, 6], [3, 7], [1, 5], [4, 8], [0, 10], [9, 11]]
                ),
            },
        )

        # Test linkage matrix making function
        self.assertTrue(
            np.allclose(
                cluster_chains.make_linkage_matx(mock_model),
                np.array(
                    [  # Mock linkage matrix
                        [2.0, 6.0, 2.10000000e04, 2.0],
                        [3.0, 7.0, 2.30000000e04, 3.0],
                        [1.0, 5.0, 2.66260000e04, 2.0],
                        [4.0, 8.0, 2.66666667e04, 4.0],
                        [0.0, 10.0, 3.25000000e04, 5.0],
                        [9.0, 11.0, 3.43219000e04, 7.0],
                    ]
                ),
            )
        )

    def test_cluster_agglomerative_default_cutoff(self):
        """
        Tests the function which calculates the distance-based score used for
        agglomerative clustering.

        Also tests the make_linkage_matx() function, which is called by
        cluster_agglomerative().

        TODO: Add way to check clustering success -- compare against example
        """
        # Expected clustering labels
        MOCK_EXPECTED_MODEL.labels_ = np.array([5, 4, 6, 2, 3, 1, 0])

        test_model = cluster_chains.cluster_agglomerative(MOCK_SCORE_MATX, cutoff=None)
        # Only test whether distances, children and labels attributes return the
        # expected output
        self.assertTrue(
            np.allclose(MOCK_EXPECTED_MODEL.distances_, test_model.distances_)
        )

        self.assertTrue(
            np.allclose(MOCK_EXPECTED_MODEL.children_, test_model.children_)
        )

        self.assertTrue(np.allclose(MOCK_EXPECTED_MODEL.labels_, test_model.labels_))

    def test_cluster_agglomerative_70pc_cutoff(self):
        """
        Tests the function which calculates the distance-based score used for
        agglomerative clustering.

        Also tests the make_linkage_matx() function, which is called by
        cluster_agglomerative().

        TODO: Add way to check clustering success -- compare against example
        """
        # Expected clustering labels
        MOCK_EXPECTED_MODEL.labels_ = np.array([2, 4, 0, 0, 3, 1, 0])

        test_model = cluster_chains.cluster_agglomerative(
            MOCK_SCORE_MATX, cutoff=0.7  # or None
        )
        # Only test whether distances, children and labels attributes return the
        # expected output
        self.assertTrue(
            np.allclose(MOCK_EXPECTED_MODEL.distances_, test_model.distances_)
        )

        self.assertTrue(
            np.allclose(MOCK_EXPECTED_MODEL.children_, test_model.children_)
        )

        self.assertTrue(np.allclose(MOCK_EXPECTED_MODEL.labels_, test_model.labels_))

    if __name__ == "__main__":

        unittest.main()
