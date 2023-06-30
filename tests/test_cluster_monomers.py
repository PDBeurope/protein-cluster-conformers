"""
Unit tests for the main object which handles clustering of monomeric protein chains
"""

# Third party imports
# import os
import pathlib

import gemmi

# import numpy as np
import pandas as pd

# Import functions/classes to test
from cluster_conformers import cluster_monomers

# Import modified TestCase class
from .test_case import TestCaseModified, remove_files_in_dir

# Global variables
PATH_BASE = pathlib.Path("./tests")
PATH_SERIALISED_MATXS = PATH_BASE.joinpath("mock_serial_matxs")
PATH_SERIALISED_10_10_MATXS = PATH_SERIALISED_MATXS.joinpath("10_by_10_full")

# Initialise other paths
PATH_SAVE_OUTPUT = PATH_BASE.joinpath("test_output")
PATH_SAVE_CA = PATH_SAVE_OUTPUT.joinpath("ca_distances")
PATH_SAVE_CLUSTER_RESULTS = PATH_SAVE_OUTPUT.joinpath("cluster_results")
PATH_SAVE_DD_MATXS = PATH_SAVE_OUTPUT.joinpath("distance_differences")
PATH_SAVE_DD_MAPS = PATH_SAVE_OUTPUT.joinpath("distance_difference_maps")

PATH_TRUNCATED_MOCK_MMCIFS = PATH_BASE.joinpath("mock_data/mock_mmcifs/")

# Mock objects
LOADED_MMCIF_OBJ_TYPE = gemmi.cif.Block
CHAINS = ["A", "B"]
TEST_UNP = "A12345"

# Test input
TEST_MMCIFS_AND_CHAINS_DICT = {
    str(PATH_TRUNCATED_MOCK_MMCIFS.joinpath("3nc3_updated.cif")): CHAINS,
    str(PATH_TRUNCATED_MOCK_MMCIFS.joinpath("3nc5_updated.cif")): CHAINS,
    str(PATH_TRUNCATED_MOCK_MMCIFS.joinpath("3nc6_updated.cif")): CHAINS,
    str(PATH_TRUNCATED_MOCK_MMCIFS.joinpath("3nc7_updated.cif")): CHAINS,
}


class TestClusterMonomersClassicalForce(TestCaseModified):
    def setUp(self):
        """
        Creates a simplistic class instance of a set of truncated protein chains. The
        instance is used by all functions, which are subsequently patched where
        function calls are tested elsewhere.
        """

        # Initialise common object instance
        self.test_cluster_conformers_obj = cluster_monomers.ClusterConformations(
            unp=TEST_UNP, mmcifs_and_chains=TEST_MMCIFS_AND_CHAINS_DICT, force=True
        )

        # Create folders if not already present
        for path in (
            PATH_SAVE_CA,
            PATH_SAVE_CLUSTER_RESULTS,
            PATH_SAVE_DD_MATXS,
            PATH_SAVE_DD_MAPS,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        """
        Flush class instance after each test.
        """
        return super().tearDown()

    ####################################################################################
    # @mock.patch(
    #     "cluster_conformers.cluster_monomers.ClusterConformations.build_clustering_inputs",
    # )
    # @mock.patch(
    #     "cluster_conformers.cluster_monomers.cluster_chains.cluster_agglomerative"
    # )
    # @mock.patch("cluster_conformers.cluster_monomers.cluster_chains.make_linkage_matx")
    ####################################################################################
    def test_cluster(self):
        """
        Tests the main function of the class, which clusters conformations. Previously,
        each private method was tested separately, but this was deemed unnecessary as
        those should (for the mosst part) not be called directly.
        """

        # Remove saved files from previous tests
        remove_files_in_dir(PATH_SAVE_CLUSTER_RESULTS)

        self.test_cluster_conformers_obj.ca_distance(
            path_save=PATH_SAVE_CA,
        )

        # Run the method
        self.test_cluster_conformers_obj.cluster(
            path_save_cluster_results=PATH_SAVE_CLUSTER_RESULTS,
            path_save_dd_matx=PATH_SAVE_DD_MATXS,
        )

        # Check a Pandas DataFrame is returned
        self.assertTrue(
            type(self.test_cluster_conformers_obj.cluster_df) == pd.DataFrame,
            msg="Clustering results are not being saved as a Pandas DataFrame",
        )

        # Check minimum number of columns have been created
        for col_name in ["UNP_ACC", "PDBe_ID", "CHAIN_ID", "CONFORMER_ID"]:
            self.assertTrue(
                col_name in self.test_cluster_conformers_obj.cluster_df.columns,
                msg=f"Clustering DataFrame is missing column {col_name}",
            )

        # Check at least one row (default = first) is correctly formatted
        row_index = 0
        self.assertEqual(  # UniProt accession
            self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["UNP_ACC"],
            TEST_UNP,
            msg="UniProt accession has not been stored correctly",
        )

        self.assertTrue(  # Check PDB ID is stored as a string
            type(self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["PDBe_ID"])
            == str,
            msg="PDB ID has not been stored as a string",
        )

        self.assertTrue(  # Check chain ID stored as a single string
            type(
                self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            )
            == str
            and len(
                self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            )
            == 1,
            msg="Chain ID has not been stored as a single character",
        )

        self.assertTrue(  # Check chain ID stored as a single string
            self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            != self.test_cluster_conformers_obj.cluster_df.iloc[row_index][
                "CHAIN_ID"
            ].lower(),
            msg="Chain ID has not been stored in upper case. Important to keep this "
            "consistent for downstream analyses",
        )

        # Expected files saved from clustering process
        expected_output_files = [
            f"{TEST_UNP}_sum_based_clustering_results.csv",
            f"{TEST_UNP}_score_matrix.npz",
            f"{TEST_UNP}_label_matrix.npz",
            f"{TEST_UNP}_linkage_matrix.npz",
        ]

        # Check expected files have been saved
        for file in expected_output_files:
            path_to_expected_file = PATH_SAVE_CLUSTER_RESULTS.joinpath(file)
            self.assertIsFile(path_to_expected_file)

    def test_cluster_with_updated_entries(self):
        """
        Tests the main function of the class in the case where a several entries have
        been updated. In this case, another
        """

        # Remove saved files from previous tests
        remove_files_in_dir(PATH_SAVE_CLUSTER_RESULTS)

        # Only method truly being tested here -- the rest were checked in the previous
        # test by running them is still important to ensure the correct results are
        # generated and file saved
        self.test_cluster_conformers_obj.remove_entry_matxs(
            pdb_ids=["3nc3", "3nc7"],
            path_ca=PATH_SAVE_CA,
            path_dd=PATH_SAVE_DD_MATXS,
        )

        self.test_cluster_conformers_obj.ca_distance(
            path_save=PATH_SAVE_CA,
        )

        # Run the method
        self.test_cluster_conformers_obj.cluster(
            path_save_cluster_results=PATH_SAVE_CLUSTER_RESULTS,
            path_save_dd_matx=PATH_SAVE_DD_MATXS,
        )

        # Check a Pandas DataFrame is returned
        self.assertTrue(
            type(self.test_cluster_conformers_obj.cluster_df) == pd.DataFrame,
            msg="Clustering results are not being saved as a Pandas DataFrame",
        )

        # Check minimum number of columns have been created
        for col_name in ["UNP_ACC", "PDBe_ID", "CHAIN_ID", "CONFORMER_ID"]:
            self.assertTrue(
                col_name in self.test_cluster_conformers_obj.cluster_df.columns,
                msg=f"Clustering DataFrame is missing column {col_name}",
            )

        # Check at least one row (default = first) is correctly formatted
        row_index = 0
        self.assertEqual(  # UniProt accession
            self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["UNP_ACC"],
            TEST_UNP,
            msg="UniProt accession has not been stored correctly",
        )

        self.assertTrue(  # Check PDB ID is stored as a string
            type(self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["PDBe_ID"])
            == str,
            msg="PDB ID has not been stored as a string",
        )

        self.assertTrue(  # Check chain ID stored as a single string
            type(
                self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            )
            == str
            and len(
                self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            )
            == 1,
            msg="Chain ID has not been stored as a single character",
        )

        self.assertTrue(  # Check chain ID stored as a single string
            self.test_cluster_conformers_obj.cluster_df.iloc[row_index]["CHAIN_ID"]
            != self.test_cluster_conformers_obj.cluster_df.iloc[row_index][
                "CHAIN_ID"
            ].lower(),
            msg="Chain ID has not been stored in upper case. Important to keep this "
            "consistent for downstream analyses",
        )

        # Expected files saved from clustering process
        expected_output_files = [
            f"{TEST_UNP}_sum_based_clustering_results.csv",
            f"{TEST_UNP}_score_matrix.npz",
            f"{TEST_UNP}_label_matrix.npz",
            f"{TEST_UNP}_linkage_matrix.npz",
        ]

        # Check expected files have been saved
        for file in expected_output_files:
            path_to_expected_file = PATH_SAVE_CLUSTER_RESULTS.joinpath(file)
            self.assertIsFile(path_to_expected_file)


class TestFigureRendering(TestCaseModified):
    def setUp(self):
        """
        Creates a simplistic class instance of a set of truncated protein chains. The
        instance is used by all functions, which are subsequently patched where
        function calls are tested elsewhere.
        """

        self.path_save = PATH_SAVE_OUTPUT.joinpath("figures")
        # Create folder if not already present
        self.path_save.mkdir(parents=True, exist_ok=True)

    def tearDown(self) -> None:
        return super().tearDown()

    def test_render_dendrogram(self):
        """
        The render_dendrogram() function is now run as this current instance of the
        ClusterConformations() object now contains the prerequisite attributes after
        executing ClusterConformations().cluster()
        """

        # Input data here is a copy of clustering results for UniProt=O34926
        test_unp = "B12345"

        # Run method
        cluster_monomers.render_dendrogram(
            unp=test_unp,
            path_results=PATH_BASE.joinpath("test_inputs"),
            path_save=self.path_save,
            png=True,
            svg=True,
        )

        # Check PNG file has been created
        self.assertIsFile(
            PATH_SAVE_OUTPUT.joinpath(
                "figures", f"{test_unp}_agglomerative_dendrogram.png"
            )
        )

        # Check SVG file has been created
        self.assertIsFile(
            PATH_SAVE_OUTPUT.joinpath(
                "figures", f"{test_unp}_agglomerative_dendrogram.svg"
            )
        )
