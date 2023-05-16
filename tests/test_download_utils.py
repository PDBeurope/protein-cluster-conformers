# Third party imports
import os
import pathlib

# For testing the test script
# import sys
import unittest
from unittest import mock

import numpy as np

# Import functions to test
from cluster_conformers.utils import download_utils

# Import modified TestCase class
from .test_case import TestCaseModified

# # Modify Python path for custom imports
# sys.path[0] = str(  # Override default Python import part
#     pathlib.Path(
#         os.path.realpath(__file__)  # Path to this file
#     ).parent.parent  # Path to folder this file is in
# )


# Paths
# PATH_HOME = pathlib.Path.home()
# PATH_HOME = PATH_HOME.joinpath(
#     "EMBL-EBI", "funclan_work", "contact_map_difference", "tests"
# )
# PATH_HOME = str(PATH_HOME)
# PATH_BASE = PATH_HOME
PATH_BASE = str(pathlib.Path("./tests"))
PATH_TEST_SAVE = pathlib.Path(f"{PATH_BASE}/test_output")
MOCK_MMCIF_PATH = PATH_BASE + "/mock_data/test_mmcifs/3ns8_updated.cif"

# Mock variables and objects
MOCK_MATX = np.full((100, 100), 1.0)


class TestIOUtils(TestCaseModified):
    """
    Tests functions contained within the cluster_conformers.utils.io_utils script.
    These functions are used for basic IO operations used by worker functions to perform
    protein monomer conformational clustering.
    """

    @mock.patch("cluster_conformers.utils.download_utils.get")
    def test_fetch_benchmark(self, mock_get):
        """
        Tests function used to download all structures from the parsed benchmark dataset
        CSV.

        :param mock_get: Patched function call to API via requests.get()
        :type mock_get: MagicMock
        """
        # Location of mmCIF to mock call from API
        mock_mmcif = open("tests/test_inputs/3ns8_updated.cif", "rb").read()
        mock_get.return_value = mock.Mock(**{"content": mock_mmcif})

        # Run function
        download_utils.fetch_benchmark_mmcifs(
            path_benchmark_df="tests/test_inputs/benchmark_ds_truncated.csv",
            path_save=pathlib.Path(  # Not actually downloaded...
                "tests/test_output/downloads"
            ),
        )

        # Check if file is saved
        self.assertIsFile(path="tests/test_output/downloads/3ns8_updated.cif")

    @mock.patch("cluster_conformers.utils.download_utils.get")
    def test_download_alphafold_mmcif(self, mock_get):
        """
        Tests the function to download AlphaFold structure from URL.

        :param mock_get: Patched function call to API via requests.get()
        :type mock_get: MagicMock
        """

        # Remove file from previous tests
        os.system("rm tests/test_output/downloads/A12345/AF-A12345-F1-model_v3.cif")

        # Location of mmCIF to mock call from API
        test_af_mmcif = open("tests/test_inputs/afv3_updated.cif", "rb").read()
        mock_get.side_effect = [
            # File not found
            mock.Mock(**{"content": None, "status_code": 404}),
            # Database down or other error
            mock.Mock(**{"content": None, "status_code": 500}),
            # Connection made successfully
            mock.Mock(**{"content": test_af_mmcif, "status_code": 200}),
        ]

        # Run function -- 404 file not found
        path_to_af = download_utils.download_alphafold_mmcif(
            uniprot="A12345",
            path_save=pathlib.Path(  # Not actually downloaded...
                "test/test_output/downloads"
            ),
        )
        # Check nothing returned
        self.assertIsNone(path_to_af)
        # Check file is not saved
        self.assertIsNotFile(
            path="tests/test_output/downloads/A12345/AF-A12345-F1-model_v3.cif"
        )

        # Run function -- server error
        return_msg = """Function must be able to handle server-connection failures.
        Failed to raise anything with 500 error."""
        with self.assertRaises(ConnectionError, msg=return_msg):
            download_utils.download_alphafold_mmcif(
                uniprot="A12345",
                path_save=pathlib.Path(  # Not actually downloaded...
                    "tests/test_output/downloads"
                ),
            )
        # Check file is not saved
        self.assertIsNotFile(
            path="tests/test_output/downloads/A12345/AF-A12345-F1-model_v3.cif"
        )

        # Run function -- healthy connection, 200
        path_to_af = download_utils.download_alphafold_mmcif(
            uniprot="A12345",
            path_save=pathlib.Path(  # Not actually downloaded...
                "tests/test_output/downloads"
            ),
        )
        # Check file is saved
        self.assertIsFile(
            path="tests/test_output/downloads/A12345/AF-A12345-F1-model_v3.cif"
        )
        # Check expected string is returned
        self.assertEqual(
            path_to_af,
            pathlib.Path(
                "tests/test_output/downloads/A12345/AF-A12345-F1-model_v3.cif"
            ),
        )


# Run unit tests on call of script
if __name__ == "__main__":

    unittest.main()
