# Third party imports
import os
import pathlib

# For testing the test script
# import sys
import unittest

import gemmi
import numpy as np
from matplotlib import pyplot as plt

# Import functions to test
from cluster_conformers.utils import io_utils

# Import modified TestCase class
from .test_case import TestCaseModified, remove_files_in_dir

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
MOCK_MMCIF_PATH = PATH_BASE + "/mock_data/test_mmcifs/3nc3_updated.cif"

# Mock variables and objects
MOCK_MATX = np.full((100, 100), 1.0)


class TestIOUtils(TestCaseModified):
    """
    Tests functions contained within the cluster_conformers.utils.io_utils script.
    These functions are used for basic IO operations used by worker functions to perform
    protein monomer conformational clustering.
    """

    def test_load_mmcif(self):
        """
        Checks if a Gemmi Block object is loaded.
        """
        # Run the function on uncompressed file
        test_cif = io_utils.load_mmcif(
            pathlib.Path("tests/mock_data/mock_mmcifs/mock_structure.cif")
        )

        # Check returned object is a Gemmi Block
        self.assertIsInstance(
            test_cif,
            gemmi.cif.Block,
            msg="Returned object is not a Gemmi Block object.",
        )

        # Run the function on compressed file
        test_cif = io_utils.load_mmcif(
            pathlib.Path("tests/mock_data/mock_mmcifs/mock_structure_compressed.cif.gz")
        )

        # Check returned object is a Gemmi Block
        self.assertIsInstance(
            test_cif,
            gemmi.cif.Block,
            msg="Returned object is not a Gemmi Block object.",
        )

        # Run the function on a non-mmCIF file: should raise error
        return_msg = """The parsed file was in .txt format and should not have been read
                        by the function."""

        with self.assertRaises(ValueError, msg=return_msg):
            io_utils.load_mmcif(
                pathlib.Path("tests/mock_data/mock_mmcifs/mock_structure_non_mmcif.txt")
            )

    def test_save_matrix_numeric(self):
        """
        Check the save_matrix() function saves a symmetric, numeric matrix to a
        compressed .npz file.
        """
        # Remove output from previous test
        remove_files_in_dir(PATH_TEST_SAVE)

        # Location to save data produced in test
        path_save = PATH_TEST_SAVE.joinpath("test_matrix_numeric")

        # Run function -- not Python2 compatible
        io_utils.save_matrix(MOCK_MATX, path_save)

        # Check matrix has been saved
        self.assertIsFile(str(path_save) + ".npz")

    def test_save_matrix_text(self):
        """
        Check the save_matrix() function saves a symmetric, test-based matrix to a
        compressed .npz file.
        """
        # Remove output from previous test
        remove_files_in_dir(PATH_TEST_SAVE)

        # Location to save data produced in test
        path_save = PATH_TEST_SAVE.joinpath("test_matrix_text")

        # Create a test matrix containing strings only
        mock_text_matx = np.full((3, 3), "A")

        # Run function -- not Python2 compatible
        io_utils.save_matrix(mock_text_matx, path_save)

        # Check matrix has been saved
        self.assertIsFile(str(path_save) + ".npz")

    def test_load_matrix(self):
        """
        Test function which simply loads in a matrix, given a path. Should be able to
        accept path as a string or pathlib object.
        """
        # Location of mock matrix
        path_mock_matx = f"{PATH_BASE}/mock_data/test_matrix.npz"

        # Check function loads in the correct matrix, given string path
        self.assertTrue(
            np.array_equal(io_utils.load_matrix(path_mock_matx), MOCK_MATX),
            msg="Matrices not being loaded correctly.",
        )

        # Check function loads in the correct matrix, given pathlib object
        self.assertTrue(
            np.array_equal(
                io_utils.load_matrix(pathlib.Path(path_mock_matx)), MOCK_MATX
            ),
            msg="Matrix not loaded correctly given pathlib.Path object.",
        )

    def test_get_fnames(self):
        """
        Checks if the get_fnames() function retrieves the correct list of file names,
        given a test directory as a string and as a pathlib object.
        """

        # Expected output from function
        expected_output_list = [
            "mock_matx_1.npy",
            "mock_matx_2.npy",
            "mock_matx_3.npy",
            "mock_matx_4.npy",
            "mock_matx_5.npy",
        ]

        # Check function works with string path
        self.assertListEqual(
            io_utils.get_fnames(f"{PATH_BASE}/mock_data/mock_matrices/"),
            expected_output_list,
            msg="Cannot get file names using string as path.",
        )

        # Check function works with pathlib object path
        self.assertListEqual(
            io_utils.get_fnames(pathlib.Path(f"{PATH_BASE}/mock_data/mock_matrices/")),
            expected_output_list,
            msg="Cannot get file names using pahtlib object as path.",
        )

    # def test_rename_file(self):

    #     path_to_rename_file = "test/test_outputs/renaming/file_last_renamed.txt"

    #     time_last_modified = os.path.getmtime(
    #         path_to_rename_file
    #     )

    #     io_utils.rename_file(
    #         path_file=pathlib.Path(path_to_rename_file),
    #         new_fname="file_last_renamed.txt"
    #     )

    #     print(time_last_modified)
    #     print(os.path.getmtime(path_to_rename_file))

    #     self.assertTrue(
    #         time_last_modified < os.path.getmtime(path_to_rename_file)
    #     )

    def test_save_figure(self):
        """
        Checks if the save_figure() function correctly saves a figure.
        """
        # Remove output test files from previous tests, before running this test
        remove_files_in_dir(PATH_TEST_SAVE)

        # Check general plt.scatter case
        # Data -- actual values are not important
        x_data = range(0, 11)
        y_data = range(0, 11)

        # Make new plot, using the plt.subplots() method and save -- run function
        _, ax = plt.subplots(1, 1)
        ax.scatter(x_data, y_data)
        io_utils.save_figure(
            PATH_TEST_SAVE, save_fname="simple_scatter_test", png=True, svg=True
        )  # Function should close figure

        # Check if PNG version saved
        self.assertIsFile(PATH_TEST_SAVE.joinpath("simple_scatter_test.png"))

        # Check if SVG version is saved
        self.assertIsFile(PATH_TEST_SAVE.joinpath("simple_scatter_test.svg"))

    def test_serial_dump(self):
        """
        Test for function which should perform serialisation and storing of an object to
        hard drive location.
        """

        # Object and path to serialise
        test_list = [1, 2, 3, 4, 5]
        path_to_serialised_file = "tests/test_output/test_serialised_list.pickle"

        # Remove file from previous tests
        os.system(f"rm {path_to_serialised_file}")

        # Run function -- test
        io_utils.serial_dump(test_list, path_to_serialised_file)

        # Check file created
        self.assertIsFile(path_to_serialised_file)

    def test_serial_load(self):
        """
        Tests function responsible for loading in a serialised file stored on hard
        drive.
        """

        # Path to input serialised list (same as object above but in a diff location)
        path_to_serialised_file = "tests/test_inputs/test_serialised_list.pickle"

        # Run function -- test
        test_serialised_list = io_utils.serial_load(path_to_serialised_file)

        # Check correctly loaded into memory
        self.assertListEqual(
            test_serialised_list,
            [1, 2, 3, 4, 5],
            msg="Function not loading serialised file correctly.",
        )


# Run unit tests on call of script
if __name__ == "__main__":

    unittest.main()
