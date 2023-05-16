"""
Some simplifications needed
"""

# import os

# Third party imports
import pathlib

# For testing the test script
# import sys
import unittest
from unittest.mock import MagicMock, patch

# import numpy as np
from gemmi import cif
from numpy import nan

# Import functions to test
from cluster_conformers.utils import parsing_utils

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
MOCK_MMCIF_PATH = PATH_BASE + "/mock_data/mock_mmcifs/mock_structure.cif"
# ""
MOCK_MMCIF = cif.read_file(MOCK_MMCIF_PATH).sole_block()

# Mock objects
mock_table = [  # Test table with HETAMS, different chains and non-CAs
    ["ATOM", "CA", "1.0", "2.0", "3.0", "A", "1"],
    ["ATOM", "C", "4.0", "5.0", "6.0", "A", "1"],
    ["ATOM", "CA", "10.0", "20.0", "30.0", "A", "2"],
    ["ATOM", "N", "7.0", "8.0", "9.0", "A", "1"],
    ["ATOM", "CA", "100.0", "200.0", "300.0", "A", "3"],
    ["ATOM", "CA", "-100.0", "-200.0", "-300.0", "A", "4"],
    ["HETATM", "C", "-100.0", "-200.0", "-300.0", "A", "4"],
    ["ATOM", "CA", "1.0", "2.0", "3.0", "B", "1"],
    ["ATOM", "CA", "10.0", "20.0", "30.0", "B", "2"],
    ["ATOM", "CA", "100.0", "200.0", "300.0", "B", "3"],
    ["ATOM", "CA", "-100.0", "-200.0", "-300.0", "B", "4"],
]


class TestParsingUtils(TestCaseModified):
    """
    Object for testing all functions in cluster_conformers.utils.parsing_utils.py .

    :param TestCaseModified: _description_
    :type TestCaseModified: _type_
    """

    def test_extract_table(self):
        """
        Tests for the function responsible for extracting row-column information from a
        parsed updated mmCIF file (Gemmi object in this case).
        """

        # All search terms in _atom_site loop of updated mmCIF as of 16/09/2022
        search_terms = [
            "group_PDB",
            "id",
            "type_symbol",
            "label_atom_id",
            "label_alt_id",
            "label_comp_id",
            "label_asym_id",
            "label_entity_id",
            "label_seq_id",
            "pdbx_PDB_ins_code",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
            "pdbx_formal_charge",
            "auth_seq_id",
            "auth_comp_id",
            "auth_asym_id",
            "auth_atom_id",
            "pdbx_PDB_model_num",
            "pdbx_label_index",
            "pdbx_sifts_xref_db_name",
            "pdbx_sifts_xref_db_acc",
            "pdbx_sifts_xref_db_num",
            "pdbx_sifts_xref_db_res",
        ]

        # Check all search terms can be queried individually
        for term in search_terms:
            test_table = parsing_utils.extract_table(
                mmcif=MOCK_MMCIF, search_list=[term]
            )
            self.assertIsInstance(
                test_table,
                cif.Table,
                msg=f"'{term}' loop could not be extracted from updated mmCIF file.",
            )
            self.assertTrue(
                len(test_table) != 0,
                msg="cif.Table object was returned but object empty.",
            )
            self.assertEqual(
                len(test_table),
                795,
                msg=f"Returned cif.Table is {len(test_table)} rows in length, not 795",
            )

        # Check all search terms can be queried simultaneously
        test_mega_searched_table = parsing_utils.extract_table(
            mmcif=MOCK_MMCIF, search_list=search_terms
        )
        self.assertEqual(  # Correct number of columns extracted
            test_mega_searched_table.width(),
            26,
            msg="Should have 26 extracted columns, not "
            f"{test_mega_searched_table.width()}",
        )
        self.assertEqual(  # Correct number of rows extracted
            len(test_mega_searched_table),
            795,
            msg="Should have 795 extracted rows, not "
            f"{len(test_mega_searched_table)}",
        )

        # Check search for only terms relevant to analyses here works
        # test_select_terms = parsing_utils.extract_table(
        #     mmcif=MOCK_MMCIF,
        #     search_list=[
        #         "group_PDB",
        #         "label_atom_id",
        #         "Cartn_x",
        #         "Cartn_y",
        #         "Cartn_z",
        #         "label_asym_id",
        #         "pdbx_sifts_xref_db_num",
        #     ]
        # )

    def test_fill_missing_unps(self):
        """
        Test the function that inserts missing UniProt residue indices and fills
        pertaining coordinate lists with np.NaN values.

        Inputs to tested function:
            3 lists of equal length
        Returned from tested function:
            3 lists of equal length. Length equals max UniProt residue index (in first
            parsed list)
        """

        mock_input_dict = {
            "unp_res_ids": {3, 4, 5, 7, 9, 10},
            "cartn_x": [10.0, 1.41, 3.33, 49.1, 100.1, 72.02],
            "cartn_y": [30.37, 44.16, 81.86, 53.95, 63.09, 21.7],
            "cartn_z": [72.46, 8.15, 31.35, 89.18, 23.32, 26.19],
        }

        # Run test
        test_output_dict = parsing_utils.fill_missing_unps(
            structure_coords=mock_input_dict
        )

        self.assertDictEqual(
            test_output_dict,
            {
                "unp_res_ids": {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                "cartn_x": [nan, nan, 10.0, 1.41, 3.33, nan, 49.1, nan, 100.1, 72.02],
                "cartn_y": [
                    nan,
                    nan,
                    30.37,
                    44.16,
                    81.86,
                    nan,
                    53.95,
                    nan,
                    63.09,
                    21.7,
                ],
                "cartn_z": [
                    nan,
                    nan,
                    72.46,
                    8.15,
                    31.35,
                    nan,
                    89.18,
                    nan,
                    23.32,
                    26.19,
                ],
            },
            msg="Missing UniProt residue entries are not being filled or not being "
            "filled correctly. ",
        )

    @patch(
        "cluster_conformers.utils.parsing_utils.extract_table",
        side_effect=[
            # So long as the patched function returns a N*7-dimensional iterable, the
            # remainder of the function should work as of 20/9/22
            [  # Test table without any problematic values
                ["ATOM", "CA", "1.0", "2.0", "3.0", "A", "1"],
                ["ATOM", "CA", "10.0", "20.0", "30.0", "A", "2"],
                ["ATOM", "CA", "100.0", "200.0", "300.0", "A", "3"],
                ["ATOM", "CA", "-100.0", "-200.0", "-300.0", "A", "4"],
            ],
            [],  # Test where the length of the table is zero
            # Test whether excludable are in fact excluded
            mock_table,
            mock_table,
        ],
    )
    def test_parse_mmcif(self, mock_extract_table):
        """
        Test for the parse_mmcif() function. The function extract_table(), which is
        tested above, is called within this function and should therefore be patched.

        :param mock_extract_table: Mock N*7 iterable, where N=number of rows and
            7=number of columns.
        :type mock_extract_table: Iterable
        :raises TypeError: Should raise TypeError when emtpy table parsed.
        """

        # Null mock -- fed into patched function so doesn't need to have any value
        mock_mmcif = MagicMock(
            spec=cif.Block, return_value=None  # Function to call this obj is patched
        )

        # Expected output from parsed mmCIF.
        expected_dict = {
            "unp_res_ids": {1, 2, 3, 4},
            "cartn_x": [1.0, 10.0, 100.0, -100.0],
            "cartn_y": [2.0, 20.0, 200.0, -200.0],
            "cartn_z": [3.0, 30.0, 300.0, -300.0],
        }

        # Test function on a 'perfect' situation -- where there are no exceptions
        self.assertDictEqual(
            parsing_utils.parse_mmcif(mock_mmcif, "A"),
            expected_dict,
            msg="Function cannot handle a situation in which there are no exceptions "
            "to mmCIF parsing.",
        )

        # Check error message is returned if the parsed table is empty
        return_msg = "Parsing in an empty cif.Table should raise TypeError"
        with self.assertRaises(TypeError, msg=return_msg):
            parsing_utils.parse_mmcif(mock_mmcif, "A")

        # Check whether a specific chain can be pulled out of a table containing
        # multiple chains
        self.assertDictEqual(
            parsing_utils.parse_mmcif(mock_mmcif, "B"),
            expected_dict,
            msg="Function cannot pull out information pertaining to a single chain "
            "from a multi-chain structure.",
        )

        # Check whether a specific chain can be pulled out of a table containing
        # multiple chains and non-CA information
        self.assertDictEqual(
            parsing_utils.parse_mmcif(mock_mmcif, "A"),
            expected_dict,
            msg="Function cannot pull out information pertaining to a single chain "
            "from a multi-chain structure, where other information must be ignored.",
        )


# Run unit tests on call of script
if __name__ == "__main__":

    unittest.main()
