"""
Contains functions for unit testing scripts
"""

import os
import pathlib
from unittest import TestCase


class TestCaseModified(TestCase):
    """
    Additional methods to expand the TestCase class

    Modifies the TestCase class from unittest, so therefore TestCaseBase is called by
    the master testing class (next).

    Modified from:
    https://stackoverflow.com/questions/59121161/python-unittest-how-to-assert-the-
    existence-of-a-file-or-folder-and-print-the-p
    """

    def assertIsFile(self, path):
        """
        Tests for the presence of a saved file.
        """
        if not pathlib.Path(path).resolve().is_file():
            raise FileNotFoundError(f"File does not exist: {str(path)}")

    def assertIsNotFile(self, path):
        if pathlib.Path(path).resolve().is_file():
            raise FileExistsError(f"File does exist: {str(path)}")


def remove_files_in_dir(path):
    """
    Utils function used at the start of test methods to remove files saved from previous
    tests.

    To save copies of previous tests, copy the folder's contents to a new location.
    """

    os.system(f"rm {str(path)}/*")
