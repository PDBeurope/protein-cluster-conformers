"""
Tools for managing IO operations for files required and created during clustering. Most
read operations pertain to the updated mmCIF file format and most write operations
relate to matrices and figures.
"""

from logging import getLogger
from os import rename
from pathlib import PosixPath
from pickle import dump, load
from subprocess import check_output

# Third party imports
from gemmi import cif
from matplotlib import pyplot as plt
from numpy import around, float32
from numpy import load as np_load
from numpy import ndarray, savez_compressed, triu

# Load logger
logger = getLogger(__name__)


def load_mmcif(path: PosixPath) -> cif.Block:
    """Returns mmCIF Gemmi object given path to file.

    :param path: Load (updated) mmCIF file into memory.
    :type path: pathlib.Path
    :raises ValueError: File parsed was not mmCIF or gzip mmCIF
    :return: Contents of mmCIF as Gemmi Block object
    :rtype: gemmi.cif.Block
    """
    path = str(path)
    if path[-6:] == "cif.gz":
        return cif.read(path).sole_block()
    elif path[-3:] == "cif":
        return cif.read_file(path).sole_block()
    else:
        raise ValueError


def save_matrix(matrix: ndarray, path: "PosixPath|str", label: str = None) -> None:
    """
    Saves a matrix in a compressed but non-numerically transformed format. Should be
    used to save matrices that contain mix data type elements, are non-symmetric, not
    N*N or contain numeric values which must have large precision preserved (e.g.
    double, float64)

    :param matrix: Matrix object to save
    :type matrix: ndarray
    :param path: Path to save file (include file name)
    :type path: PosixPath|str
    :param label: Identification label in case save is unsuccessful, defaults to None
    :type label: str, optional
    """

    try:
        savez_compressed(path, matrix)

    except Exception:
        logger.error(f"Could not save matrix: {label}", exc_info=True)


def save_compressed_matrix(
    matrix: ndarray, path: "PosixPath|str", label: str = None
) -> None:
    """
    Parse matrix and path to save it. Saves matrix as .npz file. Compression techniques
    used are specific to N*N-square, symmetric matrices. Parsing in anything else could
    result in data loss.

    :param matrix: Any N*N-square, symmetric matrix.
    :type matrix: np.ndarray
    :param path: Path to save the matrix.
    :type path: pathlib.Path|str
    """

    try:
        matrix = triu(matrix)  # Sparse array
        matrix = around(  # Round the array before saving
            matrix, 1  # Save to this many decimal places
        )

        matrix = float32(matrix)

    except Exception:
        logger.warning(
            f"Could not perform space-saving transformation on matrix ({label}), but "
            "saving anyway",
            exc_info=True,
        )

    savez_compressed(path, matrix)


def save_string_based_matx(
    matrix: ndarray, path: "PosixPath|str", label: str = None
) -> None:
    """
    Saves a string-containing matrix in an uncompressed format.

    :param matrix: Matrix filled with string elements.
    :type matrix: ndarray
    :param path: Path to save matrix.
    :type path: PosixPath|str
    :param label: Optional label, defaults to None
    :type label: str, optional
    """

    try:
        savez_compressed(path, matrix)
    except Exception:
        logger.error(f"Could not save matrix: {label}", exc_info=True)


def load_matrix(path: "PosixPath|str") -> ndarray:
    """
    Loads matrix (numpy array) into variable.

    :param path: Path to saved matrix (including file name)
    :type path: Path|str
    :return: Matrix containing in .pnz file.
    :rtype: np.ndarray
    """
    npz_obj = np_load(path, allow_pickle=True)
    matx = npz_obj["arr_0"]

    npz_obj.close()

    return matx


def get_fnames(path: "PosixPath|str") -> "list[str]":
    """
    Given a path to a dir containing files, returns a list of file names in the dir
    (including the path to their location prefixed) as strings.

    :param path: Path to folder containing files to return.
    :type path: Path|str
    :return: List of files contained in parsed path.
    :rtype: list[str]
    """

    ls = f"ls {path}"

    # Args to ensure correct formatnig of string elements
    args = {"shell": True, "encoding": "utf-8"}

    return check_output([ls], **args).splitlines()


def rename_file(path_file: PosixPath, new_fname: str) -> None:
    """
    Given a path to a file and a new file name (new_fname), the function renames the
    file. This function was made for renaming AlphaFold mmCIF structure files, to
    retain compatibility with other mmCIF parsing functions in the codebase.

    :param path_file: Path to original file name.
    :type path_file: pathlib.Path
    :param new_fname: New file name (parent path inherited from path_file).
    :type new_fname: str
    """

    rename(path_file, path_file.parent.joinpath(new_fname))


def save_figure(
    path_save: PosixPath, save_fname: str = "tmp", png: bool = True, svg: bool = False
) -> None:
    """
    Saves a matplotlib figure. Figure object does not need to be parsed in as an
    argument.

    :param path_save: Path to save figure.
    :type path_save: PosixPath
    :param save_fname: File name (saved in parent path defined by path_save), defaults
        to "tmp"
    :type save_fname: str, optional
    :param png: Save image in PNG format, defaults to True
    :type png: bool, optional
    :param svg: Save image in SVG format, defaults to False
    :type svg: bool, optional
    """
    default_dpi = 200

    if png:
        save_fig_dir = path_save.joinpath(f"{save_fname}.png")
        plt.savefig(save_fig_dir, dpi=default_dpi, bbox_inches="tight")
        logger.info(f"Figure saved as {save_fig_dir}")

    if svg:
        save_fig_dir = path_save.joinpath(f"{save_fname}.svg")
        plt.savefig(save_fig_dir, bbox_inches="tight")
        logger.info(f"Figure saved as {save_fig_dir}")

    # plt.close()


def serial_dump(object, path_fname: PosixPath) -> None:
    """
    Function for dumping a parsed object as a serialised object file.

    :param object: Python object to serialise
    :type object: Any
    :param path_fname: Path (including file name) to save serialised object
    :type path_fname: PosixPath
    """

    file = open(path_fname, "wb")
    dump(obj=object, file=file)  # pickle.dump
    file.close()


def serial_load(path_fname: PosixPath):
    """
    Loads a serialised object stored at the given path. Returns whatever object was in
    the serialised file object.

    :param path_fname: Path to serialised Python object (including file name).
    :type path_fname: PosixPath
    :return: Python object stored in serialised file.
    :rtype: Any
    """

    file = open(path_fname, "rb")
    object = load(file)  # pickle.load
    file.close()

    return object


def load_matrix_from_tri_upper(path: "PosixPath|str") -> ndarray:
    """
    Function for loading compressed, symmetric matrices in which only the upper half was
    saved.
    """
    matx_upper_tri = load_matrix(path=path)

    # square_symmetric_matx = fill_diagonal(
    #     matx_upper_tri.T + matx_upper_tri,
    #     diag(matx_upper_tri)
    # )

    return matx_upper_tri.T + matx_upper_tri
