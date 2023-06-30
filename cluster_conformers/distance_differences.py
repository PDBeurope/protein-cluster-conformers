"""
Calculates the difference in internal CA distances between all CA distance matrices in
the set of pre-calculated CA distance matrices defined as being in different
conformations.

A table containing conformation IDs must be parsed in order to extract the PDBe and
chain IDs of structures pertaining to the same UniProt accession and reqion but solved
with different conformations.

This script shares functions with intracluster_distance_difference.py, which calculates
the same distance difference matrices but for structures within the same conformation.

Generates distance difference maps as 2D histograms to visualise the residue-residue
differences between protein structures, within and between conformational states.
"""

from pathlib import PosixPath

# Standard package imports
import seaborn as sns
from matplotlib import colormaps
from matplotlib import pyplot as plt
from numpy import ndarray, where
import logging

# Custom module imports
from .utils import io_utils, linear_algebra_utils

logger = logging.getLogger(__name__)

# # To suppress Matplotlib debug messages
# import os
# os.environ["QT_LOGGING_RULES"] = "qt5ct.debug=false"
# import warnings
# warnings.filterwarnings("ignore")


def generate_matx_diff(
    ca_matx1: ndarray,
    ca_matx2: ndarray,
    path_save: "PosixPath|str" = None,
) -> "ndarray|None":
    """
    Executes all the functions needed for calculating and saving the difference matrix
    between two matrices. This function does not assume each matrix is identically
    sized but does assume that each is square and all residues account for in the
    matrices map in 1-1, with only the final rows/columns of the larger matrix being
    removed.

    :param ca_matx1: First matrix
    :type ca_matx1: ndarray
    :param ca_matx2: Second matrix
    :type ca_matx2: ndarray
    :param path_save: Path to save location, defaults to None
    :type path_save: PosixPath|str, optional
    :param return_matx: Whether to return the difference matrix, defaults to False
    :type return_matx: bool, optional
    :return: Returns difference matrix if `return_matx=True`, else the distance matrix
        is simply saved to `path_save` (if parsded in) and nothing is returned.
    :rtype: ndarray|None
    """

    # Remove the final n rows/cols from the larger of the two matrices.
    ca_matx1, ca_matx2 = linear_algebra_utils.trim_to_smallest(ca_matx1, ca_matx2)

    # Calculate difference matrix
    diff_matx = linear_algebra_utils.matx_subtract(ca_matx1, ca_matx2)

    # Set all values < 3 (cutoff) Angstroms to zero
    cutoff = 3.0
    diff_matx = where(diff_matx < cutoff, 0, diff_matx)

    if path_save:
        io_utils.save_compressed_matrix(diff_matx, path_save)

    # Return if parsed on function call
    return diff_matx


def make_heatmap_kwargs():
    """
    Generates key word arguments for the distance difference map (rendered as a 2D
    histogram). Used several times throughout scripts so changes here will propagate.
    """

    dist_diff_cmap = colormaps["viridis"].copy()
    dist_diff_cmap.set_bad("grey", alpha=0.5)

    heatmap_kwargs = {
        "vmin": 0,
        "cbar": False,
        "square": True,
        "edgecolors": "black",
        "clip_on": True,
        "cmap": dist_diff_cmap,
    }

    return heatmap_kwargs


def plot_2d_hist(dd_matx, axes, heatmap_kwargs):
    """
    Plots parsed distance difference matrix as a 2D histogram.
    """

    sns.heatmap(dd_matx, ax=axes[0], **heatmap_kwargs)

    # Custom formatting
    axes[0].set_ylabel("Sequence ID (UniProt)", fontweight="demi")
    axes[0].set_xlabel("Sequence ID (UniProt)", fontweight="demi")

    unp_range = range(0, len(dd_matx), 50)
    axes[0].set_xticks(unp_range)
    axes[0].set_yticks(unp_range)
    axes[0].set_xticklabels(unp_range)
    axes[0].set_yticklabels(unp_range)


def add_colour_bar(fig, axes):
    """
    Adds colourbar to a figure, based on the heatmap kwargs stored to the figure from
    the plot_2d_hist() function. Figure must have been created with two axes in one row.
    """

    fig.colorbar(axes[0].collections[0], cax=axes[1])


def get_conformer_id(df, pdb, chain):
    """
    Used to extract the conformer ID from a Pandas DataFrame, given a PDB and chain ID.
    DataFrame must have columns: CHAIN_ID and PDBe_ID.

    Used by format_title() to generate the histogram's title.
    """

    conformer_id = df["CONFORMER_ID"][
        (df["PDBe_ID"] == pdb) & (df["CHAIN_ID"] == chain)
    ].iloc[0]

    return conformer_id


def format_title(unp, df, pdbe1, pdbe2, chain1, chain2):
    """
    Creates a title for the distance difference histogram.
    """

    conf_id_1 = get_conformer_id(df, pdbe1, chain1)
    conf_id_2 = get_conformer_id(df, pdbe2, chain2)

    # Title
    return f"{unp} ({conf_id_1}-{conf_id_2}): {pdbe1}-{chain1} to {pdbe2}-{chain2}"


# This will be the benchmark version
def format_2d_hist(axes, title):
    """
    Applies aesthetic formatting features to the distance difference plot
    """

    axes[0].set_title(title, fontweight="bold")


def find_largest_distance_from_matxs(
    path_matxs: PosixPath, dd_matxs_fnames: "list[str]"
) -> float:
    """
    Finds the largest distance in all distance difference matrices. Used to set the
    maximum value for the colour bar in the distance difference maps.
    """

    max_distance = 0

    for i in dd_matxs_fnames:
        path_dd_matx = path_matxs.joinpath(i)
        dd_matx = io_utils.load_matrix_from_tri_upper(path_dd_matx)

        # Update max distance if larger than current max
        max_distance = linear_algebra_utils.find_max(dd_matx, starting_max=max_distance)

    return max_distance


def make_dd_maps(
    path_matxs: PosixPath,
    path_save_maps: PosixPath,
    force: bool = False,
) -> None:
    """
    Plot all unique distance difference matrices as 2D histograms (heatmaps). Must
    have the distance difference matrices using the cluster() method on a
    ClusterConformers() object first.

    :param path_save: Path to save CA distance-diiference heatmaps.
    :type path_save: pathlib.Path
    """

    # Load in the distance difference matrices to find the maximum distance
    dd_matxs_fnames = io_utils.get_fnames(path_matxs)

    # Find the maximum distance in all distance difference matrices
    logger.info("Finding maximum distance")
    max_distance = find_largest_distance_from_matxs(path_matxs, dd_matxs_fnames)
    logger.info(f"Maximum distance: {max_distance}")

    # Assign max distance to heatmap kwargs
    heatmap_kwargs = make_heatmap_kwargs()
    heatmap_kwargs["vmax"] = max_distance

    # Make directory to save maps
    path_save_maps.mkdir(parents=True, exist_ok=True)

    # Plot distance difference maps
    logger.info("Rendering distance difference maps...")
    for dd_matx_fname in dd_matxs_fnames:

        # Only render if not already rendered or force=True
        if not (
            path_save_maps.joinpath(f"{dd_matx_fname[:-4]}.png").exists() and not force
        ):

            # Initialise plot area
            fig, axes = plt.subplots(
                # figsize=(10,3),
                ncols=2,  # matrix | colour_bar
                nrows=1,
                gridspec_kw={"width_ratios": [4, 0.2]},
                tight_layout=True,
            )

            # Plot heatmap
            plot_2d_hist(
                io_utils.load_matrix_from_tri_upper(path_matxs.joinpath(dd_matx_fname)),
                axes,
                heatmap_kwargs,
            )
            logger.debug(f"Rendered {dd_matx_fname}")

            # Format plot
            add_colour_bar(fig, axes)
            format_2d_hist(axes, " ".join(dd_matx_fname.split("_"))[:-4])
            logger.debug(f"Formatted {dd_matx_fname}")

            # Save plot
            io_utils.save_figure(
                path_save_maps,
                save_fname=f"{dd_matx_fname[:-4]}",
                png=True,
                svg=False,  # Can be very expensive
            )
            logger.debug(f"Saved {dd_matx_fname}")

            # Close plot
            plt.close(fig)

        else:
            # Already rendered
            logger.debug(f"Skipping {dd_matx_fname}")

    logger.info("Rendering done.")
