#!/usr/bin/python3

"""
This script is a wrapper to run UPGMA agglomerative clustering on a set of monomeric
chains, parsed in by the user. These chains should all have the same UniProt accession
and have at least partial sequence overlap. Although the latter condition is not
essential, a lack of sufficient overlap will not provide informative clustering results.

Functions defined here handle argument parsing, such as setting up the input to the
cluster_monomers.ClusterConformations() class instance and make decisions on which
methods to execute.
"""

# Third party imports
import argparse
from pathlib import Path, PosixPath
import sys

# Custom imports
import cluster_conformers.cluster_monomers as cluster_monomers
from cluster_conformers.utils import logging_utils
from cluster_conformers.distance_differences import make_dd_maps


def extract_image_format(image_args: str):
    """
    Image formats are parsed in from the command line by the user, along with paths to
    their save location. This function separates the path from specified image format(s)
    and returns them as easily-handleable booleans.
    """

    png_bool = False
    svg_bool = False

    for arg in image_args:
        if arg == "png":
            png_bool = True
        if arg == "svg":
            svg_bool = True
    path_image = Path(image_args[0])

    return path_image, png_bool, svg_bool


def extract_structure_format(args_mmcif):
    """
    Takes arguments collected from the create_parser() function (below) and creates a
    dictionary acceptable by the cluster_monomers.ClusterConformations() object.

    Output example for 'structures' object:
        "/path/to/updated/mmcif/1atp_updated.cif" : ['A', 'B'],
        "/path/to/updated/mmcif/2adp_updated.cif" : ['C', 'D', 'E'],
        ...
        "/path/to/updated/mmcif/9amp_updated.cif" : ['A', 'B', ... 'Z']

    """

    # Add parsed list of mmCIFs to dictionary
    structures = {  # str : list
        # "/path/1atp_updated.cif" : ['A', 'B'],
        # "/path/2adp_updated.cif" : ['C', 'D', 'E'],
        # ...
        # "/path/9amp_updated.cif" : ['A', 'B', ... 'Z']
    }

    if args_mmcif:
        try:
            for i in args_mmcif:
                structures[i[0]] = i[1:]
        except Exception:
            raise IndexError("Must parse in chain ID(s)")
    else:
        raise NameError("Must parse in path to one or more mmCIF file(s)")

    return structures


def create_parser(input_args=None):
    """
    Collects command-line arguments from the user and parses them into a dictionary,
    ready for feeding into cluster_monomers.ClusterConformations(), and an 'arguments'
    object, which is used to make decisions on which methods in
    cluster_monomers.ClusterConformations() to run.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", "--verbose", help="Increase verbosity", default=False, action="store_true"
    )

    parser.add_argument(
        "-u", "--uniprot", help="UniProt accession", type=str, required=True
    )

    parser.add_argument(
        "-m",
        "--mmcif",
        nargs="+",
        action="append",
        help="Enter list of paths to mmCIFs that overlap a given UniProt segment",
        # type=pathlib.Path
        required=True,
    )

    parser.add_argument(
        "-s",
        "--path_clusters",
        help="Path to save clustering results",
        type=PosixPath,
    )

    parser.add_argument(
        "-c",
        "--path_ca",
        help="Path to save CA distance matrices",
        type=PosixPath,
        required=True,
    )

    parser.add_argument(
        "-d",
        "--path_dd",
        help="Path to save distance difference matrices",
        type=PosixPath,
        default=None,
    )

    parser.add_argument(
        "-g",
        "--path_dendrogram",
        nargs="+",
        help="Path to save dendrogram of clustering results",
        type=str,
    )

    parser.add_argument(
        "-w",
        "--path_swarm",
        help="Path to save swarm plot of scores",
        nargs="+",
        type=str,
    )

    parser.add_argument(
        "-o",
        "--path_histogram",
        help="Path to save histograms of distance difference maps",
        type=PosixPath,
    )

    parser.add_argument(
        "-a",
        "--path_alpha_fold",
        help="Path to save AlphaFold Database structure",
        type=PosixPath,
        default=None,
    )

    parser.add_argument(
        "-n", "--nproc", help="Max number of threads to utilise", type=int, default=1
    )

    parser.add_argument(
        "-f",
        "--force",
        help="Force overwrite of existing matrix files",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "-i",
        "--updated_entries",
        help="List of updated entries, matrices will be regenerated",
        # action="append",
        nargs="+",
        type=str,
        default=None,
    )

    args = parser.parse_args(input_args)

    # Add parsed list of mmCIFs to dictionary
    structures = extract_structure_format(args.mmcif)

    return args, structures


def main():
    """
    Wrapper to run cluster_monomers.ClusterConformations() on a set of parsed mmCIFs.
    """
    args, structures = create_parser(sys.argv[1:])

    # Initialise logger
    logging_utils.init_logger(verbose=args.verbose)

    # Create object for clustering
    unp_cluster = cluster_monomers.ClusterConformations(
        unp=args.uniprot,
        mmcifs_and_chains=structures,
        path_save_alphafold=args.path_alpha_fold,
        nproc=args.nproc,
        force=args.force,
    )

    # Remove any existing matrices for updated entries
    if args.updated_entries:
        unp_cluster.remove_entry_matxs(
            pdb_ids=args.updated_entries,
            path_ca=args.path_ca,
            path_dd=args.path_dd,
        )

    # Generate CA distance matrices and save
    unp_cluster.ca_distance(args.path_ca)

    # Perform agglomerative clustering and save results
    if args.path_clusters:

        unp_cluster.cluster(
            path_save_dd_matx=args.path_dd,
            path_save_cluster_results=args.path_clusters,
        )

    elif bool(args.path_clusters):
        raise NameError(
            "Must parse both path to save distance difference matrices "
            "and clustering results. Use -d </path/to/distance/difference/matrices/> "
            "-s </path/to/save/clustering/results/>"
        )
    else:
        pass

    # Render and save distance difference maps
    if args.path_histogram:

        make_dd_maps(
            path_matxs=args.path_dd,
            path_save_maps=args.path_histogram,
            force=True,
        )

    # Parsing in options for saving dendrogram
    if args.path_dendrogram:

        path_save, png_bool, svg_bool = extract_image_format(args.path_dendrogram)

        cluster_monomers.render_dendrogram(
            unp=args.uniprot,
            path_results=path_save,
            path_save=path_save,
            png=png_bool,
            svg=svg_bool,
        )


if __name__ == "__main__":

    main()
