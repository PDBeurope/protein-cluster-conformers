#!/usr/bin/python3

"""
Wrapper for running the conformational state clustering on the benchmark dataset.
"""

# Third party imports
import argparse
import logging
import pathlib
import sys

import pandas as pd
from matplotlib import pyplot as plt

# Custom imports
from cluster_conformers import cluster_monomers
from cluster_conformers.utils import download_utils

# Global variables
PATH_BASE = pathlib.Path("benchmark_data")

PATH_BENCHMARK = PATH_BASE.joinpath("benchmark_monomeric_open_closed_conformers.csv")

PATH_MMCIFS = PATH_BASE.joinpath("updated_mmcifs")

PATH_SAVE_CA = PATH_BASE.joinpath("all_uniprots", "ca_distances")

PATH_SAVE_DD_MATXS = PATH_BASE.joinpath(
    "all_uniprots", "distance_differences", "all_conformers"
)

PATH_SAVE_CLUSTER_RESULTS = PATH_BASE.joinpath("all_uniprots", "clustering_results")

PATH_SAVE_DD_MAPS = PATH_BASE.joinpath("all_uniprots", "distance_difference_maps")

PATH_ALPHAFOLD_MMCIFS = PATH_BASE.joinpath("alphafold_mmcifs")


def benchmark_parser(args):
    """
    Collect and parse command line arguements for this script.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s",
        "--single-acc",
        default=None,
        type=str,
        help="UniProt accession code. Only cluster this code from the benchmark "
        " dataset. ",
    )

    return parser.parse_args(args)


def benchmark_cluster(benchmark_df, unp):
    """
    Initialises clustering object and runs clustering methods for all chains per UniProt
    segment described in the benchmark dataset CSV file.
    """

    # Prepare input for ClusterConformations() based on info in benchmark CSV
    structures_dict = {}
    df_unp = benchmark_df[benchmark_df["UNP_ACC"] == unp]
    for pdb in df_unp["PDBe_ID"].unique():
        chains = list(df_unp["label_asym_id"][df_unp["PDBe_ID"] == pdb])
        path_mmcif = PATH_MMCIFS.joinpath(f"{pdb}_updated.cif")
        structures_dict[str(path_mmcif)] = chains

    # Prepare chains in clustering object
    unp_cluster = cluster_monomers.ClusterConformations(
        unp=unp,
        mmcifs_and_chains=structures_dict,
        path_save_alphafold=PATH_ALPHAFOLD_MMCIFS,  # Include AFDB structures
    )

    # Generate CA matrices
    unp_cluster.ca_distance(PATH_SAVE_CA)

    # Cluster
    unp_cluster.cluster(
        path_save_dd_matx=PATH_SAVE_DD_MATXS,
        path_save_cluster_results=PATH_SAVE_CLUSTER_RESULTS,
    )

    # Post-processing
    # unp_cluster.select_representatives()

    """
    To render CA distance difference maps, dendrograms and swarm plots, uncomment
    the lines between:

    here ....
    """

    # unp_cluster.make_dd_maps(PATH_SAVE_DD_MAPS)

    png_bool = True
    svg_bool = True

    unp_cluster.make_dendrogram(
        PATH_SAVE_CLUSTER_RESULTS.joinpath("dendrograms"), png=png_bool, svg=svg_bool
    )

    unp_cluster.make_swarmplot(
        PATH_SAVE_CLUSTER_RESULTS.joinpath("swarm_plots"), png=png_bool, svg=svg_bool
    )

    """
    .... and here
    """

    # Flush objects
    del unp_cluster
    plt.clf()


def main(unp=None):
    """
    Wrapper function for clustering

    E.g. :
     - `python cluster_benchmark.py` for whole dataset
     - `python cluster_benchmark.py -s <uniprot>` for specific segment in dataset
    """
    benchmark_df = pd.read_csv(PATH_BENCHMARK)

    download_utils.fetch_benchmark_mmcifs(PATH_BENCHMARK, PATH_MMCIFS)

    if unp:
        benchmark_cluster(benchmark_df, unp)
    else:
        for unp in benchmark_df["UNP_ACC"].unique():
            benchmark_cluster(benchmark_df, unp)


if __name__ == "__main__":
    """
    Run from the command line
    """

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
        datefmt="%m-%d %H:%M",
    )

    args = benchmark_parser(sys.argv[1:])

    main(args.single_acc)
