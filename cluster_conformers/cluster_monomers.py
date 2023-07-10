"""
The class here performs clustering on a set of parsed chains.
"""

# Third party imports
from logging import getLogger
from multiprocessing import Pool
from pathlib import PosixPath

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import use as mpl_use
from pandas import DataFrame

# Custom imports -- peptide_analysis
from . import cluster_chains, distance_differences

# Custom imports -- utils
from .utils import (
    appearance_utils,
    download_utils,
    io_utils,
    linear_algebra_utils,
    parsing_utils,
)

# Global variable
CLUSTERING_CUTOFF_PC = 0.7
logger = getLogger(__name__)
mpl_use("AGG")  # Lighter-weight Matplotlib backend


class ClusterConformations:
    """
    Object to collect clustering information from a parsed list of paths to updated
    mmCIF files.

    The object contains methods to create CA distance matrices, distance difference
    matrices, a CSV file containing cluster assignments for parsed mmCIFs, dendrogram
    of clustered chains, plots (as 2D histograms) of distance difference maps, and swarm
    plots of scores for chain-chain comparisons.

    The object relies on functions used throughout the codebase to perform conformation
    clustering on any given list of mmCIF files.

    TODO:
     - Implement cluster variable measure per residue
     - Add function to plot cluster variability and between-cluster variability features
    """

    # Constructor
    def __init__(
        self,
        unp: str,
        mmcifs_and_chains: "dict[str, list[str] ]",
        path_save_alphafold: PosixPath = None,
        nproc: int = 10,  # max 10 chosen based on benchmarking results
        force: bool = False,
    ) -> None:
        """
        Constructor -- setup object.
        Loads mmCIFs, PDBe IDs, chains IDs, and paths parsed into script.

        :param unp: UniProt accession
        :type unp: str
        :param mmcifs_and_chains: Dictionary of paths to updated mmCIF files (keys) and
            the desired chain (values) as label_asym_ids in a list.
        :type mmcifs_and_chains: dict[str, list[str] ]
        :param path_save_alphafold: Path to save downloaded AlphaFold structre, acting
            as a flag to include it in clustering or not, defaults to None
        :type path_save_alphafold: PosixPath, optional
        :param nproc: Number of processes to use for multiprocessing, defaults to 10
        :type nproc: int, optional
        :param force: Force re-generation of all matrices, defaults to False
        :type force: bool, optional
        """

        self.unp = unp  # UniProt accession

        self.mmcif_paths = {  # Paths to mmCIFs
            # "1atp" : path to mmCIF as str,
            # "2adp" : path to mmCIF as str,
            # ...
        }
        self.chains = {  # Unique lists of desired chains
            # "1atp" : ["A", "B", ...],
            # "2adp" : ["A", "B", ...],
            # ...
        }
        self.chains_all = []  # Redundant list of chain IDs

        # Cannot be set as it requires preservation of order
        # TODO: Consider sorting this list, it may save space and runtime later on
        self.pdbe_ids = []  # Unique list of PDBe names. Ordered

        self.pdbe_chain_ids = []  # String of pdbID_chainID. Ordered, unique

        # Track progress
        logger.info(f"Loading mmCIF files for {self.unp}")

        # Load mmCIF files into object
        for mmcif_path, chains in mmcifs_and_chains.items():

            # Only load mmCIFs
            if (mmcif_path[-3:] == "cif") or (mmcif_path[-6:] == "cif.gz"):

                # Store PDBe ID
                mmcif_fname = mmcif_path.split("/")[-1]
                pdbe_id = mmcif_fname[:4]
                self.pdbe_ids.append(pdbe_id)

                # Store list of chains
                self.chains[pdbe_id] = chains

                self.mmcif_paths[pdbe_id] = mmcif_path
                logger.debug(f"Stored {mmcif_path} to dictionary.")

                # Record order of pdbe-chain pairs in two order-depended lists
                for chain in chains:
                    self.chains_all.append(chain)
                    self.pdbe_chain_ids.append(f"{pdbe_id}_{chain}")
                    logger.debug(f"Stored {pdbe_id}_{chain} to ordered list.")

        # Handle pre-clustering AlphaFold file information
        if path_save_alphafold:

            # Download and save
            afdb_path = download_utils.download_alphafold_mmcif(
                self.unp, path_save_alphafold
            )

            # Only add AlphaFold structure info if not 404 error
            if afdb_path:

                # Rename for clustering compatibility with archive naming scheme
                af_prefix = "afv3"
                af_fname = f"{af_prefix}_updated.cif"
                io_utils.rename_file(afdb_path, af_fname)

                # Load AlphaFold structure
                path_alphafold = afdb_path.parent.joinpath(af_fname)
                afdb_structure = io_utils.load_mmcif(path_alphafold)

                # Store object for clustering in later methods
                self.mmcif_paths[af_prefix] = str(path_alphafold)
                self.pdbe_ids.append(af_prefix)

                # Record chain info for clustering in later methods
                af_chain = "A"
                self.chains[af_prefix] = [af_chain]
                self.chains_all.append(af_chain)

                # Parse AlphaFold structure for extracting chain info for superpose.py
                afdb_mmcif = parsing_utils.parse_mmcif(afdb_structure, af_chain)
                # Storing the start-end UniProt residue indices for protein-superpose
                self.af_unp_range = (
                    afdb_mmcif["unp_res_ids"][0],
                    afdb_mmcif["unp_res_ids"][-1],
                )

        # Number of threads for multiprocessing. Only use 1 if few unique chains
        if len(self.pdbe_chain_ids) > 20:  # 20 chains = 190 dd_matrices
            self.nproc = nproc
        else:
            self.nproc = 1
        logger.debug(f"Using {self.nproc} threads")

        # Regenerate matrices if already present, defaul=False
        self.force = force
        logger.debug(f"Force regeneration of matrices: {self.force}")

        # Conversion to Numpy array for performance in other functions
        self.chains_all = np.asarray(self.chains_all)

        # Convert to Numpy array for performance improvement
        self.pdbe_chain_ids = np.asarray(self.pdbe_chain_ids)

    def remove_entry_matxs(
        self, pdb_ids: "set[str]", path_ca: PosixPath, path_dd: PosixPath
    ):
        """
        Function to remove all CA and distance difference matrices for a given set of
        PDB IDs. This method should be called when existing PDB entries are updated in
        the archive, rendering their CA and distance difference matrices invalid.
        Removing them here will ensure that they are re-generated during the rest of the
        pipeline.
        """

        # Remove CA and DD matrices
        for pdb_id in pdb_ids:
            # Define expression for deletion
            path_ca_matx = path_ca.glob(f"{pdb_id}*.npz")
            path_dd_matx = path_dd.glob(f"{self.unp}*{pdb_id}*.npz")

            # Delete files
            logger.debug(
                "The following PDB entries have been updated: {}".format(pdb_ids)
            )
            for path in [path_ca_matx, path_dd_matx]:
                for file in path:
                    if file.exists():
                        logger.debug(f"Removing {file}")
                        file.unlink()

    def _generate_ca_matx(self, pdbe_chain_id: str) -> "tuple(dict)":
        """
        Method for calculating and saving the CA distance matrix for a given PDB-chain
        ID string.

        :return: tuple of dictionaries to send to fields in the ClusterConformers object
        :rtype: _type_
        """
        # Path to save (or already saved) CA matrix
        path_save_this_ca_matx = self.path_save_base_ca.joinpath(
            f"{pdbe_chain_id}_ca_distance_matrix"
        )

        # Path to save serialist list of
        path_save_this_unp = self.path_save_unps.joinpath(f"{pdbe_chain_id}.pickle")

        # Check if the file already exists
        if not path_save_this_ca_matx.with_suffix(".npz").exists() or self.force:

            pdbe_chain_id_split = pdbe_chain_id.split("_")
            # Extract x, y, z, and UNP index info from mmCIF/chain
            mmcif = io_utils.load_mmcif(self.mmcif_paths[pdbe_chain_id_split[0]])
            xyz_unp_dict = parsing_utils.parse_mmcif(mmcif, pdbe_chain_id_split[1])

            # Serialise the unimputed array of UniProt residue indices
            io_utils.serial_dump(xyz_unp_dict["unp_res_ids"], path_save_this_unp)

            # Make square CA matrix
            xyz_unp_dict = parsing_utils.fill_missing_unps(xyz_unp_dict)
            ca_matx = linear_algebra_utils.generate_ca_matx(
                xyz_unp_dict["cartn_x"],
                xyz_unp_dict["cartn_y"],
                xyz_unp_dict["cartn_z"],
            )

            # Write matrix file if specified
            io_utils.save_compressed_matrix(ca_matx, path_save_this_ca_matx)
            logger.debug(
                f"Generated CA matrix for {pdbe_chain_id}, saved to {path_save_this_ca_matx}"
            )

        else:
            logger.debug(
                f"CA matrix for {pdbe_chain_id} already exists, skipping generation"
            )

        return {  # UniProt residue ID list
            pdbe_chain_id: path_save_this_unp
        }, {  # Path to saved CA matrix
            pdbe_chain_id: path_save_this_ca_matx.with_suffix(".npz")
        }

    def ca_distance(self, path_save: PosixPath = None) -> None:
        """
        Calculates a pairwise CA distance matrix for all parsed mmCIF/chains. Stores
        a reference to matrices in dictionary and saves matrix files if path provided.
        Can be called as a public method to generate CA distance matrices.

        :param path_save: Path to save CA distance matrices, defaults to None
        :type path_save: PosixPath, optional
        """

        self.ca_matxs = {  # CA distance matrices. Ordered
            # "1atp_A" : path (as pathlib.PosixPath) to serilised np.ndarray(...) file,
            # "2adp_B" : path (as pathlib.PosixPath) to serilised np.ndarray(...) file,
            # ...
        }

        self.unp_res_ids = {
            # "1atp_A" : path (as pathlib.PosixPath) to serilised np.array(),
            # "2adp_B" : path (as pathlib.PosixPath) to serilised np.array(),
            # ...
        }

        # Dir to save the raw UniProt residue IDs as 1D np.array()s
        self.path_save_unps = path_save.joinpath("unp_residue_ids")
        self.path_save_unps.mkdir(exist_ok=True)

        self.path_save_base_ca = path_save

        # Track progress
        logger.info("Calculating CA distance matrices...")

        # # Generate CA matrices from loaded mmCIF files if force=True parsed
        pool = Pool(processes=self.nproc)
        try:
            results = pool.map(self._generate_ca_matx, self.pdbe_chain_ids)

            for i in results:
                self.unp_res_ids.update(i[0])
                self.ca_matxs.update(i[1])
        finally:
            pool.close()  # Marks the pool as closed.
            pool.join()  # Waits for workers to exit.

    def _dd_matx_to_score(self, label: str) -> "dict[str: float]":
        """
        Function to either generate de novo or retrieve existing distance-difference
        matrix, calculate its score and then return the result as a key-value paired
        dictionary. The value is the score and the key is a unique label attributed to
        the two chains contributing to the distance difference matrix.

        :return: _description_
        :rtype: _type_
        """

        label_list = label.split("_")
        id_A = f"{label_list[1]}_{label_list[2]}"
        id_B = f"{label_list[-2]}_{label_list[-1]}"

        # Check if distance difference matrix already exists. If it does, skip
        # and load the saved file.
        # TODO: Make sure all data pertaining to updated models are removed

        key = f"{label_list[1]}_{label_list[2]}_{label_list[3]}_{label_list[4]}_{label_list[5]}"
        dd_matx_file = self.path_save_base_dd.joinpath(f"{self.unp}_{key}")

        # Note: "dd_matx" = "distance difference matrix"
        if not dd_matx_file.with_suffix(".npz").exists() or self.force:

            dd_matx = distance_differences.generate_matx_diff(
                io_utils.load_matrix(self.ca_matxs[id_A]),
                io_utils.load_matrix(self.ca_matxs[id_B]),
                path_save=dd_matx_file,
            )
            logger.debug(
                f"Generated distance difference matrix from scratch for {key}."
                f"Saved to {dd_matx_file.with_suffix('.npz')}"
            )

        else:
            # TODO: This needs checking to see if it'll work
            dd_matx = io_utils.load_matrix(dd_matx_file.with_suffix(".npz"))
            logger.debug(
                f"Loaded distance difference matrix for {key} from existing file: "
                f"{dd_matx_file.with_suffix('.npz')}"
            )

        # Calculate score
        score = linear_algebra_utils.calc_score(
            dd_matx,
            io_utils.serial_load(self.unp_res_ids[id_A]),
            io_utils.serial_load(self.unp_res_ids[id_B]),
        )

        del dd_matx

        logger.debug(f"Score for {key} is {score}")

        return {key: score}

    def build_clustering_inputs(
        self,
        path_save_dd_matx: PosixPath,
    ) -> "tuple[ np.ndarray[any, float], np.ndarray[any, str] ]":
        """
        Constructs the requisite data structures for parsing into the
        cluster_agglomerative() function.

        :param path_save_dd_matx: Path to save distance difference matrices.
        :type path_save_dd_matx: PosixPath
        :return: Score matrix and label matrix, corresponding to the chain-chain
            comparisons in the score matrix.
        :rtype: tuple[ ndarray[any, float], ndarray[any, str] ]
        """

        label_matx_linear = []
        label_matx = []  # Matrix of labels corresponding to score_matx.

        self.label_score_reference = {
            # "1atp_A_to_2adp_B" : float of the score
            # "1atp_A_to_3amp_C" : float of the score
            # ...
        }

        # Converted to attribute for access by _dd_matx_to_score()
        self.path_save_base_dd = path_save_dd_matx

        index_A = 0

        # Iter rows
        for id_A in self.pdbe_chain_ids:
            label_row = []  # Structure comparison labels along matrix row

            index_B = 0
            # Fill the row
            for id_B in self.pdbe_chain_ids:

                if index_A < index_B:

                    # Store label for distance difference matrix
                    dd_name = f"{self.unp}_{id_A}_to_{id_B}"
                    label_matx_linear.append(dd_name)

                    # Store label for score
                    score_label = f"{id_A}_to_{id_B}"
                    label_row.append(score_label)

                    logger.debug(
                        f"Added {dd_name} to linear label matrix and {score_label} to "
                        "score label matrix"
                    )

                else:
                    label_row.append(0)
                index_B += 1

            # Label matrices to list
            label_matx.append(label_row)
            index_A += 1

        # Run distance difference matrix calculations in parallel
        pool = Pool(processes=self.nproc)
        try:
            results = pool.map(self._dd_matx_to_score, label_matx_linear)

            for i in results:
                self.label_score_reference.update(i)
        finally:
            pool.close()  # Marks the pool as closed.
            pool.join()  # Waits for workers to exit.

        # Format results from thread pool into score matrix
        score_matx = np.asarray(
            [
                [
                    self.label_score_reference[label] if label != 0 else 0
                    for label in row
                ]
                for row in label_matx
            ]
        )

        return np.maximum(score_matx, score_matx.transpose()), np.asarray(label_matx)

    def cluster(
        self,
        path_save_cluster_results: PosixPath = None,
        path_save_dd_matx: PosixPath = None,
    ) -> None:
        """
        Clusters the chains parsed by sum-based score. Clustering results are stored to
        the object instance under the following accessible attributes:
        - self.cluster_df : pandas.DataFrame of per-chain clustering results in tabular
        format
        - self.cluster_dict : dictionary of per-chain clustering results, predominently
        used by class in `protein-superpose`

        Also writes distance difference matrices and clustering results if paths parsed
        as arguments.

        Along with ca_distance(), this is the only other public method of the class.

        :param path_save_dd_matx: Path to save calculated CA distance-difference
            matrices, defaults to None
        :type path_save_dd_matx: PosixPath, optional
        :param path_save_cluster_results: PosixPath to save clustering results,
            defaults to None
        :type path_save_cluster_results: PosixPath, optional
        """

        logger.info("Generating distance difference matrices...")
        self.score_matx, self.label_matx = self.build_clustering_inputs(
            path_save_dd_matx
        )

        # If all scores==0, place all chains into a single cluster
        if not np.all(self.score_matx == 0):
            logger.info("Non-zero score matrix")

            # Begin clustering
            logger.info("Clustering structures into conformational states...")
            self.model = cluster_chains.cluster_agglomerative(
                self.score_matx, cutoff=CLUSTERING_CUTOFF_PC
            )

            self.linkage_matx = cluster_chains.make_linkage_matx(self.model)

            # Cluster labels
            cluster_labels = self.model.labels_

        else:
            logger.info("All zero score matrix")

            # Cluster labels
            cluster_labels = np.full(self.chains_all.shape[0], 0)
            self.linkage_matx = np.array([0])

        # Store clustering results to object as table
        unps_redundant = np.full(self.chains_all.shape[0], self.unp)
        pdbes = self.pdbe_chain_ids.astype("S4")

        # Output clustering results table
        self.cluster_df = DataFrame(
            {
                "UNP_ACC": unps_redundant,
                "PDBe_ID": pdbes.astype("<U6"),
                "CHAIN_ID": self.chains_all,
                "CONFORMER_ID": cluster_labels,
            }
        )

        # Add clustering results to dictionary for parsing into process_clusters in
        # protein-superpose
        self.cluster_dict = {}
        for conformer_id in np.flip(self.cluster_df["CONFORMER_ID"].unique()):
            # Conformer ID order flipped to preserve compatibility with prot.-superpose
            tmp_df = self.cluster_df[self.cluster_df["CONFORMER_ID"] == conformer_id]

            self.cluster_dict[conformer_id] = list(
                tmp_df["PDBe_ID"] + tmp_df["CHAIN_ID"]
            )

        # Write out clustering results if path specified
        if path_save_cluster_results:

            # Save clustering results
            path_save_all_conf = path_save_cluster_results.joinpath(
                f"{self.unp}_sum_based_clustering_results.csv"
            )
            self.cluster_df.to_csv(path_save_all_conf, index=False)

            # Save score dictionary
            io_utils.serial_dump(
                self.label_score_reference,
                path_save_cluster_results.joinpath(
                    f"{self.unp}_chain_label_scores_dict.pickle"
                ),
            )

            # Save score matrix
            io_utils.save_compressed_matrix(
                matrix=self.score_matx,
                path=path_save_cluster_results.joinpath(f"{self.unp}_score_matrix"),
                label="Score matrix",
            )

            # Save label matrix
            io_utils.save_string_based_matx(
                matrix=self.label_matx,
                path=path_save_cluster_results.joinpath(f"{self.unp}_label_matrix"),
                label="Label matrix",
            )

            # Save linkage matrix in untransformed manor
            io_utils.save_matrix(
                self.linkage_matx,  # Linkage
                path_save_cluster_results.joinpath(
                    f"{self.unp}_linkage_matrix"  # Linkage
                ),
                label="Linkage matrix",
            )

            # Save list of labels corresponding to matrix columns
            io_utils.serial_dump(
                self.pdbe_chain_ids,
                path_save_cluster_results.joinpath(
                    f"{self.unp}_linkage_matx_label_list.pickle"
                ),
            )

        logger.info("Clustering done.")

    def select_representatives(self) -> None:
        """
        Method to select representative structures from the set of predicted clusters.

        TODO:
         - Find runtime improvements, currently very slow
         - Must be either saved in DataFrame format or merged with the exhisting
           clustering CSV
         - Make toggleable.
        """

        self.representatives = {}

        for conformer_id in self.cluster_df["CONFORMER_ID"].unique():
            # Show info to the user
            logger.info(
                "Identifying representative structures for predicted conformer "
                f"{conformer_id}"
            )

            # Select rows in predicted conformation
            chains_in_conformer = self.cluster_df[
                self.cluster_df["CONFORMER_ID"] == conformer_id
            ]

            # Convert to iterable of PDB ID - chain. E.g. 1pdb_A
            pdbe_chain_ids = (
                chains_in_conformer["PDBe_ID"] + "_" + chains_in_conformer["CHAIN_ID"]
            )

            # Create list of CA matrices
            ca_matxs_in_cluster = [self.ca_matxs.get(key) for key in pdbe_chain_ids]

            # Select the smallest max index position to make truncation
            max_index_list = []
            for ca_matx in ca_matxs_in_cluster:
                max_index_list.append(io_utils.serial_load(ca_matx).shape[0])
            truncation_index = min(max_index_list)  # Set truncation value

            # Find all median values
            # Relatively slow -- a performance improvement might be needed here
            logger.info("Generating matrix of element-wise medians")
            # progress = appearance_utils.ProgressBar(truncation_index**2)
            col_index = 0
            median_matrix = []
            while col_index < truncation_index:
                row_index = 0
                median_matrix_row = []

                while row_index < truncation_index:
                    median_at_index = np.median(
                        [
                            io_utils.serial_load(x)[col_index][row_index]
                            for x in ca_matxs_in_cluster
                        ]
                    )
                    median_matrix_row.append(median_at_index)
                    # progress.update()
                    row_index += 1

                median_matrix.append(median_matrix_row)
                col_index += 1

            median_matrix = np.asarray(median_matrix_row)

            # Find number of medians per CA matrix
            median_counts_per_chain = {}
            max_num_medians = 0
            for id in pdbe_chain_ids:

                ca_matx = io_utils.serial_load(self.ca_matxs[id])
                ca_matx_truncated = linear_algebra_utils.matx_trim(
                    ca_matx, truncation_index
                )

                # Find number of medians and store
                num_medians = np.count_nonzero(ca_matx_truncated == median_matrix)
                median_counts_per_chain[id] = num_medians

                # Keep track of the highest number of medians in a given CA matrix
                if num_medians > max_num_medians:
                    max_num_medians = num_medians

            representative_chains = {}
            median_set = False
            for id, median in median_counts_per_chain.items():
                if (median == max_num_medians) and (not median_set):
                    representative_chains[id] = True
                    median_set = True
                else:
                    representative_chains[id] = False

            # Assign representative chains dict to field
            self.representatives.update(representative_chains)

            # Debug
            for key, value in representative_chains.items():
                if value:
                    logger.info(f"Representative for conformer {conformer_id} = {key}")


def render_dendrogram(
    unp: str,
    path_results: PosixPath,
    path_save: PosixPath = None,
    png: bool = False,
    svg: bool = False,
) -> None:
    """
    Plot hierachical dendrogram from clustering results. Must have a linkage matrix and
    ordered labels object already stored. The easiest way to get these is to run the
    ClusterConformations() object first and point to its output folder.

    :param path_save: Path to save rendered dendrogram image.
    :type path_save: PosixPath
    :param png: Save dendrogram image in PNG format, defaults to False
    :type png: bool, optional
    :param svg: Save dendrogram image in SVG format, defaults to False
    :type svg: bool, optional
    """

    # Set matplotlib global formatting
    appearance_utils.init_plot_appearance()

    try:
        linkage_matx = io_utils.load_matrix(
            path_results.joinpath(f"{unp}_linkage_matrix.npz")
        )

        pdbe_chain_ids = io_utils.serial_load(  # TODO: Turn str into global var
            path_results.joinpath(f"{unp}_linkage_matx_label_list.pickle")
        )

    except OSError:
        logger.error(
            "Linkage matrix and/or label list not found. Please run clustering first."
        )
        return

    if not np.array_equal(linkage_matx, np.array([0])):

        fig, ax = plt.subplots(1, 1)

        logger.info("Rendering dendogram")
        cluster_chains.plot_dendrogram(
            unp,
            ax,
            linkage_matx,
            CLUSTERING_CUTOFF_PC,
            labels=pdbe_chain_ids,
            leaf_rotation=90,
        )  # p=3

        io_utils.save_figure(
            path_save,
            save_fname=f"{unp}_agglomerative_dendrogram",
            png=png,
            svg=svg,
        )

        # plt.clf()
        plt.close(fig=fig)

    else:
        logger.info("Single cluster for segment. Not rendering dendrogram.")
