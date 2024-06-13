"""
Functions for data transformations and selections from updated mmCIF files loaded as
Gemmi objects.
"""

# Third party imports
from logging import getLogger
from gemmi import cif
import itertools
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from typing import Generator, Any


logger = getLogger(__name__)


def extract_table(mmcif: cif.Block, search_list: "list[str]") -> cif.Table:
    """
    Produces a Gemmi table based on a list of parsed column names in the _atom_site.
    loop.

    Input: Gemmi block, list of search terms
    Returned: Gemmi table

    :param mmcif: Contents of (updated) mmCIF file
    :type mmcif: gemmi.cif.Block
    :param search_list: Loop terms to extract
    :type search_list: list[str]
    :return: Values for searched loop terms
    :rtype: gemmi.cif.Table
    """

    table = mmcif.find(
        "_atom_site.",
        search_list,
    )

    return table


def fill_missing_unps(
    structure_coords: "dict[str, list[str|float]]",
) -> "dict[str, list[str|float]]":
    """
    Given a list of unique UniProt indices and separate lists of 3D Cartesian coords,
    the function inserts NaN values into each coordinate list where a gap in the UniProt
    sequence is found. Missing UniProt indices are also inserted into the unps list.
    Only gaps and N-terminal truncations are added, nothing is filled in beyond the
    highest UniProt index.

    :param structure_coords: (x, y, z) coordinates and UniProt residue indices.
    :type structure_coords: dict[float, str]
    :return: Same (x, y, z) coordinates and UniProt residue indices as input, with any
        missing coordinates filled as np.NaN
    :rtype: dict[float, str]
    """
    # Insert NaN values into Cartesian coordinates where UniProt index is missing
    complete_res_ids = set(range(1, max(structure_coords["unp_res_ids"]) + 1))

    disjoint_ids = list(complete_res_ids.difference(structure_coords["unp_res_ids"]))
    disjoint_ids.sort()

    # Add missing index to UniProt indices
    # structure_coords["unp_res_ids"] = complete_res_ids

    for missing_id in disjoint_ids:
        index = missing_id - 1

        structure_coords["unp_res_ids"].add(missing_id)

        # Add NaN into coordinates where UniProt is missing
        for coord_key in ("cartn_x", "cartn_y", "cartn_z"):
            structure_coords[coord_key].insert(index, np.NaN)

    return structure_coords


def parse_mmcif(mmcif: cif.Block, chain_id: str) -> "dict[str, list[str|float]]":
    """
    Takes a loaded updated mmCIF as a Gemmi block file and the desired author-specified
    chain ID. Returns a dictionary of four key-value pairs:
        - Cartesian x coords
        - Cartesian y coords
        - Cartesian x coords
        - UniProt residue indices

    Each key ("cartn_x", "cartn_y", "cartn_z", "unp_res_ids") is string type and each
    value is list type.

    :param mmcif: Contents of (updated) mmCIF file.
    :type mmcif: gemmi.cif.Block
    :param chain_id: Structural asymmetry ID (chain identifier).
    :type chain_id: str
    :raises TypeError: Column in updated mmCIF file is missing.
    :return: Cartesian (x, y, z) coordinates and UniProt residue indices.
    :rtype: dict[float, str]
    """

    # To be added to output dict
    cartn_x = []
    cartn_y = []
    cartn_z = []
    unp_res_nums = set()

    mmcif_table = extract_table(
        mmcif,
        [
            "group_PDB",
            "label_atom_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "label_asym_id",
            "pdbx_sifts_xref_db_num",
        ],
    )

    if len(mmcif_table) == 0:
        # Could not parse mmCIF file
        logger.error(
            f"Gemmi block {mmcif} for chain {chain_id} does not contain valid "
            "pdbx_sifts_xref_db_num (UniProt sequence ID) column"
        )

        raise TypeError(
            "Updated mmCIF file is missing SIFTs column headers needed for clustering. "
            "Please repeat clustering step for this UniProt segment once mmCIFs have "
            "been updated."
        )

    # Loop over Table and make checks
    for row in mmcif_table:

        # Check: Atom is CA in protein residue, with occupancy and in correct chain
        if (
            row[0] == "ATOM"
            and row[1] == "CA"  # Peptide atom
            and row[5] == chain_id  # Alpha carbons only
            and row[6] != "?"  # Atom from chain
            and int(row[6])  # UniProt residue ID is defined
            not in unp_res_nums  # Residue not already accounted for
        ):

            # Add Cartesian x,y,z
            cartn_x.append(float(row[2]))
            cartn_y.append(float(row[3]))
            cartn_z.append(float(row[4]))

            # Add UniProt residue number
            unp_res_nums.add(int(row[6]))

    # Add lists of Cartesian x, y, z coords and UniProt indices to dictionary.
    structure_coords = {
        "cartn_x": cartn_x,
        "cartn_y": cartn_y,
        "cartn_z": cartn_z,
        "unp_res_ids": unp_res_nums,
    }

    return structure_coords


class SequenceRanges(object):
    """
    A class handling sequence ranges: do they need to merge, etc
    """

    def __init__(self, range_region: dict):

        # This is the dictionary with each pdbid and chain with the start and end
        self.range_region = range_region

        # Dict returned to the superpose class
        self.segment_pdbchain = {}

    def get_sequence_segments(
        self,
        threshold: float = 0.99,
        linkage: str = "single",
        metric: str = "overlap",
        gap_threshold: int = 25,
    ) -> "dict[int, dict[str, list]]":
        """
        This function will take in the dictionary made by either the query, or by the
        JSON Parsing. Then, it will group the structures based on overlap of sequence
        ranges.

        :param threshold: an optional parameter defining the minimum overlap
            the distance matrix is constructed and normalised to the range [0-1]
            with 0 meaning perfect overlap, and 1 meaning no overlap
            (e.g., the default threshold value of 0.99 would mean
            that an overlap with ~1% of residues of the range would be accepted)
        :param linkage:
            type of clustering for sklearn.cluster.AgglomerativeClustering
            choice of 'single' (default), 'ward', 'complete' and 'average'
        :return: A dictionary with the number of the segment as key, and the list of all
            of the pdb chains in that specific segment as value
        :raises: Nothing

        Example input dictionary.
        Note that each chain could potentially include more than one segment
        range_region =  {
            '5re4A': [(3264, 3569)],
            '5re5A': [(3264, 3569)],
            '6vwwA': [(6452, 6798)],
            '6vwwB': [(6452, 6798)],
            '6m03A': [(3264, 3569)],
            '5rehA': [(3264, 3569)],
            ...
        }
        """

        all_chains = list(self.range_region.keys())
        n_chains = len(all_chains)

        # Create a matrix of overlap distances - each chain against each chain
        overlap_matrix = np.zeros((n_chains, n_chains), dtype=float)

        for i in range(n_chains):
            ranges_i = self.range_region[all_chains[i]]["ranges"]

            for j in range(i + 1, n_chains):
                ranges_j = self.range_region[all_chains[j]]["ranges"]
                overlap_ij = self.compute_overlap(ranges_i, ranges_j, metric=metric)
                overlap_matrix[i, j] = overlap_ij
                overlap_matrix[j, i] = overlap_ij

        # Defines the segments
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=threshold,
            metric="precomputed",
            linkage=linkage,
        ).fit(overlap_matrix)

        # Formats the segments into a dictionaries
        labels = clustering.labels_
        for i in range(n_chains):
            pdbchain = all_chains[i]

            if labels[i] + 1 not in self.segment_pdbchain:
                self.segment_pdbchain[labels[i] + 1] = {
                    "pdbchains": [pdbchain],
                    "ranges": self.range_region[pdbchain]["ranges"],
                }

            else:
                self.segment_pdbchain[labels[i] + 1]["pdbchains"].append(pdbchain)
                self.segment_pdbchain[labels[i] + 1]["ranges"].extend(
                    self.range_region[pdbchain]["ranges"]
                )

        # Merge recorded ranges for each segment
        for segment in self.segment_pdbchain:

            # Remove chains with unmodelled regions
            logger.info(f"Checking segment {segment} for chains with unmodelled gaps")

            for pdbchain in self.segment_pdbchain[segment]["pdbchains"]:
                logger.debug(f"Checking {pdbchain} for gaps")

                # Unmodelled residues
                unmodelled_residues = missing_integers(
                    self.range_region[pdbchain]["residue_ids"]
                )

                # No gaps found
                if len(unmodelled_residues) == 0:
                    logger.debug(f"No gaps found in {pdbchain}")
                    continue

                # Convert set to tuple to allow indexing
                max_gap = find_max_gap((*unmodelled_residues,))
                logger.info(f"Max gap in {pdbchain} is {max_gap}")

                # If max gap is greater than threshold, reject the chain
                if max_gap > gap_threshold:
                    logger.info(
                        f"Rejecting {pdbchain} due to gap ({max_gap}) > {gap_threshold}"
                        f" residues"
                    )
                    self.remove_pdbchain(segment, pdbchain)

            # Merge all the overlapping ranges togther into one range
            self.segment_pdbchain[segment]["ranges"] = self.merge_ranges(
                self.segment_pdbchain[segment]["ranges"]
            )

        return self.segment_pdbchain

    def remove_pdbchain(self, segment: int, pdbchain: str) -> None:
        """
        Remove a pdbchain from the segment dictionary/list records.

        :param segment: Segment number (dict key)
        :type segment: int
        :param pdbchain: pdbid and chain id (e.g. 1atpA).
        :type pdbchain: str
        """
        index_pos = self.segment_pdbchain[segment]["pdbchains"].index(pdbchain)
        self.segment_pdbchain[segment]["pdbchains"].remove(pdbchain)
        self.segment_pdbchain[segment]["ranges"].pop(index_pos)

    def compute_overlap(
        self,
        ranges_i: "list[tuple[int, int]]",
        ranges_j: "list[tuple[int, int]]",
        metric: str = "overlap",
    ) -> float:
        """
        A helper function to compute an overlap coefficient over sequence ranges
        for chains i and j.

        :param ranges_i: a list of tuples [(begin, end), ...] for the first chain
        :type ranges_i: list[tuple[int, int]]
        :param ranges_j: a list of tuples [(begin, end), ...] for the second chain
        :type ranges_j: list[tuple[int, int]]
        :param metric: a choice of overlap (default) or jaccard
            overlap: returns (1 - overlap/min(i,j))
            A perfect overlap of identical ranges would give 0
            A short range fully within a long range would also give 0
            No overlap would give a score of 1
            An overlap of 50 residues between two 100-residue strong
                ranges would give a score of 0.5
            jaccard: returns (1 - overlap/union)
            A perfect overlap of identical ranges would give 0
            A short range (say 50 residues) fully within a long range (of 100 residues)
                would give a score of 0.5
            No overlap would give a score of 1
            An overlap of 50 residues between two 100-residue strong
                ranges would give a score of 0.66
        :type metric: str
        :return: A dissimilarity score controlled by 'metric'
        :rtype: float
        """

        # For quick set operations
        set_i = set()
        set_j = set()

        # Define sets
        for r in ranges_i:
            set_i = set_i.union(range(r[0], r[1] + 1))
        for r in ranges_j:
            set_j = set_j.union(range(r[0], r[1] + 1))

        # Decide on the denominator
        if metric.lower() == "overlap":
            # Shortest set
            denominator = min(len(set_i), len(set_j))

        elif metric.lower() == "jaccard":
            # Union of sets
            denominator = len(set_i | set_j)

        else:
            raise ValueError(
                "Incorrect overlap metric provided. Should be 'overlap' or 'jaccard', "
                f"but {metric} was given."
            )

        # Intersct over denominator (1 - overlap/denominator)
        return 1 - len(set_i & set_j) / denominator

    def merge_ranges(self, ranges: "list[tuple[int, int]]") -> "list[tuple[int, int]]":
        """
        Helper function to merge sets of ranges.
        :param ranges: sets of sequence ranges for two chains i and j
        :type ranges: list[tuple[int, int]]
        :return: list of sequence ranges
        """
        merged_set = set()

        for r in ranges:
            merged_set = merged_set.union(range(r[0], r[1] + 1))

        return list(to_ranges(merged_set))


def to_ranges(iterable: "list|tuple") -> Generator[Any, Any]:
    """
    From
    https://stackoverflow.com/questions/4628333/
    :param iterable: An iterable of integers
    :type iterable: list|tuple
    :return: A generator of ranges
    :rtype: Generator[Any, Any, Any]
    """
    iterable = sorted(set(iterable))

    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]


def find_max_gap(residue_ids: "tuple|list") -> int:
    """
    Identifies the largest contiguous gap in a list of integers (e.g. residue IDs)

    :param residue_ids: Iterable of integers (e.g. tuple residue IDs). It does not need
        to be mutatble, but it must be indexable.
    :type residue_ids: tuple|list
    :return: The largest contiguous gap in the array of integers
    :rtype: int
    """
    max_gap = 0  # Key measure for recording gap size
    gap = 0
    for i, res in enumerate(residue_ids[:-1]):
        # Record incremented gap or if no gap restore to zero
        gap += 1 if res == residue_ids[i + 1] - 1 else 0

        # Record max gap
        if gap > max_gap:
            max_gap = gap

    return max_gap


def missing_integers(sequence: set) -> set:
    """
    Find the missing integers in a sequence

    :param sequence: Sequence of integers
    :type sequence: set
    :return: Set of the missing integers from the sequence
    :rtype: set
    """
    # Residue IDs NOT in common
    complete_set = set(range(min(sequence), max(sequence) + 1))

    return sequence ^ complete_set
