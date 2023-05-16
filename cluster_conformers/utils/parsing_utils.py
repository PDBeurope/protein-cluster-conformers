"""
Functions for data transformations and selections from updated mmCIF files loaded as
Gemmi objects.
"""

# Third party imports
from logging import getLogger

from gemmi import cif
from numpy import NaN

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
            structure_coords[coord_key].insert(index, NaN)

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
