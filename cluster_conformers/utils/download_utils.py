"""
This script includes functions to access AlphaFold2 structures for inclusion in monomer
superimposision via GESAMT and conformational clustering via CA distance-based score.

Functions here should be imported into superpose.py and run only upon request by user.
Modification of the clustering output should also
"""

from logging import getLogger
from pathlib import PosixPath
import requests

# Thir-party imports
from pandas import read_csv
from requests import get

# Log information
logger = getLogger(__name__)


def fetch_benchmark_mmcifs(path_benchmark_df: PosixPath, path_save: PosixPath) -> None:
    """
    Downloads all updated mmCIF files located inside the parsed benchmark dataset file.
    Should not be run frequently.

    :param path_benchmark_df: Path to benchmark dataset
    :type path_benchmark_df: PosixPath
    :param path_save: Path to save downloaded mmCIFs
    :type path_save: PosixPath
    """
    # Read in list of benchmark structures
    benchmark_df = read_csv(path_benchmark_df)
    pdbe_ids = benchmark_df["PDBe_ID"].unique()

    # Download and save
    for pdbe in pdbe_ids:
        fetch_updated_mmcif(pdbe, path_save)


def fetch_updated_mmcif(pdb_code: str, path_save: PosixPath) -> None:
    """
    Function to download updated mmCIF file from the PDBe

    :param pdb_code: PDB accession code
    :type pdb_code: str
    :param path_save: Location to save downloaded updated mmCIF
    :type path_save: Path
    """
    # Base URL to download mmCIFs
    url = "https://www.ebi.ac.uk/pdbe/entry-files/download/"

    mmcif_file_name = f"{pdb_code}_updated.cif"
    download_link = url + mmcif_file_name

    # Download
    logger.info("Downloading", mmcif_file_name)
    response = get(download_link, allow_redirects=True)
    save_to = path_save.joinpath(mmcif_file_name)
    open(save_to, "wb").write(response.content)


def download_alphafold_mmcif(uniprot: str, path_save: PosixPath) -> PosixPath:
    """
    Downloads AlphaFold Database structure to specified location.

    :param uniprot: UniProt accession
    :type uniprot: str
    :param path_save: Path to save AlphaFold structure
    :type path_save: PosixPath
    :raises ConnectionError: Query-related download error
    :return: Path to saved AlphaFold mmCIF
    :rtype: PosixPath
    """

    # Construct URL
    afdb_file_name = f"AF-{uniprot}-F1-model_v3.cif"
    afdb_url = f"https://alphafold.ebi.ac.uk/files/{afdb_file_name}"

    # Retrieve file AlphaFold Database connection
    afdb_request = get(afdb_url, stream=True)

    # Try to save file to system
    if afdb_request.status_code == 200:
        # Proceed if connection made successfully
        logger.info(f"AlphaFold file for {uniprot} retrieved successfully")

        # Make save dir
        path_save = path_save.joinpath(uniprot)
        path_save.mkdir(parents=True, exist_ok=True)

        # Save
        path_save_file = path_save.joinpath(afdb_file_name)
        with open(path_save_file, "wb") as mmcif:
            mmcif.write(afdb_request.content)
            logger.info(f"AlphaFold structure for {uniprot} saved to {path_save_file}")

    elif afdb_request.status_code == 404:
        # File was not available
        logger.warning(f"AlphaFold structure for {uniprot} not available")
        return None

    else:
        # Connection error
        logger.error(
            "AlphaFold fold could not be retrieved ( API request error "
            f"{afdb_request.status_code} )"
        )
        raise ConnectionError

    # Return file save location
    return path_save_file


def download_cluster_results(uniprot: "str", save_dir: "PosixPath") -> None:
    """
    Function to download cluster results from the PDBe's GraphAPI.

    :param uniprot:
    :type uniprot: str
    :param save_dir: _description_
    :type save_dir: str|PosixPath
    :return: _description_
    :rtype: _type_
    """

    base_api_uri = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/"
    query = f"{base_api_uri}{uniprot}"

    logger.debug(f"Querying Graph API with UniProt={query}")
    request = requests.get(query, allow_redirects=True)

    if request.status_code != 200:
        logger.debug(f"Request for {uniprot} failed with status: {request.status_code}")
        return None
    else:
        save_dir = save_dir.joinpath(f"{uniprot}_gesamt_output.csv")
        open(save_dir, "wb").write(request.content)
        logger.debug(f"Saved {uniprot} to {save_dir}")
