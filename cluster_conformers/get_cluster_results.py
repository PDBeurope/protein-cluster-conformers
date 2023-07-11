from pathlib import Path
import sys
import pandas as pd
import logging

from .utils.download_utils import download_cluster_results
from .utils.logging_utils import init_logger

logger = logging.getLogger(__name__)


if __name__ == "__main__":

    sys.argv[1:]

    init_logger(verbose=True)

    path_home = Path.home()

    path_bm_df = path_home.joinpath(
        "EMBL-EBI",
        "funclan_work",
        "static-conformer-dataset",
        "open_closed_mancur_data",
        "benchmark_monomeric_open_closed_conformers.csv",
    )

    path_save = path_home.joinpath(
        "EMBL-EBI", "funclan_work", "static-conformer-dataset", "gesamt_outputs"
    )

    bm_df = pd.read_csv(path_bm_df)

    count = 0

    df_size = len(bm_df["UNP_ACC"].unique())

    ftp_url = "http://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition"

    for unp in bm_df["UNP_ACC"].unique():

        for i in range(0, 6):
            save_dir = path_save.joinpath(f"{unp}_segment{i}_gesamt_output.txt")
            tmp_ftp_url = f"{ftp_url}/{unp[0]}/{unp}/segment{i}/Gesamt_Output.out"

            download_cluster_results(tmp_ftp_url, save_dir)

        count += 1
