import pathlib

import pandas as pd
import requests


def fetch_benchmark_mmcifs(path_benchmark_df, path_save):

    url = "https://www.ebi.ac.uk/pdbe/entry-files/download/"

    benchmark_df = pd.read_csv(path_benchmark_df)
    pdbe_ids = benchmark_df["PDBe_ID"].unique()

    for pdbe in pdbe_ids:
        mmcif_file_name = f"{pdbe}_updated.cif"
        download_link = url + mmcif_file_name

        print("Downloading", mmcif_file_name)
        r = requests.get(download_link, allow_redirects=True)
        save_to = path_save.joinpath(mmcif_file_name)
        open(save_to, "wb").write(r.content)


if __name__ == "__main__":

    path_home = pathlib.Path.home()
    path_benchmark_df = path_home.joinpath(
        "EMBL-EBI",
        "funclan_work",
        "contact_map_difference",
        "benchmark_data",
        "benchmark_monomeric_open_closed_conformers.csv",
    )

    # Path to save internal CA distances within a peptide
    path_save = path_home.joinpath(
        "EMBL-EBI",
        "funclan_work",
        "contact_map_difference",
        "benchmark_data",
        "updated_mmcifs",
    )
    run_benchmark = True

    if run_benchmark:
        fetch_benchmark_mmcifs(path_benchmark_df, path_save)
