import os
import sys

import oracledb
import pandas as pd
import sqlalchemy
from dotenv import find_dotenv, load_dotenv


def get_label_asym_ids(pdb_id, auth_asym_id):

    query = """
        SELECT
            struct_asym_id
        FROM
            sifts_xref_residue
        WHERE
            entry_id = :1 AND
            auth_asym_id = :2
    """

    data = connection.execute(query, (pdb_id, auth_asym_id))

    return str(list(data)[0][0])


def initialise_oracle_database():
    """
    Initialise Oracle with the TNS ADMIN file

    TODO: Use of this function inside of Codon must be validated
    """

    # Load environment variables
    load_dotenv(find_dotenv())

    # Workaround for SqlAlchemy not yet handling python-oracledb
    oracledb.version = "8.3.0"
    sys.modules["cx_Oracle"] = oracledb
    oracledb.defaults.config_dir = os.environ.get("TNS_ADMIN")

    engine_path = "oracle+cx_oracle://{}:{}@{}".format(
        os.environ.get("ORACLE_USER"),
        os.environ.get("ORACLE_PASSWORD"),
        os.environ.get("ORACLE_HOST"),
    )
    oracle_engine = sqlalchemy.create_engine(engine_path)

    return oracle_engine


if __name__ == "__main__":

    # Open benchmark dataset
    df = pd.read_csv("./benchmark_data/benchmark_monomeric_open_closed_conformers.csv")
    print(df.head())
    print(df.columns)

    # Make DB connection
    engine = initialise_oracle_database()
    connection = engine.connect()

    label_asym_id_list = []
    for _, row in df.iterrows():
        # Send query
        label_asym_id = get_label_asym_ids(str(row["PDBe_ID"]), str(row["CHAIN_ID"]))
        print(str(row["PDBe_ID"]), str(row["CHAIN_ID"]), label_asym_id)
        label_asym_id_list.append(label_asym_id)

    print(label_asym_id_list)
    df["label_asym_id"] = label_asym_id_list
    print(df.columns)
    print(df.head())

    df.to_csv(
        "./benchmark_data/benchmark_monomeric_open_closed_conformers_struct_asym_id.csv"
    )
