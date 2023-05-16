"""
Uncomment when ready for production
"""

# import argparse

# import gemmi
# import numpy as np

# from cluster_conformers.peptide_analysis.utils import parsing_utils


# def iter_table(gemmi_table):
#     unique_chains = []
#     for row in gemmi_table:
#         if row[0] == "ATOM" and row[1] == "CA":
#             unique_chains += [k for k in row[2] if k not in unique_chains]

#     return unique_chains


# def print_chains(path_to_mmcifs):

#     files = parsing_utils.get_fnames(path=path_to_mmcifs)

#     for i in files:
#         path_file = path_to_mmcifs + i
#         structure = parsing_utils.load_mmcif(path_file)

#         gemmi_table = parsing_utils.extract_table(
#             structure, ["group_PDB", "label_atom_id", "auth_asym_id"]
#         )
#         unique_chains = iter_table(gemmi_table)

#         out_line = path_file
#         for chain in unique_chains:
#             out_line += " " + chain

#         print("-m", out_line, "\\")


# if __name__ == "__main__":

#     parser = argparse.ArgumentParser()

#     parser.add_argument("-d", "--directory")

#     args = parser.parse_args()

#     print_chains(args.directory)
