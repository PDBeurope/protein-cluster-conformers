================
User Guide
================

To cluster a set of protein structures, run the ``find_clusters.py`` script, located in
the root of the GitHub repository:

.. code-block:: bash

    python find_conformers.py [-h] [-v] -u UNIPROT -m MMCIF [MMCIF ...]
        [-s PATH_CLUSTERS] -c PATH_CA [-d PATH_DD]
        [-g PATH_DENDROGRAM [PATH_DENDROGRAM ...]]
        [-w PATH_SWARM [PATH_SWARM ...]] [-o PATH_HISTOGRAM]
        [-a PATH_ALPHA_FOLD]


The following parameters can be parsed:

.. code-block:: bash

    required arguments:
        -u UNIPROT, --uniprot UNIPROT
                                UniProt accession
        -m MMCIF [MMCIF ...], --mmcif MMCIF [MMCIF ...]
                                Enter list of paths to mmCIFs that overlap a given UniProt segment
    optional arguments:
        -h, --help            show this help message and exit
        -v, --verbose         Increase verbosity
        -s PATH_CLUSTERS, --path_clusters PATH_CLUSTERS
                                Path to save clustering results
        -c PATH_CA, --path_ca PATH_CA
                                Path to save CA distance matrices
        -d PATH_DD, --path_dd PATH_DD
                                Path to save distance difference matrices
        -g PATH_DENDROGRAM [PATH_DENDROGRAM ...], --path_dendrogram PATH_DENDROGRAM [PATH_DENDROGRAM ...]
                                Path to save dendrogram of clustering results
        -w PATH_SWARM [PATH_SWARM ...], --path_swarm PATH_SWARM [PATH_SWARM ...]
                                Path to save swarm plot of scores
        -o PATH_HISTOGRAM, --path_histogram PATH_HISTOGRAM
                                Path to save histograms of distance difference maps
        -a PATH_ALPHA_FOLD, --path_alpha_fold PATH_ALPHA_FOLD
                                Path to save AlphaFold Database structure



Minimal run example
===================

To only cluster a set of monomeric protein structures that share part or all of the same UniProt sequence, run:

.. code-block:: bash
    python find_clusters.py -u "A12345" \
        -m /path/to/structure_1.cif [chains] \
        -m /path/to/structure_2.cif [chains] \
        ... \
        -m /path/to/structure_N.cif [chains] \
        -s /path/to/save/clustering/results/


The paths to each structure are parsed using the `-m` flag.

Chain IDs (only `struct_asym_id` is currently recognised) should be given as
space-delimited arguments after the path. Parse in multiple structures using
consecutive  `-m` flags. The UniProt accession must be parsed using the `-u` flag.

Further examples of clustering chains from the comamnd line are given in the
repository's ``README.md``.
