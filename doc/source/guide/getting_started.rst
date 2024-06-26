================
Getting Started
================



Installation
==============

You can install the latest stable version of the package from the PDBe's PyPi page
using:

.. code:: bash

   pip install protein-cluster-conformers

**or** ...

from the source code using:

.. code:: bash

   git clone https://github.com/PDBeurope/protein-cluster-conformers.git

   cd protein-cluster-conformers

   pip install .

The package is written for >=Python3.10.



Quick Start
==============

Once installed, import the package into your code with:

.. code:: python

   import cluster_conformers


``cluster_conformers`` contains many useful tools to help you with downloading,
clustering and analysing your results. The most important object for clustering can be
imported via:

.. code:: python

   from cluster_conformers.cluster_monomers import ClusterConformations

Initialise an instance of the object with the paths to your PDB files. They should be
parsed in as a dictionary like so:

.. code:: python

   input_dictionary = {
       "updated_mmcifs/P46060/5d2m_updated.cif": ["C", "F"],
       "updated_mmcifs/P46060/2grn_updated.cif": ["B"],
   }

   # Create clustering object
   conformations = ClusterConformations(unp="P46060", mmcifs_and_chains=input_dictionary)

The ``input_dictionary`` is a dictionary where the keys are the paths to the PDB files
and the values are a list of the ``structut_asym_id`` chain IDs you would like to
include in the clustering. Only the ``*_updated.cif`` are accepted because residue
correspondence is specified in this file type.

The ``unp`` is the UniProt accession code for the protein you are clustering. ``A12345``
can be parsed if you are not clustering a known UniProt protein.

Next, generate the C-alpha distance matrices with:

.. code:: python

   conformations.ca_distance(path_save_ca_matx="ca_distance_matrices")

And then cluster the conformations with:

.. code:: python

   conformations.cluster(
       path_save_dd_matx="distance_difference_matrices",
       path_save_cluster_results="cluster_results",
   )

The ``find_conformers.py`` wrapper, demonstating the steps above, is provided in the
root of the repository. An instructional notebook (``tutorials/instructions.ipynb``) is
also included in the same GitHub repository, going into more detail on a specific
worked example.



Visualising Results
====================

The cluster results can be visualised in dendrogram form with:

.. code:: python

   from cluster_conformers.cluster_monomers import render_dendrogram

   # Render dendrogram
   render_dendrogram(
       unp="A12345",
       # Generated by ClusterConformations.cluster()
       path_results="path/to/cluster_results",
       # Path to save dendrogram
       path_save="/path/to/dendrogram",
       png=True,  # Save as png
       svg=False,  # (and/or) Save as svg
   )

You can also view the chain-chain comparisons made by the distance difference matrices
by generating heatmaps of the Angstrom differences with:

.. code:: python

   from cluster_conformers.distance_differences import make_dd_maps

   make_dd_maps(
       # Should be same as used in ClusterConformations.cluster()
       path_matxs="path/to/distance_difference_matrices",
       path_save="path/to/dd_maps",
       force=False,
   )


API Reference
==================

The same information can be accessed from the `API documentatation` section as well.

.. toctree::
   :maxdepth: 10

   ../cluster_conformers
