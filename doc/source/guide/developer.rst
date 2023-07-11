================
Developer Guide
================

To contribute to the development of this project, clone the repository and install both
the base and developer dependencies:

.. code:: bash

    git clone https://github.com/PDBeurope/protein-cluster-conformers.git

    cd protein-cluster-conformers

    pip install -r requirements.txt

    pip install -r dev-requirements.txt


If you would like to make edits to the source code and have those changes immediately
update, install the package in editable mode:

.. code:: bash

   pip install -e .

Next, install the pre-commit hooks:

.. code:: bash

    pre-commit install

and use:

.. code:: bash

    pre-commit run --all-files

to run the checks on all files before committing changes.

To run the tests, use:

.. code:: bash

    pytest --cov=cluster_conformers --cov-report=html -v
