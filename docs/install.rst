Installing TEFingerprint
========================

:Date: 2017-07-31
:Authors: Tim Millar
:Organization: Plant & Food Research


Requirements
------------

Requires Python 3.4.0 or higher (Python 2 not supported)

Python library requirements (from requirements.txt):

- numpy>=1.10.0
- pysam>=0.9.0
- biopython>=1.65

External software required for extracting informative reads:

- samtools
- bwa (bwa-mem)

Development Environment
-----------------------

A development environment can be built using the conda tool available in
the `Anaconda Python
distribution <https://www.continuum.io/downloads>`__. An `environment
file <http://conda.pydata.org/docs/using/envs.html#share-an-environment>`__
"environment.yml" is available in the root directory of this project.

::

    conda env create -f environment.yml

This will create a `virtual
environment <http://conda.pydata.org/docs/using/envs.html>`__ called
``tefingerprint``.

Activate the development environment:

::

    conda activate tefingerprint

Manual Installation
-------------------------

- Uninstall any previous version:

::

    pip uninstall tefingerprint

-  Navigate to the projects root directory then build and install:

   -  (note the ``*`` wildcard should be replaced with the current
      version number to install)

::

    python setup.py sdist
    pip install dist/tefingerprint-*.tar.gz

-  Alternatively, The project can be installed in editable mode from the
   projects root directory:

::

    pip install -e ./

GitHub Installation
-------------------

Installing a specific version e.g. 0.3.0

::

    pip install git://github.com/PlantandFoodResearch/TEFingerprint.git@v0.3.0

If using two factor authentication with ssh

::

    pip install git+ssh://git@github.com/PlantandFoodResearch/TEFingerprint.git@v0.3.0

Testing
-------

Unit testing using `pytest <http://doc.pytest.org/en/latest/>`__. To run
tests from the projects root directory:

::

    python -m pytest -v ./
