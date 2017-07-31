Installing TEFingerprint
========================

:Date: 2017-07-31
:Authors: Tim Millar
:Organization: Plant & Food Research

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

Dependencies
------------

The current development version supports Python 3.5 or higher and pysam
0.9 or higher.

Building and Installation
-------------------------

-  Activate the development environment:

::

    source activate tefingerprint

-  Uninstall any previous version:

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

Testing
-------

Unit testing using `pytest <http://doc.pytest.org/en/latest/>`__. To run
tests from the projects root directory:

::

    python -m pytest -v ./
