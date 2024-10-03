.. _install:

############
Installation
############

First install ``gfortran``...
=============================
Since COSMIC requires compilation of Fortran code, you'll need a gfortran installation. Several options exist for installing gfortran including through `homebrew <https://brew.sh/>`_ or `from source <https://gcc.gnu.org/wiki/GFortran>`_. If you have a gfortran installation that works with other code bases, chances are it will work with COSMIC too!

.. tab-set::
    :sync-group: os

    .. tab-item:: MacOS
        :sync: mac

        .. code-block:: bash

            brew install gcc

    .. tab-item:: Unix
        :sync: unix

        .. code-block:: bash

            sudo apt-get install gfortran

        .. note::

            If you're using a different package manager, you can search for gfortran in the package manager's search bar.


...then install ``COSMIC``
==========================

We recommend following the code below to create a conda environment for COSMIC. This will ensure that all dependencies are installed correctly. If you don't have conda installed, you can download it from `here <https://docs.conda.io/en/latest/miniconda.html>`_.

.. tab-set::
    :sync-group: os

    .. tab-item:: MacOS
        :sync: mac

        .. code-block:: bash

            conda create -n cosmic numpy h5py python=3.10
            source activate cosmic
            pip install cosmic-popsynth

        .. note::

            The largest hurdle for installation on MacOS is keeping your gfortran installation up to date with the linking libraries in Mac's commandlinetools. When in doubt, reinstall your gfortran library then try reinstalling COSMIC.

    .. tab-item:: Unix
        :sync: unix

        .. code-block:: bash

            conda create --name cosmic python=3.10 numpy h5py
            source activate cosmic
            pip install cosmic-popsynth


Using IPython and Jupyter with COSMIC
-------------------------------------

Please note that using the global instance of the conda jupyter-notebook
or ipython will most likely fail when trying to use COSMIC.
PLEASE explicitly install both into the COSMIC environment:

.. tab-set::

    .. tab-item:: conda

        .. code-block:: bash

            conda install jupyter ipython

    .. tab-item:: pip

        .. code-block:: bash

            pip install jupyter ipython


