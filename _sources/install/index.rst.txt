.. _install:

############
Installation
############

===========
Prequisites
===========
Since COSMIC requires compilation of Fortran code, you'll need a gfortran installation. Several options exist for installing gfortran including through `homebrew <https://brew.sh/>`_ or `from source <https://gcc.gnu.org/wiki/GFortran>`_. If you have a gfortran installation that works with other code bases, chances are it will work with COSMIC too!

MacOS
-----
.. note::

    The largest hurdle for installation on MacOS is keeping your gfortran installation up to date with the linking libraries in Mac's commandlinetools. When in doubt, reinstall your gfortran library then try reinstalling COSMIC.

.. code-block:: bash

    conda create -n cosmic numpy h5py python=3.10
    source activate cosmic
    pip install cosmic-popsynth


Unix
----
.. code-block:: bash

    conda create --name cosmic python=3.10 numpy h5py
    source activate cosmic
    pip install cosmic-popsynth

Installation Notes/FAQ
----------------------

.. note::

    USING IPYTHON OR JUPYTER-NOTEBOOKS WITH COSMIC ENVIRONMENT

    Please note that using the global instance of the conda jupyter-notebook
    or ipython will most likely fail when trying to use COSMIC.
    PLEASE explicitly install both into the COSMIC environment with either

    ``conda install jupyter ipython``

    ``pip install jupyter ipython``

.. note::

    USING COMSIC WHEN BUILT FROM SOURCE

    If you want import the fortran wrapped library
    from the GITHUB folder itself, i.e.

    ``from cosmic import _evolvebin``

    then you must build the extension locally

    ``python setup.py build_ext --inplace``
