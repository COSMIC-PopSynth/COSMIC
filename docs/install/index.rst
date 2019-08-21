.. _install:

############
Installation
############

=====
Conda
=====

MacOS
-----
.. code-block:: bash

    conda create -n cosmic gfortran_osx-64 numpy h5py python=3.7
    source activate cosmic
    pip install --upgrade cosmic-popsynth

.. note::

    IMPORTANT NOTE: On OS X 10.14 (Mojave), you need to install the development header files. These must be installed by hand. To do so, run the command 

    ``open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg``

Unix
----
.. code-block:: bash

    conda create --name cosmic python=3.7 gfortran_linux-64 numpy h5py
    source activate cosmic
    pip install --upgrade cosmic-popsynth

Installation Notes/FAQ
----------------------

.. note::

    USING IPYTHON OR JUPYTER-NOTEBOOKS WITH COSMIC ENVIRONMENT

    Please note that using the global instance of the conda jupyter-notebook
    or ipython will most likely fail when trying to use cosmic.
    PLEASE explciitly install both into the cosmic environment with either

    ``conda install jupyter ipython``

    ``pip install jupyter ipython``

.. note::

    USING COMSIC WHEN BUILT FROM SOURCE

    If you want import the fortran wrapped library
    from the GITHUB folder itself, i.e.

    ``from cosmic import _evolvebin``

    then you must build the extension locally

    ``python setup.py build_ext --inplace``
