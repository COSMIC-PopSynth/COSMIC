.. _install:

############
Installation
############

=====
Conda
=====
Since COSMIC requires compilation of Fortran code, we strongly recommend the use of virtual environments like Anaconda to mitigate any potential compiling issues with other Fortran code bases.

MacOS
-----
.. note::

    IMPORTANT NOTE: On OS X 10.16 (Big Sur), you may need to supply the following extra linker flag 

    ``LDFLAGS="-L /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib" pip install cosmic-popsynth``

.. code-block:: bash

    conda create -n cosmic gfortran_osx-64 numpy h5py python=3.7
    source activate cosmic
    pip install cosmic-popsynth


Unix
----
.. code-block:: bash

    conda create --name cosmic python=3.7 numpy h5py
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
