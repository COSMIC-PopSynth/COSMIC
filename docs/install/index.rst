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
    
    conda create --name cosmic python=3.6
    source activate cosmic
    conda install gcc numpy
    conda install -c uvcdat gfortran
    pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git

Unix
----
.. code-block:: bash
    
    conda create --name cosmic python=3.6
    source activate cosmic
    conda install numpy
    pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git
