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

   conda create -n cosmic gcc numpy h5py python=3.7
   source activate cosmic
   pip install cosmic-popsynth

Unix
----
.. code-block:: bash

   conda create --name cosmic python=3.6
   source activate cosmic
   conda install numpy
   pip install cosmic-popsynth
