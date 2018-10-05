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
   pip install cosmic-popsynth

Unix
----
.. code-block:: bash

   conda create --name cosmic python=3.6
   source activate cosmic
   conda install numpy
   pip install cosmic-popsynth
