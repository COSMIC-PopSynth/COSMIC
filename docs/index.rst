.. cosmic documentation master file, created by
   sphinx-quickstart on Thu Apr 21 14:05:08 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##################################
Welcome to cosmic's documentation!
##################################
cosmic is a rapid binary population synthesis suite with a special purpose of generating realistic Milky Way compact binary populations. 

cosmic currently implements binary evolutionary processes using BSE (`Hurley et al. (2002) <http://adsabs.harvard.edu/abs/2002MNRAS.329..897H>`_). Several modifications have been applied to BSE to account for recent updates to binary evolution especially important to compact binary formation (e.g. metallicity-dependent stellar winds or black hole natal kick strengths). For a detailed discussion of these modifications, see Breivik et al. (2018), `Rodriguez et al (2016) <http://adsabs.harvard.edu/abs/2016PhRvD..93h4029R>`_, and `Dominik et al (2013) <http://adsabs.harvard.edu/abs/2013ApJ...779...72D>`_.

*****************
Installing cosmic
*****************

The easiest method to install cosmic is using `pip <https://pip.pypa.io/en/stable/>`_ directly from the `GitHub repository <https://github.com/COSMIC-PopSynth/COSMIC.git>`_:

.. code-block:: bash

   $ pip install git+https://github.com/COSMIC-PopSynth/COSMIC.git


For more details on how to use cosmic, see :ref:`examples` and :ref:`runpop`.


*****************
Table of Contents
*****************

.. toctree::
   :maxdepth: 4

   examples/index
   runpop/index

*********************
Package documentation
*********************

Please consult these pages for more details on using cosmic:

.. toctree::
   :maxdepth: 1

   command-line/index

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
