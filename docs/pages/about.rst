.. _about:

#####
About
#####
COSMIC (Compact Object Synthesis and Monte Carlo Investigation Code) is a rapid binary population synthesis suite with a special focus of generating compact binary populations. 

COSMIC currently implements stellar evolution using SSE (`Hurley, Pols, and Tout 2000 <https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract>`_) and binary interactions using BSE (`Hurley, Tout, and Pols 2002 <https://ui.adsabs.harvard.edu/abs/2002MNRAS.329..897H/abstract>`_). Several modifications have been applied to BSE to account for recent updates to binary evolution especially important to compact binary formation (e.g. metallicity-dependent stellar winds or black hole natal kick strengths). For a detailed discussion of these modifications, see `Breivik et al. 2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv191100903B/abstract>`_.

************
Using COSMIC
************
COSMIC's primary purpose is to generate synthetic populations with an adaptive size based on how the shape of binary parameter disributions change as the number of simulated binaries increases.
This is done using the cosmic-pop exectuable which is installed with COSMIC. See :ref:`fixedpop` for more details.

COSMIC can also be used to simulate a single binary at a time, a list of multiple binaries, a grid of binaries, or a fixed population size as well as restart binaries at a mid point in their evolution. For more details on how to use COSMIC in these ways see :ref:`examples` and :ref:`runpop`.
