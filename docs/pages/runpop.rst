.. _runpop:

############################
Sampling initial populations
############################

The process to generate a synthetic binary population, is similar to the
process to evolve a single/multiple binaries by hand: first generate an
initial population, then evolve it with the :class:`~cosmic.evolve.Evolve` class.

An initialized binary population consists of a collection of binary systems
with assigned primary and secondary masses, orbital periods, eccentricities,
metallicities, and star formation histories. These parameters are randomly
sampled from observationally motivated distribution functions.

In COSMIC, the initial sample is done through an initial binary sampler which works
with the :class:`~cosmic.sample.initialbinarytable.InitialBinaryTable` class. There are two samplers available:

1. ``independent`` : initialize binaries with independent parameter
distributions for the primary mass, mass ratio, eccentricity, separation,
and binary fraction

2. ``multidim`` : initialize binaries with multidimensional parameter
distributions according to `Moe & Di Stefano 2017 <http://adsabs.harvard.edu/abs/2017ApJS..230...15M>`_

We consider both cases in the guides below.

.. toctree::
    :maxdepth: 1

    sample/independent
    sample/multidim

You can also use COSMIC to sample the initial conditions for a Globular Cluster (GC) using the ClusterMonteCarlo (CMC) software package.
Check out the guide below for more information.

.. toctree::
    :maxdepth: 1

    sample/cluster

.. tip::

    The initial binary population that is generated with these sampling methods can be quickly and 
    easily evolved using the :class:`~cosmic.evolve.Evolve` class.
    This is demonstrated in `this guide <evolve/evolve_sample.html>`_.