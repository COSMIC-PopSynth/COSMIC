**************************
ClusterMonteCarlo Sampling
**************************

New in COSMIC 3.4, you can now use COSMIC to sample initial conditions that can be used in the simulation of a Globular Cluster (GC), using the ClusterMonteCarlo (CMC) software package. To create these initial conditions, and save them in a format readable by CMC, you can do the following.

.. ipython::

    In [1]: from cosmic.sample.initialcmctable import InitialCMCTable

    In [2]: from cosmic.sample.sampler import cmc

To see the arguments necessary to call the CMC sampler use the help function:

.. ipython::

    In [3]: help(cmc.get_cmc_sampler)

.. ipython::
    :okwarning:

    In [1]: from cosmic.sample import InitialCMCTable
      
    In [2]: Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=0.2, primary_model='kroupa01', ecc_model='sana12', porb_model='sana12', qmin=-1.0, cluster_profile='plummer', met=0.014, size=40000, params='../examples/Params.ini', gamma=4, r_max=100)

    In [3]: InitialCMCTable.write(Singles, Binaries, filename="input.hdf5")

    In [4]: InitialCMCTable.write(Singles, Binaries, filename="input.fits")
