.. _examples:

#######################
Using cosmic to run BSE
#######################


cosmic can evolve binaries for several different use cases. Below 
you'll find examples to run a single binary system, multiple binary
systems or a grid of binaries.

See 'Running a Milky Way population' below for the process to simulate
a population of binaries consistent with a user-supplied star formation 
history for a single compact object population (e.g. BH-BH) or a range
of compact object populations (e.g. combinations of BH, NS, WD) as 
described in Breivik+2018.



*************
single binary
*************

Below is the process to initialize and evolve a binary that 
could have formed a GW150914-like binary. First, import the modules in cosmic
that initialize and evolve the binary.

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: from cosmic.evolve import Evolve


To initialize a single binary, populate the SingleBinary method in the
InitialBinaryTable class. Each initialized binary requires the following parameters:


* m1 : ZAMS mass of the primary star in :math:`M_{\odot}`

* m2 : ZAMS mass of the secondary star in :math:`M_{\odot}`

* porb : initial orbital period in days

* ecc : initial eccentricity

* tphysf : total evolution time of the binary in Myr

* kstar1 : initial primary stellar type, following the BSE convention

* kstar2 : initial secondary stellar type, following the BSE convention

* metallicity : metallicity fraction (e.g. :math:`Z_{\odot}=0.02`)
 
.. ipython::

    In [3]: single_binary = InitialBinaryTable.SingleBinary(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)

    In [4]: print(single_binary)

The flags for the various binary evolution prescriptions used in BSE also need to be set. 
Each flag is saved in the BSEDict dictionary. Note that the BSEDict
only needs to be specified the first time a binary is evolved with cosmic or
if you need to change the binary evolution prescriptions. 

If you are unfamiliar with these prescriptions, it is highly 
advised to either run the defaults from the cosmic install (which are consistent
with `Rodriguez+2018 <http://adsabs.harvard.edu/abs/2018PhRvL.120o1101R>`_ and `Kremer+2018 <http://adsabs.harvard.edu/abs/2018PhRvL.120s1103K>`_) or refer to `Hurley+2002 <http://adsabs.harvard.edu/abs/2002MNRAS.329..897H>`_.

.. ipython::

    In [5]: BSEDict = {'xi': 0.5, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 0, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 1.0, 'ck': -1000, 'bwind': 0.0, 'lambdaf': 1.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'nsflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': -3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [-100.0,-100.0,-100.0,-100.0,-100.0,-100.0], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.5, 'ecsn_mlow' : 1.4, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 3, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0]} 

Once the binary is initialized and the BSE model is set, the system is evolved with the 
the Evolve class, which calls the evolv2.f subroutine in the BSE source code. 

.. ipython:: 

    In [6]: bpp, bcm, initC = Evolve.evolve(initialbinarytable=single_binary, BSEDict=BSEDict)


For every evolved binary system, BSE generates two arrays, which are stored as pandas DataFrames in cosmic:

* bpp - contains binary parameters at important stages in the binary's evolution, including stellar evolutionary phase changes or mass transfer episodes.

* bcm - contains several binary parameters at user specified time steps during the binary's evolution. The default setting in cosmic is to output the final stage of the binary at the evolution time specified by the user.

You can see the different parameters included in each DataFrame using the columns attribute of the DataFrame:

.. ipython:: 

    In [7]: print(bpp.columns)

    In [8]: print(bcm.columns)

The units are broadly consistent with BSE; masses, radii are in Solar units, temperatures are in log-Kelvin, luminosities are in log-Solar luminosity, orbital periods are in log-seconds, semi-major axes are in solar radii, and times are in units of Myr.

The evol_type column in bpp indicates the evolutionary change that occured for each line. The meaning of each number is described here, :ref:`evolve-type-table`.

Each of the parameters in bpp or bcm can be accessed in the usual way for DataFrames.

.. ipython:: 

    In [9]: print(bpp.mass_1)

    In [10]: print(bpp[['mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'sep', 'evol_type']])


*****************
multiple binaries
*****************

Multiple systems can also be initialized and evolved; below is an example for systems
that could form GW150914 and GW170817 - like binaries.

.. ipython::

    In [11]: binary_set = InitialBinaryTable.MultipleBinary(m1=[100.0, 11.8], m2=[85.0, 11.1], porb=[10000.0,2211.0], ecc=[0.65,0.55], tphysf=[13700.0,13700.0], kstar1=[1,1], kstar2=[1,1], metallicity=[0.005,0.02])

    In [12]: print(binary_set)

    In [13]: bpp, bcm, initC  = Evolve.evolve(initialbinarytable=binary_set, BSEDict=BSEDict)

Note that the BSEDict did not be reinitialized since the BSE model did not change.

As before, bpp, bcm, and initC are returned as pandas DataFrames which assign an index to each binary system we evolve. We can access each binary as follows

.. ipython:: 

    In [14]: print(bpp.loc[0])

    In [15]: print(bcm.loc[0])

    In [16]: print(initC.loc[0])

    In [17]: print(bpp.loc[1])


****************
grid of binaries
****************

Sometimes it is helpful to run a grid of initial binaries to explore how
changing a single paramter affects the evolved binary. Here we evolve 
the same system that produces a GW150914-like binary, but run over several initial orbital
periods spaced evenly in log space.

.. ipython::

    In [16]: n_grid = 10 

    In [17]: binary_grid = InitialBinaryTable.MultipleBinary(m1=np.ones(n_grid)*100.0, m2=np.ones(n_grid)*85.0, porb=np.logspace(3,5,n_grid), ecc=np.ones(n_grid)*0.65, tphysf=np.ones(n_grid)*13700.0, kstar1=np.ones(n_grid), kstar2=np.ones(n_grid), metallicity=np.ones(n_grid)*0.005)

    In [18]: print(binary_grid)

    In [19]: bpp, bcm, initC  = Evolve.evolve(initialbinarytable=binary_grid, BSEDict=BSEDict)

    In [20]: print(bpp)

    In [21]: print(bcm)

