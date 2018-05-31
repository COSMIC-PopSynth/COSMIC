.. _runpop:

#############################
Evolve a Milky Way population
#############################

To generate a Milky Way like population, we need to generate an initial set of
binaries that is representative of a population born in the Milky Way. This means
we need to supply several arguments to our initial binary sampler. There are two 
distinct ways to generate the initial population:

1. assume all binary paramters are independent of one another

2. try to account for parameter dependencies

We consider both cases below. 

************************************************************************
Sampling an initial population with independently distributed parameters
************************************************************************

The 

.. ipython::

    In [1]: from cosmic.sample.initialbinarytable import InitialBinaryTable

    In [2]: IBT, sampled_mass = InitialBinaryTable.sampler('independent', primary_min=0.08, primary_max=5.0, primary_model='kroupa93', ecc_model='thermal', SFH_model='const', component_age=10000.0, size=1000)

    In [3]: print(IBT)


