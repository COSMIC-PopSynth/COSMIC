*******************
Re-running binaries
*******************

COSMIC allows you to re-run a binary from a previous evolution from a COSMIC generated initC table. This can
also allow you to modify binary physics assumptions to see how they affect the later evolution of a binary
with otherwise identical initial conditions.

Re-run an identical binary
==========================

You can re-run a binary from a COSMIC generated initC table, that's all you need! This means you
can send your favourite initial conditions around to your friends and they will be able to reproduce the
exact same evolution (assuming they use the same COSMIC version).

First, let's evolve a binary and save the initC table.

.. ipython:: python

    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.evolve import Evolve

    BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0,
               'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000,
               'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5,
               'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0,
               'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0,
               'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]],
               'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90,
               'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
               'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25,
               'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0,
               'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0],
               'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1,
               'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 1,
               'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'rtmsflag' : 0,
               'wd_mass_lim': 1}

    binary = InitialBinaryTable.InitialBinaries(
        m1=20, m2=15, porb=100, ecc=0.1, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.02
    )

    bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable=binary, BSEDict=BSEDict)

We can check some of the output for this binary so that we can see it's identical after we re-run it.

.. ipython:: python

    print(bpp)
    print(kick_info)

Now we can re-run this exact same binary by passing in the initC table we just generated instead of the
initial binary table. We also do **not** need to pass in the BSEDict again, as the physics assumptions are
stored in the initC table. In fact, you should not pass in a BSEDict when re-running from an initC since things
may get inconsistent (and confusing!).

.. ipython:: python

    bpp_rerun, bcm_rerun, initC_rerun, kick_info_rerun = Evolve.evolve(initialbinarytable=initC)

    print(bpp_rerun)
    print(kick_info_rerun)

so overall, we can see that these outputs are exactly identical, even down to the kick information.

.. warning::

    When re-running from an initC table, the BSEDict is no longer necessary since the physics assumptions are stored
    in the initC table. In fact, you should **not** pass in a BSEDict when re-running from an initC since things
    may get inconsistent (and confusing!).

Re-run with different physics
=============================

You can also re-run a binary from a COSMIC generated initC table but modify the physics assumptions.
This allows you to see how different physics assumptions affect the later evolution of a binary
with otherwise identical initial conditions.

Let's re-run the same binary as before, but this time we will modify the common envelope efficiency
parameter, alpha1, to be 10 instead of 1.0.

To do this, we modify the initC column corresponding to alpha1 before re-running the binary. We **do not**
pass in a BSEDict this time, since we are modifying the physics assumptions directly in the initC table.

.. ipython:: python

    # modify alpha1 in the initC table
    initC_modified = initC.copy()
    initC_modified['alpha1'] = 10

    # re-run the binary with modified physics
    bpp_rerun_mod, bcm_rerun_mod, initC_rerun_mod, kick_info_rerun_mod = Evolve.evolve(
        initialbinarytable=initC_modified
    )

    print(bpp_rerun_mod)
    print(kick_info_rerun_mod)

Now as we can see, this output is very different from the original evolution, since the common envelope now
proceeds quite differently with the higher efficiency parameter! The system goes from a merging binary that
produces a black hole, to an unbound binary that produces two neutron stars.

Ensuring natal kicks account for modified physics
-------------------------------------------------

When re-running a binary with modified physics, you may want to ensure that the natal kicks are re-sampled
appropriately given the new physics assumptions. For example, if you change the remnant mass prescription,
you may want to re-sample the natal kicks since the remnant masses will be different. If you change the natal
kick prescription itself, you will *definitely* want to re-sample the natal kicks.

In order to re-sample the natal kicks when re-running a binary, you need to erase the existing kick information
from the initC table before re-running. This is done by setting the kick columns to -100.0, which is the flag value
indicating that no kick has been assigned yet.

Let's try this out by re-running our system with higher alpha, but this time re-sampling the natal kicks as well with a
very small sigma value for both regular core-collapse supernovae and electron-capture supernovae.

.. ipython:: python

    # modify alpha1 in the initC table
    initC_modified_kick = initC.copy()
    initC_modified_kick['alpha1'] = 10
    initC_modified_kick['sigma'] = 1.0
    initC_modified_kick['sigmadiv'] = -1.0
    initC_modified_kick['kickflag'] = 1

    # erase existing kick information by setting to -100.0
    kick_columns = [
        'natal_kick_1', 'natal_kick_2', 'phi_1', 'phi_2',
        'theta_1', 'theta_2', 'mean_anomaly_1', 'mean_anomaly_2'
    ]
    for col in kick_columns:
        initC_modified_kick[col] = -100.0

    # re-run the binary with modified physics and re-sampled kicks
    bpp_rerun_mod_kick, bcm_rerun_mod_kick, initC_rerun_mod_kick, kick_info_rerun_mod_kick = Evolve.evolve(
        initialbinarytable=initC_modified_kick,
    )

    print(bpp_rerun_mod_kick)
    print(kick_info_rerun_mod_kick)

As one might expect, the binary no longer becomes unbound since the natal kicks are now very small!
The system instead produces two neutron stars in a bound system which eventually inspirals due to gravitational
wave emission and produces a NS-NS merger.
