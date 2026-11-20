from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.plotting import evolve_and_plot
import numpy as np
np.random.seed(5)
binary_set = InitialBinaryTable.InitialBinaries(m1=[85.543645, 11.171469], m2=[84.99784, 9.67305], porb=[446.795757, 370.758343], ecc=[0.448872, 0.370], tphysf=[13700.0, 13700.0], kstar1=[1, 1], kstar2=[1, 1], metallicity=[0.002, 0.02])
fig = evolve_and_plot(binary_set, t_min=None, t_max=[6.0, 60.0], BSEDict=default_BSEDict, SSEDict=default_SSEDict, sys_obs={})