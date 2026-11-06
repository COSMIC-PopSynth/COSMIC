from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.plotting import evolve_and_plot
single_binary = InitialBinaryTable.InitialBinaries(m1=85.543645, m2=84.99784, porb=446.795757, ecc=0.448872, tphysf=13700.0, kstar1=1, kstar2=1, metallicity=0.002)
fig = evolve_and_plot(single_binary, t_min=None, t_max=None, BSEDict=default_BSEDict, SSEDict=default_SSEDict, sys_obs={})