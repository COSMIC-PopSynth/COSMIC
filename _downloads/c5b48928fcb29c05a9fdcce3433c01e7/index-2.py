from cosmic.utils import a_from_p
from cosmic.sample.initialbinarytable import InitialBinaryTable
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
final_kstar = np.linspace(0,14,15)
colors = {'green' : '#1b9e77', 'purple' : '#d95f02', 'orange' : '#7570b3'}
final_kstar = np.linspace(0,14,15)
initC_mult, m_sin_mult, m_bin_mult, n_sin_mult, n_bin_mult = InitialBinaryTable.sampler('multidim',
                                                                                        final_kstar1=final_kstar,
                                                                                        final_kstar2=final_kstar,
                                                                                        rand_seed=2,
                                                                                        nproc=1,
                                                                                        SF_start=13700.0,
                                                                                        SF_duration=0.0,
                                                                                        met=0.02,
                                                                                        size=100000)
initC_mult['sep'] = a_from_p(p=initC_mult.porb, m1=initC_mult.mass_1, m2=initC_mult.mass_2)
