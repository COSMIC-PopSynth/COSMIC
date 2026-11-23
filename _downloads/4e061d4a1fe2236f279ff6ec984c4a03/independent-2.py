import numpy as np
import matplotlib.pyplot as plt
from cosmic.sample.initialbinarytable import InitialBinaryTable

final_kstars = np.linspace(0, 14, 15)

common_kwargs = {
    'final_kstar1' : final_kstars,
    'final_kstar2' : final_kstars,
    'binfrac_model' : 1.0,
    'primary_model' : 'kroupa01',
    'ecc_model' : 'sana12',
    'porb_model' : 'sana12',
    'SF_start' : 13700.0,
    'SF_duration' : 0.0,
    'met' : 0.02,
    'size' : 100000,
    'trim_extra_samples' : True
}

init_bin_default, _, _, _, _ = InitialBinaryTable.sampler(
    'independent',
    **common_kwargs
)

init_bin_qmin, _, _, _, _ = InitialBinaryTable.sampler(
    'independent',
    qmin=0.2,
    **common_kwargs
)

init_bin_m2min, _, _, _, _ = InitialBinaryTable.sampler(
    'independent',
    m2_min=0.08,
    **common_kwargs
)

init_bin_qpower, _, _, _, _ = InitialBinaryTable.sampler(
    'independent',
    q_power_law=3,
    **common_kwargs
)

fig, ax = plt.subplots(figsize=(8,6))
ax.hist(init_bin_default.mass_2/init_bin_default.mass_1, bins=np.linspace(0, 1, 60),
        histtype='step',
        lw=3, label='default (qmin=-1)')
ax.hist(init_bin_qmin.mass_2/init_bin_qmin.mass_1, bins=np.linspace(0, 1, 60),
        histtype='step',
        lw=3, label='qmin=0.2')
ax.hist(init_bin_m2min.mass_2/init_bin_m2min.mass_1, bins=np.linspace(0, 1, 60),
        histtype='step',
        lw=3, label='m2_min=0.08 M$_\odot$')
ax.hist(init_bin_qpower.mass_2/init_bin_qpower.mass_1, bins=np.linspace(0, 1, 60),
        histtype='step',
        lw=3, label='q_power_law=3')
ax.set_xlabel(r'Mass ratio q=M$_2$/M$_1$', size=18)
ax.set_ylabel('Counts', size=18)
ax.legend(prop={'size' : 14})
fig.tight_layout()
fig.show()