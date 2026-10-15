"""Tests of the Shariat NS law and its optional extension to BHs."""
from copy import deepcopy
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.special import ndtr
from scipy.stats import kstest
from cosmic import utils
from cosmic.evolve import Evolve, INITIAL_BINARY_TABLE_SAVE_COLUMNS
from cosmic.sample.initialbinarytable import InitialBinaryTable, INITIAL_CONDITIONS_COLUMNS_ALL

SSE = {'stellar_engine': 'sse', 'path_to_tracks': '', 'path_to_he_tracks': ''}

def config(flag=9):
    initial = pd.read_hdf(Path(__file__).parent/'data/initial_conditions_for_testing.hdf5', 'initC')
    columns = list(set(INITIAL_BINARY_TABLE_SAVE_COLUMNS)-set(INITIAL_CONDITIONS_COLUMNS_ALL))
    bse = initial[columns].to_dict(orient='index')[0]
    for key in ['bin_num', 'randomseed']:
        bse.pop(key, None)
    bse.update(kickflag=flag, ecsn=0., ecsn_mlow=0., polar_kick_angle=90.,
        qcrit_array=[0.]*16, fprimc_array=[2./21]*16,
        natal_kick_array=[[-100.,-100.,-100.,-100.,0.],[-100.,-100.,-100.,-100.,0.]])
    return bse

def population(mass=20.005, n=128):
    table = InitialBinaryTable.InitialBinaries(m1=mass, m2=0., porb=0., ecc=-1.,
        tphysf=100., kstar1=1, kstar2=15, metallicity=0.02)
    table=table.loc[table.index.repeat(n)].reset_index(drop=True)
    table['randomseed'] = -(np.arange(n)+740001)
    return table

def run(table, bse, nproc=1):
    return Evolve.evolve(initialbinarytable=table.copy(), BSEDict=deepcopy(bse),
                         SSEDict=SSE, nproc=nproc, progress=False)

def mixture_cdf(v):
    result=np.zeros_like(v, dtype=float)
    for weight, mu, sig in [(0.126,1.87,0.55),(0.874,5.62,0.71)]:
        lo,hi=ndtr((np.log([0.05,1000.])-mu)/sig)
        result += weight*np.clip((ndtr((np.log(v)-mu)/sig)-lo)/(hi-lo),0,1)
    return result

@pytest.mark.parametrize("mass,remnant_type", [(20.005,13),(30.,14)])
def test_shariat_distribution_and_directions(mass,remnant_type):
    table=population(mass=mass,n=12000)
    bse=config()
    bse.update(bhflag=3,bhsigmafrac=1.)
    bpp,_,_,info=run(table,bse)
    events=info.loc[info.star > 0]
    assert len(events)==len(table)
    assert (bpp.kstar_1==remnant_type).any()
    speed=events.natal_kick.to_numpy()
    assert ((speed>0.05)&(speed<1000.)).all()
    threshold=np.sqrt(np.log(2/0.001)/(2*len(speed)))
    assert kstest(speed,mixture_cdf).statistic < threshold
    assert kstest(np.sin(np.deg2rad(events.phi)), 'uniform', args=(-1,2)).statistic < threshold
    assert kstest(events.theta/360., 'uniform').statistic < threshold

@pytest.mark.parametrize('magnitude',[0.,0.01,123.4,1200.])
@pytest.mark.parametrize('mass',[20.005,30.])
def test_shariat_supplied_magnitude(magnitude,mass):
    bse=config()
    bse.update(bhflag=3,bhsigmafrac=0.2)
    bse['natal_kick_array'][0][0]=magnitude
    info=run(population(mass=mass,n=4),bse)[-1]
    np.testing.assert_array_equal(info.loc[info.star>0,'natal_kick'], magnitude)

@pytest.mark.parametrize('bhflag',range(5))
@pytest.mark.parametrize('bhsigmafrac',[0.2,0.5,1.])
def test_shariat_black_hole_scaling(bhflag,bhsigmafrac):
    table=population(mass=30., n=32)
    base=config(9)
    base.update(bhflag=3,bhsigmafrac=1.)
    reference=run(table,base)
    assert (reference[0].kstar_1==14).any()
    remnant_mass=reference[0].loc[reference[0].kstar_1==14,'mass_1'].iloc[0]
    raw=reference[3].loc[reference[3].star>0,'natal_kick'].to_numpy()
    current=deepcopy(base)
    current.update(bhflag=bhflag,bhsigmafrac=bhsigmafrac)
    result=run(table,current)[3]
    kick=result.loc[result.star>0,'natal_kick'].to_numpy()
    if bhflag==0:
        np.testing.assert_array_equal(kick,0.)
    elif bhflag==2:
        np.testing.assert_allclose(kick,raw*bhsigmafrac*base['mxns']/remnant_mass,rtol=1e-13)
    elif bhflag==3:
        np.testing.assert_allclose(kick,raw*bhsigmafrac,rtol=1e-13)
    else:
        # Both existing fallback flags apply the same scaling to drawn kicks.
        other=deepcopy(current)
        other['bhflag']=4 if bhflag==1 else 1
        alternate=run(table,other)[3]
        np.testing.assert_array_equal(kick,alternate.loc[alternate.star>0,'natal_kick'])
        assert (kick>=0.).all()
        assert (kick<=raw*bhsigmafrac).all()

def test_shariat_repeat_and_multiprocessing():
    table=population(n=64)
    first=run(table,config())
    for second in [run(table,config()),run(table,config(),nproc=2)]:
        for a,b in zip(first,second):
            pd.testing.assert_frame_equal(a,b,check_exact=True)

def test_shariat_option_validation():
    utils.error_check(config(), SSE)
    with pytest.raises(ValueError):
        utils.error_check(config(-9), SSE)
