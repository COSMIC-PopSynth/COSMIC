import scipy.integrate
import numpy as np

##################################################################################
# DEFINE MIN AND MAX MASS SELECTOR
##################################################################################
def mass_min_max_select(kstar_1, kstar_2):
    primary_max = 150.0
    secondary_max = 150.0

    primary_min = 0.08
    secondary_min = 0.08

    kstar = [kstar_1, kstar_2]
   
    min_mass = [primary_min, secondary_min]
    max_mass = [primary_max, secondary_max]
   
    ii = 0
    for k in kstar:
        if k > 13:
            min_mass[ii] = 15.0
        elif k > 12:
            min_mass[ii] = 8.0
            max_mass[ii] = 50.0
        elif k > 11:
            min_mass[ii] = 5.0
            max_mass[ii] = 15.0
        elif k > 10:
            min_mass[ii] = 2.0
            max_mass[ii] = 10.0
        elif k > 9:
            min_mass[ii] = 0.5
            max_mass[ii] = 5.0
        else:
            max_mass[ii] = 5.0
        ii += 1

    return min_mass[0], max_mass[0], min_mass[1], max_mass[1]


def idl_tabulate(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = scipy.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in xrange(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


def rndm(a, b, g, size):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)
