import numpy as np
import scipy
from scipy import special
import sys

LARGE_DISTANCE = 1.0e40
DBL_MIN = 2.2250738585072014e-308
DBL_EPSILON = 2.2204460492503131e-016 

#Input W0 and NSTAR
W0 = np.float(sys.argv[1])
NSTAR = np.int(sys.argv[2])

OUTFILE = 'king_w'+np.str(W0)+'_n'+np.str(NSTAR)+'.dat'
filename = open(OUTFILE,'w')
SEED = 0

debug=0

class strpar:
    pass
strpar.rt = 0.0
strpar.rk = 0.0

def f0(t, x, y, w0):
    return(y)

def f1(t,x,y,w0):
    rho0 = np.exp(w0)*scipy.special.erf(np.sqrt(w0)) - np.sqrt(4.0*w0/np.pi)*(1.0+2.0*w0/3.0)
    result = -9.0*(np.exp(x)*scipy.special.erf(np.sqrt(x)) - np.sqrt(4.0*x/np.pi)*(1.0+2.0*x/3.0))/rho0
    if np.abs(y) < DBL_EPSILON and np.abs(t) < DBL_EPSILON:
        result = result/3.0
    else:
        result = result - 2.0*y/t
    #result = result/3.0
    return(result)

def calc_rho_rho0(t,x,y,w0):
    rho = np.exp(x)*scipy.special.erf(np.sqrt(x)) - np.sqrt(4.0*x/np.pi)*(1.0+2.0*x/3.0)
    rho0 = np.exp(w0)*scipy.special.erf(np.sqrt(w0)) - np.sqrt(4.0*w0/np.pi)*(1.0+2.0*w0/3.0)
    return(rho/rho0)

def incr(k1,k2,k3,k4,k5,k6):
    return(25.0/216.0*k1 + 1408.0/2565.0*k3 + 2197.0/4104.0*k4 - 1.0/5.0*k5)

def set_masses(m, N):
    mtot = 0.0
    m[0] = 0.0
    for i in range(1,N+1):
        m[i] = 1.0/N
    m[N+1] = 0
    return(m)

def scale_pos_and_vel(m, r, vr, vt, N, strpar):
    PEtot = 0.0
    KEtot = 0.0
    U = 0.0
    MM = 1.0 #because of units, the total mass has to be 1 initially 
    for i in range(N,0,-1):
        U -= MM*(1.0/r[i] - 1.0/r[i+1])
        T = 0.5 * (vr[i]*vr[i] + vt[i]*vt[i])
        MM -= m[i]
        PEtot += 0.5*U*m[i]
        KEtot += T*m[i]
    print("Before scaling: PEtot = %f, KEtot = %f, vir rat = %f\n" % (PEtot, KEtot, KEtot/PEtot))
    #scaling position and velocity
    rfac = -PEtot*2.0
    vfac = 1.0/np.sqrt(4.0*KEtot)
    for i in range(1,N+1):
        r[i] *= rfac
        vr[i] *= vfac
        vt[i] *= vfac

    PEtot = 0.0
    KEtot = 0.0
    U = 0.0
    MM = 1.0; # because of units, the total mass has to be 1 initially
    for i in range(N,0,-1):
        U -= MM*(1.0/r[i] - 1.0/r[i+1])
        T = 0.5 * (vr[i]*vr[i] + vt[i]*vt[i])
        MM -= m[i]
        PEtot += 0.5*U*m[i]
        KEtot += T*m[i]
    print("After  scaling: PEtot = %f, KEtot = %f, vir rat = %f\n" % (PEtot, KEtot, KEtot/PEtot)) 
    strpar.rk = rfac
    strpar.rt *= rfac
    print("value of r_0 is: %e\n" % (strpar.rk))
    print("value of r_t is: %e\n" % (strpar.rt))

def write_output_file(m, r, vr, vt, N, OUTFILE, strpar):
    NOBJ = N
    NBINARY = 0
    Mclus = N
    Rvir = 1.0
    Rtid = strpar.rt
    Z = 0.02

    #####cmc_malloc_fits_data_t(&cfd);
    obj_id = []
    obj_k = []
    obj_id, obj_k, obj_m, obj_Reff, obj_r, obj_vr, obj_vt, obj_binind = ([] for i in range(8))      
    for i in range(0,N+1):
        obj_id.append(i)
        obj_k.append(0)
        obj_m.append(m[i])
        obj_Reff.append(0.0)
        obj_r.append(r[i])
        obj_vr.append(vr[i])
        obj_vt.append(vt[i])
        obj_binind.append(0)

    for i in range(1,N+1):
        print(obj_id[i], obj_k[i], obj_m[i], obj_Reff[i], obj_r[i], obj_vr[i], obj_vt[i], obj_binind[i] ,file=OUTFILE)                  

def create_random_array(X,N):
    for i in range(1,N+1):
        # Is this sufficient?
        rand = np.random.uniform(0,1)
        X[i] = rand
    X2 = np.unique(X)
    while len(X2) != len(X):
        #print('Warning: X array contains duplicates')
        diff = np.int(len(X2) - len(X))
        for i in range(0,diff):
            X2.append(np.random.uniform(0,1))
        X2 = np.unique(X)
        X = X2
    return(X)

#### Main
w0 = W0
N = NSTAR
seed = SEED
tol = 1e-10
hmax = 0.01
hmin = 1e-11
error = 0.0

array_incr = 10000
array_size = array_incr

#/* reset_rng_t113(2718234592UL); */
#reset_rng_t113(seed);
#/* t is r/r0, x is Psi/sigma^2, y is dx/dt */
h = tol
t = 0.0 
x = w0
y = 0.0 
pos = []
den = []
pot = []

i = 0
while True:
    rho_rho0 = calc_rho_rho0(t,x,y,w0)
    pos.append(t)
    den.append(rho_rho0)
    pot.append(x)
    i += 1
    if i>=array_size:
        array_size += array_incr

    k1 = h*f0(t,x,y,w0)
    if np.isnan(k1) and h>hmin:
        h *= 0.5
        continue
    l1 = h*f1(t,x,y,w0)
    
    t2 = t + h/4.0
    x2 = x + k1/4.0 
    y2 = y + l1/4.0
    k2 = h*f0(t2,x2,y2,w0)
    if np.isnan(k2) and h>hmin:
        h*=0.5
        continue
    l2 = h*f1(t2,x2,y2,w0)
    
    t3 = t + 3.0/8.0*h 
    x3 = x + 3.0/32.0*k1 + 9.0/32.0*k2 
    y3 = y + 3.0/32.0*l1 + 9.0/32.0*l2 
    k3 = h*f0(t3,x3,y3,w0)
    if np.isnan(k3) and h>hmin:
        h*=0.5
        continue
    l3 = h*f1(t3,x3,y3,w0)
    
    t4 = t + 12.0/13.0*h
    x4 = x + 1932.0/2197*k1 - 7200.0/2197*k2 + 7296.0/2197*k3
    y4 = y + 1932.0/2197*l1 - 7200.0/2197*l2 + 7296.0/2197*l3
    k4 = h*f0(t4,x4,y4,w0)
    if np.isnan(k4) and h>hmin:
        h*=0.5
        continue
    l4 = h*f1(t4,x4,y4,w0)

    t5 = t + h
    x5 = x + 439.0/216.0*k1 - 8.0*k2 + 3680.0/513.0*k3 - 845.0/4104*k4
    y5 = y + 439.0/216.0*l1 - 8.0*l2 + 3680.0/513.0*l3 - 845.0/4104*l4
    k5 = h*f0(t5,x5,y5,w0)
    if np.isnan(k5) and h>hmin:
        h*=0.5
        continue
    l5 = h*f1(t5,x5,y5,w0)

    t6 = t + h/2.0
    x6 = x - 8.0/27.0*k1 + 2.0*k2 - 3544.0/2565.0*k3 + 1859.0/4104.0*k4 - 11.0/40.0*k5
    y6 = y - 8.0/27.0*l1 + 2.0*l2 - 3544.0/2565.0*l3 + 1859.0/4104.0*l4 - 11.0/40.0*l5
    k6 = h*f0(t6,x6,y6,w0)
    if np.isnan(k6) and h>hmin:
        h*=0.5
        continue
    l6 = h*f1(t6,x6,y6,w0)
    
    error = 1.0/360.0*k1 - 128.0/4275.0*k3 - 2197.0/75240.0*k4 + 1.0/50.0*k5 + 2.0/55.0*k6
    if error>tol:
        h *= 0.5
        continue
    
    x += incr(k1, k2, k3, k4, k5, k6)
    y += incr(l1, l2, l3, l4, l5, l6)
    t += h
    if error<tol/2.0 and h<hmax/2.0:
        h *= 2.0
    if x < 0 or rho_rho0 < 0:
        break

num_point = i
strpar.rt = pos[num_point-1]
mR = np.zeros(num_point)
mR[0] = 0.0;
for i in range(1,num_point):
    mR[i] = mR[i-1] + (pos[i]-pos[i-1])*(den[i]*pos[i]*pos[i]+den[i-1]*pos[i-1]*pos[i-1])/2.0

for i in range(1,num_point):
    mR[i] /= mR[num_point-1]

r = np.zeros(N+2)
m = np.zeros(N+2)
psi = np.zeros(N+2)
vr = np.zeros(N+2)
vt = np.zeros(N+2)
X = np.zeros(N+2)

create_random_array(X, N)
r[0] = DBL_MIN   
for i in range(1,N+1):
    #XXX below is uniformly spaced points, use RNG for a more random distribution distribution in r XXX changed
    jmin = np.int(0)
    jmax = np.int(num_point)
    while(jmin != jmax):
        #loop invariant is: mR[jmin]<X[i]<mR[jmax+1]
        jtry = np.int((jmin+jmax+1.0)/2.0)
        if mR[jtry] > X[i]:
            jmax = jtry-1
        else:
            jmin = jtry
    if mR[jmin]>X[i] or mR[jmin+1]<=X[i]:
        #binary search failed!
        print("binary search failed!")
        print("i = %d, jmin = %d," & (i, jmin))
        print("mR[jmin]=%e, mR[jmin+1]=%e, X[i] = %e\n" & (mR[jmin], mR[jmin+1], X[i]))
    r[i] = pos[jmin] + (pos[jmin+1]-pos[jmin])*(X[i]-mR[jmin])/(mR[jmin+1]-mR[jmin])
    psi[i] = pot[jmin] + (pot[jmin+1]-pot[jmin])*(X[i]-mR[jmin])/(mR[jmin+1]-mR[jmin])
r[N+1] = LARGE_DISTANCE

vr[0] = 0.0
vt[0] = 0.0
for i in range(1,N+1):
    F = (np.exp(psi[i])-1.0)*2.0*psi[i]
    while True:
        X1 = np.random.uniform(0,1)
        X2 = np.random.uniform(0,1)
        v_0 = X1*np.sqrt(2.0*psi[i])
        f_0 = X2*F
        f = (np.exp(psi[i]-v_0*v_0/2.0)-1.0)*v_0*v_0
        if f_0 < f:
            break
    X3 = np.random.uniform(0,1)
    vr[i] = (1.0 - 2.0*X3) * v_0
    vt[i] = np.sqrt(v_0*v_0 - vr[i]*vr[i])
vr[N+1] = 0.0
vt[N+1] = 0.0

set_masses(m, N)
scale_pos_and_vel(m, r, vr, vt, N, strpar)
write_output_file(m, r, vr, vt, N, filename, strpar)

