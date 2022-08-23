import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.integrate import simps

kB=8.6173324E-05 # eV/K
convRyd2eV=13.6057 # eV/Ryd
convcm2Ryd=1./109737.26 

def get_index(tran, idx):
    return [el[idx] for el in tran]


def initialize_upsilon_dataframe(temperature, nlevmax, nlevmin):
    if nlevmax is None or nlevmin is None: 
        print(' must define nlevmax/nlevmin ')
        sys.exit()
    tran = [(i+1, j+1) for i in range(nlevmax) for j in range(nlevmin) if (i != j) & (i >= j)] 
    ups = pd.DataFrame({'k': get_index(tran, 0), 'i': get_index(tran, 1)})
    for t in temperature:
        ups[t] = 0
    return ups


def compute_ECS(terms, omg, ener, temperature, debug=False):
    # number of points for energy integral
    npts = 2 ** 13 + 1 
    ntran = len(omg)

    # create effective collision strength dataframe
    ups = pd.DataFrame({'k': omg['k'], 'i': omg['i'], 'aki': omg['aki']})

    for t in temperature: 
        ups[t] = 0

    print('calculating ECS...')
    for tran in range(ntran):
        i, k, eik = data_trans(tran, omg, terms)

        # check eik (if all is good, nothing should be printed)
        if eik == 0: 
            print(f">>>>>>>> tran: {tran:>6d}, {i:>3d} => {k:>3d}, eik == 0")
            continue

        # copy transition data into aux variables
        yrow = omg.loc[tran][3:]
        ntype = int(yrow.pop('type'))
        yinf = yrow.pop('inf')
        xener = np.array(ener)
        yomg = np.array(yrow)
        yomg[yomg < 0.] = 0.

        if debug: 
            print(f"tran: {tran:>6d}, {i:>3d} => {k:>3d} (type {ntype})")
        # map into BT-space
        x, y = map_BTspace(eik, xener, yomg, ntype, yinf)
        # make spline of x, y in BT-space
        x_spline, y_spline = make_spline(x, y, npts)
        # mapping back from BT-space
        enerpp, omgpp = mapback_BTspace(ntype, eik, x_spline, y_spline)

        # compute effective collision strength
        yups = calculate_integral(enerpp, omgpp, temperature)
        ups.iloc[tran,3:] = yups

        # plot and check data
        iplot = 0
        if iplot == 1:
            plot_om2ups(xener, yomg, 
                        enerpp, omgpp, 
                        x, y, 
                        x_spline, y_spline, 
                        temperature, yups)

    print('... all done!\n')
    return ups


def data_trans(tran, omg, terms):
    k, i = int(omg.loc[tran]['k']), int(omg.loc[tran]['i'])
    Aki = omg.loc[tran]['aki']
    jk, ji = terms.loc[k]['J'], terms.loc[i]['J']
    gk, gi = 2*jk+1, 2*ji+1
    ek, ei = terms.loc[k]['Energy'], terms.loc[i]['Energy']
    eik = abs(ek-ei)
    fik = 0
    if eik != 0: 
        S = 3.73491E-10*gk*Aki/eik**3
        fik = eik*S/(3.*gi)
    return i, k, eik


def interp_omg(n,x,y,xnew,kind):
    if kind=="cubic": spline=interpolate.interp1d(x, y, kind=kind)
    if kind=="akima": spline=interpolate.Akima1DInterpolator(x,y)
    ynew=spline(xnew)
    return ynew


def make_spline(x, y, npts):
    nx = len(x)
    x_spline = np.linspace(x[0],x[nx-1], npts)
    y_spline = interp_omg(nx, x, y, x_spline, 'akima')
    return x_spline, y_spline


def map_BTspace(eik, xener, yomg, ntype, yinf):
    C = np.e
    xx = xener/eik
    if ntype == 1:
        x = 1-np.log(C)/np.log(xx+C)
        y = yomg/np.log(xx+np.e)
    if ntype == 2:
        x = xx/(xx+C)
        y = yomg
    if ntype == 3:
        x = xx/(xx+C)
        y = (xx+1)**2*yomg
    if ntype == 4:
        x = 1-np.log(C)/np.log(xx+C)
        y = yomg/np.log(xx+C)
    # include first point in BT-space
    lx, ly = list(x), list(y)
    # for neutral atoms, collision strenghts are zero at threshold
    lx.insert(0,0)
    ly.insert(0,0)
    x, y = np.array(lx), np.array(ly)
    # include last point in BT-space
    lx = list(x)
    xlast = 0.99
#     if ntype==1: xlast=0.99
#     if ntype==3: xlast=0.99
    if lx[-1] > xlast: xlast = abs(lx[-1]-1.)/2.+lx[-1]
    lx.append(xlast)
    x = np.array(lx)
    if ntype == 1 or ntype == 4: ylast = abs(yinf)
    if ntype == 2 or ntype == 3: ylast = y[-2]+(xlast-x[-2])/(x[-1]-x[-2])*(y[-1]-y[-2])
    ly = list(y)
    ly.append(ylast)
    y = np.array(ly)
    return x, y


def mapback_BTspace(ntype, eik, x_new, y_new):
    C = np.e
    if ntype == 1:
        arg = np.log(C)/(1.-x_new)
        enerp = (np.exp(arg)-C)*eik
        omgp = y_new*np.log(enerp/eik+np.e)
    if ntype == 2:
        enerp = x_new/(1.-x_new)*C*eik
        omgp = y_new
    if ntype == 3:
        enerp = x_new*C/(1.-x_new)*eik
        omgp = y_new/(enerp/eik+1)**2
    if ntype == 4:
        arg = np.log(C)/(1.-x_new)
        enerp = (np.exp(arg)-C)*eik
        omgp = y_new*np.log(enerp/eik+C)
    omgp[omgp < 0.] = 0.
    # nmax=npts-1
    nmax = len([i for i in enerp if i < 10.])
    enerpp = enerp[:nmax]
    omgpp = omgp[:nmax]
    return enerpp, omgpp


def calculate_integral(enerpp, omgpp, temperature):
    yups=[]
    for t in temperature:
        kBT=kB*t
        arg=enerpp*convRyd2eV/kBT
        yint=omgpp*np.exp(-arg)
        I1=simps(yint, arg)
        yups.append(I1)
    return yups

def plot_om2ups(xener,yomg,enerpp,omgpp,x,y,x_new,y_new,T,ups):
    fig, axs = plt.subplots(1, 3)
    axs[0].plot(xener, yomg, 'o-')
    axs[0].plot(enerpp, omgpp, 'r-')
    axs[1].plot(x, y, 'o-')
    axs[1].plot(x_new, y_new, 'r-')
    axs[2].plot(T, ups, 'r-', label="present")
    axs[2].legend(loc='upper left', bbox_to_anchor=[1,1], ncol=1)
    plt.show()

