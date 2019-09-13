import numpy as np
import os
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.optimize as siopt
from matplotlib import rc
from scipy import *
from numpy import *
from pylab import *

from headerFile import *

# rc("text",usetex=True)
# rc('font', family='serif')

list_file = sys.argv[1:len(sys.argv)]

def line_fit(_h, _e, _coeff) :
    hh = np.arange( min(_h)/2, 2* max(_h), (2*max(_h) - min(_h) / 2.) / 50. )
    c0 = _e[-1] / _h[-1]**_coeff
    yy = [ c0 * _y**_coeff for _y in hh ]
    return (hh, yy, _coeff)


for f in list_file :

    hf = HeaderFile( f, skiprows=1 )

    porder = 1.

    h = hf["h"]
    nod = hf["e.nod.k"]
    nod_fit = line_fit(h, nod, hf["e.nod.k.roc"][-1])

    l2 = hf["e.l2.k"]
    l2_fit = line_fit(h, l2, hf["e.l2.k.roc"][-1])

    sm = hf["e.sm.k"]
    sm_fit = line_fit(h, sm, hf["e.sm.k.roc"][-1])

    hs = hf["e.hs.k"]
    hs_fit = line_fit(h, hs, hf["e.hs.k.roc"][-1])

    opt = hf["e.opt.k"]
    opt_fit = line_fit(h, opt, hf["e.opt.k.roc"][-1])

    plt.figure(1)
    plt.title("nodal error")
    plt.plot(np.log(h), np.log(nod), marker="^", linestyle="", label="")
    plt.plot(np.log(nod_fit[0]), np.log(nod_fit[1]), linestyle="-", marker="", label="%.2f"%nod_fit[2])
    plt.ylabel("error")
    plt.xlabel("h")
    plt.legend()

    # plt.figure(2)
    # plt.title("L2 error")
    # plt.plot(np.log(h), np.log(l2), marker="^", linestyle="", label="")
    # plt.plot(np.log(l2_fit[0]), np.log(l2_fit[1]), linestyle="-", marker="", label="%.2f"%l2_fit[2])
    # plt.ylabel("error")
    # plt.xlabel("h")
    # plt.legend()

    plt.figure(3)
    plt.title("smooth error")
    plt.plot(np.log(h), np.log(sm), marker="^", linestyle="", label="")
    plt.plot(np.log(sm_fit[0]), np.log(sm_fit[1]), linestyle="-", marker="", label="%.2f"%sm_fit[2])
    plt.ylabel("error")
    plt.xlabel("h")
    plt.legend()

    # plt.figure(4)
    # plt.title("hessian error")
    # plt.plot(np.log(h), np.log(hs), marker="^", linestyle="", label="")
    # plt.plot(np.log(hs_fit[0]), np.log(hs_fit[1]), linestyle="-", marker="", label="%.2f"%hs_fit[2])
    # plt.ylabel("error")
    # plt.xlabel("h")
    # plt.legend()

    plt.figure(5)
    plt.title("optimal error")
    plt.plot(np.log(h), np.log(opt), marker="^", linestyle="", label="")
    plt.plot(np.log(opt_fit[0]), np.log(opt_fit[1]), linestyle="-", marker="", label="%.2f"%opt_fit[2])
    plt.ylabel("error")
    plt.xlabel("h")
    plt.legend()


plt.show()



