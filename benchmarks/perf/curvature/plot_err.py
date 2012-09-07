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

# rc("text",usetex=True)
# rc('font', family='serif')

list_file = sys.argv[1:len(sys.argv)]


for f in list_file :

    data = np.loadtxt(f, skiprows=1)

    h = data[:,1]
    nod = data[:,2]
    l2 = data[:,10]
    l2_int = data[:,16]
    sm = data[:,18]
    hs_proj = data[:,20]
    hs = data[:,22]

    nod_coeff = np.polyfit( np.log(h), np.log(nod), 1)
    l2_coeff = np.polyfit( np.log(h), np.log(l2), 1)
    l2_int_coeff = np.polyfit( np.log(h), np.log(l2_int), 1)
    sm_coeff = np.polyfit( np.log(h), np.log(sm), 1)
    hs_proj_coeff = np.polyfit( np.log(h), np.log(hs_proj), 1)
    hs_coeff = np.polyfit( np.log(h), np.log(hs), 1)

    plt.figure(1)
    plt.plot(np.log(h), np.log(nod), marker="+", linestyle="-", label="nod coeff = %.2f "%nod_coeff[0]+f)

    plt.figure(2)
    plt.plot(np.log(h), np.log(l2), marker="+", linestyle="-", label="ls coeff = %.2f "%l2_coeff[0]+f)

    plt.figure(3)
    plt.plot(np.log(h), np.log(l2_int), marker="+", linestyle="-", label="l2_int coeff = %.2f "%sm_coeff[0]+f)

    plt.figure(4)
    plt.plot(np.log(h), np.log(sm), marker="+", linestyle="-", label="sm coeff = %.2f "%sm_coeff[0]+f)

    plt.figure(5)
    plt.plot(np.log(h), np.log(hs_proj), marker="+", linestyle="-", label="hs_proj coeff = %.2f "%hs_proj_coeff[0]+f)

    plt.figure(6)
    plt.plot(np.log(h), np.log(hs), marker="+", linestyle="-", label="hs coeff = %.2f "%hs_coeff[0]+f)

    del data

plt.ylabel("log(e)")
plt.xlabel("log(h)")
plt.legend()
plt.show()



