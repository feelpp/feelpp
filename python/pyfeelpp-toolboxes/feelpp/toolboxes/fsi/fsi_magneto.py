import sys


import feelpp.core as fppc

import feelpp.core.quality as quality
import feelpp.toolboxes.core as tb
import feelpp.core.interpolation as I
from feelpp.toolboxes.fsi import *
import feelpp.core.meshmover as mm
import mpi4py
mpi4py.rc.thread_level="single"
import pandas as pd
import json
#import matplotlib.pyplot as plt
import numpy as np

sys.argv = ['magneto']
e = fppc.Environment(sys.argv, opts=tb.toolboxes_options("fsi"),
                    config=fppc.globalRepository("magneto"))

fppc.Environment.setConfigFile('magneto.cfg')

#============== Control parameters =======================#
freq = 0.9
ux = lambda t : 0.005
uy = lambda t : 0.005 * np.sin(2*np.pi*freq*t)
#=========================================================#

fsi_tb = fsi(dim=2, orderU=2, orderP=1, orderGeo=1)
fsi_tb.init()
#fsi_tb.printAndSaveInfo()


#Add Torque FSI
fsi_tb.addMagnetoTorqueModelFSI()
fsi_tb.addMagnetoTroqueResModelFSI()

fsi_tb.startTimeStep()

while not fsi_tb.timeStepBase().isFinished():
 
    if fppc.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s iteration : {}\n".format(fsi_tb.time(), fsi_tb.timeStepBase().iteration()))
        print("============================================================\n")
    

    #Update control at time t
    uxt = ux(fsi_tb.time())
    uyt = uy(fsi_tb.time())
    fsi_tb.addParameterInModelProperties("uxt", uxt)
    fsi_tb.addParameterInModelProperties("uyt", uyt)
    fsi_tb.addParameterInModelProperties("uzt", 0)
    fsi_tb.updateParameterValues()


    #Solve FSI
    fsi_tb.solve()

    #Export results
    fsi_tb.exportResults()

    #Update time : t <- t+dt
    fsi_tb.updateTimeStep()
