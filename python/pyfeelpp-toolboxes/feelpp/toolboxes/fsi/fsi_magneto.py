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
import matplotlib.pyplot as plt
import numpy as np

sys.argv = ['magneto']
e = fppc.Environment(sys.argv, opts=tb.toolboxes_options("fsi"),
                    config=fppc.globalRepository("magneto"))

fppc.Environment.setConfigFile('magneto.cfg')

fsi_tb = fsi(dim=2, orderU=2, orderP=1, orderGeo=1)
fsi_tb.init()
fsi_tb.printAndSaveInfo()


fsi_tb.reset_executionTime()

#Add collision force 
fsi_tb.addMagnetoTorqueModelFSI()
fsi_tb.addMagnetoTroqueResModelFSI()

fsi_tb.startTimeStep()

while not fsi_tb.timeStepBase().isFinished():
    
    # min_etaq = quality.etaQ(fsi_tb.mesh()).min()
    
    # if min_etaq < 1.0:
    #     remesh_toolbox(f, hclose, hfar, None)
    #     f.addContactForceModel()
    #     f.addContactForceResModel()
 


        # nbr_remesh += 1
        # time_remesh.append(f.time())
    
  
    if fppc.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s iteration : {}\n".format(fsi_tb.time(), fsi_tb.timeStepBase().iteration()))
        #print("  -- mesh quality: {}s\n".format(min_etaq))
        print("============================================================\n")
    
    fsi_tb.solve()
    fsi_tb.exportResults()

    fsi_tb.updateTimeStep()
