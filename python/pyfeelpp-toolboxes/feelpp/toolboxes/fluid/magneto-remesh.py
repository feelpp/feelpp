import sys
import feelpp.core as fppc 
import feelpp.toolboxes.core as tb
#import fppc.interpolation as I
from feelpp.toolboxes.fluid import *
import json

sys.argv = ['magneto-remesh']
e = fppc.Environment(sys.argv, opts=tb.toolboxes_options("fluid"),config=fppc.globalRepository("magneto-remesh"))
fppc.Environment.setConfigFile('/nvme0/vanlandeghem/feelpp-test/toolboxes/fluid/cases/head.cfg')

f = fluid(dim=2, orderVelocity=2, orderPressure=1)
f.init()
f.printAndSaveInfo()

f.startTimeStep()

# add model
f.reset_Data()
f.addMagnetoTorque()
f.addMagnetoTorqueRes()


freq_remesh = 5 # frequence donné dans le json
freq = 0

while not f.timeStepBase().isFinished():
    
    freq += 1
    
    if freq == freq_remesh:
        # add models
        freq = 0
        f.addMagnetoTorque()
        f.addMagnetoTorqueRes()

    if fppc.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s iteration : {}\n".format(f.time(), f.timeStepBase().iteration()))
        print("============================================================\n")
    
    f.solve()
    f.exportResults()
    f.updateTimeStep()
    
f.write_Data()