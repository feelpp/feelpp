import sys 
import feelpp as feelpp
import feelpp.toolboxes.core as tb
from feelpp.toolboxes.solid import *
import math

sys.argv = ['solid-contact']
e = feelpp.Environment(sys.argv, opts=tb.toolboxes_options("solid"),config=feelpp.globalRepository("ball_newmark"))

feelpp.Environment.setConfigFile('falling_ball.cfg')

s=solid(dim=2,orderDisp=1)
s.init() 
s.printAndSaveInfo()

# Export wall
s.init_Measures()
#s.setWalls()

#Add collision force 
s.addContactForceModel()

s.startTimeStep()

it = 0
while not s.timeStepBase().isFinished():
        
    if feelpp.Environment.isMasterRank():
        print("============================================================\n")
        print("time simulation: {}s\n".format(s.time()))
        print("============================================================\n")
    

    
    s.solve()
    s.exportResults()
    s.updateTimeStep()
    
    if (it > 0):
        s.energy()
    
    
    it += 1

s.export_Measures()





