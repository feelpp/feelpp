import sys

import feelpp.core as fppc
import feelpp.toolboxes as fppt as fppt

e=fppc.Environment(sys.argv,opts=tb.toolboxes_options("hdg"))

#from pyfeelpp import discr,ts,filters
from feelpp.toolboxes.hdg import *

f=hdgpoisson(dim=2,order=1)
f.init()
f.printAndSaveInfo()
print("parameters:",len(f.modelProperties().parameters()))
print("materials:",len(f.modelProperties().materials()))
print("materials:",f.modelProperties().materials())
print("omega:",f.modelProperties().materials()["omega"])
print("omega has k?:",f.modelProperties().materials()["omega"].hasProperty("k"))
for k,v in f.modelProperties().materials():
    print(k," ", v)
f.modelProperties().materials().at("omega").setProperty("k","10")
f.modelProperties().materials().at("omega").setProperty("cond","10")
print("xx",f.modelProperties().materials().at("omega"))
    
f.solve()
f.exportResults()
