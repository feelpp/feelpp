import sys
from pyfeelpptoolboxes.hdg import *

e=core.Environment(sys.argv,opts=hdg_poisson_options())


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
