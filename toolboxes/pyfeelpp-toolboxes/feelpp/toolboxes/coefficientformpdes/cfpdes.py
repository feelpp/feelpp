import sys
import os
import feelpp
import feelpp.toolboxes as tb
import feelpp.toolboxes.cfpdes as cfpdes
import pandas as pd

cwd = os.getcwd()

e=feelpp.Environment(sys.argv,opts=tb.toolboxes_options("coefficient-form-pdes","cfpdes"))         
e.setConfigFile('thermo2dCase2.cfg')   

# f=cfpdes.simulate(dim=2)
print("Create cfpdes")
f=cfpdes.cfpdes(dim=2)

print("Init pb")
f.init()
print("params: ", f.modelProperties().parameters())

f.solve()
f.exportResults()

# os.chdir(cwd)
# df = pd.read_csv("cfpdes.heat.measures.csv", sep=",")
# print("\nInspect measures data:\n", df.head(10))

# get Params
print("Rerun after modifying some params")
f.addParameterInModelProperties("T0_top", 10)
f.updateParameterValues()
print("after change: ", f.modelProperties().parameters())

f.solve()
f.exportResults()

os.chdir(cwd)
df = pd.read_csv("cfpdes.heat.measures.csv", sep=",")
print("\nInspect measures data:\n", df.head(10))
