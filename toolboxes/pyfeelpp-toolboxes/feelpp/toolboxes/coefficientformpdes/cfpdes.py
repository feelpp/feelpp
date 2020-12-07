import sys
import feelpp
import feelpp.toolboxes as tb
import feelpp.toolboxes.cfpdes as cfpdes
import pandas as pd

e=feelpp.Environment(sys.argv,opts=tb.toolboxes_options("coefficient-form-pdes","cfpdes"))         
e.setConfigFile('thermo2dCase2.cfg')   
#cd ~/Devel/feelpp.quick/toolboxes/coefficientformpdes/bratu/square/  

f=cfpdes.simulate(dim=2)
df = pd.read_csv("cfpdes.heat.measures.csv", sep=",")
df.head(10)