import sys, feelpp, feelpp.mor
import pandas as pd

e = feelpp.Environment(sys.argv,opts=feelpp.mor.makePbdwOptions(""))
p = feelpp.mor.Pbdw(name="toolboxmor-office228-residual-np1", dbId="a57b61ca-d6ed-45f6-909b-71609bea6c2d", dbLoad=3)

csv = pd.read_csv("meraki_results2.csv")
csv.loc[:,csv.columns.str.startswith('MT')] += 273.15

sn = p.sensorNames()
temp = p.outputs(csv.loc[0,sn[4:]],sn[4:])
print(temp)
