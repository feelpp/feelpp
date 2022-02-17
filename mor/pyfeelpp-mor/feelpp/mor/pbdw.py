import sys, feelpp, feelpp.mor
import pandas as pd

e = feelpp.Environment(sys.argv,opts=feelpp.mor.makePbdwOptions(""))
p = feelpp.mor.Pbdw(name="toolboxmor-office228-residual-np1", dbId="a57b61ca-d6ed-45f6-909b-71609bea6c2d", dbLoad=3)
sn = p.sensorNames()
#print(sn)

csv = pd.read_csv("meraki_results2.csv")
csv.loc[:,csv.columns.str.startswith('MT')] += 273.15

# temp = p.outputs(csv.loc[0,sn[1:]],sn[1:])
# print(temp)

# temp = p.outputsWithout(csv.loc[0,sn[:3]+sn[4:]],[sn[3]])
# print(temp)

sn2=sn.copy()
sn.remove("MT10-04")
outputs = []
for i in range(len(sn)):
    # print(sn[i])
    # print(csv.loc[[0],sn[:i]+sn[i+1:]])
    outputs.append(p.outputsWithout(csv.loc[0,sn[:i]+sn[i+1:]],[sn[i],"MT10-04"]))
    # outputs.append(p.outputs(csv.loc[0,sn[:i]+sn[i+1:]],sn[:i]+sn[i+1:]))
    # print(p.outputNames())
    # print(temp)
out = pd.DataFrame(outputs,columns=p.outputNames(),index=pd.Index(sn, name='removed sensor'))
out.to_csv("out.csv")
print(out.loc[:,sn2])
