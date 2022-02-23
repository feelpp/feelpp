import sys, feelpp, feelpp.mor
import pandas as pd
import time

e = feelpp.Environment(sys.argv,opts=feelpp.mor.makePbdwOptions(""))
p = feelpp.mor.Pbdw(name="toolboxmor-office228-residual-np1", dbId="a57b61ca-d6ed-45f6-909b-71609bea6c2d", dbLoad=3)
sn = p.sensorNames()
#print(sn)

csv = pd.read_csv("/home/u2/hild/gits/ibat/research/ibat/ibat/cpp/pbdw/cases/bureau228/meraki_results.csv")
csv.loc[:,sn] += 273.15
# print(csv.loc[[0],sn])

# temp = p.outputs(csv.loc[0,sn[1:]],sn[1:])
# print(temp)

# temp = p.outputsWithout(csv.loc[0,sn[:3]+sn[4:]],[sn[3]])
# print(temp)

# sn2=sn.copy()
# sn.remove("MT10-04")
outputs = []
# outputs.append(p.outputs(csv.loc[0,sn]))
# for i in range(len(sn)):
#     # print(sn[i])
#     # print(csv.loc[[0],sn[:i]+sn[i+1:]])
#     s = time.time()
#     outputs.append(p.outputsWithout(csv.loc[0,sn[:i]+sn[i+1:]],[sn[i]]))
#     print(time.time()-s)
#     # outputs.append(p.outputsWithout(csv.loc[0,sn[:i]+sn[i+1:]],[sn[i],"MT10-04"]))
#     # outputs.append(p.outputs(csv.loc[0,sn[:i]+sn[i+1:]],sn[:i]+sn[i+1:]))
#     # print(p.outputNames())
#     # print(temp)
# out = pd.DataFrame(outputs,columns=p.outputNames(),index=pd.Index(sn, name='removed sensor'))

per = 168
correct=False
for j in range(len(sn)):
    outputs = []
    s = time.time()
    for i in range(-per,0):
        outputs.append(p.outputs(csv.loc[csv.shape[0]+i,sn[:j]+sn[j+1:]],sn[:j]+sn[j+1:], correct))
        # outputs.append(p.outputs(csv.loc[csv.shape[0]+i,sn]))
    print("time: "+str(time.time()-s)+"s")

    out = pd.DataFrame(outputs, columns=p.outputNames(), index=pd.Index(csv.loc[csv.shape[0]-per:,"date_format"]))
    out["Top"] = (out["average_temp"]+out["average_temp_surf"])/2
    out["average_sensor"] = out.loc[:,sn[:3]+sn[4:10]].mean(axis=1)
    out["average_data"] = csv.set_index("date_format").loc[out.index,sn[:3]+sn[4:10]].mean(axis=1)
    out = out.join(csv.set_index("date_format").loc[out.index,out.columns[:10]].add_prefix("data_"))
    out.loc[:,:] -= 273.15
    out["error_rel"] = (out["average_temp"]-out["average_data"]).abs()/out["average_data"]

    name = "out_"+sn[j]
    if not correct:
        name+="_nocorrect"
    name+=".csv"
    out.to_csv(name)
print(out.head())
