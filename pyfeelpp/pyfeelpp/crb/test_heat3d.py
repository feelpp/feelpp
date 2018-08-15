import online
online.initOnline()
p=online.loadPlugin('heat3d')
Dmu=p.parameterSpace()
S=Dmu.sampling()
S.sampling(100,'random')
for mu in S.getVector():
    r=p.run(mu)
    print("mu=",r.parameter(),"s=",r.output(), " error=",r.errorBound())

