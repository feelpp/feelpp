import feelpp.mor.online
online.initOnline()
p=online.loadPlugin('heat3d')
Dmu=p.parameterSpace()
S=Dmu.sampling()
S.sample(100,'random')
for mu in S.getVector():
    r=p.run(mu)
    print("mu=",r.parameter(),"s=",r.output(), " error=",r.errorBound())

