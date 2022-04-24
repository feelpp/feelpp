import feelpp.mor.online
online.initOnline()
p=online.loadPlugin('heat3d')
D=p.parameterSpace()
s=D.sampling()
s.sampling(3,'random')
print("s1=",s.getVector())
s.sampling(3,'random')
for v in s.getVector():
    print("v=",v)
