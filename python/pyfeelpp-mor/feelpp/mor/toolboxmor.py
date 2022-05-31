import sys, os
from feelpp.toolboxes.heat import *
from feelpp.mor import *
import feelpp

o = toolboxes_options("heat")
o.add(makeToolboxMorOptions())
sys.argv = ["toolbox-mor"]
# sys.argv = ['--config-file opusheat/opusheat-heat.cfg']
e = feelpp.Environment(sys.argv, opts=o)

casefile = 'thermal-fin/2d/thermal-fin.cfg'
feelpp.Environment.setConfigFile(casefile)
name = "thermal-fin-2d"

model_path = "$cfgdir/"+os.path.splitext(os.path.basename(casefile))[0] + ".json"
crb_model_properties = CRBModelProperties(worldComm=feelpp.Environment.worldCommPtr())
crb_model_properties.setup(model_path)

model_properties = feelpp.ModelProperties(worldComm=feelpp.Environment.worldCommPtr())
model_properties.setup(model_path)

outputs = crb_model_properties.outputs()
output_names = []
output_objs = []
for n, o in outputs:
    output_names.append(n)
    output_objs.append(o)

dim = 2
assert dim in [2,3]

heatBox = heat(dim=dim, order=1)
heatBox.init()
model = toolboxmor(name=name, dim=dim, time_dependent=False)
model.setFunctionSpaces( Vh=heatBox.spaceTemperature() )

def assembleDEIM(mu):
    for i in range(0,mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleRhs()

def assembleMDEIM(mu):
    for i in range(0,mu.size()):
        heatBox.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBox.updateParameterValues()
    return heatBox.assembleMatrix()

model.setAssembleDEIM(fct=assembleDEIM)
model.setAssembleMDEIM(fct=assembleMDEIM)
model.initModel()

heatBoxDEIM=heat(dim=dim,order=1)
meshDEIM = model.getDEIMReducedMesh()
heatBoxDEIM.setMesh(meshDEIM)
heatBoxDEIM.init()

def assembleOnlineDEIM(mu):
    for i in range(0,mu.size()):
        heatBoxDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBoxDEIM.updateParameterValues()
    return heatBoxDEIM.assembleRhs()

model.setOnlineAssembleDEIM(assembleOnlineDEIM)

heatBoxMDEIM=heat(dim=dim,order=1)
meshMDEIM = model.getMDEIMReducedMesh()
heatBoxMDEIM.setMesh(meshMDEIM)
heatBoxMDEIM.init()

def assembleOnlineMDEIM(mu):
    for i in range(0,mu.size()):
        heatBoxMDEIM.addParameterInModelProperties(mu.parameterName(i),mu(i))
    heatBoxMDEIM.updateParameterValues()
    return heatBoxMDEIM.assembleMatrix()

model.setOnlineAssembleMDEIM(assembleOnlineMDEIM)

model.postInitModel()
model.setInitialized(True)

[Aq, Fq] = model.getAffineDecomposition()

print("Aq")
print(Aq)
print("Fq")
print(Fq)

Dmu = model.parameterSpace()
mu = Dmu.element(True, False)

[betaA, betaF] = model.computeBetaQm(mu)

print("mu")
print(mu)
print("betaA")
print(betaA)
print("betaF")
print(betaF)

# crbmodel = crbmodel_toolboxmor_2d(model)
# crb = crb_toolboxmor_2d(crbmodel)
# crb.offline()

print("cool")
