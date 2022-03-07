# -*- coding: iso-8859-15 -*-
from openturns import *
from math import *

########Uncertain variables
#############################################################
#k  -> thermal conductivity -> [0.2,150]
#D   -> inflow rate          -> [0.0005,0.01]
#Q   -> heat source          -> [0,10^6]
#r -> thermal conductance  -> [0.1,100]
#ea  -> thickness            -> [0.0025,0.05]


########Output variables
#############################################################
#s1 -> mean temperature of the hottest IC
#s2 -> mean temperature of the air at the outlet

class FUNC(OpenTURNSPythonFunction):
	def __init__(self):
	  OpenTURNSPythonFunction.__init__(self, 5, 7)
	def f(self,X):
	  v = []
	  for i in range(len(X)):
	    v.append(X[i])
	  v.append(1.) # meshsize times 1e-3 (default: 1)
	  v.append(2) # integer : polynomial degree for the temperature in {1,2,3,4} (default: 2)
	  return v


inputVariablesC = Description(2)
inputVariablesC[0] = "S1"
inputVariablesC[1] = "S2"
outputVariableC = Description(1)
outputVariableC[0] = "Resu"
formulaC = Description(1)
formulaC[0] = "S1"
P_C = NumericalMathFunction(inputVariablesC, outputVariableC, formulaC)

#Log.Show(Log.DBG)

####distribs de base
unifQ = Uniform(0.,1.e6)
unifkIC = Uniform(0.2,150.)
unifr = Uniform(.1,100)
unifD = Uniform(.0005,.01)
unifea = Uniform(0.0025,0.05)

####input vector
comptemp = Description(1)
inpCollection = DistributionCollection(5)

disttemp = Distribution(unifkIC)
comptemp[0] = "kIC"
disttemp.setDescription(comptemp)
disttemp.setName("random variable kIC")
inpCollection[0] = disttemp

disttemp = Distribution(unifD)
comptemp[0] = "D"
disttemp.setDescription(comptemp)
disttemp.setName("random variable D")
inpCollection[1] = disttemp

disttemp = Distribution(unifQ)
comptemp[0] = "Q"
disttemp.setDescription(comptemp)
disttemp.setName("random variable Q")
inpCollection[2] = disttemp

disttemp = Distribution(unifr)
comptemp[0] = "r"
disttemp.setDescription(comptemp)
disttemp.setName("random variable r")
inpCollection[3] = disttemp

disttemp = Distribution(unifea)
comptemp[0] = "ea"
disttemp.setDescription(comptemp)
disttemp.setName("random variable ea")
inpCollection[4] = disttemp

######copula

R = CorrelationMatrix(inpCollection.getSize())

aCopula = IndependentCopula(inpCollection.getSize())
aCopula.setName("Copula of the random input vector")

myDist = ComposedDistribution(inpCollection,Copula(aCopula))

#input and output random vector definition

inp = RandomVector(Distribution(myDist))

eads_mfem = NumericalMathFunction("opuseadsmfem")
enlarge = NumericalMathFunction(FUNC())
modelfull = NumericalMathFunction(eads_mfem,enlarge)
model = NumericalMathFunction(P_C,modelfull)

outP = RandomVector(model,inp)

#a = NumericalPoint(2)
#b = P_C(a)
#inP = NumericalPoint(5)
#inP[0] = 10   # kIC : thermal conductivity (default: 2)
#inP[1] = 5e-3 # D : fluid flow rate (default: 5e-3)
#inP[2] = 10.e6  # Q : heat flux (default: 1e6)
#inP[3] = 100  # r : conductance (default: 100)
#inP[4] = 4e-3 # ea : length air flow channel (default: 4e-3)

#model(inP)

#Sampling 

#samp = outP.getNumericalSample(100)

sizeX = 1000
Xsample = myDist.getNumericalSample(sizeX)
modelSample = model(Xsample)

# save study

study_sample = Study()

filename = "./study_sample.xml"
study_sample.setStorageManager(XMLStorageManager(filename))

study_sample.add("inpdistribution",myDist)
#study_sample.add("inpvector",inp)
#study_sample.add("modeltotmfem",model)
#study_sample.add("outvector",outP)
study_sample.add("modeltotmfemsample", modelSample)
study_sample.add("inputsample", Xsample)

study_sample.save()


