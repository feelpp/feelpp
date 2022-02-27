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
	  v.append(1) # meshsize times 1e-3 (default: 1)
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
unifQ = Uniform(0.,10.e6)
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


#disttemp = Distribution(unifkIC)
#comptemp[0] = "k2"
#disttemp.setDescription(comptemp)
#disttemp.setName("random variable k2")
#inpCollection[0] = disttemp

disttemp = Distribution(unifr)
comptemp[0] = "r"
disttemp.setDescription(comptemp)
disttemp.setName("random variable r")
inpCollection[3] = disttemp

#disttemp = Distribution(unifr)
#comptemp[0] = "r2"
#disttemp.setDescription(comptemp)
#disttemp.setName("random variable r2")
#inpCollection[0] = disttemp


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

#eads_mfem = NumericalMathFunction("opuseadsmfem")
eads_pfem = NumericalMathFunction("opuseadspfem")
enlarge = NumericalMathFunction(FUNC())
#modelfull = NumericalMathFunction(eads_mfem,enlarge)
modelfull = NumericalMathFunction(eads_pfem,enlarge)
#reduced = NumericalMathFunction(FUNC2())
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

# Polynomial chaos expansion
##########

dim = inp.getDimension()

polyColl = PolynomialFamilyCollection(dim)

legendreFamily = LegendreFactory()
polyColl[0] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
polyColl[1] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
polyColl[2] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
polyColl[3] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
polyColl[4] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
#polyColl[5] = OrthogonalUniVariatePolynomialFamily(legendreFamily)
#polyColl[6] = OrthogonalUniVariatePolynomialFamily(legendreFamily)

multivariateBasis = OrthogonalProductPolynomialFactory(polyColl,EnumerateFunction(dim))
#distributionMu = multivariateBasis.getMeasure()

maximumConsideredTerms = 500
mostSignificant = 10
significanceFactor = 1.e-6
truncatureBasisStrategy = CleaningStrategy(OrthogonalBasis(multivariateBasis),maximumConsideredTerms, mostSignificant, significanceFactor, True)
sampleSize = 500
evaluationCoeffStrategy = LeastSquaresStrategy(LHSExperiment(sampleSize))
polynomialChaosAlgorithm = FunctionalChaosAlgorithm(model, Distribution(myDist), AdaptiveStrategy(truncatureBasisStrategy),ProjectionStrategy(evaluationCoeffStrategy))

polynomialChaosAlgorithm.run()
polynomialChaosResult = polynomialChaosAlgorithm.getResult()

coefficients = polynomialChaosResult.getCoefficients()
metaModel = polynomialChaosResult.getMetaModel()

sizeX = 100
Xsample = myDist.getNumericalSample(sizeX)
modelSample = model(Xsample)
metaModelSample = metaModel(Xsample)

sampleMixed = NumericalSample(sizeX,2)
for i in range(sizeX):
  sampleMixed[i][0] = modelSample[i][0]
  sampleMixed[i][1] = metaModelSample[i][0]

legend = str(sizeX) + " realizations"
comparisonCloud = Cloud(sampleMixed, "blue", "fsquare", legend)
graphCloud = Graph('Polynomial chaos expansion', 'model', 'metamodel', True, 'topleft')
graphCloud.addDrawable(comparisonCloud)
#Show(graphCloud)
graphCloud.draw('PCE_ModelsComparison')

S1tilde = FunctionalChaosRandomVector(polynomialChaosResult)
print "mean=", S1tilde.getMean()[0]
print "variance=", sqrt(S1tilde.getCovariance()[0,0])
for i in range(dim):
  print "Sobol index" , i,"=", S1tilde.getSobolIndex(i)

S1sample = S1tilde.getNumericalSample(100000)
S1tildequantile = []
S1quantile=[]
for i in range(10):
  S1tildequantile.append(S1sample.computeQuantile(i/10.)[0])
  S1quantile.append(modelSample.computeQuantile(i/10.)[0])

quantileMixed = NumericalSample(10,2)
for i in range(10):
  quantileMixed[i][0] = S1quantile[i]
  quantileMixed[i][1] = S1tildequantile[i]

legend = " Quantile Comparison"
comparisonquantileCloud = Cloud(quantileMixed, "blue", "fsquare", legend)
graphquantileCloud = Graph('Polynomial chaos expansion', 'model', 'metamodel', True, 'topleft')
graphquantileCloud.addDrawable(comparisonquantileCloud)
#Show(graphCloud)
graphquantileCloud.draw('PCE_QuantilesComparison')




