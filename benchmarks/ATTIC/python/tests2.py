#! /usr/bin/python

import xmlParser
from xml.dom import minidom
import error
import displayGraphics
import validationbis
import util
from pyx import *
from pyx.deco import barrow,earrow
from pyx.style import linewidth, linestyle
from pyx.graph.axis import *
from pyx.graph.axis import painter, tick

from scipy import *
from scipy.optimize import leastsq
#from scipy.io import array_import

from numpy import *

################################################################################
#				Main program
################################################################################
reponse = xmlParser.parse_xml_response("essai_xml_response2.xml")
reponse2 = xmlParser.parse_xml_response("essai_xml_response3.xml")

########################################################
#On parse l'ensemble des datas et funcs et on les stocke
########################################################


progname="laplacian"
path="./"

precision=0.01
val_param={"dim":["1"],"order":["1"],"beta":["0.1"],"nu":["0.2"],"h":["0.5"]}
val_depend={"h":["0.5","0.75","0.95"],"nu":["0.2","0.4","0.8"],"beta":["0.1","0.4","0.7","0.9","1.4"]}

data=[[0.5,0.75,1],[-1,2,3]]
func="exp(x**2)"
depend="x"
precision=0.1

#print(error.slope_calc(data,func))
#print(error.validation(data,func,depend,precision))

#print(validationbis.validation_code(reponse[0],reponse[1],reponse[2],"essai_xml2.xml",val_param,val_depend,precision,path))

#print(validationbis.validation_code2(reponse[0],reponse[1],reponse[2],"essai_xml2.xml",val_param,val_depend,precision))

validationbis.data_plot(reponse[0],reponse2[1],reponse2[2][0],"essai_xml2.xml",val_param,"beta",val_depend,precision,path)

#validationbis.data_interpol(reponse[0],reponse2[1],reponse2[2][0],"essai_xml2.xml",val_param,"beta",val_depend,precision,1,path)
#validationbis.data_interpol(reponse[0],reponse[1],reponse[2][0],"essai_xml2.xml",val_param,"h",val_depend,precision,1,path)
