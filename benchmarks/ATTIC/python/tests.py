#! /usr/bin/python

import xmlParser
import validation
import error
import displayGraphics
from pyx import color,text,unit,canvas

from numpy import log

################################################################################
#				Main program
################################################################################




reponse = xmlParser.parse_xml_response("essai_xml_response.xml")
print(validation.validation_code(reponse,"essai_xml.xml","laplacian1DP1",0.32))
#print(validation.validation_code_norm(reponse,"essai_xml.xml","laplacian1DP1","H1",0.05))
data=[[1,2,20],[1,4,6.78]]
func="x**2"
#print(error.validation_all(data,func,0.1))

#data=datas[0]
#func=funcs[0]
data2=error.changedata(data,func)
func2=error.changefunc(func)
data3=error.slope_calc(data,func)

#print(data2)
#print(func2)

##V2 affichage en utilisant nos fonctions

#h=displayGraphics.create_graph_XlogYlin(12,1,20,1,20,"h","Error")






#displayGraphics.plotPointsXlinYlog(h,data,"Data Points")
#displayGraphics.plotLinesXlinYlin(h,data[0],error.peval(data2[0],data3),color.rgb.red,"Interpolated Error")
#displayGraphics.plotFunction(h,func2,color.rgb.blue,"Theoretical Error")


#x1,x2 = h.pos(23,10)
#tbox0 =text.text(x1, x2, r"Theoretical error : \ \ \ \ \ \  ")
#tbox1  = text.text(x1, x2-20*unit.x_pt, r"\left\| . \right\|_{H^{1}} \sim "+func,[text.mathmode, text.valign.bottom])
#c=canvas.canvas()
#c.insert(tbox0)
#c.insert(tbox1)


#h.insert(c)
#h.writePDFfile("plotL2")

#print(error.validation(data,func,0.1))











###########################
#On parse pour la norme H1
###########################
#data,func = xmlParser.parse_xml_result("essai_xml.xml", prognamelist[0], dimlist[0], plist[0], normlist[1], paramlist)


#data2=error.changedata(data,func)
#func2=error.changefunc(func)
#data3=error.slope_calc(data,func)

#print(data2)
#print(func2)

#V2 affichage en utilisant nos fonctions
#c=canvas.canvas()
#h=displayGraphics.create_graph_XlogYlin(8,1,10,1,10,"h","Error")
#displayGraphics.plotPointsXlinYlog(h,data,"Data Points")
#displayGraphics.plotLinesXlinYlin(h,data[0],error.peval(data2[0],data3),color.rgb.red,"Interpolated Error")
#displayGraphics.plotFunction(h,func2,color.rgb.blue,"Theoretical Error")
#h.insert(c)
#h.writePDFfile("plotH1")


#print(error.validation(data,func,0.1))