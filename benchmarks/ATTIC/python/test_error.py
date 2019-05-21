#!/usr/bin/python

import error
import displayGraphics
from pyx import *
from pyx.deco import barrow,earrow
from pyx.style import linewidth, linestyle
from pyx.graph.axis import *
from pyx.graph.axis import painter, tick

from scipy import *
from scipy.optimize import leastsq
#from scipy.io import array_import

from numpy import *

data=[[2,3,4],[7.0,27.54,74]]
func='x**3'
#func='exp(x**2)'
data2=error.changedata(data,func)
func2=error.changefunc(func)
#funcexp2=error.changefunc(funcexp)
print(data2)
print(func2)

data3=error.slope_calc(data,func)
print(data3)

#V1 affichage sans utiliser nos fonctions
c=canvas.canvas()
g=c.insert(graph.graphxy(width=8,
                  x=graph.axis.log(min=1, max=10),
                  y=graph.axis.linear(min=1, max=10)))

g.plot([graph.data.list(array([data[0],log(data[1])]).transpose(), x=1,y=2,)],[graph.style.symbol()])

g.plot([graph.data.list(array([data[0],error.peval(data2[0],data3)]).transpose(), x=1, y=2,)],[graph.style.line([color.rgb.red])])

g.plot([graph.data.function('y(x)='+func2)],[graph.style.line([color.rgb.blue])])


g.writePDFfile("plot")

print(error.validation(data,func,0.01))


#V2 affichage en utilisant nos fonctions
h=displayGraphics.create_graph_XlogYlin(8,1,10,1,10)
displayGraphics.plotPointsXlinYlog(h,data)
displayGraphics.plotLinesXlinYlin(h,data[0],error.peval(data2[0],data3),color.rgb.red)
displayGraphics.plotFunction(h,func2,color.rgb.blue)
h.writePDFfile("plot2")
