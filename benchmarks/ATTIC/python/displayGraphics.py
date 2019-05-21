# -*- mode: python -*-
#
#  This file is part of the Feel library
#
#  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
#             Monzat Christophe <christophe.monzat@gmail.com>
#        Date: 2009-04-07
#
#   Copyright (C) 2009
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Universite Joseph Fourier (Grenoble I)
#
# \file codeCalculParaview.py
# \author Florent Vielfaure <florent.vielfaure@gmail.com>
# \author Monzat Christophe <christophe.monzat@gmail.com>
# \date 2009-04-07
#

## \package displayGraphics
#  This package contains methods to plot datas on graphs

from pyx import *
from pyx.deco import barrow,earrow
from pyx.style import linewidth, linestyle
from pyx.graph.axis import painter, tick
from pyx.graph.axis import *

import decimal
from scipy import *
#from scipy.io import array_import
from numpy import *

GRAPH_LIN_LIN=0
GRAPH_SEMI_LOG_X=1
GRAPH_SEMI_LOG_Y=2
GRAPH_LOG_LOG=3

class graphBuilder:

		def __init__(self,title,type,graphL,titleXaxis,titleYaxis):
				# Set properties of the defaulttexrunner, e.g. switch to LaTeX.
				try:
						print "[Validation] info : setting latex mode..."
						text.set(mode="latex")
				except RuntimeError:
						print "[Validation] info : latex mode already launched..."

				self.__title=title
				self.__c=canvas.canvas()
				self.__graphL=graphL
				if graphL==0:
						length=8
				else:
						length=13.5

				if type == GRAPH_LIN_LIN:
						self.__g=self.__c.insert(graph.graphxy(width=length,height=8,key=graph.key.key(),
                                                       x=graph.axis.linear(title=titleXaxis),
                                                       y=graph.axis.linear(title=titleYaxis)))
				elif type == GRAPH_SEMI_LOG_X:
						self.__g=self.__c.insert(graph.graphxy(width=length,height=8,key=graph.key.key(),
                                                       x=graph.axis.log(title=titleXaxis),
                                                       y=graph.axis.linear(title=titleYaxis)))
				elif type == GRAPH_SEMI_LOG_Y:
						self.__g=self.__c.insert(graph.graphxy(width=length,height=8,key=graph.key.key(),
                                                       x=graph.axis.linear(title=titleXaxis),
                                                       y=graph.axis.log(title=titleYaxis)))
				elif type == GRAPH_LOG_LOG:
						self.__g=self.__c.insert(graph.graphxy(width=length,height=8,key=graph.key.key(),
                                                       x=graph.axis.log(title=titleXaxis),
                                                       y=graph.axis.log(title=titleYaxis)))
				else:
						raise ValueError, "[displayGraphics] error : the type %d is not valid" % type

		def addPoints(self,data,name):
				self.__g.plot([graph.data.list(array([data[0],data[1]]).transpose(),
                                               title=name,
                                               x=1,y=2,)],
                              [graph.style.symbol()])

		def addLines(self,data,name,col):
				self.__g.plot([graph.data.list(array([data[0],data[1]]).transpose(),
                                               title=name,
                                               x=1,y=2,)],
                              [graph.style.line([col])])

		def addFunction(self,func,name,col):
				self.__g.plot(graph.data.function(func,title=name),
                              [graph.style.line([col])])

		def addTriangle(self,Tslope,Pslope):
				if self.__graphL==1:
						return
				self.__c.stroke(path.rect(8.5, 3.25, 5, 4.75))

				c1={'x':11,'y':6.8}
				c2={'x':11,'y':4.5}
				YTslope=2.0
				if float(Tslope)>0:
					XTslope=YTslope/float(Tslope)
				else:
					XTslope=0
				YPslope=2.0
				if float(Pslope)>0:
					XPslope=YPslope/float(Pslope)
				else:
					XPslope=0
				# Theoritical error triangle
				self.__c.stroke(path.line(c1['x']-XTslope/2,
								          c1['y']-YTslope/2,
								          c1['x']+XTslope/2+0.1,
								          c1['y']-YTslope/2))

				self.__c.stroke(path.line(c1['x']-XTslope/2,
								          c1['y']-YTslope/2,
								          c1['x']+XTslope/2,
								          c1['y']+YTslope/2),[color.rgb.red])

				self.__c.stroke(path.line(c1['x']+XTslope/2,
								          c1['y']+YTslope/2,
								          c1['x']+XTslope/2,
								          c1['y']-YTslope/2-0.1))

				self.__c.stroke(path.line(c1['x']-XTslope/2,
								          c1['y']-YTslope/2,
								          c1['x']-XTslope/2,
								          c1['y']-YTslope/2-0.1))

				self.__c.stroke(path.line(c1['x']+XTslope/2+0.1,
								          c1['y']+YTslope/2,
								          c1['x']+XTslope/2,
								          c1['y']+YTslope/2))

				self.__c.text(c1['x']+XTslope/2+0.1,
					          c1['y'],
					          Tslope)

				# Practical error triangle
				self.__c.stroke(path.line(c2['x']-XPslope/2,
								          c2['y']-YTslope/2,
								          c2['x']+XPslope/2+0.1,
								          c2['y']-YTslope/2))

				self.__c.stroke(path.line(c2['x']-XPslope/2,
								          c2['y']-YTslope/2,
								          c2['x']+XPslope/2,
								          c2['y']+YTslope/2),[color.rgb.blue])

				self.__c.stroke(path.line(c2['x']+XPslope/2,
								          c2['y']+YTslope/2,
								          c2['x']+XPslope/2,
								          c2['y']-YTslope/2-0.1))

				self.__c.stroke(path.line(c2['x']-XPslope/2,
								          c2['y']-YTslope/2,
								          c2['x']-XPslope/2,
								          c2['y']-YTslope/2-0.1))

				self.__c.stroke(path.line(c2['x']+XPslope/2+0.1,
								          c2['y']+YTslope/2,
								          c2['x']+XPslope/2,
								          c2['y']+YTslope/2))

				self.__c.text(c2['x']+XPslope/2+0.1,
					          c2['y'],
					          Pslope)


		def addInfo(self,theoritical,practical,valid,precision):
				if self.__graphL==1:
						return

				self.__c.stroke(path.rect(8.5, 0, 5, 2.75))

				self.__c.text( 9, 2.1,theoritical)
				self.__c.text( 9, 1.5,practical)

				tbox = 0;
				if (valid):
						tbox=text.text( 9, 1,"Validation : TRUE",[text.mathmode, text.valign.middle])
						tpath = tbox.bbox().enlarged(3*unit.x_pt).path()
						self.__c.draw(tpath, [deco.filled([color.rgb.green]), deco.stroked()])
				else:
						tbox=text.text( 9, 1,"Validation : FALSE",[text.mathmode, text.valign.middle])
						tpath = tbox.bbox().enlarged(3*unit.x_pt).path()
						self.__c.draw(tpath, [deco.filled([color.rgb.red]), deco.stroked()])
				self.__c.insert(tbox)

				self.__c.text( 9, 0.5,"Precision : "+str(precision),[text.mathmode, text.valign.middle])

		def writePDFfile(self,PDFname):
				self.__c.stroke(path.rect(0, 8.5, 13.5, 1.5))
				self.__c.text(6.75,9,self.__title,[text.size(2),text.halign.center])
				self.__c.writePDFfile(PDFname)



# Local Variables:
# indent-tabs-mode: t
# End:
