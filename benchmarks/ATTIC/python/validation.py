#! /usr/bin/python

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

## \package validation
#  This package contains tools for c++ code validation.

import xmlParser
import error
import displayGraphics
import util

import decimal
from scipy import *
#from pyx import *
#from pyx.deco import barrow,earrow
#from pyx.style import linewidth, linestyle
#from pyx.graph.axis import painter, tick
#from pyx.graph.axis import *
from pyx import color#,text,unit,canvas
import time,datetime
from numpy import log,exp


##	validation_code
#	Validate code for a set of output and param, generate graphics
#
#	@param name name of the c++ code
#	@param param set of parameters
#	@param output set of output
#	@param xmlname name of the xml file containing the result
#	@param precision precision of the validation
#	@param val_param values of the parameters specified by the user
#	@param val_depend value of the parameter depend
#
#	@return true or false for the validation
def validation_code(name,param,output,xmlname,val_param,val_depend,precision,pathD):
	result=1
	for out in output:
		if len(out.getDependencies())==0:
			print "[Validation] error : no dependencies..."
			break
		depend=out.getDependencies()[0]
		func=out.getFuncs()[0]
		data=[]
		for val_dep in val_depend[depend]:
			list_param=[]
			val=[]
			for par in param:
				if (par.getName()!=depend):
					list_param.append(par.getName())
					val.append(val_param[par.getName()][0])
				else:
					list_param.append(par.getName())
					val.append(val_dep)
					par_depend=par
			data.append(xmlParser.parse_xml_result(xmlname,name,list_param,val,out.getName()))

		data=map(float,data)
		val=map(float,val_depend[depend])
		data_temp=[val,data]

		temp=error.validation(data_temp,func,depend,precision)
		result=temp and result

		nb_decimal=len(str(precision))-1

		#Transformation des donnees pour affichage graphique
		data2=error.changedata(data_temp,func)
		func2=error.changefunc(func)
		data3 = error.fit(data2[0],data2[1])

		#On veut un certain nombre de chiffres significatifs pour les estimations
		#bo abscisse a l'origine et b1 pente
		b0=str(decimal.Decimal(str(round(exp(data3[1][0]),nb_decimal))))
		b1=str(decimal.Decimal(str(round(data3[1][1],nb_decimal))))

		#On teste si le parametre variable contient une description latex
		if (par_depend.getAttrNames().count("latex")):
				chaine=par_depend.getAttrValues()[3]
				func=func.replace(par_depend.getName(),par_depend.getAttrValues()[3])
		else:
				chaine=depend

		#On teste le nombre de log appliques a la fonction
	   	#Permet de mettre les bons titres aux axes et de rappliquer un log aux donnees si besoin
		if ((func2.count("log")-func.count("log"))==2):
				debut_log="log(log("
				fin_log="))"
				data_temp=[data_temp[0],log(data_temp[1])]
				interp="exp(chaine+**{"+b1+"})"
		else:
				debut_log="log("
				fin_log=")"
				interp=chaine+"**{"+b1+"}"

		a=float(exp(data3[1][0]))
		b=float(data3[1][1])

		graph = displayGraphics.graphBuilder(name+" (generated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S)"),3,0,"log(x)","log(y)")

		graph.addTriangle(func.split("**")[-1],b1)

		graph.addFunction("y(x)="+func.replace(par_depend.getName(),'x'),"theoritical",color.rgb.red)
		graph.addFunction("y(x)="+str(a)+"*x**"+str(b),"practical",color.rgb.blue)

		graph.addPoints(data_temp,"data")

		graph.addInfo("$"+str(out.getAttrValues()[2])+"$"+" : "+"$"
                      +str(chaine)+"$"+"$\longrightarrow$"+"$"+func.replace("**","^")+"$",
                      "Interpolated : "+"$"+str(chaine)+"$"
					  +"$\longrightarrow$"+"$"+interp.replace("**","^")+"$",
                      temp,precision)

		print("[Validation] : writing "+ pathD+out.getName()+".pdf")
		graph.writePDFfile(pathD+out.getName())
	return result


##	data_plot
#	Plot a set of datas in function of a particular parameter
#
#	@param name name of the c++ code
#	@param param set of parameters
#	@param out output we want to plot
#	@param xmlname name of the xml file containing the result
#	@param precision precision of the validation
#	@param val_param values of the parameters specified by the user
#	@param var name of the parameter we want to make vary
#	@param val_depend value of the parameter depend
#	@return the value 0 for convenience
def data_plot(name,param,out,xmlname,val_param,var,val_depend,precision,path):

	data=[]
	for val_dep in val_depend[var]:
		list_param=[]
		val=[]
		for par in param:
			if (par.getName()!=var):
				list_param.append(par.getName())
				val.append(val_param[par.getName()][0])
			else:
				list_param.append(par.getName())
				val.append(val_dep)
		data.append(xmlParser.parse_xml_result(xmlname,name,list_param,val,out.getName()))

	data=map(float,data)
	val=map(float,val_depend[var])
	data_temp=[val,data]

	diff=float(max(data_temp[1])-min(data_temp[1]))/2.0

	graph = displayGraphics.graphBuilder(name+" (generated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S)"),3,1,"log(x)","log(y)")

	graph.addPoints(data_temp,"$"+out.getAttrValues()[2]+"$")
	graph.addLines(data_temp,"Lines",color.rgb.red)

	print("[Validation] : writing "+path+out.getName()+".pdf")
	graph.writePDFfile(path+out.getName())

	return 0


##	data_interpol
#	Take a set of datas and make an interpolation in function of var
#
#	@param name name of the c++ code
#	@param param set of parameters
#	@param out output we want to plot
#	@param xmlname name of the xml file containing the result
#	@param val_param values of the parameters specified by the user
#	@param var name of the parameter we want to make vary
#	@param val_depend value of the parameter depend
#	@param type_int type of the function followed by out
#   <ul>
#			<li>1 : a*x^p</li>
#			<li>2 : exp(a*x^p)</li>
#   </ul>
#	@return the value 0 for convenience
def data_interpol(name,param,out,xmlname,val_param,var,val_depend,precision,type,path):

	# Set properties of the defaulttexrunner, e.g. switch to LaTeX.
	text.set(mode="latex")

	data=[]
	for val_dep in val_depend[var]:
		list_param=[]
		val=[]
		for par in param:
			if (par.getName()!=var):
				list_param.append(par.getName())
				val.append(val_param[par.getName()][0])
			else:
				list_param.append(par.getName())
				val.append(val_dep)
				par_depend=par
		data.append(xmlParser.parse_xml_result(xmlname,name,list_param,val,out.getName()))
	data=map(float,data)
	val=map(float,val_depend[var])
	data_temp=[val,data]

	if (type==1):
		func=var
	elif (type==2):
		func="exp("+var+")"


	#Transformation des donnees pour affichage graphique
	data2=error.changedata(data_temp,func)
	a,b = error.fit(data2[0],data2[1])
	func2=error.changefunc(func)


	#On teste si le parametre variable contient une description latex
	if (par_depend.getAttrNames().count("latex")):
		chaine="$"+par_depend.getAttrValues()[3]+"$"
		func=func.replace("nu",par_depend.getAttrValues()[3])

		if (type==1):
			func_th="a*"+par_depend.getAttrValues()[3]+"**b"
			func=str(exp(b[0]))+"*"+var+"**"+str(b[1])
		elif (type==2):
			func_th="exp(a*"+par_depend.getAttrValues()[3]+"**b)"
			func="exp("+str(exp(b[0]))+"*"+var+"**"+str(b[1])+")"
	else:
		chaine=var
		if (type==1):
			func_th="a*"+chaine+"**b"
			func=str(exp(b[0]))+"*"+var+"**"+str(b[1])
		elif (type==2):
			func_th="exp(a*"+chaine+"**b)"
			func="exp("+str(exp(b[0]))+"*"+var+"**"+str(b[1])+")"
	print(b)
	b0=str(decimal.Decimal(str(round(exp(b[0]),3))))
	b1=str(decimal.Decimal(str(round(b[1],3))))
	print(b0,b1)

	graph = displayGraphics.graphBuilder(name+" (generated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S)"),0,1,"x","y")

	graph.addPoints(data_temp,"$"+out.getAttrValues()[2]+"$")
	graph.addFunction(func.replace(var,"x"),color.rgb.red,"Interpolated function",color.rgb.red)

	#x1,x2 = g.pos(min(data_temp[0])+0.1,max(data_temp[1]))
	#tbox0 =text.text(x1, x2+50*unit.x_pt,"$"+out.getAttrValues()[2]+"$"+" : "+chaine+"$\longrightarrow$"+"$"+func_th.replace("**","^")+"$")
	#tbox1 = text.text(x1, x2+30*unit.x_pt,"a : "+b0,[text.mathmode, text.valign.bottom])
	#tbox2  = text.text(x1, x2+10*unit.x_pt,"b : "+b1,[text.mathmode, text.valign.bottom])
	#c=canvas.canvas()
	#c.insert(tbox0)
	#c.insert(tbox1)
	#c.insert(tbox2)
	#g.insert(c)

	print("[Validation] : writing "+path+out.getName()+".pdf")

	graph.writePDFfile(path+out.getName())

	return 0

# Local Variables:
# indent-tabs-mode: t
# End:
