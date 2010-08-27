#! /usr/bin/python

# -*- mode: python -*-
#
#  This file is part of the Feel library
#
#  Author(s): Florent Vielfaure <florent.vielfaure@gmail.com>
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
# \date 2009-04-07
#

## \package util
#  This package contains classes to handle parameters and outputs of the c++ codes

import sys
import string
import os
import scipy
import random

def nbrToSci(nbr):
		return "%.5e" % (nbr)

def strToSci(str):
		return "%.5e" % (float(str))

def intervalLog(Vmin,Vmax,nbVal):
		v=[]
		for i in scipy.arange(0,nbVal):
				v.append(Vmax/(2**i))
		return v

		# print "Vmin,Vmax=",Vmin, Vmax
		# return (scipy.ceil(100*scipy.exp(scipy.linspace(scipy.log(Vmin),
		# 												 scipy.log(Vmax),
		# 												 nbVal)))/100)

def intervalLin(Vmin,Vmax,nbVal):
		v=[]
		for i in scipy.arange(0,nbVal):
				v.append(Vmax/(2**i))
		return v
#		return scipy.arange(Vmin,Vmax,(Vmax-Vmin)/nbVal).tolist()[1:]

def pickInInterval(Vmin,Vmax,nbVal,nbPos,method):
		if method==0:
				return intervalLin(Vmin,Vmax,nbPos)[:nbVal]
				#return random.sample(intervalLin(Vmin,Vmax,nbPos),nbVal)
		if method==1:
				return intervalLog(Vmin,Vmax,nbPos)[:nbVal]
				#return random.sample(intervalLog(Vmin,Vmax,nbPos),nbVal)

class Parameter:
	def __init__(self,attrNames="",attrValues="",values=""):
		self.__name=""
		for ind in range(0,len(attrNames)):
			if attrNames[ind] == "name":
				self.__name=attrValues[ind]
		self.__type=0
		for ind in range(0,len(attrNames)):
			if attrNames[ind] == "type":
				if attrValues[ind] == "discrete":
					self.__type=1
				else:
					self.__type=2
		self.__cmdName=self.__name
		for ind in range(0,len(attrNames)):
			if attrNames[ind] == "cmd_name":
				self.__cmdName=attrValues[ind]
		self.__latex=""
		for ind in range(0,len(attrNames)):
			if attrNames[ind] == "latex":
				self.__latex=attrValues[ind]
		self.__values=values

	def getAttrNames(self):
		AttrNamesTab=[]
		AttrNamesTab.append("name")
		if (self.__type!=0):
			AttrNamesTab.append("type")
		if (self.__cmdName!=""):
			AttrNamesTab.append("cmdName")
		if (self.__latex!=""):
			AttrNamesTab.append("latex")
		if (self.__values!=""):
			AttrNamesTab.append("values")
		return AttrNamesTab

	def getAttrValues(self):
		AttrValuesTab=[]
		AttrValuesTab.append(self.__name)
		if (self.__type!=0):
			if (self.__type==1):
				AttrValuesTab.append("discrete")
			else:
				AttrValuesTab.append("continuous")
		if (self.__cmdName!=""):
			AttrValuesTab.append(self.__cmdName)
		if (self.__latex!=""):
			AttrValuesTab.append(self.__latex)
		if (self.__values!=""):
			AttrValuesTab.append(self.__values)
		return AttrValuesTab

	def getValues(self):
		return self.__values

	def getName(self):
		return self.__name

	def getCmdName(self):
		return self.__cmdName;

class Output(Parameter):
	def __init__(self, attrNames, attrValues, dependencies, funcs):
		Parameter.__init__(self,attrNames, attrValues)
		self.__dependencies=dependencies
		self.__funcs=funcs

	def getDependencies(self):
		return self.__dependencies

	def getFuncs(self):
		return self.__funcs

# Local Variables:
# indent-tabs-mode: t
# End:
