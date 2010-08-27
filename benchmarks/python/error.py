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

## @package error
#  This package contain functions for validating input data
#
import re,string
from scipy import *
from scipy.optimize import leastsq
from numpy import *


## changedata
#  Changes the data into log(data) or log(log(data))
#  @param data the datas to transform
#  @param func the theoretical function of the error
#
#  @return log(data) or log(log(data))
def changedata(data,func):
	if re.search('exp', func):
		if (min(data[1])<=1):
			raise ValueError, "Output\\le1"
		return log(data[0]),log(log(data[1]))
	else:
		if (min(data[1])<=0):
			raise ValueError, "Output\\le0"
		return log(data[0]),log(data[1])


## changefunc
#  Changes the function into log(func) or log(log(func))
# @param func the function to transform
#
# @return log(func) or log(log(func))
def changefunc(func):
	if re.search('exp', func):

		func2='log(log('+func+'))'
		return func2;
	else:
		func2='log('+func+')'
		return func2;

## peval
#  Returns the value of a polynom p on the point x
#  This function is essential to use the residuals function
#
#  @param x the point where we evaluate the polynom
#  @param p coefficients of the polynom
#
#  @return value of the polynom
def peval(x, p):
	return p[0]+p[1]*x

## expeval
#  Returns the value of a polynom p on the point exp(x)
#  This function is essential to use the residuals function
#
#  @param x the point where we evaluate the polynom
#  @param p coefficients of the polynom
#
#  @return value of the polynom
def expeval(x,p):
	return p[0]+p[1]*exp(x)

## fit
#  Returns the values of the least square algorithm
#
#  @param x datas used by leastsq
#  @param y datas used by leastsq
#
#  @return x first datas
#  @return plsq[0] values of the least square algorithm
def fit(x,y):
    est = [1, 1]
    plsq = leastsq(residuals, est, args=(y, x))
    return x,plsq[0]

## residuals
#  Returns the value of the residuals
#
#  @param x datas x
#  @param y datas we want to approximate by a line
#  @param p coefficients of the polynom approximating datas x
#
#  @return value of the residuals
def residuals(p, y, x):
    err = y-peval(x,p)
    return err

## slope_calc
#  Returns the slope of the line approximating data
#
#  @param data the datas we approximate
#  @param func theoretical error
#
#
#  @return the slope of the estimated error
def slope_calc(data,func):
	try:
		data1,data2=changedata(data,func)
		a,b = fit(data1,data2)
		return b
	except ValueError:
		raise

## validation of a set of datas,func
#  Returns true if the approximated slope verifies
#  abs(approximated_slope-real_slope)<precision.
#  Returns false in the other case
#
#  @param data the datas we approximate
#  @param func theoretical error
#  @param precision a real number
#
#  @return a boolean
def validation(data,func,depend,precision):
	try:
		result=slope_calc(data,func)
		pente_pratique=result[1]
		func2=changefunc(func)
		x=1
		a=eval(func2,globals(),{depend:x})
		x=2
		b=eval(func2,globals(),{depend:x})
		pente_theorique=(b-a)/log(2)
		return pente_theorique-pente_pratique<=precision
	except ValueError:
		raise

## validation of multiple sets of datas,func
#  Iterates throw the sets and return True if all the validations
#  are True, False in other cases.
#
#  @param datas the datas we approximate
#  @param funcs theoretical error
#  @param precision a real number
#
#  @return a boolean
def validation_all(datas,funcs,depends,precision):
	i=0;
	result=1
	for data in datas:
		try:
			temp=validation(data,funcs[i],depends[i],precision)
		except ValueError:
			raise
		result=result and temp
		i=i+1
	return result

# Local Variables:
# indent-tabs-mode: t
# End:
