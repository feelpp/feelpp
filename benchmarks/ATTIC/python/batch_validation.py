#!@PYTHON_EXECUTABLE@

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

from optparse import OptionParser
import codeCalcul
import util
import sys
import string


# Creation des options
parser = OptionParser()
parser.add_option("--path", dest="path", help="path of the program to execute")
parser.add_option("--program", dest="program", help="program to execute")
parser.add_option("--precision", type="float",dest="precision", help="precision of the validation", default=0.25)
parser.add_option("--mode", type="int",dest="mode", help="mode of the validation(1=range verification, 0=no range verification)", default=0)
parser.add_option("--cmdargs", dest="cmdargs", help="command line arguments for the c++ code execution", default="")
parser.add_option("--launchnbr", dest="lnbr", type="int", help="number of times the c++ code is launched", default=3)
parser.add_option("--method", dest="method", type="int", help="method for generating intervals : 0 - linear, 1 - logarithmic", default=1)
(options, args) = parser.parse_args()

print "program:",options.program
print "path:",options.path
print "mode:",options.mode


# Verification des options
if (not options.path or not options.program):
	print "[validation] error : path and/or program name is not defined"
	sys.exit(-1)

if (options.precision < 0):
	print "[validation] error : precision must be positive"
	sys.exit(-1)

# Initialisation du code de calcul
code=0
try:
	code = codeCalcul.CodeCalcul(options.program,options.path)
except IOError:
	print "[validation] error : CodeCalcul initialisation failed"
	sys.exit(-1)

# Ajout des eventuelles options
code.addOptions(string.replace(string.replace(options.cmdargs,"\ ", " "), "\"", "" ) )

# Initialisation des variables
val_param={}
val_depend={}

# Creations des valeurs nominales
for p in code.params:
	if p.getAttrValues()[1]=="discrete":
		val_param[p.getName()]=[util.strToSci(str(p.getValues().split(',')[len(p.getValues().split(','))/2]))]
	else:
			if len(p.getValues().split(':'))<3:
					val_param[p.getName()]=[util.nbrToSci((float(p.getValues().split(':')[0])+float(p.getValues().split(':')[-1]))/2)]
			else:
					val_param[p.getName()]=[util.nbrToSci(float(p.getValues().split(':')[1]))]

print "[validation] info : nominal values : "
print val_param

# Creation des valeurs variables
for o in code.outputs:
	for d in o.getDependencies():
		if not val_depend.has_key(str(d)):
			for p in code.params:
				if p.getName()==str(d):
					if p.getAttrValues()[1]=="discrete":
						val_depend[str(d)]=map(util.strToSci,p.getValues().split(','))
					else:
						val_depend[str(d)]=map(util.strToSci,util.pickInInterval(float(p.getValues().split(':')[0]),float(p.getValues().split(':')[-1]),options.lnbr,options.lnbr+1,options.method))

# les outputs n'ayant pas de dependances sont enleves
dep_outputs=[]
for o in code.outputs:
	if len(o.getDependencies())>0:
		dep_outputs.append(o)

print "[validation] info : variable values : "
print val_depend

# Lancement de la validation
try:
		valide,graphfile,imgs=code.launch(val_param,val_depend,dep_outputs,options.precision,options.mode)
except IOError:
		print "[validation] error : the test case stopped unexpectedly..."
		valide=0

print "[validation] results : "
if valide:
	print "["+options.program+"] : Ok"
else:
	print "["+options.program+"] : Failed"

sys.exit(int(not valide))

# Local Variables:
# indent-tabs-mode: t
# End:
