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

import codeCalcul
import displaySolution

class CodeCalculParaview(codeCalcul.CodeCalcul):
	def __init__(self,pname,pdir):
		codeCalcul.CodeCalcul.__init__(self,pname,pdir)

	def launch(self,params_val,depend_val,output,precision,mode):

		valide,pdfname,plop = codeCalcul.CodeCalcul.launch(self,params_val,depend_val,output,precision,mode)

		images=[]
		for i in range(0,len(depend_val)):
			for j in range(0,len(depend_val[depend_val.keys()[i]])):
				img_dir=self.base_dir+self.vname+"/"
				for k in range(0,len(self.params)):
					if self.params[k].getName() == depend_val.keys()[i]:
						img_dir+=self.params[k].getName()+"_"+depend_val[self.params[k].getName()][j]+"/"
					else:
						img_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

				case=str(img_dir+self.vname+"-1_0.case")
				#print case
				image=str(img_dir+"image.png")
				images.append(image)
				title=str(self.vname)
				legend="legende"
				displaySolution.paraviewScreenshot2D(case,image,title,legend,600,1,-45,0,0)

		return valide,pdfname,images


	def launch_plot(self,params_val,depend_val,output,precision,mode):
		valide,pdfname,plop = codeCalcul.CodeCalcul.launch_plot(self,params_val,depend_val,output,precision,mode)
		images=[]
		for i in range(0,len(depend_val)):
			for j in range(0,len(depend_val[depend_val.keys()[i]])):
				img_dir=self.base_dir+self.vname+"/"
				for k in range(0,len(self.params)):
					if self.params[k].getName() == depend_val.keys()[i]:
						img_dir+=self.params[k].getName()+"_"+depend_val[self.params[k].getName()][j]+"/"
					else:
						img_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

				case=str(img_dir+self.vname+"-1_0.case")
				#print case
				image=str(img_dir+"image.png")
				images.append(image)
				title=str(self.vname)
				legend="legende"
				displaySolution.paraviewScreenshot2D(case,image,title,legend,600,1,-45,0,0)

		return valide,pdfname,images

# Local Variables:
# indent-tabs-mode: t
# End:
