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

from xmlParser import *
import validation
import util

class CodeCalcul:
	def __init__(self,pname,pdir):
		self.name=pname
		self.exec_path=pdir
		self.xml_path=pdir
		self.other_options=""
		self.base_dir=os.environ["HOME"]+"/feel/"
		self.info()

	def info(self):
		print self.exec_path+self.name+" --capabilities"
		if not os.system(self.exec_path+self.name+" --capabilities")==0:
			raise IOError, "[CodeCalcul] error : file not found"
		self.vname, self.params, self.outputs=parse_xml_response(self.base_dir+"xml/xml_response.xml")

	def addOptions(self,opt):
		self.other_options=opt

	def verify(self,params_val,depend_val,output):
		for out in output:
			if not (out in self.outputs):
				return False
		for pv in params_val.keys():
			present=False
			for p in self.params:
				if p.getName()==pv:
					present=True
					if p.getAttrValues()[1]=="discrete":
						if not (str(params_val[pv][0]) in map(util.strToSci,p.getValues().split(','))):
							return False
					else:
						if float(str(params_val[pv][0]))<float(p.getValues().split(':')[0]) or float(str(params_val[pv][0]))>float(p.getValues().split(':')[-1]):
							return False
			if not present:
				return False
		for pv in depend_val.keys():
			present=False
			for p in self.params:
				if p.getName()==pv:
					present=True
					if p.getAttrValues()[1]=="discrete":
						for val in depend_val[pv]:
							if not (str(val) in p.getValues().split(',')):
								return False
					else:
						for val in depend_val[pv]:
							if float(str(val))<float(p.getValues().split(':')[0]) or float(str(val))>float(p.getValues().split(':')[-1]):
								return False
			if not present:
				return False
		return True

	def launch(self,params_val,depend_val,output,precision,mode):
		if (mode==1):
			if not self.verify(params_val,depend_val,output):
				print "[CodeCalcul] error : Some values are out of range"
				return
		for i in range(0,len(depend_val)):
			for j in range(0,len(depend_val[depend_val.keys()[i]])):
				lString=" "+self.other_options
				img_dir=self.base_dir+self.vname+"/"
				for k in range(0,len(self.params)):
					if self.params[k].getName() == depend_val.keys()[i]:
						img_dir+=self.params[k].getName()+"_"+depend_val[self.params[k].getName()][j]+"/"
						lString+=" --"+self.params[k].getCmdName()+"="+depend_val[self.params[k].getName()][j]
					else:
						img_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"
						lString+=" --"+self.params[k].getCmdName()+"="+params_val[self.params[k].getName()][0]
				print str(self.exec_path+self.name+lString)
				if not os.system(str(self.exec_path+self.name+lString))==0:
					raise IOError, "error in execution"

			graph_dir=self.base_dir+"graph/"+self.name+"/"
			for k in range(0,len(self.params)):
				if not (self.params[k].getName() == depend_val.keys()[i]):
					graph_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

			try:
				os.makedirs(str(graph_dir))
			except OSError:
				pass

			pdf=str(graph_dir)
			pdfname=[]
			for tmpOut in output:
				pdfname.append(pdf+tmpOut.getName()+".pdf")

			valide=validation.validation_code(str(self.vname),
											self.params,
											output,
											str(self.base_dir+"xml/xml_result.xml"),
											params_val,
											depend_val,
											precision,
											pdf)
			index=0
			for tmpPdf in pdfname:
				os.system("convert -density 150 "+tmpPdf+" "+tmpPdf+".png")
				pdfname[index]=pdfname[index]+".png"
				index+=1

		return valide,pdfname,[]


	def launch_plot(self,params_val,depend_val,output,precision,mode):
		if (mode==1):
			if not self.verify(params_val,depend_val,output):
				print "[CodeCalcul] error : Some values are wrong"
				return

		for i in range(0,len(depend_val)):
			for j in range(0,len(depend_val[depend_val.keys()[i]])):
				lString=" "+self.other_options
				img_dir=self.base_dir+self.vname+"/"
				for k in range(0,len(self.params)):
					if self.params[k].getName() == depend_val.keys()[i]:
						img_dir+=self.params[k].getName()+"_"+depend_val[self.params[k].getName()][j]+"/"
						lString+=" --"+self.params[k].getCmdName()+"="+depend_val[self.params[k].getName()][j]
					else:
						img_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"
						lString+=" --"+self.params[k].getCmdName()+"="+params_val[self.params[k].getName()][0]
				print str(self.exec_path+self.name+lString)
				if not os.system(str(self.exec_path+self.name+lString))==0:
					raise IOError, "error in execution"

			graph_dir=self.base_dir+"graph/"+self.name+"/"
			for k in range(0,len(self.params)):
				if not (self.params[k].getName() == depend_val.keys()[i]):
					graph_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

			try:
				os.makedirs(str(graph_dir))
			except OSError:
				pass

			pdf=str(graph_dir)
			pdfname=[]
			for tmpOut in output:
				pdfname.append(pdf+tmpOut.getName()+".pdf")

			for i in range(0,len(output)):
				valide=validation.data_plot(str(self.vname),
												self.params,
												output[i],
												str(self.base_dir+"xml/xml_result.xml"),
												params_val,
												depend_val.keys()[0],
												depend_val,
												precision,
												pdf)
			index=0
			for tmpPdf in pdfname:
				os.system("convert -density 150 "+tmpPdf+" "+tmpPdf+".png")
				pdfname[index]=pdfname[index]+".png"
				index+=1

		return valide,pdfname,[]


	def launchXml(self,params_val,depend_val,output,precision,mode,xmlname):
		for i in range(0,len(depend_val)):
			graph_dir=self.base_dir+"graph/"+self.name+"/"
			for k in range(0,len(self.params)):
				if not (self.params[k].getName() == depend_val.keys()[i]):
					graph_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

			try:
				os.makedirs(str(graph_dir))
			except OSError:
				pass

			pdf=str(graph_dir)
			pdfname=[]
			for tmpOut in output:
				pdfname.append(pdf+tmpOut.getName()+".pdf")

			valide=validation.validation_code(str(self.vname),
											self.params,
											output,
											str(xmlname),
											params_val,
											depend_val,
											precision,
											pdf)
			index=0
			for tmpPdf in pdfname:
				os.system("convert -density 150 "+tmpPdf+" "+tmpPdf+".png")
				pdfname[index]=pdfname[index]+".png"
				index+=1

		return valide,pdfname


	def launch_plotXml(self,params_val,depend_val,output,precision,mode,xmlname):
		for i in range(0,len(depend_val)):
			graph_dir=self.base_dir+"graph/"+self.name+"/"
			for k in range(0,len(self.params)):
				if not (self.params[k].getName() == depend_val.keys()[i]):
					graph_dir+=self.params[k].getName()+"_"+params_val[self.params[k].getName()][0]+"/"

			try:
				os.makedirs(str(graph_dir))
			except OSError:
				pass

			pdf=str(graph_dir)
			pdfname=[]
			for tmpOut in output:
				pdfname.append(pdf+tmpOut.getName()+".pdf")

			for i in range(0,len(output)):
				valide=validation.data_plot(str(self.vname),
												self.params,
												output[i],
												str(xmlname),
												params_val,
												depend_val.keys()[0],
												depend_val,
												precision,
												pdf)
			index=0
			for tmpPdf in pdfname:
				os.system("convert -density 150 "+tmpPdf+" "+tmpPdf+".png")
				pdfname[index]=pdfname[index]+".png"
				index+=1

		return valide,pdfname

# Local Variables:
# indent-tabs-mode: t
# End:

