#!/usr/bin/python

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

import sys
import string
import os
from PyQt4 import QtGui,QtCore
from view import *
import util
try:
		from codeCalculParaview import CodeCalculParaview as CodeCalcul
except ImportError:
		print "[Interface] warning : Python module for paraview not found"
		print "[Interface] warning : switch to default..."
		from codeCalcul import *

class dlgCmd(QtGui.QDialog):
		def __init__(self,parent = None,name = None,modal = 0,fl = 0):
				QtGui.QDialog.__init__(self,parent,name,modal,fl)



class app:
	def __init__(self,args):
		# application Qt
		self.qtapp=QtGui.QApplication(args)
		# creation de la fenetre principale
		self.win=Ui_Main()
		self.win.setupUi(self.win)
		# affichage de la fenetre
		self.win.show()

		# connection SLOT/SIGNAL
		self.qtapp.connect(self.win.actionOpen,QtCore.SIGNAL("activated()"),self.openprogram)
		self.qtapp.connect(self.win.actionImport,QtCore.SIGNAL("activated()"),self.importxml)
		self.qtapp.connect(self.win.buttonRun,QtCore.SIGNAL("clicked()"),self.runProg)
		self.qtapp.connect(self.win.actionRun,QtCore.SIGNAL("activated()"),self.runProg)
		self.qtapp.connect(self.win.actionNormal,QtCore.SIGNAL("activated()"),self.changeModeNormal)
		self.qtapp.connect(self.win.actionExpert,QtCore.SIGNAL("activated()"),self.changeModeExpert)
		self.qtapp.connect(self.win.actionExit, QtCore.SIGNAL("activated()"),self.qtapp,QtCore.SLOT("quit()"))
		self.qtapp.connect(self.qtapp, QtCore.SIGNAL("lastWindowClosed()"),self.qtapp,QtCore.SLOT("quit()"))
		self.qtapp.connect(self.win.actionVersion,QtCore.SIGNAL("activated()"),self.displayAbout)
		self.qtapp.connect(self.win.actionHelp,QtCore.SIGNAL("activated()"),self.displayHelp)
		self.qtapp.connect(self.win.tableParam,QtCore.SIGNAL("itemSelectionChanged()"),self.displayPossibleVals)
		self.qtapp.connect(self.win.tableParam,QtCore.SIGNAL("cellChanged(int,int)"),self.resetLaunch)

		self.mode=1 # normal
		self.launched=False

		self.qtapp.exec_()

	# definition des differents SLOTS
	def openprogram(self):
		tmp = QtGui.QFileDialog.getOpenFileName(self.win, "Open Image", "/home", "Programs (*)");
		pname=tmp.split('/')[len(tmp.split('/'))-1]
		pdir=tmp[0:len(tmp)-len(pname)]
		self.cc=CodeCalcul(str(pname),str(pdir))
		self.win.progname.setText(QtGui.QApplication.translate("Main", "Program : "+self.cc.name, None, QtGui.QApplication.UnicodeUTF8))
		self.win.buttonRun.setEnabled(True)
		self.win.actionRun.setEnabled(True)
		self.win.actionImport.setEnabled(True)

		self.win.tableParam.clear()
		self.win.tableParam.setColumnCount(2)
		self.win.tableParam.setRowCount(len(self.cc.params))

		headerItem = QtGui.QTableWidgetItem()
		headerItem.setText(QtGui.QApplication.translate("Main", "Parameter", None, QtGui.QApplication.UnicodeUTF8))
		self.win.tableParam.setHorizontalHeaderItem(0,headerItem)

		headerItem1 = QtGui.QTableWidgetItem()
		headerItem1.setText(QtGui.QApplication.translate("Main", "Value", None, QtGui.QApplication.UnicodeUTF8))
		self.win.tableParam.setHorizontalHeaderItem(1,headerItem1)

		self.win.precision.setText(QtGui.QApplication.translate("Main", "0.2", None, QtGui.QApplication.UnicodeUTF8))

		for i in range(0,len(self.cc.params)):
			newItem = QtGui.QTableWidgetItem()
			newItem.setText(self.cc.params[i].getName())
			newItem.setFlags(QtCore.Qt.ItemIsSelectable)
			self.win.tableParam.setItem(i,0,newItem)

		self.win.comboOutput.clear()
		for out in self.cc.outputs:
			self.win.comboOutput.addItem(out.getName())

		self.launched=False

	def displayPossibleVals(self):
		#print self.win.tableParam.currentRow()
		self.win.lineEdit.setText(self.cc.params[self.win.tableParam.currentRow()].getValues())

	def resetLaunch(self,i,j):
		self.launched=False

	def importxml(self):
		tmp = QtGui.QFileDialog.getOpenFileName(self.win, "Open XML file", "/home", "XML files (*.xml)");
		params_val={}
		depend_val={}
		for i in range(0,len(self.cc.params)):
			if len(self.win.tableParam.item(i,1).text().split(","))>1:
				depend_val[self.cc.params[i].getName()]=map(util.strToSci,self.win.tableParam.item(i,1).text().split(","))
			else:
				params_val[self.cc.params[i].getName()]=map(util.strToSci,self.win.tableParam.item(i,1).text().split(","))

		if len(depend_val)==0:
			QtGui.QMessageBox.critical(self.win, "Error","There is no variable parameter", QtGui.QMessageBox.Ok)
			return

		#output=[self.cc.outputs[self.win.comboOutput.currentIndex()]]
		output=self.cc.outputs

		if len(self.win.precision.text())>0:
			precision=abs(float(self.win.precision.text()))
		else:
			precision=0.2
		if depend_val.keys()[0] in output[0].getDependencies():
			valide,self.graphfile=self.cc.launchXml(params_val,depend_val,output,precision,self.mode,tmp)
		else:
			valide,self.graphfile=self.cc.launch_plotXml(params_val,depend_val,output,precision,self.mode,tmp)

		graph = QtGui.QPixmap(self.graphfile[self.win.comboOutput.currentIndex()])
		self.win.graphicsView.setFixedSize(graph.width(),graph.height())
		self.win.graphicsView.setPixmap(graph)
		self.launched=True

	def runProg(self):

		if self.launched==False:
			params_val={}
			depend_val={}
			for i in range(0,len(self.cc.params)):
				if len(self.win.tableParam.item(i,1).text().split(","))>1:
					depend_val[self.cc.params[i].getName()]=map(util.strToSci,self.win.tableParam.item(i,1).text().split(","))
				else:
					params_val[self.cc.params[i].getName()]=map(util.strToSci,self.win.tableParam.item(i,1).text().split(","))

			if len(depend_val)==0:
				QtGui.QMessageBox.critical(self.win, "Error","There is no variable parameter", QtGui.QMessageBox.Ok)
				return

			#output=[self.cc.outputs[self.win.comboOutput.currentIndex()]]
			output=self.cc.outputs

			if len(self.win.precision.text())>0:
				precision=abs(float(self.win.precision.text()))
			else:
				precision=0.2
			if depend_val.keys()[0] in output[0].getDependencies():
				valide,self.graphfile,self.imgfiles=self.cc.launch(params_val,depend_val,output,precision,self.mode)
			else:
				valide,self.graphfile,self.imgfiles=self.cc.launch_plot(params_val,depend_val,output,precision,self.mode)
		if (len(self.imgfiles)>0):
				image = QtGui.QPixmap(self.imgfiles[0])
				self.win.screenshotView.setFixedSize(image.width(),image.height())
				self.win.screenshotView.setPixmap(image)
		if (len(self.graphfile)>0):
				graph = QtGui.QPixmap(self.graphfile[self.win.comboOutput.currentIndex()])
				self.win.graphicsView.setFixedSize(graph.width(),graph.height())
				self.win.graphicsView.setPixmap(graph)
		self.launched=True

	def changeModeNormal(self):
		self.mode=1

	def changeModeExpert(self):
		QtGui.QMessageBox.information(self.win, "Expert mode","In expert mode, all values are available.\nBe careful...", QtGui.QMessageBox.Ok)
		self.mode=0

	def displayAbout(self):
		ret = QtGui.QMessageBox.information(self.win, "About","Validation Tool Version 0.1\nThe Feel Test Team", QtGui.QMessageBox.Ok)

	def displayHelp(self):
		ret = QtGui.QMessageBox.information(self.win, "About","For more information about how to use the interface, see the 2009 report of Feel-Test", QtGui.QMessageBox.Ok)

# lancement de l'application
mapp = app(sys.argv)

# Local Variables:
# indent-tabs-mode: t
# End:
