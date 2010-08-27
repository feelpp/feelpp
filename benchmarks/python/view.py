# -*- coding: utf-8 -*-

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

from PyQt4 import QtCore, QtGui

class Ui_Main(QtGui.QMainWindow):
    def setupUi(self, Main):
        Main.setObjectName("Main")
        Main.resize(QtCore.QSize(QtCore.QRect(0,0,800,600).size()).expandedTo(Main.minimumSizeHint()))
        Main.setMinimumSize(QtCore.QSize(400,300))
        Main.setAutoFillBackground(False)

        self.centralwidget = QtGui.QWidget(Main)
        self.centralwidget.setObjectName("centralwidget")

        self.gridlayout = QtGui.QGridLayout(self.centralwidget)
        self.gridlayout.setObjectName("gridlayout")

        self.gridlayout1 = QtGui.QGridLayout()
        self.gridlayout1.setObjectName("gridlayout1")

        self.tabs = QtGui.QTabWidget(self.centralwidget)

        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(100)
        sizePolicy.setVerticalStretch(100)
        sizePolicy.setHeightForWidth(self.tabs.sizePolicy().hasHeightForWidth())
        self.tabs.setSizePolicy(sizePolicy)
        self.tabs.setBaseSize(QtCore.QSize(571,521))
        self.tabs.setObjectName("tabs")

        self.tab = QtGui.QWidget()
        self.tab.setObjectName("tab")

        self.graphicsView = QtGui.QLabel(self.tab)
        self.graphicsView.setGeometry(QtCore.QRect(0,0,511,491))
        self.graphicsView.setScaledContents(True)
        self.graphicsView.setAlignment(QtCore.Qt.AlignCenter)
        self.graphicsView.setObjectName("graphicsView")
        self.tabs.addTab(self.tab,"")

        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName("tab_2")

        self.screenshotView = QtGui.QLabel(self.tab_2)
        self.screenshotView.setGeometry(QtCore.QRect(0,0,511,491))
        self.screenshotView.setScaledContents(True)
        self.screenshotView.setAlignment(QtCore.Qt.AlignCenter)
        self.screenshotView.setObjectName("screenshotView")
        self.tabs.addTab(self.tab_2,"")
        self.gridlayout1.addWidget(self.tabs,0,0,1,1)

        self.vboxlayout = QtGui.QVBoxLayout()
        self.vboxlayout.setObjectName("vboxlayout")

        self.progname = QtGui.QLabel(self.centralwidget)

        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(90)
        font.setBold(True)
        self.progname.setFont(font)
        self.progname.setObjectName("progname")
        self.vboxlayout.addWidget(self.progname)

        spacerItem = QtGui.QSpacerItem(20,40,QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
        self.vboxlayout.addItem(spacerItem)

        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.vboxlayout.addWidget(self.label_3)

        self.tableParam = QtGui.QTableWidget(self.centralwidget)
        self.tableParam.setObjectName("tableParam")
        self.tableParam.setColumnCount(2)
        self.tableParam.setRowCount(2)
        self.tableParam.horizontalHeader().resizeSection(0, 110)
        self.tableParam.horizontalHeader().resizeSection(1, 132)
        self.tableParam.verticalHeader().hide()
        self.vboxlayout.addWidget(self.tableParam)

        spacerItem1 = QtGui.QSpacerItem(20,40,QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
        self.vboxlayout.addItem(spacerItem1)

        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.vboxlayout.addWidget(self.label_2)

        self.lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.lineEdit.setReadOnly(True)
        self.lineEdit.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit.setObjectName("lineEdit")
        self.vboxlayout.addWidget(self.lineEdit)

        spacerItem2 = QtGui.QSpacerItem(20,40,QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
        self.vboxlayout.addItem(spacerItem2)

        self.label = QtGui.QLabel(self.centralwidget)

        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.vboxlayout.addWidget(self.label)


        self.comboOutput = QtGui.QComboBox(self.centralwidget)
        self.comboOutput.setObjectName("comboOutput")
        self.vboxlayout.addWidget(self.comboOutput)


        self.hboxlayout = QtGui.QHBoxLayout()
        self.hboxlayout.setObjectName("hboxlayout")

        self.label2 = QtGui.QLabel(self.centralwidget)

        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label2.setFont(font)
        self.label2.setObjectName("label2")
        self.vboxlayout.addWidget(self.label2)

        self.precision = QtGui.QLineEdit(self.centralwidget)
        self.precision.setReadOnly(False)
        self.precision.setAlignment(QtCore.Qt.AlignCenter)
        self.precision.setObjectName("precision")
        self.vboxlayout.addWidget(self.precision)

        self.vboxlayout.addLayout(self.hboxlayout)

        spacerItem3 = QtGui.QSpacerItem(20,40,QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Expanding)
        self.vboxlayout.addItem(spacerItem3)

        self.buttonRun = QtGui.QPushButton(self.centralwidget)
        self.buttonRun.setEnabled(False)
        self.buttonRun.setObjectName("buttonRun")
        self.vboxlayout.addWidget(self.buttonRun)
        self.gridlayout1.addLayout(self.vboxlayout,0,1,1,1)
        self.gridlayout.addLayout(self.gridlayout1,0,0,1,1)
        Main.setCentralWidget(self.centralwidget)

        self.menubar = QtGui.QMenuBar(Main)
        self.menubar.setGeometry(QtCore.QRect(0,0,800,34))
        self.menubar.setObjectName("menubar")

        self.menuProgram = QtGui.QMenu(self.menubar)
        self.menuProgram.setObjectName("menuProgram")

        self.menuMode = QtGui.QMenu(self.menubar)
        self.menuMode.setObjectName("menuMode")

        self.menuAbout = QtGui.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        Main.setMenuBar(self.menubar)

        self.statusbar = QtGui.QStatusBar(Main)
        self.statusbar.setObjectName("statusbar")
        Main.setStatusBar(self.statusbar)

        self.actionOpen = QtGui.QAction(Main)
        self.actionOpen.setObjectName("actionOpen")

        self.actionNormal = QtGui.QAction(Main)
        self.actionNormal.setObjectName("actionNormal")

        self.actionExpert = QtGui.QAction(Main)
        self.actionExpert.setObjectName("actionExpert")

        self.actionChange_XML_dir = QtGui.QAction(Main)
        self.actionChange_XML_dir.setObjectName("actionChange_XML_dir")

        self.actionRun = QtGui.QAction(Main)
        self.actionRun.setObjectName("actionRun")
        self.actionRun.setEnabled(False)

        self.actionImport = QtGui.QAction(Main)
        self.actionImport.setObjectName("actionImport")
        self.actionImport.setEnabled(False)

        self.actionExit = QtGui.QAction(Main)
        self.actionExit.setObjectName("actionExit")

        self.actionValidation = QtGui.QAction(Main)
        self.actionValidation.setObjectName("actionValidation")

        self.actionTests = QtGui.QAction(Main)
        self.actionTests.setObjectName("actionTests")

        self.actionHelp = QtGui.QAction(Main)
        self.actionHelp.setObjectName("actionHelp")

        self.actionVersion = QtGui.QAction(Main)
        self.actionVersion.setObjectName("actionVersion")

        self.menuProgram.addAction(self.actionOpen)
        self.menuProgram.addAction(self.actionImport)
        self.menuProgram.addAction(self.actionRun)
        self.menuProgram.addAction(self.actionExit)
        self.menuMode.addAction(self.actionNormal)
        self.menuMode.addAction(self.actionExpert)
        self.menuAbout.addAction(self.actionHelp)
        self.menuAbout.addAction(self.actionVersion)
        self.menubar.addAction(self.menuProgram.menuAction())
        self.menubar.addAction(self.menuMode.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())

        self.retranslateUi(Main)
        self.tabs.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Main)

    def retranslateUi(self, Main):
        Main.setWindowTitle(QtGui.QApplication.translate("Main", "Validation tool", None, QtGui.QApplication.UnicodeUTF8))
        self.tabs.setTabText(self.tabs.indexOf(self.tab), QtGui.QApplication.translate("Main", "Graph view", None, QtGui.QApplication.UnicodeUTF8))
        self.graphicsView.setText(QtGui.QApplication.translate("Main", "No data", None, QtGui.QApplication.UnicodeUTF8))
        self.screenshotView.setText(QtGui.QApplication.translate("Main", "No data", None, QtGui.QApplication.UnicodeUTF8))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_2), QtGui.QApplication.translate("Main", "Screenshot view", None, QtGui.QApplication.UnicodeUTF8))
        self.progname.setText(QtGui.QApplication.translate("Main", "Program : ", None, QtGui.QApplication.UnicodeUTF8))
        self.tableParam.setRowCount(0)
        self.tableParam.clear()
        self.tableParam.setColumnCount(2)
        self.tableParam.setRowCount(0)

        headerItem = QtGui.QTableWidgetItem()
        headerItem.setText(QtGui.QApplication.translate("Main", "Parameter", None, QtGui.QApplication.UnicodeUTF8))
        self.tableParam.setHorizontalHeaderItem(0,headerItem)

        headerItem1 = QtGui.QTableWidgetItem()
        headerItem1.setText(QtGui.QApplication.translate("Main", "Value", None, QtGui.QApplication.UnicodeUTF8))
        self.tableParam.setHorizontalHeaderItem(1,headerItem1)

        self.label_2.setText(QtGui.QApplication.translate("Main", "Available values for selected param. :", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("Main", "Parameters :", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("Main", "Outputs :", None, QtGui.QApplication.UnicodeUTF8))
        self.label2.setText(QtGui.QApplication.translate("Main", "Precision :", None, QtGui.QApplication.UnicodeUTF8))
        self.buttonRun.setText(QtGui.QApplication.translate("Main", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.menuProgram.setTitle(QtGui.QApplication.translate("Main", "Program", None, QtGui.QApplication.UnicodeUTF8))
        self.menuMode.setTitle(QtGui.QApplication.translate("Main", "Mode", None, QtGui.QApplication.UnicodeUTF8))
        self.menuAbout.setTitle(QtGui.QApplication.translate("Main", "About", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen.setText(QtGui.QApplication.translate("Main", "Open", None, QtGui.QApplication.UnicodeUTF8))
        self.actionImport.setText(QtGui.QApplication.translate("Main", "Import XML", None, QtGui.QApplication.UnicodeUTF8))
        self.actionNormal.setText(QtGui.QApplication.translate("Main", "Mode normal", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExpert.setText(QtGui.QApplication.translate("Main", "Mode expert", None, QtGui.QApplication.UnicodeUTF8))
        self.actionChange_XML_dir.setText(QtGui.QApplication.translate("Main", "Change XML dir", None, QtGui.QApplication.UnicodeUTF8))
        self.actionRun.setText(QtGui.QApplication.translate("Main", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setText(QtGui.QApplication.translate("Main", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionValidation.setText(QtGui.QApplication.translate("Main", "Validation", None, QtGui.QApplication.UnicodeUTF8))
        self.actionTests.setText(QtGui.QApplication.translate("Main", "Tests", None, QtGui.QApplication.UnicodeUTF8))
        self.actionHelp.setText(QtGui.QApplication.translate("Main", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.actionVersion.setText(QtGui.QApplication.translate("Main", "Version", None, QtGui.QApplication.UnicodeUTF8))

# Local Variables:
# indent-tabs-mode: t
# End:
