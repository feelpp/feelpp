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

## \package displaySolution
#  This package creates a bridge between Paraview and the graphical interface

from paraview import servermanager as sm


##	paraviewScreenshot2D
#	Creates an png image of the solution using the paraview module in python.
#
#	@param casefilename filename of the EnSight .case file.
#	@param imagefilename filename of the image to be saved.
#	@param title title of the image
#	@param legendlabel legend string
#	@param size size of the image (final image is size*size)
#	@param domain the type of representation for the domain
#   <ul>
#	   <li>0 : Points</li>
#	   <li>1 : Wireframes</li>
#	   <li>2 : Surfaces</li>
#	   <li>3 : Outline</li>
#   </ul>
#	@param elevation vertical angle of the camera
#	@param azimuth horizontal angle of the camera
#	@param roll rotation angle of the camera
def paraviewScreenshot2D(casefilename,imagefilename,title,legendlabel,size,domain,elevation,azimuth,roll):
	print "[displayGraphics] : paraviewScreenshot2D...\n"
	# Connection aux donnees et creation d'une entite de visualisation
	view,reader=connection(casefilename)
	# Recuperation de la solution
	# Sous forme de surfaces colorees representant la solution u
	func = sm.filters.ElevationFilter(Input=reader)
	# Et en relief par la meme occasion
	func = sm.filters.WarpScalar(Input=reader)
	# Creation de la representation contenant la solution
	rep = sm.CreateRepresentation(func,view)
	# Recuperation du range des valeurs a afficher
	ai = func.GetDataInformation().GetPointDataInformation().GetArrayInformation(0)
	# Un truc obscure qui permet de representer la solution
	lt = sm.rendering.PVLookupTable()
	rep.LookupTable = lt
	# Creation des couleurs pour le ElevationFilter
	rep.ColorAttributeType = 0
	rep.ColorArrayName = ai.GetName()
	lt.RGBPoints = [ai.GetComponentRange(0)[0], 0, 0, 1, ai.GetComponentRange(0)[1], 1, 0, 0]
	func.SetPropertyWithName("ScaleFactor",1/(ai.GetComponentRange(0)[1]-ai.GetComponentRange(0)[0]))
	lt.ColorSpace = 1

	# Affichage des legendes
	printLegend(title, view, 20, [0.4,0.9])
	printLegend(legendlabel, view, 14, [0.05,0.05])

	# Affichage de l'echelle a droite de l'image
	printScalarBar(view,ai,lt)

	# Mode de representation
	rep.Representation=domain

	# Changement de point de vue de la camera pour voir l'integralite du maillage
	view.ViewSize=[size,size]
	camera = view.GetActiveCamera()
	view.ResetCamera()
	pos=view.CameraPosition.GetData()
	view.CameraPosition=[pos[0],pos[1],pos[2]*1.25]
	camera.Elevation(elevation)
	camera.Roll(roll)
	camera.Azimuth(azimuth)

	view.StillRender()
	print "[displayGraphics] : saving image..."
	view.WriteImage(imagefilename, "vtkPNGWriter", 1)

def paraviewScreenshot3D(casefilename,imagefilename,title,legendlabel,size,domain,elevation,azimuth,roll):
	print "[displayGraphics] : paraviewScreenshot3D...\n"
	# Connection aux donnees et creation d'une entite de visualisation
	view,reader=connection(casefilename)
	# Recuperation de la solution
	# Sous forme de surfaces colorees representant la solution u
	func = sm.filters.ElevationFilter(Input=reader)
	# Et en relief par la meme occasion
	func = sm.filters.Contour(Input=reader)
	# Creation de la representation contenant la solution
	rep = sm.CreateRepresentation(func,view)
	# Recuperation du range des valeurs a afficher
	ai = func.GetDataInformation().GetPointDataInformation().GetArrayInformation(0)
	# Un truc obscure qui permet de representer la solution
	lt = sm.rendering.PVLookupTable()
	rep.LookupTable = lt

	func.ComputeScalars=1
	func.ComputeNormals=0
	func.ContourValues=[ai.GetComponentRange(0)[0],ai.GetComponentRange(0)[0]/2,0,ai.GetComponentRange(0)[1]/2,ai.GetComponentRange(0)[1]]
	func.SelectInputScalars=ai.GetName()

	print ai
	print ai.GetComponentRange(0)[0]
	print ai.GetComponentRange(0)[1]

	# Creation des couleurs pour le ElevationFilter
	rep.ColorAttributeType = 0
	rep.ColorArrayName = ai.GetName()
	lt.RGBPoints = [ai.GetComponentRange(0)[0], 0, 0, 1, ai.GetComponentRange(0)[1], 1, 0, 0]

	lt.ColorSpace = 1

	# Affichage des legendes
	printLegend(title, view, 20, [0.4,0.9])
	printLegend(legendlabel, view, 14, [0.05,0.05])

	# Affichage de l'echelle a droite de l'image
	printScalarBar(view,ai,lt)

	# Mode de representation
	rep.Representation=domain

	# Changement de point de vue de la camera pour voir l'integralite du maillage
	view.ViewSize=[size,size]
	camera = view.GetActiveCamera()
	view.ResetCamera()
	pos=view.CameraPosition.GetData()
	view.CameraPosition=[pos[0],pos[1],pos[2]*1.25]
	camera.Elevation(elevation)
	camera.Roll(roll)
	camera.Azimuth(azimuth)

	view.StillRender()
	print "[displayGraphics] : saving image..."
	view.WriteImage(imagefilename, "vtkPNGWriter", 1)

def connection(filename):
	print "[displayGraphics] : setting paraview connection..."
	# Connection au servermanager
	if not sm.ActiveConnection:
		connection = sm.Connect()
	# Lecture du fichier .case
	print "[displayGraphics] : reading file..."
	reader=sm.sources.ensight(CaseFileName = filename,ByteOrder = 1)
	# Ajout du contenu du fichier au servermanager
	reader.UpdatePipelineInformation()
	# Creation de la fenetre de rendu
	return sm.CreateRenderView(),reader


def printLegend(string, view, size, pos):
	legend = sm.sources.TextSource()
	legend.Text = string
	legendrep = sm.rendering.TextSourceRepresentation()
	legendrep.TextScaleMode=0
	legendrep.FontSize=size
	legendrep.Bold=1
	legendrep.FontFamily=1
	legendrep.Input = legend
	legendrep.Position=pos
	legendrep.UpdateVTKObjects()
	view.Representations.append(legendrep)

def printScalarBar(view,ai,lt):
	bar = sm.rendering.ScalarBarWidgetRepresentation(registrationGroup='scalar_bars', registrationName="ScalarBarWidgetRepresentation1")
	bar.LabelFontFamily=1
	bar.TitleFontFamily=1
	bar.Title =ai.GetName()
	bar.LookupTable = lt
	view.Representations.append(bar)


# Local Variables:
# indent-tabs-mode: t
# End:
