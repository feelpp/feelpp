#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
###theStudy = salome.myStudy

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

Bitter = geompy.ImportSTEP("BE-02-01M.stp", True, True)

Channel0 = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(Channel0, [412, 409])

Channel1 = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(Channel1, [419, 266])

TieRods = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(TieRods, [400, 404])

top = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(top, [13])

bottom = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(bottom, [140])

V0 = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(V0, [417, 415])

V1 = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
geompy.UnionIDs(V1, [3, 406])

Bord = geompy.CreateGroup(Bitter, geompy.ShapeType["FACE"])
Faces = geompy.ExtractShapes(Bitter, geompy.ShapeType["FACE"], True)
print("Faces: %d" % len(Faces))
geompy.UnionList(Bord, Faces)

CoolingHoles = geompy.CutListOfGroups([Bord], [Channel0, Channel1, TieRods, top, bottom, V0, V1])

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

geompy.addToStudy( Bitter, 'Bitter' )
geompy.addToStudyInFather( Bitter, Channel0, 'Channel0' )
geompy.addToStudyInFather( Bitter, Channel1, 'Channel1' )
geompy.addToStudyInFather( Bitter, TieRods, 'TieRods' )
geompy.addToStudyInFather( Bitter, top, 'top' )
geompy.addToStudyInFather( Bitter, bottom, 'bottom' )
geompy.addToStudyInFather( Bitter, V0, 'V0' )
geompy.addToStudyInFather( Bitter, V1, 'V1' )
#geompy.addToStudyInFather( Bitter, Bord, 'Bord' )
geompy.addToStudyInFather( Bitter, CoolingHoles, 'CoolingHoles' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Bitter)
MG_CADSurf = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf)
MG_CADSurf_Parameters_1 = MG_CADSurf.Parameters()
MG_CADSurf_Parameters_1.SetPhySize( 11.9696 )
MG_CADSurf_Parameters_1.SetMinSize( 0.119696 )
MG_CADSurf_Parameters_1.SetMaxSize( 23.9391 )
MG_CADSurf_Parameters_1.SetChordalError( 5.98478 )

MG_Tetra_1 = Mesh_1.Tetrahedron(algo=smeshBuilder.MG_Tetra)
MG_Tetra_Parameters_1 = MG_Tetra_1.Parameters()

Bitter_1 = Mesh_1.GroupOnGeom(Bitter,'Bitter',SMESH.VOLUME)
Channel0_1 = Mesh_1.GroupOnGeom(Channel0,'Channel0',SMESH.FACE)
Channel1_1 = Mesh_1.GroupOnGeom(Channel1,'Channel1',SMESH.FACE)
TieRods_1 = Mesh_1.GroupOnGeom(TieRods,'TieRods',SMESH.FACE)
top_1 = Mesh_1.GroupOnGeom(top,'top',SMESH.FACE)
bottom_1 = Mesh_1.GroupOnGeom(bottom,'bottom',SMESH.FACE)
V0_1 = Mesh_1.GroupOnGeom(V0,'V0',SMESH.FACE)
V1_1 = Mesh_1.GroupOnGeom(V1,'V1',SMESH.FACE)
CoolingHoles_1 = Mesh_1.GroupOnGeom(CoolingHoles,'CoolingHoles',SMESH.FACE)

isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(MG_CADSurf.GetAlgorithm(), 'MG-CADSurf')
smesh.SetName(MG_Tetra_Parameters_1, 'MG-Tetra Parameters_1')
smesh.SetName(MG_CADSurf_Parameters_1, 'MG-CADSurf Parameters_1')

smesh.SetName(Channel0_1, 'Channel0')
smesh.SetName(Channel1_1, 'Channel1')
smesh.SetName(TieRods_1, 'TieRods')
smesh.SetName(top_1, 'top')
smesh.SetName(bottom_1, 'bottom')
smesh.SetName(V0_1, 'V0')
smesh.SetName(V1_1, 'V1')
smesh.SetName(Bitter_1, 'Bitter')

try:
  Mesh_1.ExportMED('Bitter.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')



if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
