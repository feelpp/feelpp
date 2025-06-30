import sys


import feelpp.core as fppc

import feelpp.core.quality as quality
import feelpp.toolboxes.core as tb
import feelpp.core.interpolation as I
from feelpp.toolboxes.fluid import *
from feelpp.toolboxes.cfpdes import *
import feelpp.core.meshmover as mm
import mpi4py
mpi4py.rc.thread_level="single"
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np

## Utils 

def remesh_toolbox(f, hclose, hfar, required_facets, required_elts, parent_mesh, cst):
    
        n_required_elts_before=fppc.nelements(fppc.markedelements(f.mesh(),required_elts))
        n_required_facets_before=fppc.nelements(fppc.markedfaces(f.mesh(),required_facets))
        print(" . [before remesh]   n required elts: {}".format(n_required_elts_before))
        print(" . [before remesh] n required facets: {}".format(n_required_facets_before))

        new_mesh, cpt = fppc.remesh(
            mesh=f.mesh(), metric="gradedls({},{})".format(hclose, hfar), required_elts=required_elts, required_facets=required_facets, params='{"remesh":{ "verbose":-1}}')
        
        print(" . [after remesh]  n remeshes: {}".format(cpt))
        n_required_elts_after=fppc.nelements(fppc.markedelements(new_mesh,required_elts))
        n_required_facets_after=fppc.nelements(fppc.markedfaces(new_mesh,required_facets))
        print(" . [after remesh]  n required elts: {}".format(n_required_elts_after))
        print(" . [after remesh] n required facets: {}".format(n_required_facets_after))
        f.applyRemesh(f.mesh(),new_mesh)

##

## Genral parameters
folder = "shaposwimmer"
sys.argv = [folder]
e = fppc.Environment(
                sys.argv, opts= fppc.backend_options("Iv")
                                .add(tb.toolboxes_options("fluid", "fluid"))
                                .add(tb.toolboxes_options("fluid", "pfluid"))
                                .add(tb.toolboxes_options("fluid", "dfluid"))
                                .add(tb.toolboxes_options("coefficient-form-pdes", "expansion")),
                config=fppc.localRepository(folder))

cst = 0.98
hfar = 0.005
hclose = 0.005
qual = 0.4  #lower bound for the quality of the mesh





mesh1 = fppc.load(fppc.mesh(dim=3,realdim=3), "fluidandswimmer.geo" , 0.05)
mesh2 = fppc.createSubmesh(mesh=mesh1, range=fppc.markedelements(mesh1, "Fluid"))


fppc.Environment.setConfigFile('test2.cfg')
exporter1 = fppc.exporter(mesh=mesh1, name="fluidandswimmer", geo="change")
exporter2 = fppc.exporter(mesh=mesh2, name="fluid", geo="change")


#interpolation between mesh1 and mesh2
Pchv2_mesh1 = fppc.functionSpace(mesh=mesh1, space = "Pchv", order=2)
Pchv2_mesh2 = fppc.functionSpace(mesh=mesh2, space = "Pchv", order=2)

interp_swimmer_to_laplacian = I.interpolator(domain = Pchv2_mesh1, image = Pchv2_mesh2,  range = fppc.elements(mesh2)) 

ud = Pchv2_mesh1.element()
ud.on(range=fppc.elements(mesh1), expr=fppc.expr("{1,1,1}", row=3))

print("Ud interpolation")
ud_interp = interp_swimmer_to_laplacian.interpolate(ud)
print("Ud save")

w = Pchv2_mesh2.element()
w.on(range=fppc.elements(mesh2), expr=fppc.expr("{1,1,1}", row=3))

v = w - ud_interp
l2_v = fppc.normL2(range=fppc.elements(mesh2), expr=v)

print(l2_v)

