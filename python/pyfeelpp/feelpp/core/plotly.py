from ._plot import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

colors = ["#1f77b4","#bcbd22", "#8c564b", "#2ca02c","#9467bd","#ff7f0e", "#d62728", "#e377c2", "#17becf"]

_plotly={
        'plotly(2)':ExporterPlotly2,
        'plotly(3)':ExporterPlotly3,
    }

def plotly(mesh):
    key='plotly('+str(mesh.dimension())+')'
    e=_plotly[key](mesh)
    e.plotMesh = lambda **kwargs : plotMesh(e,**kwargs)
    e.plotMesh2d = lambda z=None: plotMesh2d(e,z=z)
    e.plotMesh3d = lambda z=None, intensity=None : plotMesh3d(e,z=z,intensity=intensity )
    e.plotFill = lambda z=None : plotFill(e,z=z)
    e.xRange = lambda : xRange(e)
    e.yRange = lambda : yRange(e)
    return e

def xRange(e):
    return np.amin(e.getPoints()[:,0]), np.amax(e.getPoints()[:,0])

def yRange(e):
    return np.amin(e.getPoints()[:,1]), np.amax(e.getPoints()[:,1])

def plotMesh(e, color=colors[0],legendgroup=0, idx=0):
    mesh_traces = []
    p = e.getPoints()
    
    for indices in e.getConnectivity():
        mesh_traces.append(
            go.Scatter(
                x=p[indices,0],
                y=p[indices,1],
                mode="lines",
                line=dict(color=color, width=1),
                legendgroup="Ω"+str(legendgroup),
                name = "Ω"+str((legendgroup+1 if legendgroup < 4 else 0)),
                legendrank=(legendgroup+1 if legendgroup < 4 else 0),
                showlegend=legendgroup == idx,
            )
        )
        idx+=1
    return mesh_traces


def plotMesh3d(e, z=None, intensity=None, color=colors[0],legendgroup=0, idx=0):
    p = e.getPoints()
    c = e.getConnectivity()
    if z is None:
        z = p[:,0]+p[:,1]
    if intensity is None:
        intensity = z
    return go.Mesh3d(
             x=p[:,0],
             y=p[:,1],
             z=z,
             i = c[:,0],
             j = c[:,1],
             k = c[:,2],
             intensity=intensity,
             showscale=True
         )
def plotFill(e, z=None, color=colors[0],legendgroup=0, idx=0):
    p = e.getPoints()
    c = e.getConnectivity()
    if z is None:
        z = p[:,0]+p[:,1]
    camera = dict(
        eye=dict(x=0, y=0, z=1),
        up=dict(x=0, y=1, z=0),
        projection=dict(type="orthographic"),
    )
    return go.Mesh3d(
             x=p[:,0],
             y=p[:,1],
             z=0*p[:,0],
             i = c[:,0],
             j = c[:,1],
             k = c[:,2],
             intensity=z,
             showscale=True
         ), camera
    
