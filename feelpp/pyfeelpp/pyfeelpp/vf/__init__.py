import pyfeelpp.core as core
from _vf import *


def expr(e,filename=None,worldComm=None,directory=None):
    if filename is None:
        filename=""
    if directory is None:
        directory=""
    if worldComm is None:
        worldComm=core.Environment.worldComm()
    return expr_(e,filename=filename,dir=directory,worldComm=worldComm)
