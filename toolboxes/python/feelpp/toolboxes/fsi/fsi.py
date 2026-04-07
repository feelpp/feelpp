import sys
import feelpp.core as fppc
from feelpp.toolboxes.core import toolboxes_options
from feelpp.toolboxes.fsi import fsi

e = fppc.Environment(sys.argv, opts=toolboxes_options("fsi"))



f = fsi(dim=2, orderU=2, orderP=1, orderGeo=1, worldComm=e.worldCommPtr())
f.init()
f.printAndSaveInfo()
f.solve()
f.exportResults()
