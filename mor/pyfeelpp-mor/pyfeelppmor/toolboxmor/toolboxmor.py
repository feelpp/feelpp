import sys
from pyfeelpptoolboxes.heat import *
from pyfeelppmor.toolboxmor import *
from pyfeelpp import discr

o=toolboxes_options("heat")
o.add(makeToolboxMorOptions())
e=core.Environment(sys.argv,opts=o)

heatBox=heat(dim=2,order=1)
heatBox.init()
model = toolboxmor(2)
vh = discr.functionSpace( heatBox.mesh(), "Pch",1)
model.setFunctionSpaces( Vh=vh)

# def assembleDEIM():
#     heatBox.updateParameterValues()
#     return heatBox.rhs()

# model.setAssembleDEIM(assembleDEIM)
# xh = f.spaceTemperature()
# f.solve()
# f.exportResults()
