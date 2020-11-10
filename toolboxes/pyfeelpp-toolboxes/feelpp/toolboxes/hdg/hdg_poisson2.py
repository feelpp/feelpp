import sys
from pyfeelpptoolboxes.hdg import *

e=core.Environment(sys.argv,opts=hdg_poisson_options())


p1=hdgpoisson(dim=2,order=1)
p2=hdgpoisson(dim=2,order=1)

print("solving p1...")
p1.init()
p1.printAndSaveInfo()
p1.solve()
p1.exportResults()

print("solving p2...")
p2.init()
p2.printAndSaveInfo()
p2.solve()
p2.exportResults()


