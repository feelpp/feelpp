#! /usr/bin/env python
# author(s): Guillaume Dolle
#
# This example shows some basic operations to manipulate Feel++ journals (JSON
# files) generated automatically. Several partial journal files have been gathered
# in the ./json folder. Note that the format might change and is application based!
#
# This example is for demonstration purpose only! You should add better checks
# for proper comparison!

import os
import json
from decimal import Decimal
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt

j=[]
for filename in os.listdir("./json"):
    with open( os.getcwd()+ "/json/" + filename ) as fp:
        j.append( json.load(fp) )

#pprint(j[0])
#pprint( j[0]["database"])

# Display several informations
for i in j:
    print( "start time: " + i["journal"]["time"]["gm"] )
    print( "number of processors: " + i["environment"]["mpi"]["number_processors"] )
    print( "average mesh size: " + i["mesh"]["mesh-3"]["h_average"] )
    print( "---")


mesh_size = []
time_a_solve = []
for i in j:
    mesh_size.append( Decimal(i["mesh"]["mesh-3"]["h_average"]) )
    for timer in i["timers"]:
        if( i["timers"][timer]["message"] == "a.solve" ):
            time_a_solve.append( Decimal( i["timers"][timer]["total"] ) )

plt.bar(mesh_size, time_a_solve, width=0.008 )
plt.show()
