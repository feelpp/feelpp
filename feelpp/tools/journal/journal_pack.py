#! /usr/bin/env python
# author(s): Guillaume Dolle
#
# This example shows some basic operations to manipulate Feel++ journals (JSON
# files) generated automatically. Several partial journal files have been gathered
# in the ./json folder. Note that the format might change and is application based!
#
# This example is for demonstration purpose only! You should add better checks
# for proper comparison!
#
# How to pack/merge several journal files into a static json database. 
# Note that this is not equivalent with the mongodb.

import os
import argparse
import json
from decimal import Decimal
import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt

# Command line arguments
parser = argparse.ArgumentParser(description='Pack journal .json files into a unique one (database)')
parser.add_argument('-o','--output' , metavar="outfile", default="pack.json", type=str, help='output file (.json)' )
parser.add_argument('-i','--inputs' ,metavar="infiles", help='input journal files (.json)')
args = parser.parse_args()

outfile = args.output
infiles = args.inputs

packname = 'pack.json'

j=[]
if(not infiles ):
    for filename in os.listdir("./json"):
        with open( os.getcwd()+ "/json/" + filename ) as fd:
            j.append( json.load(fd) )

jdb={}
for i in j:
    # We use the same uuid for the database.
    db_uuid = i["environment"]["run"]["uuid"]
    jdb[db_uuid]=j[0]
    print(db_uuid)

with open( outfile, 'w') as fd:
    json.dump(jdb,fd, sort_keys=True, indent=4,  ensure_ascii=False)

print("------------")
print("json files packed into '" + outfile + "' file")
