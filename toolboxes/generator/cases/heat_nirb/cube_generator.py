from liquid import Environment
from liquid import FileSystemLoader

import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--Nv", help="number of fins in vertical direction [default=3]", type=str, default="3")
parser.add_argument("--Nh", help="number of fins in horizontal direction [default=3]", type=str, default="3")
parser.add_argument("--Nd", help="number of fins in depth direction [default=3]", type=str, default="3")
parser.add_argument("--L", help="width of a cube [default=1]", type=str, default="1")
parser.add_argument("--h", help="height of a cube [default=1]", type=str, default="1")
parser.add_argument("--d", help="depth of a cube (for 3D only)", type=str, default="1")
parser.add_argument("--dim", help="dimension of the case (2 or 3) [default=2]", type=str, default="2")
parser.add_argument("--odir", help="output directory", type=str, default=".")


args = parser.parse_args()

if not args.odir[-1] == "/":
    args.odir += "/"

if not os.path.isdir(args.odir):
    os.mkdir(args.odir)



if args.dim not in ["2","3"]:
    raise ValueError("dimension must be 2 or 3")


DIRPATH = os.path.dirname(__file__)
if len(DIRPATH) != 0:
    DIRPATH += '/'

templates = ["cube_templates/"]

meshfile = "cube.geo"
modelfile = "heat-cube.json"

env = Environment(loader=FileSystemLoader(DIRPATH + templates[0]))
templateGeo = env.get_template("cube.geo")
templateCfg = env.get_template("heat-cube.cfg")
templateJson = env.get_template("heat-cube.json")



renderGeo = ""
renderCfg = ""

if args.dim == "2":
    ElementShape = "Rectangle"
    ElementArgs = "{ s*L, r*h, 0, L, h, 0}"
    eltDim = "Surface"
    eltDimM1 = "Curve"
    diffVal = 1
    fourierVal = 4

else:
    ElementShape = "Box"
    ElementArgs = "{ s*L, r*h, 0, L, h, d}"
    eltDim = "Volume"
    eltDimM1 = "Surface"
    diffVal = 3
    fourierVal = 1
     

    

renderGeo = templateGeo.render(
    Nv = args.Nv,
    Nh = args.Nh,
    L = args.L,
    h = args.h,
    d = args.d,
    ElementShape = ElementShape,
    ElementArgs = ElementArgs, 
    eltDim = eltDim,
    eltDimM1 = eltDimM1,
    diffVal = diffVal,
    dim = args.dim,
    fourierVal = fourierVal
)

renderCfg = templateCfg.render(
    dim = args.dim
)

fins = list(range(2,int(args.Nv)*int(args.Nh)+1))

renderJson = templateJson.render(
    fins = fins,
    dim = args.dim,
    geofilename = meshfile,
    modelfilename = modelfile
)


fileGeo = open(args.odir+"cube.geo", "w")
a = fileGeo.write(renderGeo)
fileGeo.close()

fileCfg = open(args.odir+"heat-cube.cfg", "w")
a = fileCfg.write(renderCfg)
fileCfg.close()


fileJson = open(args.odir+"heat-cube.json", "w")
a = fileJson.write(renderJson)
fileJson.close()
