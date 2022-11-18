from liquid import Environment
from liquid import FileSystemLoader

import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--Nv", help="number of fins in vertical direction [default=3]", type=str, default="3")
parser.add_argument("--Nh", help="number of fins in horizontal direction [default=3]", type=str, default="3")
parser.add_argument("--L", help="width of a fin [default=1]", type=str, default="1")
parser.add_argument("--t", help="thickness of a fin [default=1]", type=str, default="1")
parser.add_argument("--dim", help="dimension of the case (2 or 3) [default=2]", type=str, default="2")
parser.add_argument("--cylinder", help="shape of fin and post (0=boxes, 1=box/cylinders, 2=cylinders) [default=0]", type=int, default=0)
parser.add_argument("--odir", help="output directory", type=str, default=".")


args = parser.parse_args()

if not args.odir[-1] == "/":
    args.odir += "/"

if not os.path.isdir(args.odir):
    os.mkdir(args.odir)



if args.dim not in ["2","3"]:
    raise ValueError("dimension must be 2 or 3")

if args.dim == "3" and args.cylinder not in [0, 1, 2]:
    raise ValueError("cylinder must be 0, 1 or 2")

print(__file__)

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
    ElementArgs = "{ -L + (s-1)*L, -t + (r-1)*t, 0, L, t, 0}"
    eltDim = "Surface"
    eltDimM1 = "Curve"
    step = 1
    diffVal = 1

else:
    print(f"Dim {args.dim} not implemented")

    exit()
     

    

renderGeo = templateGeo.render(
    Nv = args.Nv,
    Nh = args.Nh,
    L = args.L,
    t = args.t,
    ElementShape = ElementShape,
    ElementArgs = ElementArgs, 
    step = step,
    eltDim = eltDim,
    eltDimM1 = eltDimM1,
    diffVal = diffVal,
)

renderCfg = templateCfg.render(
    dim = args.dim
)

fins = list(range(1,int(args.Nv)*int(args.Nh)+1))

renderJson = templateJson.render(
    fins = fins,
    dim = args.dim
    geofilename = geofilename
    modelfilename = modelfilename
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
