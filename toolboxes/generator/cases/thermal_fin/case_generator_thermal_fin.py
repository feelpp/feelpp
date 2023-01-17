from liquid import Environment
from liquid import FileSystemLoader

import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--N", help="number of fins [default=4]", type=str, default="4")
parser.add_argument("--L", help="width of a fin [default=2.5]", type=str, default="2.5")
parser.add_argument("--d", help="distance between two fins [default=0.75]", type=str, default="0.75")
parser.add_argument("--t", help="thickness of a fin [default=0.25]", type=str, default="0.25")
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

DIRPATH = os.path.dirname(__file__)
if len(DIRPATH) != 0:
    DIRPATH += '/'


env = Environment(loader=FileSystemLoader(DIRPATH + "templates/"))
templateGeo = env.get_template("fin.geo")
templateCfg = env.get_template("thermal-fin.cfg")
templateJson = env.get_template("thermal-fin.json")



renderGeo = ""
renderCfg = ""

if args.dim == "2":
    PostShape = "Rectangle"
    PostArgs = "{0, 0, 0, 1, N*(d+t)+t, 0}"
    FinShape = "Rectangle"
    FinArgs = "{-L, r*(d+t), 0, 2*L+1, t, 0}"
    eltDim = "Surface"
    eltDimM1 = "Curve"
    step = 2
    physicalArg = "{r,r+1}"
    diffVal = 4

else:
    eltDim = "Volume"
    eltDimM1 = "Surface"
    
    if args.cylinder <= 1:
        PostShape = "Box"
        PostArgs = "{0, 0, 0, 1, 1, N*(d+t)+t}"
    else:
        PostShape = "Cylinder"
        PostArgs = "{0.5, 0.5, 0, 0, 0, N*(d+t)+t, 0.5, 2*Pi}"

    step = 1
    physicalArg = "{ r }"
    
    if args.cylinder >= 1:
        FinShape = "Cylinder"
        FinArgs = "{0.5, 0.5, r*(d+t), 0, 0, t, L, 2*Pi}"
        
    else:
        FinShape = "Box"
        FinArgs = "{-L, -L, r*(d+t), 2*L+1, 2*L+1, t}"

    diffVal = [5, 5, 17][args.cylinder]
    
renderGeo = templateGeo.render(
    N = args.N,
    L = args.L,
    t = args.t,
    d = args.d,
    PostShape = PostShape,
    PostArgs = PostArgs,
    FinShape = FinShape,
    FinArgs = FinArgs,
    step = step,
    physicalArg = physicalArg,
    eltDim = eltDim,
    eltDimM1 = eltDimM1,
    diffVal = diffVal,
)

renderCfg = templateCfg.render(
    dim = args.dim
)

fins = list(range(1,int(args.N)+1))

renderJson = templateJson.render(
    fins = fins,
    dim = args.dim
)


fileGeo = open(args.odir+"fin.geo", "w")
a = fileGeo.write(renderGeo)
fileGeo.close()

fileCfg = open(args.odir+"thermal-fin.cfg", "w")
a = fileCfg.write(renderCfg)
fileCfg.close()


fileJson = open(args.odir+"thermal-fin.json", "w")
a = fileJson.write(renderJson)
fileJson.close()
