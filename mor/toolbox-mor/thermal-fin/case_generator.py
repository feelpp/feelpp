from liquid import Environment, template
from liquid import FileSystemLoader

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--N", help="number of fins [default=4]", type=str, default="4")
parser.add_argument("--L", help="width of a fin [default=2.5]", type=str, default="2.5")
parser.add_argument("--d", help="distance between two fins [default=0.5]", type=str, default="0.5")
parser.add_argument("--t", help="thickness of a fin [default=0.25]", type=str, default="0.25")
parser.add_argument("--P", help="tickness of the thermal fin (only for 3D case) [default=1]", type=str, default="1")
parser.add_argument("--dim", help="dimension of the case (2 or 3) [default=2]", type=str, default="2")
parser.add_argument("--cylinder", help="shape of fin and post (0=boxes, 1=box/cylinders, 2=cylinders) [default=0]", type=int, default=0)

args = parser.parse_args()



if args.dim not in ["2","3"]:
    raise ValueError("dimension must be 2 or 3")

if args.dim == "3" and args.cylinder not in [0, 1, 2]:
    raise ValueError("cylinder must be 0, 1 or 2")



env = Environment(loader=FileSystemLoader("templates/"))
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
    diffVal = 4
    Elt = "N"

else:
    eltDim = "Volume"
    eltDimM1 = "Surface"
    
    if args.cylinder <= 1:
        PostShape = "Box"
        PostArgs = ["{0, 0, 0, 1, P, N*(d+t)+t}", "{0, 0, 0, 1, t, N*(d+t)+t}"][args.cylinder]
    else:
        PostShape = "Cylinder"
        PostArgs = "{t/2, t/2, 0, 0, 0, N*(d+t)+t, t/2, 2*Pi}"
    
    if args.cylinder >= 1:
        FinShape = "Cylinder"
        FinArgs = "{-L, t/2, r*(d+t)+t/2, 2*L+1, 0, 0, t/2, 2*Pi}"
        Elt = "(N+1)"
    else:
        FinShape = "Box"
        FinArgs = "{-L, 0, r*(d+t), 2*L+1, P, t}"
        Elt = "N"

    diffVal = [5, 5, 33][args.cylinder]
    
renderGeo = templateGeo.render(
    P = args.P,
    R = args.R,
    N = args.N,
    L = args.L,
    t = args.t,
    d = args.d,
    PostShape = PostShape,
    PostArgs = PostArgs,
    FinShape = FinShape,
    FinArgs = FinArgs,
    eltDim = eltDim,
    eltDimM1 = eltDimM1,
    diffVal = diffVal,
    Elt = Elt
)

renderCfg = templateCfg.render(
    dim = args.dim
)

fins = list(range(1,int(args.N)+1))

renderJson = templateJson.render(
    fins = fins,
    dim = args.dim
)


fileGeo = open("fin.geo", "w")
a = fileGeo.write(renderGeo)
fileGeo.close()

fileCfg = open("thermal-fin.cfg", "w")
a = fileCfg.write(renderCfg)
fileCfg.close()


fileJson = open("thermal-fin.json", "w")
a = fileJson.write(renderJson)
fileJson.close()
