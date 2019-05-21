from pyx import *
import errorplot
import slopes
from style import graph_style

text.set(mode="latex")
text.preamble(r"\usepackage{amssymb}")

revision = "4497"

xmin = 1
xmax = 1e6
ymin = 1e-4
ymax = 1e4

c = canvas.canvas()
g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
                        x=graph.axis.log(title="nDof",min=xmin,max=xmax),
                        y=graph.axis.log(title="$time$",min=ymin,max=ymax),
                        key=graph.key.key(pos="tl", dist=0.05)))


s1=slopes.getSlope(g, "$2", [7], "P1_timings_" + revision + ".txt")

g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="1/$2", y=4,title="stabilisation assembly"), graph_style)
g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="$2", y=5,title="stokes assembly"), graph_style)
g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="$2", y=6,title="global assembly"), graph_style)
g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="$2", y=7,title="rhs assembly"), graph_style)
g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="$2", y=8,title="boundary conditions"), graph_style)
g.plot(graph.data.file("P1_timings_" + revision + ".txt", x="$2", y=9,title="calculate $||u - u_h ||_{L^2(\Omega)}$"), graph_style)

c.writePDFfile("timings/P1_timings_" + revision + "")










xmin = 1
xmax = 1e6
ymin = 1e-4
ymax = 1e4

c = canvas.canvas()
g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
                        x=graph.axis.log(title="nDof",min=xmin,max=xmax),
                        y=graph.axis.log(title="$time$",min=ymin,max=ymax),
                        key=graph.key.key(pos="tl", dist=0.05)))


#s1=slopes.getSlope(g, "$2", [7], "P2_timings_" + revision + ".txt")

g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=4,title="stabilisation assembly"), graph_style)
g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=5,title="stokes assembly"), graph_style)
g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=6,title="global assembly"), graph_style)
g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=7,title="rhs assembly"), graph_style)
g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=8,title="boundary conditions"), graph_style)
g.plot(graph.data.file("P2_timings_" + revision + ".txt", x="$2", y=9,title="calculate $||u - u_h ||_{L^2(\Omega)}$"), graph_style)

c.writePDFfile("timings/P2_timings_" + revision + "")















xmin = 1
xmax = 1e6
ymin = 1e-4
ymax = 1e4

xmin_slopes = 0.1
xmax_slopes = 0.8

c = canvas.canvas()
g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
                        x=graph.axis.log(title="nDof",min=xmin,max=xmax),
                        y=graph.axis.log(title="$time$",min=ymin,max=ymax),
                        key=graph.key.key(pos="tl", dist=0.03)))


#s1=slopes.getSlope(g, "$2", [7], "P5_timings_" + revision + ".txt")

g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=4,title="stabilisation assembly"), graph_style)
g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=5,title="stokes assembly"), graph_style)
g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=6,title="global assembly"), graph_style)
g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=7,title="rhs assembly"), graph_style)
g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=8,title="boundary conditions"), graph_style)
g.plot(graph.data.file("P5_timings_" + revision + ".txt", x="$2", y=9,title="calculate $||u - u_h ||_{L^2(\Omega)}$"), graph_style)

c.writePDFfile("timings/P5_timings_" + revision + "")
