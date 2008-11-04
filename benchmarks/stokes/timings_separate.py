from pyx import *
#import errorplot
#import slopes
from style import graph_style


graph_style = [graph.style.symbol(graph.style.symbol.changesquare, 0.1, [deco.filled]),
               graph.style.line([style.linestyle.dashed, style.linewidth.thin])]

#orders = ["1", "2", "5"]
orders = ["1"]
for i in range(3):
	order_str = orders[i]

	xmin = 1e0
	xmax = 1e6
	ymin = 1e-2
	ymax = 100

	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=4,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=4,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=4,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_stabilisation")













	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=5,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=5,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=5,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_stokes")






	ymin = 1e-4
	ymax = 1e1

	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=6,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=6,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=6,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_global_assembly")










	ymin = 1e-3
	ymax = 1e3

	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=7,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=7,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=7,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_rhs")











	ymin = 1e-3
	ymax = 1

	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=8,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=8,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=8,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_boundary_conditions")











	ymin = 1e-3
	ymax = 1e2

	c = canvas.canvas()
	g = c.insert(graph.graphxy(width=8, ratio=4./4.5,
				x=graph.axis.log(title="nDof",min=xmin,max=xmax),
				y=graph.axis.log(title="$time$",min=ymin,max=ymax),
				key=graph.key.key(pos="tl", dist=0.05)))

	#g.plot(graph.data.file("P" + order_str + "_timings_0.txt", x="$2", y=9,title="revision 0"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4362.txt", x="$2", y=9,title="revision 4362"), graph_style)
	g.plot(graph.data.file("P" + order_str + "_timings_4497.txt", x="$2", y=9,title="revision 4497"), graph_style)

	c.writePDFfile("timings/P" + order_str + "_timings_L2_error_velocity")
