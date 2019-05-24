#! /usr/bin/python
#prend en argument la dim,d=2,3 puis le type d'element Pk,k=1 ou 2 pour d=2 et k=1 pour d=3

import os
import sys
from pyx import *
from pyx.deco import barrow,earrow
from pyx.style import linewidth, linestyle
from pyx.graph.axis import painter, tick
from pyx.graph.axis import *
from scipy import *
from scipy.optimize import leastsq
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--prefix", dest="prefix", help="prefix of the program")
parser.add_option("--program", dest="program", help="program to execute")
parser.add_option("--path", dest="path", help="path of the program to execute")
parser.add_option("--ntests", type="int",dest="ntests", help="number of tests", default=3)
parser.add_option("--hsize", type="float",dest="hsize", help="starting h size", default=.2)
parser.add_option("--dim", type="int",dest="dim", help="dimension of problem")
parser.add_option("--order", type="int",dest="order", help="order of the discretization")
(options, args) = parser.parse_args()


unit.set(xscale=2)
text.set(mode="latex")

def peval(x, p):
	return p[0]-p[1]*x**p[2];

def fit(x,y):
    p0 = [1, 1, 1]
    plsq = leastsq(residuals, p0, args=(y, x))
    return x,plsq[0]

def residuals(p, y, x):
    p1,p2,p3 = p
    err = y-p1+p2*x
    return err


mu=1
beta=1
c=options.order*options.order
#lesh1D = [1./4,1./8,1./16,1./32,1./64]
lesh1D = [1./2, 1./4, 1./8, 1./10, 1./20]
#lesh2D = [1./2,1./4,1./8,1./16,1./32]
lesh2D = [1./2,1./4,1./8,1./12,1./16]
#lesh3D = [1./2,1./4,1./8,1./10,1./12]
lesh3D = [ .25, .125, 0.1, .08, .05]
backend= '--backend=petsc'
print "cleanup resulTest.txt"
os.system('rm ~/feel/'+options.prefix+'/TestConv/'+str(options.dim)+'D/P'+str(options.order)+'/nu_1/beta_1/resultTest.txt')

for i in range(1,options.ntests+1) :
    print 'Execute '+str(options.program)+' avec h='+str(i)
    os.system(options.path+'/'+options.program+' --hsize='+str(options.hsize/i)+' '+backend)

print 'Recuperation des resultats'
os.system('cp ~/feel/'+options.prefix+'/TestConv/'+str(options.dim)+'D/P'+str(options.order)+'/nu_1/beta_1/resultTest.txt .')

data=io.array_import.read_array('resultTest.txt')

print 'Calcul de la droite des monidres carree'
NbDeTest=len(data[:,0])
xf1,yf1 = fit( log(data[:,0]),log(data[:,1]))
d=[[x,peval(x,yf1)] for x in xf1]
slope1=yf1[1]

xf2,yf2 = fit( log(data[:,0]),log(data[:,2]))
d=[[x,peval(x,yf2)] for x in xf2]
slope2=yf2[1]


print 'Creation du graph'
g = graph.graphxy(width=29,x=graph.axis.log(title="$h$"),
                  y=graph.axis.log(),key=graph.key.key())


g.plot(graph.data.file("resultTest.txt", x=1, y=2,title="Norme H1"),styles=[graph.style.line([color.rgb.red])])
g.plot(graph.data.file("resultTest.txt", x=1, y=3,title="Norme L2"),[graph.style.line([color.rgb.green])])
#g.plot(graph.data.function("y(x)=exp(x)"))
g.text(g.width/2, g.height + 0.2, "Test \ de \ Convergence : Laplacian "+str(options.dim)+"d P"+str(options.order)+" : \mu = "+str(mu)+" , beta = "+str(beta),[text.mathmode,text.halign.center, text.valign.bottom, text.size.Large])
g.finish()

x1, y1 = g.pos((data[NbDeTest/2+1,0]),(data[NbDeTest/2+1,1]))
x2, y2 = g.pos((data[NbDeTest/2-1,0]),(data[NbDeTest/2+1,1]))
x3, y3 = g.pos((data[NbDeTest/2-1,0]),(data[NbDeTest/2-1,1]))
g.stroke(path.line(x1-.5, y1, x2+.5, y2), [linestyle.dashed])
g.stroke(path.line(x2, y2-.5, x3, y3+.5), [linestyle.dashed])
g.stroke(path.line(x1, y1-.5, x1, y1+.5), [linestyle.dashed])
g.stroke(path.line(x3-.5, y3, x3+.5, y3), [linestyle.dashed])
g.stroke(path.line(x1,y1-.25,x2,y2-.25 ), [barrow.normal, earrow.normal])
g.stroke(path.line(x2+.25,y2,x3+.25,y3 ), [barrow.normal, earrow.normal])
g.text((x1+x2)/2,y1-.5, r"$a1$", [text.vshift.middlezero])
g.text(x2+.5,(y2+y3)/2, r"$b1$", [text.vshift.middlezero])
g.text(x1, y3, r"$b1/a1 = $"+str(abs(round(slope1,3))), [text.halign.center])



x1, y1 = g.pos((data[NbDeTest/2+1,0]),(data[NbDeTest/2+1,2]))
x2, y2 = g.pos((data[NbDeTest/2-1,0]),(data[NbDeTest/2+1,2]))
x3, y3 = g.pos((data[NbDeTest/2-1,0]),(data[NbDeTest/2-1,2]))
g.stroke(path.line(x1-.5, y1, x2+.5, y2), [linestyle.dashed])
g.stroke(path.line(x2, y2-.5, x3, y3+.5), [linestyle.dashed])
g.stroke(path.line(x1, y1-.5, x1, y1+.5), [linestyle.dashed])
g.stroke(path.line(x3-.5, y3, x3+.5, y3), [linestyle.dashed])
g.stroke(path.line(x1,y1-.25,x2,y2-.25 ), [barrow.normal, earrow.normal])
g.stroke(path.line(x2+.25,y2,x3+.25,y3 ), [barrow.normal, earrow.normal])
g.text((x1+x2)/2,y1-.5, r"$a2$", [text.vshift.middlezero])
g.text(x2+.5,(y2+y3)/2, r"$b2$", [text.vshift.middlezero])
g.text(x1, y3, r"$b2/a2 = $"+str(abs(round(slope2,3))), [text.halign.center])

# fenetre parametre

T1 ,T2 =g.pos(data[0,0], data[0,2])

PH1=options.order
PL2=PH1+1

tbox0 =text.text(T1, T2+0*unit.x_pt, r" Hypothesis on the slope :")
tbox1  = text.text(T1+5*3*unit.x_pt, T2-20*unit.x_pt, r"\left\| . \right\|_{H^{1}} \geq "+str(PH1),[text.mathmode, text.valign.bottom])
tbox2 = text.text(T1+5*3*unit.x_pt, T2+-40*unit.x_pt, r"\left\| . \right\|_{L^{2}}\geq"+str(PL2),[text.mathmode, text.valign.bottom])
tbox3 =text.text(T1, T2-60*unit.x_pt, r" Test Validation :")
tbox4T  = text.text(T1+5*3*unit.x_pt, T2-80*unit.x_pt, r"OK",[text.mathmode, text.valign.bottom])
tbox4F = text.text(T1+5*3*unit.x_pt, T2-80*unit.x_pt, r"PROBLEM",[text.mathmode, text.valign.bottom])


tpath1 = tbox1.bbox().enlarged(2*unit.x_pt).path()
tpath2 = tbox2.bbox().enlarged(2*unit.x_pt).path()
tpath4T = tbox4T.bbox().enlarged(2*unit.x_pt).path()
tpath4F = tbox4F.bbox().enlarged(2*unit.x_pt).path()


c = canvas.canvas()
c.draw(tpath1, [deco.filled([color.gray(0.6)]), deco.stroked()])
c.draw(tpath2, [deco.filled([color.gray(0.6)]), deco.stroked()])
c.insert(tbox0)
c.insert(tbox1)
c.insert(tbox2)
c.insert(tbox3)

MargeErreur=0.1

ok=0
if ( ((abs(slope1)+MargeErreur)>=PH1) & ((abs(slope2)+MargeErreur)>=PL2) ) :
	c.draw(tpath4T, [deco.filled([color.cmyk.Green]), deco.stroked()])
	c.insert(tbox4T)
	print 'Test OK'
else :
	c.draw(tpath4F, [deco.filled([color.cmyk.Red]), deco.stroked()])
	c.insert(tbox4F)
	print 'Test FAILED'
	ok=1


g.insert(c)

g.writePDFfile('convergence_'+str(options.dim)+'D_P'+str(options.order))
g.writeEPSfile('convergence_'+str(options.dim)+'D_P'+str(options.order))

os.system('rm resultTest.txt')

sys.exit(ok)
