#! /usr/bin/python

# -*- mode: python -*-
#
#  This file is part of the Feel library
#
#  Author(s): Goncalo Pena <gpena@mat.uc.pt>
#        Date: 2010-01-09
#
#   Copyright (C) 2010 University of Coimbra
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#   University of Coimbra
#
# \file simplex.py
# \author Goncalo Pena <gpena@mat.uc.pt>
# \date 2010-01-09
#

# This file contains the functions to draw the contours of a simplex,
# draw the equispaced pointset up to any order,
# and draw the numbering of the points (in the pointset) by subentity


from pyx import *
import math
from scipy import *



def border(c, x0, line_size):
    a, b = x0
    result = []

    #border of triangle
    tri = path.path( path.moveto(a, b), path.lineto(a+line_size, b),
                     path.lineto(a, b+line_size),
                     path.closepath())

    c.stroke(tri, [style.linewidth(0.04)])

    return result


def Vertex(c, line_size, origin, local_id, circle_radius):

    result = []

    x, y = origin

    if local_id == 1:
	x = x + line_size
    else:
	if local_id == 2:
	    y = y + line_size

    c.stroke(path.circle(x, y, circle_radius), [style.linewidth.thin, deco.filled([color.rgb.red])])

    return result


def Edge(c, line_size, origin, local_id, order, circle_radius):
    # plots the points in an edge of a triangle
    result = []

    h = 1.0*line_size/order

    x0 = origin[0]
    x1 = origin[1]

    if local_id == 0:
	x0 = x0 + line_size
	v = [-1,1]
    else:
	if local_id == 1:
	    x1 = x1 + line_size
    	    v=[0,-1]
        else:
	    v=[1,0]

    for i in range(order-1):
	p0 = x0 + (i+1)*h*v[0]
	p1 = x1 + (i+1)*h*v[1]

	c.stroke(path.circle(p0, p1, circle_radius), [style.linewidth.thin, deco.filled([color.rgb.blue])])


    return result


def Face(c, line_size, origin, order, circle_radius):

    result = []

    if order >= 3:
    	h = 1.0*line_size/(order)

    	for i in range(order):
           for j in range(order-i-2):
	       p0 = origin[0] + 1.0*h*( (i+1)  )
	       p1 = origin[1] + 1.0*h*( (j+1) )

	       c.stroke(path.circle(p0, p1, circle_radius), [style.linewidth.thin, deco.filled([color.rgb.green])])


    return result




def numberVertex(c, line_size, origin, local_id, circle_radius, shift):

    result = []

    x = origin[0] - shift
    y = origin[1] - 3.0*shift

    if local_id == 0:
	c.text(x, y, local_id, [text.halign.boxright])

    if local_id == 1:
	x = x + line_size + 2.0*shift
        c.text(x, y, local_id, [text.halign.boxleft])

    if local_id == 2:
	y = y + line_size + 4.0*shift
        c.text(x, y, local_id, [text.halign.boxright])

    return result



def numberEdge(c, line_size, origin, local_id, order, circle_radius, shift):

    result = []

    x = origin[0] - shift
    y = origin[1]

    h = 1.0*line_size/order


    if local_id == 0:
	x = x + line_size + 3.0*shift
	v = [-1,1]
	first = 3
    else:
	if local_id == 1:
	    y = y + line_size - 1.0*shift
	    x = x - 1.0*shift
    	    v=[0,-1]
            first = order+2
        else:
	    x = x + 1.0*shift
	    y = y - 4.0*shift
	    v=[1,0]
	    first = 2*order+1

    for i in range(order-1):
	p0 = x + (i+1)*h*v[0]
	p1 = y + (i+1)*h*v[1]

	if local_id == 0:
	    c.text(p0, p1, i+first, [text.halign.boxleft])
	else:
	    if local_id == 1:
		c.text(p0, p1, i+first, [text.halign.boxright])
	    else:
		c.text(p0, p1, i+first, [text.halign.boxcenter])


    return result





def numberFace(c, line_size, origin, order, circle_radius, shift):

    result = []

    x = origin[0]
    y = origin[1]



    if order >= 3:
    	h = 1.0*line_size/(order)

    p = 3 + 3*(order-1)

    for i in range(order):
        for j in range(order-i-2):
	    p0 = origin[0] + 1.0*h*( j+1 ) - 2.0*shift
	    p1 = origin[1] + 1.0*h*( i+1 ) - 1.0*shift

	    c.text(p0, p1, p, [text.halign.boxright])
	    p = p+1


    return result

def n_side_points(order):
    return order-2

def npoints(order):
    n = n_side_points(order)
    return n*(n+1)/2


def identity(order, i):
    return i

def base(order, index):
    N = n_side_points(order)

    for i in range(N):
	for j in range(i+1):
	    first = i + j*(j-1)/2 + j*(N - j);
	    last = i + (i-j)*(i-j-1)/2 + (i-j)*(N - (i-j))

	    if first == index:
		return last

def height(order, index):
    N = n_side_points(order)

    for i in range(N):
	first = i*N - i*(i-1)/2
	last = (i+1)*N - (i+1)*i/2 - 1

	for j in range(N-i):
	    if index == first+j:
		return last-j

def anticlockwise(order, index):
    p = height(order, index)
    return base(order, p)


def hypotenuse(order, index):
    p = base(order, index)
    return anticlockwise(order, p)

def clockwise(order, index):
    p = base(order, index)
    return height(order, p)



def addIndices(c, index, line_size, origin, order, circle_radius, shift):

    result = []

    x = origin[0]
    y = origin[1]

    if order >= 3:
    	h = 1.0*line_size/(order)

    p = 0
    index_shift = 0
    #3 + 3*(order-1)


    for i in range(order):
        for j in range(order-i-2):
	    p0 = origin[0] + 1.0*h*( j+1 ) - 2.0*shift
	    p1 = origin[1] + 1.0*h*( i+1 ) - 1.0*shift

	    if index == 0:
	    	c.text(p0, p1, index_shift + identity(order, p), [text.halign.boxright])

	    if index == 1:
	    	c.text(p0, p1, index_shift + base(order, p), [text.halign.boxright])

	    if index == 2:
	    	c.text(p0, p1, index_shift + height(order, p), [text.halign.boxright])

	    if index == 3:
	    	c.text(p0, p1, index_shift + hypotenuse(order, p), [text.halign.boxright])

	    if index == 4:
	    	c.text(p0, p1, index_shift + anticlockwise(order, p), [text.halign.boxright])

	    if index == 5:
	    	c.text(p0, p1, index_shift + clockwise(order, p), [text.halign.boxright])

	    p = p+1
    return result


def addLegend(c, index, x0, hpos, vpos):
    result = []

    p0 = x0[0] + hpos
    p1 = x0[1] - vpos

    if index == 0:
    	c.text(p0, p1, "$\sigma_0$", [text.halign.boxcenter])
    if index == 1:
    	c.text(p0, p1, "$\sigma_1$", [text.halign.boxcenter])
    if index == 2:
    	c.text(p0, p1, "$\sigma_2$", [text.halign.boxcenter])
    if index == 3:
    	c.text(p0, p1, "$\sigma_3$", [text.halign.boxcenter])
    if index == 4:
    	c.text(p0, p1, "$\sigma_4$", [text.halign.boxcenter])
    if index == 5:
    	c.text(p0, p1, "$\sigma_5$", [text.halign.boxcenter])

    return result

# Local Variables:
# indent-tabs-mode: t
# End:
