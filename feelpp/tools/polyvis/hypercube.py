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
# \file hypercube.py
# \author Goncalo Pena <gpena@mat.uc.pt>
# \date 2010-01-09
#

# This file contains the functions to draw the contours of a hypercube,
# draw the equispaced pointset up to any order,
# and draw the numbering of the points (in the pointset) by subentity

from pyx import *
import math
from scipy import *


def line(c, x0, x1):
    a, b = x0
    d, e = x1
    result = []

    sqa = path.path( path.moveto(a, b), path.lineto(d,e), path.closepath())

    c.stroke(sqa, [style.linewidth(0.04)])

    return result


def border(c, x0, line_size):
    a, b = x0
    result = []

    #border of square

    sqa = path.path( path.moveto(a, b), path.lineto(a+line_size, b),
                     path.lineto(a+line_size, b+line_size), path.lineto(a, b+line_size),
                     path.closepath())

    c.stroke(sqa, [style.linewidth(0.04)])

    return result


def borderRect(c, x0, line_size1, line_size2):
    a, b = x0
    result = []

    #border of square

    sqa = path.path( path.moveto(a, b), path.lineto(a+line_size1, b),
                     path.lineto(a+line_size1, b+line_size2), path.lineto(a, b+line_size2),
                     path.closepath())

    c.stroke(sqa, [style.linewidth(0.04)])

    return result


def Vertex(c, line_size, origin, local_id, circle_radius):

    result = []

    x = origin[0]
    y = origin[1]

    if local_id == 1:
	x = x + line_size
    else:
	if local_id == 2:
	    x = x + line_size
	    y = y + line_size
	else:
	    if local_id == 3:
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
	v = [1,0]
    else:
	if local_id == 1:
	    x0 = x0 + line_size
    	    v=[0,1]
        else:
	    if local_id == 2:
		x0 = x0 + line_size
		x1 = x1 + line_size
	        v=[-1,0]
	    else:
		x1 = x1 + line_size
	        v=[0,-1]


    for i in range(order-1):
	p0 = x0 + (i+1)*h*v[0]
	p1 = x1 + (i+1)*h*v[1]

	c.stroke(path.circle(p0, p1, circle_radius), [style.linewidth.thin, deco.filled([color.rgb.blue])])


    return result


def Face(c, line_size, origin, order, circle_radius):

    result = []

    if order >= 2:
    	h = 1.0*line_size/(order)

    	for i in range(order-1):
           for j in range(order-1):
	       p0 = origin[0] + 1.0*h*( (i+1)  )
	       p1 = origin[1] + 1.0*h*( (j+1) )

	       c.stroke(path.circle(p0, p1, circle_radius), [style.linewidth.thin, deco.filled([color.rgb.green])])

    return result


def numberLine(c, line_size, origin, order, circle_radius, shift):

    result = []

    local_id = 0

    x = origin[0]
    y = origin[1]

    h = 1.0*line_size/order

    x = origin[0] - shift
    y = origin[1] - 4.0*shift
    v = [1,0]
    first = 2

    for i in range(order-1):
        p0 = x + (i+1)*h*v[0]
        p1 = y + (i+1)*h*v[1]

        c.text(p0, p1, i+first, [text.halign.boxleft])

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
	x = x + line_size + 4.0*shift
	y = y + line_size + 4.0*shift
        c.text(x, y, local_id, [text.halign.boxright])

    if local_id == 3:
	y = y + line_size + 4.0*shift
        c.text(x, y, local_id, [text.halign.boxright])

    return result

def numberEdge(c, line_size, origin, local_id, order, circle_radius, shift):

    result = []

    x = origin[0]
    y = origin[1]

    h = 1.0*line_size/order

    if local_id == 0:
	x = x - shift
	y = y - 4.0*shift
	v = [1,0]
	first = 4
    else:
	if local_id == 1:
	    x = x + line_size + 2.0*shift
	    y = y - shift
    	    v=[0,1]
	    first = order + 3
        else:
	    if local_id == 2:
		x = x + line_size
		y = y + line_size + 2.0*shift
	        v=[-1,0]
		first = 2*order +2
	    else:
		x = x - 2.0*shift
		y = y + line_size - shift
	        v=[0,-1]
		first = 3*order +1


    for i in range(order-1):
	p0 = x + (i+1)*h*v[0]
	p1 = y + (i+1)*h*v[1]

	if local_id == 0:
	    c.text(p0, p1, i+first, [text.halign.boxleft])
	else:
	    if local_id == 1:
		c.text(p0, p1, i+first, [text.halign.boxleft])
	    else:
		if local_id == 2:
		    c.text(p0, p1, i+first, [text.halign.boxcenter])
		else:
		    c.text(p0, p1, i+first, [text.halign.boxright])


    return result

def numberFace(c, line_size, origin, order, circle_radius, shift):

    result = []

    x = origin[0]
    y = origin[1]

    if order >= 2:
    	h = 1.0*line_size/(order)


        p = 4 + 4*(order-1)

        for i in range(order-1):
            for j in range(order-1):
	        p0 = origin[0] + 1.0*h*( i+1 ) - 2.0*shift
	        p1 = origin[1] + 1.0*h*( j+1 ) - 1.0*shift

	        c.text(p0, p1, p, [text.halign.boxright])
	        p = p+1


    return result


def n_side_points(order):
    return order-1

def identity(order, i):
    return i

def base(order, index):
    N = n_side_points(order)

    p = 0

    result = index

    for i in range(N):
	k = N-1-i

	first = k * N

	for j in range(N):
	    if index == p:
		result = first+j

   	    p = p+1

    return result


def once_anticlockwise(order, index):
    N = n_side_points(order)

    p = 0

    for i in range(N):
	k = N-1-i

	for j in range(N):
	    if index == p:
		return k + N*j

   	    p = p+1


def second_diagonal(order, index):
    p = base(order, index)
    return once_anticlockwise(order, p)


def height(order, index):
    p = second_diagonal(order, index)
    return once_anticlockwise(order, p)

def twice_clockwise(order, index):
    p = base(order, index)
    return height(order, p)

def diagonal(order, index):
    p = height(order, index)
    return once_anticlockwise(order, p)

def once_clockwise(order, index):
    p = once_anticlockwise(order, index)
    return twice_clockwise(order, p)

def addIndices(c, index, line_size, origin, order, circle_radius, shift):

    result = []

    x = origin[0]
    y = origin[1]

    if order >= 2:
    	h = 1.0*line_size/(order)

    p = 0
    index_shift = 0 #4 + 4*(order-1)


    for i in range(order-1):
        for j in range(order-1):
	    p0 = origin[0] + 1.0*h*( i+1 ) - 2.0*shift
	    p1 = origin[1] + 1.0*h*( j+1 ) - 1.0*shift

	    if index == 0:
	    	c.text(p0, p1, index_shift + identity(order, p), [text.halign.boxright])

	    if index == 1:
	    	c.text(p0, p1, index_shift + second_diagonal(order, p), [text.halign.boxright])

	    if index == 2:
	    	c.text(p0, p1, index_shift + base(order, p), [text.halign.boxright])

	    if index == 3:
	    	c.text(p0, p1, index_shift + once_anticlockwise(order, p), [text.halign.boxright])

	    if index == 4:
	    	c.text(p0, p1, index_shift + diagonal(order, p), [text.halign.boxright])

	    if index == 5:
	    	c.text(p0, p1, index_shift + twice_clockwise(order, p), [text.halign.boxright])
	    if index == 6:
	    	c.text(p0, p1, index_shift + once_clockwise(order, p), [text.halign.boxright])
	    if index == 7:
	    	c.text(p0, p1, index_shift + height(order, p), [text.halign.boxright])
	    p = p+1
    return result

















def addLegend(c, index, x0, hpos, vpos):
    result = []

    p0 = x0[0] + hpos
    p1 = x0[1] - vpos

    if index == 0:
    	c.text(p0, p1, "$\omega_0$", [text.halign.boxcenter])
    if index == 1:
    	c.text(p0, p1, "$\omega_1$", [text.halign.boxcenter])
    if index == 2:
    	c.text(p0, p1, "$\omega_2$", [text.halign.boxcenter])
    if index == 3:
    	c.text(p0, p1, "$\omega_3$", [text.halign.boxcenter])
    if index == 4:
    	c.text(p0, p1, "$\omega_4$", [text.halign.boxcenter])
    if index == 5:
    	c.text(p0, p1, "$\omega_5$", [text.halign.boxcenter])
    if index == 6:
    	c.text(p0, p1, "$\omega_6$", [text.halign.boxcenter])
    if index == 7:
    	c.text(p0, p1, "$\omega_7$", [text.halign.boxcenter])


    return result

# Local Variables:
# indent-tabs-mode: t
# End:

