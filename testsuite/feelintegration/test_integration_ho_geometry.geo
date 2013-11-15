// -*- mode: c++ -*-
//
//  This file is part of the Feel library
//
//  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
//       Date: 2008-09-15
//
//  Copyright (C) 2008 Université de Grenoble et CNRS, Laboratoire Jean  Kuntzmann
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//
//
//  \file test_integration_ho.geo
//  \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
//  \date 2008-09-15
//
/*
h4=h;
Point(1) = { 0.0, 0.0  , 0.0 , h4};
Point(2) = { 1.0, 0.0  , 0.0 , h4};
//Point(3) = { 1.0, Exp(1.0)  , 0.0 , h};
Point(3) = { 1.0, 1.0  , 0.0 , h};
Point(4) = { 0.0, 1.0  , 0.0 , h};
Point(5) = { 0.5, 1.5  , 0.0 , h4};
Point(6) = { 0.5, 1.0  , 0.0 , h4};

Line(20)  = {1 ,4};
Line(21)  = {3 ,2};
Line(22)  = {2 ,1};

N=3;
//N=2;
amp=.25;

pts[0]=4;
pts[N]=3;



For t In {1:N}

  If (N > 1 && t < N)
    pts[t]=newp;
//Point(pts[t]) = { 0+t/N, 1+amp*Sin(Pi*t/N), 0., h };
//Point(pts[t]) = { 0+t/N, Exp(t/N)*(1+amp*Sin(2*Pi*t/N)), 0., h };
Point(pts[t]) = { 0+t/N, 1+t/N*(1-t/N)*(10-t/N), 0., h };
  EndIf
  //  lines[t-1]=newreg;
  //  Line(lines[t-1])={pts[t-1],pts[t]};
  //Printf("segment %g, pts[%g]=%g, pts[%g]=%g, line %g",t,t-1,pts[t-1],t,pts[t],lines[t-1]);
EndFor


Spline(16)={pts[]};
Line Loop(15) = {20,16,21,22};


Circle(16)={4,6,5};
Circle(17)={5,6,3};
Line Loop(15) = {20,16,17,21,22};


Plane Surface(30) = {-15};

Physical Line(1) = {15};
Physical Surface(1) = {30};
*/

h4=h;
Point(1) = { 0.0, 0.0  , 0.0 , h4};
Point(2) = { 1.0, 0.0  , 0.0 , h4};
Point(3) = { 0.0, 1.0  , 0.0 , h};
Point(4) = {-1.0, 0.0  , 0.0 , h};
Point(5) = { 0.0,-1.0  , 0.0 , h4};


Circle(20)  = {2 ,1, 3};
Circle(21)  = {3 ,1, 4};
Circle(22)  = {4 ,1, 5};
Circle(23)  = {5 ,1, 2};

Line Loop(24) = {20,21,22,23};

//Line(31) = {1, 2};
//Line(32) = {3, 1};
//Line Loop(24) = {31,20,32};
Plane Surface(30) = {24};

Physical Line(1) = {24};
Physical Surface(1) = {30};

View "comments" {
  // 10 pixels from the left and 15 pixels from the top of the graphic
  // window:
  T2(10,15,0){StrCat("File created on ", Today)};

  T2(10,-30,0){Sprintf("h=%.2f, Order=%.0f, Version=%.0f",h,Mesh.ElementOrder,Mesh.MshFileVersion)};
  // 10 pixels from the left and 10 pixels from the bottom of the
  // graphic window:
  T2(10,-10,0){"Copyright (C) Université de Grenoble et CNRS, Laboratoire Jean  Kuntzmann"};

};
