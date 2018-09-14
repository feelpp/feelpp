/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-01-09

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file geom.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-01-09
 */
#include <sstream>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
namespace Feel
{
gmsh_ptrtype
oneelement_geometry_ref()
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=2;\n"
          <<"Point(1) = {-1,-1,0,h};\n"
          <<"Point(2) = {1,-1,0,h};\n"
          <<"Point(3) = {-1,1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n"
          <<"Transfinite Line{1} = 1;\n"
          <<"Transfinite Line{2} = 1;\n"
          <<"Transfinite Line{3} = 1;\n"
          <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;
    nameStr << "one-elt-ref";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}


gmsh_ptrtype
oneelement_geometry_real()
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=2;\n"
          <<"Point(1) = {-2,-2,0,h};\n"
          <<"Point(2) = {2,-2,0,h};\n"
          <<"Point(3) = {-2,2,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n"
          <<"Transfinite Line{1} = 1;\n"
          <<"Transfinite Line{2} = 1;\n"
          <<"Transfinite Line{3} = 1;\n"
          <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;
    nameStr << "one-elt-real";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

}
