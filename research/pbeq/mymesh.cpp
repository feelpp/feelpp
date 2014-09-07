/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2007-07-17

  Copyright (C) 2007 Unil

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file mesh.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2007-07-17
 */

#include "mymesh.hpp"

namespace Feel
{



std::string
myMesh::getDescription( ) const
{
    std::ostringstream ostr;

    ostr <<"Mesh.MshFileVersion = 1;"
         <<"Function mySphere \n"
         <<"  // In the following commands we use the reserved variable name \n"
         <<"  // `newp', which automatically selects a new point number. This \n"
         <<"  // number is chosen as the highest current point number, plus \n"
         <<"  // one. (Note that, analogously to `newp', the variables `newc', \n"
         <<"  // `news', `newv' and `newreg' select the highest number amongst \n"
         <<"  // currently defined curves, surfaces, volumes and `any entities \n"
         <<"  // other than points', respectively.) \n"
         <<" \n"
         <<"  p1 = newp; Point(p1) = {x,  y,  z,  lcar} ; \n"
         <<"  p2 = newp; Point(p2) = {x+r,y,  z,  lcar} ; \n"
         <<"  p3 = newp; Point(p3) = {x,  y+r,z,  lcar} ; \n"
         <<"  p4 = newp; Point(p4) = {x,  y,  z+r,lcar} ; \n"
         <<"  p5 = newp; Point(p5) = {x-r,y,  z,  lcar} ; \n"
         <<"  p6 = newp; Point(p6) = {x,  y-r,z,  lcar} ; \n"
         <<"  p7 = newp; Point(p7) = {x,  y,  z-r,lcar} ; \n"
         <<" \n"
         <<"  c1 = newreg; Circle(c1) = {p2,p1,p7}; \n"
         <<"  c2 = newreg; Circle(c2) = {p7,p1,p5}; \n"
         <<"  c3 = newreg; Circle(c3) = {p5,p1,p4}; \n"
         <<"  c4 = newreg; Circle(c4) = {p4,p1,p2}; \n"
         <<"  c5 = newreg; Circle(c5) = {p2,p1,p3}; \n"
         <<"  c6 = newreg; Circle(c6) = {p3,p1,p5}; \n"
         <<"  c7 = newreg; Circle(c7) = {p5,p1,p6}; \n"
         <<"  c8 = newreg; Circle(c8) = {p6,p1,p2}; \n"
         <<"  c9 = newreg; Circle(c9) = {p7,p1,p3}; \n"
         <<"  c10 = newreg; Circle(c10) = {p3,p1,p4}; \n"
         <<"  c11 = newreg; Circle(c11) = {p4,p1,p6}; \n"
         <<"  c12 = newreg; Circle(c12) = {p6,p1,p7}; \n"
         <<" \n"
         <<"  // We need non-plane surfaces to define the spherical holes. Here we \n"
         <<"  // use ruled surfaces, which can have 3 or 4 sides: \n"
         <<" \n"
         <<"  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1}; \n"
         <<"  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2}; \n"
         <<"  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3}; \n"
         <<"  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4}; \n"
         <<"  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5}; \n"
         <<"  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6}; \n"
         <<"  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7}; \n"
         <<"  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8}; \n"
         <<" \n"
         <<"  // We then store the surface loops identification numbers in list \n"
         <<"  // for later reference (we will need these to define the final \n"
         <<"  // volume): \n"
         <<" \n"
         <<"  theloops[t] = newreg ; \n"
         <<" \n"
         <<"  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; \n"
         <<" \n"
         <<"  thesphere = newreg ; \n"
         <<"  Volume(thesphere) = theloops[t] ; \n"
         <<" \n"
         <<"Return \n";


    ostr << " x = 0 ; y = 0 ; z = 0 ; "
         << "lcar = " << this->h() * M_farfactor << ";\n"
         << "t=2; r = 1.8*" << M_farBnd << " ;\n"
         << "Call mySphere ;\n"
         << "// We define a physical volume :\n"
         << "//Physical Volume (4) = thesphere ;\n"
         << "\n"
         << "lcar = " << this->h() << ";\n"
         << "t=1; r = 1.8" << " ;\n"
         << "Call mySphere ;\n"
         << "// We define a physical volume :\n"
         << "//Physical Volume (3) = thesphere ;\n"
         << "\n"
         << " // Point Characteristic Description\n"
         //<< getPointCharacteristicDescription()
         << "\n"
         << "v1=newv; Volume(v1) = {theloops[1]} ;\n"
         << "v2=newv; Volume(v2) = {theloops[2],theloops[1]} ;\n"
         << "\n"
         << "Physical Volume (1) = {v1,v2};\n";

    return ostr.str();

}


std::string
myMesh::getPointCharacteristicDescription( ) const
{
    int i;
    std::ostringstream ostr;
    std::list<std::vector<double> >::const_iterator it( M_ptChar.begin() );

    if ( it == M_ptChar.end() ) return ostr.str();

    ostr << "fact = " << M_factor << ";\n";
    ostr << "p1 = newp; Point(p1) = { ";

    for ( i=0; i< nDim; i++ )
        ostr << ( *it )[i] << ", ";

    ostr << " fact};\n";

    for ( it++ ; it != M_ptChar.end() ; it++ )
    {
        ostr << "p = newp; Point(p) = { ";

        for ( i=0; i< nDim; i++ )
            ostr << ( *it )[i] << ", ";

        ostr << " fact};\n";
    }

    ostr << "Attractor Point {p1:p} = { h/2, h*fact,h};\n";

    return ostr.str();
}

} // namespace Feel
