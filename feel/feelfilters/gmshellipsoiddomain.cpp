/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-08-02

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file gmshellipsoiddomain.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-08-02
 */
#include <feel/feelcore/feel.hpp>

#include <feel/feelfilters/gmshellipsoiddomain.hpp>

namespace Feel
{
GmshEllipsoidDomain::GmshEllipsoidDomain( int Dim, int Order, DomainType dt )
    :
    super( Dim, Order )
{
    switch ( dt )
    {
    case GMSH_REAL_DOMAIN:
    {

        if ( this->dimension() >= 1 )
            this->M_I[0] = std::make_pair( -1, 1 );

        if ( this->dimension() >= 2 )
            this->M_I[1] = std::make_pair( -1, 1 );

        if ( this->dimension() >= 3 )
            this->M_I[2] = std::make_pair( -1, 1 );
    }
    break;

    case GMSH_REFERENCE_DOMAIN:
    {
        if ( this->dimension() >= 1 )
            this->M_I[0] = std::make_pair( -1, 1 );

        if ( this->dimension() >= 2 )
            this->M_I[1] = std::make_pair( -1, 1 );

        if ( this->dimension() >= 3 )
            this->M_I[2] = std::make_pair( -1, 1 );
    }
    break;
    }
}

std::string
GmshEllipsoidDomain::getDescription() const
{
    switch ( this->dimension() )
    {
    case 1:
        return this->getDescription1D();

    case 2:
        return this->getDescription2D();

    case 3:
        return this->getDescription3D();

    default:
        return std::string();
    }
}
std::string
GmshEllipsoidDomain::getDescription1D() const
{
    std::ostringstream ostr;
    ostr << this->preamble() << "\n";

    ostr << "Point(1) = {" << this->M_I[0].first << ",0,0,h};\n"
         << "Point(2) = {" << this->M_I[0].second << ",0,0,h};\n";

    if ( this->addMidPoint() )
    {
        ostr << "Point(3) = {" << ( this->M_I[0].second+this->M_I[0].first )/2 << ",0,0,h};\n"
             << "Line(1) = {1,3};\n"
             << "Line(2) = {3,2};\n";

        if ( this->usePhysicalNames() == false )
        {
            ostr    << "Physical Point(1) = {1};\n"
                    << "Physical Point(3) = {2};\n"
                    << "Physical Point(2) = {3};\n"
                    << "Physical Line(1) = {1};\n"
                    << "Physical Line(2) = {2};\n";
        }

        else
        {
            ostr    << "Physical Point(\"Dirichlet\") = {1};\n"
                    << "Physical Point(\"Neumann\") = {2};\n"
                    << "Physical Point(3) = {3};\n"
                    << "Physical Line(\"Mat1\") = {1};\n"
                    << "Physical Line(\"Mat2\") = {2};\n";
        }
    }

    else
    {
        if ( this->usePhysicalNames() == false )
        {
            ostr << "Line(1) = {1,2};\n"
                 << "Physical Point(1) = {1};\n"
                 << "Physical Point(3) = {2};\n"
                 << "Physical Line(1) = {1};\n";
        }

        else
        {
            ostr << "Line(1) = {1,2};\n"
                 << "Physical Point(\"Dirichlet\") = {1};\n"
                 << "Physical Point(\"Neumann\") = {2};\n"
                 << "Physical Line(\"Mat1\") = {1};\n";
        }
    }

    return ostr.str();
}
// 2D
std::string
GmshEllipsoidDomain::getDescription2D() const
{
    std::ostringstream ostr;
    ostr << this->preamble() << "\n";

    double xcenter = ( this->xmax()+this->xmin() )/2;
    double ycenter = ( this->ymax()+this->ymin() )/2;
    double zcenter = 0;
    ostr << "Point(1) = {" << xcenter << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(2) = {" << this->xmax() << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(3) = {" << xcenter << "," << this->ymax() << "," << zcenter << ",h};\n"
         << "Circle(1) = {2,1,3};\n"
         << "Point(4) = {" << this->xmin() << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(5) = {" << xcenter << "," << this->ymin() << "," << zcenter << ",h};\n"
         << "Circle(2) = {3,1,4};\n"
         << "Circle(3) = {4,1,5};\n"
         << "Circle(4) = {5,1,2};\n"
         << "Line Loop(5) = {2, 3, 4, 1};\n"
         << "Ruled Surface(6) = {5};\n";

    // add center point in the mesh
    ostr << "Point{ 1 } In Surface{ 6 };\n";
    ostr << "Physical Point(\"center\") = {1};\n";

    if ( this->usePhysicalNames() == false )
    {
        ostr << "Physical Line(7) = {1, 4, 3, 2};\n"
             << "Physical Surface(8) = {6};\n";
    }

    else
    {
        ostr << "Physical Line(\"Neumann\") = {1, 4};\n"
             << "Physical Line(\"Dirichlet\") = {3, 2};\n"
             << "Physical Surface(\"Mat_1\") = {6};\n";
    }



    return ostr.str();
}
// 3D
std::string
GmshEllipsoidDomain::getDescription3D() const
{
    std::ostringstream ostr;
    ostr << this->preamble() << "\n";

    double xcenter = ( this->xmax()+this->xmin() )/2;
    double ycenter = ( this->ymax()+this->ymin() )/2;
    double zcenter = ( this->zmax()+this->zmin() )/2;
    ostr << "Point(1) = {" << xcenter << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(2) = {" << this->xmax() << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(3) = {" << xcenter << "," << this->ymax() << "," << zcenter << ",h};\n"
         << "Circle(1) = {2,1,3};\n"
         << "Point(4) = {" << this->xmin() << "," << ycenter << "," << zcenter << ",h};\n"
         << "Point(5) = {" << xcenter << "," << this->ymin() << "," << zcenter << ",h};\n"
         << "Circle(2) = {3,1,4};\n"
         << "Circle(3) = {4,1,5};\n"
         << "Circle(4) = {5,1,2};\n"
         << "Point(6) = {" << xcenter << "," << ycenter << "," << this->zmin() << ",h};\n"
         << "Point(7) = {" << xcenter << "," << ycenter << "," << this->zmax() << ",h};\n"
         << "Circle(5) = {3,1,6};\n"
         << "Circle(6) = {6,1,5};\n"
         << "Circle(7) = {5,1,7};\n"
         << "Circle(8) = {7,1,3};\n"
         << "Circle(9) = {2,1,7};\n"
         << "Circle(10) = {7,1,4};\n"
         << "Circle(11) = {4,1,6};\n"
         << "Circle(12) = {6,1,2};\n"
         << "Line Loop(13) = {2,8,-10};\n"
         << "Ruled Surface(14) = {13};\n"
         << "Line Loop(15) = {10,3,7};\n"
         << "Ruled Surface(16) = {15};\n"
         << "Line Loop(17) = {-8,-9,1};\n"
         << "Ruled Surface(18) = {17};\n"
         << "Line Loop(19) = {-11,-2,5};\n"
         << "Ruled Surface(20) = {19};\n"
         << "Line Loop(21) = {-5,-12,-1};\n"
         << "Ruled Surface(22) = {21};\n"
         << "Line Loop(23) = {-3,11,6};\n"
         << "Ruled Surface(24) = {23};\n"
         << "Line Loop(25) = {-7,4,9};\n"
         << "Ruled Surface(26) = {25};\n"
         << "Line Loop(27) = {-4,12,-6};\n"
         << "Ruled Surface(28) = {27};\n"
         << "Surface Loop(29) = {28,26,16,14,20,24,22,18};\n"
         << "Volume(30) = {29};\n";

    if ( this->usePhysicalNames() == false )
    {
        ostr << "Physical Surface(1) = {28,26,16,14,20,24,22,18};\n"
             << "Physical Volume(2) = {30};\n";
    }

    else
    {
        ostr << "Physical Surface(\"Neumann\") = {28,26,16,14};\n"
             << "Physical Surface(\"Dirichlet\") = {20,24,22,18};\n"
             << "Physical Volume(\"Mat_1\") = {30};\n";
    }

    return ostr.str();
}


}
