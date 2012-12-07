/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-10-11

  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
   \file gmshsimplexdomain.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-10-11
 */
#ifndef __GMSHSIMPLEXDOMAIN_CPP
#define __GMSHSIMPLEXDOMAIN_CPP

#include <feel/feelcore/feel.hpp>

#include <feel/feelfilters/gmshsimplexdomain.hpp>


namespace Feel
{

GmshSimplexDomain::GmshSimplexDomain( int dim, int order, DomainType dt )
    :
    super( dim,order ),
    M_descr()
{
    switch ( dt )
    {
    case GMSH_REAL_DOMAIN:
    {

        if ( this->dimension() >= 1 )
            this->M_I[0] = std::make_pair( 0, 1 );

        if ( this->dimension() >= 2 )
            this->M_I[1] = std::make_pair( 0, 1 );

        if ( this->dimension() >= 3 )
            this->M_I[2] = std::make_pair( 0, 1 );
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
GmshSimplexDomain::getDescription() const
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
GmshSimplexDomain::getDescription1D() const
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
                    << "Physical Line(\"Mat1\") = {1};\n"
                    << "Physical Line(\"Mat2\") = {2};\n";
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
                 << "Physical Line(\"Mat1\") = {1};\n";
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
GmshSimplexDomain::getDescription2D() const
{
    std::ostringstream ostr;
    ostr << this->preamble() << "\n";
    ostr << "Point(1) = {" << this->M_I[0].first << "," << this->M_I[1].first << ",0,h};\n"
         << "Point(2) = {" << this->M_I[0].second << "," << this->M_I[1].first << ",0,h};\n"
         << "Point(3) = {" << this->M_I[0].first << "," << this->M_I[1].second << ",0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,1};\n";

    //if ( this->isReference() )
    if ( this->h() >= ( this->M_I[0].second-this->M_I[0].first ) )
        ostr <<"Transfinite Line{1} = 1;\n"
             <<"Transfinite Line{2} = 1;\n"
             <<"Transfinite Line{3} = 1;\n";

    ostr << "Line Loop(4) = {3,1,2};\n";

    if ( this->usePhysicalNames() == false )
    {
        ostr << "Plane Surface(5) = {4};\n"
             << "Physical Line(6) = {1};\n"
             << "Physical Line(7) = {2};\n"
             << "Physical Line(8) = {3};\n"
             << "Physical Surface(9) = {5};\n";
    }

    else
    {
        ostr << "Plane Surface(5) = {4};\n"
             << "Physical Line(\"Dirichlet\") = {1};\n"
             << "Physical Line(\"Neumann\") = {2,3};\n"
             << "Physical Surface(\"Mat1\") = {5};\n";
    }

    return ostr.str();
}
// 3D
std::string
GmshSimplexDomain::getDescription3D() const
{
    std::ostringstream ostr;
    ostr << this->preamble() << "\n";
    ostr << "Point(1) = {" << this->M_I[0].first << "," << this->M_I[1].first << "," << this->M_I[2].first << ",h};\n"
         << "Point(2) = {" << this->M_I[0].second << "," << this->M_I[1].first << "," << this->M_I[2].first << ",h};\n"
         << "Point(3) = {" << this->M_I[0].first << "," << this->M_I[1].second << "," << this->M_I[2].first << ",h};\n"
         << "Point(4) = {" << this->M_I[0].first << "," << this->M_I[1].first << "," << this->M_I[2].second << ",h};\n"
         << "Line(1) = {1,2};" << "\n"
         << "Line(2) = {2,3};" << "\n"
         << "Line(3) = {3,1};" << "\n"
         << "Line(4) = {1,4};" << "\n"
         << "Line(5) = {2,4};" << "\n"
         << "Line(6) = {3,4};" << "\n"
         << "Line Loop(4) = {3,1,2};" << "\n"
         << "Plane Surface(5) = {4};" << "\n"
         << "Line Loop(10) = {6, -4, -3};" << "\n"
         << "Plane Surface(11) = {10};" << "\n"
         << "Line Loop(12) = {6, -5, 2};" << "\n"
         << "Plane Surface(13) = {12};" << "\n"
         << "Line Loop(14) = {4, -5, -1};" << "\n"
         << "Plane Surface(15) = {14};" << "\n";

    //if ( this->isReference() )
    if ( this->h() >= ( this->M_I[0].second-this->M_I[0].first ) )
    {
        ostr << "Transfinite Line(1) = 1;\n"
             << "Transfinite Line(2) = 1;\n"
             << "Transfinite Line(3) = 1;\n"
             << "Transfinite Line(4) = 1;\n"
             << "Transfinite Line(5) = 1;\n"
             << "Transfinite Line(6) = 1;\n";

        ostr << "Transfinite Surface(5)=1;\n"
             << "Transfinite Surface(11)=1;\n"
             << "Transfinite Surface(13)=1;\n"
             << "Transfinite Surface(15)=1;\n";
    }

    ostr << "Surface Loop(20) = {11, 13, 15, 5};" << "\n"
         << "Volume(21) = {20};" << "\n"
         << "" << "\n";

    if ( this->usePhysicalNames() == false )
    {
        ostr << "Physical Surface(16) = {11};" << "\n"
             << "Physical Surface(17) = {15};" << "\n"
             << "Physical Surface(18) = {5};" << "\n"
             << "Physical Surface(19) = {13};" << "\n"
             << "Physical Volume(22) = {21};" << "\n";
    }

    else
    {
        ostr << "Physical Surface(\"Neumann\") = {11,15,13};" << "\n"
             << "Physical Surface(\"Dirichlet\") = {5};" << "\n"
             << "Physical Volume(\"Mat1\") = {21};" << "\n";
    }

    return ostr.str();
}


} // Feel
#endif // __GMSHSIMPLEXDOMAIN_CPP
