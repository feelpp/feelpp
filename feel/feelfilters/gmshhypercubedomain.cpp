/* -*- Mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-12-29

  Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble)

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
   \file gmshhypercubedomain.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-29
 */
#ifndef __GMSHTENSORIZEDDOMAIN_HPP
#define __GMSHTENSORIZEDDOMAIN_HPP 1

#include <feel/feelfilters/gmshhypercubedomain.hpp>

namespace Feel
{
GmshHypercubeDomain::GmshHypercubeDomain( int dim, int order  )
    :
    super( dim, order ),
    M_rdim( dim ),
    M_use_hypercube( false )
{}

GmshHypercubeDomain::GmshHypercubeDomain( int dim, int order, int rdim, bool use_hypercube  )
    :
    super( dim, order ),
    M_rdim( rdim ),
    M_use_hypercube( use_hypercube )
{}

GmshHypercubeDomain::GmshHypercubeDomain( GmshHypercubeDomain const& d )
    :
    super( d ),
    M_rdim( d.M_rdim ),
    M_use_hypercube( d.M_use_hypercube )
{}
GmshHypercubeDomain::~GmshHypercubeDomain()
{}
std::string
GmshHypercubeDomain::getDescription() const
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
GmshHypercubeDomain::getDescription1D() const
{
    std::ostringstream ostr;
    ostr << this->preamble();
    ostr << "Point(1) = {" << this->M_I[0].first << ",";

    if ( M_rdim == this->dimension() + 1 )
        ostr << this->M_I[1].first;

    else
        ostr << 0;

    ostr << ",0,h};\n"
         << "Point(2) = {" << this->M_I[0].second << ",";

    if ( M_rdim == this->dimension() + 1 )
        ostr << this->M_I[1].second;

    else
        ostr << 0;

    ostr << ",0,h};\n";

    if ( this->addMidPoint() )
    {
        ostr << "Point(3) = {" << ( this->M_I[0].second+this->M_I[0].first )/2 << ",";

        if ( M_rdim == this->dimension() + 1 )
            ostr << ( this->M_I[1].second+this->M_I[1].first )/2;

        else
            ostr << 0;

        ostr << ",0,h};\n"
             << "Line(1) = {1,3};\n"
             << "Line(2) = {3,2};\n";

        if ( this->usePhysicalNames() == false )
        {
            ostr  << "Physical Point(1) = {1};\n"
                  << "Physical Point(3) = {2};\n"
                  << "Physical Point(2) = {3};\n"
                  << "Physical Line(\"Mat1\") = {1};\n"
                  << "Physical Line(\"Mat2\") = {2};\n";
        }

        else
        {
            ostr  << "Physical Point(\"Dirichlet\") = {1};\n"
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
GmshHypercubeDomain::getDescription2D() const
{
    std::ostringstream ostr;
    ostr << this->preamble();


    ostr << "xmin=" << this->M_I[0].first << ";\n"
         << "xmax=" << this->M_I[0].second << ";\n"
         << "ymin=" << this->M_I[1].first << ";\n"
         << "ymax=" << this->M_I[1].second << ";\n";

    ostr << "Point(1) = {xmin,ymin,0.0,h};\n"
         << "Point(2) = {xmax,ymin+"<<this->shear()<<",0.0,h};\n"
         << "Point(3) = {xmax+" << this->shear() << ",ymax,0.0,h};\n"
         << "Point(4) = {xmin+" << this->shear() << ",ymax+"<<this->shear()<<",0.0,h};\n"
         << "Line(1) = {4,1};\n"
         << "Line(2) = {1,2};\n"
         << "Line(3) = {2,3};\n"
         << "Line(4) = {3,4};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Plane Surface(6) = {5};\n";

    if ( !this->periodic().empty() )
    {
        std::for_each( this->periodic().begin(),
                      this->periodic().end(),
                       [&ostr]( std::pair<int,std::pair<int,int> > const& p )
                       {
                           CHECK( p.first == 1 ) << "Invalid periodic entity dimension : " << p.first << ",  should be 1";
                           ostr << "Periodic Line ( " << p.second.first << " ) = { " << p.second.second << " };\n";
                      } );
    }

    if ( this->subStructuring() == true )
    {
        ostr << "Physical Point(\"CrossPoints\") = {1,2,3,4};\n";
        ostr << "Physical Line(\"NORTH\") = {4};\n"
             << "Physical Line(\"WEST\") = {1};\n"
             << "Physical Line(\"SOUTH\") = {2};\n"
             << "Physical Line(\"EAST\") = {3};\n"
             << "Physical Surface(\"Omega\") = {6};\n";

    }
    else if ( this->usePhysicalNames() == false )
    {
        ostr << "Physical Line(1) = {1};\n"
             << "Physical Line(2) = {2};\n"
             << "Physical Line(3) = {3};\n"
             << "Physical Line(4) = {4};\n"
             << "Physical Surface(6) = {6};\n";
    }
    else
    {
        ostr << "Physical Line(\"Dirichlet\") = {1,3};\n"
             << "Physical Line(\"Neumann\") = {2,4};\n"
             << "Physical Surface(\"Mat1\") = {6};\n";
    }

    if ( this->structuredMesh() || M_use_hypercube )
    {
        //if ( this->structuredMesh() == 1 || M_use_hypercube )
        if ( this->structuredMesh() == 1 )
            ostr << "nx = (xmax-xmin)/h;\n"
                 << "ny = (ymax-ymin)/h;\n";
        else if ( this->structuredMesh() == 2 )
            ostr << "nx = 1/h;\n"
                 << "ny = 1/h;\n";
        else if ( this->structuredMesh() == 3 )
            ostr << "nx=" << M_nx << ";\n"
                 << "ny=" << M_ny << ";\n";


        ostr << "\n"
             << "Transfinite Line {1,3} = ny + 1 Using Progression 1.0;\n"
             << "Transfinite Line {2,4} = nx + 1 Using Progression 1.0;\n"
             << "\n"
             << "Transfinite Surface {6} = {1,2,3,4};\n";

    }
    if ( M_use_hypercube )
    {
        ostr << "Recombine Surface {6};\n";
    }

    return ostr.str();
}
// 3D
std::string
GmshHypercubeDomain::getDescription3D() const
{
    std::ostringstream ostr;
    ostr << this->preamble();

    ostr << "xmin=" << this->M_I[0].first << ";\n"
         << "xmax=" << this->M_I[0].second << ";\n"
         << "ymin=" << this->M_I[1].first << ";\n"
         << "ymax=" << this->M_I[1].second << ";\n"
         << "zmin=" << this->M_I[2].first << ";\n"
         << "zmax=" << this->M_I[2].second << ";\n"
         << "Point(1) = {xmin,ymin,zmin,h};\n"
         << "Point(2) = {xmax,ymin,zmin,h};\n"
         << "Point(3) = {xmax,ymax,zmin,h};\n"
         << "Point(4) = {xmin,ymax,zmin,h};\n"
         << "Line(1) = {2,3};\n"
         << "Line(2) = {3,4};\n"
         << "Line(3) = {4,1};\n"
         << "Line(4) = {1,2};\n"
         << "Line Loop(5) = {2,3,4,1};\n"
         << "Plane Surface(6) = {5};\n"
         << "\n"
         << "Extrude Surface {6, {0,0,zmax-zmin} } {\n"
         << "  Layers { {(zmax-zmin)/h}, {1.0} };\n";



    if ( M_use_hypercube )
        ostr << "  Recombine;\n";

    ostr << "};\n";

    // if one wants a mesh which is regular using transfinite, remove
    // the 'if' however tre is a bug in gmsh which prevent of doing
    // thatat momentsee // https://github.com/feelpp/feelpp/issues/147
    // for more details

    if ( this->structuredMesh() || M_use_hypercube )
    {
        if ( this->structuredMesh() == 1 )
            ostr << "nx = (xmax-xmin)/h;\n"
                 << "ny = (ymax-ymin)/h;\n"
                 << "nz = (zmax-zmin)/h;\n";
        else if ( this->structuredMesh() == 2 )
            ostr << "nx = 1/h;\n"
                 << "ny = 1/h;\n"
                 << "nz = 1/h;\n";
        else if ( this->structuredMesh() == 3 )
            ostr << "nx=" << M_nx << ";\n"
                 << "ny=" << M_ny << ";\n"
                 << "nz=" << M_nz << ";\n";
    }


    if ( M_use_hypercube )
        ostr << "Transfinite Line {4,10,2,8} = nx + 1  Using Progression 1;\n"
             << "Transfinite Line {3,9,1,11} = ny + 1 Using Progression 1;\n"
             << "Transfinite Line {14,18,13,22} = nz + 1 Using Progression 1;\n"
             << "\n"
             << "Transfinite Surface {6} = {3,2,1,4};\n"
             << "Transfinite Surface {23} = {14,2,1,10};\n"
             << "Transfinite Surface {19} = {6,10,1,4};\n"
             << "Transfinite Surface {15} = {5,3,4,6};\n"
             << "Transfinite Surface {27} = {5,14,2,3};\n"
             << "Transfinite Surface {28} = {6,10,14,5};\n";
    if ( this->usePhysicalNames() == false && this->subStructuring() == false )
    {
        ostr << "Physical Line(1) = {1};\n"
             << "Physical Line(2) = {2};\n"
             << "Physical Line(3) = {3};\n"
             << "Physical Line(4) = {4};\n";

        ostr << "Physical Surface(6) = {6};\n"
             << "Physical Surface(15) = {15};\n"
             << "Physical Surface(19) = {19};\n"
             << "Physical Surface(23) = {23};\n"
             << "Physical Surface(27) = {27};\n"
             << "Physical Surface(28) = {28};\n"
             << "Physical Volume(30) = {1};\n";
    }
    else if ( this->subStructuring() == true )
    {
        ostr << "Physical Point(\"CrossPoints\") = {1,2,3,4,5,6,10,14};\n";
        ostr << "Physical Line(\"WireBasket\") = {1,2,3,4,8,9,10,11,13,14,18,22};\n";
        ostr << "Physical Surface(\"TOP\") = {6};\n"
             << "Physical Surface(\"NORTH\") = {15};\n"
             << "Physical Surface(\"WEST\") = {19};\n"
             << "Physical Surface(\"SOUTH\") = {23};\n"
             << "Physical Surface(\"EAST\") = {27};\n"
             << "Physical Surface(\"BOTTOM\") = {28};\n"
             << "Physical Volume(\"Omega\") = {1};\n";
            //<< "Physical Volume(30) = {1};\n";
    }
    else
    {
        ostr << "Physical Surface(\"Neumann\") = {6,19,27,28};\n"
             << "Physical Surface(\"Dirichlet\") = {15,23};\n"
             << "Physical Volume(\"Mat1\") = {1};\n";
    }

    if ( M_use_hypercube )
    {
        ostr << "\n"
             << "Transfinite Line {4,10,2,8} = nx + 1  Using Progression 1;\n"
             << "Transfinite Line {3,9,1,11} = ny + 1 Using Progression 1;\n"
             << "Transfinite Line {14,18,13,22} = nz + 1 Using Progression 1;\n"
             << "\n"
             << "//Transfinite Surface {23} = {14,2,1,10};\n"
             << "//Transfinite Surface {19} = {6,10,1,4};\n"
             << "//Transfinite Surface {15} = {5,3,4,6};\n"
             << "Transfinite Surface {6} = {3,2,1,4};\n"
             << "//Transfinite Surface {27} = {5,14,2,3};\n"
             << "//Transfinite Surface {28} = {6,10,14,5};\n"
             << "Recombine Surface {27,23,6,19,15,28};\n";
    }

    return ostr.str();
}


} // Feel

#endif // __GMSHTENSORIZEDDOMAIN_HPP
