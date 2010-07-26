/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-12-29

  Copyright (C) 2006 Université Joseph Fourier (Grenoble)

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
   \file gmshtensorizeddomain.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-12-29
 */
#ifndef __GMSHTENSORIZEDDOMAIN_HPP
#define __GMSHTENSORIZEDDOMAIN_HPP 1

#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifemesh/geoentity.hpp>

namespace Life
{
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
const uint16_type GmshTensorizedDomain<Dim, Order, RDim, Entity>::nDim;
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
const uint16_type GmshTensorizedDomain<Dim, Order, RDim, Entity>::nOrder;
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
const uint16_type GmshTensorizedDomain<Dim, Order, RDim, Entity>::nRealDim;


template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
std::string
GmshTensorizedDomain<Dim, Order, RDim, Entity>::getDescription( mpl::int_<1>,  mpl::bool_<false> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {" << _M_I[0].first << ",";
    if ( nRealDim == nDim + 1 )
        ostr << _M_I[1].first;
    else
        ostr << 0;
    ostr << ",0,h};\n"
         << "Point(2) = {" << _M_I[0].second << ",";
    if ( nRealDim == nDim + 1 )
        ostr << _M_I[1].second;
    else
        ostr << 0;
    ostr << ",0,h};\n";
    if ( this->addMidPoint() )
    {
        ostr << "Point(3) = {" << (_M_I[0].second+_M_I[0].first)/2 << ",";
        if ( nRealDim == nDim + 1 )
            ostr << (_M_I[1].second+_M_I[1].first)/2;
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
                  << "Physical Line(1) = {1};\n"
                  << "Physical Line(2) = {2};\n";
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
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
std::string
GmshTensorizedDomain<Dim, Order, RDim, Entity>::getDescription( mpl::int_<2>,  mpl::bool_<false> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "a=" << _M_I[0].first << ";\n"
         << "b=" << _M_I[0].second << ";\n"
         << "c=" << _M_I[1].first << ";\n"
         << "d=" << _M_I[1].second << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {a,c,0.0,h};\n"
         << "Point(2) = {b,c,0.0,h};\n"
         << "Point(3) = {b,d,0.0,h};\n"
         << "Point(4) = {a,d,0.0,h};\n"
         << "Line(1) = {4,1};\n"
         << "Line(2) = {1,2};\n"
         << "Line(3) = {2,3};\n"
         << "Line(4) = {3,4};\n"
         << "Line Loop(5) = {1,2,3,4};\n"
         << "Plane Surface(6) = {5};\n";
    if ( this->usePhysicalNames() == false )
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
    return ostr.str();
}
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
std::string
GmshTensorizedDomain<Dim, Order, RDim, Entity>::getDescription( mpl::int_<2>,  mpl::bool_<true> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << getDescription( mpl::int_<2>(), mpl::bool_<false>() )
         << "nx = 1/h;\n"
         << "ny = 1/h;\n"
         << "\n"
         << "Transfinite Line {1,3} = ny + 1 Using Progression 1.0;\n"
         << "Transfinite Line {2,4} = nx + 1 Using Progression 1.0;\n"
         << "\n"
         << "Transfinite Surface {6} = {1,2,3,4};\n"
         << "Recombine Surface {6};\n";
    return ostr.str();
}
// 3D
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
std::string
GmshTensorizedDomain<Dim, Order, RDim, Entity>::getDescription( mpl::int_<3>,  mpl::bool_<false>, bool do_recombine ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "a=" << _M_I[0].first << ";\n"
         << "b=" << _M_I[0].second << ";\n"
         << "c=" << _M_I[1].first << ";\n"
         << "d=" << _M_I[1].second << ";\n"
         << "e=" << _M_I[2].first << ";\n"
         << "f=" << _M_I[2].second << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {a,c,e,h};\n"
         << "Point(2) = {b,c,e,h};\n"
         << "Point(3) = {b,d,e,h};\n"
         << "Point(4) = {a,d,e,h};\n"
         << "Line(1) = {1,4};\n"
         << "Line(2) = {4,3};\n"
         << "Line(3) = {3,2};\n"
         << "Line(4) = {2,1};\n"
         << "Line Loop(5) = {3,4,1,2};\n"
         << "Plane Surface(6) = {5};\n"
         << "\n"
         << "Extrude Surface {6, {0,0,f-e} } {\n"
         << "  Layers { {(f-e)/h}, {1.0} };\n";
    if ( do_recombine )
        ostr << "  Recombine;\n";

    ostr << "};\n"
         << "Physical Line(1) = {1};\n"
         << "Physical Line(2) = {2};\n"
         << "Physical Line(3) = {3};\n"
         << "Physical Line(4) = {4};\n";
    if ( this->usePhysicalNames() == false )
    {
        ostr << "Physical Surface(6) = {6};\n"
             << "Physical Surface(15) = {15};\n"
             << "Physical Surface(19) = {19};\n"
             << "Physical Surface(23) = {23};\n"
             << "Physical Surface(27) = {27};\n"
             << "Physical Surface(28) = {28};\n"
             << "Physical Volume(30) = {1};\n";
    }
    else
    {
        ostr << "Physical Surface(\"Neumann\") = {6,19,27,28};\n"
             << "Physical Surface(\"Dirichlet\") = {15,23};\n"
             << "Physical Volume(\"Mat1\") = {1};\n";
    }
    return ostr.str();
}
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
std::string
GmshTensorizedDomain<Dim, Order, RDim, Entity>::getDescription( mpl::int_<3>,  mpl::bool_<true> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << getDescription( mpl::int_<3>(), mpl::bool_<false>(), true )
         << "nx = 1/h;\n"
         << "ny = 1/h;\n"
         << "nz = 1/h;\n"
         << "\n"
         << "Transfinite Line {4,10,2,8} = nx Using Progression 1;\n"
         << "Transfinite Line {3,9,1,11} = ny Using Progression 1;\n"
         << "Transfinite Line {14,18,13,22} = nz Using Progression 1;\n"
         << "\n"
         << "//Transfinite Surface {23} = {1,10,14,2};\n"
         << "//Transfinite Surface {19} = {1,10,6,4};\n"
         << "//Transfinite Surface {15} = {4,6,5,3};\n"
         << "Transfinite Surface {6} = {2,1,4,3};\n"
         << "//Transfinite Surface {27} = {14,2,3,5};\n"
         << "//Transfinite Surface {28} = {6,10,14,5};\n"
         << "Recombine Surface {27,23,6,19,15,28};\n";
    return ostr.str();

}

#if defined( LIFE_INSTANTIATION_MODE )

// Instantiations
// 1D
template class GmshTensorizedDomain<1,1,1,Simplex>;
template class GmshTensorizedDomain<1,1,2,Simplex>;
template class GmshTensorizedDomain<1,1,1,SimplexProduct>;
template class GmshTensorizedDomain<1,2,1,Simplex>;
template class GmshTensorizedDomain<1,2,1,SimplexProduct>;
// 2D
template class GmshTensorizedDomain<2,1,2,Simplex>;
template class GmshTensorizedDomain<2,1,2,SimplexProduct>;
template class GmshTensorizedDomain<2,2,2,Simplex>;
template class GmshTensorizedDomain<2,2,2,SimplexProduct>;
// 2D3D
template class GmshTensorizedDomain<2,1,3,Simplex>;
// 3D
template class GmshTensorizedDomain<3,1,3,Simplex>;
template class GmshTensorizedDomain<3,1,3,SimplexProduct>;
template class GmshTensorizedDomain<3,2,3,Simplex>;
template class GmshTensorizedDomain<3,2,3,SimplexProduct>;
#endif // LIFE_INSTANTIATION_MODE

}

#endif // __GMSHTENSORIZEDDOMAIN_HPP
