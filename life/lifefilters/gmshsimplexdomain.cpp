/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-11
 */
#ifndef __GMSHSIMPLEXDOMAIN_CPP
#define __GMSHSIMPLEXDOMAIN_CPP

#include <life/lifefilters/gmshsimplexdomain.hpp>
#include <life/lifemesh/geoentity.hpp>

namespace Life
{
template<int Dim, int Order>
const uint16_type GmshSimplexDomain<Dim, Order>::nDim;
template<int Dim, int Order>
const uint16_type GmshSimplexDomain<Dim, Order>::nOrder;

template<int Dim, int Order>
std::string
GmshSimplexDomain<Dim, Order>::getDescription( mpl::int_<1> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {" << _M_I[0].first << ",0,0,h};\n"
         << "Point(2) = {" << _M_I[0].second << ",0,0,h};\n"
         << "Point(3) = {" << (_M_I[0].second+_M_I[0].first)/2 << ",0,0,h};\n"
         << "Line(1) = {1,3};\n"
         << "Line(2) = {3,2};\n"
         << "Physical Point(1) = {1};\n"
         << "Physical Point(2) = {2};\n"
         << "Physical Point(3) = {3};\n"
         << "Physical Line(1) = {1};\n"
         << "Physical Line(2) = {2};\n";
    return ostr.str();
}
// 2D
template<int Dim, int Order>
std::string
GmshSimplexDomain<Dim, Order>::getDescription( mpl::int_<2> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {-1,-1,0,h};\n"
         << "Point(2) = {1,-1,0,h};\n"
         << "Point(3) = {-1,1,0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,1};\n"
         << "Line Loop(4) = {3,1,2};\n"
         << "Plane Surface(5) = {4};\n"
         << "Physical Line(6) = {1};\n"
         << "Physical Line(7) = {2};\n"
         << "Physical Line(8) = {3};\n"
         << "Physical Surface(9) = {5};\n";
    return ostr.str();
}
// 3D
template<int Dim, int Order>
std::string
GmshSimplexDomain<Dim, Order>::getDescription( mpl::int_<3> ) const
{
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
         << "h=" << _M_h << ";\n"
         << "Point(1) = {-1,-1,-1,h};" << "\n"
         << "Point(2) = {1,-1,-1,h};" << "\n"
         << "Point(3) = {-1,1,-1,h};" << "\n"
         << "Point(4) = {-1,-1,1,h};" << "\n"
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
         << "Plane Surface(15) = {14};" << "\n"
         << "Surface Loop(20) = {11, 13, 15, 5};" << "\n"
         << "Volume(21) = {20};" << "\n"
         << "" << "\n"
         << "Physical Surface(16) = {11};" << "\n"
         << "Physical Surface(17) = {15};" << "\n"
         << "Physical Surface(18) = {5};" << "\n"
         << "Physical Surface(19) = {13};" << "\n"
         << "Physical Volume(22) = {21};" << "\n";
    return ostr.str();
}

#if defined( LIFE_INSTANTIATION_MODE )

// Instantiations
// 1D
template class GmshSimplexDomain<1,1>;
template class GmshSimplexDomain<1,2>;
// 2D
template class GmshSimplexDomain<2,1>;
template class GmshSimplexDomain<2,2>;
// 3D
template class GmshSimplexDomain<3,1>;
template class GmshSimplexDomain<3,2>;
#endif // LIFE_INSTANTIATION_MODE

}
#endif // __GMSHSIMPLEXDOMAIN_CPP
