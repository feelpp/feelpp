/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-02-20

  Copyright (C) 2006 EPFL

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
   \file simplex.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-02-20
 */
#include <life/lifemesh/simplex.hpp>

namespace Life
{
/// \cond DETAIL
namespace details
{

/**
 * \class line
 */
const uint16_type line::__e2p[2] = { 0, 1 };

/**
 * \class triangle

           2
           |\
           | \
         1 |  \ 0
           |   \
           |    \
           |     \
         0 ------- 1
              2
*/
template<uint16_type Order>
const uint16_type triangle<Order>::__e2p_order1[6] =
    {
        1, 2, // edge 0
        2, 0, // edge 1
        0, 1  // edge 2
    };
template<uint16_type Order>
const uint16_type triangle<Order>::__e2p_order2[9] =
    {
        1, 2, 3, // edge 0
        2, 0, 4, // edge 1
        0, 1, 5  // edge 2
    };
template<uint16_type Order>
const uint16_type triangle<Order>::__e2p_order3[12] =
    {
        1, 2, 3, 4, // edge 0
        2, 0, 5, 6, // edge 1
        0, 1, 7, 8  // edge 2
    };
template<uint16_type Order>
const uint16_type triangle<Order>::__e2p_order4[15] =
    {
        1, 2, 3, 4, 5, // edge 0
        2, 0, 6, 7, 8, // edge 1
        0, 1, 9,10,11  // edge 2
    };
template<uint16_type Order>
const uint16_type triangle<Order>::__e2p_order5[18] =
    {
        1, 2, 3, 4, 5, 6, // edge 0
        2, 0, 7, 8, 9,10, // edge 1
        0, 1,11,12,13,14  // edge 2
    };


template<uint16_type Order>
const uint16_type triangle<Order>::__f2p[3] =
    {
        0, // point 0
        1, // point 1
        2  // point 2
    };

template<uint16_type Order>
const uint16_type triangle<Order>::__f2e[3] =
    {
        0, // edge 0
        1, // edge 1
        2  // edge 2
    };

/**
   \class tetra

           3
           |\\
           | \\
           |  \\
           |   \\
           |    \\
           |     \.2
           |    . \\
           |  .    \\
           |.       \!
         0 ----------1
*/
template<uint16_type Order >
const uint16_type tetra<Order>::__e2p_order1[12] =
{
    1, 2,    // edge 0
    2, 0,    // edge 1
    0, 1,    // edge 2
    0, 3,    // edge 3
    1, 3,    // edge 4
    2, 3     // edge 5
};

template<uint16_type Order >
const uint16_type tetra<Order>::__e2p_order2[18] =
{
    1, 2, 4,   // edge 0
    2, 0, 5,   // edge 1
    0, 1, 6,   // edge 2
    0, 3, 7,   // edge 3
    1, 3, 8,   // edge 4
    2, 3, 9    // edge 5
};
#if 0 // GMSH ORDERING
    {
        1, 2, 9,   // edge 0
        2, 0, 8,   // edge 1
        0, 1, 6,   // edge 2
        0, 3, 4,   // edge 3
        1, 3, 5,   // edge 4
        2, 3, 7    // edge 5
    };
#endif
template<uint16_type Order >
const uint16_type tetra<Order>::__e2p_order3[24] =
{
    1, 2, 4, 5,   // edge 0
    2, 0, 6, 7,   // edge 1
    0, 1, 8, 9,   // edge 2
    0, 3,10,11,   // edge 3
    1, 3,12,13,   // edge 4
    2, 3,14,15    // edge 5
};
#if 0 // GMSH ORDERING
    {
        1, 2, 15, 14,   // edge 0
        2, 0, 12, 13,   // edge 1
        0, 1,  9,  8,   // edge 2
        0, 3,  5,  4,   // edge 3
        1, 3,  7,  6,   // edge 4
        2, 3, 11, 10    // edge 5
    };
#endif
template<uint16_type Order >
const uint16_type tetra<Order>::__e2p_order4[30] =
    {
        1, 2, 21, 20, 19,  // edge 0
        2, 0, 17, 17, 18,  // edge 1
        0, 1, 12, 11, 10,  // edge 2
        0, 3,  6,  5,  4,  // edge 3
        1, 3,  9,  8,  7,  // edge 4
        2, 3, 15, 14, 13  // edge 5
    };

template<uint16_type Order >
const uint16_type tetra<Order>::__e2p_order5[36] =
    {
        1, 2, 11, 10,  9,  8,  // edge 0
        2, 0, 20, 21, 22, 23,  // edge 1
        0, 1, 15, 14, 13, 12,  // edge 2
        0, 3,  7,  6,  5,  4,  // edge 3
        1, 3, 11, 10,  9,  8,  // edge 4
        2, 3, 19, 18, 17, 16   // edge 5
    };


template<uint16_type Order >
const uint16_type tetra<Order>::__f2p[12] =
    {
#if 0 // Old numerotation
        1, 3, 2, // face 0
        2, 3, 0, // face 1
        3, 1, 0, // face 2
        0, 1, 2  // face 3
#else // New numerotation needed to construct multidomain high order
      // modal basis
        1, 2, 3, // face 0
        0, 2, 3, // face 1
        0, 1, 3, // face 2
        0, 1, 2  // face 3
#endif
    };

template<uint16_type Order >
const uint16_type tetra<Order>::__f2e[12] =
    {
        0, 5, 4, // face 0
        1, 3, 5, // face 1
        2, 4, 3, // face 2
        2, 0, 1  // face 3
    };

template<uint16_type Order >
const int16_type tetra<Order>::__f2e_orientation[12] =
    {
        1, 1,-1, // face 0
        1, 1,-1, // face 1
        1, 1,-1, // face 2
        1, 1, 1  // face 3
    };

} // details
/// \endcond

} // Life


template class Life::details::triangle<1>;
template class Life::details::triangle<2>;
template class Life::details::triangle<3>;
template class Life::details::triangle<4>;
template class Life::details::triangle<5>;

template class Life::details::tetra<1>;
template class Life::details::tetra<2>;
template class Life::details::tetra<3>;
template class Life::details::tetra<4>;
template class Life::details::tetra<5>;
