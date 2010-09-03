/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

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
   \file hypercube.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-02-20
 */
#include <feel/feelmesh/hypercube.hpp>

namespace Feel
{

/// \cond DETAIL
namespace details
{

/**
 * \class quad

\code
          3   2   2
           |-----|
           |     |
         3 |     | 1
           |     |
         0 ------- 1
              0
\endcode
*/
template<uint16_type Order >
const uint16_type quad<Order>::__e2p_order1[8] =
    {
        0, 1, // points in edge 0
        1, 2, // points in edge 1
        2, 3, // points in edge 2
        3, 0  // points in edge 3
    };

template<uint16_type Order >
const uint16_type quad<Order>::__e2p_order2[12] =
    {
        0, 1, 4, // points in edge 0
        1, 2, 5, // points in edge 1
        2, 3, 6, // points in edge 2
        3, 0, 7  // points in edge 3
    };

template<uint16_type Order >
const uint16_type quad<Order>::__e2p_order3[16] =
    {
        0, 1,  4,  5, // points in edge 0
        1, 2,  6,  7, // points in edge 1
        2, 3,  8,  9, // points in edge 2
        3, 0, 10, 11  // points in edge 3
    };


template<uint16_type Order >
const uint16_type quad<Order>::__f2p[4] =
    {
        0, // point 0
        1, // point 1
        2, // point 2
        3  // point 3
    };
template<uint16_type Order >
const uint16_type quad<Order>::__f2e[4] =
    {
        0, // edge 0
        1, // edge 1
        2, // edge 2
        3  // edge 3
    };


/**
   \class hexa

\code
        7-------6
       /.      /|
      / .     / |
     4_______5  |
     |  .    |  |
     |  3....|..2
     | .     | /
     |.      |/
     0_______1
\endcode
*/


// edge to point relation
template<uint16_type Order >
const uint16_type hexa<Order>::__e2p_order1[24] =
    {
        0, 1,    // edge 0
        1, 2,    // edge 1
        2, 3,    // edge 2
        3, 0,    // edge 3
        1, 5,    // edge 4
        5, 4,    // edge 5
        4, 0,    // edge 6
        2, 6,    // edge 7
        6, 5,    // edge 8
        3, 7,    // edge 9
        7, 6,    // edge 10
        4, 7     // edge 11
    };

template<uint16_type Order >
const uint16_type hexa<Order>::__e2p_order2[36] =
    {
        0, 1, 8,    // edge 0
        1, 2, 9,    // edge 1
        2, 3, 10,   // edge 2
        3, 0, 11,   // edge 3
        1, 5, 12,   // edge 4
        5, 4, 13,   // edge 5
        4, 0, 14,   // edge 6
        2, 6, 15,   // edge 7
        6, 5, 16,   // edge 8
        3, 7, 17,   // edge 9
        7, 6, 18,   // edge 10
        4, 7, 19    // edge 11
    };


template<uint16_type Order >
const uint16_type hexa<Order>::__e2p_order3[48] =
    {
        0, 1, 8, 9,     // edge 0
        1, 2, 10, 11,   // edge 1
        2, 3, 12, 13,   // edge 2
        3, 0, 14, 15,   // edge 3
        1, 5, 16, 17,   // edge 4
        5, 4, 18, 19,   // edge 5
        4, 0, 20, 21,   // edge 6
        2, 6, 22, 23,   // edge 7
        6, 5, 24, 25,   // edge 8
        3, 7, 26, 27,   // edge 9
        7, 6, 28, 29,   // edge 10
        4, 7, 30, 31    // edge 11
    };


// face to point relation
template<uint16_type Order >
const uint16_type hexa<Order>::__f2p_order1[24] =
    {
        0, 1, 2, 3, // face 0
        0, 1, 5, 4, // face 1
        1, 2, 6, 5, // face 2
        2, 3, 7, 6, // face 3
        3, 0, 4, 7, // face 4
        4, 5, 6, 7  // face 5
    };

template<uint16_type Order >
const uint16_type hexa<Order>::__f2p_order2[54] =
    {
        0, 1, 2, 3,  8, 9,  10, 11, 20, // face 0
        0, 1, 5, 4,  8, 12, 13, 14, 21, // face 1
        1, 2, 6, 5,  9, 15, 16, 12, 22, // face 2
        2, 3, 7, 6, 10, 17, 18, 15, 23, // face 3
        3, 0, 4, 7, 11, 14, 19, 17, 24, // face 4
        4, 5, 6, 7, 13, 16, 18, 19, 25  // face 5
    };


template<uint16_type Order >
const uint16_type hexa<Order>::__f2p_order3[96] =
    {
        0, 1, 2, 3,  8,  9, 10, 11, 12, 13, 14, 15, 32, 33, 34, 35, // face 0
        0, 1, 5, 4,  8,  9, 16, 17, 18, 19, 20, 21, 36, 37, 38, 39, // face 1
        1, 2, 6, 5, 10, 11, 22, 23, 24, 25, 17, 16, 40, 41, 42, 43, // face 2
        2, 3, 7, 6, 12, 13, 26, 27, 28, 29, 23, 22, 44, 45, 46, 47, // face 3
        3, 0, 4, 7, 14, 15, 21, 20, 30, 31, 27, 26, 48, 49, 50, 51, // face 4
        4, 5, 6, 7, 19, 18, 25, 24, 29, 28, 31, 30, 52, 53, 54, 55  // face 5
    };

// face to edge relation
template<uint16_type Order >
const uint16_type hexa<Order>::__f2e[24] =
    {
        0,  1,  2,  3, // face 0
        0,  4,  5,  6, // face 1
        1,  7,  8,  4, // face 2
        2,  9, 10,  7, // face 3
        3,  6, 11,  9, // face 4
        5,  8, 10, 11  // face 5
    };

// edge permutation in face
template<uint16_type Order >
const int16_type hexa<Order>::__f2e_permutation[24] =
    {
        1,  1,  1,  1, // face 0
        1,  1,  1,  1, // face 1
        1,  1,  1, -1, // face 2
        1, -1,  1, -1, // face 3
        1, -1,  1, -1, // face 4
       -1, -1, -1, -1  // face 5
    };

} // details
/// \endcond

} // Feel

template class Feel::details::quad<1>;
template class Feel::details::quad<2>;
template class Feel::details::quad<3>;
template class Feel::details::quad<4>;
template class Feel::details::quad<5>;

template class Feel::details::hexa<1>;
template class Feel::details::hexa<2>;
template class Feel::details::hexa<3>;
template class Feel::details::hexa<4>;
template class Feel::details::hexa<5>;
