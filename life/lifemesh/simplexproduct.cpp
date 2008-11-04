/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-02-20

  Copyright (C) 2006 EPFL

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
   \file simplexproduct.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-02-20
 */
#include <life/lifemesh/simplexproduct.hpp>

namespace Life
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
const uint16_type quad::__e2p[8] =
    {
        0, 1, // points in edge 0
        1, 2, // points in edge 1
        2, 3, // points in edge 2
        3, 0  // points in edge 3
    };
const uint16_type quad::__f2p[4] =
    {
        0, // point 0
        1, // point 1
        2, // point 2
        3  // point 3
    };
const uint16_type quad::__f2e[4] =
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
const uint16_type hexa::__e2p[24] =
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
// face to point relation
const uint16_type hexa::__f2p[24] =
    {
        0, 1, 2, 3, // face 0
        0, 1, 5, 4, // face 1
        1, 2, 6, 5, // face 2
        2, 3, 7, 6, // face 3
        3, 0, 4, 7, // face 4
        4, 5, 6, 7  // face 5
    };

// face to edge relation
const uint16_type hexa::__f2e[24] =
    {
        0,  1,  2,  3, // face 0
        0,  4,  5,  6, // face 1
        1,  7,  8,  4, // face 2
        2,  9, 10,  7, // face 3
        3,  6, 11,  9, // face 4
        5,  8, 10, 11  // face 5
    };

// edge permutation in face
const int16_type hexa::__f2e_permutation[24] =
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

} // Life
