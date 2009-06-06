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
const uint16_type triangle::__e2p[6] =
    {
        1, 2, // edge 0
        2, 0, // edge 1
        0, 1  // edge 2
    };
const uint16_type triangle::__f2p[3] =
    {
        0, // point 0
        1, // point 1
        2  // point 2
    };
const uint16_type triangle::__f2e[3] =
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
const uint16_type tetra::__e2p[12] =
    {
        1, 2,    // edge 0
        2, 0,    // edge 1
        0, 1,    // edge 2
        0, 3,    // edge 3
        1, 3,    // edge 4
        2, 3     // edge 5
    };
const uint16_type tetra::__f2p[12] =
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
const uint16_type tetra::__f2e[12] =
    {
        0, 5, 4, // face 0
        1, 3, 5, // face 1
        2, 4, 3, // face 2
        2, 0, 1  // face 3
    };
const int16_type tetra::__f2e_orientation[12] =
    {
        1, 1,-1, // face 0
        1, 1,-1, // face 1
        1, 1,-1, // face 2
        1, 1, 1  // face 3
    };

} // details
/// \endcond

} // Life
