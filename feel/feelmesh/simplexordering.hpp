/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-11-28

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
   \file SimplexOrdering.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-11-28
 */
#ifndef __SimplexOrdering_H
#define __SimplexOrdering_H 1

#include <feel/feelmesh/lineordering.hpp>

namespace Feel
{
/// \cond DETAIL
namespace details
{
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

struct triangle
{
    static constexpr uint16_type f2e( uint16_type /*f*/, uint16_type e )
    {
        return __f2e[e];
    }
    // f2eLoc : (num face,num edge glob)->num edge loc in face
    static constexpr uint16_type f2eLoc( uint16_type /*f*/, uint16_type e )
    {
        return __f2e[e];
    }
    static constexpr uint16_type f2p( uint16_type o, uint16_type /*f*/, uint16_type p )
    {
        return __f2p[p];
    }

    static constexpr uint16_type e2p( uint16_type o, uint16_type e, uint16_type p )
    {
        return 
            (o==1)?e2p( e,p,boost::mpl::int_<1>() ):
            ((o==2)?e2p( e,p,boost::mpl::int_<2>() ):
             ((o==3)?e2p( e,p,boost::mpl::int_<3>() ):
              ((o==4)?e2p( e,p,boost::mpl::int_<4>() ):
               ((o==5)?e2p( e,p,boost::mpl::int_<5>() ): e2p( e,p,boost::mpl::int_<1>() )))));
    }

    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<1> )
    {
        return __e2p_order1[2*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<2> )
    {
        return __e2p_order2[3*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<3> )
    {
        return __e2p_order3[4*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<4> )
    {
        return __e2p_order4[5*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<5> )
        {
            return __e2p_order5[6*e+p];
        }


    static constexpr uint16_type __e2p_order1[6] =
    {
        1, 2, // edge 0
        2, 0, // edge 1
        0, 1  // edge 2
    };

    static constexpr  uint16_type __e2p_order2[9] =
    {
        1, 2, 3, // edge 0
        2, 0, 4, // edge 1
        0, 1, 5  // edge 2
    };
    
    static constexpr uint16_type __e2p_order3[12] =
    {
        1, 2, 3, 4, // edge 0
        2, 0, 5, 6, // edge 1
        0, 1, 7, 8  // edge 2
    };
    
    static constexpr uint16_type __e2p_order4[15] =
    {
        1, 2, 3, 4, 5, // edge 0
        2, 0, 6, 7, 8, // edge 1
        0, 1, 9,10,11  // edge 2
    };
    
    static constexpr uint16_type __e2p_order5[18] =
    {
        1, 2, 3, 4, 5, 6, // edge 0
        2, 0, 7, 8, 9,10, // edge 1
        0, 1,11,12,13,14  // edge 2
    };
    
    static constexpr uint16_type __f2p[21] =
    {
        0, // point 0 - vertex
        1, // point 1 - vertex
        2, // point 2 - vertex
        3, 4, 5, 6, 7, 8, 9, // edge
        10, 11, 12, 13, 14,  // edge
        15, 16, 17, 18, 19, 20 // interior points
        
    };


    static constexpr uint16_type __f2e[3] =
    {
        0, // edge 0
        1, // edge 1
        2  // edge 2
    };

    std::vector<uint16_type> entity( uint16_type /*topo_dim*/, uint16_type id ) const
    {
        std::vector<uint16_type> __entity = { __e2p_order1[2*id], __e2p_order1[2*id+1] };
        return std::move(__entity);
    }
};
constexpr uint16_type triangle::__e2p_order1[6];
constexpr uint16_type triangle::__e2p_order2[9];
constexpr uint16_type triangle::__e2p_order3[12];
constexpr uint16_type triangle::__e2p_order4[15];
constexpr uint16_type triangle::__e2p_order5[18];

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
struct tetra
{
    static constexpr uint16_type f2e( uint16_type f, uint16_type e )
    {
        return __f2e[3*f+e];
    }

    // f2eLoc : (num face,num edge glob)->num edge loc in face
    static constexpr uint16_type f2eLoc( uint16_type f, uint16_type e )
    {
        return __f2eLoc[6*f+e];
    }

    static int16_type f2eOrientation( uint16_type f, uint16_type e )
    {
        return __f2e_orientation[3*f+e];
    }



    static constexpr uint16_type f2p( uint16_type o, uint16_type e, uint16_type p )
    {
        return 
            (o==1)?f2p( e,p,boost::mpl::int_<1>() ):
            ((o==2)?f2p( e,p,boost::mpl::int_<2>() ):
             ((o==3)?f2p( e,p,boost::mpl::int_<3>() ):
              ((o==4)?f2p( e,p,boost::mpl::int_<4>() ):f2p( e,p,boost::mpl::int_<1>() ))));
               
    }
    static constexpr uint16_type f2p( uint16_type f, uint16_type p, boost::mpl::int_<1> )
    {
        return __f2p_order1[3*f+p];
    }
    static constexpr uint16_type f2p( uint16_type f, uint16_type p, boost::mpl::int_<2> )
    {
        return __f2p_order2[6*f+p];
    }
    static constexpr uint16_type f2p( uint16_type f, uint16_type p, boost::mpl::int_<3> )
    {
        return __f2p_order3[10*f+p];
    }
    static constexpr uint16_type f2p( uint16_type f, uint16_type p, boost::mpl::int_<4> )
    {
        return __f2p_order4[15*f+p];
    }

    static constexpr uint16_type e2p( uint16_type o, uint16_type e, uint16_type p )
    {
        return 
            (o==1)?e2p( e,p,boost::mpl::int_<1>() ):
            ((o==2)?e2p( e,p,boost::mpl::int_<2>() ):
             ((o==3)?e2p( e,p,boost::mpl::int_<3>() ):
              ((o==4)?e2p( e,p,boost::mpl::int_<4>() ):
               ((o==5)?e2p( e,p,boost::mpl::int_<5>() ):e2p( e,p,boost::mpl::int_<1>() ) ))));
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<1> )
    {
        return __e2p_order1[2*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<2> )
    {
        return __e2p_order2[3*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<3> )
    {
        return __e2p_order3[4*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<4> )
    {
        return __e2p_order4[5*e+p];
    }
    static constexpr uint16_type e2p( uint16_type e, uint16_type p,boost::mpl::int_<5> )
    {
        return __e2p_order5[6*e+p];
    }


    std::vector<uint16_type> entity( uint16_type topo_dim, uint16_type id ) const
    {
        std::vector<uint16_type> __entity( topo_dim+1 );

        if ( topo_dim == 1 )
        {
            __entity[0] = __e2p_order1[2*id];
            __entity[1] = __e2p_order1[2*id+1];

        }

        else if ( topo_dim == 2 )
        {
            __entity[0] = __f2p_order1[3*id];
            __entity[1] = __f2p_order1[3*id+1];
            __entity[2] = __f2p_order1[3*id+2];
        }

        return __entity;
    }

static constexpr uint16_type __e2p_order1[12] =
{
    1, 2,    // edge 0
    2, 0,    // edge 1
    0, 1,    // edge 2
    0, 3,    // edge 3
    1, 3,    // edge 4
    2, 3     // edge 5
};

static constexpr uint16_type __e2p_order2[18] =
{
    1, 2, 4,   // edge 0
    2, 0, 5,   // edge 1
    0, 1, 6,   // edge 2
    0, 3, 7,   // edge 3
    1, 3, 8,   // edge 4
    2, 3, 9    // edge 5
};

static constexpr uint16_type __e2p_order3[24] =
{
    1, 2, 4, 5,   // edge 0
    2, 0, 6, 7,   // edge 1
    0, 1, 8, 9,   // edge 2
    0, 3,10,11,   // edge 3
    1, 3,12,13,   // edge 4
    2, 3,14,15    // edge 5
};

static constexpr uint16_type __e2p_order4[30] =
{
    1, 2,  4,  5,  6,  // edge 0
    2, 0,  7,  8,  9,  // edge 1
    0, 1, 10, 11, 12,  // edge 2
    0, 3, 13, 14, 15,  // edge 3
    1, 3, 16, 17, 18,  // edge 4
    2, 3, 19, 20, 21   // edge 5
};


static constexpr uint16_type __e2p_order5[36] =
{
    1, 2, 4, 5, 6, 7,      // edge 0
    2, 0, 8, 9, 10, 11,    // edge 1
    0, 1, 12, 13, 14, 15,  // edge 2
    0, 3, 16, 17, 18, 19,  // edge 3
    1, 3, 20, 21, 22, 23,  // edge 4
    2, 3, 25, 26, 27, 28   // edge 5
};


static constexpr uint16_type __f2p_order1[12] =
{
    1, 2, 3, // face 0
    0, 2, 3, // face 1
    0, 1, 3, // face 2
    0, 1, 2  // face 3
};

static constexpr uint16_type __f2p_order2[24] =
{
    1, 2, 3, 9, 8, 4, // face 0
    0, 2, 3, 9, 7, 5, // face 1
    0, 1, 3, 8, 7, 6, // face 2
    0, 1, 2, 4, 5, 6  // face 3
};

static constexpr uint16_type __f2p_order3[40] =
{
    1, 2, 3,  14, 15, 13, 12, 4, 5, 16, // face 0
    0, 2, 3,  14, 15, 11, 10, 7, 6, 17, // face 1
    0, 1, 3,  12, 13, 11, 10, 8, 9, 18, // face 2
    0, 1, 2,   4,  5,  6,  7, 8, 9, 19  // face 3
};

static constexpr uint16_type __f2p_order4[60] =
{
    1, 2, 3, 19, 20, 21, 18, 17, 16,  4,  5,  6, 22, 23, 24, // face 0
    0, 2, 3, 19, 20, 21, 15, 14, 13,  9,  8,  7, 25, 26, 27, // face 1
    0, 1, 3, 16, 17, 18, 15, 14, 13, 10, 11, 12, 28, 29, 30, // face 2
    0, 1, 2,  4,  5,  6,  7,  8,  9, 10, 11, 12, 31, 32, 33  // face 3
};

static constexpr uint16_type __f2e[12] =
{
    5, 4, 0, // face 0
    5, 3, 1, // face 1
    4, 3, 2, // face 2
    0, 1, 2  // face 3
};

//99 for bad value
static constexpr uint16_type  __f2eLoc[24] =
{
    2 , 99, 99, 99,  1,  0, // face 0
    99,  2, 99,  1, 99,  0, // face 1
    99, 99,  2,  1,  0, 99, // face 2
    0 ,  1,  2, 99, 99, 99  // face 3
};

static constexpr int16_type __f2e_orientation[12] =
{
    1, -1, 1, // face 0
    1, -1,-1, // face 1
    1, -1, 1, // face 2
    1,  1, 1  // face 3
};

};

constexpr uint16_type tetra::__e2p_order1[12];
constexpr uint16_type tetra::__e2p_order2[18];
constexpr uint16_type tetra::__e2p_order3[24];
constexpr uint16_type tetra::__e2p_order4[30];
constexpr uint16_type tetra::__e2p_order5[36];
constexpr uint16_type tetra::__f2p_order1[12];
constexpr uint16_type tetra::__f2p_order2[24];
constexpr uint16_type tetra::__f2p_order3[40];
constexpr uint16_type tetra::__f2p_order4[60];
constexpr uint16_type tetra::__f2e[12];

} // details
/// \endcond

}
#endif /* __SimplexOrdering_H */
