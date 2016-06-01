/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Jun 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_MESHSTRUCTURED_HPP
#define FEELPP_MESHSTRUCTURED_HPP 1

namespace Feel {

/**
 *
 */
class MeshStructured: public Mesh<Hypercube<2>>
{

  public:
    using super = Mesh<Hypercube<2>>;
    using point_type = super::point_type;
    using element_type = super::element_type;
    using node_type = super::node_type;
    MeshStructured() = default;
    MeshStructured( MeshStructured const& ) = default;
    MeshStructured( MeshStructured && ) = default;
    MeshStructured& operator=()( MeshStructured const& ) = default;
    MeshStructured& operator=()( MeshStructured && ) = default;

    MeshStructured( int nx, int ny );
};

MeshStructured::MeshStructured( int nx, int ny,
                                double pixelsize,
                                WorldComm const& wc )
    :
    super( wc ),
    M_nx( nx ),
    M_ny( ny ),
    M_pixelsize( pixelsize )
{
    // origin at (0,0)
    node_type coords( 2 );
    rank_type partId = this->worldComm().localRank();
    for( int i = 0; i <= M_nx; ++i )
        for( int j = 0; j <= M_ny; ++j )
        {
            // point
            int ptid = (M_ny+1)*i+j;
            coords[0]=pixelsize*i;
            coords[1]=pixelsize*j;
            point_type pt( ptid, coords, false );
            pt.setProcessIdInPartition( partId );
            pt.setProcessId( partId );
            mesh->addPoint( pt );
        }

    int eltid = 0;
    rank_type partId = mesh->worldComm().localRank();
    element_type e;

    for( int i = 0; i < M_nx; ++i )
        for( int j = 0; j < M_ny; ++j )
        {
            int eid = (M_ny)*i+j;
            e.setId( eid );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            //
            int ptid[4] = { (M_ny+1)*(i+1)+j,   // 0
                            (M_ny+1)*(i+1)+j+1, // 1
                            (M_ny+1)*i+j+1,     // 2
                            (M_ny+1)*i+j        // 3
            };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
            mesh->addElement( e, true );
        }
}
}

#endif
