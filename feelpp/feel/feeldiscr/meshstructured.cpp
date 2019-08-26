/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 04 Jan 2017

 Copyright (C) 2017 Feel++ Consortium

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
#include <feel/feelconfig.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <feel/feeldiscr/meshstructured.hpp>

typedef boost::geometry::model::d2::point_xy<double> poly_point_type;
typedef boost::geometry::model::segment<poly_point_type> segment_type;
typedef boost::geometry::model::polygon<poly_point_type> polygon_type;

namespace Feel
{

MeshStructured::MeshStructured( int nx, int ny, double pixelsize,
                                std::optional<holo3_image<float>> const& cx,
                                std::optional<holo3_image<float>> const& cy,
                                worldcomm_ptr_t const& wc = Environment::worldCommPtr(),
                                bool withCoord,
                                std::string pathPoly, bool withPoly )
    : super( wc ),
      M_nx( nx ),
      M_ny( ny ),
      M_cx( cx ),
      M_cy( cy ),
      M_pixelsize( pixelsize )
{
    if ( !M_cx && !M_cy )
        this->setStructureProperty( "00110" );
    else
        this->setStructureProperty( "00010" );
    VLOG( 1 ) << "nx x ny = " << nx << " x " << ny << "\t" << nx * ny << std::endl;

    rank_type nProc = wc->localSize();
    rank_type partId = wc->localRank();

    // compute parallel distribution with column scheme partitioning
    size_type nTotalPointsByCol = M_ny + ( nProc - 1 );
    size_type nPtByColByProc = ( M_ny + ( nProc - 1 ) ) / nProc;
    if ( nTotalPointsByCol <= nPtByColByProc * ( nProc - 1 ) )
        --nPtByColByProc;
    std::vector<size_type> nPtByCol( nProc, nPtByColByProc );
    if ( nPtByColByProc * nProc < nTotalPointsByCol )
    {
        size_type nPointToDistribute = nTotalPointsByCol - nPtByColByProc * nProc;
        for ( size_type k = 0; k < nPointToDistribute; ++k )
        {
            if ( k < nProc )
                ++nPtByCol[nProc - 1 - k];
        }
    }
    std::vector<size_type> startColPtId( nProc, 0 );
    for ( rank_type p = 1; p < nProc; ++p )
        startColPtId[p] = startColPtId[p - 1] + ( nPtByCol[p - 1] - 1 );

#if 0
    if ( this->worldComm().isMasterRank() )
    {
        std::cout<< "startColPtId : ";
        for (int p=0;p<nProc;++p )
            std::cout << " " << startColPtId[p] << "(" << nPtByCol[p] << ") ";
        std::cout<< "\n";
    }
#endif
    std::unordered_map<int, boost::tuple<int, rank_type>> mapGhostElt;
    std::unordered_map<int, int> idStructuredMeshToFeelMesh;
    node_type coords( 2 );

    // polygon build
    std::vector<polygon_type> list_poly;
    //std::cout << list_poly.size() << std::endl;
    polygon_type poly;
    bool inPoly = false;
    bool testInPoly = false;
    if ( withPoly )
    {
        polygon_type p;
        std::ifstream flux( pathPoly );
        int nbPolygons;
        int nbVertices;
        double coordX, coordY;
        flux >> nbPolygons;
        //std::cout<< "nP " << nbPolygons <<  std::endl;
        for ( int j = 0; j < nbPolygons; j++ )
        {
            flux >> nbVertices;
            //std::cout<< "nV: " << nbVertices <<  std::endl;
            for ( int i = 0; i < nbVertices; i++ )
            {
                flux >> coordX;
                flux >> coordY;
                p.outer().push_back( poly_point_type( coordX * pixelsize, coordY * pixelsize ) );
            }
            p.outer().push_back( poly_point_type( p.outer().at( 0 ).x(), p.outer().at( 0 ).y() ) );

            list_poly.push_back( p );
            p.clear();
        }
    }
    //std::cout<< "poly build end"<< std::endl;

    // active points
    for ( int j = startColPtId[partId]; j < startColPtId[partId] + nPtByCol[partId]; ++j )
    {
        for ( int i = 0; i < M_nx; ++i )
        {
            int ptid = (M_ny)*i + j;
            if ( withCoord && M_cx && M_cy )
            {
                coords[0] = (*M_cy)( j, i );
                coords[1] = (*M_cx)( j, i );
            }
            else
            {
                coords[0] = M_pixelsize * j;
                coords[1] = M_pixelsize * i;
            }
            if ( withPoly )
            {
                poly_point_type p( coords[0], coords[1] );
                for ( int k = 0; k < list_poly.size(); k++ )
                {
                    poly = list_poly[k];
                    testInPoly = boost::geometry::within( p, poly );
                    inPoly = ( inPoly ) || ( testInPoly );
                }
            }
            if ( ( !withPoly ) || ( !inPoly ) )
            {
                point_type pt( ptid, coords );
                pt.setProcessId( partId );
                pt.setProcessIdInPartition( partId );
                this->addPoint( pt );
            }
            inPoly = false;
        }
    }

    // ghost points (not belong to active partition)
    std::vector<size_type> ghostPointColDesc;
    if ( partId > 0 )
        ghostPointColDesc.push_back( startColPtId[partId] - 1 );
    if ( partId < ( nProc - 1 ) )
        ghostPointColDesc.push_back( startColPtId[partId] + nPtByCol[partId] );
    for ( int j : ghostPointColDesc )
    {
        for ( int i = 0; i < M_nx; ++i )
        {
            int ptid = (M_ny)*i + j;
            if ( withCoord && M_cx && M_cy )
            {
                coords[0] = (*M_cy)( j, i );
                coords[1] = (*M_cx)( j, i );
            }
            else
            {
                coords[0] = M_pixelsize * j;
                coords[1] = M_pixelsize * i;
            }
            point_type pt( ptid, coords );
            pt.setProcessId( invalid_rank_type_value );
            pt.setProcessIdInPartition( partId );
            this->addPoint( pt );
        }
    }

    // active elements
    size_type startColPtIdForElt = startColPtId[partId];
    size_type stopColPtIdForElt = startColPtId[partId] + ( nPtByCol[partId] - 1 );
    for ( int j = startColPtIdForElt; j < stopColPtIdForElt; ++j )
    {
        rank_type neighborProcessId = invalid_rank_type_value;
        if ( partId > 0 )
            if ( j == startColPtIdForElt )
                neighborProcessId = partId - 1;
        if ( partId < ( nProc - 1 ) )
            if ( j == ( stopColPtIdForElt - 1 ) )
                neighborProcessId = partId + 1;

        for ( int i = 0; i < M_nx - 1; ++i )
        {
            if ( withPoly )
            {
                poly_point_type p1;
                poly_point_type p2;
                poly_point_type p3;
                poly_point_type p4;
                if ( withCoord && M_cx && M_cy )
                {
                    p1 = poly_point_type( (*M_cy)( j, i + 1 ), (*M_cx)( j, i + 1 ) );
                    p2 = poly_point_type( (*M_cy)( j + 1, i + 1 ), (*M_cx)( j + 1, i + 1 ) );
                    p3 = poly_point_type( (*M_cy)( j + 1, i ), (*M_cx)( j + 1, i ) );
                    p4 = poly_point_type( (*M_cy)( j, i ), (*M_cx)( j, i ) );
                }
                else
                {
                    p1 = poly_point_type( M_pixelsize * j, M_pixelsize * ( i + 1 ) );
                    p2 = poly_point_type( M_pixelsize * ( j + 1 ), M_pixelsize * ( i + 1 ) );
                    p3 = poly_point_type( M_pixelsize * ( j + 1 ), M_pixelsize * i );
                    p4 = poly_point_type( M_pixelsize * j, M_pixelsize * i );
                }

                for ( int k = 0; k < list_poly.size(); k++ )
                {
                    poly = list_poly[k];
                    bool testInPoly1 = boost::geometry::within( p1, poly );
                    bool testInPoly2 = boost::geometry::within( p2, poly );
                    bool testInPoly3 = boost::geometry::within( p3, poly );
                    bool testInPoly4 = boost::geometry::within( p4, poly );
                    inPoly = ( inPoly ) || ( testInPoly1 ) || ( testInPoly2 ) || ( testInPoly3 ) || ( testInPoly4 );
                }
            }
            if ( ( !withPoly ) || ( !inPoly ) )
            {
                size_type eid = ( M_ny - 1 ) * i + j; // StructuredMesh Id
                element_type e;
                e.setMarker( 1, 1 );
                e.setProcessIdInPartition( partId );
                e.setProcessId( partId );
                size_type ptid[4] = {( M_ny ) * ( i + 1 ) + j,     // 0
                                     ( M_ny ) * ( i + 1 ) + j + 1, // 1
                                     (M_ny)*i + j + 1,             // 2
                                     (M_ny)*i + j};                // 3

                for ( uint16_type k = 0; k < 4; ++k )
                    e.setPoint( k, this->point( ptid[k] ) );
                if ( neighborProcessId != invalid_rank_type_value )
                    e.addNeighborPartitionId( neighborProcessId );

                auto [eit,inserted] = this->addElement( e, true ); // e.id() is defined by Feel++
                idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eit->second.id() ) );
            }
            inPoly = false;
        }
    }

    // ghost elements (touch active partition with a point)
    std::vector<std::pair<size_type, rank_type>> ghostEltColDesc;
    if ( partId > 0 )
        ghostEltColDesc.push_back( std::make_pair( startColPtId[partId] - 1, partId - 1 ) );
    if ( partId < ( nProc - 1 ) )
        ghostEltColDesc.push_back( std::make_pair( startColPtId[partId] + nPtByCol[partId] - 1, partId + 1 ) );
    for ( auto const& ghostEltCol : ghostEltColDesc )
    {
        int j = ghostEltCol.first;
        rank_type partIdGhost = ghostEltCol.second;
        for ( int i = 0; i < M_nx - 1; ++i )
        {
            if ( withPoly )
            {
                poly_point_type p1;
                poly_point_type p2;
                poly_point_type p3;
                poly_point_type p4;
                if ( withCoord && M_cx && M_cy )
                {
                    p1 = poly_point_type( (*M_cy)( j, i + 1 ), (*M_cx)( j, i + 1 ) );
                    p2 = poly_point_type( (*M_cy)( j + 1, i + 1 ), (*M_cx)( j + 1, i + 1 ) );
                    p3 = poly_point_type( (*M_cy)( j + 1, i ), (*M_cx)( j + 1, i ) );
                    p4 = poly_point_type( (*M_cy)( j, i ), (*M_cx)( j, i ) );
                }
                else
                {
                    p1 = poly_point_type( M_pixelsize * j, M_pixelsize * ( i + 1 ) );
                    p2 = poly_point_type( M_pixelsize * ( j + 1 ), M_pixelsize * ( i + 1 ) );
                    p3 = poly_point_type( M_pixelsize * ( j + 1 ), M_pixelsize * i );
                    p4 = poly_point_type( M_pixelsize * j, M_pixelsize * i );
                }
                for ( int k = 0; k < list_poly.size(); k++ )
                {
                    poly = list_poly[k];
                    bool testInPoly1 = boost::geometry::within( p1, poly );
                    bool testInPoly2 = boost::geometry::within( p2, poly );
                    bool testInPoly3 = boost::geometry::within( p3, poly );
                    bool testInPoly4 = boost::geometry::within( p4, poly );
                    inPoly = ( inPoly ) || ( testInPoly1 ) || ( testInPoly2 ) || ( testInPoly3 ) || ( testInPoly4 );
                }
            }
            if ( ( !withPoly ) || ( !inPoly ) )
            {
                size_type eid = ( M_ny - 1 ) * i + j; // StructuredMesh Id
                element_type e;
                e.setMarker( 1, 1 );
                e.setProcessIdInPartition( partId );
                e.setProcessId( partIdGhost );

                size_type ptid[4] = {( M_ny ) * ( i + 1 ) + j,     // 0
                                     ( M_ny ) * ( i + 1 ) + j + 1, // 1
                                     (M_ny)*i + j + 1,             // 2
                                     (M_ny)*i + j};                // 3

                for ( uint16_type k = 0; k < 4; ++k )
                {
                    CHECK( this->hasPoint( ptid[k] ) ) << "mesh doesnt have this point id : " << ptid[k];
                    e.setPoint( k, this->point( ptid[k] ) );
                }

                auto [eit,inserted] = this->addElement( e, true ); // e.id() is defined by Feel++
                idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eit->second.id() ) );
                mapGhostElt.insert( std::make_pair( eid, boost::make_tuple( idStructuredMeshToFeelMesh[eid], partIdGhost ) ) );
            }
            inPoly = false;
        }
    }

    std::vector<int> nbMsgToRecv( nProc, 0 );
    // ghost elts on left
    if ( partId > 0 )
        nbMsgToRecv[partId - 1] = 1;
    // ghost elts on right
    if ( partId < ( nProc - 1 ) )
        nbMsgToRecv[partId + 1] = 1;

    this->updateGhostCellInfoByUsingNonBlockingComm( idStructuredMeshToFeelMesh, mapGhostElt, nbMsgToRecv );
}

#if 0
MeshStructured::MeshStructured( int nx, int ny, double pixelsize, worldcomm_ptr_t const& wc )
    :
    super( wc ),
    M_nx( nx ),
    M_ny( ny ),
    M_pixelsize( pixelsize )
{
    tic();
    Feel::cout << "nx x ny = " << nx << " x " << ny << "\t" << nx*ny << std::endl;
    // origin at (0,0)
    node_type coords( 2 );
    rank_type partId = wc->localRank();
    std::vector<int> nbMsgToRecv(wc->godSize(),0);
    //std::vector<int> nbMsgToRecv((M_nx-1)*(M_ny-1),0);
    int procSize = wc->ogodSize();
    // cx[rank] = first point Id of partition
    int cx[procSize+1];
    rank_type id;
    std::vector<int>::iterator lowProc;

    GmshOrdering<element_type> ordering;

    for (int tmpProc=0;tmpProc<=procSize;tmpProc++)
    {
        cx[tmpProc]=(tmpProc)*(M_ny-1)/procSize;
        Feel::cout << "cx[" << tmpProc << "] = " << cx[tmpProc] << "\t";
    }
    Feel::cout << "\n";

    /*
     * Nodes creations
     * NO ghost points to avoid costly tests at this moment
     */
#pragma omp parallel
#if 1
    // First column
   if(partId > 0)
   {
       int j = cx[partId];
       std::cout << "j = " << j << std::endl;
       ghosts.clear();
       ghosts.push_back(partId-1);
       for( int i = 0; i < M_nx; ++i )
       {
            // point
            int ptid = (M_ny)*i+j;
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid,
                           coords,
                           i == 0 || i == M_nx-1 || j%((M_ny-1)/procSize) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            //pt.setMarker(1);
            // Actual rank ID
            pt.setProcessId( partId  );
            // NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            //pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
    }
#endif
    int start = (partId == 0) ?  0  : cx[partId]+1;
    int stop = (partId == (procSize-1)) ?  M_ny : cx[partId+1];
    for( int j = start ; j < stop; ++j )
    {
        std::cout << "j = " << j << std::endl;
        ghosts.clear();
        if (j == cx[partId] && j != 0)
            ghosts.push_back(partId-1);
        else if (j == cx[partId+1]-1 && j != M_ny-1)
            ghosts.push_back(partId+1);

        for( int i = 0; i < M_nx; ++i )
        {
            // point
            int ptid = (M_ny)*i+j;
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid,
                           coords,
                           i == 0 || i == M_nx-1 || j%((M_ny-1)/procSize) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            //pt.setMarker(1);
            //// Actual rank ID
            pt.setProcessId( partId );
            //// NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            //pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
    }
#if 1
    // Last column
   if(partId < (procSize-1))
   {
       int j = cx[partId+1];
       //std::cout << "j = " << j << std::endl;
       ghosts.clear();
       ghosts.push_back(partId+1);
       for( int i = 0; i < M_nx; ++i )
       {
            // point
            int ptid = (M_ny)*i+j;
            //std::cout << "ptId = " << ptid << std::endl;
            coords[0]=M_pixelsize*j;
            coords[1]=M_pixelsize*i;
            point_type pt( ptid,
                           coords,
                           i == 0 || i == M_nx-1 || j%((M_ny-1)/procSize) == 0  ); // Is on boundary ?
                           //i == 0 || i == M_nx-1 || j == 0 || j == M_ny-1 ); // Is on boundary ?
            //pt.setMarker(1);
            // Actual rank ID
            pt.setProcessId( partId );
            // NO GHOST POINTS I SAID
            pt.setProcessIdInPartition( partId );
            //pt.setNeighborPartitionIds( ghosts );
            this->addPoint( pt );
        }
   }
#endif
   toc("MeshStructured Nodes", FLAGS_v>0);
   tic();

   /*
    * Generates edges
    */
#if 0
    int feid = 0;
    int nelemG=(M_ny-1)*(M_nx-1);
    int nelemL=(cx[partId+1]-cx[partId])*(M_nx-1);
    if ( partId > 0 )
    {
        for( int i = 0; i < M_nx-2; ++i )
        {
            face_type f;
            int eid =nelemL+(M_nx-1)*(2*partId-1)+i+1;
            //f.setId(feid);
            f.setId(eid);
            f.setProcessIdInPartition( partId );
            f.setProcessId( partId );

            int ptid[2] = { (M_ny)*i+cx[partId],   // 0
                            (M_ny)*(i+1)+cx[partId] // 1
                          };

            for( int k = 0; k < 2; ++k )
            {
                f.setPoint( k, this->point( ptid[k]  ) );
            }

            //f.setOnBoundary( i==0 || i==(M_nx-1) || j%((M_ny-1)/procSize)==0 );
            face_iterator fit;
            bool inserted;
            boost::tie( fit, inserted )=this->addFace( f);
            __idGmshToFeel.insert( std::make_pair(feid,f.id()));
            feid++;
            //mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( f.id(), partId-1 )));
            }

    }
    if ( partId < procSize-1 )
    {
        for( int i = 0; i < M_nx-2; ++i )
        {
            face_type f;
            int eid=nelemL+(M_nx-1)*(2*partId)+i+1;
            //f.setId(feid);
            f.setId(eid);
            f.setProcessIdInPartition( partId );
            f.setProcessId( partId );

            int ptid[2] = { (M_ny)*i+cx[partId+1],   // 0
                            (M_ny)*(i+1)+cx[partId+1] // 1
                          };

            for( int k = 0; k < 2; ++k )
            {
                f.setPoint( k, this->point( ptid[k]  ) );
            }

            //f.setOnBoundary( i==0 || i==(M_nx-1) || j%((M_ny-1)/procSize)==0 );
            face_iterator fit;
            bool inserted;
            boost::tie( fit, inserted )=this->addFace( f);

            __idGmshToFeel.insert( std::make_pair(feid,f.id()));
            feid++;
            //mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( f.id(), partId+1 )));
            }
        }
#endif
int eid=0;
//int eid =(cx[partId+1]-cx[partId])*(M_nx-1);
//int eid =(M_ny-1)*(M_nx-1);
//int nelemG =(M_ny-1)*(M_nx-1);
/*
//#if 0
    tic();
    int nelemG=(M_ny-1)*(M_nx-1);
    int eid=nelemG+1;
    for( int j = cx[partId]; j < cx[partId+1]; j++ )
    {
        for( int i = 0; i < M_nx-1; i++ )
        {
            face_type f;
            f.setId(eid);
            f.setProcessIdInPartition( partId );
            f.setProcessId( partId );

             int ptid[2] = { (M_ny)*i+j,   // 0
                            (M_ny)*i+j+1, // 1
                           };

            for( int k = 0; k < 2; ++k )
            {
                f.setPoint( k, this->point( ptid[k]  ) );
            }

            f.setOnBoundary( i==0 || i==(M_nx-1) || j%((M_ny-1))==0 );
            this->addFace( f );
            eid++;

            face_type f2;
            f2.setId(eid);
            f2.setProcessIdInPartition( partId );
            f2.setProcessId( partId );
            int ptid2[2] = { (M_ny)*i+j,   // 0
                            (M_ny)*(i+1)+j, // 1
                           };

            for( int k = 0; k < 2; ++k )
            {
                f2.setPoint( k, this->point( ptid2[k]  ) );
            }

            f2.setOnBoundary( i==0 || i==(M_nx-1) || j%((M_ny-1))==0 );
            this->addFace( f2 );
            eid++;
        }
        face_type f;
        f.setId(eid);
        f.setProcessIdInPartition( partId );
        f.setProcessId( partId );

        int ptid[2] = { (M_ny)*(M_nx-1)+j,   // 0
                        (M_ny)*(M_nx-1)+j+1, // 1
                      };

        for( int k = 0; k < 2; ++k )
        {
            f.setPoint( k, this->point( ptid[k]  ) );
        }

        f.setOnBoundary( true );
        this->addFace( f );
        eid++;
    }


    // Edges
    int final = partId;
    //if( partId < (procSize-1) )
    //    final = partId+1;
    if ( partId < (procSize-1))
    {
    for( int i = 0; i < M_nx-1; ++i )
    {
           face_type f;
           //int eid=nelemG+(M_nx-1)*(2*partId)+i+1;
           f.setId(eid);
           f.setProcessIdInPartition( partId );
           f.setProcessId( final );
           int ptid[2] = { (M_ny)*i+cx[partId+1],   // 0
                            (M_ny)*(i+1)+cx[partId+1], // 1
                         };

           for( int k = 0; k < 2; ++k )
           {
               f.setPoint( k, this->point( ptid[k]  ) );
           }
           //f.setMarker(1);
           f.setOnBoundary( true );
           this->addFace( f );
           __idGmshToFeel.insert( std::make_pair(eid,f.id()));
           eid ++;

     }
     }

     if ( partId > 0)
    {
    for( int i = 0; i < M_nx-1; ++i )
    {
           face_type f;
           //int eid =nelemG+(M_nx-1)*(2*partId-1)+i+1;
           f.setId(eid);
           f.setProcessIdInPartition( partId );
           f.setProcessId( final-1 );
           int ptid[2] = { (M_ny)*i+cx[partId],   // 0
                            (M_ny)*(i+1)+cx[partId], // 1
                         };

           for( int k = 0; k < 2; ++k )
           {
               f.setPoint( k, this->point( ptid[k]  ) );
           }
           //f.setMarker(1);
           f.setOnBoundary( true );
           this->addFace( f );
           __idGmshToFeel.insert( std::make_pair(eid,f.id()));
           std::cout << eid << ":" << f.id() << std::endl;
           //mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( nelemG+(M_nx-1)*(2*partId)+i+1, partId-1 )));
           eid ++;

     }
     }
*/


    toc("MeshStructured Edges", FLAGS_v>0);
//#endif
    tic();
    int eltid = 0;
    /*
     * Generates elements whithout ghosts
     */
    tic();
    for( int j = cx[partId]; j < cx[partId+1]; ++j )
    {
        ghosts.clear();
        if (j == cx[partId] && j != 0)
            ghosts.push_back(partId-1);
        else if (j == cx[partId+1]-1 && j != M_ny-2)
        {
            std::cout << " I'm proc " << partId << " and I have " << partId+1 << "as neigh at " << j << std::endl;
            ghosts.push_back(partId+1);
        }
        for( int i = 0; i < M_nx-1; ++i )
        {
            element_type e;
            int eid = (M_ny-1)*i+j; // GMSH Id
            //int eid =((M_ny-1)/procSize)*i+j%((M_ny-1)/procSize)+M_nx-1;
            e.setId( eid ); // eid == gmsh id
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId );
            //e.setMarker(1);
            //int l_j = j-cx[procId];
            //int l_i = i;
            //int l_eid = (cx[procId+1]-cx[procId])*l_i+l_j;
            //std::cout << eid << "\t" << partId << "\t" << partId << "\n";
            e.setNeighborPartitionIds( ghosts );

            // points definition
            int ptid[4] = { (M_ny)*(i+1)+j,   // 0
                            (M_ny)*(i+1)+j+1, // 1
                            (M_ny)*i+j+1,     // 2
                            (M_ny)*i+j        // 3
                          };
           /* std::cout << "eid = " << eid << std::endl;
            for(auto it : ptid)
                std::cout << it << "\t" << std::endl;
            std::cout << std::endl;
      */

            for( int k = 0; k < 4; ++k )
            {
                //this->points().modify( this->pointIterator( ptid[k] ), Feel::detail::UpdateProcessId(e.processId()) );
                e.setPoint( k, this->point( ptid[k]  ) );
            }
                //e.setPoint( k, this->point( ptid[k]  ) );
            //e.setOnBoundary( i==0 || j==0 || i==(M_nx-2) || j==(M_ny-2) );
            e.setOnBoundary( i==0 || i==(M_nx-2) || j%((M_ny-2)/procSize)==0 );
                // true -> id = global elt id
                // false -> id = local elt id
            this->addElement( e, true ); // e.id() is modified to Feel++ Id
            //std::cout << partId << "\t" << i << " : " << j << "\t" << e.id() << std::endl;
            __idGmshToFeel.insert( std::make_pair(eid,e.id())); // e.id() == Feel++ Id
            //__idGmshToFeel.insert( std::make_pair(eid,((M_ny-1)/procSize)*i+j%((M_ny-1)/procSize)+M_nx-1)); // e.id() == Feel++ Id
        }
    }
    toc("setElementInPart", FLAGS_v>0);
    tic();

#if 1
/* Ghost management */
    int j_inf =cx[partId];
    int j_sup =cx[partId+1];
    bool is_first = (partId == 0); /* proc 0 ? */
    bool is_last = (partId == (procSize-1));
    if(! is_last ) // look at the right
    {
        std::cout << "I'm proc " << partId << " and go in !is_last\n";
        /*
         * Indicate I am a ghost to neighborhood partition
         */
        ghosts.clear();
        ghosts.push_back(partId+1);
        for( int i = 0; i < M_nx; ++i )
        {
            int ptId = (M_ny)*i+j_sup+1;
            coords[0]=M_pixelsize*j_sup+1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                       i == 0 || i == M_nx-1  ); // Is on boundary ?

            pt.setProcessIdInPartition( partId );
            pt.setProcessId( invalid_rank_type_value );
            //pt.setMarker(1);
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        }
        for( int i = 0; i < M_nx-1; ++i )
        {
            element_type e;
            int eid = (M_ny-1)*i+j_sup;
            int eidF = i;
            //int eidF =(cx[partId]-cx[partId-1])*(M_nx-1)+i+M_nx-1 ;
            e.setId( eidF );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId+1 );
            //e.setMarker(1);
            //e.setIdInOtherPartitions(partId+1,i+M_nx-1);
            //e.setNeighborPartitionIds( ghosts );

            int ptid[4] = { (M_ny)*(i+1)+j_sup,   // 0
                            (M_ny)*(i+1)+j_sup+1, // 1
                            (M_ny)*i+j_sup+1,     // 2
                            (M_ny)*i+j_sup        // 3
                          };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );

            e.setOnBoundary( true );
            nbMsgToRecv[partId+1]+=1;
            //std::cout << partId << std::endl;
            //__idGmshToFeel.insert( std::make_pair(eid,eid));
            // true -> id = global elt id
            // false -> id = local elt id
            /*std::cout << "*******************************\n";
            std::cout << "id \t"    << e.id() << std::endl;
            std::cout << "pid\t"    << e.processId() << std::endl;
            std::cout << "pidinp\t" << e.pidInPartition() << std::endl;
            std::cout << "refd\t"   << e.refDim() << std::endl;
            std::cout << "np\t"     << e.nPoints() << std::endl;
            auto tete =  e.idInOthersPartitions();
            std::cout << tete.size() << std::endl;
            for (auto it : tete )
                std::cout << it.first << "\t" << it.second << std::endl;*/
            auto el=this->addElement( e, false );
            /*std::cout <<"id \t"    <<  toto.id() << std::endl;
            std::cout <<"pid\t"    <<  toto.processId() << std::endl;
            std::cout <<"pidinp\t" <<  toto.pidInPartition() << std::endl;
            std::cout <<"refd\t"   <<  toto.refDim() << std::endl;
            std::cout <<"np\t"     <<  toto.nPoints() << std::endl;
            auto titi =  toto.idInOthersPartitions();
            std::cout << titi.size() << std::endl;
            for (auto it : titi )
                std::cout << it.first << "\t" << it.second << std::endl;
            */
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple( el.id(), partId+1 )));
            /*mapGhostElt.insert( std::make_pair( eid,boost::make_tuple(
            (nx-1)*cx[partId+1]+i*(cx[partId+1]-cx[partId]),
            partId+1 ) ) );*/
            __idGmshToFeel.insert( std::make_pair(eid,el.id()));
            //__idGmshToFeel.insert( std::make_pair(eid,(nx-1)*cx[partId+1]+i*(cx[partId+1]-cx[partId])));


        }
    } // is_last
    ghosts.clear();
    if(! is_first) // look at the left
    {
        ghosts.push_back(partId-1);
        for( int i = 0; i < M_nx; ++i )
        {
            int ptId = (M_ny)*i+j_inf-1;
            coords[0]=M_pixelsize*j_inf-1;
            coords[1]=M_pixelsize*i;
            point_type pt( ptId, coords,
                    i == 0 || i == M_nx-1  ); // Is on boundary ?
            pt.setProcessIdInPartition( partId );
            pt.setProcessId( invalid_rank_type_value );
            //pt.setMarker(1);
            //pt.setNeighborPartitionIds( ghosts );

            this->addPoint( pt );
        }

        for( int i = 0; i < M_nx-1; ++i )
        {
            element_type e;
            int eid = (M_ny-1)*i+j_inf-1;
            int eidF =(cx[partId]-cx[partId-1]-1)*(M_nx-1)+i ;
            //int eidF = i ;
            e.setId( eidF );
            e.setProcessIdInPartition( partId );
            e.setProcessId( partId-1 );
            //e.setMarker(1);
            //e.setIdInOtherPartitions(partId-1,(cx[partId]-cx[partId-1])*(M_nx-1)+i);
            //e.setNeighborPartitionIds( ghosts );


            int ptid[4] = { (M_ny)*(i+1)+j_inf-1,   // 0
                            (M_ny)*(i+1)+j_inf, // 1
                            (M_ny)*i+j_inf,     // 2
                            (M_ny)*i+j_inf-1        // 3
            };

            for( int k = 0; k < 4; ++k )
                e.setPoint( k, this->point( ptid[k]  ) );
            e.setOnBoundary( true );
            nbMsgToRecv[partId-1]+=1;
            // true -> id = global elt id
            // false -> id = local elt id
            /*std::cout << "*******************************\n";
            std::cout << "id \t"    << e.id() << std::endl;
            std::cout << "pid\t"    << e.processId() << std::endl;
            std::cout << "pidinp\t" << e.pidInPartition() << std::endl;
            std::cout << "refd\t"   << e.refDim() << std::endl;
            std::cout << "np\t"     << e.nPoints() << std::endl;
            auto tete =  e.idInOthersPartitions();
            std::cout << tete.size() << std::endl;
            for (auto it : tete )
                std::cout << it.first << "\t" << it.second << std::endl;*/
            //this->addElement( e, true );
            auto el=this->addElement( e, false );
            /*std::cout <<"id \t"    <<  toto.id() << std::endl;
            std::cout <<"pid\t"    <<  toto.processId() << std::endl;
            std::cout <<"pidinp\t" <<  toto.pidInPartition() << std::endl;
            std::cout <<"refd\t"   <<  toto.refDim() << std::endl;
            std::cout <<"np\t"     <<  toto.nPoints() << std::endl;
            auto titi =  toto.idInOthersPartitions();
            std::cout << titi.size() << std::endl;
            for (auto it : titi )
                std::cout << it.first << "\t" << it.second << std::endl;
            */
            mapGhostElt.insert( std::make_pair( eid,boost::make_tuple(el.id(),partId-1)));

            /*mapGhostElt.insert( std::make_pair( eid,boost::make_tuple(
            (nx-1)*cx[partId-1]+(i+1)*(cx[partId]-cx[partId-1])-1,
            partId-1 ) ) );*/
            __idGmshToFeel.insert( std::make_pair(eid,el.id()));
            //__idGmshToFeel.insert( std::make_pair(eid,(nx-1)*cx[partId-1]+(i+1)*(cx[partId]-cx[partId-1])-1));

        }
    } // is_first

    for (auto it=mapGhostElt.begin(); it!=mapGhostElt.end(); ++it)
        std::cout << "test : " << it->first << " ( " << it->second.get<0>()<< "," << it->second.get<1>() << " )" << std::endl;
    this->worldComm().barrier();
    // this->addMarkerName("omega",1,2);
    updateGhostCellInfoByUsingNonBlockingComm(this,__idGmshToFeel,mapGhostElt,nbMsgToRecv);
    toc("Ghost Management", FLAGS_v>0);
    toc("MeshStructured Elements", FLAGS_v>0);
#endif
}
#endif

//template<typename MeshType>
void MeshStructured::updateGhostCellInfoByUsingNonBlockingComm( std::unordered_map<int, int> const& idStructuredMeshToFeelMesh,
                                                                std::unordered_map<int, boost::tuple<int, rank_type>> const& mapGhostElt,
                                                                std::vector<int> const& nbMsgToRecv )
{
    DVLOG( 1 ) << "updateGhostCellInfoNonBlockingComm : start on rank " << this->worldComm().localRank() << "\n";

    const int nProc = this->worldComm().localSize();

    //std::cout << nProc << std::endl;

    //-----------------------------------------------------------//
    // compute size of container to send
    std::unordered_map<rank_type, int> nDataInVecToSend;
    auto it_map = mapGhostElt.begin();
    auto const en_map = mapGhostElt.end();
    for ( ; it_map != en_map; ++it_map )
    {
        const rank_type idProc = it_map->second.template get<1>();
        if ( nDataInVecToSend.find( idProc ) == nDataInVecToSend.end() )
            nDataInVecToSend[idProc] = 0;
        nDataInVecToSend[idProc]++;
    }
    //-----------------------------------------------------------//
    // init and resize the container to send
    std::unordered_map<rank_type, std::vector<int>> dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend != enNDataInVecToSend; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //-----------------------------------------------------------//
    // prepare container to send
    std::unordered_map<rank_type, std::unordered_map<int, int>> memoryMsgToSend;
    std::unordered_map<rank_type, int> nDataInVecToSendBis;
    it_map = mapGhostElt.begin();
    for ( ; it_map != en_map; ++it_map )
    {
        const int idGmsh = it_map->first;
        const int idFeel = it_map->second.template get<0>();
        const rank_type idProc = it_map->second.template get<1>();

        if ( nDataInVecToSendBis.find( idProc ) == nDataInVecToSendBis.end() )
            nDataInVecToSendBis[idProc] = 0;
        // save request
        memoryMsgToSend[idProc][nDataInVecToSendBis[idProc]] = idFeel;
        // update container
        dataToSend[idProc][nDataInVecToSendBis[idProc]] = idGmsh;
        // update counter
        nDataInVecToSendBis[idProc]++;
        // std::cout << idProc << std::endl;
    }
    //-----------------------------------------------------------//
    // counter of request
    int nbRequest = 0;
    for ( rank_type proc = 0; proc < nProc; ++proc )
    {
        if ( dataToSend.find( proc ) != dataToSend.end() )
            ++nbRequest;
        if ( nbMsgToRecv[proc] > 0 )
            ++nbRequest;
    }
    if ( nbRequest == 0 ) return;

    mpi::request* reqs = new mpi::request[nbRequest];
    int cptRequest = 0;
    //-----------------------------------------------------------//
    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend != enDataToSend; ++itDataToSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToSend->first, 0, itDataToSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // first recv
    std::unordered_map<rank_type, std::vector<int>> dataToRecv;
    for ( rank_type proc = 0; proc < nProc; ++proc )
    {
        if ( nbMsgToRecv[proc] > 0 )
        {
            reqs[cptRequest] = this->worldComm().localComm().irecv( proc, 0, dataToRecv[proc] );
            ++cptRequest;
        }
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    //-----------------------------------------------------------//
    // build the container to ReSend
    std::unordered_map<rank_type, std::vector<int>> dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv != enDataRecv; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[idProc].resize( nDataRecv );
        //store the idFeel corresponding
        for ( int k = 0; k < nDataRecv; ++k )
        {
            dataToReSend[idProc][k] = idStructuredMeshToFeelMesh.find( itDataRecv->second[k] )->second;
        }
    }
    //-----------------------------------------------------------//
    // send respond to the request
    cptRequest = 0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend != enDataToReSend; ++itDataToReSend )
    {
        reqs[cptRequest] = this->worldComm().localComm().isend( itDataToReSend->first, 0, itDataToReSend->second );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // recv the initial request
    std::unordered_map<rank_type, std::vector<int>> finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend != enDataToSend; ++itDataToSend )
    {
        const rank_type idProc = itDataToSend->first;
        reqs[cptRequest] = this->worldComm().localComm().irecv( idProc, 0, finalDataToRecv[idProc] );
        ++cptRequest;
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    // delete reqs because finish comm
    delete[] reqs;
    //-----------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();

    for ( ; itFinalDataToRecv != enFinalDataToRecv; ++itFinalDataToRecv )
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        //std::cout << idProc << ":" << nDataRecv << std::endl;
        for ( int k = 0; k < nDataRecv; ++k )
        {
            /* std::cout << "I want element " << memoryMsgToSend[idProc][k] << ": " << idProc << std::endl;*/
            auto& eltToUpdate = this->elementIterator( memoryMsgToSend[idProc][k] /*,idProc*/ )->second;
#if 0
            std::cout << "k = " << k << std::endl;
            std::cout << "itFinalDataToRecv->second[k]  " << itFinalDataToRecv->second[k]  << std::endl;
            std::cout << "eltToUpdate->id()             " << eltToUpdate->id()             << std::endl;
            std::cout << "eltToUpdate->processId()      " << eltToUpdate->processId()      << std::endl;
            std::cout << "eltToUpdate->pidInPartition() " << eltToUpdate->pidInPartition() << std::endl;
            std::cout << "eltToUpdate->refDim()         " << eltToUpdate->refDim()         << std::endl;
            std::cout << "eltToUpdate->nPoints()        " << eltToUpdate->nPoints()        << std::endl;
#endif
            eltToUpdate.setIdInOtherPartitions( idProc, itFinalDataToRecv->second[k] );
        }
    }
    //-----------------------------------------------------------//
    DVLOG( 1 ) << "updateGhostCellInfoNonBlockingComm : finish on rank " << this->worldComm().localRank() << "\n";
}

} // namespace Feel
