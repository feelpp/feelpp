/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2012-01-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file vectorblock.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2012-01-18
 */

#include <feel/feelalg/vectorblock.hpp>



namespace Feel
{


template <typename T>
VectorBlockBase<T>::VectorBlockBase( vf::BlocksBase<vector_ptrtype> const & blockVec,
                                     backend_type &backend,
                                     bool copy_values )
    :
    M_vec()
{
    auto nRow = blockVec.nRow();

#if 0
    size_type _size = 0;
    for ( uint i=0; i<nRow; ++i )
        _size += blockVec( i,0 )->size();

    M_vec = backend.newVector( _size,_size );
#else

    DataMap dm(blockVec(0,0)->map().worldComm());
    const int myrank = dm.worldComm().globalRank();
    const int worldsize = dm.worldComm().globalSize();
    for (int proc = 0 ; proc < worldsize ; ++proc)
    {
        size_type firstDofGlobalCluster=0;
        for ( int p=0; p<proc; ++p )
            for ( uint16_type i=0 ; i<nRow; ++i )
                firstDofGlobalCluster += blockVec(i,0)->map().nLocalDofWithoutGhost( p );
        dm.setFirstDofGlobalCluster( proc, firstDofGlobalCluster );

        size_type sizeWithoutGhost=0, sizeWithGhost=0, sizeGlobalCluster=0;
        for ( uint16_type i=0 ; i<nRow; ++i)
        {
            sizeWithoutGhost += blockVec(i,0)->map().nLocalDofWithoutGhost( proc );
            sizeWithGhost += blockVec(i,0)->map().nLocalDofWithGhost( proc );
            sizeGlobalCluster += blockVec(i,0)->map().nDof();
        }
        dm.setNLocalDofWithoutGhost( proc, sizeWithoutGhost );
        dm.setNLocalDofWithGhost( proc, sizeWithGhost );
        dm.setFirstDof( proc, 0 );
        dm.setLastDof( proc, (sizeWithGhost == 0)?0:sizeWithGhost-1 );
        dm.setLastDofGlobalCluster(proc,  (sizeWithoutGhost ==0)? firstDofGlobalCluster : ( firstDofGlobalCluster +sizeWithoutGhost-1 ));
        if ( proc==myrank )
            dm.setNDof( sizeGlobalCluster );
    }

    dm.resizeMapGlobalProcessToGlobalCluster( dm.nLocalDofWithGhost(myrank) );
    dm.resizeMapGlobalClusterToGlobalProcess( dm.nLocalDofWithoutGhost(myrank) );
    const size_type firstDofGC = dm.firstDofGlobalCluster(myrank);
    size_type start_i = firstDofGC;
    size_type nLocalDofStart = dm.firstDof();
    for ( uint16_type i=0 ; i<nRow; ++i)
    {
        const size_type firstBlockDofGC =  blockVec(i,0)->map().firstDofGlobalCluster(myrank);
        for (size_type gdof = blockVec(i,0)->map().firstDof(myrank) ; gdof < blockVec(i,0)->map().nLocalDofWithGhost(myrank) ; ++gdof )
        {
            const size_type localDof = nLocalDofStart+gdof;
            size_type gdofGC = blockVec(i,0)->map().mapGlobalProcessToGlobalCluster(gdof);
            if ( blockVec(i,0)->map().dofGlobalClusterIsOnProc( gdofGC ) )
            {
                const size_type globalDof = start_i+(gdofGC-firstBlockDofGC);
                dm.setMapGlobalProcessToGlobalCluster( localDof, globalDof );
                dm.setMapGlobalClusterToGlobalProcess( globalDof-firstDofGC ,localDof );
            }
            else
            {
                const int realproc = blockVec(i,0)->map().procOnGlobalCluster(gdofGC);
                size_type nDofStart=dm.firstDofGlobalCluster(realproc);
                for ( uint16_type k=0; k<i; ++k )
                    nDofStart += blockVec(k,0)->map().nLocalDofWithoutGhost( realproc );
                const size_type globDof = nDofStart+(gdofGC- blockVec(i,0)->map().firstDofGlobalCluster(realproc));
                dm.setMapGlobalProcessToGlobalCluster( localDof, globDof );
            }
        }
        nLocalDofStart += blockVec(i,0)->map().nLocalDofWithGhost( myrank );

        start_i += blockVec(i,0)->map().nLocalDofWithoutGhost( myrank );
    }

    //dm.showMeMapGlobalProcessToGlobalCluster();

    M_vec = backend.newVector( dm );
#endif
    M_vec->zero();

    if ( copy_values )
    {
        size_type start_i=0;

        for ( uint i=0; i<nRow; ++i )
        {
            this->updateBlockVec( blockVec( i,0 ),start_i );
            start_i += blockVec( i,0 )->size();
        }
    }
}

template <typename T>
void
VectorBlockBase<T>::updateBlockVec( vector_ptrtype const& m, size_type start_i )
{
    const size_type size = m->localSize();

    for ( uint i=0; i<size; ++i )
        M_vec->set( start_i+i,m->operator()( i ) );
}


template class VectorBlockBase<double>;

} // Feel


