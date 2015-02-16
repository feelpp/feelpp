/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-10-29

  Copyright (C) 2007-2011 Universite Joseph Fourier (Grenoble I)

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
   \file graphcsr.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-10-29
 */
#include <boost/timer.hpp>

#include <feel/feelcore/application.hpp>
#include <feel/feelalg/graphcsr.hpp>

namespace Feel
{
GraphCSR::GraphCSR( size_type n,
                    size_type first_row_entry_on_proc,
                    size_type last_row_entry_on_proc,
                    size_type first_col_entry_on_proc,
                    size_type last_col_entry_on_proc,
                    WorldComm const& worldcomm )
    :
    M_is_closed( false ),
    M_worldComm( worldcomm ),
    M_first_row_entry_on_proc( this->worldComm().globalSize(), first_row_entry_on_proc ),
    M_last_row_entry_on_proc( this->worldComm().globalSize(),last_row_entry_on_proc ),
    M_first_col_entry_on_proc( this->worldComm().globalSize(),first_col_entry_on_proc ),
    M_last_col_entry_on_proc( this->worldComm().globalSize(),last_col_entry_on_proc ),
    M_max_nnz( 0 ),
    M_n_total_nz( n, 0 ),
    M_n_nz( n, 0 ),
    M_n_oz( n, 0 ),
    M_storage(),
    M_mapRow( new DataMap(worldcomm) ),
    M_mapCol( new DataMap(worldcomm) )
{
    const int myrank = this->worldComm().globalRank();
    const int worldsize = this->worldComm().globalSize();
    const size_type _size1 = last_row_entry_on_proc - first_row_entry_on_proc+1;
    const size_type _size2 = last_col_entry_on_proc - first_col_entry_on_proc+1;

    for (int proc = 0 ; proc < worldsize ; ++proc)
    {
        M_mapRow->setFirstDofGlobalCluster( proc,first_row_entry_on_proc );
        M_mapRow->setLastDofGlobalCluster(proc, last_row_entry_on_proc );
        M_mapCol->setFirstDofGlobalCluster( proc,first_col_entry_on_proc );
        M_mapCol->setLastDofGlobalCluster(proc, last_col_entry_on_proc );


        M_mapRow->setNLocalDofWithoutGhost( proc, _size1 );
        M_mapRow->setNLocalDofWithGhost( proc, _size1 );
        M_mapCol->setNLocalDofWithoutGhost( proc, _size2 );
        M_mapCol->setNLocalDofWithGhost( proc, _size2 );

        M_mapRow->setFirstDof( proc, 0 );
        M_mapCol->setFirstDof( proc, 0 );
        if (_size2==0)
            M_mapCol->setLastDof( proc, 0 );
        else
            M_mapCol->setLastDof( proc, _size2-1 );

        if ( _size1==0 )
            M_mapRow->setLastDof( proc, 0 );
        else
            M_mapRow->setLastDof( proc, _size1-1 );

        if ( proc==myrank )
        {
            M_mapRow->setNDof( _size1 );
            M_mapCol->setNDof( _size2 );
        }

    }
    M_mapRow->resizeMapGlobalProcessToGlobalCluster( M_mapRow->nLocalDofWithGhost(myrank ) );
    M_mapRow->resizeMapGlobalClusterToGlobalProcess( M_mapRow->nLocalDofWithoutGhost(myrank) );
    M_mapCol->resizeMapGlobalProcessToGlobalCluster( M_mapCol->nLocalDofWithGhost(myrank) );
    M_mapCol->resizeMapGlobalClusterToGlobalProcess( M_mapCol->nLocalDofWithoutGhost(myrank) );

    for ( size_type i=0;i<_size1;++i )
    {
        M_mapRow->setMapGlobalProcessToGlobalCluster( i,first_row_entry_on_proc+i );
        M_mapRow->setMapGlobalClusterToGlobalProcess( first_row_entry_on_proc+i,i );
    }

    for( size_type j=0;j<_size2;++j )
    {
        M_mapCol->setMapGlobalProcessToGlobalCluster( j,first_col_entry_on_proc+j );
        M_mapCol->setMapGlobalClusterToGlobalProcess( first_col_entry_on_proc+j,j );
    }

}


GraphCSR::GraphCSR( boost::shared_ptr<DataMap> const& mapRow,
                    boost::shared_ptr<DataMap> const& mapCol )
    :
    M_is_closed( false ),
    M_worldComm( mapRow->worldComm() ),
    M_first_row_entry_on_proc( mapRow->firstDofGlobalClusterWorld() ),
    M_last_row_entry_on_proc( mapRow->lastDofGlobalClusterWorld() ),
    M_first_col_entry_on_proc( mapCol->firstDofGlobalClusterWorld() ),
    M_last_col_entry_on_proc( mapCol->lastDofGlobalClusterWorld() ),
    M_max_nnz( 0 ),
    M_n_total_nz( 0 ),
    M_n_nz( 0 ),
    M_n_oz( 0 ),
    M_storage(),
    M_mapRow(mapRow),
    M_mapCol(mapCol)
{}

GraphCSR::GraphCSR( DataMap const& mapRow,
                    DataMap const& mapCol )
    :
    M_is_closed( false ),
    M_worldComm( mapRow.worldComm() ),
    M_first_row_entry_on_proc( mapRow.firstDofGlobalClusterWorld() ),
    M_last_row_entry_on_proc( mapRow.lastDofGlobalClusterWorld() ),
    M_first_col_entry_on_proc( mapCol.firstDofGlobalClusterWorld() ),
    M_last_col_entry_on_proc( mapCol.lastDofGlobalClusterWorld() ),
    M_max_nnz( 0 ),
    M_n_total_nz( 0 ),
    M_n_nz( 0 ),
    M_n_oz( 0 ),
    M_storage(),
    M_mapRow( new DataMap(mapRow) ),
    M_mapCol( new DataMap(mapCol) )
{}




GraphCSR::GraphCSR( vf::BlocksBase<self_ptrtype> const & blockSet,
                    bool diagIsNonZero, bool close )
    :
    M_is_closed( false ),
    M_worldComm( blockSet(0,0)->worldComm()/* Environment::worldComm()*/ ),
    M_first_row_entry_on_proc( this->worldComm().globalSize(),0 ),
    M_last_row_entry_on_proc( this->worldComm().globalSize(),0 ),
    M_first_col_entry_on_proc( this->worldComm().globalSize(),0 ),
    M_last_col_entry_on_proc( this->worldComm().globalSize(),0 ),
    M_max_nnz( 0 ),
    M_n_total_nz( /*n*/0, 0 ),
    M_n_nz( /*n*/0, 0 ),
    M_n_oz( /*n*/0, 0 ),
    M_storage(),
    M_mapRow( new DataMap(this->worldComm()) ),
    M_mapCol( new DataMap(this->worldComm()) )
{
    DVLOG(2) << " GraphCSR constructor from block graph\n";

    auto const nRow = blockSet.nRow();
    auto const nCol = blockSet.nCol();

    const int myrank = this->worldComm().globalRank();
    const int worldsize = this->worldComm().globalSize();

    this->updateDataMap( blockSet );

    DVLOG(2) << "myrank " << myrank
             << " this->firstRowEntryOnProc() " << this->firstRowEntryOnProc()
             << " this->firstColEntryOnProc() " << this->firstColEntryOnProc()
             << " this->lastRowEntryOnProc() " << this->lastRowEntryOnProc()
             << " this->lastColEntryOnProc() " << this->lastColEntryOnProc() << "\n";

    size_type start_i=this->firstRowEntryOnProc();
    size_type start_j=this->firstColEntryOnProc();
    for ( int i=0; i<nRow; ++i )
    {
        start_j=this->firstColEntryOnProc();

        for ( int j=0; j<nCol; ++j )
        {
            if ( blockSet(i,j)->empty() )
            {
                start_j += blockSet(i,j)->mapCol().nLocalDofWithoutGhost( myrank );
                continue;
            }

            //blockSet(i,j)->close();
            if (this->worldComm().globalSize()==1)
                this->mergeBlockGraph( blockSet(i,j),start_i,start_j );
            else
                this->mergeBlockGraphMPI( blockSet(i,j),blockSet,i,j,start_i,start_j );

            start_j += blockSet(i,j)->mapCol().nLocalDofWithoutGhost( myrank );
        }

        start_i += blockSet(i,0)->mapRow().nLocalDofWithoutGhost( myrank );
    }

    //M_mapCol=M_mapRow;
    //M_mapCol.setMapGlobalClusterToGlobalProcess(M_mapRow.mapGlobalClusterToGlobalProcess());
    //M_mapCol.setMapGlobalProcessToGlobalCluster(M_mapRow.mapGlobalProcessToGlobalCluster());

    if ( diagIsNonZero ) this->addMissingZeroEntriesDiagonal();

    if ( close ) this->close();
}



GraphCSR::GraphCSR( GraphCSR const & g )
    :
    M_is_closed( g.M_is_closed ),
    M_worldComm( g.M_worldComm ),
    M_first_row_entry_on_proc( g.M_first_row_entry_on_proc ),
    M_last_row_entry_on_proc( g.M_last_row_entry_on_proc ),
    M_first_col_entry_on_proc( g.M_first_col_entry_on_proc ),
    M_last_col_entry_on_proc( g.M_last_col_entry_on_proc ),
    M_max_nnz( g.M_max_nnz ),
    M_n_total_nz( g.M_n_total_nz ),
    M_n_nz( g.M_n_nz ),
    M_n_oz( g.M_n_oz ),
    M_storage( g.M_storage ),
    M_graphT( g.M_graphT ),
    M_mapRow(g.M_mapRow),
    M_mapCol(g.M_mapCol)
{
}

GraphCSR::~GraphCSR()
{
    M_storage.clear();
}

GraphCSR&
GraphCSR::operator=( GraphCSR const& g )
{
    if ( this != &g )
    {
        M_first_row_entry_on_proc = g.M_first_row_entry_on_proc;
        M_last_row_entry_on_proc = g.M_last_row_entry_on_proc;
        M_first_col_entry_on_proc = g.M_first_col_entry_on_proc;
        M_last_col_entry_on_proc = g.M_last_col_entry_on_proc;
        M_max_nnz = g.M_max_nnz;
        M_n_total_nz= g.M_n_total_nz;
        M_n_nz = g.M_n_nz;
        M_n_oz = g.M_n_oz;
        M_storage = g.M_storage;
        M_graphT = g.M_graphT;
        M_is_closed = g.M_is_closed;
        M_worldComm = g.M_worldComm;
        M_mapRow = g.M_mapRow;
        M_mapCol = g.M_mapCol;
    }

    return *this;
}



void
GraphCSR::updateDataMap( vf::BlocksBase<self_ptrtype> const & blockSet )
{
    const int myrank = this->worldComm().globalRank();
    const int worldsize = this->worldComm().globalSize();
    auto const nRow = blockSet.nRow();
    auto const nCol = blockSet.nCol();

    for (int proc = 0 ; proc < worldsize ; ++proc)
    {
        //-------------------------------------------------------------//
        size_type size1=0, size2=0;
        for ( int p=0; p<proc; ++p )
        {
            for ( uint i=0; i<nCol; ++i )
                size2 += blockSet(0,i)->mapCol().nLocalDofWithoutGhost( p );
            for ( uint i=0; i<nRow; ++i )
                size1 += blockSet(i,0)->mapRow().nLocalDofWithoutGhost( p );
        }
        M_first_col_entry_on_proc[proc] = size2;
        M_first_row_entry_on_proc[proc] = size1;
        M_mapCol->setFirstDofGlobalCluster( proc, size2 );
        M_mapRow->setFirstDofGlobalCluster( proc, size1 );
        //-------------------------------------------------------------//
        M_mapCol->setFirstDof( proc, 0 );
        M_mapRow->setFirstDof( proc, 0 );
        size_type _size2WithoutGhost=0,_size2WithGhost=0,_globalClusterSize2=0;
        for ( uint i=0; i<nCol; ++i )
        {
            _size2WithoutGhost += blockSet(0,i)->mapCol().nLocalDofWithoutGhost( proc );
            _size2WithGhost += blockSet(0,i)->mapCol().nLocalDofWithGhost( proc );
            _globalClusterSize2 += blockSet(0,i)->mapCol().nDof();
        }
        size_type _size1WithoutGhost=0,_size1WithGhost=0,_globalClusterSize1=0;
        for ( uint i=0; i<nRow; ++i )
        {
            _size1WithoutGhost += blockSet(i,0)->mapRow().nLocalDofWithoutGhost( proc );
            _size1WithGhost += blockSet(i,0)->mapRow().nLocalDofWithGhost( proc );
            _globalClusterSize1 += blockSet(i,0)->mapRow().nDof();
        }
        //-------------------------------------------------------------//
        M_mapCol->setNLocalDofWithoutGhost( proc, _size2WithoutGhost );
        M_mapCol->setNLocalDofWithGhost( proc, _size2WithGhost );
        if ( _size2WithGhost==0 )
            M_mapCol->setLastDof( proc, 0 );
        else
            M_mapCol->setLastDof( proc, _size2WithGhost-1 );
        if ( proc==myrank )
            M_mapCol->setNDof( _globalClusterSize2 );
        //-------------------------------------------------------------//
        M_mapRow->setNLocalDofWithoutGhost( proc, _size1WithoutGhost );
        M_mapRow->setNLocalDofWithGhost( proc, _size1WithGhost );
        if ( _size1WithGhost==0 )
            M_mapRow->setLastDof( proc, 0 );
        else
            M_mapRow->setLastDof( proc, _size1WithGhost-1 );

        if ( proc==myrank )
            M_mapRow->setNDof( _globalClusterSize1 );
        //-------------------------------------------------------------//
        if ( M_mapRow->nLocalDofWithoutGhost(proc) == 0 )
            M_mapRow->setLastDofGlobalCluster(proc,  M_mapRow->firstDofGlobalCluster(proc) );
        else
            M_mapRow->setLastDofGlobalCluster(proc,  M_mapRow->firstDofGlobalCluster(proc) + M_mapRow->nLocalDofWithoutGhost(proc)-1 );

        if ( M_mapCol->nLocalDofWithoutGhost(proc) == 0 )
            M_mapCol->setLastDofGlobalCluster(proc,  M_mapCol->firstDofGlobalCluster(proc) );
        else
            M_mapCol->setLastDofGlobalCluster(proc,  M_mapCol->firstDofGlobalCluster(proc) + M_mapCol->nLocalDofWithoutGhost(proc)-1 );
        //-------------------------------------------------------------//
        M_last_row_entry_on_proc[proc] = M_mapRow->lastDofGlobalCluster(proc);
        M_last_col_entry_on_proc[proc] = M_mapCol->lastDofGlobalCluster(proc);
    }

    M_mapRow->resizeMapGlobalProcessToGlobalCluster( M_mapRow->nLocalDofWithGhost(myrank) );
    M_mapRow->resizeMapGlobalClusterToGlobalProcess( M_mapRow->nLocalDofWithoutGhost(myrank) );
    M_mapCol->resizeMapGlobalProcessToGlobalCluster( M_mapCol->nLocalDofWithGhost(myrank) );
    M_mapCol->resizeMapGlobalClusterToGlobalProcess( M_mapCol->nLocalDofWithoutGhost(myrank) );

    // update DataMapRow : setMapGlobalProcessToGlobalCluster setMapGlobalClusterToGlobalProcess
    const size_type firstDofGCRow = M_mapRow->firstDofGlobalCluster(myrank);
    size_type start_ii = firstDofGCRow;
    size_type nLocalDofStartRow = M_mapRow->firstDof();
    for ( uint16_type i=0 ; i<nRow; ++i)
    {
        auto const& mapRowOnBlock = blockSet(i,0)->mapRow();

        for ( rank_type procIdNeigh : mapRowOnBlock.neighborSubdomains() )
            M_mapRow->addNeighborSubdomain(procIdNeigh);

        const size_type firstBlockDofGC =  mapRowOnBlock.firstDofGlobalCluster(myrank);

        for (size_type gdof = mapRowOnBlock.firstDof(myrank) ; gdof < mapRowOnBlock.nLocalDofWithGhost(myrank) ; ++gdof )
        {
            const size_type localDofRow = nLocalDofStartRow+gdof;
            const size_type gdofGC = mapRowOnBlock.mapGlobalProcessToGlobalCluster(gdof);
            if ( mapRowOnBlock.dofGlobalClusterIsOnProc( gdofGC ) )
            {
                const size_type globalDofRow = start_ii+(gdofGC-firstBlockDofGC);
                M_mapRow->setMapGlobalProcessToGlobalCluster( localDofRow, globalDofRow );
                M_mapRow->setMapGlobalClusterToGlobalProcess( globalDofRow-firstDofGCRow,localDofRow );
            }
            else
            {
                const int realproc = mapRowOnBlock.procOnGlobalCluster(gdofGC);
                const size_type nDofStartRow = M_mapRow->firstDofGlobalCluster(realproc) + this->nLocalDofWithoutGhostOnProcStartRow( blockSet, realproc, i, 0);
                const size_type globalDofRow = nDofStartRow+(gdofGC- mapRowOnBlock.firstDofGlobalCluster(realproc));
                M_mapRow->setMapGlobalProcessToGlobalCluster( localDofRow, globalDofRow );
            }
         }

        for ( auto const& activeDofShared : mapRowOnBlock.activeDofSharedOnCluster() )
        {
            M_mapRow->setActiveDofSharedOnCluster( nLocalDofStartRow + activeDofShared.first, activeDofShared.second );
        }

        nLocalDofStartRow += mapRowOnBlock.nLocalDofWithGhost( myrank );
        start_ii += mapRowOnBlock.nLocalDofWithoutGhost( myrank );
    }

    // update DataMapCol : setMapGlobalProcessToGlobalCluster setMapGlobalClusterToGlobalProcess
    const size_type firstDofGCCol = M_mapCol->firstDofGlobalCluster(myrank);
    size_type start_jj = firstDofGCCol;
    size_type nLocalDofStartCol = M_mapCol->firstDof();
    for ( uint16_type j=0 ; j<nCol; ++j)
    {
        auto const& mapColOnBlock = blockSet(0,j)->mapCol();

        for ( rank_type procIdNeigh : mapColOnBlock.neighborSubdomains() )
            M_mapCol->addNeighborSubdomain(procIdNeigh);

        const size_type firstBlockDofGC =  mapColOnBlock.firstDofGlobalCluster( myrank );

        for (size_type gdof = mapColOnBlock.firstDof(myrank) ; gdof < mapColOnBlock.nLocalDofWithGhost(myrank) ; ++gdof )
        {
            const size_type localDofCol = nLocalDofStartCol+gdof;
            const size_type gdofGC = mapColOnBlock.mapGlobalProcessToGlobalCluster(gdof);
            if ( mapColOnBlock.dofGlobalClusterIsOnProc( gdofGC ) )
            {
                const size_type globalDofCol = start_jj+(gdofGC-firstBlockDofGC);
                M_mapCol->setMapGlobalProcessToGlobalCluster( localDofCol,globalDofCol );
                M_mapCol->setMapGlobalClusterToGlobalProcess( globalDofCol-firstDofGCCol,localDofCol );
            }
            else
            {
                const int realproc = mapColOnBlock.procOnGlobalCluster(gdofGC);
                const size_type nDofStartCol = M_mapCol->firstDofGlobalCluster(realproc) + this->nLocalDofWithoutGhostOnProcStartCol( blockSet, realproc, 0, j);
                const size_type globalDofCol = nDofStartCol+(gdofGC- mapColOnBlock.firstDofGlobalCluster(realproc));
                M_mapCol->setMapGlobalProcessToGlobalCluster( localDofCol, globalDofCol );
            }
         }

        for ( auto const& activeDofShared : mapColOnBlock.activeDofSharedOnCluster() )
        {
            M_mapCol->setActiveDofSharedOnCluster( nLocalDofStartCol + activeDofShared.first, activeDofShared.second );
        }

        nLocalDofStartCol += mapColOnBlock.nLocalDofWithGhost( myrank );
        start_jj += mapColOnBlock.nLocalDofWithoutGhost( myrank );
    }

    // index split ( only do for row )
    bool computeIndexSplit = true;
    if ( computeIndexSplit )
    {
        //const uint16_type nRow = blockgraph.nRow();
        boost::shared_ptr<IndexSplit> indexSplit( new IndexSplit() );
        const size_type firstDofGC = this->mapRow().firstDofGlobalCluster();
        for ( uint16_type i=0; i<nRow; ++i )
        {
            indexSplit->addSplit( firstDofGC, blockSet(i,0)->mapRow().indexSplit() );
        }
        //indexSplit->showMe();
        this->mapRowPtr()->setIndexSplit( indexSplit );
    }



}

void
GraphCSR::mergeBlockGraph( self_ptrtype const& g,
                           size_type start_i, size_type start_j )
{

    auto it = g->begin();
    auto en = g->end();
    for ( ; it != en; ++it )
    {
        int theglobalrow = start_i + it->first;
        row_type & row = this->row( theglobalrow );

        row.get<0>() = it->second.get<0>();//rank

        int thelocalrow = start_i + it->second.get<1>();
        row.get<1>() = thelocalrow;

        auto nbDof = it->second.get<2>().size();

        if ( nbDof>0 )
        {
            // Get the row of the sparsity pattern
#if 0
            std::vector<size_type> ivec(  it->second.template get<2>().begin(),  it->second.template get<2>().end() );
            std::for_each( ivec.begin(), ivec.end(), boost::phoenix::arg_names::arg1 += start_j );
#else
            std::vector<size_type> ivec( nbDof );
            auto itDof=it->second.get<2>().begin();

            for ( int i=0; i<( int )nbDof; ++i,++itDof )
                ivec[i]=*itDof+start_j;
#endif
            //std::set<size_type> iout( ivec.size()+ M_graph->row(theglobalrow).template get<2>().size() );
            //std::set<size_type> iout( ivec.begin(), ivec.end() );
            //iout.insert( globGraph->row(theglobalrow).template get<2>().begin(),
            //             globGraph->row(theglobalrow).template get<2>().end() );
            //globGraph->row(theglobalrow).template get<2>() = iout;
            //globGraph->row(theglobalrow).template get<2>().insert(ivec.begin(), ivec.end());
            row.get<2>().insert( ivec.begin(), ivec.end() );
        }
    }

}


void
GraphCSR::mergeBlockGraphMPI( self_ptrtype const& g,vf::BlocksBase<self_ptrtype> const & blockSet, int i, int j,
                              size_type start_i, size_type start_j )
{
    const size_type nLocalDofWithGhostOnProcStartRow = this->nLocalDofWithGhostOnProcStartRow( blockSet,g->worldComm().globalRank(), i,j );

    auto it = g->begin();
    auto const en = g->end();
    for ( ; it != en; ++it )
    {
        if ( boost::get<2>( it->second ).empty() ) continue;

        const int proc = it->second.get<0>();
        int theglobalrow = start_i + (it->first - g->firstRowEntryOnProc());
        const int thelocalrow = nLocalDofWithGhostOnProcStartRow + it->second.get<1>();

        if (it->second.get<0>()!=g->worldComm().globalRank() )
        {
            //continue;
            const size_type realrowStart = this->mapRow().firstDofGlobalCluster( proc )
                + this->nLocalDofWithoutGhostOnProcStartRow( blockSet, proc, i, j );
            theglobalrow = realrowStart+(it->first-g->mapRow().firstDofGlobalCluster( proc ));


            DCHECK( this->mapRow().searchGlobalProcessDof(theglobalrow).get<0>() )
                << " my rank " << g->worldComm().globalRank()
                << " does not contain this ghost dof " << theglobalrow
                << "in DataMapRow\n";
#if 0
            bool find=false;
            size_type gDofProcess = 0;
            boost::tie(find,gDofProcess) = this->mapRow().searchGlobalProcessDof(theglobalrow);
            if (!find) { std::cout << "STRANGE(continue) "<< std::endl; continue; }
#endif
            DCHECK(M_mapRow->mapGlobalProcessToGlobalCluster(thelocalrow) == theglobalrow)
                << " my rank " << g->worldComm().globalRank()
                << " thelocalrow " << thelocalrow
                << " M_mapRow. " << M_mapRow->mapGlobalProcessToGlobalCluster(thelocalrow)
                << " theglobalrow" << theglobalrow << "\n";
#if 0
            M_mapRow->setMapGlobalProcessToGlobalCluster(thelocalrow, theglobalrow);
#endif
            }
#if 0
        else
            {
                M_mapRow->setMapGlobalProcessToGlobalCluster(thelocalrow, theglobalrow);
                M_mapRow->setMapGlobalClusterToGlobalProcess( theglobalrow,thelocalrow );
            }
#endif

        DVLOG(2) << "rank " << this->worldComm().rank() << "update from : "
                  << " it->first " << it->first << " it->second.get<1>() " << it->second.get<1>()
                  << " into : theglobalrow " << theglobalrow << " thelocalrow " << thelocalrow << "\n";


        row_type & row = this->row( theglobalrow );
        row.get<0>() = proc;
        row.get<1>() = thelocalrow;

        std::set<size_type>& row1_entries = row.get<2>();
        std::set<size_type> const& row2_entries = boost::get<2>( it->second );
        if ( !row2_entries.empty() )
        {
            for ( auto itcol = row2_entries.begin(), encol = row2_entries.end() ; itcol!=encol; ++itcol )
            {
                if (g->mapCol().dofGlobalClusterIsOnProc(*itcol))
                {
                    const size_type dofcol = start_j + ( *itcol - g->firstColEntryOnProc() );
                    row1_entries.insert( dofcol );
                }
                else
                {
                    const int realproc = g->mapCol().procOnGlobalCluster(*itcol);
                    const size_type realcolStart = this->mapCol().firstDofGlobalCluster(realproc)
                        + this->nLocalDofWithoutGhostOnProcStartCol( blockSet, realproc, i, j );
                    const size_type dofColGC = realcolStart+ (*itcol-g->mapCol().firstDofGlobalCluster(realproc));
                    row1_entries.insert( dofColGC );
                }
            }
        }

    } // for ( ; it != en; ++it )


}

size_type
GraphCSR::nLocalDofWithoutGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int i=0; i<rowIndex; ++i )
        nDofStart += blockSet(i,colIndex)->mapRow().nLocalDofWithoutGhost( proc );

    return nDofStart;
}

size_type
GraphCSR::nLocalDofWithoutGhostOnProcStartCol( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int j=0; j<colIndex ; ++j )
        nDofStart += blockSet(rowIndex,j)->mapCol().nLocalDofWithoutGhost( proc );

    return nDofStart;
}

size_type
GraphCSR::nLocalDofWithGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int i=0; i<rowIndex ; ++i )
        nDofStart += blockSet(i,colIndex)->mapRow().nLocalDofWithGhost( proc );

    return nDofStart;
}

size_type
GraphCSR::nLocalDofWithGhostOnProcStartCol( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int j=0; j<colIndex ; ++j )
        nDofStart += blockSet(rowIndex,j)->mapCol().nLocalDofWithGhost( proc );

    return nDofStart;
}



void
GraphCSR::zero()
{
#if 0
    const size_type firstRow = this->firstRowEntryOnProc();
    const size_type rowEnd = firstRow + std::min(this->mapRow().nLocalDofWithoutGhost(),this->mapCol().nLocalDofWithoutGhost());
    for ( size_type i = firstRow ; i<rowEnd ; ++i )
    {
        if ( this->storage().find( i )!=this->end() )
        {
            row_type& row = this->row( i );
            row.get<0>() = this->worldComm().globalRank();//proc
            row.get<1>() = this->mapRow().mapGlobalClusterToGlobalProcess( i- firstRow );//local index
            row.get<2>().clear(); //all is zero
        }
    }
#else
    if ( !this->empty() )
    {
        M_storage.clear();
        M_is_closed = false;
        this->close();
    }

#endif
}

GraphCSR::self_ptrtype
GraphCSR::transpose( bool doClose )
{
    DVLOG(2) << " GraphCSR compute transpose graph\n";

    if ( M_graphT ) return M_graphT;

    //this->close();

    M_graphT = self_ptrtype( new self_type( this->mapColPtr(), this->mapRowPtr() ) );

    const int myrank = M_graphT->mapRow().worldComm().globalRank();
    const size_type myfirstGC = M_graphT->mapRow().firstDofGlobalCluster(myrank);
    for ( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
    {
        // Get the row of the sparsity pattern
        row_type const& irow = it->second;
        // num line
        size_type globalindex = it->first;

        for ( auto colit = boost::get<2>( irow ).begin(), colen=boost::get<2>( irow ).end() ; colit!=colen ; ++colit )
            {
                bool hasEntry = M_graphT->storage().find( *colit )!=M_graphT->end();

                if ( M_graphT->mapRow().dofGlobalClusterIsOnProc( *colit ) )
                {
                    self_type::row_type& row = M_graphT->row( *colit );
                    if ( !hasEntry )
                    {
                        row.get<0>()=myrank;
                        row.get<1>()= M_graphT->mapRow().mapGlobalClusterToGlobalProcess( *colit - myfirstGC );
                    }
                    row.get<2>().insert( globalindex );
                }
                else
                {
                    const int realproc = M_graphT->mapRow().procOnGlobalCluster(*colit);

                    bool find=false;
                    size_type gDofProcess = 0;
                    boost::tie(find,gDofProcess) = M_graphT->mapRow().searchGlobalProcessDof(*colit);
                    // only if find
                    if (find)
                    {
                        self_type::row_type& row = M_graphT->row( *colit );
                        if ( !hasEntry )
                        {
                            row.get<0>()=realproc;
                            row.get<1>()= gDofProcess;
                        }
                        row.get<2>().insert( globalindex );
                    }
                }
            } // for ( auto colit ... )
    }


    if ( doClose ) M_graphT->close();

    return M_graphT;
}

void
GraphCSR::addMissingZeroEntriesDiagonal()
{
    //if ( this->mapRow().worldComm().size() >1 && this->mapRow().nDof()!=this->mapCol().nDof() )
    if ( this->mapRow().nDof()!=this->mapCol().nDof() )
        return;

    DVLOG(2) << " GraphCSR addMissingZeroEntriesDiagonal in graph\n";

    const size_type firstRow = this->firstRowEntryOnProc();
    size_type firstCol = this->firstColEntryOnProc();
    const size_type rowEnd = firstRow+std::min(this->mapRow().nLocalDofWithoutGhost(),this->mapCol().nLocalDofWithoutGhost());
    for ( size_type i = firstRow ; i<rowEnd ; ++i,++firstCol )
    {
        if ( this->storage().find( firstRow )!=this->end() )
        {
            this->row( i ).get<2>().insert( i/*firstCol*/ ); // insert col on diag
        }
        else // if (this->mapCol().dofGlobalClusterIsOnProc( i/*firstCol*/ ) ) //&& this->mapCol().dofGlobalClusterIsOnProc( firstCol ))
        {
            this->row( i ).get<0>() = this->mapRow().worldComm().globalRank();//proc
            this->row( i ).get<1>() = this->mapRow().mapGlobalClusterToGlobalProcess( i-firstRow );//local index
            this->row( i ).get<2>().insert(i/*firstCol*/); // insert col on diag
        }
    }
}


void
GraphCSR::close()
{

    if ( M_is_closed )
    {
        //std::cout << "already closed graph " << this << "...\n";
        return ;
    }

    M_is_closed = true;
    //return;
    //std::cout << "closing graph " << this << "...\n";
    boost::timer ti;
    DVLOG(2) << "[close] nrows=" << this->size() << "\n";
    DVLOG(2) << "[close] firstRowEntryOnProc()=" << this->firstRowEntryOnProc() << "\n";
    DVLOG(2) << "[close] lastRowEntryOnProc()=" << this->lastRowEntryOnProc() << "\n";
    DVLOG(2) << "[close] firstColEntryOnProc()=" << this->firstColEntryOnProc() << "\n";
    DVLOG(2) << "[close] lastColEntryOnProc()=" << this->lastColEntryOnProc() << "\n";
    DVLOG(2) << "[close] M_n_total_nz=" << M_n_total_nz.size() << "\n";
    DVLOG(2) << "[close] M_storage size=" << M_storage.size() << "\n";
    DVLOG(2) << "[close] nrows=" << this->size() << "\n";

    M_n_total_nz.resize( this->mapRow().nLocalDofWithGhost(),0 );//this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1 );
    M_n_nz.resize( this->mapRow().nLocalDofWithGhost(),0 );//this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1 );
    M_n_oz.resize( this->mapRow().nLocalDofWithGhost(),0 );//this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1 );

    const int proc_id = this->worldComm().globalRank();
    const int nProc = this->worldComm().globalSize();

    std::vector<size_type> vecDofCol;//(1,0);

    std::vector< std::vector<size_type> > vecToSend( nProc );
    std::vector< std::vector<size_type> > vecToRecv( nProc );
    std::vector< std::vector<size_type> > vecToSend_nElt( nProc );
    std::vector< std::vector<size_type> > vecToRecv_nElt( nProc );

    std::vector< std::list<std::vector<size_type> > > memory_graphMPI( nProc );
    std::vector<size_type> memory_n_send(this->worldComm().globalSize() );

    for ( int proc=0 ; proc<nProc ; ++proc )
    {
        vecToSend[proc].clear();
        vecToRecv[proc].clear();
        vecToSend_nElt[proc].clear();
        vecToRecv_nElt[proc].clear();
    }

    std::fill( M_n_nz.begin(), M_n_nz.end(), 0 );
    std::fill( M_n_oz.begin(), M_n_oz.end(), 0 );

    size_type sum_nz = 0;

    M_max_nnz = 0;
    VLOG(2) << "Closing graph...";
    for ( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
    {
        // Get the row of the sparsity pattern
        row_type const& irow = it->second;

        if ( ( int )boost::get<0>( irow ) == this->worldComm().globalRank() )
        {
            size_type globalindex = it->first;
            size_type localindex = irow.get<1>();
            size_type vec_size = irow.get<2>().size();

            FEELPP_ASSERT( globalindex >= firstRowEntryOnProc() )
                ( globalindex <= lastRowEntryOnProc() )
                ( globalindex )( firstRowEntryOnProc() )
                ( lastRowEntryOnProc() ).error ( "invalid local/global index" );
            FEELPP_ASSERT( globalindex >= 0 )( globalindex < M_n_total_nz.size() )
                ( globalindex )
                ( M_n_total_nz.size() ).error ( "invalid local/global index for M_n_total_nz" );
            M_n_total_nz[localindex] = vec_size;
            sum_nz += vec_size;

            for ( auto vecit = boost::get<2>( irow ).begin(), vecen = boost::get<2>( irow ).end(); vecit != vecen; ++vecit )
            {
                if ( ( *vecit < firstColEntryOnProc() ) ||
                     ( *vecit > lastColEntryOnProc() ) )
                {
                    // entry is off block-diagonal
                    ++M_n_oz[localindex];
                }
                else
                {
                    // entry is in block-diagonal
                    ++M_n_nz[localindex];
                }
            }


#if !defined( NDEBUG )
            DVLOG(2) << "M_total_nz [  " << localindex << "]=" << M_n_total_nz[localindex] << "\n";

            DVLOG(2) << "M_nz [  " << localindex << "]=" << M_n_nz[localindex] << "\n";
            DVLOG(2) << "M_oz [  " << localindex << "]=" << M_n_oz[localindex] << "\n";
#endif // NDEBUG

            M_max_nnz = std::max( M_n_total_nz[localindex], M_max_nnz );
        }

        else
        {
            // Get the row of the sparsity pattern
            const auto dofOnGlobalCluster = it->first;
            row_type const& irow = it->second;
            const auto procOnGlobalCluster = irow.get<0>();
            const auto dofOnProc = irow.get<1>();

            vecDofCol.resize( irow.get<2>().size()+2 );
            vecDofCol[0]=dofOnGlobalCluster;
            vecDofCol[1]=dofOnProc;
            auto icol_it = irow.get<2>().begin();
            auto icol_en = irow.get<2>().end();
            for ( int i=0; icol_it!=icol_en ; ++i,++icol_it )
            {
                vecDofCol[i+2] = *icol_it;
            }

#if 0
            std::cout << "/n I am proc " << this->worldComm().globalRank()
                      << " god Rank " << this->worldComm().godRank()
                      << " global Rank " << this->worldComm().globalRank()
                      << " I send to " << procOnGlobalCluster
                      << " with dofGlobCluster " << dofOnGlobalCluster
                      << " with dofLoc " << dofOnProc
                      << " size of ivec " << ivec.size()
                      << std::endl;
#endif

            memory_graphMPI[procOnGlobalCluster].push_back(vecDofCol);//boost::make_tuple( irow.get<1>(),vecDofCol));
            memory_n_send[procOnGlobalCluster]+=vecDofCol.size();
        }

    }
    VLOG(2) << "Closing graph done.";

    if ( nProc > 1 )
    {
#if 1
        std::vector< boost::tuple< std::vector<size_type>,std::vector<size_type> > > vecToSendBase( nProc );
        std::vector< boost::tuple< std::vector<size_type>,std::vector<size_type> > > vecToRecvBase( nProc );

        const rank_type myrank = M_mapRow->worldComm().rank();

        // init data to send and counter of request
        int nbRequestToSend = 0, nbRequestToRecv = 0;
        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {
            if ( procIdNeigh < myrank )
            {
                vecToSend[procIdNeigh].resize(memory_n_send[procIdNeigh]);
                vecToSend_nElt[procIdNeigh].resize( memory_graphMPI[procIdNeigh].size() );
                auto vtsit = vecToSend[procIdNeigh].begin();
                auto it_mem = memory_graphMPI[procIdNeigh].begin();
                auto const en_mem = memory_graphMPI[procIdNeigh].end();
                for ( int cpt = 0 ; it_mem !=en_mem ; ++it_mem)
                {
                    vtsit = std::copy( it_mem->begin(), it_mem->end(), vtsit );
                    vecToSend_nElt[procIdNeigh][cpt] = it_mem->size();
                    ++cpt;
                }
                vecToSendBase[procIdNeigh] = boost::make_tuple( vecToSend[procIdNeigh],vecToSend_nElt[procIdNeigh] );
                ++nbRequestToSend;
            }
            else
            {
                ++nbRequestToRecv;
            }
        }

        // do isend/irecv
        int nbRequest = nbRequestToSend+nbRequestToRecv;
        mpi::request * reqs = new mpi::request[nbRequest];
        int cptRequest=0;
        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {
            if ( procIdNeigh < myrank )
            {
                reqs[cptRequest] = this->worldComm().globalComm().isend( procIdNeigh, 0, vecToSendBase[procIdNeigh] );
                ++cptRequest;
            }
            else
            {
                reqs[cptRequest] = this->worldComm().globalComm().irecv( procIdNeigh, 0, vecToRecvBase[procIdNeigh] );
                ++cptRequest;
            }
        }
        // wait all requests
        mpi::wait_all(reqs, reqs + nbRequest);
        // delete reqs because finish comm
        delete [] reqs;

        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {
            if ( procIdNeigh > myrank )
            {
                vecToRecv[procIdNeigh] = vecToRecvBase[procIdNeigh].get<0>();
                vecToRecv_nElt[procIdNeigh] = vecToRecvBase[procIdNeigh].get<1>();
            }
        }
        //------------------------------------------------------
        //this->worldComm().globalComm().barrier();
        //------------------------------------------------------
        for ( int proc=0; proc<nProc; ++proc )
        {
            if (vecToRecv[proc].size()>0 )
            {
                for (int cpt=0,istart=0;cpt<vecToRecv_nElt[proc].size();++cpt)
                {
                    const auto nRecvElt = vecToRecv_nElt[proc][cpt];

                    size_type globalindex = vecToRecv[proc][istart];

                    DCHECK( this->mapRow().dofGlobalClusterIsOnProc( globalindex ) ) << " GlobalCluster dofGlobalClusterIsOnProc Is not on proc "
                                                                                     << " with globalindex " << globalindex << "\n";

                    size_type localindex = this->mapRow().mapGlobalClusterToGlobalProcess( globalindex - this->mapRow().firstDofGlobalCluster()  );

                    bool hasEntry = this->storage().find( globalindex )!=this->end();
                    row_type& row = this->row( globalindex );
                    if ( !hasEntry )
                    {
                        row.get<0>()=proc_id;
                        row.get<1>()=localindex;
                    }

                    for ( int k=2;k<nRecvElt;++k)
                    {
                        bool isInserted = false;
                        if ( row.get<2>().find( vecToRecv[proc][istart+k] ) == row.get<2>().end() )
                        {
                            row.get<2>().insert( vecToRecv[proc][istart+k] );
                            isInserted = true;
                        }

                        if ( !this->mapCol().dofGlobalClusterIsOnProc(vecToRecv[proc][istart+k]) )
                        {
                            if ( isInserted )
                            {
                                ++M_n_oz[localindex];
                                ++sum_nz;
                            }
                        }

                        else
                        {
                            if ( isInserted )
                            {
                                ++M_n_nz[localindex];
                                ++sum_nz;
                            }
                        }
                    }

                    istart += nRecvElt;
                }
            }
        }
    } //if ( nProc > 1 )
#endif // MPI_MODE



    size_type nRowLoc = this->lastRowEntryOnProc()-this->firstRowEntryOnProc()+1;
    if ( nRowLoc>1 || ( sum_nz>0 )/* nRowLoc==1 && this->worldComm().globalRank()==4)*/ )
    {
        M_ia.resize( nRowLoc+1,0 );
        M_ja.resize( /*sum_n_nz*/sum_nz, 0. );
        //M_a.resize(  /*sum_n_nz*/sum_nz, 0. );
        size_type col_cursor = 0;
        auto jait = M_ja.begin();

        for ( size_type i = 0 ; i< nRowLoc; ++i )
        {
            if ( M_storage.find( this->firstRowEntryOnProc()+i )!=M_storage.end() )
            {
                row_type const& irow = this->row( this->firstRowEntryOnProc()+i );
                //size_type localindex = boost::get<1>( irow );
                M_ia[i/*localindex*/] = col_cursor;
                jait = std::copy( boost::get<2>( irow ).begin(), boost::get<2>( irow ).end(), jait );

                col_cursor+=boost::get<2>( irow ).size();
            }

            else
            {
                M_ia[i] = col_cursor;
            }
        }

        M_ia[nRowLoc] = /*sum_n_nz*/sum_nz;
    }
    else
    {
        M_ia.resize( 1,0 );
        M_ja.resize( 0 );
        M_a.resize(  0 );

    }
} // close


void
GraphCSR::showMe( std::ostream& __out ) const
{
    __out << std::endl;
    this->worldComm().globalComm().barrier();

    for ( int proc = 0; proc<this->worldComm().globalSize(); ++proc )
    {
        if ( proc==this->worldComm().globalRank() )
        {
            __out << "--------------------------------------------------------------" << std::endl;
            __out << "-------------Graph (on proc " << proc << ")------------------------------"<< std::endl;
            __out << "first_row_entry_on_proc " << this->firstRowEntryOnProc() << std::endl;
            __out << "last_row_entry_on_proc " << this->lastRowEntryOnProc() << std::endl;
            __out << "first_col_entry_on_proc " << this->firstColEntryOnProc() << std::endl;
            __out << "last_col_entry_on_proc " << this->lastColEntryOnProc() << std::endl;
            __out << "max_nnz " << M_max_nnz << std::endl;

            for ( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
            {
                // Get the row of the sparsity pattern
                row_type const& row = it->second;

                __out << " proc " << row.get<0>()
                      << " globalindex " << it->first
                      << " localindex " << row.get<1>();
#if 1
                if ( it->first>=this->firstRowEntryOnProc() && it->first<=this->lastRowEntryOnProc() )
                {
                    __out << "(nz " << M_n_nz[row.get<1>()]
                          << " oz " << M_n_oz[row.get<1>()]
                          << ") : ";

                    //size_type vec_size = boost::get<2>( irow ).size();
                    //size_type globalindex = it->first;
                    //size_type localindex = boost::get<1>( irow );
                    for ( auto it = row.get<2>().begin(), en= row.get<2>().end() ; it!=en ; ++it )
                        __out << *it << " ";
                }
#endif
                __out << std::endl;
            }

            __out << "--------------------------------------------------------------" << std::endl;

        } // if (proc==this->comm().rank())

        this->worldComm().globalComm().barrier();

    } // for (int proc = 0;proc<nProc;++proc)

}

void
GraphCSR::printPython( std::string const& nameFile ) const
{

#if 0
    std::cout << "first_row_entry_on_proc " << this->firstRowEntryOnProc() << std::endl;
    std::cout << "last_row_entry_on_proc " << this->lastRowEntryOnProc() << std::endl;
    std::cout << "first_col_entry_on_proc " << this->firstColEntryOnProc() << std::endl;
    std::cout << "last_col_entry_on_proc " << this->lastColEntryOnProc() << std::endl;
    std::cout << "max_nnz " << M_max_nnz << std::endl;
#endif

    //compute first_row_entry last_row_entry on global graph
    std::vector<size_type> first_row_entry( this->worldComm().globalSize() );
    std::vector<size_type> last_row_entry( this->worldComm().globalSize() );
    std::vector<size_type> first_col_entry( this->worldComm().globalSize() );
    std::vector<size_type> last_col_entry( this->worldComm().globalSize() );
    mpi::all_gather( this->worldComm().globalComm(),
                     (this->mapRow().nLocalDofWithoutGhost()>0)? this->firstRowEntryOnProc() : 0,
                     first_row_entry );
    mpi::all_gather( this->worldComm().globalComm(),
                     (this->mapRow().nLocalDofWithoutGhost()>0)? this->lastRowEntryOnProc() : 0,
                     last_row_entry );
    mpi::all_gather( this->worldComm().globalComm(),
                     (this->mapCol().nLocalDofWithoutGhost()>0)? this->firstColEntryOnProc() : 0,
                     first_col_entry );
    mpi::all_gather( this->worldComm().globalComm(),
                     (this->mapCol().nLocalDofWithoutGhost()>0)? this->lastColEntryOnProc() : 0,
                     last_col_entry );
    size_type thefirstRow = *std::min_element( first_row_entry.begin(),first_row_entry.end() );
    size_type thelastRow = *std::max_element( last_row_entry.begin(),last_row_entry.end() );
    size_type thefirstCol = *std::min_element( first_col_entry.begin(),first_col_entry.end() );
    size_type thelastCol = *std::max_element( last_col_entry.begin(),last_col_entry.end() );


    //std::ofstream graphFile(nameFile, std::ios::out /*| std::ios::app*/);
    std::ofstream graphFile;//(nameFile, std::ios::out

    // start file : init
    if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        graphFile.open( nameFile.c_str(), std::ios::out );

        graphFile << "import numpy" << std::endl
                  << "from scipy.sparse import * " << std::endl
                  << "from scipy import * " << std::endl
                  << "from pylab import * " << std::endl;

        graphFile << "nRow=" << thelastRow-thefirstRow+1 << std::endl
                  << "nCol=" << thelastCol-thefirstCol+1 << std::endl;

        graphFile << "mattt = array([" << std::endl;

        graphFile.close();
    }

    // synchro
    this->worldComm().barrier();

    // big part : graph
    for ( int proc = 0; proc<this->worldComm().globalSize(); ++proc )
    {
        if ( proc==this->worldComm().globalRank() )
        {
            graphFile.open( nameFile.c_str(), std::ios::out | std::ios::app );

            if (M_storage.size() > 0)
                {
                    for ( auto it = M_storage.begin(), en = --M_storage.end() ; it != en; ++it )
                        {
                            auto const& row = it->second;

                            if ( ( int )row.get<0>()==proc )
                                for ( auto it2 = row.get<2>().begin(), en2= row.get<2>().end() ; it2!=en2 ; ++it2 )
                                    graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";// << std::endl;
                        }
                    auto it = --M_storage.end();
                    auto const& row = it->second;

                    if ( ( int )row.get<0>()==proc )
                        {
                            if ( row.get<2>().size()>0 )
                                {
                                    if ( row.get<2>().size()>1 )
                                        {
                                            for ( auto it2 = row.get<2>().begin(), en2= --row.get<2>().end() ; it2!=en2 ; ++it2 )
                                                graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";
                                        }

                                    auto it2 = --row.get<2>().end();

                                    if ( proc==this->worldComm().globalSize()-1 || M_storage.size()==1 )
                                        graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ] ])" << std::endl;

                                    else
                                        graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ], " << std::endl;
                                }

                            else
                                {
                                    /*???*/
                                }
                        }
                }

            graphFile.close();
        }

        this->worldComm().barrier();
    }

    //endfile
    if ( this->worldComm().globalRank() == this->worldComm().masterRank() )
    {
        graphFile.open( nameFile.c_str(), std::ios::out | std::ios::app );

        graphFile << "row = array(mattt[:,0],dtype=int);" << std::endl
                  << "col = array(mattt[:,1],dtype=int);" << std::endl
                  << "data = array(mattt[:,2]);" << std::endl
                  << "A = csr_matrix( (data,(row,col)), shape=(nRow,nCol) );" << std::endl
                  << "fig = plt.figure();" << std::endl
                  << "matplotlib.pyplot.spy(A,precision=1e-8,aspect='equal');" << std::endl
                  << "plt.show();" << std::endl;

        graphFile.close();

    }


} // printPython

void
BlocksBaseGraphCSR::close()
{
    if ( this->isClosed() ) return;

    std::vector<boost::shared_ptr<DataMap> > dataMapRowRef(this->nRow());
    std::vector<boost::shared_ptr<DataMap> > dataMapColRef(this->nCol());

    // search a reference row datamap foreach row
    for ( index_type i=0 ; i<this->nRow() ;++i)
    {
        // search a data row avalaible
        bool findDataMapRow=false;
        for ( index_type j=0 ; j<this->nCol() && !findDataMapRow ;++j)
        {
            if ( this->operator()(i,j) )
            {
                dataMapRowRef[i] = this->operator()(i,j)->mapRowPtr();
                findDataMapRow=true;
            }
        }
        CHECK ( findDataMapRow ) << "not find at least one DataMap initialized in row "<< i << "\n";
    }

    // search reference col datamap foreach col
    for ( index_type j=0 ; j<this->nCol() ;++j)
    {
        // search a data col avalaible
        bool findDataMapCol=false;
        for ( index_type i=0 ; i<this->nRow() && !findDataMapCol ;++i)
        {
            if ( this->operator()(i,j) )
            {
                dataMapColRef[j] = this->operator()(i,j)->mapColPtr();
                findDataMapCol=true;
            }
        }
        CHECK ( findDataMapCol ) << "not find at least one DataMap initialized in col "<< j << "\n";
    }

    // if a block is missing then completed with zero graph
    for ( index_type i=0 ; i<this->nRow() ;++i)
    {
        for ( index_type j=0 ; j<this->nCol() ;++j)
        {
            if ( this->operator()(i,j) ) continue;

            DVLOG(1) << "add zero graph in block ("<<i<<","<<j<<")\n";
            typedef graph_ptrtype::element_type graph_type;
            graph_ptrtype zerograph( new graph_type( dataMapRowRef[i], dataMapColRef[j] ) );
            zerograph->zero();
            this->operator()(i,j) = zerograph;
        }
    }

    M_isClosed = true;

}


} // Feel
