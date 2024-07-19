/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

#include <feel/feelcore/application.hpp>
#include <feel/feelalg/graphcsr.hpp>

namespace Feel
{
GraphCSR::GraphCSR( size_type n,
                    size_type first_row_entry_on_proc,
                    size_type last_row_entry_on_proc,
                    size_type first_col_entry_on_proc,
                    size_type last_col_entry_on_proc,
                    worldcomm_ptr_t const& worldcomm )
    :
    super( worldcomm ),
    M_is_closed( false ),
    M_max_nnz( 0 ),
    M_n_total_nz( n, 0 ),
    M_n_nz( n, 0 ),
    M_n_oz( n, 0 ),
    M_storage(),
    M_mapRow( new datamap_type(this->worldCommPtr()) ),
    M_mapCol( new datamap_type(this->worldCommPtr()) )
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
    M_mapCol->resizeMapGlobalProcessToGlobalCluster( M_mapCol->nLocalDofWithGhost(myrank) );

    for ( size_type i=0;i<_size1;++i )
    {
        M_mapRow->setMapGlobalProcessToGlobalCluster( i,first_row_entry_on_proc+i );
    }

    for( size_type j=0;j<_size2;++j )
    {
        M_mapCol->setMapGlobalProcessToGlobalCluster( j,first_col_entry_on_proc+j );
    }

}


GraphCSR::GraphCSR( datamap_ptrtype const& mapRow,
                    datamap_ptrtype const& mapCol )
    :
    super( mapRow->worldCommPtr() ),
    M_is_closed( false ),
    M_max_nnz( 0 ),
    M_n_total_nz( 0 ),
    M_n_nz( 0 ),
    M_n_oz( 0 ),
    M_storage(),
    M_mapRow(mapRow),
    M_mapCol(mapCol)
{}

GraphCSR::GraphCSR( datamap_type const& mapRow,
                    datamap_type const& mapCol )
    :
    super( mapRow.worldCommPtr() ),
    M_is_closed( false ),
    M_max_nnz( 0 ),
    M_n_total_nz( 0 ),
    M_n_nz( 0 ),
    M_n_oz( 0 ),
    M_storage(),
    M_mapRow( new datamap_type(mapRow) ),
    M_mapCol( new datamap_type(mapCol) )
{}




GraphCSR::GraphCSR( vf::BlocksBase<self_ptrtype> const & blockSet,
                    bool diagIsNonZero, bool close )
    :
    super( blockSet(0,0)->worldCommPtr() ), 
    M_is_closed( false ),
    M_max_nnz( 0 ),
    M_n_total_nz( /*n*/0, 0 ),
    M_n_nz( /*n*/0, 0 ),
    M_n_oz( /*n*/0, 0 ),
    M_storage(),
    M_mapRow(),
    M_mapCol()
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
    //M_mapCol.setMapGlobalProcessToGlobalCluster(M_mapRow.mapGlobalProcessToGlobalCluster());

    if ( diagIsNonZero ) this->addMissingZeroEntriesDiagonal();

    if ( close ) this->close();
}



GraphCSR::GraphCSR( GraphCSR const & g )
    :
    super( g ),
    M_is_closed( g.M_is_closed ),
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
        super::operator=( g );
        M_max_nnz = g.M_max_nnz;
        M_n_total_nz= g.M_n_total_nz;
        M_n_nz = g.M_n_nz;
        M_n_oz = g.M_n_oz;
        M_storage = g.M_storage;
        M_graphT = g.M_graphT;
        M_is_closed = g.M_is_closed;
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
    if ( nRow==1 && nCol==1)
    {
        M_mapRow = blockSet(0,0)->mapRowPtr();
        M_mapCol = blockSet(0,0)->mapColPtr();
    }
    else
    {
        std::vector<datamap_ptrtype > listofdmRow, listofdmCol;
        for ( uint i=0; i<nRow; ++i )
            listofdmRow.push_back( blockSet(i,0)->mapRowPtr() );
        for ( uint i=0; i<nCol; ++i )
            listofdmCol.push_back( blockSet(0,i)->mapColPtr() );
        M_mapRow.reset( new datamap_type( listofdmRow, this->worldCommPtr() ) );
        M_mapCol.reset( new datamap_type( listofdmCol, this->worldCommPtr() ) );
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
    const size_type nLocalDofWithoutGhostOnProcStartRow = this->nLocalDofWithoutGhostOnProcStartRow( blockSet,g->worldComm().globalRank(), i,j );
    rank_type myrank=g->worldComm().globalRank();
    auto it = g->begin();
    auto const en = g->end();
    for ( ; it != en; ++it )
    {
        if ( boost::get<2>( it->second ).empty() ) continue;

        const int proc = it->second.get<0>();
        int theglobalrow = start_i + (it->first - g->firstRowEntryOnProc());
        int thelocalrow = nLocalDofWithoutGhostOnProcStartRow + it->second.get<1>();

        if (it->second.get<0>()!=g->worldComm().globalRank() )
        {
            const size_type realrowStart = this->mapRow().firstDofGlobalCluster( proc )
                + this->nLocalDofWithoutGhostOnProcStartRow( blockSet, proc, i, j );
            theglobalrow = realrowStart+(it->first-g->mapRow().firstDofGlobalCluster( proc ));

            thelocalrow = this->mapRow().nLocalDofWithoutGhost()
                + ( nLocalDofWithGhostOnProcStartRow - nLocalDofWithoutGhostOnProcStartRow)
                + (it->second.get<1>() - g->mapRow().nLocalDofWithoutGhost());

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
            }

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

typename GraphCSR::size_type
GraphCSR::nLocalDofWithoutGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int i=0; i<rowIndex; ++i )
        nDofStart += blockSet(i,colIndex)->mapRow().nLocalDofWithoutGhost( proc );

    return nDofStart;
}

typename GraphCSR::size_type
GraphCSR::nLocalDofWithoutGhostOnProcStartCol( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int j=0; j<colIndex ; ++j )
        nDofStart += blockSet(rowIndex,j)->mapCol().nLocalDofWithoutGhost( proc );

    return nDofStart;
}

typename GraphCSR::size_type
GraphCSR::nLocalDofWithGhostOnProcStartRow( vf::BlocksBase<self_ptrtype> const & blockSet, int proc, int rowIndex, int colIndex ) const
{
    size_type nDofStart=0;
    for ( int i=0; i<rowIndex ; ++i )
        nDofStart += blockSet(i,colIndex)->mapRow().nLocalDofWithGhost( proc );

    return nDofStart;
}

typename GraphCSR::size_type
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
    if ( !this->empty() )
    {
        M_storage.clear();
        M_is_closed = false;
        this->close();
    }
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
                        row.get<1>()= *colit - myfirstGC;
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


GraphCSR::self_ptrtype
GraphCSR::createSubGraph( std::vector<size_type> const& _rows, std::vector<size_type> const& _cols,
                          datamap_ptrtype const& _subMapRow, datamap_ptrtype const& _subMapCol,
                          bool useSameDataMap, bool checkAndFixRange ) const
{
    auto const& theMapRow = this->mapRow();
    auto const& theMapCol = this->mapCol();

    // update maybe input index set
    std::vector<size_type> rows = ( checkAndFixRange )?
        this->mapRowPtr()->buildIndexSetWithParallelMissingDof( _rows ) : _rows;
    std::vector<size_type> cols = ( checkAndFixRange && !useSameDataMap )?
        this->mapColPtr()->buildIndexSetWithParallelMissingDof( _cols ) : ( useSameDataMap )? rows : _cols;

    // build subdatamp if not given
    datamap_ptrtype subMapRow = (!_subMapRow)? this->mapRowPtr()->createSubDataMap( rows, false ) : _subMapRow;
    datamap_ptrtype subMapCol;
    if ( !_subMapCol )
        subMapCol = ( useSameDataMap )? subMapRow : this->mapColPtr()->createSubDataMap( cols, false );
    else
        subMapCol = _subMapCol;


    const size_type firstDofGCrow = theMapRow.firstDofGlobalCluster();
    const size_type firstDofGCcol = theMapCol.firstDofGlobalCluster();

    // convert input idExtract into std::set ( preserve ordering and remove duplicated dofs )
    std::set<size_type> idExtractRow;
    idExtractRow.insert( rows.begin(), rows.end() );
    CHECK( idExtractRow.size() == subMapRow->nLocalDofWithGhost() ) << "invalid size " << idExtractRow.size() << " vs "<< subMapRow->nLocalDofWithGhost();
    std::set<size_type> idExtractCol;
    idExtractCol.insert( cols.begin(), cols.end() );
    CHECK( idExtractCol.size() == subMapCol->nLocalDofWithGhost() ) << "invalid size " << idExtractCol.size() << " vs "<< subMapCol->nLocalDofWithGhost();

    // container ( with optimized access ) for global dof active on proc ( pair : globalProcessDof, globalClusterDof )
    std::vector<std::pair<size_type,size_type> > relationDofActiveRow( theMapRow.nLocalDofWithoutGhost(),std::make_pair(invalid_v<size_type>,invalid_v<size_type>) );
    // container ( map ) for non active global dof ( pair : globalProcessDof, globalClusterDof )
    std::map<size_type,std::pair<size_type,size_type> > relationDofNonActiveRow;
    size_type subDofGPCounter = 0;
    for ( size_type id : idExtractRow )
    {
        const size_type gdof = theMapRow.mapGlobalProcessToGlobalCluster( id );
        const size_type subDofGC = subMapRow->mapGlobalProcessToGlobalCluster( subDofGPCounter );
        if ( theMapRow.dofGlobalClusterIsOnProc( gdof ) )
            relationDofActiveRow[ gdof-firstDofGCrow] = std::make_pair(subDofGPCounter,subDofGC);
        else
            relationDofNonActiveRow[ gdof ] = std::make_pair(subDofGPCounter,subDofGC);
        ++subDofGPCounter;
    }

    // container ( with optimized access ) for global dof active on proc
    std::vector<size_type> relationDofActiveCol( theMapCol.nLocalDofWithoutGhost(),invalid_v<size_type> );
    // container ( map ) for non active global dof
    std::map<size_type,size_type> relationDofNonActiveCol;
    subDofGPCounter = 0;
    for ( size_type id : idExtractCol )
    {
        const size_type gdof = theMapCol.mapGlobalProcessToGlobalCluster( id );
        const size_type subDofGC = subMapCol->mapGlobalProcessToGlobalCluster( subDofGPCounter );
        if ( theMapCol.dofGlobalClusterIsOnProc( gdof ) )
            relationDofActiveCol[ gdof-firstDofGCcol] = subDofGC;
        else
            relationDofNonActiveCol[ gdof ] = subDofGC;
        ++subDofGPCounter;
    }

    // build sub graph
    self_ptrtype subgraph( new self_type(subMapRow, subMapCol) );

    for ( auto it = M_storage.begin(), en = M_storage.end() ; it != en; ++it )
    {
        size_type gdofRow = it->first;

        size_type subDofRowGP = invalid_v<size_type>, subDofRowGC = invalid_v<size_type>;
        if ( theMapRow.dofGlobalClusterIsOnProc( gdofRow ) )
        {
            subDofRowGP = relationDofActiveRow[gdofRow-firstDofGCrow].first;
            subDofRowGC = relationDofActiveRow[gdofRow-firstDofGCrow].second;
        }
        else
        {
            auto itFindDof = relationDofNonActiveRow.find( gdofRow );
            if ( itFindDof != relationDofNonActiveRow.end() )
            {
                subDofRowGP = itFindDof->second.first;
                subDofRowGC = itFindDof->second.second;
            }
        }
        // ignore this row if dof not present in extraction
        if ( subDofRowGC == invalid_v<size_type> ) continue;

        // Get the row of the sparsity pattern
        row_type const& datarow = it->second;

        row_type& subdatarow = subgraph->row(subDofRowGC);
        subdatarow.get<0>() = datarow.get<0>(); // same proc
        subdatarow.get<1>() = subDofRowGP; // global process dof

        // iterate on each column dof
        for ( size_type dofColGC : datarow.get<2>() )
        {
            size_type subDofColGC = invalid_v<size_type>;
            if ( theMapCol.dofGlobalClusterIsOnProc( dofColGC ) )
            {
                subDofColGC = relationDofActiveCol[dofColGC-firstDofGCcol];
            }
            else
            {
                auto itFindDof = relationDofNonActiveCol.find( dofColGC );
                if ( itFindDof != relationDofNonActiveCol.end() )
                {
                    subDofColGC = itFindDof->second;
                }
            }
            // add column indice if dof present in extraction
            if ( subDofColGC != invalid_v<size_type> )
                subdatarow.get<2>().insert( subDofColGC );
        }
    }

    subgraph->close();

    return subgraph;
}


void
GraphCSR::addMissingZeroEntriesDiagonal()
{
    //if ( this->mapRow().worldComm().size() >1 && this->mapRow().nDof()!=this->mapCol().nDof() )
    if ( this->mapRow().nDof()!=this->mapCol().nDof() ) // only square matrix
        return;

    DVLOG(2) << " GraphCSR addMissingZeroEntriesDiagonal in graph\n";

    const size_type firstRow = this->firstRowEntryOnProc();
    //size_type firstCol = this->firstColEntryOnProc();
    //const size_type rowEnd = firstRow+std::min(this->mapRow().nLocalDofWithoutGhost(),this->mapCol().nLocalDofWithoutGhost());
    const size_type rowEnd = firstRow + this->mapRow().nLocalDofWithoutGhost();
    rank_type procId = this->mapRow().worldComm().globalRank();
    for ( size_type i = firstRow ; i<rowEnd ; ++i/*,++firstCol*/ )
    {
        auto & row = this->row( i );
        //! init if empty
        if ( row.get<2>().empty() )
        {
            row.get<0>() = procId;
            row.get<1>() = i-firstRow; //local index (just a shift because active dof)
        }
        //! insert diagonal entry
        row.get<2>().insert( i ); // insert col on diag
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

    std::fill( M_n_total_nz.begin(), M_n_total_nz.end(), 0 );
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
                if ( ( M_mapCol->nLocalDofWithoutGhost() == 0 ) ||
                     ( *vecit < firstColEntryOnProc() ) || ( *vecit > lastColEntryOnProc() ) )
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
        VLOG(2) << "Closing graph parallel work start.";
#if 1
        std::vector< boost::tuple< std::vector<size_type>,std::vector<size_type> > > vecToSendBase( nProc );
        std::vector< boost::tuple< std::vector<size_type>,std::vector<size_type> > > vecToRecvBase( nProc );

        const rank_type myrank = M_mapRow->worldComm().rank();

        // init data to send and counter of request
        int nbRequestToSend = 0, nbRequestToRecv = 0;
        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {
#if 0
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
#else
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
                ++nbRequestToRecv;
#endif
        }

        // do isend/irecv
        int nbRequest = nbRequestToSend+nbRequestToRecv;
        mpi::request * reqs = new mpi::request[nbRequest];
        int cptRequest=0;
        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {
#if 0
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
#else
            reqs[cptRequest++] = this->worldComm().globalComm().isend( procIdNeigh, 0, vecToSendBase[procIdNeigh] );
            reqs[cptRequest++] = this->worldComm().globalComm().irecv( procIdNeigh, 0, vecToRecvBase[procIdNeigh] );
#endif
        }
        VLOG(2) << "Closing graph parallel wait all process " << nbRequest;
        // wait all requests
        mpi::wait_all(reqs, reqs + nbRequest);
        // delete reqs because finish comm
        if ( nbRequest )
            delete [] reqs;

        for ( rank_type procIdNeigh : M_mapRow->neighborSubdomains() )
        {

#if 0
            if ( procIdNeigh > myrank )
#endif
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

                    size_type localindex = globalindex - this->mapRow().firstDofGlobalCluster();

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
        VLOG(2) << "Closing graph parallel work done ";
    } //if ( nProc > 1 )
#endif // MPI_MODE


#if 0
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
#else
    VLOG(2) << "Closing graph almost done";
    size_type nRowLoc = this->mapRow().nLocalDofWithoutGhost();
    VLOG(2) << "Closing graph almost done " << nRowLoc << "," << sum_nz;

    M_ia.resize( nRowLoc+1,0 );
    M_ja.resize( sum_nz, 0. );

    //M_a.resize(  /*sum_n_nz*/sum_nz, 0. );
    size_type col_cursor = 0;
    auto jait = M_ja.begin();

    for ( size_type i = 0 ; i< nRowLoc; ++i )
    {
        if ( M_storage.find( this->firstRowEntryOnProc()+i ) != M_storage.end() )
        {
            row_type const& irow = this->row( this->firstRowEntryOnProc()+i );
            M_ia[i] = col_cursor;
            jait = std::copy( boost::get<2>( irow ).begin(), boost::get<2>( irow ).end(), jait );
            col_cursor+=boost::get<2>( irow ).size();
        }
        else
        {
            M_ia[i] = col_cursor;
        }
    }

    M_ia[nRowLoc] = sum_nz;
    VLOG(2) << "Closing graph done";
#endif
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
                    for ( auto it = M_storage.begin(), en = std::prev(M_storage.end()) ; it != en; ++it )
                        {
                            auto const& row = it->second;

                            if ( ( int )row.get<0>()==proc )
                                for ( auto it2 = row.get<2>().begin(), en2= row.get<2>().end() ; it2!=en2 ; ++it2 )
                                    graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";// << std::endl;
                        }
                    auto it = std::prev(M_storage.end());
                    auto const& row = it->second;

                    if ( ( int )row.get<0>()==proc )
                        {
                            if ( row.get<2>().size()>0 )
                                {
                                    if ( row.get<2>().size()>1 )
                                        {
                                            for ( auto it2 = row.get<2>().begin(), en2= std::prev(row.get<2>().end()) ; it2!=en2 ; ++it2 )
                                                graphFile << "[" << it->first << " , " << *it2 << " , 1.0 ],";
                                        }

                                    auto it2 = std::prev(row.get<2>().end());

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

#endif
} // printPython

void
BlocksBaseGraphCSR::close()
{
    if ( this->isClosed() ) return;

    std::vector<datamap_ptrtype<>> dataMapRowRef(this->nRow());
    std::vector<datamap_ptrtype<>> dataMapColRef(this->nCol());

    // search a reference row datamap foreach row
    for ( index_type i=0 ; i<this->nRow() ;++i)
    {
        // search a data row available
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
        // search a data col available
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
