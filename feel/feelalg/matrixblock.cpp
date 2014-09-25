/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2008-01-03

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
   \file matrixblock.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-06-10
 */

#include <feel/feelalg/matrixblock.hpp>



namespace Feel
{
namespace detail {
    
    template<typename T>
    typename Backend<T>::ptrtype backend( T t  ) {}

    template<>
    Backend<double>::ptrtype backend( double t  ) { return Feel::backend(); }

    template<>
    Backend<std::complex<double>>::ptrtype backend( std::complex<double> t  ) { return Feel::cbackend(); }
    
        

}

template <typename T>
void
BlocksBaseSparseMatrix<T>::close()
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

    // if a block is missing then completed with zero matrix
    for ( index_type i=0 ; i<this->nRow() ;++i)
    {
        for ( index_type j=0 ; j<this->nCol() ;++j)
        {
            if ( this->operator()(i,j) ) continue;

            DVLOG(1) << "add zero matrix in block ("<<i<<","<<j<<")\n";
            
            this->operator()(i,j) = Feel::detail::backend(T(0))->newZeroMatrix(dataMapRowRef[i], dataMapColRef[j]);
        }
    }

    M_isClosed = true;

}


template <typename T>
MatrixBlockBase<T>::MatrixBlockBase( vf::BlocksBase<matrix_ptrtype> const & blockSet,
                                     backend_ptrtype backend,
                                     bool copy_values,
                                     bool diag_is_nonzero )
    :
    super(),
    M_backend(backend),
    M_mat()
{
    const uint16_type nRow = blockSet.nRow();
    const uint16_type nCol = blockSet.nCol();

    BlocksBaseGraphCSR blockGraph(nRow,nCol);
    for ( uint16_type i=0; i<nRow; ++i )
        for ( uint16_type j=0; j<nCol; ++j )
            blockGraph(i,j) = blockSet(i,j)->graph();

    M_mat = M_backend->newBlockMatrix(_block=blockGraph);

    if ( copy_values )
    {
        auto const& mapRow = M_mat->mapRow();
        auto const& mapCol = M_mat->mapCol();
        const int worldsize = mapRow.worldComm().globalSize();
        const int myrank = mapRow.worldComm().globalRank();

        std::vector<size_type> start_i = M_mat->mapRow().firstDofGlobalClusterWorld();
        for ( uint16_type i=0; i<nRow; ++i )
        {
            std::vector<size_type> start_j = M_mat->mapCol().firstDofGlobalClusterWorld();

            for ( uint16_type j=0; j<nCol; ++j )
            {
                blockSet(i,j)->close();
                this->updateBlockMat( blockSet(i,j),start_i,start_j );

                for ( int proc=0 ;proc < worldsize;++proc )
                    start_j[proc] += blockSet(i,j)->mapCol().nLocalDofWithoutGhost(proc);
            }

            for ( int proc=0 ;proc < worldsize;++proc )
                start_i[proc] += blockSet(i,0)->mapRow().nLocalDofWithoutGhost(proc);
        }
    }

}


template <typename T>
MatrixBlockBase<T>::MatrixBlockBase( vf::BlocksBase<graph_ptrtype> const & blockgraph, 
                                     backend_ptrtype backend,
                                     bool diag_is_nonzero )
    :
    super(),
    M_backend(backend),
    M_mat()
{
    graph_ptrtype graph( new graph_type( blockgraph,diag_is_nonzero,true) );
    //graph->showMe();
    //graph->mapRow().showMeMapGlobalProcessToGlobalCluster();
    //graph->mapCol().showMeMapGlobalProcessToGlobalCluster();
    //graph->printPython("GraphG.py");

    size_type properties = NON_HERMITIAN;
    M_mat = M_backend->newMatrix( graph->mapColPtr(),  graph->mapRowPtr(), properties, false );
    M_mat->init( graph->mapRow().nDof(), graph->mapCol().nDof(),
                 graph->mapRow().nLocalDofWithoutGhost(), graph->mapCol().nLocalDofWithoutGhost(),
                 graph );

    M_mat->zero();

    M_mat->setIndexSplit( graph->mapRow().indexSplit() );

#if 0
    bool computeIndexSplit = true;
    if ( computeIndexSplit )
    {
#if 0
        const uint16_type nRow = blockgraph.nRow();
        const uint16_type nCol = blockgraph.nCol();
        // index container for field split preconditioner
        std::vector< std::vector<size_type> > indexSplit( nRow );
        size_type startIS = M_mat->mapRow().firstDofGlobalCluster();
        for ( uint16_type i=0; i<nRow; ++i )
        {
            auto const& mapRowBlock = blockgraph(i,0)->mapRow();
            const size_type nLocDofWithoutGhostBlock = mapRowBlock.nLocalDofWithoutGhost();
            const size_type nLocDofWithGhostBlock = mapRowBlock.nLocalDofWithGhost();
            const size_type firstDofGCBlock = mapRowBlock.firstDofGlobalCluster();

            indexSplit[i].resize( nLocDofWithoutGhostBlock );

            for ( size_type l = 0; l< nLocDofWithGhostBlock ; ++l )
            {
                if ( mapRowBlock.dofGlobalProcessIsGhost(l) ) continue;

                const size_type globalDof = mapRowBlock.mapGlobalProcessToGlobalCluster(l);
                indexSplit[i][globalDof - firstDofGCBlock ] = startIS + (globalDof - firstDofGCBlock);
            }

            startIS += nLocDofWithoutGhostBlock;
        }
#else
        const uint16_type nRow = blockgraph.nRow();
        indexsplit_ptrtype indexSplit( new IndexSplit() );
        const size_type firstDofGC = graph->mapRow().firstDofGlobalCluster();
        for ( uint16_type i=0; i<nRow; ++i )
        {
            indexSplit->addSplit( firstDofGC, blockgraph(i,0)->mapRow().indexSplit() );
        }
        //indexSplit.showMe();
        graph->mapRowPtr()->setIndexSplit( indexSplit );
#endif

        // update
        M_mat->setIndexSplit( indexSplit );
    }
#endif

}



template <typename T>
void
MatrixBlockBase<T>::init ( const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           const size_type nnz,
                           const size_type /*noz*/ )
{
    if ( ( m==0 ) || ( n==0 ) )
        return;

    {
        // Clear initialized matrices
        if ( this->isInitialized() )
            this->clear();

        this->setInitialized( true );
    }

    M_mat->init( m,n,m_l,n_l,nnz );

}

template <typename T>
void
MatrixBlockBase<T>::init ( const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           graph_ptrtype const& graph )
{
    this->setGraph( graph );

    {
        // Clear initialized matrices
        if ( this->isInitialized() )
            this->clear();

        this->setInitialized(  true );
    }

    M_mat->init( m,n,m_l,n_l,graph );
}


template <typename T>
void
MatrixBlockBase<T>::zero ()
{
    M_mat->zero();
}

template <typename T>
void
MatrixBlockBase<T>::zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{
    M_mat->zero();
}

template <typename T>
void
MatrixBlockBase<T>::clear ()
{
    M_mat->clear();
}

template <typename T>
inline
void
MatrixBlockBase<T>::close () const
{
    M_mat->close();
}

template <typename T>
inline
size_type
MatrixBlockBase<T>::size1 () const
{
    return M_mat->size1();
}

template <typename T>
inline
size_type
MatrixBlockBase<T>::size2 () const
{
    return M_mat->size2();
}

template <typename T>
inline
size_type
MatrixBlockBase<T>::rowStart () const
{
    return M_mat->rowStart();
}

template <typename T>
inline
size_type
MatrixBlockBase<T>::rowStop () const
{
    return M_mat->rowStop();
}

template <typename T>
inline
void
MatrixBlockBase<T>::set ( const size_type i,
                          const size_type j,
                          const value_type& value )
{
    M_mat->set( i,j,value );
}

template <typename T>
inline
void
MatrixBlockBase<T>::add ( const size_type i,
                          const size_type j,
                          const value_type& value )
{
    M_mat->add( i,j,value );
}

template <typename T>
inline
bool
MatrixBlockBase<T>::closed() const
{
    return M_mat->closed();
}

template <typename T>
void
MatrixBlockBase<T> ::addMatrix( const ublas::matrix<value_type>& dm,
                                const std::vector<size_type>& rows,
                                const std::vector<size_type>& cols )
{
    M_mat->addMatrix( dm,rows,cols );
}

template <typename T>
void
MatrixBlockBase<T>::addMatrix ( int* rows, int nrows,
                                int* cols, int ncols,
                                value_type* data )
{
    M_mat->addMatrix( rows,nrows,cols,ncols,data );
}

template <typename T>
void
MatrixBlockBase<T>::printMatlab ( const std::string name ) const
{
    M_mat->printMatlab( name );
}

template <typename T>
inline
void
MatrixBlockBase<T>::addMatrix ( const value_type a_in, MatrixSparse<value_type> &X_in )
{
    M_mat->addMatrix( a_in,X_in );
}

template <typename T>
void
MatrixBlockBase<T>::scale( value_type const a )
{
    M_mat->scale( a );
}

template <typename T>
typename MatrixBlockBase<T>::real_type
MatrixBlockBase<T>::energy( Vector<value_type> const& __v,
                            Vector<value_type> const& __u,
                            bool _transpose ) const
{
    return M_mat->energy( __v,__u,_transpose );
}

template <typename T>
inline
typename MatrixBlockBase<T>::real_type
MatrixBlockBase<T>::l1Norm() const
{
    return M_mat->l1Norm();
}

template <typename T>
inline
typename MatrixBlockBase<T>::real_type
MatrixBlockBase<T>::linftyNorm() const
{
    return M_mat->linftyNorm();
}

template <typename T>
inline
typename MatrixBlockBase<T>::value_type
MatrixBlockBase<T>::operator () ( const size_type i,
                                  const size_type j ) const
{
    return ( *M_mat )( i,j );
}


template <typename T>
MatrixBlockBase<T> &
MatrixBlockBase<T>::operator = ( MatrixSparse<value_type> const& M )
{
    *M_mat = M;
    return *this;
}

template <typename T>
void
MatrixBlockBase<T>::zeroRows( std::vector<int> const& rows, Vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context )
{
    M_mat->zeroRows( rows,values,rhs,on_context );
}

template <typename T>
void
MatrixBlockBase<T>::diagonal( Vector<value_type>& out ) const
{
    M_mat->diagonal( out );
}

template <typename T>
void
MatrixBlockBase<T>::transpose( MatrixSparse<value_type>& Mt, size_type options ) const
{
    M_mat->transpose( Mt, options );
}

template <typename T>
void
MatrixBlockBase<T>::updateBlockMat( boost::shared_ptr<MatrixSparse<value_type> > m, std::vector<size_type> start_i, std::vector<size_type> start_j )
{
    M_mat->updateBlockMat( m,start_i,start_j );
}

template class BlocksBaseSparseMatrix<double>;
template class MatrixBlockBase<double>;
template class BlocksBaseSparseMatrix<std::complex<double>>;
template class MatrixBlockBase<std::complex<double>>;

} // Feel
