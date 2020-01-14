/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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


template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::localize( vector_ptrtype const& vb, size_type _start_i )
{
    if ( !vb->closed() )
        vb->close();

    auto const& dm = vb->map();
    int currentDataBaseId = _start_i;
    for ( uint16_type i=0; i<this->nRow(); ++i )
    {
        //size_type nBlockRow = this->operator()( i,0 )->localSize();
        auto const& dmb = this->operator()( i,0 )->map();
        int nDataBase = dmb.numberOfDofIdToContainerId();
        for ( int tag=0;tag<nDataBase;++tag,++currentDataBaseId )
        {
            auto const& dofIdToContainerIdBlock = dmb.dofIdToContainerId(tag);
            auto const& dofIdToContainerIdVec = dm.dofIdToContainerId(currentDataBaseId);
            for (int k=0;k<dofIdToContainerIdBlock.size();++k)
            {
                this->operator()( i,0 )->set( dofIdToContainerIdBlock[k], vb->operator()( dofIdToContainerIdVec[k] ) );
            }
        }
        this->operator()( i,0 )->close();
    }
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::localize()
{
    if ( !this->vector() ) return;

    this->localize( this->vector() );
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::buildVector( backend_ptrtype backend )
{
    M_vector = backend->newBlockVector( _block=*this );
}

template <typename T, typename SizeT>
typename BlocksBaseVector<T, SizeT>::vector_ptrtype&
BlocksBaseVector<T, SizeT>::vectorMonolithic()
{
    CHECK( M_vector ) << "vector not initialized";
    std::shared_ptr< VectorBlockBase<T> > vcast = std::dynamic_pointer_cast< VectorBlockBase<T> >( M_vector );
    if ( vcast )
        return vcast->getVector();
    else
        return M_vector;
}
template <typename T, typename SizeT>
typename BlocksBaseVector<T, SizeT>::vector_ptrtype const&
BlocksBaseVector<T, SizeT>::vectorMonolithic() const
{
    CHECK( M_vector ) << "vector not initialized";
    std::shared_ptr< VectorBlockBase<T> const > vcast = std::dynamic_pointer_cast< VectorBlockBase<T> const >( M_vector );
    if ( vcast )
        return vcast->getVector();
    else
        return M_vector;
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::setVector( vector_type & vec, vector_type const& subvec, int dtId, bool closeVector ) const
{
    auto const& dmVec = vec.map();
    auto const& dmSubVec = subvec.map();
    for ( int tag=0 ; tag<dmSubVec.numberOfDofIdToContainerId() ; ++tag )
    {
        auto const& basisGpToContainerGpSubVec = dmSubVec.dofIdToContainerId( tag );
        CHECK( dtId+tag < dmVec.numberOfDofIdToContainerId() ) << "error "<<dtId+tag << " vs " << dmVec.numberOfDofIdToContainerId();
        auto const& basisGpToContainerGpVec = dmVec.dofIdToContainerId( dtId+tag );
        CHECK( basisGpToContainerGpSubVec.size() == basisGpToContainerGpVec.size() ) << " aii " << basisGpToContainerGpSubVec.size() << " vs " << basisGpToContainerGpVec.size();
        for ( int k=0;k<basisGpToContainerGpSubVec.size();++k )
            vec.set( basisGpToContainerGpVec[k], subvec( basisGpToContainerGpSubVec[k] ) );
            //vec( basisGpToContainerGpVec[k] ) = subvec( basisGpToContainerGpSubVec[k] );
    }
    if ( closeVector )
        vec.close();
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::setSubVector( vector_type & subvec, vector_type const& vec , int idStart ) const
{
    auto const& dmVec = vec.map();
    auto const& dmSubVec = subvec.map();
    int basisIndexSubVec = idStart;
    for ( int tag=0 ; tag<dmSubVec.numberOfDofIdToContainerId() ; ++tag )
    {
        auto const& basisGpToContainerGpSubVec = dmSubVec.dofIdToContainerId( tag );
        CHECK( basisIndexSubVec+tag < dmVec.numberOfDofIdToContainerId() ) << "error "<<basisIndexSubVec+tag<< " vs " << dmVec.numberOfDofIdToContainerId();
        auto const& basisGpToContainerGpVec = dmVec.dofIdToContainerId( basisIndexSubVec+tag );
        CHECK( basisGpToContainerGpSubVec.size() == basisGpToContainerGpVec.size() ) << " error " << basisGpToContainerGpSubVec.size() << " vs " << basisGpToContainerGpVec.size();
        for ( int k=0;k<basisGpToContainerGpSubVec.size();++k )
            subvec.set( basisGpToContainerGpSubVec[k], vec( basisGpToContainerGpVec[k] ) );
    }
    subvec.close();
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::updateVectorFromSubVectors( vector_type & vec ) const
{
    int nBlock = this->nRow();
    for ( int k = 0, dtId = 0 ; k<nBlock ;++k )
    {
        vector_ptrtype subvec = this->operator()(k);
        this->setVector( vec, *subvec, dtId, false );
        dtId += subvec->map().numberOfDofIdToContainerId();
    }
    vec.close();
}

template <typename T, typename SizeT>
void
BlocksBaseVector<T, SizeT>::updateVectorFromSubVectors()
{
    if ( !M_vector )
        return;
    this->updateVectorFromSubVectors( *M_vector );
}

template class BlocksBaseVector<double,uint32_type>;


template <typename T, typename SizeT>
VectorBlockBase<T,SizeT>::VectorBlockBase( BlocksBaseVector<T> const & blockVec,
                                     backend_type &backend,
                                     bool copy_values )
    :
    M_backend( backend.shared_from_this() ),
    M_vec()
    
{
    auto nRow = blockVec.nRow();

    std::shared_ptr<DataMap<>> dm;
    if ( nRow == 1 )
    {
        dm = blockVec(0,0)->mapPtr();
    }
    else
    {
        std::vector<std::shared_ptr<DataMap<>> > listofdm;
        for ( uint16_type i=0 ; i<nRow; ++i )
            listofdm.push_back( blockVec(i,0)->mapPtr() );
        dm.reset( new DataMap<>( listofdm, blockVec(0,0)->map().worldCommPtr() ) );
    }
    M_vec = backend.newVector( dm );
    this->setMap( M_vec->mapPtr() );

    if ( copy_values )
    {
        blockVec.updateVectorFromSubVectors( *M_vec );
    }

}

template <typename T, typename SizeT>
void
VectorBlockBase<T,SizeT>::updateBlockVec( vector_ptrtype const& m, size_type start_i )
{
    auto const& dmb = m->map();
    auto const& dm = M_vec->map();

    // int startTagIdInVec=0;
    // for ( int i=0; i<start_i; ++i )
    //     startTagIdInVec += this->operator[]( i,0 )->map().indexSplit()->size();

    int nTag = dmb.numberOfDofIdToContainerId();
    for ( int tag=0;tag<nTag;++tag )
    {
        auto const& dofIdToContainerIdBlock = dmb.dofIdToContainerId(tag);
        auto const& dofIdToContainerIdVec = dm.dofIdToContainerId(start_i+tag);
        CHECK( dofIdToContainerIdBlock.size() == dofIdToContainerIdVec.size() ) << "incompatibility with size";
        for (int k=0;k<dofIdToContainerIdBlock.size();++k)
            M_vec->set( dofIdToContainerIdVec[k],m->operator()( dofIdToContainerIdBlock[k] ) );
    }
    //this->setMap( M_vec->mapPtr() );
}

template <typename T, typename SizeT>
void
VectorBlockBase<T,SizeT>::addVector ( int* rows, int nrows, value_type* data, size_type K, size_type K2 )
{
    M_vec->addVector( rows, nrows, data, K, K2 );
}
template class VectorBlockBase<double,uint32_type>;

} // Feel


