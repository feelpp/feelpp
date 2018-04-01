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


template <typename T>
void
BlocksBaseVector<T>::localize( vector_ptrtype const& vb, size_type _start_i )
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

template <typename T>
void
BlocksBaseVector<T>::localize()
{
    if ( !this->vector() ) return;

    this->localize( this->vector() );
}

template <typename T>
void
BlocksBaseVector<T>::buildVector( backend_ptrtype backend )
{
    M_vector = backend->newBlockVector( _block=*this );
}

template<typename T>
typename BlocksBaseVector<T>::vector_ptrtype&
BlocksBaseVector<T>::vectorMonolithic()
{
    boost::shared_ptr< VectorBlockBase<T> > vcast = boost::dynamic_pointer_cast< VectorBlockBase<T> >( M_vector );
    return vcast->getVector();
}
template<typename T>
typename BlocksBaseVector<T>::vector_ptrtype const&
BlocksBaseVector<T>::vectorMonolithic() const
{
    boost::shared_ptr< VectorBlockBase<T> const > vcast = boost::dynamic_pointer_cast< VectorBlockBase<T> const >( M_vector );
    return vcast->getVector();
}

template <typename T>
void
BlocksBaseVector<T>::setVector( vector_type & vec, vector_type const& subvec, int blockId, bool closeVector ) const
{
    auto const& dmVec = vec.map();
    auto const& dmSubVec = subvec.map();
    for ( int tag=0 ; tag<dmSubVec.numberOfDofIdToContainerId() ; ++tag )
    {
        auto const& basisGpToContainerGpSubVec = dmSubVec.dofIdToContainerId( tag );
        CHECK( blockId+tag < dmVec.numberOfDofIdToContainerId() ) << "error "<<blockId+tag << " vs " << dmVec.numberOfDofIdToContainerId();
        auto const& basisGpToContainerGpVec = dmVec.dofIdToContainerId( blockId+tag );
        CHECK( basisGpToContainerGpSubVec.size() == basisGpToContainerGpVec.size() ) << " aii " << basisGpToContainerGpSubVec.size() << " vs " << basisGpToContainerGpVec.size();
        for ( int k=0;k<basisGpToContainerGpSubVec.size();++k )
            vec.set( basisGpToContainerGpVec[k], subvec( basisGpToContainerGpSubVec[k] ) );
            //vec( basisGpToContainerGpVec[k] ) = subvec( basisGpToContainerGpSubVec[k] );
    }
    if ( closeVector )
        vec.close();
}

template <typename T>
void
BlocksBaseVector<T>::setSubVector( vector_type & subvec, vector_type const& vec , int idStart ) const
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

template <typename T>
void
BlocksBaseVector<T>::updateVectorFromSubVectors()
{
    if ( !M_vector )
        return;
    int nBlock = this->nRow();
    for ( int k = 0 ; k<nBlock ;++k )
        this->setVector( *M_vector, *(this->operator()(k)), k, false );
    M_vector->close();
}


template class BlocksBaseVector<double>;


template <typename T>
VectorBlockBase<T>::VectorBlockBase( vf::BlocksBase<vector_ptrtype> const & blockVec,
                                     backend_type &backend,
                                     bool copy_values )
    :
    M_vec()
{
    auto nRow = blockVec.nRow();

    boost::shared_ptr<DataMap> dm;
    if ( nRow == 1 )
    {
        dm = blockVec(0,0)->mapPtr();
    }
    else
    {
        std::vector<boost::shared_ptr<DataMap> > listofdm;
        for ( uint16_type i=0 ; i<nRow; ++i )
            listofdm.push_back( blockVec(i,0)->mapPtr() );
        dm.reset( new DataMap( listofdm, blockVec(0,0)->map().worldComm() ) );
    }
    M_vec = backend.newVector( dm );

    M_vec->zero();

    if ( copy_values )
    {
        size_type start_i = 0;//M_vec->map().firstDof();
        for ( int i=0; i<nRow; ++i )
        {
            blockVec( i,0 )->close(); // not good but necessary here (TODO)
            this->updateBlockVec( blockVec( i,0 ), start_i );
            start_i += blockVec( i,0 )->map().numberOfDofIdToContainerId();
            //start_i += blockVec( i,0 )->map().nLocalDofWithGhost();
        }
    }
    this->setMap( M_vec->mapPtr() );
}

template <typename T>
void
VectorBlockBase<T>::updateBlockVec( vector_ptrtype const& m, size_type start_i )
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
    this->setMap( M_vec->mapPtr() );
}

template<typename T>
void
VectorBlockBase<T>::addVector ( int* rows, int nrows, value_type* data, size_type K, size_type K2 )
{
    M_vec->addVector( rows, nrows, data, K, K2 );
}
template class VectorBlockBase<double>;

} // Feel


