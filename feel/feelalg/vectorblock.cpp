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
    vb->close();

    size_type _start_iloc=0;
    for ( uint16_type i=0; i<this->nRow(); ++i )
    {
        size_type nBlockRow = this->operator()( i,0 )->localSize();

        if ( !this->vector() )
        {
            for ( size_type k=0; k<nBlockRow; ++k )
            {
                this->operator()( i,0 )->set( k, vb->operator()( _start_i+k ) );
            }
        }
        else
        {
            for ( size_type k=0; k<nBlockRow; ++k )
            {
                const T val = vb->operator()( _start_i+k );
                this->operator()( i,0 )->set( k, val );
                this->vector()->set( _start_iloc+k, val );
            }
        }

        this->operator()( i,0 )->close();
        _start_i += nBlockRow;
        _start_iloc += nBlockRow;
    }
}

template <typename T>
void
BlocksBaseVector<T>::localize()
{
    if ( !this->vector() ) return;

    this->vector()->close();
#if 0
    size_type _start_i=0;
    for ( uint16_type i=0; i<this->nRow(); ++i )
    {
        size_type nBlockRow = this->operator()( i,0 )->localSize();

            for ( size_type k=0; k<nBlockRow; ++k )
            {
                this->operator()( i,0 )->set( k, this->vector()->operator()( _start_i+k ) );
            }

        this->operator()( i,0 )->close();
        _start_i += nBlockRow;
    }
#else
    auto const& dm = this->vector()->map();
    int currentDataBaseId = 0;
    for ( uint16_type i=0; i<this->nRow(); ++i )
    {
        //size_type nBlockRow = this->operator()( i,0 )->localSize();
        auto const& dmb = this->operator()( i,0 )->map();
        int nDataBase = dmb.nBasisGp();
        for ( int tag=0;tag<nDataBase;++tag,++currentDataBaseId )
        {
            auto const& basisGpToCompositeGpBlock = dmb.basisGpToCompositeGp(tag);
            auto const& basisGpToCompositeGpVec = dm.basisGpToCompositeGp(currentDataBaseId);
            for (int k=0;k<basisGpToCompositeGpBlock.size();++k)
            {
                this->operator()( i,0 )->set( basisGpToCompositeGpBlock[k], this->vector()->operator()( basisGpToCompositeGpVec[k] ) );
            }
        }
        this->operator()( i,0 )->close();
    }


#endif
}

template <typename T>
void
BlocksBaseVector<T>::buildVector( backend_ptrtype backend )
{
    M_vector = backend->newBlockVector( _block=*this );
}


template <typename T>
void
BlocksBaseVector<T>::setVector( vector_type & vec, vector_type const& subvec , int blockId ) const
{
    auto const& dmVec = vec.map();
    auto const& dmSubVec = subvec.map();
    for ( int tag=0 ; tag<dmSubVec.nBasisGp() ; ++tag )
    {
        auto const& basisGpToContainerGpSubVec = dmSubVec.basisGpToCompositeGp( tag );
        CHECK( blockId+tag < dmVec.nBasisGp() ) << "error "<<blockId+tag << " vs " << dmVec.nBasisGp();
        auto const& basisGpToContainerGpVec = dmVec.basisGpToCompositeGp( blockId+tag );
        CHECK( basisGpToContainerGpSubVec.size() == basisGpToContainerGpVec.size() ) << " aii " << basisGpToContainerGpSubVec.size() << " vs " << basisGpToContainerGpVec.size();
        for ( int k=0;k<basisGpToContainerGpSubVec.size();++k )
            vec( basisGpToContainerGpVec[k] ) = subvec( basisGpToContainerGpSubVec[k] );
    }
}

template <typename T>
void
BlocksBaseVector<T>::setSubVector( vector_type & subvec, vector_type const& vec , int idStart ) const
{
    auto const& dmVec = vec.map();
    auto const& dmSubVec = subvec.map();
    int basisIndexSubVec = idStart;
    for ( int tag=0 ; tag<dmSubVec.nBasisGp() ; ++tag )
    {
        auto const& basisGpToContainerGpSubVec = dmSubVec.basisGpToCompositeGp( tag );
        CHECK( basisIndexSubVec+tag < dmVec.nBasisGp() ) << "error "<<basisIndexSubVec+tag<< " vs " << dmVec.nBasisGp();
        auto const& basisGpToContainerGpVec = dmVec.basisGpToCompositeGp( basisIndexSubVec+tag );
        CHECK( basisGpToContainerGpSubVec.size() == basisGpToContainerGpVec.size() ) << " error " << basisGpToContainerGpSubVec.size() << " vs " << basisGpToContainerGpVec.size();
        for ( int k=0;k<basisGpToContainerGpSubVec.size();++k )
            subvec( basisGpToContainerGpSubVec[k] ) = vec( basisGpToContainerGpVec[k] );
    }
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
            start_i += blockVec( i,0 )->map().nBasisGp();
            //start_i += blockVec( i,0 )->map().nLocalDofWithGhost();
        }
    }
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

    int nTag = dmb.nBasisGp();
    for ( int tag=0;tag<nTag;++tag )
    {
        auto const& basisGpToCompositeGpBlock = dmb.basisGpToCompositeGp(tag);
        auto const& basisGpToCompositeGpVec = dm.basisGpToCompositeGp(start_i+tag);
        CHECK( basisGpToCompositeGpBlock.size() == basisGpToCompositeGpVec.size() ) << "incompatibility with size";
        for (int k=0;k<basisGpToCompositeGpBlock.size();++k)
            M_vec->set( basisGpToCompositeGpVec[k],m->operator()( basisGpToCompositeGpBlock[k] ) );
    }
}

template class VectorBlockBase<double>;

} // Feel


