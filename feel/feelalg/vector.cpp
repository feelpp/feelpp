/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-03-20

  Copyright (C) 2008, 2009 Universite Joseph Fourier (Grenoble I)

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
   \file vector.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-03-20
 */
#include <feel/feelalg/matrixshell.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel
{
template <typename T>
Vector<T>::Vector( WorldComm const& _worldComm )
    : M_is_closed( false ),
      M_is_initialized( false ),
      M_map( new datamap_type( _worldComm ) )
{
}

/*template <typename T>
Vector<T>::Vector ( datamap_type const& dm ) :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map ( new datamap_type(dm) )
{}
    */
template <typename T>
Vector<T>::Vector( datamap_ptrtype const& dm )
    : M_is_closed( false ),
      M_is_initialized( false ),
      M_map( dm )
{
}

template <typename T>
Vector<T>::Vector( const size_type n, WorldComm const& _worldComm )
    : M_is_closed( false ),
      M_is_initialized( false ),
      M_map( new datamap_type( n, n, _worldComm ) )
{
}

template <typename T>
Vector<T>::Vector( const size_type n,
                   const size_type n_local,
                   WorldComm const& _worldComm )
    : M_is_closed( false ),
      M_is_initialized( false ),
      M_map( new datamap_type( n, n_local, _worldComm ) )
{
}

template <typename T>
Vector<T>::Vector( Vector const& v )
    : M_is_closed( v.M_is_closed ),
      M_is_initialized( v.M_is_initialized ),
      M_map( v.M_map )
{
}
template <typename T>

Vector<T>::~Vector()
{
    clear();
}

template <typename T>

Vector<T>& Vector<T>::operator=( const T )
{
    //  error();

    return *this;
}

template <typename T>
void Vector<T>::init( const size_type n,
                      const size_type nl,
                      const bool fast )
{
    boost::ignore_unused_variable_warning( fast );
    //deja fait dans le constructeur!!
    //M_map=DataMap( n, nl );
}

template <typename T>
void Vector<T>::init( const size_type n,
                      const bool fast )
{
    this->init( n, n, fast );
}

template <typename T>

Vector<T>& Vector<T>::operator=( const Vector<T>& v )
{
    if ( this != &v )
    {
        if ( !M_map->isCompatible( v.map() ) )
            M_map = v.mapPtr();

        for ( size_type i = 0; i < this->map().nLocalDofWithGhost(); ++i )
        {
            this->set( i, v( v.firstLocalIndex() + i ) );
        }
        this->close();
    }

    return *this;
}

template <typename T>

Vector<T>& Vector<T>::operator=( const std::vector<T>& v )
{
    M_map.reset( new datamap_type( v.size(), 0 ) );
    return *this;
}

template <typename T>

void Vector<T>::clear()
{
    M_is_closed = false;
    M_is_initialized = false;
}
template <typename T>
void Vector<T>::addVector( const Vector<T>& V_in,
                           const MatrixShell<T>& A_in )
{
    A_in.multVector( V_in, *this );
}

template <typename T>
void Vector<T>::addVector( const boost::shared_ptr<Vector<T>>& V_in,
                           const boost::shared_ptr<MatrixShell<T>>& A_in )
{
    A_in->multVector( *V_in, *this );
}

template <typename T>
int Vector<T>::reciprocal()
{
    LOG( WARNING ) << "Invalid call to reciprocal. Not implement in Vector base class";
    return 0;
}

#if 0
// Full specialization of the print() member for complex
// variables.  This must precede the non-specialized
// version, at least according to icc v7.1
template <>

void Vector<Complex>::print( std::ostream& os ) const
{
    assert ( this->initialized() );
    os << "Size\tglobal =  " << this->size()
       << "\t\tlocal =  " << this->local_size() << std::endl;

    // std::complex<>::operator<<() is defined, but use this form
    os << "#\tReal part\t\tImaginary part" << std::endl;

    for ( size_type i=this->first_local_index(); i<this->last_local_index(); i++ )
        os << i << "\t"
           << ( *this )( i ).real() << "\t\t"
           << ( *this )( i ).imag() << std::endl;
}

#endif

template <>
int Vector<float>::compare( const Vector<float>& other_vector,
                            const real_type threshold ) const
{
    FEELPP_ASSERT( this->isInitialized() )
        .error( "vector not initialized" );
    FEELPP_ASSERT( other_vector.isInitialized() )
        .error( "vector not initialized" );
    FEELPP_ASSERT( this->firstLocalIndex() == other_vector.firstLocalIndex() )
        .error( "" );
    FEELPP_ASSERT( this->lastLocalIndex() == other_vector.lastLocalIndex() )
        .error( "" );

    int rvalue = -1;
    size_type i = firstLocalIndex();

    do
    {
        if ( std::abs( ( *this )(i)-other_vector( i ) ) > threshold )
            rvalue = i;

        else
            i++;
    } while ( rvalue == -1 && i < lastLocalIndex() );

    return rvalue;
}

// Full specialization for double datatypes
template <>
int Vector<double>::compare( const Vector<double>& other_vector,
                             const real_type threshold ) const
{
    FEELPP_ASSERT( this->isInitialized() )
        .error( "vector not initialized" );
    FEELPP_ASSERT( other_vector.isInitialized() )
        .error( "vector not initialized" );
    FEELPP_ASSERT( this->firstLocalIndex() == other_vector.firstLocalIndex() )
        .error( "" );
    FEELPP_ASSERT( this->lastLocalIndex() == other_vector.lastLocalIndex() )
        .error( "" );

    int rvalue = -1;
    size_type i = firstLocalIndex();

    do
    {
        if ( std::abs( ( *this )(i)-other_vector( i ) ) > threshold )
            rvalue = i;

        else
            i++;
    } while ( rvalue == -1 && i < lastLocalIndex() );

    return rvalue;
}

#if 0
// Full specialization for long double datatypes
template <>
int
Vector<long double>::compare ( const Vector<long double> &other_vector,
                               const real_type threshold ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
    FEELPP_ASSERT ( other_vector.isInitialized() ).error( "vector not initialized" );
    FEELPP_ASSERT ( this->firstLocalIndex() == other_vector.firstLocalIndex() ).error( "" );
    FEELPP_ASSERT ( this->lastLocalIndex()  == other_vector.lastLocalIndex() ).error( "" );

    int rvalue     = -1;
    size_type i = firstLocalIndex();

    do
    {
        if ( std::abs( ( *this )( i ) - other_vector( i ) ) > threshold )
            rvalue = i;

        else
            i++;
    }
    while ( rvalue==-1 && i<lastLocalIndex() );

    return rvalue;
}
#endif
// Full specialization for Complex datatypes
template <>
int Vector<std::complex<double>>::compare( const Vector<std::complex<double>>& other_vector,
                                           const real_type threshold ) const
{
    CHECK( this->isInitialized() ) << "vector not initialized";
    CHECK( other_vector.isInitialized() ) << "other vector not initialized";
    CHECK( this->firstLocalIndex() == other_vector.firstLocalIndex() ) << "invalid index";
    CHECK( this->lastLocalIndex() == other_vector.lastLocalIndex() ) << "invalid index";

    int rvalue = -1;
    size_type i = firstLocalIndex();

    do
    {
        if ( ( std::abs( ( *this )( i ).real() - other_vector( i ).real() ) > threshold ) ||
             ( std::abs( ( *this )( i ).imag() - other_vector( i ).imag() ) > threshold ) )
            rvalue = i;

        else
            i++;
    } while ( rvalue == -1 && i < this->lastLocalIndex() );

    return rvalue;
}

template <typename T>

void Vector<T>::print( std::ostream& os ) const
{
    FEELPP_ASSERT( this->isInitialized() )
        .error( "vector not initialized" );
    os << "Size\tglobal =  " << this->size()
       << "\t\tlocal =  " << this->localSize() << std::endl;

    os << "#\tValue" << std::endl;

    for ( size_type i = this->firstLocalIndex(); i < this->lastLocalIndex(); i++ )
        os << i << "\t" << ( *this )( i ) << std::endl;
}
template <typename T>
void Vector<T>::localize( Vector<T> const& v )
{
}

template class Vector<double>;
template class Vector<std::complex<double>>;
//template class Vector<long double>;

namespace detail
{
template <typename T>
struct syncOperatorEqual : syncOperator<T>
{
    typedef syncOperator<T> super_type;
    typedef typename super_type::storage_ghostdof_type storage_ghostdof_type;
    syncOperatorEqual( bool hasOperator )
        : super_type(),
          M_hasOperator( hasOperator )
    {
    }
    syncOperatorEqual( std::map<size_type, std::set<rank_type>> const& m )
        : super_type( m ),
          M_hasOperator( true )
    {
    }
    virtual T operator()( size_type gcdof, rank_type activeProcId, T activeDofValue, storage_ghostdof_type const& ghostDofs ) const
    {
        auto itFindGCDof = this->activeDofClusterUsedByProc().find( gcdof );
        if ( itFindGCDof != this->activeDofClusterUsedByProc().end() )
        {
            if ( itFindGCDof->second.find( activeProcId ) != itFindGCDof->second.end() )
                return activeDofValue;

            for ( auto const& ghostDof : ghostDofs )
            {
                rank_type ghostProcId = ghostDof.first;
                if ( itFindGCDof->second.find( ghostProcId ) != itFindGCDof->second.end() )
                {
                    T dofValue = ghostDof.second;
                    return dofValue;
                }
            }
        }
        return activeDofValue;
    }
    virtual bool hasOperator() const { return M_hasOperator; }
  private:
    bool M_hasOperator;
};

template <typename T>
struct syncOperatorPlus : syncOperator<T>
{
    typedef syncOperator<T> super_type;
    typedef typename super_type::storage_ghostdof_type storage_ghostdof_type;
    syncOperatorPlus()
        : super_type()
    {
    }
    syncOperatorPlus( std::map<size_type, std::set<rank_type>> const& m )
        : super_type( m )
    {
    }
    virtual T operator()( size_type gcdof, rank_type activeProcId, T activeDofValue, storage_ghostdof_type const& ghostDofs ) const
    {
        T res = 0;
        auto itFindGCDof = this->activeDofClusterUsedByProc().find( gcdof );
        if ( itFindGCDof == this->activeDofClusterUsedByProc().end() )
        {
            res += activeDofValue;
            for ( auto const& ghostDof : ghostDofs )
            {
                T dofValue = ghostDof.second;
                res += dofValue;
            }
            return res;
        }
        else
        {
            if ( itFindGCDof->second.find( activeProcId ) != itFindGCDof->second.end() )
                res += activeDofValue;
            for ( auto const& ghostDof : ghostDofs )
            {
                rank_type ghostProcId = ghostDof.first;
                if ( itFindGCDof->second.find( ghostProcId ) != itFindGCDof->second.end() )
                {
                    T dofValue = ghostDof.second;
                    res += dofValue;
                }
            }
            return res;
        }
    }
    virtual bool hasOperator() const { return true; }
};
// BinaryFuncType = 0 -> min, BinaryFuncType=1 -> max
template <typename T, int BinaryFuncType>
struct syncOperatorBinaryFunc : syncOperator<T>
{
    typedef syncOperator<T> super_type;
    typedef typename super_type::storage_ghostdof_type storage_ghostdof_type;
    syncOperatorBinaryFunc()
        : super_type()
    {
    }
    syncOperatorBinaryFunc( std::map<size_type, std::set<rank_type>> const& m )
        : super_type( m )
    {
    }
    virtual T operator()( size_type gcdof, rank_type activeProcId, T activeDofValue, storage_ghostdof_type const& ghostDofs ) const
    {
        auto itFindGCDof = this->activeDofClusterUsedByProc().find( gcdof );
        if ( itFindGCDof == this->activeDofClusterUsedByProc().end() )
        {
            T res = activeDofValue;
            for ( auto const& ghostDof : ghostDofs )
            {
                T dofValue = ghostDof.second;
                res = this->applyBinaryFunc( res, dofValue, mpl::int_<BinaryFuncType>() );
            }
            return res;
        }
        else
        {
            T res = 0;
            bool init = false;
            if ( itFindGCDof->second.find( activeProcId ) != itFindGCDof->second.end() )
            {
                res = activeDofValue;
                init = true;
            }
            for ( auto const& ghostDof : ghostDofs )
            {
                rank_type ghostProcId = ghostDof.first;
                if ( itFindGCDof->second.find( ghostProcId ) != itFindGCDof->second.end() )
                {
                    T dofValue = ghostDof.second;
                    if ( !init )
                    {
                        res = dofValue;
                        init = true;
                    }
                    else
                        res = this->applyBinaryFunc( res, dofValue, mpl::int_<BinaryFuncType>() );
                }
            }
            return res;
        }
    }

    static T applyBinaryFunc( T const& val1, T const& val2, mpl::int_<0> /**/ ) { return std::min( val1, val2 ); }
    static T applyBinaryFunc( T const& val1, T const& val2, mpl::int_<1> /**/ ) { return std::max( val1, val2 ); }

    virtual bool hasOperator() const { return true; }
};
}

template <typename T>
void sync( Vector<T>& v, std::string const& opSyncStr )
{
    if ( opSyncStr == "=" )
        sync( v, detail::syncOperatorEqual<T>( false ) );
    else if ( opSyncStr == "+" )
        sync( v, detail::syncOperatorPlus<T>() );
    else if ( opSyncStr == "min" )
        sync( v, detail::syncOperatorBinaryFunc<T, 0>() );
    else if ( opSyncStr == "max" )
        sync( v, detail::syncOperatorBinaryFunc<T, 1>() );
}

template <typename T>
void sync( Vector<T>& v, std::string const& opSyncStr, std::set<size_type> const& dofGlobalProcessPresent )
{
    auto activeDofData = v.mapPtr()->activeDofClusterUsedByProc( dofGlobalProcessPresent );
    if ( opSyncStr == "=" )
        sync( v, detail::syncOperatorEqual<T>( activeDofData ) );
    else if ( opSyncStr == "+" )
        sync( v, detail::syncOperatorPlus<T>( activeDofData ) );
    else if ( opSyncStr == "min" )
        sync( v, detail::syncOperatorBinaryFunc<T, 0>( activeDofData ) );
    else if ( opSyncStr == "max" )
        sync( v, detail::syncOperatorBinaryFunc<T, 1>( activeDofData ) );
}

template <typename T>
void sync( Vector<T>& v, detail::syncOperator<T> const& opSync )
{
    auto dataMap = v.mapPtr();

    // if sequential return
    if ( !dataMap || dataMap->worldComm().localSize() == 1 )
        return;

    // prepare mpi com
    rank_type currentProcId = dataMap->worldComm().localRank();
    int nbRequest = 2 * dataMap->neighborSubdomains().size();
    mpi::request* reqs = new mpi::request[nbRequest];
    // apply isend/irecv
    int cptRequest = 0;

    if ( opSync.hasOperator() )
    {
        // init data to send : values of ghost dof is send to the unique active dof
        std::map<rank_type, std::vector<boost::tuple<size_type, T>>> dataToSend, dataToRecv;
        for ( rank_type p : dataMap->neighborSubdomains() )
            dataToSend[p].clear();
        for ( size_type id = 0; id < dataMap->nLocalDofWithGhost(); ++id )
        {
            if ( dataMap->dofGlobalProcessIsGhost( id ) )
            {
                size_type gcdof = dataMap->mapGlobalProcessToGlobalCluster( id );
                T val = v( id );
                rank_type procIdFinded = dataMap->procOnGlobalCluster( gcdof );
                CHECK( procIdFinded != invalid_rank_type_value ) << " proc not find for gcdof : " << gcdof;
                dataToSend[procIdFinded].push_back( boost::make_tuple( gcdof, val ) );
            }
        }
        // apply isend/irecv
        cptRequest = 0;
        for ( rank_type p : dataMap->neighborSubdomains() )
        {
            CHECK( dataToSend.find( p ) != dataToSend.end() ) << " no data to send to proc " << p << "\n";
            reqs[cptRequest++] = dataMap->worldComm().localComm().isend( p, 0, dataToSend.find( p )->second );
            reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( p, 0, dataToRecv[p] );
        }
        // wait all requests
        mpi::wait_all( reqs, reqs + nbRequest );
        // update value of active dofs (repect to the sync operator)
        std::map<size_type, std::set<std::pair<rank_type, T>>> ghostDofValues;
        for ( auto const& dataR : dataToRecv )
        {
            rank_type theproc = dataR.first;
            for ( auto const& dataRfromproc : dataR.second )
            {
                size_type gcdof = boost::get<0>( dataRfromproc );
                T valRecv = boost::get<1>( dataRfromproc );
                ghostDofValues[gcdof].insert( std::make_pair( theproc, valRecv ) );
            }
        }
        for ( auto const& ghostDofVal : ghostDofValues )
        {
            size_type gcdof = ghostDofVal.first;
            size_type gpdof = gcdof - dataMap->firstDofGlobalCluster();
#if !defined( NDEBUG )
            auto resSearchDof = dataMap->searchGlobalProcessDof( gcdof );
            CHECK( boost::get<0>( resSearchDof ) ) << "dof not found";
            size_type gpdof2 = boost::get<1>( resSearchDof );
            CHECK( gpdof == gpdof2 ) << "dof id must be the same : " << gpdof << " vs " << gpdof2;
            CHECK( dataMap->activeDofSharedOnCluster().find( gpdof ) != dataMap->activeDofSharedOnCluster().end() ) << "not an active dof shared";
#endif
            T valCurrent = v( gpdof );
            v.set( gpdof, opSync( gcdof, currentProcId, valCurrent, ghostDofVal.second ) );
        }
    }

    // init data to re-send : update values of active dof in ghost dof associated
    std::map<rank_type, std::vector<boost::tuple<size_type, T>>> dataToReSend, dataToReRecv;
    for ( rank_type p : dataMap->neighborSubdomains() )
        dataToReSend[p].clear();
    for ( auto const& dofActive : dataMap->activeDofSharedOnCluster() )
    {
        size_type gpdof = dofActive.first;
        CHECK( !dataMap->dofGlobalProcessIsGhost( gpdof ) ) << "must be active";

        T val = v( gpdof );
        size_type gcdof = dataMap->mapGlobalProcessToGlobalCluster( gpdof );
        for ( rank_type pNeighborId : dofActive.second )
        {
            if ( pNeighborId != dataMap->worldComm().localRank() ) // normaly this check is useless
                dataToReSend[pNeighborId].push_back( boost::make_tuple( gcdof, val ) );
        }
    }
    // prepare mpi com
    cptRequest = 0;
    for ( rank_type p : dataMap->neighborSubdomains() )
    {
        CHECK( dataToReSend.find( p ) != dataToReSend.end() ) << " no data to send to proc " << p << "\n";
        reqs[cptRequest++] = dataMap->worldComm().localComm().isend( p, 0, dataToReSend.find( p )->second );
        reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( p, 0, dataToReRecv[p] );
    }
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    delete[] reqs;
    // update values of ghost dofs
    for ( auto const& dataR : dataToReRecv )
    {
        rank_type theproc = dataR.first;
        for ( auto const& dataRfromproc : dataR.second )
        {
            size_type gcdof = boost::get<0>( dataRfromproc );
            T valRecv = boost::get<1>( dataRfromproc );
            auto resSearchDof = dataMap->searchGlobalProcessDof( gcdof );
            DCHECK( boost::get<0>( resSearchDof ) ) << "dof not found";
            size_type gpdof = boost::get<1>( resSearchDof );
            DCHECK( dataMap->dofGlobalProcessIsGhost( gpdof ) ) << "dof is not ghost : " << gcdof << " and " << gpdof;
            v.set( gpdof, valRecv );
        }
    }
}

template void sync<double>( Vector<double>& v, std::string const& opSyncStr, std::set<size_type> const& dofGlobalProcessPresent );
template void sync<double>( Vector<double>& v, std::string const& opSyncStr );
template void sync<double>( Vector<double>& v, detail::syncOperator<double> const& opSync );
}
