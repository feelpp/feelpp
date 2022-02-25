/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2021-10-29

  Copyright (C) 2021 Feel++ Consortium
  Copyright (C) 2021 INRIA

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

#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelcore/hdf5.hpp>

#define FEELPP_INSTANTIATE_VECTORUBLAS 1

namespace Feel {

namespace detail {

template< typename T >
VectorUblasBase<T>::VectorUblasBase( size_type s ):
    super_type( s, Environment::worldCommSeqPtr() )
{
    this->init( s, s, false ); 
}

template< typename T >
VectorUblasBase<T>::VectorUblasBase( const datamap_ptrtype & dm ):
    super_type( dm )
{
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false ); 
}

template< typename T >
VectorUblasBase<T>::VectorUblasBase( size_type s, size_type n_local ):
    super_type( s, n_local, Environment::worldCommSeqPtr() )
{
    this->init( this->size(), this->localSize(), false ); 
}

template< typename T >
void VectorUblasBase<T>::init( const size_type n, const size_type n_local, const bool fast )
{}

template< typename T >
void VectorUblasBase<T>::init( const size_type n, const bool fast )
{
    this->init( n, n, fast ); 
}

template< typename T >
void VectorUblasBase<T>::init( const datamap_ptrtype & dm )
{
    super_type::init( dm );
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false ); 
}

template< typename T >
Vector<T> & VectorUblasBase<T>::operator=( const Vector<T> & v )
{
    if ( this == &v )
        return *this;

    if ( !this->map().isCompatible( v.map() ) )
    {
        this->setMap( v.mapPtr() );
        this->resize( this->map().nLocalDofWithGhost() );
    }

#if !defined(NDEBUG)
    this->checkInvariants();
    FEELPP_ASSERT( this->localSize() == v.localSize() &&
                   this->map().nLocalDofWithoutGhost() == v.map().nLocalDofWithoutGhost() )
    ( this->localSize() )( this->map().nLocalDofWithoutGhost() )
    ( this->size() )
    ( v.localSize() )( v.map().nLocalDofWithoutGhost() ).warn( "may be vector invalid copy" );
#endif

    this->set( v );
    return *this;
}

template< typename T >
VectorUblasBase<T> & VectorUblasBase<T>::operator=( const VectorUblasBase<T> & v )
{
    if ( this == &v )
        return *this;

    if ( !this->map().isCompatible( v.map() ) )
    {
        this->setMap( v.mapPtr() );
        this->resize( this->map().nLocalDofWithGhost() );
    }

#if !defined(NDEBUG)
    this->checkInvariants();
    FEELPP_ASSERT( this->localSize() == v.localSize() &&
                   this->map().nLocalDofWithoutGhost() == v.map().nLocalDofWithoutGhost() )
    ( this->localSize() )( this->map().nLocalDofWithoutGhost() )
    ( this->size() )
    ( v.localSize() )( v.map().nLocalDofWithoutGhost() ).warn( "may be vector invalid copy" );
#endif

    this->setVector( v );
    return *this;
}

template< typename T >
void VectorUblasBase<T>::set( const size_type i, const value_type & value )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    this->operator()( i ) = value;
}

template< typename T >
void VectorUblasBase<T>::set( const Vector<T> & v )
{
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->setVector( *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->setVector( *vecPetsc );
#endif
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) = v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::setVector( int * i, int n, value_type * v )
{
    for( int j = 0; j < n; ++j )
        this->operator()( i[j] ) = v[j];
}

template< typename T >
void VectorUblasBase<T>::add( const size_type i, const value_type & value )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    this->operator()( i ) += value;
}

template< typename T >
void VectorUblasBase<T>::add( const value_type & value )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += value;
}

template< typename T >
void VectorUblasBase<T>::add( const Vector<T> & v )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->addVector( *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->addVector( *vecPetsc );
#endif
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::add( const value_type & a, const Vector<T> & v )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->maddVector( a, *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->maddVector( a, *vecPetsc );
#endif
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += a * v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::sub( const Vector<T> & v )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->subVector( *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->subVector( *vecPetsc );
#endif
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) -= v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::sub( const value_type & a, const Vector<T> & v )
{
#ifndef NDEBUG
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->msubVector( a, *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->msubVector( a, *vecPetsc );
#endif
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) -= a * v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::addVector( int * i, int n, value_type * v, size_type K, size_type K2 )
{
    for( int j = 0; j < n; ++j )
        this->operator()( i[j] ) += v[j];
}

template< typename T >
void VectorUblasBase<T>::addVector( const std::vector<value_type> & v, const std::vector<size_type> & dof_ids )
{
    FEELPP_ASSERT( v.size() == dof_ids.size() ).error( "invalid dof indices" );

    for ( size_type j = 0; j < v.size(); ++j )
        this->operator()( dof_ids[j] ) += v[j];
}

template< typename T >
void VectorUblasBase<T>::addVector( const Vector<value_type> & v, const std::vector<size_type> & dof_ids )
{
    FEELPP_ASSERT( v.size() == dof_ids.size() ).error( "invalid dof indices" );

    for ( size_type j = 0; j < v.size(); ++j )
        this->operator()( dof_ids[j] ) += v(j);
}

template< typename T >
void VectorUblasBase<T>::addVector( const ublas::vector<value_type> & v, const std::vector<size_type> & dof_ids )
{
    FEELPP_ASSERT( v.size() == dof_ids.size() ).error( "invalid dof indices" );

    for ( size_type j = 0; j < v.size(); ++j )
        this->operator()( dof_ids[j] ) += v(j);
}

template< typename T >
typename VectorUblasBase<T>::value_type VectorUblasBase<T>::dot( const Vector<T> & v ) const
{
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->dotVector( *vecUblas );
#if FEELPP_HAS_PETSC
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->dotVector( *vecPetsc );
#endif
    // Default
    value_type localRes = 0;
    for( size_type i = 0; i < this->map().nLocalDofWithoutGhost(); ++i )
        localRes += this->operator()( i ) * v( i );
    value_type globalRes = localRes;
#ifdef FEELPP_HAS_MPI
    if( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localRes, globalRes, std::plus<value_type>() );
#endif
    return globalRes;
}

namespace helpers {
#if defined(FEELPP_HAS_TBB)
template<typename VectorType>
struct Sqrt
{
    VectorType M_in;
    VectorType& M_out;
    Sqrt( VectorType const& _in, VectorType& _out )
        :
        M_in( _in ), M_out( _out )
    { }
    void operator() ( const tbb::blocked_range<size_t>& r ) const
    {
        for ( size_t i = r.begin(); i != r.end(); ++i )
        {
            M_out[i] = math::sqrt( M_in[i] );
        }
    }
};
#endif // FEELPP_HAS_TBB
}

template< typename T >
std::unique_ptr<VectorUblasBase<T>>
VectorUblasBase<T>::sqrt() const
{
    std::unique_ptr<VectorUblasBase<T>> sqrtV( this->emptyPtr() );
    sqrtV->init( this->mapPtr() );

#if defined( FEELPP_HAS_TBB )
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, this->localSize() ),
                       helpers::Sqrt<this_type>( *this, *sqrtV ) );
#else
    for ( size_type i = 0; i < this->localSize(); ++i )
        sqrtV->operator()(i) = math::sqrt( this->operator()( i ) );
#endif // FEELPP_HAS_TBB

    return sqrtV;
}

template< typename T >
std::unique_ptr<VectorUblasBase<T>>
VectorUblasBase<T>::pow( int n ) const
{
    std::unique_ptr<VectorUblasBase<T>> powV( this->emptyPtr() );
    powV->init( this->mapPtr() );

    for ( size_type i = 0; i < this->localSize(); ++i )
        powV->operator()(i) = math::pow( this->operator()( i ), n );

    return powV;
}

template< typename T >
void VectorUblasBase<T>::printMatlab( const std::string filename, bool renumber ) const
{
    std::string name = filename;
    std::string separator = " , ";

    // check on the file name
    int i = filename.find( "." );

    if ( i <= 0 )
        name = filename + ".m";

    else
    {
        if ( ( unsigned int ) i != filename.size() - 2 ||
                filename[ i + 1 ] != 'm' )
        {
            DVLOG(2) << "[VectorUblas::printMatlab] adding .m extension to given file name '"
                          << filename << "'\n";
            name = filename + ".m";
        }
    }


    ublas::vector<value_type> v_local;
    this->localizeToOneProcessor( v_local, 0 );

    if ( this->comm().rank() == 0 )
    {
        std::ofstream file_out( name.c_str() );

        FEELPP_ASSERT( file_out )( filename ).error( "[VectorUblas::printMatlab] ERROR: File cannot be opened for writing." );

        file_out << "var_"<<filename<<" = [ ";
        file_out.precision( 16 );
        file_out.setf( std::ios::scientific );

        for ( size_type i = 0; i < this->size(); ++i )
        {
            file_out << v_local[i] << separator << std::endl;
        }

        file_out << "];" << std::endl;
    }
}

#ifdef FEELPP_HAS_HDF5
template< typename T >
void VectorUblasBase<T>::saveHDF5( const std::string & filename, const std::string & tableName, bool appendMode ) const
{
    bool useTransposedStorage = true;
    const int dimsComp0 = (useTransposedStorage)? 1 : 0;
    const int dimsComp1 = (useTransposedStorage)? 0 : 1;
    DataMap<> const& dm = this->map();
    std::vector<uint> sizeValues(1, dm.nLocalDofWithoutGhost() );
    hsize_t dims[2];
    dims[dimsComp0] = this->comm().localSize();dims[dimsComp1] = 1;
    hsize_t dims2[2];

    dims2[dimsComp0] = sizeValues.size();dims2[dimsComp1] = 1;
    hsize_t offset[2];
    offset[dimsComp0] = this->comm().localRank(); offset[dimsComp1] = 0;

    hsize_t dimsElt[2];
    dimsElt[dimsComp0] = dm.nDof();
    dimsElt[dimsComp1] = 1;

    hsize_t dimsElt2[2];
    dimsElt2[dimsComp0] = dm.nLocalDofWithoutGhost();
    dimsElt2[dimsComp1] = 1;
    hsize_t offsetElt[2];
    size_type offsetCount = 0;
    for (rank_type p=0 ; p < this->comm().localRank() ; ++p )
        offsetCount += dm.nLocalDofWithoutGhost( p );

    offsetElt[dimsComp0] = offsetCount;
    offsetElt[dimsComp1] = 0;

    std::vector<double> dataStorage( dm.nLocalDofWithoutGhost() );

    HDF5 hdf5;

    /* If appendMode is true, then we want to append the table to the hdf5 file */
    /* To do so we open the hdf5 as existing and allow read/write in it */
    if(appendMode)
    {
        hdf5.openFile( filename, this->comm().localComm(), true, true );
    }
    else
    {
        hdf5.openFile( filename, this->comm().localComm(), false );
    }
    if ( false )
    {
        // create size tab
        hdf5.createTable( "size", H5T_NATIVE_UINT, dims );
        hdf5.write( "size", H5T_NATIVE_UINT, dims2, offset, sizeValues.data() );
        hdf5.closeTable( "size" );
    }

    for ( size_type k=0;k<dataStorage.size();++k )
        dataStorage[k] = this->operator()( k );

    // create double tab
    hdf5.createTable( tableName, H5T_NATIVE_DOUBLE, dimsElt );
    if ( dm.nLocalDofWithoutGhost() > 0 )
        hdf5.write( tableName, H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, dataStorage.data()/*&(M_vec[0])*/ );
    hdf5.closeTable( tableName );

    hdf5.closeFile();
}

template< typename T >
void VectorUblasBase<T>::loadHDF5( const std::string & filename, const std::string & tableName )
{
    bool useTransposedStorage = true;
    const int dimsComp0 = (useTransposedStorage)? 1 : 0;
    const int dimsComp1 = (useTransposedStorage)? 0 : 1;
    DataMap<> const& dm = this->map();
    if ( dm.nDof() == 0 )
        return;

    // create maybe a mpi sub comm if some process has no dof
    std::vector<rank_type> rankInGroup;
    for (rank_type p=0 ; p < this->comm().localSize() ; ++p )
        if ( dm.nLocalDofWithoutGhost( p ) > 0 )
            rankInGroup.push_back( p );
    std::shared_ptr<boost::mpi::communicator> subComm;
    if ( rankInGroup.size() != this->comm().size() )
    {
        boost::mpi::group mpiGroup = this->comm().group().include( rankInGroup.begin(), rankInGroup.end() );
        subComm = std::make_shared<boost::mpi::communicator>( this->comm(), mpiGroup );
    }

    if ( dm.nLocalDofWithoutGhost() > 0 )
    {
        hsize_t dimsElt2[2];
        dimsElt2[dimsComp0] = dm.nLocalDofWithoutGhost();
        dimsElt2[dimsComp1] = 1;
        hsize_t offsetElt[2];
        size_type offsetCount = 0;
        for (rank_type p=0 ; p < this->comm().localRank() ; ++p )
            offsetCount += dm.nLocalDofWithoutGhost( p );

        offsetElt[dimsComp0] = offsetCount;
        offsetElt[dimsComp1] = 0;

        std::vector<double> dataStorage( dm.nLocalDofWithoutGhost() );

        HDF5 hdf5;
        hdf5.openFile( filename, (subComm)? *subComm : this->comm().comm(), true );

        hsize_t dimsGlob[2];
        hdf5.openTable( tableName, dimsGlob );
        CHECK( dimsGlob[dimsComp0] == dm.nDof() ) << "invalid table dimension";
        CHECK( dimsGlob[dimsComp1] == 1 ) << "invalid table dimension";

        hdf5.read( tableName, H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, dataStorage.data()/*&(M_vec[0])*/ );

        hdf5.closeTable( tableName );

        hdf5.closeFile();

        for ( size_type k=0;k<dataStorage.size();++k )
            this->set( k, dataStorage[k] );
    }
    sync( *this );
}
#endif

template< typename T >
void VectorUblasBase<T>::localizeToOneProcessor( ublas::vector<T> & v_local, const size_type pid ) const
{
    checkInvariants();

    v_local.resize( this->size() );
    std::fill( v_local.begin(), v_local.end(), 0. );

    ublas::vector<value_type> v_tmp( this->size() );
    std::fill( v_tmp.begin(), v_tmp.end(), 0. );

    for( size_type i = 0; i < this->localSize(); ++i )
        v_tmp[i+this->firstLocalIndex()] = this->operator()( this->firstLocalIndex()+i );

#ifdef FEELPP_HAS_MPI
    if( this->comm().size() > 1 )
    {
        MPI_Reduce( &v_tmp[0], &v_local[0], v_local.size(),
                    MPI_DOUBLE, MPI_SUM, pid, this->comm() );
    }
    else
    {
        std::copy( v_tmp.begin(), v_tmp.end(), v_local.begin() );
    }
#else
    FEELPP_ASSERT( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );
    FEELPP_ASSERT( pid == 0 )( pid ).error( "invalid pid in non MPI mode" );
#endif
}

template< typename T >
void VectorUblasBase<T>::localizeToOneProcessor( std::vector<T> & v_local, const size_type pid ) const
{
    ublas::vector<T> ublasV;
    this->localizeToOneProcessor( ublasV, pid );
    v_local.resize( ublasV.size() );
    std::copy( ublasV.begin(), ublasV.end(), v_local.begin() );
}

/*-----------------------------------------------------------------------------*/

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::init( const size_type n, const size_type n_local, const bool fast ) 
{
    FEELPP_ASSERT ( n_local <= n )
    ( n_local )( n )
    ( this->comm().rank() )
    ( this->comm().size() ).error( "Invalid local vector size" );

    // Clear the data structures if already initialized
    if ( this->isInitialized() )
        this->clear();

    super_type::init( n, n_local, fast );

    // Initialize data structures
    this->resize( this->localSize() );

    // Set the initialized flag
    this->M_is_initialized = true;

    DVLOG(2) << "        global size = " << n << "\n";
    DVLOG(2) << "        global size = " << n_local << "\n";
    DVLOG(2) << "        global size = " << this->size() << "\n";
    DVLOG(2) << "        local  size = " << this->localSize() << "\n";
    DVLOG(2) << "  first local index = " << this->firstLocalIndex() << "\n";
    DVLOG(2) << "   last local index = " << this->lastLocalIndex() << "\n";

    // Zero the components unless directed otherwise
    if ( !fast )
        this->zero();
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::self_type *
VectorUblasContiguousGhosts<T, Storage>::clonePtr() const
{
    return new VectorUblasContiguousGhosts<T, Storage>( *this );
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::base_type *
VectorUblasContiguousGhosts<T, Storage>::emptyPtr() const
{
    if constexpr ( is_vector_proxy ) // should not happen
        return new VectorUblasContiguousGhosts<T, typename Storage::vector_type>( );
    else
        return new VectorUblasContiguousGhosts<T, Storage>( );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::resize( size_type n )
{
    if constexpr ( std::is_same_v<storage_type, vector_storage_type> )
        M_vec.resize( n );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::clear( )
{
    if constexpr ( std::is_same_v<storage_type, vector_storage_type> )
        M_vec.resize( 0 );
}
        
template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::value_type VectorUblasContiguousGhosts<T, Storage>::operator()( size_type i ) const
{
    return M_vec.operator()( i-this->firstLocalIndex() );
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::value_type& VectorUblasContiguousGhosts<T, Storage>::operator()( size_type i )
{
    return M_vec.operator()( i-this->firstLocalIndex() );
}

template< typename T, typename Storage > 
typename VectorUblasContiguousGhosts<T, Storage>::size_type
VectorUblasContiguousGhosts<T, Storage>::startActive() const
{
    if constexpr ( is_vector_proxy )
        return M_vec.start(); // should not happen but not dangerous
    else
        return 0;
}

template< typename T, typename Storage > 
typename VectorUblasContiguousGhosts<T, Storage>::size_type
VectorUblasContiguousGhosts<T, Storage>::startGhost() const
{
    if constexpr ( is_vector_proxy )
        CHECK( false ) << "should not happen";
    else
        return this->map().nLocalDofWithoutGhost();
}

template< typename T, typename Storage > 
void VectorUblasContiguousGhosts<T, Storage>::setConstant( value_type v )
{
    M_vec = ublas::scalar_vector<value_type>( M_vec.size(), v );
}  

template< typename T, typename Storage > 
void VectorUblasContiguousGhosts<T, Storage>::setZero()
{
    M_vec = ublas::zero_vector<value_type>( M_vec.size() );
}
  

template< typename T, typename Storage > 
void VectorUblasContiguousGhosts<T, Storage>::scale( const value_type factor )
{
    M_vec.operator*=( factor );
}
 
template< typename T, typename Storage > 
typename VectorUblasContiguousGhosts<T, Storage>::real_type VectorUblasContiguousGhosts<T, Storage>::min( bool parallel ) const
{
    checkInvariants();

    size_type nActiveDof = this->map().nLocalDofWithoutGhost();
    size_type nGhostDof = this->map().nLocalGhosts();
    real_type local_min = ( nActiveDof > 0 ) ?
        *std::min_element( M_vec.begin(), ( nGhostDof == 0 ) ? M_vec.end() : M_vec.begin()+nActiveDof ) :
        std::numeric_limits<real_type>::max();
    real_type global_min = local_min;

#ifdef FEELPP_HAS_MPI
    if ( parallel && this->comm().size() > 1 )
    {
        MPI_Allreduce( &local_min, &global_min, 1,
                MPI_DOUBLE, MPI_MIN, this->comm() );
    }
#endif
    return global_min;
}

template< typename T, typename Storage > 
typename VectorUblasContiguousGhosts<T, Storage>::real_type VectorUblasContiguousGhosts<T, Storage>::max( bool parallel ) const
{
    checkInvariants();

    size_type nActiveDof = this->map().nLocalDofWithoutGhost();
    size_type nGhostDof = this->map().nLocalGhosts();
    real_type local_max = ( nActiveDof > 0 ) ?
        *std::max_element( M_vec.begin(), ( nGhostDof == 0 ) ? M_vec.end() : M_vec.begin()+nActiveDof ) :
        std::numeric_limits<real_type>::min();
    real_type global_max = local_max;

#ifdef FEELPP_HAS_MPI
    if ( parallel && this->comm().size() > 1 )
    {
        MPI_Allreduce( &local_max, &global_max, 1,
                MPI_DOUBLE, MPI_MAX, this->comm() );
    }
#endif
    return global_max;
}
        
template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::real_type VectorUblasContiguousGhosts<T, Storage>::l1Norm() const
{
    checkInvariants();

    real_type local_l1 = 0;
    if ( this->comm().size() == 1 )
        local_l1 = ublas::norm_1( M_vec );
    else
        local_l1 = ublas::norm_1( ublas::project( M_vec, ublas::range( 0, this->map().nLocalDofWithoutGhost() ) ) );

    real_type global_l1 = local_l1;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
    {
        mpi::all_reduce( this->comm(), local_l1, global_l1, std::plus<real_type>() );
    }
#endif
    return global_l1;
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::real_type VectorUblasContiguousGhosts<T, Storage>::l2Norm() const
{
    checkInvariants();

    real_type local_norm2 = 0;
    if ( this->comm().size() == 1 )
    {
        local_norm2 = ublas::inner_prod( M_vec, M_vec );
    }
    else
    {
        auto vecActiveDof = ublas::project( M_vec, ublas::range( 0, this->map().nLocalDofWithoutGhost() ) );
        local_norm2 = ublas::inner_prod( vecActiveDof, vecActiveDof );
    }
    real_type global_norm2 = local_norm2;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), local_norm2, global_norm2, std::plus<real_type>() );
#endif

    return math::sqrt( global_norm2 );
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::real_type VectorUblasContiguousGhosts<T, Storage>::linftyNorm() const
{
    checkInvariants();

    real_type local_norminf = 0;
    if ( this->comm().size() == 1 )
    {
        local_norminf = ublas::norm_inf( M_vec );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        local_norminf = ublas::norm_inf( ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) );
    }

    real_type global_norminf = local_norminf;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
    {
        mpi::all_reduce( this->comm(), local_norminf, global_norminf, mpi::maximum<real_type>() );
    }
#endif
    return global_norminf;
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::value_type VectorUblasContiguousGhosts<T, Storage>::sum() const
{
    checkInvariants();

    value_type local_sum = 0;
    if ( this->comm().size() == 1 )
    {
        local_sum = ublas::sum( M_vec );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        local_sum = ublas::sum( ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) );
    }
    value_type global_sum = local_sum;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), local_sum, global_sum, std::plus<value_type>() );
#endif

    return global_sum;
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::checkInvariants() const
{
    DCHECK( this->isInitialized() ) <<  "vector not initialized" ;
    DCHECK( this->localSize() <= this->size() ) << "vector invalid size: " << this->size() << "," << this->localSize();
    DCHECK( this->localSize() == ( M_vec.size() ) ) << "vector invalid size: " 
        << M_vec.size() << ","
        << this->localSize();
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::setVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec = v.vec();
    std::visit( [this]( auto && _v ) { this->M_vec = *_v; }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::setVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec = v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec = *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) = v.vec();
            std::visit( [this, nLocalDofWithoutGhost]( auto && _v ) { ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) = *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) = v.vecNonContiguousGhosts();
            std::visit( [this, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) = *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::addVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec += v.vec();
    std::visit( [this]( auto && _v ) { this->M_vec += *_v; }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::addVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec += *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) += v.vec();
            std::visit( [this, nLocalDofWithoutGhost]( auto && _v ) { ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) += *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) += v.vecNonContiguousGhosts();
            std::visit( [this, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) += *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec += a * v.vec();
    std::visit( [this, a]( auto && _v ) { this->M_vec += a * (*_v); }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += a * v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec += a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) += a * v.vec();
            std::visit( [this, a, nLocalDofWithoutGhost]( auto && _v ) { ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) += a * (*_v); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) += a * v.vecNonContiguousGhosts();
            std::visit( [this, a, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) += a * (*_v); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::subVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec -= v.vec();
    std::visit( [this]( auto && _v ) { this->M_vec -= *_v; }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::subVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec -= *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) -= v.vec();
            std::visit( [this, nLocalDofWithoutGhost]( auto && _v ) { ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) -= *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) -= v.vecNonContiguousGhosts();
            std::visit( [this, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) -= *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec -= a * v.vec();
    std::visit( [this, a]( auto && _v ) { this->M_vec -= a * (*_v); }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= a * v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec -= a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) -= a * v.vec();
            std::visit( [this, a, nLocalDofWithoutGhost]( auto && _v ) { ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) -= a * (*_v); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) -= a * v.vecNonContiguousGhosts();
            std::visit( [this, a, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) -= a * (*_v); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::mulVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    //M_vec *= v.vec();
    // [noalias] 
    // We use "noalias" assignment to prevent unneeded copy. Despite the obvious M_vec aliasing, this should
    // be safe as both element_prod and assignment work on a coefficient-wise level.
    std::visit( [this]( auto && _v ) { ublas::noalias(this->M_vec) = ublas::element_prod( this->M_vec, *_v ); }, v.vec() );
}

template< typename T, typename Storage >
void VectorUblasContiguousGhosts<T, Storage>::mulVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec *= v.vec();
        // c.f. [noalias]
        std::visit( [this]( auto && _v ) { ublas::noalias(this->M_vec) = ublas::element_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ) *= v.vec();
            // c.f. [noalias]
            std::visit( [this, nLocalDofWithoutGhost]( auto && _v ) { auto activeRange = ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ); ublas::noalias( activeRange ) = ublas::element_prod( activeRange, *_v ); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //ublas::project( M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ) *= v.vecNonContiguousGhosts();
            // c.f. [noalias]
            std::visit( [this, nLocalDofWithoutGhost, nLocalGhosts]( auto && _v ) { auto ghostRange = ublas::project( this->M_vec, ublas::range( nLocalDofWithoutGhost, nLocalDofWithoutGhost+nLocalGhosts ) ); ublas::noalias( ghostRange ) = ublas::element_prod( ghostRange, *_v ); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::value_type 
VectorUblasContiguousGhosts<T, Storage>::dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const
{
    value_type localResult = 0;
    if ( this->comm().localSize() == 1 )
    {
        //localResult = ublas::inner_prod( M_vec, v.vec() );
        localResult = std::visit( [this]( auto && _v ) -> value_type { return ublas::inner_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        //localResult = ublas::inner_prod( ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ),
                                         //ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) ) );
        localResult = std::visit( [this, nLocalDofWithoutGhost, nLocalDofWithoutGhostOther]( auto && _v ) -> value_type { 
                return ublas::inner_prod( ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ),
                                          ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ) );
                }, v.vec() );
    }

    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif
    return globalResult;
}

template< typename T, typename Storage >
typename VectorUblasContiguousGhosts<T, Storage>::value_type 
VectorUblasContiguousGhosts<T, Storage>::dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const
{
    value_type localResult = 0;
    if ( this->comm().localSize() == 1 )
    {
        //localResult = ublas::inner_prod( M_vec, v.vec() );
        localResult = std::visit( [this]( auto && _v ) -> value_type { return ublas::inner_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        //localResult = ublas::inner_prod( ublas::project( M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ),
                                         //v.vec() );
        localResult = std::visit( [this, nLocalDofWithoutGhost]( auto && _v ) -> value_type {
                return ublas::inner_prod( ublas::project( this->M_vec, ublas::range( 0, nLocalDofWithoutGhost ) ), *_v );
                }, v.vec() );

    }

    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif
    return globalResult;
}   

#if FEELPP_HAS_PETSC
template< typename T, typename Storage > 
void 
VectorUblasContiguousGhosts<T, Storage>::setVector( const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->operator=( v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) = v( i );
    }
}  

template< typename T, typename Storage > 
void 
VectorUblasContiguousGhosts<T, Storage>::addVector( const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( 1., v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) += v( i );
    }
}  
template< typename T, typename Storage > 
void 
VectorUblasContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( a, v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) += a * v( i );
    }
}  
template< typename T, typename Storage > 
void 
VectorUblasContiguousGhosts<T, Storage>::subVector( const VectorPetsc<T> & v )
{
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( -1., v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) -= v( i );
    }
}
template< typename T, typename Storage > 
void 
VectorUblasContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( -a, v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) -= a * v( i );
    }
}  
template< typename T, typename Storage > 
typename VectorUblasContiguousGhosts<T, Storage>::value_type 
VectorUblasContiguousGhosts<T, Storage>::dotVector( const VectorPetsc<T> & v ) const
{
    if constexpr ( !is_vector_slice )
    {
        return this->vectorPetsc()->dot( v );
    }
    else
    {
        value_type localRes = 0;
        for( size_type i = 0; i < this->map().nLocalDofWithoutGhost(); ++i )
            localRes += this->operator()( i ) * v( i );
        value_type globalRes = localRes;
#ifdef FEELPP_HAS_MPI
        if( this->comm().size() > 1 )
            mpi::all_reduce( this->comm(), localRes, globalRes, std::plus<value_type>() );
#endif
        return globalRes;
    }
}
#endif // FEELPP_HAS_PETSC

template< typename T, typename Storage > 
VectorUblasRange<T, Storage> * VectorUblasContiguousGhosts<T, Storage>::rangeImpl( const range_type & rangeActive, const range_type & rangeGhost )
{
    if constexpr ( is_vector_range ) // should not happen, but still implemented
    {
        size_type rActiveStart = this->startActive() + rangeActive.start();
        size_type rActiveSize = std::min( this->map().nLocalDofWithoutGhost(), rangeActive.size() );
        range_type rActive = range_type( rActiveStart, rActiveStart + rActiveSize );

        size_type rGhostStart = this->startGhost() + rangeGhost.start();
        size_type rGhostSize = std::min( M_vec.size() - this->map().nLocalDofWithoutGhost(), rangeGhost.size() );
        range_type rGhost = range_type( rGhostStart, rGhostStart + rGhostSize );

        return new VectorUblasRange<T, typename Storage::vector_type>( M_vec.data().expression(), rActive, rGhost, this->mapPtr() );
    }
    else if constexpr ( is_vector_slice ) // should not happen, but still implemented
    {
        slice_type sActive = slice_type( this->startActive() + rangeActive.start() * M_vec.stride(), M_vec.stride(), std::min( this->map().nLocalDofWithoutGhost(), rangeActive.size() ) );
        slice_type sGhost = slice_type( this->startGhost() + rangeGhost.start() * M_vec.stride(), M_vec.stride(), std::min( M_vec.size() - this->map().nLocalDofWithoutGhost(), rangeGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, sGhost, this->mapPtr() );
    }
    else
    {
        size_type rGhostStart = this->startGhost() + rangeGhost.start();
        size_type rGhostSize = rangeGhost.size();
        range_type rGhost = range_type( rGhostStart, rGhostStart + rGhostSize );
        return new VectorUblasRange<T, Storage>( M_vec, rangeActive, rGhost, this->mapPtr() );
    }
}

template< typename T, typename Storage >
VectorUblasSlice<T, Storage> * VectorUblasContiguousGhosts<T, Storage>::sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost )
{
    if constexpr ( is_vector_range ) // should not happen, but still implemented
    {
        slice_type sActive = slice_type( this->startActive() + sliceActive.start(), sliceActive.stride(), std::min( this->map().nLocalDofWithoutGhost(), sliceActive.size() ) );
        slice_type sGhost = slice_type( this->startGhost() + sliceGhost.start(), sliceGhost.stride(), std::min( M_vec.size() - this->map().nLocalDofWithoutGhost(), sliceGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, sGhost, this->mapPtr() );
    }
    else if constexpr ( is_vector_slice ) // should not happen, but still implemented
    {
        slice_type sActive = slice_type( this->startActive() + sliceActive.start() * M_vec.stride(), M_vec.stride() * sliceActive.stride(), std::min( this->map().nLocalDofWithoutGhost(), sliceActive.size() ) );
        slice_type sGhost = slice_type( this->startGhost() + sliceGhost.start() * M_vec.stride(), M_vec.stride() * sliceGhost.stride(), std::min( M_vec.size() - this->map().nLocalDofWithoutGhost(), sliceGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, sGhost, this->mapPtr() );
    }
    else
    {
        slice_type sGhost = slice_type( this->startGhost() + sliceGhost.start(), sliceGhost.stride(), sliceGhost.size() );
        return new VectorUblasSlice<T, Storage>( M_vec, sliceActive, sGhost, this->mapPtr() );
    }
}

#if FEELPP_HAS_PETSC
template< typename T, typename Storage > 
VectorPetsc<T> * VectorUblasContiguousGhosts<T, Storage>::vectorPetscImpl() const
{
    const T * beginV = ( this->map().nLocalDofWithGhost() > 0 ) ? std::addressof( * M_vec.begin() ) : nullptr;
    if ( this->comm().size() > 1 )
        return new VectorPetscMPI<T>( beginV, this->mapPtr() );
    else
        return new VectorPetsc<T>( beginV, this->mapPtr() );
}
#endif

/*-----------------------------------------------------------------------------*/

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::init( const size_type n, const size_type n_local, const bool fast ) 
{
    FEELPP_ASSERT ( n_local <= n )
    ( n_local )( n )
    ( this->comm().rank() )
    ( this->comm().size() ).error( "Invalid local vector size" );

    // Clear the data structures if already initialized
    if ( this->isInitialized() )
        this->clear();

    super_type::init( n, n_local, fast );

    // Initialize data structures
    this->resize( this->localSize() );

    // Set the initialized flag
    this->M_is_initialized = true;

    DVLOG(2) << "        global size = " << n << "\n";
    DVLOG(2) << "        global size = " << n_local << "\n";
    DVLOG(2) << "        global size = " << this->size() << "\n";
    DVLOG(2) << "        local  size = " << this->localSize() << "\n";
    DVLOG(2) << "  first local index = " << this->firstLocalIndex() << "\n";
    DVLOG(2) << "   last local index = " << this->lastLocalIndex() << "\n";

    // Zero the components unless directed otherwise
    if ( !fast )
        this->zero();
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::self_type *
VectorUblasNonContiguousGhosts<T, Storage>::clonePtr() const
{
    return new VectorUblasNonContiguousGhosts<T, Storage>( *this );
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::base_type *
VectorUblasNonContiguousGhosts<T, Storage>::emptyPtr() const
{
    if constexpr ( is_vector_proxy )
        return new VectorUblasNonContiguousGhosts<T, typename Storage::vector_type>( );
    else
        return new VectorUblasNonContiguousGhosts<T, Storage>( );
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::resize( size_type n )
{
    if constexpr ( std::is_same_v<storage_type, vector_storage_type> )
        M_vec.resize( n );
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::clear( )
{
    if constexpr ( std::is_same_v<storage_type, vector_storage_type> )
    {
        M_vec.resize( 0 );
        M_vecNonContiguousGhosts.resize( 0 );
    }
}
        
template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type VectorUblasNonContiguousGhosts<T, Storage>::operator()( size_type i ) const
{
    const size_type nLocalActiveDof = this->map().nLocalDofWithoutGhost();
    if ( i < nLocalActiveDof )
        return M_vec.operator()( i );
    else
        return M_vecNonContiguousGhosts.operator()( i-nLocalActiveDof );
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type& VectorUblasNonContiguousGhosts<T, Storage>::operator()( size_type i )
{
    const size_type nLocalActiveDof = this->map().nLocalDofWithoutGhost();
    if ( i < nLocalActiveDof )
        return M_vec.operator()( i );
    else
        return M_vecNonContiguousGhosts.operator()( i-nLocalActiveDof );
}

template< typename T, typename Storage > 
typename VectorUblasNonContiguousGhosts<T, Storage>::size_type
VectorUblasNonContiguousGhosts<T, Storage>::startActive() const
{
    if constexpr ( is_vector_proxy )
        return M_vec.start();
    else
        return 0;
}

template< typename T, typename Storage > 
typename VectorUblasNonContiguousGhosts<T, Storage>::size_type
VectorUblasNonContiguousGhosts<T, Storage>::startGhost() const
{
    if constexpr ( is_vector_proxy )
        return M_vecNonContiguousGhosts.start();
    else
        return 0;
}
        
template< typename T, typename Storage > 
void VectorUblasNonContiguousGhosts<T, Storage>::setConstant( value_type v )
{
    M_vec = ublas::scalar_vector<value_type>( M_vec.size(), v );
    M_vecNonContiguousGhosts = ublas::scalar_vector<value_type>( M_vecNonContiguousGhosts.size(), v );
}
  

template< typename T, typename Storage > 
void VectorUblasNonContiguousGhosts<T, Storage>::setZero()
{
    M_vec = ublas::zero_vector<value_type>( M_vec.size() );
    M_vecNonContiguousGhosts = ublas::zero_vector<value_type>( M_vecNonContiguousGhosts.size() );
}
  

template< typename T, typename Storage > 
void VectorUblasNonContiguousGhosts<T, Storage>::scale( const value_type factor )
{
    M_vec.operator*=( factor );
    M_vecNonContiguousGhosts.operator*=( factor );
}
 
template< typename T, typename Storage > 
typename VectorUblasNonContiguousGhosts<T, Storage>::real_type VectorUblasNonContiguousGhosts<T, Storage>::min( bool parallel ) const
{
    checkInvariants();

    size_type nActiveDof = this->map().nLocalDofWithoutGhost();
    real_type local_min = ( nActiveDof > 0 ) ?
        *std::min_element( M_vec.begin(), M_vec.end() ) :
        std::numeric_limits<real_type>::max();
    real_type global_min = local_min;

#ifdef FEELPP_HAS_MPI
    if ( parallel && this->comm().size() > 1 )
    {
        MPI_Allreduce( &local_min, &global_min, 1,
                MPI_DOUBLE, MPI_MIN, this->comm() );
    }
#endif
    return global_min;
}

template< typename T, typename Storage > 
typename VectorUblasNonContiguousGhosts<T, Storage>::real_type VectorUblasNonContiguousGhosts<T, Storage>::max( bool parallel ) const
{
    checkInvariants();

    size_type nActiveDof = this->map().nLocalDofWithoutGhost();
    real_type local_max = ( nActiveDof > 0 ) ?
        *std::max_element( M_vec.begin(), M_vec.end() ) :
        std::numeric_limits<real_type>::min();
    real_type global_max = local_max;

#ifdef FEELPP_HAS_MPI
    if ( parallel && this->comm().size() > 1 )
    {
        MPI_Allreduce( &local_max, &global_max, 1,
                MPI_DOUBLE, MPI_MAX, this->comm() );
    }
#endif
    return global_max;
}
        
template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::real_type VectorUblasNonContiguousGhosts<T, Storage>::l1Norm() const
{
    checkInvariants();

    real_type local_l1 = ublas::norm_1( M_vec );
    real_type global_l1 = local_l1;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
    {
        mpi::all_reduce( this->comm(), local_l1, global_l1, std::plus<real_type>() );
    }
#endif
    return global_l1;
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::real_type VectorUblasNonContiguousGhosts<T, Storage>::l2Norm() const
{
    checkInvariants();

    real_type local_norm2 = ublas::inner_prod( M_vec, M_vec );
    real_type global_norm2 = local_norm2;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), local_norm2, global_norm2, std::plus<real_type>() );
#endif

    return math::sqrt( global_norm2 );
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::real_type VectorUblasNonContiguousGhosts<T, Storage>::linftyNorm() const
{
    checkInvariants();

    real_type local_norminf = ublas::norm_inf( M_vec );
    real_type global_norminf = local_norminf;

#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
    {
        mpi::all_reduce( this->comm(), local_norminf, global_norminf, mpi::maximum<real_type>() );
    }
#endif
    return global_norminf;
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type VectorUblasNonContiguousGhosts<T, Storage>::sum() const
{
    checkInvariants();

    value_type local_sum = ublas::sum( M_vec );
    value_type global_sum = local_sum;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), local_sum, global_sum, std::plus<value_type>() );
#endif

    return global_sum;
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::checkInvariants() const
{
    DCHECK( this->isInitialized() ) <<  "vector not initialized" ;
    DCHECK( this->localSize() <= this->size() ) << "vector invalid size: " << this->size() << "," << this->localSize();
    DCHECK( this->localSize() == ( M_vec.size() + M_vecNonContiguousGhosts.size() ) ) << "vector invalid size: " 
        << M_vec.size() << ","
        << M_vecNonContiguousGhosts.size() << ","
        << this->localSize();
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::setVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec = v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec = *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() = ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther]( auto && _v) { this->M_vec = ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() = ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { this->M_vecNonContiguousGhosts = ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::setVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec = v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec = *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() = v.vec();
            std::visit( [this]( auto && _v ) { this->M_vec = *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() = v.vecNonContiguousGhosts();
            std::visit( [this]( auto && _v ) { this->M_vecNonContiguousGhosts = *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::addVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec += *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() += ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther]( auto && _v) { this->M_vec += ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() += ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { this->M_vecNonContiguousGhosts += ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::addVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec += *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() += v.vec();
            std::visit( [this]( auto && _v ) { this->M_vec += *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() += v.vecNonContiguousGhosts();
            std::visit( [this]( auto && _v ) { this->M_vecNonContiguousGhosts += *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += a*v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec += a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() += a * ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            std::visit( [this, a, nLocalDofWithoutGhostOther]( auto && _v) { this->M_vec += a * ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() += a * ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, a, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { this->M_vecNonContiguousGhosts += a * ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec += a * v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec += a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() += a * v.vec();
            std::visit( [this, a]( auto && _v ) { this->M_vec += a * (*_v); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() += a * v.vecNonContiguousGhosts();
            std::visit( [this, a]( auto && _v ) { this->M_vecNonContiguousGhosts += a * (*_v); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::subVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec -= *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() -= ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther]( auto && _v) { this->M_vec -= ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() -= ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { this->M_vecNonContiguousGhosts -= ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::subVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= v.vec();
        std::visit( [this]( auto && _v ) { this->M_vec -= *_v; }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() -= v.vec();
            std::visit( [this]( auto && _v ) { this->M_vec -= *_v; }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() -= v.vecNonContiguousGhosts();
            std::visit( [this]( auto && _v ) { this->M_vecNonContiguousGhosts -= *_v; }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= a*v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec -= a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() -= a * ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            std::visit( [this, a, nLocalDofWithoutGhostOther]( auto && _v) { this->M_vec -= a * ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() -= a * ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, a, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { this->M_vecNonContiguousGhosts -= a * ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec -= a * v.vec();
        std::visit( [this, a]( auto && _v ) { this->M_vec -= a * (*_v); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() -= a * v.vec();
            std::visit( [this, a]( auto && _v ) { this->M_vec -= a * (*_v); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() -= a * v.vecNonContiguousGhosts();
            std::visit( [this, a]( auto && _v ) { this->M_vecNonContiguousGhosts -= a * (*_v); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::mulVector( const VectorUblasContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec *= v.vec();
        // c.f. [noalias]
        std::visit( [this]( auto && _v ) { ublas::noalias( this->M_vec ) = ublas::element_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        size_type nLocalGhostsOther = v.map().nLocalGhosts();
        if ( nLocalDofWithoutGhostOther > 0 )
        {
            //this->vec() *= ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostOther ) );
            // c.f. [noalias]
            std::visit( [this, nLocalDofWithoutGhostOther]( auto && _v) { ublas::noalias( this->M_vec ) = ublas::element_prod( this->M_vec, ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ) ); }, v.vec() );
        }
        if ( nLocalGhostsOther > 0 )
        {
            //this->vecNonContiguousGhosts() *= ublas::project( v.vec(), ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) );
            std::visit( [this, nLocalDofWithoutGhostOther, nLocalGhostsOther]( auto && _v ) { ublas::noalias( this->M_vecNonContiguousGhosts ) = ublas::element_prod( this->M_vecNonContiguousGhosts, ublas::project( *_v, ublas::range( nLocalDofWithoutGhostOther, nLocalDofWithoutGhostOther+nLocalGhostsOther ) ) ); }, v.vec() );
        }
    }
}

template< typename T, typename Storage >
void VectorUblasNonContiguousGhosts<T, Storage>::mulVector( const VectorUblasNonContiguousGhostsBase<T> & v )
{
    if ( this->comm().localSize() == 1 )
    {
        //M_vec *= v.vec();
        // c.f. [noalias]
        std::visit( [this]( auto && _v ) { ublas::noalias( this->M_vec ) = ublas::element_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhost = this->map().nLocalDofWithoutGhost();
        size_type nLocalGhosts = this->map().nLocalGhosts();
        if ( nLocalDofWithoutGhost > 0 )
        {
            //this->vec() *= v.vec();
            // c.f. [noalias]
            std::visit( [this]( auto && _v ) { ublas::noalias( this->M_vec ) = ublas::element_prod( this->M_vec, *_v ); }, v.vec() );
        }
        if ( nLocalGhosts > 0 )
        {
            //this->vecNonContiguousGhosts() *= v.vecNonContiguousGhosts();
            // c.f. [noalias]
            std::visit( [this]( auto && _v ) { ublas::noalias( this->M_vecNonContiguousGhosts ) = ublas::element_prod( this->M_vecNonContiguousGhosts, *_v ); }, v.vecNonContiguousGhosts() );
        }
    }
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type 
VectorUblasNonContiguousGhosts<T, Storage>::dotVector( const VectorUblasContiguousGhostsBase<T> & v ) const
{
    value_type localResult = 0;
    if ( this->comm().localSize() == 1 )
    {
        //localResult = ublas::inner_prod( M_vec, v.vec() );
        localResult = std::visit( [this]( auto && _v ) -> value_type { return ublas::inner_prod( this->M_vec, *_v ); }, v.vec() );
    }
    else
    {
        size_type nLocalDofWithoutGhostOther = v.map().nLocalDofWithoutGhost();
        //localResult = ublas::inner_prod( this->vec(),
                                         //ublas::project( v.vec(), ublas::range( 0, nLocalDofWithoutGhostV ) ) );
        localResult = std::visit( [this, nLocalDofWithoutGhostOther]( auto && _v ) -> value_type { 
                return ublas::inner_prod( this->M_vec,
                                          ublas::project( *_v, ublas::range( 0, nLocalDofWithoutGhostOther ) ) );
                }, v.vec() );
    }

    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif
    return globalResult;
}

template< typename T, typename Storage >
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type 
VectorUblasNonContiguousGhosts<T, Storage>::dotVector( const VectorUblasNonContiguousGhostsBase<T> & v ) const
{
    //value_type localResult = ublas::inner_prod( M_vec, v.vec() );
    value_type localResult = std::visit( [this]( auto && _v ) -> value_type { return ublas::inner_prod( this->M_vec, *_v ); }, v.vec() );

    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif
    return globalResult;
}   

#if FEELPP_HAS_PETSC
template< typename T, typename Storage > 
void 
VectorUblasNonContiguousGhosts<T, Storage>::setVector( const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->operator=( v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) = v( i );
    }
}  

template< typename T, typename Storage > 
void 
VectorUblasNonContiguousGhosts<T, Storage>::addVector( const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( 1., v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) += v( i );
    }
}  
template< typename T, typename Storage > 
void 
VectorUblasNonContiguousGhosts<T, Storage>::maddVector( const value_type & a, const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( a, v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) += a * v( i );
    }
}  
template< typename T, typename Storage > 
void 
VectorUblasNonContiguousGhosts<T, Storage>::subVector( const VectorPetsc<T> & v )
{
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( -1., v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) -= v( i );
    }
}
template< typename T, typename Storage > 
void 
VectorUblasNonContiguousGhosts<T, Storage>::msubVector( const value_type & a, const VectorPetsc<T> & v )
{ 
    if constexpr ( !is_vector_slice )
    {
        this->vectorPetsc()->add( -a, v );
    }
    else
    {
        for( size_type i = 0; i < this->localSize(); ++i )
            this->operator()( i ) -= a * v( i );
    }
}  
template< typename T, typename Storage > 
typename VectorUblasNonContiguousGhosts<T, Storage>::value_type 
VectorUblasNonContiguousGhosts<T, Storage>::dotVector( const VectorPetsc<T> & v ) const
{
    if constexpr ( !is_vector_slice )
    {
        return this->vectorPetsc()->dot( v );
    }
    else
    {
        value_type localRes = 0;
        for( size_type i = 0; i < this->map().nLocalDofWithoutGhost(); ++i )
            localRes += this->operator()( i ) * v( i );
        value_type globalRes = localRes;
#ifdef FEELPP_HAS_MPI
        if( this->comm().size() > 1 )
            mpi::all_reduce( this->comm(), localRes, globalRes, std::plus<value_type>() );
#endif
        return globalRes;
    }
}
#endif // FEELPP_HAS_PETSC

template< typename T, typename Storage > 
VectorUblasBase<T> * VectorUblasNonContiguousGhosts<T, Storage>::rangeImpl( const range_type & rangeActive, const range_type & rangeGhost )
{
    if constexpr ( is_vector_range )
    {
        size_type rActiveStart = M_vec.start() + rangeActive.start();
        size_type rActiveSize = std::min( M_vec.size(), rangeActive.size() );
        range_type rActive = range_type( rActiveStart, rActiveStart + rActiveSize );

        size_type rGhostStart = M_vecNonContiguousGhosts.start() + rangeGhost.start();
        size_type rGhostSize = std::min( M_vecNonContiguousGhosts.size(), rangeGhost.size() );
        range_type rGhost = range_type( rGhostStart, rGhostStart + rGhostSize );

        return new VectorUblasRange<T, typename Storage::vector_type>( M_vec.data().expression(), rActive, M_vecNonContiguousGhosts.data().expression(), rGhost, this->mapPtr() );
    }
    else if constexpr ( is_vector_slice )
    {
        slice_type sActive = slice_type( M_vec.start() + rangeActive.start() * M_vec.stride(), M_vec.stride(), std::min( M_vec.size(), rangeActive.size() ) );
        slice_type sGhost = slice_type( M_vecNonContiguousGhosts.start() + rangeGhost.start() * M_vecNonContiguousGhosts.stride(), M_vecNonContiguousGhosts.stride(), std::min( M_vecNonContiguousGhosts.size(), rangeGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, M_vecNonContiguousGhosts.data().expression(), sGhost, this->mapPtr() );
    }
    else
    {
        return new VectorUblasRange<T, Storage>( M_vec, rangeActive, M_vecNonContiguousGhosts, rangeGhost, this->mapPtr() );
    }
}

template< typename T, typename Storage > 
VectorUblasBase<T> * VectorUblasNonContiguousGhosts<T, Storage>::sliceImpl( const slice_type & sliceActive, const slice_type & sliceGhost )
{
    if constexpr ( is_vector_range )
    {
        slice_type sActive = slice_type( M_vec.start() + sliceActive.start(), sliceActive.stride(), std::min( M_vec.size(), sliceActive.size() ) );
        slice_type sGhost = slice_type( M_vecNonContiguousGhosts.start() + sliceGhost.start(), sliceGhost.stride(), std::min( M_vecNonContiguousGhosts.size(), sliceGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, M_vecNonContiguousGhosts.data().expression(), sGhost, this->mapPtr() );
    }
    else if constexpr ( is_vector_slice )
    {
        slice_type sActive = slice_type( M_vec.start() + sliceActive.start() * M_vec.stride(), M_vec.stride() * sliceActive.stride(), std::min( M_vec.size(), sliceActive.size() ) );
        slice_type sGhost = slice_type( M_vecNonContiguousGhosts.start() + sliceGhost.start() * M_vecNonContiguousGhosts.stride(), M_vecNonContiguousGhosts.stride() * sliceGhost.stride(), std::min( M_vecNonContiguousGhosts.size(), sliceGhost.size() ) );
        return new VectorUblasSlice<T, typename Storage::vector_type>( M_vec.data().expression(), sActive, M_vecNonContiguousGhosts.data().expression(), sGhost, this->mapPtr() );
    }
    else
    {
        return new VectorUblasSlice<T, Storage>( M_vec, sliceActive, M_vecNonContiguousGhosts, sliceGhost, this->mapPtr() );
    }
}

#if FEELPP_HAS_PETSC
template< typename T, typename Storage > 
VectorPetsc<T> * VectorUblasNonContiguousGhosts<T, Storage>::vectorPetscImpl() const
{
    if ( this->comm().size() > 1 )
    {
        const T * beginA = ( this->map().nLocalDofWithoutGhost() > 0 ) ? std::addressof( * M_vec.begin() ) : nullptr;
        const T * beginG = ( this->map().nLocalGhosts() > 0 ) ? std::addressof( * M_vecNonContiguousGhosts.begin() ) : nullptr;
        return new VectorPetscMPIRange<T>( beginA, beginG, this->mapPtr() );
    }
    else
        return new VectorPetsc<T>( std::addressof( * M_vec.begin() ), this->mapPtr() );
}
#endif

//template< typename T, typename Storage >
//VectorUblasRange<T, Storage>::VectorUblasRange( VectorUblasContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost ):
    //super_type( ublas::vector_range( v.M_vec, rangeActive ), ublas::vector_range( v.M_vec, rangeGhost ), v.mapPtr() )
//{
//}

//template< typename T, typename Storage >
//VectorUblasRange<T, Storage>::VectorUblasRange( VectorUblasNonContiguousGhosts<T, Storage> & v, const range_type & rangeActive, const range_type & rangeGhost ):
    //super_type( ublas::vector_range( v.M_vec, rangeActive ), ublas::vector_range( v.M_vecNonContiguousGhosts, rangeGhost ), v.mapPtr() )
//{
//}

template< typename T, typename Storage >
VectorUblasRange<T, Storage>::VectorUblasRange( Storage & v, const range_type & range, const datamap_ptrtype & dm ):
    super_type( ublas::vector_range( v, range ), ublas::vector_range( v, range_type(0,0) ), dm )
{
}

template< typename T, typename Storage >
VectorUblasRange<T, Storage>::VectorUblasRange( Storage & v, const range_type & rangeActive, const range_type & rangeGhost, const datamap_ptrtype & dm ):
    super_type( ublas::vector_range( v, rangeActive ), ublas::vector_range( v, rangeGhost ), dm )
{
}

template< typename T, typename Storage >
VectorUblasRange<T, Storage>::VectorUblasRange( Storage & vActive, const range_type & rangeActive, Storage & vGhost, const range_type & rangeGhost, const datamap_ptrtype & dm ):
    super_type( ublas::vector_range( vActive, rangeActive ), ublas::vector_range( vGhost, rangeGhost ), dm )
{
}

//template< typename T, typename Storage >
//VectorUblasSlice<T, Storage>::VectorUblasSlice( VectorUblasContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost ):
    //super_type( ublas::vector_slice( v.M_vec, sliceActive ), ublas::vector_slice( v.M_vec, sliceGhost ), v.mapPtr() )
//{
//}

//template< typename T, typename Storage >
//VectorUblasSlice<T, Storage>::VectorUblasSlice( VectorUblasNonContiguousGhosts<T, Storage> & v, const slice_type & sliceActive, const slice_type & sliceGhost ):
    //super_type( ublas::vector_slice( v.M_vec, sliceActive ), ublas::vector_slice( v.M_vecNonContiguousGhosts, sliceGhost ), v.mapPtr() )
//{
//}

template< typename T, typename Storage >
VectorUblasSlice<T, Storage>::VectorUblasSlice( Storage & v, const slice_type & slice, const datamap_ptrtype & dm ):
    super_type( ublas::vector_slice( v, slice ), ublas::vector_slice( v, slice_type(0,1,0) ), dm )
{
}

template< typename T, typename Storage >
VectorUblasSlice<T, Storage>::VectorUblasSlice( Storage & v, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm ):
    super_type( ublas::vector_slice( v, sliceActive ), ublas::vector_slice( v, sliceGhost ), dm )
{
}

template< typename T, typename Storage >
VectorUblasSlice<T, Storage>::VectorUblasSlice( Storage & vActive, slice_type const& sliceActive,
                Storage & vGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm ):
    super_type( ublas::vector_slice( vActive, sliceActive ), ublas::vector_slice( vGhost, sliceGhost ), dm )
{
}

} // namespace detail

/*-----------------------------------------------------------------------------*/
template< typename T >
VectorUblas<T>::VectorUblas():
    super_type(),
    M_vectorImpl( new detail::VectorUblasContiguousGhosts<T>( this->mapPtr() ) )
{
    //CHECK(false) << "empty ctor: TODO";
}
template< typename T >
VectorUblas<T>::VectorUblas( size_type s ):
    super_type( s, Environment::worldCommSeqPtr() ),
    M_vectorImpl( new detail::VectorUblasContiguousGhosts<T>( this->mapPtr() ) )
{
    this->init( s, s, false );
}
template< typename T >
VectorUblas<T>::VectorUblas( datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( new detail::VectorUblasContiguousGhosts<T>( dm ) )
{
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false );
}
template< typename T >
VectorUblas<T>::VectorUblas( size_type s, size_type n_local ):
    super_type( s, n_local, Environment::worldCommSeqPtr() ),
    M_vectorImpl( new detail::VectorUblasContiguousGhosts<T>( this->mapPtr() ) )
{
    this->init( this->size(), this->localSize(), false );
}
template< typename T >
VectorUblas<T>::VectorUblas( VectorUblas<T>& v, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( v.M_vectorImpl->range( rangeActive, rangeGhost ) )
{
    M_vectorImpl->setMap( dm );
}
template< typename T >
VectorUblas<T>::VectorUblas( VectorUblas<T>& v, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( v.M_vectorImpl->slice( sliceActive, sliceGhost ) )
{
    M_vectorImpl->setMap( dm );
}
template< typename T >
VectorUblas<T>::VectorUblas( ublas::vector<T>& v, range_type const& range, datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( new detail::VectorUblasRange<T, ublas::vector<T>>( v, range, dm ) )
{
}
template< typename T >
VectorUblas<T>::VectorUblas( ublas::vector<T>& vActive, slice_type const& sliceActive,
        ublas::vector<value_type>& vGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( new detail::VectorUblasSlice<T, ublas::vector<T>>( vActive, sliceActive, vGhost, sliceGhost, dm ) )
{
}

template< typename T >
VectorUblas<T>::VectorUblas( size_type nActiveDof, value_type * arrayActiveDof,
        size_type nGhostDof, value_type * arrayGhostDof,
        datamap_ptrtype const& dm ):
    super_type( dm ),
    M_vectorImpl( new detail::VectorUblasNonContiguousGhosts<T, ublas::vector<T, Feel::detail::shallow_array_adaptor<T>>>( 
                ublas::vector<T, Feel::detail::shallow_array_adaptor<T>>( Feel::detail::shallow_array_adaptor<T>( nActiveDof, arrayActiveDof ) ), 
                ublas::vector<T, Feel::detail::shallow_array_adaptor<T>>( Feel::detail::shallow_array_adaptor<T>( nGhostDof, arrayGhostDof ) ), 
                dm ) )
{
}

template< typename T >
Vector<T> & VectorUblas<T>::operator=( const Vector<T> & v )
{
    if ( this == &v )
        return *this;

    if ( !this->map().isCompatible( v.map() ) )
    {
        this->setMap( v.mapPtr() );
        this->resize( this->map().nLocalDofWithGhost() );
    }

#if !defined(NDEBUG)
    //checkInvariants();

    FEELPP_ASSERT( this->localSize() == v.localSize() &&
                   this->map().nLocalDofWithoutGhost() == v.map().nLocalDofWithoutGhost() )
    ( this->localSize() )( this->map().nLocalDofWithoutGhost() )
    ( this->vectorImpl() )
    ( v.localSize() )( v.map().nLocalDofWithoutGhost() ).warn( "may be vector invalid copy" );
#endif

    //return M_vectorImpl->operator=( v );
    // do not call M_vectorImpl->operator= since map setting and resizing has already been performed
    // however, optimize for VectorUblas which can use setVector directly instead of for loop
    const VectorUblas<T> * vecUblas = dynamic_cast<const VectorUblas<T> *>( &v );
    if( vecUblas )
        M_vectorImpl->setVector( *( vecUblas->M_vectorImpl ) );
    else
        M_vectorImpl->set( v );
    return *this;
}

template< typename T >
VectorUblas<T> & VectorUblas<T>::operator=( const VectorUblas<T> & v )
{
    if ( this == &v )
        return *this;

    if ( !this->map().isCompatible( v.map() ) )
    {
        this->setMap( v.mapPtr() );
        this->resize( this->map().nLocalDofWithGhost() );
    }

#if !defined(NDEBUG)
    //checkInvariants();

    FEELPP_ASSERT( this->localSize() == v.localSize() &&
                   this->map().nLocalDofWithoutGhost() == v.map().nLocalDofWithoutGhost() )
    ( this->localSize() )( this->map().nLocalDofWithoutGhost() )
    ( this->vectorImpl() )
    ( v.localSize() )( v.map().nLocalDofWithoutGhost() ).warn( "may be vector invalid copy" );
#endif

    //return M_vectorImpl->operator=( v );
    // do not call M_vectorImpl->operator= since map setting and resizing has already been performed
    M_vectorImpl->setVector( *( v.M_vectorImpl ) );
    return *this;
}

/*-----------------------------------------------------------------------------*/
// Explicit instantiations
template class VectorUblas<double>;
namespace detail {
template class VectorUblasBase< double >;
template class VectorUblasContiguousGhosts< double, ublas::vector<double> >;
template class VectorUblasNonContiguousGhosts< double, ublas::vector<double> >;
template class VectorUblasRange< double, ublas::vector<double> >;
template class VectorUblasSlice< double, ublas::vector<double> >;
}


} // namespace Feel
