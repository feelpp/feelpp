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
void VectorUblasBase<T>::set( const size_type i, const value_type & value )
{
#ifndef(NDEBUG)
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
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->setVector( *vecPetsc );
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
#ifndef(NDEBUG)
    checkInvariants();
#endif
    this->operator()( i ) += value;
}

template< typename T >
void VectorUblasBase<T>::add( const value_type & value )
{
#ifndef(NDEBUG)
    checkInvariants();
#endif
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += value;
}

template< typename T >
void VectorUblasBase<T>::add( const Vector<T> & v )
{
#ifndef(NDEBUG)
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->addVector( *vecUblas );
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->addVector( *vecPetsc );
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::add( const value_type & a, const Vector<T> & v )
{
#ifndef(NDEBUG)
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->maddVector( a, *vecUblas );
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->maddVector( a, *vecPetsc );
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += a * v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::sub( const Vector<T> & v )
{
#ifndef(NDEBUG)
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->subVector( *vecUblas );
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->subVector( *vecPetsc );
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) -= v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::sub( const value_type & a, const Vector<T> & v )
{
#ifndef(NDEBUG)
    checkInvariants();
#endif
    // Ublas case
    const VectorUblasBase<T> * vecUblas = dynamic_cast<const VectorUblasBase<T> *>( &v );
    if( vecUblas )
        return this->msubVector( a, *vecUblas );
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->msubVector( a, *vecPetsc );
    // Default
    for( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) -= a * v( i );
    return;
}

template< typename T >
void VectorUblasBase<T>::addVector( int * i, int n, value_type * v, size_type K = 0, size_type K2 = invalid_v<size_type> )
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
    // Petsc case
    const VectorPetsc<T> * vecPetsc = dynamic_cast<const VectorPetsc<T> *>( &v );
    if( vecPetsc )
        return this->dotVector( *vecPetsc );
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


template<typename T>
void VectorUblasBase<T>::printMatlab( const std::string & filename, bool renumber ) const
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
template<typename T>
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

template<typename T,>
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

template<typename T>
void VectorUblasBase<T,Storage>::localizeToOneProcessor( ublas::vector<T> & v_local, const size_type pid ) const
{
    checkInvariant();

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
    FEELPP_ASSERT( pid == 0  )( pid ).error( "invalid pid in non MPI mode" );
#endif
}

template<typename T>
void VectorUblasBase<T,Storage>::localizeToOneProcessor( std::vector<T> & v_local, const size_type pid ) const
{
    ublas::vector<T> ublasV;
    this->localizeToOneProcessor( ublasV, pid );
    v_local.resize( ublasV.size() );
    std::copy( ublasV.begin(), ublasV.end(), v_local.begin() );
}


} // namespace detail

} // namespace Feel
