/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-06-20

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file vectoreigen.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-06-20
 */
/**
 * \sa vectoreigen.hpp
 * FEELPP_INSTANTIATE_VECTOREIGEN is never defined except in vectoreigen.cpp
 * where we do the instantiate. This allows to reduce the VectorEigen
 * instantiation to the strict minimum
 */
#define FEELPP_INSTANTIATE_VECTOREIGEN 1

#include <Eigen/Core>


#include <feel/feelalg/vectoreigen.hpp>

#if defined( FEELPP_HAS_TBB )
#include <parallel_for.h>
#include <blocked_range.h>
#endif // FEELPP_HAS_TBB

namespace Feel
{
template <typename T>
VectorEigen<T>::VectorEigen()
    :
    super1(),
    M_vec()
{
}

template <typename T>
VectorEigen<T>::VectorEigen( size_type __s, WorldComm const& _worldComm )
    :
    super1( __s, _worldComm ),
    M_vec( __s )
{
    this->init( __s, __s, false );
}

template <typename T>
VectorEigen<T>::VectorEigen( datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( dm->nDof() )
{
    //this->init( dm.nGlobalElements(), dm.nMyElements(), false );
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false );
}

template <typename T>
VectorEigen<T>::VectorEigen( size_type __s, size_type __n_local, WorldComm const& _worldComm  )
    :
    super1( __s, __n_local, _worldComm ),
    M_vec( __s )
{
    this->init( this->size(), this->localSize(), false );
}

template <typename T>
VectorEigen<T>::VectorEigen( VectorEigen const & m )
    :
    super1( m ),
    M_vec( m.M_vec )
{
    DVLOG(2) << "[VectorEigen] copy constructor with range: size:" << this->size() << ", start:" << this->start() << "\n";
    DVLOG(2) << "[VectorEigen] copy constructor with range: size:" << this->vec().size() << "\n";
}

template <typename T>
VectorEigen<T>::~VectorEigen()
{
    this->clear();
}

template <typename T>
typename VectorEigen<T>::clone_ptrtype
VectorEigen<T>::clone () const
{
    return clone_ptrtype( new VectorEigen<T>( *this ) );
}

template <typename T>
void
VectorEigen<T>::resize( size_type s, bool preserve )
{
    M_vec.conservativeResize( s );
}

template <typename T>
size_type
VectorEigen<T>::start( ) const
{
    return 0;
}

template <typename T>
Vector<T> &
VectorEigen<T>::operator= ( const Vector<value_type> &V )
{
    checkInvariant();
    FEELPP_ASSERT( this->localSize() == V.localSize() )( this->localSize() )( V.localSize() ).warn ( "invalid vector size" );
    FEELPP_ASSERT( this->firstLocalIndex() == V.firstLocalIndex() &&
                   this->lastLocalIndex() == V.lastLocalIndex() )
    ( this->firstLocalIndex() )( this->lastLocalIndex() )
    ( this->vec().size() )
    ( V.firstLocalIndex() )( V.lastLocalIndex() ).warn( "may be vector invalid  copy" );

    for ( size_type i = 0; i < this->localSize(); ++i )
    {
        M_vec.operator()( i ) = V( V.firstLocalIndex() + i );
        //M_vec.operator()( i ) = V(  i );

    }

    return *this;
}

template <typename T>
void
VectorEigen<T>::init ( const size_type n,
                               const size_type n_local,
                               const bool      fast )
{
    FEELPP_ASSERT ( n_local <= n )
    ( n_local )( n )
    ( this->comm().rank() )
    ( this->comm().size() ).error( "Invalid local vector size" );

    // Clear the data structures if already initialized
    if ( this->isInitialized() )
        this->clear();

    super1::init( n, n_local, fast );

    // Initialize data structures
    M_vec.resize( this->localSize() );

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

template<typename T>
void
VectorEigen<T>::init( datamap_ptrtype const& dm )
{
    super1::init( dm );
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false );
}


template<typename T>
void
VectorEigen<T>::init ( const size_type n,
                               const bool      fast )
{
    this->init( n,n,fast );
}

template<typename T>
void
VectorEigen<T>::clear()
{
    M_vec.resize( 0 );
}

template<typename T>
void
VectorEigen<T>::close() const
{
}

template<typename T>
void
VectorEigen<T>::printMatlab( const std::string filename, bool renumber  ) const
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
            DVLOG(2) << "[VectorEigen::printMatlab] adding .m extension to given file name '"
                          << filename << "'\n";
            name = filename + ".m";
        }
    }


    vector_type v_local;
    this->localizeToOneProcessor ( v_local, 0 );

    if ( this->comm().rank() == 0 )
    {
        std::ofstream file_out( name.c_str() );

        FEELPP_ASSERT( file_out )( filename ).error( "[VectorEigen::printMatlab] ERROR: File cannot be opened for writing." );

				std::string varName = "var_" + filename.substr(0,filename.find("."));
        file_out << varName <<" = [ ";
        file_out.precision( 16 );
        file_out.setf( std::ios::scientific );

        for ( size_type i = 0; i < this->size(); ++i )
        {
            file_out << v_local[i] << separator << std::endl;
        }

        file_out << "];" << std::endl;
    }
}


template <typename T>
void
VectorEigen<T>::localize ( Vector<T>& v_local_in ) const

{
    checkInvariant();

    VectorEigen<T>* v_local = dynamic_cast<VectorEigen<T>*>( &v_local_in );
    FEELPP_ASSERT( v_local != 0 ).error ( "dynamic_cast failed: invalid vector object" );

#if 0
    v_local->firstLocalIndex() = 0;

    v_local->size() =
        v_local->localSize() = this->size();

    v_local->lastLocalIndex() = this->size()-1;

    v_local->M_is_initialized =
        v_local->M_is_closed = true;
#else
    DataMap dm( this->size(), this->size() );
    v_local->init( this->size(), this->size() );

#endif // 0

    // Call localize on the vector's values.  This will help
    // prevent code duplication
    localize ( v_local->M_vec );

#ifndef FEELPP_HAS_MPI

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );

#endif
}



template <typename T >
void VectorEigen<T>::localize ( Vector<T>& v_local_in,
                                        const std::vector<size_type>& ) const
{
    checkInvariant();

    // We don't support the send list.  Call the less efficient localize(v_local_in)
    localize ( v_local_in );
}



template <typename T>
void
VectorEigen<T>::localize ( const size_type first_local_idx,
                                   const size_type last_local_idx,
                                   const std::vector<size_type>& send_list )
{
    // Only good for serial vectors
    FEELPP_ASSERT ( this->size() == this->localSize() )( this->size() )( this->localSize() ).error( "invalid local/global size" );
    FEELPP_ASSERT ( last_local_idx > first_local_idx )( last_local_idx )( first_local_idx ).error( "invalid first/last local indices" );
    FEELPP_ASSERT ( send_list.size() <= this->size() )( send_list.size() )( this->size() ).error( "invalid send list size" );
    FEELPP_ASSERT ( last_local_idx < this->size() )( last_local_idx )( this->size() ).error( "invalid last local index" );
    Feel::detail::ignore_unused_variable_warning( send_list );

    const size_type size       = this->size();
    const size_type local_size = ( last_local_idx - first_local_idx + 1 );

    // Don't bother for serial cases
    if ( ( first_local_idx == 0 ) &&
            ( local_size == size ) )
        return;

#if 0

    // Build a parallel vector, initialize it with the local
    // parts of (*this)
    VectorEigen<T> parallel_vec;

    parallel_vec.init ( size, local_size );

    // Copy part of *this into the parallel_vec
    for ( size_type i=first_local_idx; i<=last_local_idx; i++ )
        parallel_vec.operator()( i ) = this->operator()( i );

    // localize like normal
    parallel_vec.localize ( *this, send_list );
#endif
}

// template <typename T>
// void
// VectorEigen<T>::localize ( eigen::vector_range<vector_type >& /*v_local*/ ) const
// {
// }

// template <typename T>
// void
// VectorEigen<T>::localize ( eigen::vector_slice<vector_type >& /*v_local*/ ) const
// {
// }

template <typename T>
void
VectorEigen<T>::localize ( vector_type& v_local ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    vector_type v_local_in( this->size() );

#ifdef FEELPP_HAS_MPI

    if ( this->comm().size() > 1 )
    {
        v_local.setConstant( this->size(), 0. );
        v_local_in.setConstant( this->size(), 0. );

        for ( size_type i=0; i< this->localSize(); i++ )
        {
            v_local_in[i+this->firstLocalIndex()] = M_vec.operator[]( i );
        }

        boost::mpi::all_reduce(this->comm(), &v_local_in[0], v_local.size(), &v_local[0], std::plus<value_type>());
        DVLOG(2) << "[VectorEigen::localize] Allreduce size = " << v_local.size() << "\n";

    }

    else
    {
        FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );
        v_local = M_vec;
    }

#else

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );

#endif
}



template <typename T>
void
VectorEigen<T>::localizeToOneProcessor ( vector_type& v_local,
        const size_type pid ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    v_local.setConstant( this->size(), 0. );

    vector_type v_tmp( this->size() );
    v_tmp.setConstant( this->size(), 0. );

    for ( size_type i=0; i< this->localSize(); i++ )
        v_tmp[i+this->firstLocalIndex()] = this->operator()( this->firstLocalIndex()+i );

#ifdef FEELPP_HAS_MPI

    if ( this->comm().size() > 1 )
    {
        boost::mpi::reduce(this->comm(), &v_tmp[0], v_local.size(), &v_local[0], std::plus<value_type>(), pid);
    }

    else
    {
        v_local = v_tmp;
    }

#else

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );
    FEELPP_ASSERT ( pid == 0  )( pid ).error( "invalid pid in non MPI mode" );

#endif
}
template <typename T>
void
VectorEigen<T>::localizeToOneProcessor ( std::vector<value_type>& v_local,
        const size_type proc_id ) const
{
    vector_type eigenvector;
    localizeToOneProcessor( eigenvector, proc_id );
    v_local.resize( eigenvector.size() );
    Eigen::Map<vector_type> v( v_local.data(), v_local.size() );
    v = eigenvector;
}

template <typename T>
void
VectorEigen<T>::checkInvariant() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
    FEELPP_ASSERT ( this->localSize() <= this->size() )
    ( this->size() )( this->localSize() ).error( "vector invalid size" );
    FEELPP_ASSERT ( M_vec.size() == this->localSize() )
    ( M_vec.size() )( this->localSize() ).error( "vector invalid size" );
    FEELPP_ASSERT ( ( this->lastLocalIndex() - this->firstLocalIndex() ) == this->localSize() )
    ( this->size() )
    ( this->lastLocalIndex() )
    ( this->firstLocalIndex() )
    ( this->localSize() ).error( "vector invalid size" );
}

namespace detail
{
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
} //detail
template <typename T>
typename VectorEigen<T>::this_type
VectorEigen<T>::sqrt() const
{
    this_type _tmp( this->mapPtr() );

#if defined( FEELPP_HAS_TBB )
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, this->localSize() ),
                       detail::Sqrt<this_type>( *this, _tmp ) );
#else

    for ( int i = 0; i < ( int )this->localSize(); ++i )
        _tmp[i] = math::sqrt( this->operator[]( i ) );

#endif // FEELPP_HAS_TBB

    return _tmp;
}

template <typename T>
typename VectorEigen<T>::this_type
VectorEigen<T>::pow( int n ) const
{
    this_type _out( this->mapPtr() );

    _out.M_vec = M_vec.array().pow(n);

    return _out;
}

template <typename T>
void
VectorEigen<T>::insert ( const ublas::vector<T>& V,
                         const std::vector<size_type>& dof_indices )
{
    FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
}

template <typename T>
void VectorEigen<T>::addVector ( const Vector<value_type>& V_in,
                                 const MatrixSparse<value_type>& A_in )
{
    //FEELPP_ASSERT( 0 ).error( "invalid call, not implemented yet" );
    const VectorEigen<T>* V = dynamic_cast<const VectorEigen<T>*>( &V_in );
    //const MatrixEigenDense<T>* A = dynamic_cast<const MatrixEigenDense<T>*>( &A_in );
    const MatrixEigenDense<T>* A = dynamic_cast<const MatrixEigenDense<T>*>( &A_in );

    CHECK ( A != 0 ) << "Invalid Eigen matrix\n";

    M_vec += A->mat()*V->vec();
}
//
// instantiation
//
template class VectorEigen<double>;
template class VectorEigen<std::complex<double>>;
} // Feel
