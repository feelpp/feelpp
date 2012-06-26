/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-06-20
 */
/**
 * \sa vectoreigen.hpp
 * FEELPP_INSTANTIATE_VECTOREIGEN is never defined except in vectoreigen.cpp
 * where we do the instantiate. This allows to reduce the VectorEigen
 * instantiation to the strict minimum
 */
#define FEELPP_INSTANTIATE_VECTOREIGEN 1

#include <boost/numeric/eigen/io.hpp>

#include <feel/feelalg/vectoreigen.hpp>

#if defined( FEELPP_HAS_TBB )
#include <parallel_for.h>
#include <blocked_range.h>
#endif // FEELPP_HAS_TBB

namespace Feel
{
/// \cond detail
namespace detail
{
template<typename T>
struct fake
{

};
template<>
struct fake<eigen::vector<double> >: public eigen::vector<double>
{
    fake( eigen::vector<double> v, eigen::range  r )
        :
        eigen::vector<double>( v )
    {
        boost::ignore_unused_variable_warning( r );
    }
    fake( eigen::vector<double> v, eigen::slice  r )
        :
        eigen::vector<double>( v )
    {
        boost::ignore_unused_variable_warning( r );
    }
};


#if 0
template<>
struct fake<eigen::vector<long double> >: public eigen::vector<long double>
{
    fake( eigen::vector<long double> v, eigen::range  r )
        :
        eigen::vector<long double>( v )
    {}
};
#endif //  0
template<>
struct fake<eigen::vector_range<eigen::vector<double> > >: public eigen::vector_range<eigen::vector<double> >
{
    fake( eigen::vector<double>& v, eigen::range const& r )
        :
        eigen::vector_range<eigen::vector<double> >( v, r )
    {
    }
    fake( eigen::vector<double>& v, eigen::slice const& r )
        :
        eigen::vector_range<eigen::vector<double> >( v, eigen::range( r.start(), r.size() ) )
    {
    }
};

template<>
struct fake<eigen::vector_slice<eigen::vector<double> > >: public eigen::vector_slice<eigen::vector<double> >
{
    fake( eigen::vector<double>& v, eigen::slice const& r )
        :
        eigen::vector_slice<eigen::vector<double> >( eigen::project( v, r ) )
    {
    }
    fake( eigen::vector<double>& v, eigen::range const& r )
        :
        eigen::vector_slice<eigen::vector<double> >( v, eigen::slice( r.start(),1,r.size() ) )
    {
    }

};

#if 0
template<>
struct fake<eigen::vector_range<eigen::vector<long double> > >: public eigen::vector_range<eigen::vector<long double> >
{
    fake( eigen::vector<long double>& v, eigen::range const& r )
        :
        eigen::vector_range<eigen::vector<long double> >( v, r )
    {}

};
#endif

template<typename T>
void resize( eigen::vector<T>& v, size_type s, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );

    v.resize( s );
}

template<typename T>
void resize( eigen::vector_range<eigen::vector<T> >& /*v*/, size_type /*s*/, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );
}
template<typename T>
void resize( eigen::vector_slice<eigen::vector<T> >& /*v*/, size_type /*s*/, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );
}
template<typename T>
size_type start( eigen::vector<T> const& /*v*/ )
{
    return 0;
}
template<typename T>
size_type start( eigen::vector_range<eigen::vector<T> > const& v )
{
    return v.start();
}

template<typename T>
size_type start( eigen::vector_slice<eigen::vector<T> > const& v )
{
    return v.start();
}

}
/// \endcond detail

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen()
    :
    super1(),
    _M_vec( detail::fake<Storage>( *new eigen::vector<value_type>, eigen::range() ) ),
    M_global_values_updated( false ),
    M_global_values()
{
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( size_type __s )
    :
    super1( __s ),
    _M_vec( detail::fake<Storage>( *new eigen::vector<value_type>, eigen::range() ) ),
    M_global_values_updated( false ),
    M_global_values( __s )
{
    this->init( __s, __s, false );
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( DataMap const& dm )
    :
    super1( dm ),
    _M_vec( detail::fake<Storage>( *new eigen::vector<value_type>, eigen::range() ) ),
    M_global_values_updated( false ),
    M_global_values( dm.nGlobalElements() )
{
    //this->init( dm.nGlobalElements(), dm.nMyElements(), false );
    this->init( dm.nDof(), dm.nLocalDofWithGhost(), false );
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( size_type __s, size_type __n_local )
    :
    super1( __s, __n_local ),
    _M_vec( detail::fake<Storage>( *new eigen::vector<value_type>(), eigen::range() ) ),
    M_global_values_updated( false ),
    M_global_values( this->size() )
{
    this->init( this->size(), this->localSize(), false );
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( VectorEigen const & m )
    :
    super1( m ),
    _M_vec( m._M_vec ),
    M_global_values_updated( m.M_global_values_updated ),
    M_global_values( m.M_global_values )
{
    Debug( 5600 ) << "[VectorEigen] copy constructor with range: size:" << this->size() << ", start:" << this->start() << "\n";
    Debug( 5600 ) << "[VectorEigen] copy constructor with range: size:" << this->vec().size() << "\n";
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( VectorEigen<value_type>& m, range_type const& range, DataMap const& dm )
    :
    super1( dm ),
    _M_vec( detail::fake<Storage>( m.vec(), range ) ),
    M_global_values_updated( false ),
    M_global_values( range.size() )
{
    Debug( 5600 ) << "[VectorEigen] constructor with range: size:" << range.size() << ", start:" << range.start() << "\n";
    Debug( 5600 ) << "[VectorEigen] constructor with range: size:" << _M_vec.size() << "\n";
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( eigen::vector<value_type>& m, range_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    _M_vec( detail::fake<Storage>( m, range ) ),
    M_global_values_updated( false ),
    M_global_values( range.size() )
{
    //this->init( m.size(), m.size(), false );
}
template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( VectorEigen<value_type>& m, slice_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    _M_vec( detail::fake<Storage>( m.vec(), range ) ),
    M_global_values_updated( false ),
    M_global_values( range.size() )
{
    Debug( 5600 ) << "[VectorEigen] constructor with range: size:" << range.size() << ", start:" << range.start() << "\n";
    Debug( 5600 ) << "[VectorEigen] constructor with range: size:" << _M_vec.size() << "\n";
    this->init( invalid_size_type_value, _M_vec.size(), true );

}

template <typename T, typename Storage>
VectorEigen<T,Storage>::VectorEigen( eigen::vector<value_type>& m, slice_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    _M_vec( detail::fake<Storage>( m, range ) ),
    M_global_values_updated( false ),
    M_global_values( range.size() )
{
    //this->init( m.size(), m.size(), false );
}

template <typename T, typename Storage>
VectorEigen<T,Storage>::~VectorEigen()
{
    this->clear();
}

template <typename T, typename Storage>
void
VectorEigen<T,Storage>::resize( size_type s, bool preserve )
{
    detail::resize( _M_vec, s, preserve );
}

template <typename T, typename Storage>
size_type
VectorEigen<T,Storage>::start( ) const
{
    return detail::start( _M_vec );
}

template <typename T, typename Storage>
Vector<T> &
VectorEigen<T,Storage>::operator= ( const Vector<value_type> &V )
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
        _M_vec.operator()( i ) = V( V.firstLocalIndex() + i );
        //_M_vec.operator()( i ) = V(  i );

    }

    this->outdateGlobalValues();

    return *this;
}

template <typename T, typename Storage>
void
VectorEigen<T,Storage>::init ( const size_type n,
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

    M_global_values_updated = false;
    M_global_values.resize( this->size() );

    // Initialize data structures
    detail::resize( _M_vec, this->localSize() );

    // Set the initialized flag
    this->M_is_initialized = true;

    Debug( 5600 ) << "        global size = " << n << "\n";
    Debug( 5600 ) << "        global size = " << n_local << "\n";
    Debug( 5600 ) << "        global size = " << this->size() << "\n";
    Debug( 5600 ) << "        local  size = " << this->localSize() << "\n";
    Debug( 5600 ) << "  first local index = " << this->firstLocalIndex() << "\n";
    Debug( 5600 ) << "   last local index = " << this->lastLocalIndex() << "\n";


    // Zero the components unless directed otherwise
    if ( !fast )
        this->zero();

}

template<typename T, typename Storage>
void
VectorEigen<T,Storage>::init( DataMap const& dm )
{
    super1::init( dm );
    this->init( dm.nDof(), dm.nLocalDofWithGhost(), false );
}


template<typename T, typename Storage>
void
VectorEigen<T,Storage>::init ( const size_type n,
                               const bool      fast )
{
    this->init( n,n,fast );
}

template<typename T, typename Storage>
void
VectorEigen<T,Storage>::clear()
{
    detail::resize( _M_vec, 0, false );
}

template<typename T, typename Storage>
void
VectorEigen<T,Storage>::close() const
{
}

template<typename T, typename Storage>
void
VectorEigen<T,Storage>::printMatlab( const std::string filename ) const
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
            Debug( 5600 ) << "[VectorEigen::printMatlab] adding .m extension to given file name '"
                          << filename << "'\n";
            name = filename + ".m";
        }
    }


    eigen::vector<value_type> v_local;
    this->localizeToOneProcessor ( v_local, 0 );

    if ( this->comm().rank() == 0 )
    {
        std::ofstream file_out( name.c_str() );

        FEELPP_ASSERT( file_out )( filename ).error( "[VectorEigen::printMatlab] ERROR: File cannot be opened for writing." );

        file_out << "F = [ ";
        file_out.precision( 16 );
        file_out.setf( std::ios::scientific );

        for ( size_type i = 0; i < this->size(); ++i )
        {
            file_out << v_local[i] << separator << std::endl;
        }

        file_out << "];" << std::endl;
    }
}


template <typename T, typename Storage>
void
VectorEigen<T,Storage>::localize ( Vector<T>& v_local_in ) const

{
    checkInvariant();

    VectorEigen<T,Storage>* v_local = dynamic_cast<VectorEigen<T,Storage>*>( &v_local_in );
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
    localize ( v_local->_M_vec );

#ifndef FEELPP_HAS_MPI

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );

#endif
}



template <typename T, typename Storage >
void VectorEigen<T,Storage>::localize ( Vector<T>& v_local_in,
                                        const std::vector<size_type>& ) const
{
    checkInvariant();

    // We don't support the send list.  Call the less efficient localize(v_local_in)
    localize ( v_local_in );
}



template <typename T,typename Storage>
void
VectorEigen<T,Storage>::localize ( const size_type first_local_idx,
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
    VectorEigen<T,Storage> parallel_vec;

    parallel_vec.init ( size, local_size );

    // Copy part of *this into the parallel_vec
    for ( size_type i=first_local_idx; i<=last_local_idx; i++ )
        parallel_vec.operator()( i ) = this->operator()( i );

    // localize like normal
    parallel_vec.localize ( *this, send_list );
#endif
}

template <typename T, typename Storage>
void
VectorEigen<T, Storage>::localize ( eigen::vector_range<eigen::vector<value_type> >& /*v_local*/ ) const
{
}

template <typename T, typename Storage>
void
VectorEigen<T, Storage>::localize ( eigen::vector_slice<eigen::vector<value_type> >& /*v_local*/ ) const
{
}

template <typename T, typename Storage>
void
VectorEigen<T, Storage>::localize ( eigen::vector<value_type>& v_local ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    eigen::vector<value_type> v_local_in( this->size() );

#ifdef FEELPP_HAS_MPI

    if ( this->comm().size() > 1 )
    {
        std::fill ( v_local.begin(), v_local.end(), value_type( 0. ) );
        std::fill ( v_local_in.begin(), v_local_in.end(), value_type( 0. ) );

        for ( size_type i=0; i< this->localSize(); i++ )
        {
            v_local_in[i+this->firstLocalIndex()] = _M_vec.operator[]( i );
        }

        MPI_Allreduce ( &v_local_in[0], &v_local[0], v_local.size(),
                        MPI_DOUBLE, MPI_SUM, this->comm() );
        Debug( 5600 ) << "[VectorEigen::localize] Allreduce size = " << v_local.size() << "\n";

    }

    else
    {
        FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );
        std::copy( this->begin(), this->end(), v_local.begin() );
    }

#else

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );

#endif
}



template <typename T, typename Storage>
void
VectorEigen<T,Storage>::localizeToOneProcessor ( eigen::vector<value_type>& v_local,
        const size_type pid ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    std::fill ( v_local.begin(), v_local.end(), 0. );

    eigen::vector<value_type> v_tmp( this->size() );
    std::fill ( v_tmp.begin(), v_tmp.end(), 0. );

    for ( size_type i=0; i< this->localSize(); i++ )
        v_tmp[i+this->firstLocalIndex()] = this->operator()( this->firstLocalIndex()+i );

#ifdef FEELPP_HAS_MPI

    if ( this->comm().size() > 1 )
    {
        MPI_Reduce ( &v_tmp[0], &v_local[0], v_local.size(),
                     MPI_DOUBLE, MPI_SUM, pid, this->comm() );
    }

    else
    {
        std::copy( v_tmp.begin(), v_tmp.end(), v_local.begin() );
    }

#else

    FEELPP_ASSERT ( this->localSize() == this->size() )( this->localSize() )( this->size() ).error( "invalid size in non MPI mode" );
    FEELPP_ASSERT ( pid == 0  )( pid ).error( "invalid pid in non MPI mode" );

#endif
}
template <typename T, typename Storage>
void
VectorEigen<T,Storage>::localizeToOneProcessor ( std::vector<value_type>& v_local,
        const size_type proc_id ) const
{
    eigen::vector<T> eigenvector;
    localizeToOneProcessor( eigenvector, proc_id );
    v_local.resize( eigenvector.size() );
    std::copy( eigenvector.begin(), eigenvector.end(), v_local.begin() );
}

template <typename T, typename Storage>
void
VectorEigen<T,Storage>::checkInvariant() const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
    FEELPP_ASSERT ( this->localSize() <= this->size() )
    ( this->size() )( this->localSize() ).error( "vector invalid size" );
    FEELPP_ASSERT ( _M_vec.size() == this->localSize() )
    ( _M_vec.size() )( this->localSize() ).error( "vector invalid size" );
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
template <typename T, typename Storage>
typename VectorEigen<T,Storage>::this_type
VectorEigen<T,Storage>::sqrt() const
{
    this_type _tmp( this->map() );

#if defined( FEELPP_HAS_TBB )
    tbb::parallel_for( tbb::blocked_range<size_t>( 0, this->localSize() ),
                       detail::Sqrt<this_type>( *this, _tmp ) );
#else

    for ( int i = 0; i < ( int )this->localSize(); ++i )
        _tmp[i] = math::sqrt( this->operator[]( i ) );

#endif // FEELPP_HAS_TBB

    return _tmp;
}

template <typename T, typename Storage>
typename VectorEigen<T,Storage>::this_type
VectorEigen<T,Storage>::pow( int n ) const
{
    this_type _out( this->map() );

    for ( int i = 0; i < ( int )this->localSize(); ++i )
        _out[i] = math::pow( this->operator[]( i ), n );

    return _out;
}

//
// instantiation
//
template class VectorEigen<double,eigen::vector<double> >;
template class VectorEigen<double,eigen::vector_range<eigen::vector<double> > >;
template class VectorEigen<double,eigen::vector_slice<eigen::vector<double> > >;
//template class VectorEigen<long double,eigen::vector<long double> >;
//template class VectorEigen<long double,eigen::vector_range<eigen::vector<long double> > >;


} // Feel
