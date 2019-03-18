/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-03-20

  Copyright (C) 2008, 2009, 2010 Université Joseph Fourier (Grenoble I)

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
   \file vectorublas.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-03-20
 */
/**
 * \sa vectorublas.hpp
 * FEELPP_INSTANTIATE_VECTORUBLAS is never defined except in vectorublas.cpp
 * where we do the instantiate. This allows to reduce the VectorUblas
 * instantiation to the strict minimum
 */
#define FEELPP_INSTANTIATE_VECTORUBLAS 1

#include <boost/numeric/ublas/io.hpp>

#include <feel/feelalg/vectorublas.hpp>

#if defined( FEELPP_HAS_TBB )
#include <parallel_for.h>
#include <blocked_range.h>
#endif // FEELPP_HAS_TBB

#ifdef FEELPP_HAS_HDF5
#include <feel/feelcore/hdf5.hpp>
#endif

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
struct fake<ublas::vector<double> >: public ublas::vector<double>
{
    fake( ublas::vector<double> v, ublas::range  r )
        :
        ublas::vector<double>( v )
    {
        boost::ignore_unused_variable_warning( r );
    }
    fake( ublas::vector<double> v, ublas::slice  r )
        :
        ublas::vector<double>( v )
    {
        boost::ignore_unused_variable_warning( r );
    }
    fake( size_type /*nDof*/, value_type* /*arrayDof*/ )
        :
        ublas::vector<double>( 0 )
    {
        CHECK( false ) << "not allowed";
    }

};


#if 0
template<>
struct fake<ublas::vector<long double> >: public ublas::vector<long double>
{
    fake( ublas::vector<long double> v, ublas::range  r )
        :
        ublas::vector<long double>( v )
    {}
};
#endif //  0
template<>
struct fake<ublas::vector_range<ublas::vector<double> > >: public ublas::vector_range<ublas::vector<double> >
{
    fake( ublas::vector<double>& v, ublas::range const& r )
        :
        ublas::vector_range<ublas::vector<double> >( v, r )
    {
    }
    fake( ublas::vector<double>& v, ublas::slice const& r )
        :
        ublas::vector_range<ublas::vector<double> >( v, ublas::range( r.start(), r.size() ) )
    {
    }
    fake( size_type /*nDof*/, value_type* /*arrayDof*/ )
        :
        ublas::vector_range<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::range() )
    {
        CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::range const& /*r*/ )
        :
        ublas::vector_range<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::range() )
    {
        CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::slice const& /*r*/ )
        :
        ublas::vector_range<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::range() )
    {
        CHECK( false ) << "not allowed";
    }

};

template<>
struct fake<ublas::vector_slice<ublas::vector<double> > >: public ublas::vector_slice<ublas::vector<double> >
{
    fake( ublas::vector<double>& v, ublas::slice const& r )
        :
        ublas::vector_slice<ublas::vector<double> >( ublas::project( v, r ) )
    {
    }
    fake( ublas::vector<double>& v, ublas::range const& r )
        :
        ublas::vector_slice<ublas::vector<double> >( v, ublas::slice( r.start(),1,r.size() ) )
    {
    }
    fake( size_type /*nDof*/, value_type* /*arrayDof*/ )
        :
        ublas::vector_slice<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::slice() )
    {
        CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::range const& /*r*/ )
        :
        ublas::vector_slice<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::slice() )
    {
        CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::slice const& /*r*/ )
        :
        ublas::vector_slice<ublas::vector<double> >( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::slice() )
    {
        CHECK( false ) << "not allowed";
    }
};

template<>
struct fake<ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > >: public ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >
{
    fake( ublas::vector<double>& /*v*/, ublas::slice const& /*r*/ )
        :
        ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >( Feel::detail::shallow_array_adaptor<double>( 0, nullptr ) )
    {
        CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double>& /*v*/, ublas::range const& /*r*/ )
        :
        ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >( Feel::detail::shallow_array_adaptor<double>( 0, nullptr ) )
    {
        CHECK( false ) << "not allowed";
    }
    fake( size_type nDof, value_type* arrayDof )
        :
        ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >( Feel::detail::shallow_array_adaptor<double>( nDof, arrayDof ) )
    {}
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::range const& /*r*/ )
        :
        ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >( Feel::detail::shallow_array_adaptor<double>( 0, nullptr ) )
     {
         CHECK( false ) << "not allowed"; //CHECK
     }
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::slice const& /*r*/ )
        :
        ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >( Feel::detail::shallow_array_adaptor<double>( 0, nullptr ) )
     {
         CHECK( false ) << "not allowed"; //CHECK
     }

};

template<>
struct fake<ublas::vector_range<ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > > > :
        public ublas::vector_range<ublas::vector<double, Feel::detail::shallow_array_adaptor<double>> >
{
    typedef ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > vector_type;
    typedef ublas::vector_range<vector_type > super_type;

    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& v, ublas::range const& r )
        :
        super_type( v, r )
    {}

    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& /*v*/, ublas::slice const& /*r*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::range() )
    {
        CHECK( false ) << "not allowed";
    }

    fake( ublas::vector<double>& /*v*/, ublas::range const& /*r*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::range() )
    {
         CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double>& /*v*/, ublas::slice const& /*r*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::range() )
    {
         CHECK( false ) << "not allowed";
    }

    fake( size_type /*nDof*/, value_type* /*arrayDof*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::range() )
    {
        CHECK( false ) << "not allowed";
    }
};

template<>
struct fake<ublas::vector_slice<ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > > > :
        public ublas::vector_slice<ublas::vector<double, Feel::detail::shallow_array_adaptor<double>> >
{
    typedef ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > vector_type;
    typedef ublas::vector_slice<vector_type > super_type;

    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& v, ublas::slice const& r )
        :
        super_type( v, r ) // super_type( ublas::project( v, r ) ) )
    {}
    fake( ublas::vector<double, Feel::detail::shallow_array_adaptor<double> >& v, ublas::range const& r )
        :
        super_type( v, ublas::slice( r.start(),1,r.size() ) )
    {}

    fake( ublas::vector<double>& /*v*/, ublas::range const& /*r*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::slice() )
    {
         CHECK( false ) << "not allowed";
    }
    fake( ublas::vector<double>& /*v*/, ublas::slice const& /*r*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::slice() )
    {
         CHECK( false ) << "not allowed";
    }

    fake( size_type /*nDof*/, value_type* /*arrayDof*/ )
        :
        super_type( *std::shared_ptr<vector_type>( new vector_type(0) ), ublas::slice() )
    {
        CHECK( false ) << "not allowed";
    }
};


#if 0
template<>
struct fake<ublas::vector_range<ublas::vector<long double> > >: public ublas::vector_range<ublas::vector<long double> >
{
    fake( ublas::vector<long double>& v, ublas::range const& r )
        :
        ublas::vector_range<ublas::vector<long double> >( v, r )
    {}

};
#endif

template<typename T>
void resize( ublas::vector<T>& v, size_type s, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );

    v.resize( s );
}

template<typename T>
void resize( ublas::vector_range<ublas::vector<T> >& /*v*/, size_type /*s*/, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );
}
template<typename T>
void resize( ublas::vector_slice<ublas::vector<T> >& /*v*/, size_type /*s*/, bool preserve = true )
{
    boost::ignore_unused_variable_warning( preserve );
}
template<typename T>
void resize( ublas::vector<T,Feel::detail::shallow_array_adaptor<T> >& /*v*/, size_type /*s*/, bool preserve = true )
{
}
template<typename T>
void resize( ublas::vector_range<ublas::vector<T,Feel::detail::shallow_array_adaptor<T> > >& /*v*/, size_type /*s*/, bool preserve = true )
{
}
template<typename T>
void resize( ublas::vector_slice<ublas::vector<T,Feel::detail::shallow_array_adaptor<T> > >& /*v*/, size_type /*s*/, bool preserve = true )
{
}

template<typename T>
size_type start( ublas::vector<T> const& /*v*/ )
{
    return 0;
}
template<typename T>
size_type start( ublas::vector_range<ublas::vector<T> > const& v )
{
    return v.start();
}

template<typename T>
size_type start( ublas::vector_slice<ublas::vector<T> > const& v )
{
    return v.start();
}
template<typename T>
size_type start( ublas::vector<T,Feel::detail::shallow_array_adaptor<T>> const& /*v*/ )
{
    return 0;
}
template<typename T>
size_type start( ublas::vector_range<ublas::vector<T,Feel::detail::shallow_array_adaptor<T> > > const& v )
{
    return v.start();
}
template<typename T>
size_type start( ublas::vector_slice<ublas::vector<T,Feel::detail::shallow_array_adaptor<T> > > const& v )
{
    return v.start();
}

}
/// \endcond detail

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas()
    :
    super1(),
    M_vec( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) )
{
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( size_type __s )
    :
    super1( __s, Environment::worldCommSeqPtr() ),
    M_vec( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) )
{
    this->init( __s, __s, false );
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( *std::shared_ptr<ublas::vector<value_type>>( new ublas::vector<value_type> ), ublas::range() ) )
{
    //this->init( dm.nGlobalElements(), dm.nMyElements(), false );
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false );
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( size_type __s, size_type __n_local )
    :
    super1( __s, __n_local, Environment::worldCommSeqPtr() ),
    M_vec( detail::fake<Storage>( *new ublas::vector<value_type>(), ublas::range() ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( *std::make_shared<ublas::vector<value_type>>(), ublas::range() ) )
{
    this->init( this->size(), this->localSize(), false );
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas const & m )
    :
    super1( m ),
    M_vec( m.M_vec ),
    M_vecNonContiguousGhosts( m.M_vecNonContiguousGhosts )
{
    DVLOG(2) << "[VectorUblas] copy constructor with range: size:" << this->size() << ", start:" << this->start() << "\n";
    DVLOG(2) << "[VectorUblas] copy constructor with range: size:" << this->vec().size() << "\n";
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas<value_type>& m, range_type const& range, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vec(), ublas::range(0,0) ) )
{
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << range.size() << ", start:" << range.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << M_vec.size() << "\n";
}
template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas<value_type>& m, range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), rangeActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vec(), rangeGhost ) )
{
    CHECK( is_range_vector ) << "is not a range vector";
    DVLOG(2) << "[VectorUblas] constructor with active range: size:" << rangeActive.size() << ", start:" << rangeActive.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with ghost range: size:" << rangeGhost.size() << ", start:" << rangeGhost.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << M_vec.size() << "\n";
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m,
                                     range_type const& rangeActive, range_type const& rangeGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), rangeActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vecNonContiguousGhosts(), rangeGhost ) )
{
    CHECK( is_range_vector ) << "is not a range vector";
    CHECK( is_shallow_array_adaptor_vector ) << "is_shallow_array_adaptor_vector";
}


template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas<value_type>& m, slice_type const& range, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vec(), ublas::range(0,0) ) )
{
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << range.size() << ", start:" << range.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << M_vec.size() << "\n";
}
template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas<value_type>& m, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), sliceActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vec(), sliceGhost ) )
{
    DVLOG(2) << "[VectorUblas] constructor with active range: size:" << sliceActive.size() << ", start:" << sliceActive.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with ghost range: size:" << sliceGhost.size() << ", start:" << sliceGhost.start() << "\n";
}
template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( typename VectorUblas<value_type>::shallow_array_adaptor::type& m, slice_type const& sliceActive, slice_type const& sliceGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m.vec(), sliceActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vecNonContiguousGhosts(), sliceGhost ) )
{
    DVLOG(2) << "[VectorUblas] constructor with active range: size:" << sliceActive.size() << ", start:" << sliceActive.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with ghost range: size:" << sliceGhost.size() << ", start:" << sliceGhost.start() << "\n";
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( ublas::vector<value_type>& m, range_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    M_vec( detail::fake<Storage>( m, range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m, range_type(0,0) ) )
{
    //this->init( m.size(), m.size(), false );
}
template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( ublas::vector<value_type>& m, range_type const& range, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m, range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m, range_type(0,0) ) )
{}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( VectorUblas<value_type>& m, slice_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    M_vec( detail::fake<Storage>( m.vec(), range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m.vec(), slice_type(0,1,0) ) )
{
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << range.size() << ", start:" << range.start() << "\n";
    DVLOG(2) << "[VectorUblas] constructor with range: size:" << M_vec.size() << "\n";
    this->init( invalid_size_type_value, M_vec.size(), true );

}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( ublas::vector<value_type>& m, slice_type const& range, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( m, range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m, slice_type(0,1,0) ) )
{}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( ublas::vector<value_type>& m, slice_type const& range )
    :
    super1( invalid_size_type_value, range.size() ),
    M_vec( detail::fake<Storage>( m, range ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( m, slice_type(0,1,0) ) )
{
    //this->init( m.size(), m.size(), false );
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( ublas::vector<value_type>& mActive, slice_type const& sliceActive,
                                     ublas::vector<value_type>& mGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( mActive, sliceActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( mGhost, sliceGhost ) )
{
    CHECK( is_slice_vector ) << "is not a slice vector";
}

template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( typename this_type::shallow_array_adaptor::subtype& mActive, slice_type const& sliceActive,
                                     typename this_type::shallow_array_adaptor::subtype& mGhost, slice_type const& sliceGhost, datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( mActive, sliceActive ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( mGhost, sliceGhost ) )
{
    CHECK( is_slice_vector ) << "is not a slice vector";
}


template <typename T, typename Storage>
VectorUblas<T,Storage>::VectorUblas( size_type nActiveDof, value_type* arrayActiveDof,
                                     size_type nGhostDof, value_type* arrayGhostDof,
                                     datamap_ptrtype const& dm )
    :
    super1( dm ),
    M_vec( detail::fake<Storage>( nActiveDof, arrayActiveDof ) ),
    M_vecNonContiguousGhosts( detail::fake<Storage>( nGhostDof, arrayGhostDof ) )
{}


template <typename T, typename Storage>
VectorUblas<T,Storage>::~VectorUblas()
{
    this->clear();
}

template <typename T, typename Storage>
void
VectorUblas<T,Storage>::resize( size_type s, bool preserve )
{
    detail::resize( M_vec, s, preserve );
}

template <typename T, typename Storage>
size_type
VectorUblas<T,Storage>::start( ) const
{
    return detail::start( M_vec );
}

template <typename T, typename Storage>
size_type
VectorUblas<T,Storage>::startNonContiguousGhosts( ) const
{
    return detail::start( M_vecNonContiguousGhosts );
}

template <typename T, typename Storage>
Vector<T> &
VectorUblas<T,Storage>::operator= ( const Vector<value_type> &v )
{
    if ( this == &v )
        return *this;

    if ( !this->map().isCompatible( v.map() ) )
    {
        this->setMap( v.mapPtr() );
        this->resize( this->map().nLocalDofWithGhost() );
    }

#if !defined(NDEBUG)
    checkInvariant();

    FEELPP_ASSERT( this->localSize() == v.localSize() &&
                   this->map().nLocalDofWithoutGhost() == v.map().nLocalDofWithoutGhost() )
    ( this->localSize() )( this->map().nLocalDofWithoutGhost() )
    ( this->vec().size() )
    ( v.localSize() )( v.map().nLocalDofWithoutGhost() ).warn( "may be vector invalid  copy" );
#endif

    typedef VectorUblas<T> the_vector_ublas_type;
    typedef typename the_vector_ublas_type::range::type the_vector_ublas_range_type;
    typedef typename the_vector_ublas_type::slice::type the_vector_ublas_slice_type;
    typedef typename VectorUblas<T>::shallow_array_adaptor::type the_vector_ublas_extarray_type;
    typedef typename the_vector_ublas_extarray_type::range::type the_vector_ublas_extarray_range_type;
    typedef typename the_vector_ublas_extarray_type::slice::type the_vector_ublas_extarray_slice_type;

    const the_vector_ublas_type * vecUblas = dynamic_cast<the_vector_ublas_type const*>( &v );
    if ( vecUblas )
    {
        this->assignWithUblasImpl( *vecUblas );
        return *this;
    }
    const the_vector_ublas_range_type * vecUblasRange = dynamic_cast<the_vector_ublas_range_type const*>( &v );
    if ( vecUblasRange )
    {
        this->assignWithUblasImpl( *vecUblasRange );
        return *this;
    }
    const the_vector_ublas_slice_type * vecUblasSlice = dynamic_cast<the_vector_ublas_slice_type const*>( &v );
    if ( vecUblasSlice )
    {
        this->assignWithUblasImpl( *vecUblasSlice );
        return *this;
    }
    const the_vector_ublas_extarray_type * vecUblasExtArray = dynamic_cast<the_vector_ublas_extarray_type const*>( &v );
    if ( vecUblasExtArray )
    {
        this->assignWithUblasImpl( *vecUblasExtArray );
        return *this;
    }
    const the_vector_ublas_extarray_range_type * vecUblasExtArrayRange = dynamic_cast<the_vector_ublas_extarray_range_type const*>( &v );
    if ( vecUblasExtArrayRange )
    {
        this->assignWithUblasImpl( *vecUblasExtArrayRange );
        return *this;
    }
    const the_vector_ublas_extarray_slice_type * vecUblasExtArraySlice = dynamic_cast<the_vector_ublas_extarray_slice_type const*>( &v );
    if ( vecUblasExtArraySlice )
    {
        this->assignWithUblasImpl( *vecUblasExtArraySlice );
        return *this;
    }
#if FEELPP_HAS_PETSC
    const VectorPetsc<T> * vecPetsc = dynamic_cast<VectorPetsc<T> const*>( &v );
    if ( vecPetsc && !is_slice_vector )
    {
        //toPETSc( *this ).operator=( *vecPetsc );
        toPETScPtr( *this )->operator=( *vecPetsc );
        return *this;
    }
#endif

    // default operator=
    for ( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) = v(  i );

    return *this;
}
template <typename T, typename Storage>
Vector<T>&
VectorUblas<T,Storage>::operator= ( const this_type &V )
{
    return this->operator=( *dynamic_cast< Vector<value_type> const* >( &V ) );
}

template <typename T, typename Storage>
void
VectorUblas<T,Storage>::init ( const size_type n,
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
    detail::resize( M_vec, this->localSize() );

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

template<typename T, typename Storage>
void
VectorUblas<T,Storage>::init( datamap_ptrtype const& dm )
{
    super1::init( dm );
    this->init( dm->nDof(), dm->nLocalDofWithGhost(), false );
}


template<typename T, typename Storage>
void
VectorUblas<T,Storage>::init ( const size_type n,
                               const bool      fast )
{
    this->init( n,n,fast );
}

template<typename T, typename Storage>
void
VectorUblas<T,Storage>::clear()
{
    if ( !is_extarray_vector )
    {
        detail::resize( M_vec, 0, false );
        if ( has_non_contiguous_ghosts )
            detail::resize( M_vecNonContiguousGhosts, 0, false );
    }
}

template<typename T, typename Storage>
void
VectorUblas<T,Storage>::close() const
{
}

template<typename T, typename Storage>
void
VectorUblas<T,Storage>::add( const T& a, const Vector<T>& v )
{
    checkInvariant();

    typedef VectorUblas<T> the_vector_ublas_type;
    typedef typename the_vector_ublas_type::range::type the_vector_ublas_range_type;
    typedef typename the_vector_ublas_type::slice::type the_vector_ublas_slice_type;
    typedef typename VectorUblas<T>::shallow_array_adaptor::type the_vector_ublas_extarray_type;
    typedef typename the_vector_ublas_extarray_type::range::type the_vector_ublas_extarray_range_type;
    typedef typename the_vector_ublas_extarray_type::slice::type the_vector_ublas_extarray_slice_type;

    const the_vector_ublas_type * vecUblas = dynamic_cast<the_vector_ublas_type const*>( &v );
    if ( vecUblas )
    {
        this->addWithUblasImpl( a,*vecUblas );
        return;
    }
    const the_vector_ublas_range_type * vecUblasRange = dynamic_cast<the_vector_ublas_range_type const*>( &v );
    if ( vecUblasRange )
    {
        this->addWithUblasImpl( a,*vecUblasRange );
        return;
    }
    const the_vector_ublas_slice_type * vecUblasSlice = dynamic_cast<the_vector_ublas_slice_type const*>( &v );
    if ( vecUblasSlice )
    {
        this->addWithUblasImpl( a,*vecUblasSlice );
        return;
    }
    const the_vector_ublas_extarray_type * vecUblasExtArray = dynamic_cast<the_vector_ublas_extarray_type const*>( &v );
    if ( vecUblasExtArray )
    {
        this->addWithUblasImpl( a,*vecUblasExtArray );
        return;
    }
    const the_vector_ublas_extarray_range_type * vecUblasExtArrayRange = dynamic_cast<the_vector_ublas_extarray_range_type const*>( &v );
    if ( vecUblasExtArrayRange )
    {
        this->addWithUblasImpl( a,*vecUblasExtArrayRange );
        return;
    }
    const the_vector_ublas_extarray_slice_type * vecUblasExtArraySlice = dynamic_cast<the_vector_ublas_extarray_slice_type const*>( &v );
    if ( vecUblasExtArraySlice )
    {
        this->addWithUblasImpl( a,*vecUblasExtArraySlice );
        return;
    }
    
#if FEELPP_HAS_PETSC
    const VectorPetsc<T> * vecPetsc = dynamic_cast<VectorPetsc<T> const*>( &v );
    if ( vecPetsc && !is_slice_vector )
    {
        //toPETSc( *this ).add( a,*vecPetsc );
        toPETScPtr( *this )->add( a,*vecPetsc );
        return;
    }
#endif

    // default add operator
    for ( size_type i = 0; i < this->localSize(); ++i )
        this->operator()( i ) += a*v(  i );

    return;
}

template<typename T, typename Storage>
void
VectorUblas<T,Storage>::printMatlab( const std::string filename, bool renumber ) const
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
    this->localizeToOneProcessor ( v_local, 0 );

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


template <typename T, typename Storage>
void
VectorUblas<T,Storage>::localize ( Vector<T>& v_local_in ) const

{
    checkInvariant();

    VectorUblas<T,Storage>* v_local = dynamic_cast<VectorUblas<T,Storage>*>( &v_local_in );
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



template <typename T, typename Storage >
void VectorUblas<T,Storage>::localize ( Vector<T>& v_local_in,
                                        const std::vector<size_type>& ) const
{
    checkInvariant();

    // We don't support the send list.  Call the less efficient localize(v_local_in)
    localize ( v_local_in );
}



template <typename T,typename Storage>
void
VectorUblas<T,Storage>::localize ( const size_type first_local_idx,
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
    VectorUblas<T,Storage> parallel_vec;

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
VectorUblas<T, Storage>::localize ( ublas::vector_range<ublas::vector<value_type> >& /*v_local*/ ) const
{
}

template <typename T, typename Storage>
void
VectorUblas<T, Storage>::localize ( ublas::vector_slice<ublas::vector<value_type> >& /*v_local*/ ) const
{
}

template <typename T, typename Storage>
void
VectorUblas<T, Storage>::localize ( ublas::vector<value_type>& v_local ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    ublas::vector<value_type> v_local_in( this->size() );

#ifdef FEELPP_HAS_MPI

    if ( this->comm().size() > 1 )
    {
        std::fill ( v_local.begin(), v_local.end(), value_type( 0. ) );
        std::fill ( v_local_in.begin(), v_local_in.end(), value_type( 0. ) );

        for ( size_type i=0; i< this->localSize(); i++ )
        {
            v_local_in[i+this->firstLocalIndex()] = M_vec.operator[]( i );
        }

        MPI_Allreduce ( &v_local_in[0], &v_local[0], v_local.size(),
                        MPI_DOUBLE, MPI_SUM, this->comm() );
        DVLOG(2) << "[VectorUblas::localize] Allreduce size = " << v_local.size() << "\n";

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
VectorUblas<T,Storage>::localizeToOneProcessor ( ublas::vector<value_type>& v_local,
        const size_type pid ) const
{
    checkInvariant();

    v_local.resize( this->size() );
    std::fill ( v_local.begin(), v_local.end(), 0. );

    ublas::vector<value_type> v_tmp( this->size() );
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
VectorUblas<T,Storage>::localizeToOneProcessor ( std::vector<value_type>& v_local,
        const size_type proc_id ) const
{
    ublas::vector<T> ublasvector;
    localizeToOneProcessor( ublasvector, proc_id );
    v_local.resize( ublasvector.size() );
    std::copy( ublasvector.begin(), ublasvector.end(), v_local.begin() );
}

template <typename T, typename Storage>
void
VectorUblas<T,Storage>::checkInvariant() const
{
    DCHECK ( this->isInitialized() ) <<  "vector not initialized" ;
    DCHECK ( this->localSize() <= this->size() ) << "vector invalid size: " << this->size() << "," << this->localSize();
    DCHECK ( this->localSize() == (M_vec.size()+M_vecNonContiguousGhosts.size() ) ) << "vector invalid size: " << M_vec.size() << ","
                                                                                      << M_vecNonContiguousGhosts.size() << ","
                                                                                      << this->localSize();
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
typename VectorUblas<T,Storage>::this_type
VectorUblas<T,Storage>::sqrt() const
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

template <typename T, typename Storage>
typename VectorUblas<T,Storage>::this_type
VectorUblas<T,Storage>::pow( int n ) const
{
    this_type _out( this->mapPtr() );

    for ( int i = 0; i < ( int )this->localSize(); ++i )
        _out[i] = math::pow( this->operator[]( i ), n );

    return _out;
}

template<typename T, typename Storage>
typename VectorUblas<T,Storage>::value_type
VectorUblas<T,Storage>::dot( Vector<T> const& v ) const
{
    typedef VectorUblas<T> the_vector_ublas_type;
    typedef typename the_vector_ublas_type::range::type the_vector_ublas_range_type;
    typedef typename the_vector_ublas_type::slice::type the_vector_ublas_slice_type;
    typedef typename VectorUblas<T>::shallow_array_adaptor::type the_vector_ublas_extarray_type;
    typedef typename the_vector_ublas_extarray_type::range::type the_vector_ublas_extarray_range_type;
    typedef typename the_vector_ublas_extarray_type::slice::type the_vector_ublas_extarray_slice_type;

    const the_vector_ublas_type * vecUblas = dynamic_cast<the_vector_ublas_type const*>( &v );
    if ( vecUblas )
    {
        return this->dotWithUblasImpl( *vecUblas );
    }
    const the_vector_ublas_range_type * vecUblasRange = dynamic_cast<the_vector_ublas_range_type const*>( &v );
    if ( vecUblasRange )
    {
        return this->dotWithUblasImpl( *vecUblasRange );
    }
    const the_vector_ublas_slice_type * vecUblasSlice = dynamic_cast<the_vector_ublas_slice_type const*>( &v );
    if ( vecUblasSlice )
    {
        return this->dotWithUblasImpl( *vecUblasSlice );
    }
    const the_vector_ublas_extarray_type * vecUblasExtArray = dynamic_cast<the_vector_ublas_extarray_type const*>( &v );
    if ( vecUblasExtArray )
    {
        return this->dotWithUblasImpl( *vecUblasExtArray );
    }
    const the_vector_ublas_extarray_range_type * vecUblasExtArrayRange = dynamic_cast<the_vector_ublas_extarray_range_type const*>( &v );
    if ( vecUblasExtArrayRange )
    {
        return this->dotWithUblasImpl( *vecUblasExtArrayRange );
    }
    const the_vector_ublas_extarray_slice_type * vecUblasExtArraySlice = dynamic_cast<the_vector_ublas_extarray_slice_type const*>( &v );
    if ( vecUblasExtArraySlice )
    {
        return this->dotWithUblasImpl( *vecUblasExtArraySlice );
    }

#if FEELPP_HAS_PETSC
    const VectorPetsc<T> * vecPetsc = dynamic_cast<VectorPetsc<T> const*>( &v );
    if ( vecPetsc && !is_slice_vector )
    {
        return toPETScPtr( *this )->dot( *vecPetsc );
    }
#endif

    // default dot operator
    value_type localResult = 0;
    for ( size_type k = 0; k < this->map().nLocalDofWithoutGhost(); ++k )
    {
        localResult += this->operator()(k)*v(k);
    }
    value_type globalResult = localResult;
#ifdef FEELPP_HAS_MPI
    if ( this->comm().size() > 1 )
        mpi::all_reduce( this->comm(), localResult, globalResult, std::plus<value_type>() );
#endif
    return globalResult;
}

#ifdef FEELPP_HAS_HDF5

#if 0
template<typename T, typename Storage>
void
VectorUblas<T,Storage>::ioHDF5( bool isLoad, std::string const& filename )
{
    bool useTransposedStorage = true;
    const int dimsComp0 = (useTransposedStorage)? 1 : 0;
    const int dimsComp1 = (useTransposedStorage)? 0 : 1;
    std::vector<uint> sizeValues(1, this->map().nLocalDofWithGhost() );
    hsize_t dims[2];
    dims[dimsComp0] = this->comm().localSize();dims[dimsComp1] = 1;
    hsize_t dims2[2];

    dims2[dimsComp0] = sizeValues.size();dims2[dimsComp1] = 1;
    hsize_t offset[2];
    offset[dimsComp0] = this->comm().localRank(); offset[dimsComp1] = 0;

    size_type nLocalDofWithGhostTotal = 0;
    for (rank_type p=0 ; p < this->comm().localSize() ; ++p )
        nLocalDofWithGhostTotal += this->map().nLocalDofWithGhost( p );

    hsize_t dimsElt[2];
    dimsElt[dimsComp0] = nLocalDofWithGhostTotal;//this->map().nDof();
    dimsElt[dimsComp1] = 1;

    hsize_t dimsElt2[2];
    dimsElt2[dimsComp0] = this->map().nLocalDofWithGhost();
    dimsElt2[dimsComp1] = 1;
    hsize_t offsetElt[2];
    size_type offsetCount = 0;
    for (rank_type p=0 ; p < this->comm().localRank() ; ++p )
        offsetCount += this->map().nLocalDofWithGhost( p );

    offsetElt[dimsComp0] = offsetCount;
    offsetElt[dimsComp1] = 0;

    HDF5 hdf5;
    hdf5.openFile( filename, this->comm().localComm(), isLoad );

    if ( isLoad )
    {
        if ( false )
        {
            std::vector<uint> sizeValuesReload( 1 );
            hdf5.openTable( "size",dims );
            hdf5.read( "size", H5T_NATIVE_UINT, dims2, offset, sizeValuesReload.data() );
            hdf5.closeTable( "size" );
            CHECK( sizeValuesReload[0] == this->map().nLocalDofWithGhost() ) << "error : must be equal "  << sizeValuesReload[0] << " " << this->map().nLocalDofWithGhost();
        }

        hdf5.openTable( "element",dimsElt );
        if ( this->map().nLocalDofWithGhost() > 0 )
            hdf5.read( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, &(M_vec[0]) );
        else
        {
            double uselessValue=0;
            hdf5.read( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, &uselessValue );
        }
        hdf5.closeTable( "element" );
    }
    else // save
    {
        if ( false )
        {
            // create size tab
            hdf5.createTable( "size", H5T_NATIVE_UINT, dims );
            hdf5.write( "size", H5T_NATIVE_UINT, dims2, offset, sizeValues.data() );
            hdf5.closeTable( "size" );
        }

        // create double tab
        hdf5.createTable( "element", H5T_NATIVE_DOUBLE, dimsElt );
        //hdf5.write( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, M_vec.data().begin()) );
        if ( this->map().nLocalDofWithGhost() > 0 )
            hdf5.write( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, &(M_vec[0]) );
        hdf5.closeTable( "element" );
    }
    hdf5.closeFile();

}
#else
template<typename T, typename Storage>
void
VectorUblas<T,Storage>::ioHDF5( bool isLoad, std::string const& filename, std::string tableName, bool appendMode )
{
    bool useTransposedStorage = true;
    const int dimsComp0 = (useTransposedStorage)? 1 : 0;
    const int dimsComp1 = (useTransposedStorage)? 0 : 1;
    DataMap const& dm = this->map();
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

    if ( isLoad )
    {
        hdf5.openFile( filename, this->comm().localComm(), isLoad );

        if ( false )
        {
            std::vector<uint> sizeValuesReload( 1 );
            hdf5.openTable( "size",dims );
            hdf5.read( "size", H5T_NATIVE_UINT, dims2, offset, sizeValuesReload.data() );
            hdf5.closeTable( "size" );
            CHECK( sizeValuesReload[0] == dm.nLocalDofWithoutGhost() ) << "error : must be equal "  << sizeValuesReload[0] << " " << dm.nLocalDofWithoutGhost();
        }

        hdf5.openTable( tableName, dimsElt );
        if ( dm.nLocalDofWithoutGhost() > 0 )
            hdf5.read( tableName, H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, dataStorage.data()/*&(M_vec[0])*/ );
        else
        {
            double uselessValue=0;
            hdf5.read( tableName, H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, &uselessValue );
        }
        hdf5.closeTable( tableName );

        for ( size_type k=0;k<dataStorage.size();++k )
            this->set( k, dataStorage[k] );
        sync( *this );
    }
    else // save
    {
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
    }
    hdf5.closeFile();

}

#endif
template<typename T, typename Storage>
void
VectorUblas<T,Storage>::saveHDF5( std::string const& filename, std::string tableName, bool appendMode )
{
    this->ioHDF5( false, filename, tableName, appendMode);
}
template<typename T, typename Storage>
void
VectorUblas<T,Storage>::loadHDF5( std::string const& filename, std::string tableName )
{
    this->ioHDF5( true, filename, tableName );
}
#endif



template<typename T, typename Storage>
typename VectorUblas<T,Storage>::shallow_array_adaptor::type
VectorUblas<T,Storage>::createView( Vector<T> const& vec )
{
    typedef typename VectorUblas<T,Storage>::shallow_array_adaptor::type ret_view_type;
    datamap_ptrtype dmVec = vec.mapPtr();
#if FEELPP_HAS_PETSC
    VectorPetsc<value_type> * vecPetsc = const_cast< VectorPetsc<value_type> *>( dynamic_cast< VectorPetsc<value_type> const*>( &vec ) );
    if ( vecPetsc )
    {
        //VectorPetscMPI<value_type> * vecPetsc = const_cast< VectorPetscMPI<value_type> *>( dynamic_cast< VectorPetscMPI<value_type> const*>( &vec ) );
        size_type nActiveDof = dmVec->nLocalDofWithoutGhost();
        value_type* arrayActiveDof = (nActiveDof>0)? std::addressof( (*vecPetsc)( 0  ) ) : nullptr;
        size_type nGhostDof = dmVec->nLocalGhosts();
        value_type* arrayGhostDof = (nGhostDof>0)? std::addressof( (*vecPetsc)( nActiveDof ) ) : nullptr;
        ret_view_type u( nActiveDof, arrayActiveDof, nGhostDof, arrayGhostDof, dmVec );
        return u;
    }
#endif

    CHECK( false ) << "fails in vector cast (only implemented for Petsc vector)";
    ret_view_type u;
    return u;
}

//
// instantiation
//
template class VectorUblas<double,ublas::vector<double> >;
template class VectorUblas<double,ublas::vector_range<ublas::vector<double> > >;
template class VectorUblas<double,ublas::vector_slice<ublas::vector<double> > >;
template class VectorUblas<double,ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > >;
template class VectorUblas<double,ublas::vector_range<ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > > >;
template class VectorUblas<double,ublas::vector_slice<ublas::vector<double, Feel::detail::shallow_array_adaptor<double> > > >;
//template class VectorUblas<long double,ublas::vector<long double> >;
//template class VectorUblas<long double,ublas::vector_range<ublas::vector<long double> > >;


} // Feel
