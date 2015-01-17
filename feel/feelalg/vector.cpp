/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
Vector<T>::Vector( WorldComm const& _worldComm ) :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map ( new datamap_type( _worldComm ) )
{}


    /*template <typename T>
Vector<T>::Vector ( datamap_type const& dm ) :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map ( new datamap_type(dm) )
{}
    */
template <typename T>
Vector<T>::Vector( datamap_ptrtype const& dm ) :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map ( dm )
{}


template <typename T>
Vector<T>::Vector ( const size_type n, WorldComm const& _worldComm )
    :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map( new datamap_type(n, n, _worldComm) )
{}



template <typename T>
Vector<T>::Vector ( const size_type n,
                    const size_type n_local,
                    WorldComm const& _worldComm )
    :
    M_is_closed( false ),
    M_is_initialized( false ),
    M_map( new datamap_type(n, n_local, _worldComm) )
{}

template <typename T>
Vector<T>::Vector ( Vector const& v )
    :
    M_is_closed( v.M_is_closed ),
    M_is_initialized( v.M_is_initialized ),
    M_map( v.M_map )
{
}
template <typename T>

Vector<T>::~Vector ()
{
    clear ();
}



template <typename T>

Vector<T> & Vector<T>::operator= ( const T )
{
    //  error();

    return *this;
}


template <typename T>
void
Vector<T>::init ( const size_type n,
                  const size_type nl,
                  const bool fast )
{
    boost::ignore_unused_variable_warning( fast );
    //deja fait dans le constructeur!!
    //M_map=DataMap( n, nl );
}

template <typename T>
void
Vector<T>::init ( const size_type n,
                  const bool fast )
{
    this->init( n, n, fast );
}

template <typename T>

Vector<T> & Vector<T>::operator= ( const Vector<T>& v )
{
    if ( this != &v )
    {
        M_map = v.mapPtr();

        for ( size_type i = 0; i < this->map().nLocalDofWithGhost(); ++i )
        {
            this->set( i,  v( v.firstLocalIndex() + i ) );
        }
    }

    return *this;
}



template <typename T>

Vector<T> & Vector<T>::operator= ( const std::vector<T>& v )
{
    M_map.reset( new datamap_type( v.size(), 0 ) );
    return *this;
}



template <typename T>

void Vector<T>::clear ()
{
    M_is_closed      = false;
    M_is_initialized = false;
}
template<typename T>
void
Vector<T>::addVector ( const Vector<T>& V_in,
                       const MatrixShell<T>& A_in )
{
    A_in.multVector( V_in, *this );
}

template<typename T>
void
Vector<T>::addVector ( const boost::shared_ptr<Vector<T> >& V_in,
                       const boost::shared_ptr<MatrixShell<T> >& A_in )
{
    A_in->multVector( *V_in, *this );
}

template<typename T>
int
Vector<T>::reciprocal ()
{
    LOG(WARNING) << "Invalid call to reciprocal. Not implement in Vector base class";
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
int
Vector<float>::compare ( const Vector<float> &other_vector,
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

// Full specialization for double datatypes
template <>
int
Vector<double>::compare ( const Vector<double> &other_vector,
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
int
Vector<std::complex<double>>::compare ( const Vector<std::complex<double>> &other_vector,
                                        const real_type threshold ) const
{
    CHECK ( this->isInitialized() ) << "vector not initialized";
    CHECK ( other_vector.isInitialized() ) << "other vector not initialized";
    CHECK ( this->firstLocalIndex() == other_vector.firstLocalIndex() ) << "invalid index";
    CHECK ( this->lastLocalIndex()  == other_vector.lastLocalIndex() ) << "invalid index";

    int rvalue     = -1;
    size_type i = firstLocalIndex();

    do
    {
        if ( ( std::abs( ( *this )( i ).real() - other_vector( i ).real() ) > threshold ) ||
                ( std::abs( ( *this )( i ).imag() - other_vector( i ).imag() ) > threshold ) )
            rvalue = i;

        else
            i++;
    }
    while ( rvalue==-1 && i<this->lastLocalIndex() );

    return rvalue;
}


template <typename T>

void Vector<T>::print( std::ostream& os ) const
{
    FEELPP_ASSERT ( this->isInitialized() ).error( "vector not initialized" );
    os << "Size\tglobal =  " << this->size()
       << "\t\tlocal =  " << this->localSize() << std::endl;

    os << "#\tValue" << std::endl;

    for ( size_type i=this->firstLocalIndex(); i<this->lastLocalIndex(); i++ )
        os << i << "\t" << ( *this )( i ) << std::endl;
}
template <typename T>
void Vector<T>::localize( Vector<T> const& v )
{
}

template class Vector<double>;
template class Vector<std::complex<double>>;
//template class Vector<long double>;

}
