/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2007-08-14

   Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file vectorepetra.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-08-14
*/
#include <feel/feelalg/vectorepetra.hpp>
#include <feel/feelalg/matrixepetra.hpp>

namespace Feel
{
#if defined(FEELPP_HAS_TRILINOS_EPETRA)

template<typename T>
VectorEpetra<T>::VectorEpetra( )
    :
    super(),
#ifdef FEELPP_HAS_MPI
    M_emap( Epetra_BlockMap( -1, 0, 0, Epetra_MpiComm( super::comm() ) ) ),
    M_vec( M_emap ) // false (zerout)?
#else
    M_emap( Epetra_BlockMap( -1, 0, 0, Epetra_SerialComm ) ),
    M_vec( M_emap )
#endif
{
}

template<typename T>
VectorEpetra<T>::VectorEpetra ( Epetra_BlockMap const& emap )
    :
    super( emap.NumGlobalElements(), emap.NumMyElements() ),
    M_emap( emap ),
    M_vec( M_emap )
    //M_destroy_vec_on_exit( true )
{
    this->init( M_emap, true );
    checkInvariants();
}

template<typename T>
VectorEpetra<T>::VectorEpetra ( Epetra_FEVector const * v )
    :
    super( v->Map().NumGlobalElements(), v->Map().NumMyElements() ),
    M_emap( v->Map() ),
    M_vec( *v )
{
    this->init( M_emap, true );
    checkInvariants();
}

template<typename T>
VectorEpetra<T>::VectorEpetra ( Epetra_Vector const * v )
    :
    super( v->Map().NumGlobalElements(), v->Map().NumMyElements() ),
    M_emap( v->Map() ),
    M_vec ( M_emap )
{
    this->init( M_emap, true );
    //double** V;
    //v->ExtractView(&V);
    //printf("first val : %f\n",V[0][0]);


    M_vec.Update( 1.0,*v,1.0 );

    //for( size_type i = 0; i < this->localSize(); ++i )
    //    {
    //        this->set( i,  V[ 0 ][i] );
    //    }
    checkInvariants();
}
template<typename T>
VectorEpetra<T>::VectorEpetra ( VectorEpetra const& v )
    :
    super( v ),
    M_emap( v.Map() ),
    M_vec( v.vec() )

{
    this->M_is_initialized = true;
    checkInvariants();
}

template<typename T>
VectorEpetra<T>::~VectorEpetra ()
{

}

template<typename T>
void
VectorEpetra<T>::init ( Epetra_BlockMap const& emap, const bool fast )
{
    // Clear initialized vectors
    if ( this->isInitialized() )
        this->clear();


    //create an MPI-enabled vector
#ifdef FEELPP_HAS_MPI
    {
        M_vec = Epetra_FEVector( emap );
    }
    // otherwise create a sequential vector if on only 1 processor
#else
    {
        M_vec = Epetra_SerialDenseVector();
    }
#endif

    this->M_is_initialized = true;

    if ( fast == false )
        this->zero ();
}

template<typename T>
void
VectorEpetra<T>::init( const size_type N, const size_type n_local, const bool fast )
{
    //int ierr=0;
    //int epetra_n=static_cast<int>(n);
    int epetra_n_local=static_cast<int>( n_local );

    // Clear initialized vectors
    if ( this->isInitialized() )
        this->clear();

    FEELPP_ASSERT( n_local < N )( n_local )( N ).error( "invalid local size" );

    /*
    M_emap = Epetra_Map( -1, epetra_n_local, 0, Epetra_MpiComm(M_comm) );
    M_vec = Epetra_FEVector(M_emap);
    */

    this->M_is_initialized = true;

    if ( fast == false )
        this->zero ();
}

template<typename T>   //actually all should be double
void
VectorEpetra<T>::set ( const value_type& value )
{
    value_type epetra_value = static_cast<value_type>( value );

    M_vec.PutScalar( epetra_value );
}

template<typename T>
void
VectorEpetra<T>::set ( size_type i, const value_type& value )
{
    FEELPP_ASSERT( i<size() )( i )( size() ).error( "invalid index" );

    int i_val = static_cast<int>( i );
    value_type epetra_value = static_cast<value_type>( value );

    //M_vec[i_val] = epetra_value;
    M_vec.ReplaceGlobalValues( 1, &i_val, &epetra_value );

    DVLOG(2) << "[Vector] Replacing value in row " << i << " for " << value << "\n";

}

template<typename T>
void
VectorEpetra<T>::add ( const size_type i, const value_type& value )
{
    FEELPP_ASSERT( i<size() )( i )( size() ).error( "invalid index" );

    int i_val = static_cast<int>( i );
    value_type epetra_value = static_cast<value_type>( value );
    int ierr;
    //ierr= M_vec.SumIntoGlobalValues(1,&epetra_value,&i_val);
    ierr= M_vec.SumIntoGlobalValues( 1,&i_val, &epetra_value ); //indices are in global index space

    if ( ierr != 0 )
    {
        VLOG(1) << "ERRORCODE SumIntoGlobalValues VECTOR: " << ierr <<  " in V(" << i_val << ") for value "<< epetra_value << "." << "\n";
    }
}
template<typename T>
void
VectorEpetra<T>::addVector ( int* i, int n, value_type* v )
{
    //FEELPP_ASSERT(n<=size())( n )( size() ).error( "invalid index array size" );

    int ierr;
    //indices are in global index space
    ierr= M_vec.SumIntoGlobalValues( n,i,v );

    if ( ierr != 0 )
    {
        VLOG(1) << "ERRORCODE SumIntoGlobalValues VECTOR: " << ierr <<  " in V \n";
    }
}
template<typename T>
void
VectorEpetra<T>::checkInvariants () const
{
    FEELPP_ASSERT( this->localSize() == M_emap.NumMyElements() )
    ( this->localSize() )
    ( M_emap.NumMyElements() ).error( "Invalid VectorEpetra (Local size)" );
    FEELPP_ASSERT( this->size() == M_emap.NumGlobalElements() )
    ( this->size() )
    ( M_emap.NumGlobalElements() ).error( "Invalid VectorEpetra (Global size)" );

    FEELPP_ASSERT( this->firstLocalIndex() == M_emap.MinLID() )
    ( this->firstLocalIndex() )
    ( M_emap.MinLID() ).error( "Invalid VectorEpetra (Min global index)" );

    FEELPP_ASSERT( this->lastLocalIndex() == M_emap.MaxLID()+1 )
    ( this->lastLocalIndex() )
    ( M_emap.MaxLID() ).error( "Invalid VectorEpetra (Min global index)" );

    FEELPP_ASSERT( this->localSize() == M_vec.Map().NumMyElements() )
    ( this->localSize() )
    ( M_vec.Map().NumMyElements() ).error( "Invalid VectorEpetra (Local size)" );
    FEELPP_ASSERT( this->size() == M_vec.Map().NumGlobalElements() )
    ( this->size() )
    ( M_vec.Map().NumGlobalElements() ).error( "Invalid VectorEpetra (Global size)" );

    FEELPP_ASSERT( this->firstLocalIndex() == M_vec.Map().MinLID() )
    ( this->firstLocalIndex() )
    ( M_vec.Map().MinLID() ).error( "Invalid VectorEpetra (Min global index)" );

    FEELPP_ASSERT( this->lastLocalIndex() == M_vec.Map().MaxLID()+1 )
    ( this->lastLocalIndex() )
    ( M_vec.Map().MaxLID() ).error( "Invalid VectorEpetra (Min global index)" );

}
template<typename T>
void
VectorEpetra<T>::printMatlab ( const std::string name, bool renumber ) const
{
    FEELPP_ASSERT ( this->closed() ).warn( "epetra vector not closed" );
#if 0

    if ( !this->closed() )
        const_cast<VectorEpetra<T>*>( this )->close();

#endif

		std::string vectorName = "var_"+name;
    VLOG(1) << "[printMatlab] print vector in matlab file " << name << "\n";
    EpetraExt::MultiVectorToMatlabFile( name.c_str(), M_vec, vectorName.c_str() );
}
template<typename T>
void
VectorEpetra<T>::addVector ( const Vector<T>& _v,
                             const MatrixSparse<T>& _M )
{
    epetra_vector_type res( this->Map() );

    dynamic_cast<MatrixEpetra const&>( _M ).multiply( false, _v, res );
    this->add( 1., res );
}

//
// Debug stream
//
template<typename T>
DebugStream&
operator<<( DebugStream& __os, VectorEpetra<T> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.vec();

        __os << __str.str() << "\n";
    }

    return __os;
}

template<typename T>
NdebugStream&
operator<<( NdebugStream& __os, VectorEpetra<T> const& __n )
{
    return __os;
}

DebugStream&
operator<<( DebugStream& __os, Epetra_Map const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}

NdebugStream&
operator<<( NdebugStream& __os, Epetra_Map const& __n )
{
    return __os;
}


//
// Explicit instantiation
//
template class VectorEpetra<double>;

#endif // VectorEpetra
}
