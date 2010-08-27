/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2007-08-14

   Copyright (C) 2007, 2009 Universit� Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-08-14
*/
#include <feel/feelalg/vectorepetra.hpp>

namespace Feel
{
#if defined(HAVE_TRILINOS_EPETRA)
template<typename T>
void
VectorEpetra<T>::init ( Epetra_BlockMap const& emap, const bool fast )
{
    // Clear initialized vectors
    if (this->isInitialized())
        this->clear();


	//create an MPI-enabled vector
#ifdef HAVE_MPI
    {
        _M_vec = Epetra_FEVector(emap);
    }
	// otherwise create a sequential vector if on only 1 processor
#else
    {
        _M_vec = Epetra_SerialDenseVector();
    }
#endif

    this->M_is_initialized = true;

    if (fast == false)
        this->zero ();
}

template<typename T>
void
VectorEpetra<T>::init(const size_type N, const size_type n_local, const bool fast)
{
    //int ierr=0;
    //int epetra_n=static_cast<int>(n);
    int epetra_n_local=static_cast<int>(n_local);

    // Clear initialized vectors
    if (this->isInitialized())
        this->clear();

    FEEL_ASSERT(n_local < N)( n_local )( N ).error( "invalid local size" );

    /*
    _M_emap = Epetra_Map( -1, epetra_n_local, 0, Epetra_MpiComm(M_comm) );
    _M_vec = Epetra_FEVector(_M_emap);
    */

    this->M_is_initialized = true;

    if (fast == false)
        this->zero ();
}

template<typename T>   //actually all should be double
void
VectorEpetra<T>::set ( const value_type& value)
{
    value_type epetra_value = static_cast<value_type>(value);

    _M_vec.PutScalar(epetra_value);
}

template<typename T>
void
VectorEpetra<T>::set ( size_type i, const value_type& value)
{
    FEEL_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int i_val = static_cast<int>(i);
    value_type epetra_value = static_cast<value_type>(value);

    //_M_vec[i_val] = epetra_value;
    _M_vec.ReplaceGlobalValues(1, &i_val, &epetra_value);

    Debug(10010) << "[Vector] Replacing value in row " << i << " for " << value << "\n";

}

template<typename T>
void
VectorEpetra<T>::add (const size_type i, const value_type& value)
{
    FEEL_ASSERT(i<size())( i )( size() ).error( "invalid index" );

    int i_val = static_cast<int>(i);
    value_type epetra_value = static_cast<value_type>(value);
    int ierr;
    //ierr= _M_vec.SumIntoGlobalValues(1,&epetra_value,&i_val);
    ierr= _M_vec.SumIntoGlobalValues(1,&i_val, &epetra_value);//indices are in global index space

    if (ierr != 0)
        {
            Debug() << "ERRORCODE SumIntoGlobalValues VECTOR: " << ierr <<  " in V(" << i_val << ") for value "<< epetra_value << "." << "\n";
        }
}
template<typename T>
void
VectorEpetra<T>::addVector ( int* i, int n, value_type* v )
{
    FEEL_ASSERT(n<size())( n )( size() ).error( "invalid index array size" );

    int ierr;
    //indices are in global index space
    ierr= _M_vec.SumIntoGlobalValues(n,i,v);

    if (ierr != 0)
    {
        Debug() << "ERRORCODE SumIntoGlobalValues VECTOR: " << ierr <<  " in V \n";
    }
}

template<typename T>
void
VectorEpetra<T>::printMatlab (const std::string name) const
{
    FEEL_ASSERT (this->closed()).warn("epetra vector not closed");
#if 0
    if ( !this->closed() )
        const_cast<VectorEpetra<T>*>(this)->close();
#endif

    Debug() << "[printMatlab] print vector in matlab file " << name << "\n";
    EpetraExt::MultiVectorToMatlabFile( name.c_str(), _M_vec);
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

