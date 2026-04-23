/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Chabannes Vincent <vincent.chabannes@feelpp.org>
       Date: 2016-03-29

  Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

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
   \file topetsc.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-03-29
 */

#include <feel/feelalg/topetsc.hpp>

namespace Feel
{
namespace detail
{
template<typename T>
std::shared_ptr<VectorPetsc<T>>
toOwnedPETScCopy( Vector<T> const& vec )
{
    std::shared_ptr<VectorPetsc<T>> vec_petscClone;
    if ( vec.comm().size() > 1 )
        vec_petscClone = std::make_shared<VectorPetscMPI<T>>( vec.mapPtr() );
    else
        vec_petscClone = std::make_shared<VectorPetsc<T>>( vec.mapPtr() );

    if ( auto const* vec_ublas = dynamic_cast<VectorUblas<T> const*>( &vec ) )
    {
        int ierr = 0;
        if ( vec.comm().size() > 1 )
        {
            Vec lx = nullptr;
            PetscScalar* valuesOut = nullptr;
            ierr = VecGhostGetLocalForm( vec_petscClone->vec(), &lx );
            CHKERRABORT( vec.comm(), ierr );
            ierr = VecGetArray( lx, &valuesOut );
            CHKERRABORT( vec.comm(), ierr );
            for ( typename Vector<T>::size_type k = 0; k < vec_ublas->localSize(); ++k )
                valuesOut[k] = static_cast<PetscScalar>( (*vec_ublas)( k ) );
            ierr = VecRestoreArray( lx, &valuesOut );
            CHKERRABORT( vec.comm(), ierr );
            ierr = VecGhostRestoreLocalForm( vec_petscClone->vec(), &lx );
            CHKERRABORT( vec.comm(), ierr );
        }
        else
        {
            PetscScalar* valuesOut = nullptr;
            ierr = VecGetArray( vec_petscClone->vec(), &valuesOut );
            CHKERRABORT( vec.comm(), ierr );
            for ( typename Vector<T>::size_type k = 0; k < vec_ublas->localSize(); ++k )
                valuesOut[k] = static_cast<PetscScalar>( (*vec_ublas)( k ) );
            ierr = VecRestoreArray( vec_petscClone->vec(), &valuesOut );
            CHKERRABORT( vec.comm(), ierr );
        }
        return vec_petscClone;
    }

    *vec_petscClone = vec;
    return vec_petscClone;
}

template<typename T>
std::pair<VectorPetsc<T> *, std::shared_ptr<VectorPetsc<T> > >
toPETScPairPtr( Vector<T> & vec )
{
    if ( !vec.closed() )
        const_cast<Vector<T>*>( &vec )->close();

    // vector need if the cases ublas or vector copy are used
    std::shared_ptr<VectorPetsc<T>> vec_petscClone;

    // petsc vector
    VectorPetsc<T> * vec_petsc = dynamic_cast<VectorPetsc<T> *>( &vec );
    if ( vec_petsc )
        return std::make_pair(vec_petsc, vec_petscClone);

    VectorPetsc<T> * vec_petscUsed = nullptr;
    // ublas vector
    const VectorUblas<T> * vec_ublas = dynamic_cast<const VectorUblas<T> *>( &vec );
    if( vec_ublas )
    {
        vec_petscClone = toPETScPtr( *vec_ublas );
        vec_petscUsed = &(*vec_petscClone);
        return std::make_pair( vec_petscUsed, vec_petscClone );
    }

    return std::make_pair( vec_petscUsed,vec_petscClone );
}

template<typename T>
std::pair<const VectorPetsc<T> *, std::shared_ptr<VectorPetsc<T> > >
toPETScPairPtr( Vector<T> const& vec, bool allowCopy )
{
    if ( !vec.closed() )
        const_cast<Vector<T>*>( &vec )->close();

    // vector need if the cases ublas or vector copy are used
    std::shared_ptr<VectorPetsc<T>> vec_petscClone;

    // petsc vector
    VectorPetsc<T> const* vec_petsc = dynamic_cast<VectorPetsc<T> const*>( &vec );
    if ( vec_petsc )
    {
        return std::make_pair(vec_petsc, vec_petscClone);
    }

    const VectorPetsc<T> * vec_petscUsed = nullptr;
    // ublas vector
    //typedef VectorUblas<T> vector_ublas_type;
    //typedef typename vector_ublas_type::range::type vector_ublas_range_type;
    //typedef typename vector_ublas_type::shallow_array_adaptor::type vector_ublas_extarray_type;
    //typedef typename vector_ublas_extarray_type::range::type vector_ublas_extarray_range_type;
    //const vector_ublas_type * vec_ublas = dynamic_cast<vector_ublas_type const*>( &vec );
    //const vector_ublas_range_type * vec_ublasRange = dynamic_cast<vector_ublas_range_type const*>( &vec );
    //const vector_ublas_extarray_type * vec_ublasExtArray = dynamic_cast<vector_ublas_extarray_type const*>( &vec );
    //const vector_ublas_extarray_range_type * vec_ublasExtArrayRange = dynamic_cast<vector_ublas_extarray_range_type const*>( &vec );
    //bool vecIsUblasVarients = vec_ublas || vec_ublasRange || vec_ublasExtArray || vec_ublasExtArrayRange;
    //if ( vecIsUblasVarients )
    //{
        //if ( vec_ublas )
            //vec_petscClone = toPETScPtr( *(const_cast<vector_ublas_type *>( vec_ublas ) ) );
        //else if ( vec_ublasRange )
            //vec_petscClone = toPETScPtr( *(const_cast<vector_ublas_range_type *>( vec_ublasRange ) ) );
        //else if ( vec_ublasExtArray )
            //vec_petscClone = toPETScPtr( *(const_cast<vector_ublas_extarray_type *>( vec_ublasExtArray ) ) );
        //else if ( vec_ublasExtArrayRange )
            //vec_petscClone = toPETScPtr( *(const_cast<vector_ublas_extarray_range_type *>( vec_ublasExtArrayRange ) ) );
        //vec_petscUsed = &(*vec_petscClone);
        //return std::make_pair( vec_petscUsed,vec_petscClone );
    //}
    const VectorUblas<T> * vec_ublas = dynamic_cast<const VectorUblas<T> *>( &vec );
    if( vec_ublas )
    {
        vec_petscClone = toOwnedPETScCopy( vec );
        vec_petscUsed = &(*vec_petscClone);
        return std::make_pair( vec_petscUsed, vec_petscClone );
    }
    // create a new vector and copy values
    if ( allowCopy )
    {
        vec_petscClone = toOwnedPETScCopy( vec );
        vec_petscUsed = &(*vec_petscClone);
        return std::make_pair( vec_petscUsed,vec_petscClone );
    }

    return std::make_pair( vec_petscUsed,vec_petscClone );
}

// instantiation
template
std::pair<VectorPetsc<double> *, std::shared_ptr<VectorPetsc<double> > >
toPETScPairPtr( Vector<double> & vec );
template
std::pair<const VectorPetsc<double> *, std::shared_ptr<VectorPetsc<double> > >
toPETScPairPtr( Vector<double> const& vec, bool allowCopy );
template
std::shared_ptr<VectorPetsc<double>>
toOwnedPETScCopy( Vector<double> const& vec );

} // namespace detail

}// namespace Feel
