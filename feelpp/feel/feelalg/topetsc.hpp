/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 31 Jan 2016

 Copyright (C) 2016 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_TOPETSC_HPP
#define FEELPP_TOPETSC_HPP 1

#include <boost/smart_ptr/shared_ptr.hpp>
#include <feel/feelalg/backendpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>


namespace Feel {

/**
 * returns a BackendPetsc shared_ptr from a Backend.
 * if the cast fails then it is a null shared_ptr
 *
 * here is a sample code:
 * @code
 * auto b = backend(); // default backend is petsc type
 * if ( auto bp = toPETSc( b ) ) // bp is BackendPetsc
 * {
 *   // use BackendPetsc interface
 * }
 * else // bp is a null shared_ptr and evaluates to false in the if above
 * {
 *   // it was not a BackendPetsc
 * }
 *
 * @endcode
 */
template<typename T, typename SizeT = uint32_type>
inline std::shared_ptr<BackendPetsc<T,SizeT>>
toPETSc( std::shared_ptr<Backend<T,SizeT>> const& b )
{
    return std::dynamic_pointer_cast<BackendPetsc<T,SizeT>>( b );
}

/**
 * returns a VectorPetsc shared_ptr from a Vector if it was a VectorPetsc
 * if the cast fails then it is a null shared_ptr
 *
 * here is a sample code:
 * @code
 * auto b = backend(); // default backend is petsc type
 * auto v = b->newVector();
 * if ( auto vp = toPETSc( v ) ) // bp is VectorPetsc
 * {
 *   // use VectorPetsc interface
 * }
 * else // bp is a null shared_ptr and evaluates to false in the if above
 * {
 *   // it was not a VectorPetsc
 * }
 *
 * @endcode
 */
template<typename T>
inline std::shared_ptr<VectorPetsc<T>>
toPETSc( std::shared_ptr<Vector<T>> const& v )
{
    return std::dynamic_pointer_cast<VectorPetsc<T>>( v );
}


/**
 * returns a MatrixPetsc shared_ptr from a Matrix if it was a MatrixPetsc
 * if the cast fails then it is a null shared_ptr
 *
 * here is a sample code:
 * @code
 * auto b = backend(); // default backend is petsc type
 * auto v = b->newMatrix();
 * if ( auto vp = toPETSc( v ) ) // bp is MatrixPetsc
 * {
 *   // use MatrixPetsc interface
 * }
 * else // bp is a null shared_ptr and evaluates to false in the if above
 * {
 *   // it was not a MatrixPetsc
 * }
 *
 * @endcode
 */
template<typename T>
inline std::shared_ptr<MatrixPetsc<T>>
toPETSc( std::shared_ptr<MatrixSparse<T>> const& m )
{
    return std::dynamic_pointer_cast<MatrixPetsc<T>>( m );
}

/**
 * returns a PreconditionerPetsc shared_ptr from a Preconditioner if it was a PreconditionerPetsc
 * if the cast fails then it is a null shared_ptr
 *
 * here is a sample code:
 * @code
 * auto b = backend(); // default backend is petsc type
  * if ( auto vp = toPETSc( b->preconditioner() ) ) // vp is PreconditionerPetsc
 * {
 *   // use PreconditionerPetsc interface
 * }
 * else // bp is a null shared_ptr and evaluates to false in the if above
 * {
 *   // it was not a PreconditionerPetsc
 * }
 *
 * @endcode
 */
template<typename T>
inline std::shared_ptr<PreconditionerPetsc<T>>
toPETSc( std::shared_ptr<Preconditioner<T>> const& p )
{
    return std::dynamic_pointer_cast<PreconditionerPetsc<T>>( p );
}


namespace detail
{

/**
 * returns a pair of  (c++ pointer VectorPetsc, shared_ptr VectorPetsc )
 * from input vector vec (with ref).
 * - vec is a VectorPetsc, return the pointer of &vec and null shared_ptr
 * - vec is a VectorUblas , return the pointer of petsc view and view in shared_ptr
 * - otherwise return nullptr and null shared_ptr
 * the shared_ptr is returned in order to not destroy the temporary object (ublas cases)
 */
template<typename T>
std::pair<VectorPetsc<T> *, std::shared_ptr<VectorPetsc<T> > >
toPETScPairPtr( Vector<T> & vec );

/**
 * returns a pair of  (c++ pointer VectorPetsc, shared_ptr VectorPetsc )
 * from input vector vec (with const ref).
 * - vec is a VectorPetsc, return the pointer of &vec and null shared_ptr
 * - vec is a ublas vector, return the pointer of petsc view and view in shared_ptr
 * - vec is another vector type and allowCopy is true, return a new VectorPetsc with values copied
 * - otherwise return nullptr and null shared_ptr
 * the shared_ptr is returned in order to not destroy the temporary object (ublas or copy cases)
 */
template<typename T>
std::pair<const VectorPetsc<T> *, std::shared_ptr<VectorPetsc<T> > >
toPETScPairPtr( Vector<T> const& vec, bool allowCopy = false );

} // namespace detail


}

#endif
