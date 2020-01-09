/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Chabannes Vincent <vincent.chabannes@feelpp.org>
       Date: 2015-03-17

  Copyright (C) 2007-2015 Université Joseph Fourier (Grenoble I)

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
   \file nullspace.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-03-17
 */

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/nullspace.hpp>

namespace Feel
{

template <typename T, typename SizeT>
NullSpace<T,SizeT>::NullSpace( std::vector<vector_ptrtype> const& vecBasis, backend_ptrtype const& mybackend )
    :
    M_basisVector( vecBasis.size() )
{
    this->updateForUse( vecBasis, mybackend );
}

//template<typename T> class Backend;

template <typename T, typename SizeT>
NullSpace<T,SizeT>::NullSpace( backend_ptrtype const& mybackend, NullSpace const& ns, bool redoOrthonormalization )
    :
    M_basisVector( ns.size() )
{
    // copy vectors
    for ( int k=0;k<ns.size();++k )
    {
        M_basisVector[k] = mybackend->newVector( ns.basisVector(k).mapPtr() );
        *M_basisVector[k] = ns.basisVector(k);
        M_basisVector[k]->close();
    }
    // we assume that nullspace has an orthonormal basis, so do nothing if not forced
    if ( redoOrthonormalization )
        this->orthonormalizeBasisVector();
}

template <typename T, typename SizeT>
NullSpace<T,SizeT>::NullSpace( NullSpace const& ns )
    :
    M_basisVector( ns.M_basisVector )
{}

template <typename T, typename SizeT>
void
NullSpace<T,SizeT>::close()
{
    if ( M_basisVector.size() )
    {
        LOG(INFO) << "Close NullSpace";
        // orthonormalise basis if necessary
        if ( !this->hasOrthonormalBasisVectors() )
        {
            LOG(INFO) << "NullSpace: orthonormalize basis functions";
            this->orthonormalizeBasisVector();
        }
    }
}


template <typename T, typename SizeT>
void
NullSpace<T,SizeT>::updateForUse( std::vector<vector_ptrtype> const& vecBasis, backend_ptrtype const& mybackend )
{
    // copy vectors
    for ( int k=0;k<vecBasis.size();++k )
    {
        M_basisVector[k] = mybackend->newVector( vecBasis[k]->mapPtr() );
        *M_basisVector[k] = *vecBasis[k];
        M_basisVector[k]->close();
    }
    // we assume that nullspace has an orthonormal basis, so do nothing if not forced
    // orthonormalise basis if necessary
    if ( !this->hasOrthonormalBasisVectors() )
        this->orthonormalizeBasisVector();
}

template <typename T, typename SizeT>
void
NullSpace<T,SizeT>::orthonormalizeBasisVector()
{
    int dimNullSpace = this->size();
    std::vector<double> dots(dimNullSpace);
    for ( int k = 0 ; k<dimNullSpace ; ++k )
    {
        auto const& myvec = this->basisVector(k);
        for ( int k2 = 0 ; k2 < k ; ++k2 )
        {
            auto const& myvec2 = this->basisVector(k2);
            dots[k2] = dot(myvec,myvec2);
        }
        for (int k2=0; k2<k; k2++)
            dots[k2] *= -1.;

        for ( int k2 = 0 ; k2<k ; ++k2 )
            this->basisVectorPtr(k)->add( dots[k2],this->basisVector(k2) );

        // normalize
        double normL2 = this->basisVector(k).l2Norm();
        if ( normL2 > 1e-9 )
            this->basisVectorPtr(k)->scale(1.0/normL2 );
    }
    CHECK( this->hasOrthonormalBasisVectors() ) << "basis vectors are not orthonormal\n";
}

template <typename T, typename SizeT>
bool
NullSpace<T,SizeT>::hasOrthonormalBasisVectors()
{
    int dimNullSpace = this->size();
    for ( int k = 0 ; k<dimNullSpace ; ++k )
    {
        auto const& myvec = this->basisVector(k);
        for ( int k2 = 0 ; k2<dimNullSpace ; ++k2 )
        {
            auto const& myvec2 = this->basisVector(k2);
            if ( k == k2 )
            {
                double dotEval = myvec.l2Norm();
                if ( math::abs(dotEval-1.0) > 1e-11 ) return false;
            }
            else
            {
                double dotEval = dot(myvec,myvec2);
                if ( math::abs(dotEval) > 1e-11 ) return false;
            }
        }
    }
    return true;
}


template class NullSpace<double, uint32_type>;

} // namespace Feel
