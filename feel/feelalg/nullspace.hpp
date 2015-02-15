/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Chabannes Vincent <vincent.chabannes@feelpp.org>
       Date: 2015-02-13

  Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)

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
   \file nullspace.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-02-13
 */

#ifndef NullSpace_H
#define NullSpace_H 1

#include <feel/feeldiscr/functionspacebase.hpp>

namespace Feel
{

template <typename T=double>
class NullSpace
{
public :
    typedef T value_type;
    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    NullSpace() {}

    template <typename EltType>
    NullSpace( std::vector<EltType> const& vecBasis )
        :
        M_basisVector( vecBasis.size() )
    {
        this->updateForUse( vecBasis );
    }

    template <typename EltType>
    NullSpace( std::initializer_list<EltType> const& listBasis )
        :
        M_basisVector( listBasis.size() )
    {
        this->updateForUse( listBasis );
    }

    NullSpace( backend_ptrtype const& mybackend, NullSpace const& ns, bool redoOrthonormalization = false )
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

    NullSpace( NullSpace const& ns )
        :
        M_basisVector( ns.M_basisVector )
        {}


    int size() const { return M_basisVector.size(); }
    vector_ptrtype const& basisVectorPtr(int k) const { return M_basisVector[k]; }
    vector_type const& basisVector(int k) const { return *M_basisVector[k]; }

private :

    template <typename EltType>
    void updateForUse( std::vector<EltType> const& vecBasis )
    {
        this->updateForUse( vecBasis, mpl::bool_<boost::is_base_of<FunctionSpaceBase::ElementBase,EltType>::value >() );
    }
    template <typename EltType>
    void updateForUse( std::vector<EltType> const& vecBasis, mpl::true_ )
    {
        // copy vectors
        for ( int k=0;k<vecBasis.size();++k )
        {
            M_basisVector[k] = vecBasis[k].functionSpace()->elementPtr();
            *M_basisVector[k] = vecBasis[k];
        }
        // orthonormalise basis if necessary
        if ( !this->hasOrthonormalBasisVectors() )
            this->orthonormalizeBasisVector();
    }
    template <typename EltType>
    void updateForUse( std::vector<EltType> const& vec, mpl::false_ )
    {
        CHECK( false ) << "case not implemented : vector not define with a function space";
    }

    template <typename EltType>
    void updateForUse( std::initializer_list<EltType> const& listBasis )
    {
        this->updateForUse( listBasis, mpl::bool_<boost::is_base_of<FunctionSpaceBase::ElementBase,EltType>::value >() );
    }
    template <typename EltType>
    void updateForUse( std::initializer_list<EltType> const& listBasis, mpl::true_ )
    {
        M_basisVector.resize( listBasis.size() );
        // copy vectors
        int k=0;
        for (  auto const& vec : listBasis )
        {
            M_basisVector[k] = vec.functionSpace()->elementPtr();
            *M_basisVector[k] = vec;
            ++k;
        }
        // orthonormalise basis if necessary
        if ( !this->hasOrthonormalBasisVectors() )
            this->orthonormalizeBasisVector();
    }
    template <typename EltType>
    void updateForUse( std::initializer_list<EltType> const& listVec, mpl::false_ )
    {
        CHECK( false ) << "case not implemented : vector not define with a function space";
    }

    void orthonormalizeBasisVector()
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

    bool hasOrthonormalBasisVectors()
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


private :
    std::vector<vector_ptrtype> M_basisVector;
};

} // namespace Feel

#endif /* NullSpace_H */
