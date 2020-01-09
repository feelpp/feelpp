/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#ifndef FEELPP_NULLSPACE_HPP
#define FEELPP_NULLSPACE_HPP 1

#include <feel/feeldiscr/functionspacebase.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel
{
template<typename T,typename SizeT> class Backend;

template <typename T=double, typename SizeT = uint32_type>
class FEELPP_EXPORT NullSpace
{
public :
    typedef T value_type;
    using size_type = SizeT;
    typedef Vector<value_type,size_type> vector_type;
    typedef std::shared_ptr<vector_type> vector_ptrtype;
    typedef Backend<value_type,size_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

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

    NullSpace( std::vector<vector_ptrtype> const& vecBasis, backend_ptrtype const& mybackend );

    NullSpace( backend_ptrtype const& mybackend, std::shared_ptr<NullSpace> const& ns, bool redoOrthonormalization = false )
        :
        NullSpace( mybackend, *ns, redoOrthonormalization )
        {
        }
    NullSpace( backend_ptrtype const& mybackend, NullSpace const& ns, bool redoOrthonormalization = false );

    NullSpace( NullSpace const& ns );

    int size() const { return M_basisVector.size(); }
    vector_ptrtype const& basisVectorPtr(int k) const { return M_basisVector[k]; }
    vector_type const& basisVector(int k) const { return *M_basisVector[k]; }
    void push_back( vector_ptrtype v ) { M_basisVector.push_back( v ); }
    void close();
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

    void updateForUse( std::vector<vector_ptrtype> const& vecBasis, backend_ptrtype const& mybackend );

    void orthonormalizeBasisVector();

    bool hasOrthonormalBasisVectors();


private :
    std::vector<vector_ptrtype> M_basisVector;
};

template<typename SpaceType, typename ExprType>
void nullspace( NullSpace<double,typename SpaceType::size_type>& K, std::shared_ptr<SpaceType> Xh, ExprType&& e )
{
    auto u = Xh->elementPtr(e);
    K.push_back( u );
    K.close();
}

template<typename SpaceType, typename ExprType, typename... Args>
void nullspace( NullSpace<double,typename SpaceType::size_type>& K, std::shared_ptr<SpaceType> Xh, ExprType&& e, Args... args )
{
    auto u = Xh->elementPtr(e);
    K.push_back( u );
    nullspace( K, Xh, args...);
}

template<typename SpaceType, typename ExprType, typename... Args>
NullSpace<double,typename SpaceType::size_type> nullspace( std::shared_ptr<SpaceType> Xh, ExprType&& e,  Args... args )
{
    NullSpace<double,typename SpaceType::size_type> K;
    nullspace( K, Xh, e, args...);
    return K;
}

template<typename SpaceType, typename ExprType, typename... Args>
std::shared_ptr<NullSpace<double,typename SpaceType::size_type>>  nullspace_ptr( std::shared_ptr<SpaceType> Xh, ExprType&& e, Args... args )
{
    std::shared_ptr<NullSpace<double,typename SpaceType::size_type>> K( std::make_shared<NullSpace<double,typename SpaceType::size_type>>() );
    nullspace( *K, Xh, e, args...);
    return K;
}

template<typename SizeT>
inline std::shared_ptr<NullSpace<double,SizeT>>
toBackend( std::shared_ptr<Backend<double,SizeT>> const& b, std::shared_ptr<NullSpace<double,SizeT>> const& Kspace, bool redoOrtho = false )
{
    return std::make_shared<NullSpace<double,SizeT>>( b, Kspace, redoOrtho );
}

} // namespace Feel

#endif /* NullSpace_H */
