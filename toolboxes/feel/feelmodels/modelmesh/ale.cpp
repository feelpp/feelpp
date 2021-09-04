/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Copyright (C) 2010 University of Coimbra

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
 \file ale.cpp
 \author Goncalo Pena <gpena@mat.uc.pt>
 \date 2010-10-12
 */

//#include <boost/preprocessor/comparison/greater_equal.hpp>

#include <feel/feelmodels/modelmesh/ale.hpp>
#include <feel/feelmodels/modelmesh/ale_impl.hpp>


//#include <feel/feelfilters/gmsh.hpp>


namespace Feel
{
namespace FeelModels
{

template < class Convex, int Order >
ALE<Convex,Order>::ALE( /*mesh_ptrtype mesh,*/ std::string prefix, worldcomm_ptr_t const& worldcomm,
                        ModelBaseRepository const& modelRep )
    :
    super_type( prefix/*prefixvm(prefix,"alemesh")*/, worldcomm,"",modelRep )
    //super_type( mesh,prefix,worldcomm,moveGhostEltFromExtendedStencil ),
{
    M_flagSet["fixed"].clear();
    M_flagSet["moving"].clear();
    M_flagSet["free"].clear();
}

template< class Convex, int Order >
typename ALE<Convex,Order>::self_ptrtype
ALE<Convex,Order>::build( mesh_ptrtype mesh, std::string prefix,
                          worldcomm_ptr_t const& worldcomm,
                          ModelBaseRepository const& modelRep)
{
    return self_ptrtype( new ALE_IMPL::ALE<Convex,Order>(mesh,prefix,worldcomm,modelRep ) );
}

template< class Convex, int Order >
typename ALE<Convex,Order>::self_ptrtype
ALE<Convex,Order>::build( mesh_ptrtype mesh, range_elements_type const& rangeElt, std::string prefix,
                          worldcomm_ptr_t const& worldcomm,
                          ModelBaseRepository const& modelRep)
{
    return self_ptrtype( new ALE_IMPL::ALE<Convex,Order>(mesh,rangeElt,prefix,worldcomm,modelRep ) );
}



template< class Convex, int Order >
typename ALE<Convex,Order>::flagSet_type const&
ALE<Convex,Order>::flagSet() const { return M_flagSet; }

template < class Convex, int Order >
void
ALE<Convex,Order>::addBoundaryFlags( std::string str, flag_type flag )
{
    if ( str == "fixed" )
        M_flagSet["fixed"].push_back(flag);
    else if ( str == "moving" )
        M_flagSet["moving"].push_back(flag);
    else if ( str == "free" )
        M_flagSet["free"].push_back(flag);
    else
        CHECK( false ) << "invalid flag type" << str << " with flag name " << flag;
}
template < class Convex, int Order >
void
ALE<Convex,Order>::addBoundaryFlags( flagSet_type flags )
{
    M_flagSet = flags;
}
template < class Convex, int Order >
void
ALE<Convex,Order>::clearFlagSet()
{
    //M_flagSet.clear();
    M_flagSet["fixed"].clear();
    M_flagSet["moving"].clear();
    M_flagSet["free"].clear();
}
template < class Convex, int Order >
std::vector<flag_type> const&
ALE<Convex,Order>::flagSet(std::string key) const
{
    CHECK( M_flagSet.find(key) != M_flagSet.end() ) << "the flag type " << key << " is unknown \n";
    return M_flagSet.find(key)->second;
}
template < class Convex, int Order >
flag_type
ALE<Convex,Order>::flagSet(std::string key, int k) const
{
    CHECK( M_flagSet.find(key) != M_flagSet.end() ) << "the flag type " << key << " is unknown \n";
    CHECK( M_flagSet.find(key)->second.size() > k ) << "the key " << k << " must be <  " <<  M_flagSet.find(key)->second.size() << "\n";
    return M_flagSet.find(key)->second.at(k);
}
template < class Convex, int Order >
void
ALE<Convex,Order>::setFlagSet( flagSet_type const & fl )
{
    M_flagSet=fl;
}

} // namespace FeelModels
} // namespace Feel

