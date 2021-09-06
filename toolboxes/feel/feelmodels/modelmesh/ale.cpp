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


#include <feel/feelmodels/modelmesh/ale.hpp>
#include <feel/feelmodels/modelmesh/ale_impl.hpp>

namespace Feel
{
namespace FeelModels
{

template < class Convex, int Order >
ALE<Convex,Order>::ALE( std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                        ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldcomm,"",modelRep )
{
    M_bcToMarkers["fixed"].clear();
    M_bcToMarkers["moving"].clear();
    M_bcToMarkers["free"].clear();
}

template< class Convex, int Order >
typename ALE<Convex,Order>::self_ptrtype
ALE<Convex,Order>::build( mesh_ptrtype mesh, std::string const& prefix,
                          worldcomm_ptr_t const& worldcomm,
                          ModelBaseRepository const& modelRep)
{
    return self_ptrtype( new ALE_IMPL::ALE<Convex,Order>(mesh,prefix,worldcomm,modelRep ) );
}

template< class Convex, int Order >
typename ALE<Convex,Order>::self_ptrtype
ALE<Convex,Order>::build( mesh_ptrtype mesh, range_elements_type const& rangeElt, std::string const& prefix,
                          worldcomm_ptr_t const& worldcomm,
                          ModelBaseRepository const& modelRep)
{
    return self_ptrtype( new ALE_IMPL::ALE<Convex,Order>(mesh,rangeElt,prefix,worldcomm,modelRep ) );
}


template < class Convex, int Order >
void
ALE<Convex,Order>::addMarkerInBoundaryCondition( std::string const& bc, std::string const& marker )
{
    CHECK( bc == "fixed" || bc == "moving" || bc == "free" ) << "invalid bc : " << bc;
    M_bcToMarkers[bc].insert( marker );
}

#if 0
template < class Convex, int Order >
void
ALE<Convex,Order>::clearFlagSet()
{
    //M_flagSet.clear();
    M_flagSet["fixed"].clear();
    M_flagSet["moving"].clear();
    M_flagSet["free"].clear();
}
#endif
} // namespace FeelModels
} // namespace Feel

