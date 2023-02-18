/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 07 juil. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#pragma once

#include <type_traits>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * @brief create an isotropic metric field based on a levelset gradation
 * 
 * @tparam SpaceT function space type
 * @param Xh function space
 * @param facets ranges of facets to compute the levelset function based on the distance function to these facets
 * @param hclose mesh size close to the facets
 * @param hfar mesh size far away from the facets
 * @param cst real between 0 and 1 providing  where the gradation is constant
 * @return Pch_element_type<MeshT, 1> continuous piecewise P1 function representing the metric
 */
template <typename SpaceT, std::enable_if_t<std::is_base_of_v<FunctionSpaceBase, SpaceT>,SpaceT>* = nullptr>
typename SpaceT::element_type
gradedfromls( std::shared_ptr<SpaceT> const& Xh, faces_reference_wrapper_t<typename SpaceT::mesh_ptrtype> const& facets, double hclose, double hfar, double cst = 0. )
{
    using namespace vf;
    auto d = distanceToRange( _space=Xh, _range=facets );
    d *= 1./d.max(); // max of d is 1
    auto g = Xh->element();
    CHECK( hclose > 0 && hfar > 0 ) << fmt::format("hclose({}) and/or hfar({}) must be strictly positive", hclose, hfar);
    LOG_IF(WARNING, ( (cst < 0) || (cst > 1) ) ) << fmt::format("[gradedfromls] invalid constant percentage {} (should be > 0 and < 1)",cst);

    cst = (cst < 0)?0:(cst > 1)?1:cst;
    g.on(_range=elements(Xh->mesh()),_expr=(idv(d) > cst)*(hclose*(1-idv(d))+hfar*idv(d))+(idv(d)<=cst)*hclose );
    return g;
}

/**
 * @brief create an isotropic metric field based on a scalar expression
 * 
 * @tparam SpaceT function space type
 * @tparam ExprT expression type of the metric
 * @param Xh function space
 * @param e expression
 * @return Pch_element_type<MeshT, 1> continuous piecewise P1 function representing the metric
 */
template <typename SpaceT, typename ExprT,std::enable_if_t<std::is_base_of_v<FunctionSpaceBase, SpaceT>&& std::is_base_of_v<ExprBase,ExprT>,SpaceT>* = nullptr>
typename SpaceT::element_type
expr( std::shared_ptr<SpaceT> const& Xh, ExprT const& e )
{
    auto g = Xh->element();
    g.on( _range = elements( Xh->mesh() ), _expr = e );
    return g;
}

/**
 * @brief create an isotropic metric field based on a scalar expression
 * 
 * @tparam MeshT MEsh type
 * @tparam ExprT expression type of the metric
 * @param mesh mesh
 * @param e expression
 * @return Pch_element_type<MeshT, 1> continuous piecewise P1 function representing the metric
 */
template <typename MeshT, typename ExprT, std::enable_if_t<std::is_base_of_v<MeshBase<typename MeshT::index_type>, MeshT> && std::is_base_of_v<ExprBase,ExprT>,MeshT>* = nullptr>
Pch_element_type<MeshT, 1>
expr( std::shared_ptr<MeshT> const& mesh, ExprT const& e )
{
    auto Xh = Pch<1>( mesh );
    auto g = Xh->element();
    g.on( _range = elements( Xh->mesh() ), _expr = e );
    return g;
}
}
