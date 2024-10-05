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
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>
#include <regex>
#include <type_traits>

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
template <typename SpaceT, typename RangeType, 
          std::enable_if_t<std::is_base_of_v<FunctionSpaceBase, SpaceT> && 
                           is_range_v<decay_type<RangeType>> &&
                           decay_type<RangeType>::entities() == MESH_FACES,SpaceT>* = nullptr>
typename SpaceT::element_type
gradedfromls( std::shared_ptr<SpaceT> const& Xh, RangeType&& facets, double hclose, double hfar, double cst = 0. )
{
    using namespace vf;
    tic();
    tic();
    auto d = distanceToRange( _space=Xh, _range=std::forward<RangeType>(facets) );
    toc("distanceToRange",no_display);
    d *= 1./d.max(); // max of d is 1
    auto g = Xh->element();
    CHECK( hclose > 0 && hfar > 0 ) << fmt::format("hclose({}) and/or hfar({}) must be strictly positive", hclose, hfar);
    LOG_IF(WARNING, ( (cst < 0) || (cst > 1) ) ) << fmt::format("[gradedfromls] invalid constant percentage {} (should be > 0 and < 1)",cst);

    cst = (cst < 0)?0:(cst > 1)?1:cst;
    g.on(_range=elements(Xh->mesh()),_expr=(idv(d) > cst)*(hclose*(1-idv(d))+hfar*idv(d))+(idv(d)<=cst)*hclose );
    toc("gradedfromls",no_display);
    return g;
}

/**
 * @brief Create a Graded Levelset metric
 * 
 * @tparam SpacePtrType P1 function space type
 * @tparam MetricSetter functor that takes a P1 function and set the metric
 * @param P1h the function space in which the metric is defined
 * @param metric_expr the expression of the metric
 * @param metric_setter the setter functor
 * @param default_markers the default markers to use if the metric_expr does not contain markers
 * @return true if the metric has been set
 * @return false otherwise
 */
template<typename SpacePtrType, typename MetricSetter>
bool
createGradedLevelsetMetric( SpacePtrType const& P1h, std::string const& metric_expr, MetricSetter metric_setter, std::vector<std::string> default_markers = {} )
{
    std::vector<std::string> markers = default_markers;
    double hclose = 0.05;
    double hfar = 0.1;
    double treshold = 0.;
    if ( std::smatch match;
             std::regex_search( metric_expr, match, std::regex("gradedls\\(\\[(.*?)\\],\\s*([0-9.]+),\\s*([0-9.]+),\\s*([0-9.]+)") ) )
        {
            std::string themarkers = match[1];
            hclose = std::stof(match[2]);
            hfar = std::stof(match[3]);
            treshold = std::stof(match[4]);
            std::regex marker_regex( "\\s*,\\s*" );
            std::vector<std::string> marker_vector{
                std::sregex_token_iterator( themarkers.begin(), themarkers.end(), marker_regex, -1 ),
                std::sregex_token_iterator() };
            if ( !marker_vector.empty() )
                markers = marker_vector;            
        }
        else if ( std::smatch sm;
             std::regex_match( metric_expr, sm, std::regex( "gradedls\\((.*),(.*)\\)" ) ) )
        {
            hclose = std::stod( sm[1] );
            hfar = std::stod( sm[2] );
        }
        else if ( std::smatch sm;
                  std::regex_match( metric_expr, sm, std::regex( "gradedls\\((.*),(.*),(.*)\\)" ) ) )
        {
            hclose = std::stod( sm[1] );
            hfar = std::stod( sm[2] );
            treshold = std::stod( sm[3] );
        }
        else 
            return false;
    LOG( INFO ) << fmt::format( "[remesh] gradedfromls markers=[{}] hclose={} hfar={} treshold={}", 
                                markers, hclose, hfar, treshold ) << std::endl;
    metric_setter( gradedfromls( P1h, markedfaces( P1h->mesh(), markers ), hclose, hfar, treshold ) ); 
    return true;
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
expr( std::shared_ptr<SpaceT> const& Xh, ExprT && e )
{
    auto g = Xh->element();
    g.on( _range = elements( Xh->mesh() ), _expr = std::forward<ExprT>(e) );
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
expr( std::shared_ptr<MeshT> const& mesh, ExprT && e )
{
    auto Xh = Pch<1>( mesh );
    auto g = Xh->element();
    g.on( _range = elements( Xh->mesh() ), _expr = std::forward<ExprT>(e) );
    return g;
}

} // namespace Feel
