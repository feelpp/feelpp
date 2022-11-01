/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
       Date: 2021-12-09

       Copyright (C) Université de Strasbourg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#if !defined( FEELPP_MESH_MESHADAPTATIONMETRIC_HPP )
#define FEELPP_MESH_MESHADAPTATIONMETRIC_HPP 1
#include <regex>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelmesh/metric.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

template<typename MeshT>
class MeshAdaptationMetricGradedLevelSet
{
public:
    using mesh_t = MeshT;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;

    MeshAdaptationMetricGradedLevelSet() = default;
    MeshAdaptationMetricGradedLevelSet( nl::json const& j )
        : j_( j ) {}
    MeshAdaptationMetricGradedLevelSet( nl::json&& j )
        : j_( std::move( j ) ) {}

    void setProperties( nl::json const& j ) 
    {
        j_.merge_patch( j );
    }
    auto operator()( mesh_ptr_t const& m ) const
    {
        std::string t{ j_["metric/type"_json_pointer] };
        std::string e{ j_["metric/expr"_json_pointer] };
        if ( t == "gradedls")
        {
            std::vector<std::string> facet_markers = j_["metric/markers"_json_pointer];
            auto P1h = Pch<1>( m );
            
            if ( std::smatch sm;
                std::regex_match( e, sm, std::regex( "gradedls\\((.*),(.*)\\)" ) ) )
            {
                double hclose = std::stod( sm[1] );
                double hfar = std::stod( sm[2] );
                LOG( INFO ) << fmt::format( "[remesh] gradedfromls hclose={} hfar={}", hclose, hfar, 0. ) << std::endl;
                return idv( gradedfromls( P1h, markedfaces( m, facet_markers ), hclose, hfar, 0. ) );
            }
            else if ( std::smatch sm;
                      std::regex_match( e, sm, std::regex( "gradedls\\((.*),(.*),(.*)\\)" ) ) )
            {
                double hclose = std::stod( sm[1] );
                double hfar = std::stod( sm[2] );
                double cst = std::stod( sm[3] );
                LOG( INFO ) << fmt::format( "[remesh] gradedfromls hclose={} hfar={} percent={}", hclose, hfar, cst ) << std::endl;
                return idv( gradedfromls( P1h, markedfaces( m, facet_markers ), hclose, hfar, cst ) );
            }
        }
    throw std::logic_error(fmt::format("invalid json data {} assocated to MeshAdaptationMetricGradedLevelset",j_));
    }

private:
  nl::json j_ = R"({
                    "metric" : { 
                       "type":"gradedls"
                    }
                })"_json;
};

template <typename MeshT>
class MeshAdaptationMetricExpr
{
  public:
    using mesh_t = MeshT;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;

    MeshAdaptationMetricExpr() = default;
    MeshAdaptationMetricExpr( nl::json const& j )
        : j_( j ) {}
    MeshAdaptationMetricExpr( nl::json&& j )
        : j_( std::move( j ) ) {}

    void setProperties( nl::json const& j )
    {
        j_.merge_patch( j );
    }
    auto operator()( mesh_ptr_t const& m ) const
    {
        std::string t{ j_["metric/type"_json_pointer] };
        std::string e{ j_["metric/expr"_json_pointer] };
        if ( t == "gradedls" )
        {
            return expr( e );
        }
         throw std::logic_error(fmt::format("invalid json data {} assocated to MeshAdaptationMetricExpr",j_));
    }

  private:
    nl::json j_ = R"({
                    "metric" : { 
                       "type":"expr",
                       "expr":"h:h"
                    }
                })"_json;
};

} // namespace Feel


#endif