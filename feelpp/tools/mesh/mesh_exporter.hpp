/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 22 Nov 2019

 Copyright (C) 2019 Feel++ Consortium

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


#include <map>

#include <feel/feelmesh/dump.hpp>
#include <boost/hana/pair.hpp>
#include <boost/hana/tuple.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/ginac.hpp>

namespace Feel
{
template <typename mesh_t>
void doExport(std::shared_ptr<mesh_t> mesh = {} )
{
    if ( !mesh )
        mesh = loadMesh( _mesh = new mesh_t );
    dump( mesh, "export mesh" );
    auto ex = exporter( _mesh = mesh );

    auto getfns = [=]( std::string const& sexpr )
    {
        std::vector<std::string> fields;
        boost::split( fields, sexpr, boost::is_any_of( "|" ), boost::token_compress_on );
        CHECK( fields.size() > 1 ) << "bad expression format";
        std::set<std::string> reps;
        for ( auto it = fields.begin() + 2; it != fields.end(); ++it )
            reps.insert( *it );
        return std::tuple{ fields[0], fields[1], reps };
    };

    for ( auto const& sexpr : vsoption( _name = "scalar_expr" ) )
    {
        auto const& [strname, strexpr, reps] = getfns( sexpr );
        ex->add( strname, expr( strexpr ), reps );
    }
    if ( boption( "export_pid" ) )
        ex->addRegions();
# if 0
    for ( auto const& sexpr : vsoption( _name = "vectorial_expr" ) )
    {
        auto const& [strname, strexpr, reps] = getfns( sexpr );
        ex->add( strname, expr<dimension_v<mesh_t>, 1>( strexpr ), reps );
    }
    ex->add( "P", P(), "nodal" );
    ex->add( "facesmarker", semarker(), faces( mesh ), "nodal" );
    if constexpr ( dimension_v < mesh_t > > 2 )
        ex->add( "edgesmarker", semarker(), edges( mesh ), "nodal" );
    ex->add( "pointmarker", semarker(), points( mesh ), "nodal" );

    using namespace std::string_literals;
    auto exprs = hana::make_tuple( hana::pair{ "detJ"s, detJ() }, hana::pair{ "marker"s, emarker() }, hana::pair{ "pid"s, epid() }, hana::pair{ "h"s, h() } );
    hana::for_each( exprs,
                    [&]( const auto& x )
                    { ex->add( hana::first( x ), hana::second( x ), "element" ); } );
#endif
    ex->save();
}
}