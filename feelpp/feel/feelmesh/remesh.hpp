/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 6 Oct 2021

 Copyright (C) 2021 Feel++ Consortium

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

#include <fmt/core.h>
#include <regex>
#include <feel/feelmesh/remesher.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelmesh/metric.hpp>
#include <feel/feelfilters/partitionio.hpp>
namespace Feel
{
/**
 * @brief given a metric and required elements/facets, remesh a mesh
 * 
 * @tparam MeshT type of the mesh to be remeshed
 * @param r mesh to be remeshed
 * @param metric_expr metric expression
 * @param req_elts required elements set
 * @param req_facets  required facets set
 * @param parent parent mesh (possibly nullptr)
 * @return std::shared_ptr<MeshT> 
 */
template <typename MeshT>
std::tuple<std::shared_ptr<MeshT>,int>
remesh( std::shared_ptr<MeshT> const& r,
        std::string const& metric_expr,
        std::vector<std::string> const& req_elts,
        std::vector<std::string> const& req_facets,
        std::shared_ptr<MeshT> const& parent )
{
    static int cpt = 0;
    using mesh_t = MeshT;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    if ( r->worldComm().localSize() == 1 )
    {
        auto R = std::make_shared<Remesh<mesh_t>>( r, req_elts, req_facets );
        auto P1h = Pch<1>( r );

        if ( std::smatch sm; 
             std::regex_match( metric_expr, sm, std::regex( "gradedls\\((.*),(.*)\\)" ) )  )
        {
            double hclose = std::stod( sm[1] );
            double hfar = std::stod( sm[2] );
            LOG(INFO) << fmt::format( "[remesh] gradedfromls hclose={} hfar={}", hclose, hfar, 0. ) << std::endl;
            R->setMetric( gradedfromls( P1h, markedfaces( r, req_facets ), hclose, hfar, 0. ) );
        }
        else if ( std::smatch sm; 
             std::regex_match( metric_expr, sm, std::regex( "gradedls\\((.*),(.*),(.*)\\)" ) )  )
        {
            double hclose = std::stod( sm[1] );
            double hfar = std::stod( sm[2] );
            double cst = std::stod( sm[3] );
            LOG(INFO) << fmt::format( "[remesh] gradedfromls hclose={} hfar={} percent={}", hclose, hfar, cst ) << std::endl;
            R->setMetric( gradedfromls( P1h, markedfaces( r, req_facets ), hclose, hfar, cst ) );
        }
        else
        {
            R->setMetric( expr( P1h, expr( metric_expr ) ) );
        }

        auto m = R->execute();
        return std::tuple{ m, ++cpt };
    }
    else // parallel
    {
        std::string imeshParaName = fmt::format( "themesh{}_i.json", cpt );
        std::string omeshParaName = fmt::format( "themesh{}_o.json", cpt );
        std::string imeshParaPath = ( fs::path( Environment::appRepository() ) / imeshParaName ).string();
        std::string omeshParaPath = ( fs::path( Environment::appRepository() ) / imeshParaName ).string();
        int nPartition = r->worldComm().localSize();
        r->saveHDF5( imeshParaPath );
        r->worldComm().barrier();
        if ( r->worldComm().isMasterRank() )
        {
            auto meshSeq = loadMesh( _mesh = new mesh_t( Environment::worldCommSeqPtr() ), _savehdf5 = 0,
                                     _filename = imeshParaPath,
                                     //_update=update_,  TODO test FACE_MINIMAL
                                     _straighten = false );
            auto P1h = Pch<1>( meshSeq );
            auto R = std::make_shared<Remesh<mesh_t>>( meshSeq, req_elts, req_facets );

            if ( std::smatch sm; std::regex_match( metric_expr, sm, std::regex( "gradedls\\((.*),(.*)\\)" ) ) )
            {
                double hclose = std::stod( sm[1] );
                double hfar = std::stod( sm[2] );
                std::cout << fmt::format( "hclose={} hfar={}", hclose, hfar ) << std::endl;
                R->setMetric( gradedfromls( P1h, markedfaces( meshSeq, req_facets ), hclose, hfar ) );
            }
            else
            {
                R->setMetric( expr( P1h, expr( metric_expr ) ) );
            }
            auto m = R->execute();
            using io_t = PartitionIO<mesh_t>;
            io_t io( omeshParaPath );
            io.write( partitionMesh( m, nPartition /*, partitionByRange, partconfig*/ ) );
        }
        r->worldComm().barrier();
        auto meshAdaptPara = loadMesh( _mesh = new mesh_t, _savehdf5 = 0,
                                       _filename = omeshParaPath /*, _update=update_*/ );
        return std::tuple{ meshAdaptPara, ++cpt };
    }
}
}
