
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 4 May 2021

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

#include <boost/preprocessor/cat.hpp>
#include <mmg/libmmg.h>
#include <parmmg/libparmmg.h>
#include <variant>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
namespace Feel
{
Feel::po::options_description
remesh_options( std::string const& prefix )
{
    Feel::po::options_description mmgoptions( "Remesh Options" );
    // clang-format off
    mmgoptions.add_options()
        ( prefixvm( prefix, "remesh.save").c_str(),Feel::po::value<std::string>()->default_value(""), "name of the filename prefix to save")
        ( prefixvm( prefix, "remesh.verbose").c_str(),Feel::po::value<int>(), "[-1..10], Tune level of verbosity")
        ( prefixvm( prefix, "remesh.debug").c_str(),Feel::po::value<int>(), "Turn on/off debug mode")

        ( prefixvm( prefix, "remesh.hmin").c_str(),Feel::po::value<double>(), "[val], Minimal mesh size")
        ( prefixvm( prefix, "remesh.hmax").c_str(),Feel::po::value<double>(), "[val], Maximal mesh size")
        ( prefixvm( prefix, "remesh.hsiz").c_str(),Feel::po::value<double>(), "[val], Constant mesh size")

        ( prefixvm( prefix, "remesh.nosizreq").c_str(),Feel::po::value<int>(), "[0/1], Allow/avoid overwritten of sizes at required points (advanced usage)")
        ( prefixvm( prefix, "remesh.hgradreq").c_str(),Feel::po::value<double>(), "[val], Control gradation on required entites (advanced usage)")
        ( prefixvm( prefix, "remesh.hausd").c_str(),Feel::po::value<double>(), "[val], Control global Hausdorff distance (on all the boundary surfaces of the mesh)")
        ;
    // clang-format on
    return mmgoptions;
}

template<int topoDim, int realDim>
void setMMGOptions( nl::json const& params, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh_, MMG5_pSol sol )
{
    if ( !params.contains( "remesh") )
        return;

    if ( std::holds_alternative<MMG5_pMesh>( mesh_ ) )
    {
        auto mesh = std::get<MMG5_pMesh>( mesh_ );
        if constexpr ( topoDim == 3 )
        {
            if ( params["remesh"].contains("verbose") )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_verbose, params["/remesh/verbose"_json_pointer].get<int>() ); /*!< [-1..10], Tune level of verbosity */
            if ( params["remesh"].contains("debug") )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_debug, params["/remesh/debug"_json_pointer].get<int>() ); /*!< [1/0], Turn on/off debug mode */

            if ( params["remesh"].contains( "hmin" ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hmin, params["/remesh/hmin"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "hmax" ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hmax, params["/remesh/hmax"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "hsiz" ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hsiz, params["/remesh/hsiz"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "nosizreq" ) )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_nosizreq, params["/remesh/nosizreq"_json_pointer].get<int>() );
            if ( params["remesh"].contains( "hgradreq" ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hgradreq, params["/remesh/hgradreq"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "opnbdy" ) )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_opnbdy, params["/remesh/opnbdy"_json_pointer].get<int>() );
            if ( params["remesh"].contains( "angle" ) )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_angle, params["/remesh/angle"_json_pointer].get<int>() );
        }
        else if constexpr ( topoDim == 2 && realDim == 3 )
        {
#if 0            
            if ( countoption( prefixvm(prefix,"remesh.verbose") ) )
                MMGS_Set_iparameter( mesh, sol, MMGS_IPARAM_verbose, ioption( prefixvm(prefix,"remesh.verbose") ) ); /*!< [-1..10], Tune level of verbosity */
            if ( countoption(  prefixvm(prefix,"remesh.debug" )) )
                MMGS_Set_iparameter( mesh, sol, MMGS_IPARAM_debug, ioption(  prefixvm(prefix,"remesh.debug") ) ); /*!< [1/0], Turn on/off debug mode */
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hmin" ) ) )
                MMGS_Set_dparameter( mesh, sol, MMGS_DPARAM_hmin, doption( _name = "remesh.hmin", _prefix = prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hmax" ) ) )
                MMGS_Set_dparameter( mesh, sol, MMGS_DPARAM_hmax, doption( _name = "remesh.hmax", _prefix = prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hsiz" ) ) )
                MMGS_Set_dparameter( mesh, sol, MMGS_DPARAM_hsiz, doption( _name = "remesh.hsiz", _prefix = prefix ) );
#endif                
        }
        else if constexpr ( topoDim == 2 && realDim == 2 )
        {
            if ( params["remesh"].contains( "verbose" ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_verbose, params["/remesh/verbose"_json_pointer].get<int>() ); /*!< [-1..10], Tune level of verbosity */
            if ( params["remesh"].contains( "debug" ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_debug, params["/remesh/debug"_json_pointer].get<int>() ); /*!< [1/0], Turn on/off debug mode */

            if ( params["remesh"].contains( "hmin" ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hmin, params["/remesh/hmin"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "hmax" ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hmax, params["/remesh/hmax"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "hsiz" ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hsiz, params["/remesh/hsiz"_json_pointer].get<double>() );
            if ( params["remesh"].contains( "nosizreq" ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_nosizreq, params["/remesh/nosizreq"_json_pointer].get<int>() );
            if ( params["remesh"].contains( "hgradreq" ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hgradreq, params["/remesh/hgradreq"_json_pointer].get<int>() );
            if ( params["remesh"].contains( "opnbdy" ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_opnbdy, params["/remesh/opnbdy"_json_pointer].get<int>() );
            if ( params["remesh"].contains( "angle" ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_angle, params["/remesh/angle"_json_pointer].get<int>() );
        }
    }
        
}

template void setMMGOptions<2,2>( nl::json const& j, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh, MMG5_pSol sol );
template void setMMGOptions<2,3>( nl::json const& j, std::variant<MMG5_pMesh,  PMMG_pParMesh> mesh, MMG5_pSol sol );
template void setMMGOptions<3,3>( nl::json const& j, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh, MMG5_pSol sol );
}