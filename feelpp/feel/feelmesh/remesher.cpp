
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
void setMMGOptions( std::string const& prefix, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh_, MMG5_pSol sol )
{
#if 0    
    if ( std::holds_alternative<MMG5_pMesh>( mesh_ ) )
    {
        auto mesh = std::get<MMG5_pMesh>( mesh_ );
        if constexpr ( topoDim == 3 )
        {
            if ( countoption( prefixvm( prefix, "remesh.verbose" ) ) )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_verbose, ioption( prefixvm( prefix, "remesh.verbose" ) ) ); /*!< [-1..10], Tune level of verbosity */
            if ( countoption( prefixvm( prefix, "remesh.debug" ) ) )
                MMG3D_Set_iparameter( mesh, sol, MMG3D_IPARAM_debug, ioption( prefixvm( prefix, "remesh.debug" ) ) ); /*!< [1/0], Turn on/off debug mode */

            if ( Environment::vm().count( prefixvm( prefix, "remesh.hmin" ) ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hmin, doption( _name = "remesh.hmin", _prefix = prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hmax" ) ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hmax, doption( _name = "remesh.hmax", _prefix = prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hsiz" ) ) )
                MMG3D_Set_dparameter( mesh, sol, MMG3D_DPARAM_hsiz, doption( _name = "remesh.hsiz", _prefix = prefix ) );
              
            //MMG3D_Set_iparameter( M_mmg_mesh, M_mmg_sol, MMG3D_IPARAM_verbose, value<MmgOption::Verbose>() );
            //MMG3D_Set_iparameter( M_mmg_mesh, M_mmg_sol, MMG3D_IPARAM_mem, value<MmgOption::Mem>() );
    
        }
        else if constexpr ( topoDim == 2 && realDim == 3 )
        {
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
        }
        else if constexpr ( topoDim == 2 && realDim == 2 )
        {
            if ( countoption( prefixvm( prefix, "remesh.verbose" ) ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_verbose, ioption( prefixvm( prefix, "remesh.verbose" ) ) ); /*!< [-1..10], Tune level of verbosity */
            if ( countoption( prefixvm( prefix, "remesh.debug" ) ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_debug, ioption( prefixvm( prefix, "remesh.debug" ) ) ); /*!< [1/0], Turn on/off debug mode */

            if (  Environment::vm().count( prefixvm(prefix, "remesh.hmin" ) ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hmin, doption( _name="remesh.hmin", _prefix=prefix ) );
            if (  Environment::vm().count( prefixvm(prefix, "remesh.hmax" ) ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hmax, doption( _name="remesh.hmax", _prefix=prefix ) );
            if (  Environment::vm().count( prefixvm(prefix, "remesh.hsiz" ) ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hsiz, doption( _name="remesh.hsiz", _prefix=prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.nosizreq" ) ) )
                MMG2D_Set_iparameter( mesh, sol, MMG2D_IPARAM_nosizreq, ioption( _name="remesh.nosizreq", _prefix=prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hgradreq" ) ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hgradreq, doption( _name = "remesh.hgradreq", _prefix = prefix ) );
            if ( Environment::vm().count( prefixvm( prefix, "remesh.hausd" ) ) )
                MMG2D_Set_dparameter( mesh, sol, MMG2D_DPARAM_hausd, doption( _name = "remesh.hausd", _prefix = prefix ) );

        }
    }
#endif
        
}

template void setMMGOptions<2,2>( std::string const& prefix, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh, MMG5_pSol sol );
template void setMMGOptions<2,3>( std::string const& prefix, std::variant<MMG5_pMesh,  PMMG_pParMesh> mesh, MMG5_pSol sol );
template void setMMGOptions<3,3>( std::string const& prefix, std::variant<MMG5_pMesh, PMMG_pParMesh> mesh, MMG5_pSol sol );
}