/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-17

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file codegen_fluidmec.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-17
 */


#ifndef FEELPP_CODEGEN_FLUIDMECHANICS_HPP
#define FEELPP_CODEGEN_FLUIDMECHANICS_HPP 1

#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>

#undef FLUIDMECHANICS
#undef FLUIDMECHANICS0
#undef FLUIDMECHANICS1
#undef FLUIDMECHANICS2
#include "feelmodelscoreconfig.h"
#undef NUMFLUID
#if defined( FLUIDMECHANICS )
#define NUMFLUID /**/
#endif
#if defined( FLUIDMECHANICS0 )
#define NUMFLUID 0 /**/
#endif
#if defined( FLUIDMECHANICS1 )
#define NUMFLUID 1 /**/
#endif
#if defined( FLUIDMECHANICS2 )
#define NUMFLUID 2 /**/
#endif

#define FLUIDMECHANICS_CLASS_NAME BOOST_PP_CAT(FluidMechanics,NUMFLUID)

#include "bctool.hpp"
#undef FLUIDMECHANICS_BC
#undef FLUIDMECHANICS_VOLUME_FORCE
#include "fluid.bc"

namespace Feel
{

namespace FeelModels
    {

        class FLUIDMECHANICS_CLASS_NAME : public FLUIDMECHANICSBASE_CLASS_TYPE,
                                          public boost::enable_shared_from_this< FLUIDMECHANICS_CLASS_NAME >
        {
        public:
            typedef FLUIDMECHANICSBASE_CLASS_TYPE super_type;

            typedef FLUIDMECHANICS_CLASS_NAME self_type;
            typedef boost::shared_ptr<self_type> self_ptrtype;
            //___________________________________________________________________________________//
            //___________________________________________________________________________________//
            //___________________________________________________________________________________//

            //typedef decltype( FLUIDMECHANICS_BC( self_ptrtype() ) ) bcdef_type;
            //typedef mpl::bool_<bcdef_type::hasThisType<cl::paroi_mobile>::value> hasBcParoiMobile_type;

            //___________________________________________________________________________________//
            // constructor
            FLUIDMECHANICS_CLASS_NAME( bool __isStationary,
                                       std::string prefix,
                                       WorldComm const& _worldComm=WorldComm(),
                                       bool __buildMesh=true,
                                       std::string subPrefix="",
                                       std::string appliShortRepository=option(_name="exporter.directory").as<std::string>() );
            //___________________________________________________________________________________//
            // constructor by copy
#if 0
            FLUIDMECHANICS_CLASS_NAME( self_type const & M )
                :
                super_type(M)
            {
                std::cout <<"\nWARNING constructeur copie\n";
            }
#endif
            //___________________________________________________________________________________//

            // load config files
            void loadConfigBCFile();
            void loadConfigMeshFile( std::string const& geofilename );

            // update for use
            void init(bool buildMethodNum=true);

            //___________________________________________________________________________________//
            // assembly using bc
            void updateSourceTermResidual( vector_ptrtype& R ) const;
            void updateCLDirichlet(vector_ptrtype& U) const;
            void updateBCStrongDirichletJacobian(sparse_matrix_ptrtype& J) const;
            void updateBCStrongDirichletResidual(vector_ptrtype& R) const;
            void updateBCDirichletLagMultResidual( vector_ptrtype& R ) const;
            void updateBCDirichletNitscheResidual( vector_ptrtype& R ) const;
            void updateBCNeumannResidual( vector_ptrtype& R ) const;
            void updateBCPressureResidual( vector_ptrtype& R ) const;

            void updateSourceTermLinearPDE( vector_ptrtype& F, bool BuildCstPart ) const;
            void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
            void updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const;
            void updateBCDirichletNitscheLinearPDE( vector_ptrtype& F ) const;
            void updateBCNeumannLinearPDE( vector_ptrtype& F ) const;
            void updateBCPressureLinearPDE( vector_ptrtype& F ) const;

            typedef boost::function<void ( vector_ptrtype& F, bool buildCstPart )> updateSourceTermLinearPDE_function_type;
            updateSourceTermLinearPDE_function_type M_overwritemethod_updateSourceTermLinearPDE;
            typedef boost::function<void ( vector_ptrtype& R )> updateSourceTermResidual_function_type;
            updateSourceTermResidual_function_type M_overwritemethod_updateSourceTermResidual;

        }; // FluidMechanics

        //---------------------------------------------------------------------------------------------------------//
        //---------------------------------------------------------------------------------------------------------//


    } // namespace FeelModels

} // namespace Feel






#endif /* FEELPP_CODEGEN_FLUIDMECHANICS_HPP */

