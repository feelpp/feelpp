/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2014-06-04

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

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
   \file codegen_thermodyn.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */


#ifndef __THERMODYNAMICS_H
#define __THERMODYNAMICS_H 1

#include <feel/feelmodels/thermodyn/thermodynbase.hpp>

#undef THERMODYNAMICS
#undef THERMODYNAMICS0
#undef THERMODYNAMICS1
#undef THERMODYNAMICS2
#include "feelmodelscoreconfig.h"
#undef NUM_THERMO
#if defined( THERMODYNAMICS )
#define NUM_THERMO /**/
#endif
#if defined( THERMODYNAMICS0 )
#define NUM_THERMO 0 /**/
#endif
#if defined( THERMODYNAMICS1 )
#define NUM_THERMO 1 /**/
#endif
#if defined( THERMODYNAMICS2 )
#define NUM_THERMO 2 /**/
#endif
#define THERMODYNAMICS_CLASS_NAME BOOST_PP_CAT(ThermoDynamics,NUM_THERMO)

#include "bctool.hpp"
#undef THERMODYNAMICS_BC
#undef THERMODYNAMICS_VOLUME_FORCE
#include "thermodyn.bc"

namespace Feel
{
namespace FeelModels
{

class THERMODYNAMICS_CLASS_NAME : public THERMODYNAMICSBASE_CLASS_TYPE,
                                  public boost::enable_shared_from_this< THERMODYNAMICS_CLASS_NAME >
    {
    public:
        typedef THERMODYNAMICSBASE_CLASS_TYPE super_type;

        typedef THERMODYNAMICS_CLASS_NAME self_type;
        typedef boost::shared_ptr<self_type> self_ptrtype;

        // mesh
        typedef super_type::mesh_type mesh_type;
        typedef super_type::mesh_ptrtype mesh_ptrtype;

        //___________________________________________________________________________________//
        // constructor
        THERMODYNAMICS_CLASS_NAME( std::string prefix,
                                   bool __buildMesh=true,
                                   WorldComm const& _worldComm=Environment::worldComm(),
                                   std::string subPrefix="",
                                   std::string appliShortRepository=soption(_name="exporter.directory") );

        // load config files
        void loadConfigBCFile();
        void loadConfigMeshFile( std::string const& geofilename );

        // update for use
        void init( bool buildMethodNum = true );

        //___________________________________________________________________________________//
        // assembly using bc
        void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const;
        void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
        void updateSourceTermLinearPDE(vector_ptrtype& F, bool buildCstPart) const;

    };

} // namespace FeelModels
} // namespace Feel

#endif /* __THERMODYNAMICS_H */
