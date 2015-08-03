/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
   \file codegen_solidmec.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-17
 */

#ifndef FEELPP_CODEGEN_SOLIDMECHANICS_HPP
#define FEELPP_CODEGEN_SOLIDMECHANICS_HPP 1

#include <feel/feelmodels/solid/solidmecbase.hpp>

#include <feel/feelvf/vf.hpp>


#undef SOLIDMECHANICS
#undef SOLIDMECHANICS0
#undef SOLIDMECHANICS1
#undef SOLIDMECHANICS2
#include "feelmodelscoreconfig.h"
#undef NUMSOLID
#if defined( SOLIDMECHANICS )
#define NUMSOLID /**/
#endif
#if defined( SOLIDMECHANICS0 )
#define NUMSOLID 0 /**/
#endif
#if defined( SOLIDMECHANICS1 )
#define NUMSOLID 1 /**/
#endif
#if defined( SOLIDMECHANICS2 )
#define NUMSOLID 2 /**/
#endif

#define SOLIDMECHANICS_CLASS_NAME BOOST_PP_CAT(SolidMechanics,NUMSOLID)

#include "bctool.hpp"
#undef SOLIDMECHANICS_BC
#undef SOLIDMECHANICS_VOLUME_FORCE
#include "solid.bc"

namespace Feel
{
namespace FeelModels
{

class SOLIDMECHANICS_CLASS_NAME : public SOLIDMECHANICSBASE_CLASS_TYPE,
                                  public boost::enable_shared_from_this< SOLIDMECHANICS_CLASS_NAME >
{
public:
    typedef SOLIDMECHANICSBASE_CLASS_TYPE super_type;

    typedef SOLIDMECHANICS_CLASS_NAME self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    //___________________________________________________________________________________//

    SOLIDMECHANICS_CLASS_NAME( std::string prefix,
                               bool __buildMesh=true,
                               WorldComm const& _worldComm=Environment::worldComm(),
                               std::string subPrefix="",
                               std::string appliShortRepository=soption(_name="exporter.directory") );


    //___________________________________________________________________________________//
    // load config files
    void loadConfigBCFile();
    void loadConfigMeshFile( std::string const& geofilename );
    void loadConfigMeshFile1dReduced( std::string const& geofilename );

    // update for use
    void init( bool buildMethodNum = true );

    //___________________________________________________________________________________//
    // assembly using bc
    void updateNewtonInitialGuess( vector_ptrtype& U ) const;

    void updateBCDirichletStrongResidual( vector_ptrtype& R ) const;
    void updateBCNeumannResidual( vector_ptrtype& R ) const;
    void updateBCRobinResidual( element_displacement_type const& u, vector_ptrtype& R ) const;
    void updateBCFollowerPressureResidual(element_displacement_type const& u, vector_ptrtype& R ) const;
    void updateSourceTermResidual( vector_ptrtype& R ) const;

    void updateBCDirichletStrongJacobian(sparse_matrix_ptrtype& J) const;
    void updateBCFollowerPressureJacobian(element_displacement_type const& u, sparse_matrix_ptrtype& J) const;
    void updateBCRobinJacobian( sparse_matrix_ptrtype& J) const;

    void updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const;
    void updateBCNeumannLinearPDE( vector_ptrtype& F ) const;
    void updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;
    void updateSourceTermLinearPDE( vector_ptrtype& F ) const;

}; // SolidMechanics



} // FeelModels
} // Feel
#endif /* FEELPP_CODEGEN_SOLIDMECHANICS_HPP */
