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
   \file thermodynamics.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#include <feel/feelmodels/thermodyn/thermodynamics.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

namespace Feel
{
namespace FeelModels
{

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::ThermoDynamics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    :
    super_type( prefix, worldComm, buildMesh, subPrefix, rootRepository )
{
    this->log("ThermoDynamics","constructor", "start" );

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"ThermoDynamics.info") );
    //-----------------------------------------------------------------------------//
    // load info from .bc file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh, space, exporter,...
    if ( buildMesh ) this->build();
    //-----------------------------------------------------------------------------//
    this->log("ThermoDynamics","constructor", "finish");
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
{
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();

    M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Dirichlet" );
    for( auto const& d : M_bcDirichlet )
        this->addMarkerDirichletBC("elimination", marker(d) );
    M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Neumann" );
    for( auto const& d : M_bcNeumann )
        this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));

    M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::loadConfigMeshFile(std::string const& geofilename)
{
    CHECK( false ) << "not allow";
}


THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::init(bool buildMethodNum)
{
    super_type::init( buildMethodNum, this->shared_from_this() );
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::solve()
{
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    M_bcDirichlet.setParameterValues( paramValues );
    M_bcNeumann.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
    super_type::solve();
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( M_bcDirichlet.empty() ) return;

    this->log("ThermoDynamics","updateBCStrongDirichletLinearPDE","start" );

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = *this->fieldTemperature();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix(),
                                              _colstart=this->colStartInMatrix() );

    for( auto const& d : M_bcDirichlet )
    {
        bilinearForm_PatternCoupled +=
            on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                _element=u,_rhs=F,_expr=expression(d) );
    }

    this->log("ThermoDynamics","updateBCStrongDirichletLinearPDE","finish" );
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const
{
    if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        this->M_overwritemethod_updateSourceTermLinearPDE(F,buildCstPart);
        return;
    }

    if ( M_volumicForcesProperties.empty() ) return;

    if ( !buildCstPart )
    {
        auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=F,
                                   _rowstart=this->rowStartInVector() );
        auto const& v = *this->fieldTemperature();

        for( auto const& d : M_volumicForcesProperties )
        {
            if ( marker(d).empty() )
                myLinearForm +=
                    integrate( _range=elements(this->mesh()),
                               _expr= expression(d)*id(v),
                               _geomap=this->geomap() );
            else
                myLinearForm +=
                    integrate( _range=markedelements(this->mesh(),marker(d)),
                               _expr= expression(d)*id(v),
                               _geomap=this->geomap() );
        }
    }
}

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
    if ( M_bcNeumann.empty() ) return;

    if ( !buildCstPart )
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = *this->fieldTemperature();
        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );
        for( auto const& d : M_bcNeumann )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                           _expr= expression(d)*id(v),
                           _geomap=this->geomap() );
        }
    }
}


} // end namespace FeelModels
} // end namespace Feel
