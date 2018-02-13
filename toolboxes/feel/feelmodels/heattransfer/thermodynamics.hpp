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
   \file thermodyn.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-06-04
 */

#ifndef FEELPP_THERMODYNAMICS_HPP
#define FEELPP_THERMODYNAMICS_HPP 1

#include <feel/feelmodels/thermodyn/thermodynbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisTemperatureType>
class ThermoDynamics : public ThermoDynamicsBase<ConvexType,BasisTemperatureType>,
                       public boost::enable_shared_from_this< ThermoDynamics<ConvexType,BasisTemperatureType> >
{

public:
    typedef ThermoDynamicsBase<ConvexType,BasisTemperatureType> super_type;

    typedef ThermoDynamics<ConvexType,BasisTemperatureType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    ThermoDynamics( std::string const& prefix,
                    bool buildMesh = true,
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& appliShortRepository = "" );

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

#endif // INCLUDE_THERMODYNAMICS_HPP
