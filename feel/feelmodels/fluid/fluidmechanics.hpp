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
 \file fluidmechanics.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#ifndef FEELPP_FLUIDMECHANICS_HPP
#define FEELPP_FLUIDMECHANICS_HPP 1

#include <feel/feelmodels/fluid/fluidmecbase.hpp>

namespace Feel
{
namespace FeelModels
{
template <typename ConvexType, typename BasisVelocityType,
          typename BasisPressureType = Lagrange<( BasisVelocityType::nOrder > 1 ) ? ( BasisVelocityType::nOrder - 1 ) : BasisVelocityType::nOrder, Scalar, Continuous, PointSetFekete>,
          typename BasisDVType = Lagrange<0, Scalar, Discontinuous /*,PointSetFekete*/>,
          bool UsePeriodicity = false>
class FluidMechanics : public FluidMechanicsBase<ConvexType, BasisVelocityType, BasisPressureType, BasisDVType, UsePeriodicity>,
                       public boost::enable_shared_from_this<FluidMechanics<ConvexType, BasisVelocityType, BasisPressureType, BasisDVType, UsePeriodicity>>
{
  public:
    typedef FluidMechanicsBase<ConvexType, BasisVelocityType, BasisPressureType, BasisDVType, UsePeriodicity> super_type;
    typedef FluidMechanics<ConvexType, BasisVelocityType, BasisPressureType, BasisDVType, UsePeriodicity> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    using element_velocity_type = typename super_type::element_fluid_velocity_type;

    //___________________________________________________________________________________//
    // constructor
    FluidMechanics( std::string const& prefix,
                    bool buildMesh = true,
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    FluidMechanics( self_type const& FM ) = default;
    //___________________________________________________________________________________//
    static self_ptrtype New( std::string const& prefix,
                             bool buildMesh = true,
                             WorldComm const& worldComm = Environment::worldComm(),
                             std::string const& subPrefix = "",
                             std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    //___________________________________________________________________________________//
    // load config files
    void loadConfigBCFile();
    void loadConfigPostProcess();
    void loadConfigMeshFile( std::string const& geofilename );
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    void solve();
    //___________________________________________________________________________________//
    // assembly using bc
    void updateSourceTermResidual( vector_ptrtype& R ) const;
    void updateInitialNewtonSolutionBCDirichlet( vector_ptrtype& U ) const;
    void updateBCStrongDirichletJacobian( sparse_matrix_ptrtype& J, vector_ptrtype& RBis ) const;
    void updateBCStrongDirichletResidual( vector_ptrtype& R ) const;
    void updateBCDirichletLagMultResidual( vector_ptrtype& R ) const;
    void updateBCDirichletNitscheResidual( vector_ptrtype& R ) const;
    void updateBCNeumannResidual( vector_ptrtype& R ) const;
    void updateBCPressureResidual( vector_ptrtype& R ) const;

    void updateSourceTermLinearPDE( vector_ptrtype& F, bool BuildCstPart ) const;
    void updateBCStrongDirichletLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;
    void updateBCDirichletLagMultLinearPDE( vector_ptrtype& F ) const;
    void updateBCDirichletNitscheLinearPDE( vector_ptrtype& F ) const;
    void updateBCNeumannLinearPDE( vector_ptrtype& F ) const;
    void updateBCPressureLinearPDE( vector_ptrtype& F ) const;

    void updateInHousePreconditionerPCD( sparse_matrix_ptrtype const& mat, vector_ptrtype const& vecSol ) const;

    //___________________________________________________________________________________//

    bool hasDirichletBC() const
    {
        return ( !M_bcDirichlet.empty() ||
                 !M_bcDirichletComponents.find( Component::X )->second.empty() ||
                 !M_bcDirichletComponents.find( Component::Y )->second.empty() ||
                 !M_bcDirichletComponents.find( Component::Z )->second.empty() );
    }

  private:
    map_vector_field<super_type::nDim, 1, 2> M_bcDirichlet;
    std::map<ComponentType, map_scalar_field<2>> M_bcDirichletComponents;
    map_scalar_field<2> M_bcNeumannScalar, M_bcPressure;
    map_vector_field<super_type::nDim, 1, 2> M_bcNeumannVectorial;
    map_matrix_field<super_type::nDim, super_type::nDim, 2> M_bcNeumannTensor2;

    map_vector_field<super_type::nDim, 1, 2> M_volumicForcesProperties;

}; // FluidMechanics

} // namespace FeelModels
} // namespace Feel

#endif /* FEELPP_FLUIDMECHANICS_HPP */
