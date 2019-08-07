/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2018-03-06

  Copyright (C) 2018 Feel++ Consortium

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

#ifndef FEELPP_TOOLBOXES_HEATFLUID_HPP
#define FEELPP_TOOLBOXES_HEATFLUID_HPP 1

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename HeatType, typename FluidType>
class HeatFluid : public ModelNumerical,
                  public std::enable_shared_from_this< HeatFluid<HeatType,FluidType> >
{

public:
    typedef ModelNumerical super_type;
    typedef HeatFluid<HeatType,FluidType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef HeatType heat_model_type;
    typedef std::shared_ptr<heat_model_type> heat_model_ptrtype;

    typedef FluidType fluid_model_type;
    typedef std::shared_ptr<fluid_model_type> fluid_model_ptrtype;

    // mesh
    typedef typename heat_model_type::mesh_type mesh_heat_type;
    typedef typename fluid_model_type::mesh_type mesh_fluid_type;
    typedef mesh_fluid_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;
    // exporter
    typedef Exporter<mesh_type,mesh_type::nOrder> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    HeatFluid( std::string const& prefix,
               std::string const& keyword = "heat-fluid",
               worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
               std::string const& subPrefix = "",
               ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"HeatFluidMesh.path"); }
    std::shared_ptr<std::ostringstream> getInfo() const override;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initPostProcess();
public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    void updateParameterValues();

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }

    heat_model_ptrtype const& heatModel() const { return M_heatModel; }
    heat_model_ptrtype heatModel() { return M_heatModel; }

    fluid_model_ptrtype const& fluidModel() const { return M_fluidModel; }
    fluid_model_ptrtype fluidModel() { return M_fluidModel; }

    std::map<std::string, elements_reference_wrapper_t<mesh_type> > const& rangeMeshElementsByMaterial() const { return M_rangeMeshElementsByMaterial; }

    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }

    //___________________________________________________________________________________//

    std::shared_ptr<TSBase> timeStepBase() { return this->heatModel()->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->heatModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateInHousePreconditioner( DataUpdateLinear & data ) const override;
    void updateInHousePreconditioner( DataUpdateJacobian & data ) const override;

    void postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const override;

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

    void updateLinearFluidSolver( DataUpdateLinear & data ) const;
private :

    heat_model_ptrtype M_heatModel;
    fluid_model_ptrtype M_fluidModel;

    //bool M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    //elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;
    // materials range
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;

    // physical parameter
    bool M_useNaturalConvection;
    double M_BoussinesqRefTemperature;
    vector_field_expression<nDim,1,2> M_gravityForce;

    // solver
    //std::string M_solverName;
    bool M_useSemiImplicitTimeScheme;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

    // post-process
    export_ptrtype M_exporter;
    std::set<std::string> M_postProcessFieldExportedHeatt, M_postProcessFieldExportedFluid;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_HEATFLUID_HPP
