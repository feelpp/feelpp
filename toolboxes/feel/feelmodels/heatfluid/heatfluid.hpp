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

#include <feel/feelmodels/heattransfer/heattransfer.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename HeatTransferType, typename FluidType>
class HeatFluid : public ModelNumerical,
                  public boost::enable_shared_from_this< HeatFluid<HeatTransferType,FluidType> >
{

public:
    typedef ModelNumerical super_type;
    typedef HeatFluid<HeatTransferType,FluidType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef HeatTransferType heattransfer_model_type;
    typedef boost::shared_ptr<heattransfer_model_type> heattransfer_model_ptrtype;

    typedef FluidType fluid_model_type;
    typedef boost::shared_ptr<fluid_model_type> fluid_model_ptrtype;

    // mesh
    typedef typename heattransfer_model_type::mesh_type mesh_heattransfer_type;
    typedef typename fluid_model_type::mesh_type mesh_fluid_type;
    typedef mesh_fluid_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    static const uint16_type nDim = mesh_type::nDim;
    // exporter
    typedef Exporter<mesh_type,mesh_type::nOrder> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    HeatFluid( std::string const& prefix,
               bool buildMesh = true,
               WorldComm const& _worldComm = Environment::worldComm(),
               std::string const& subPrefix = "",
               ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"HeatFluidMesh.path"); }
    boost::shared_ptr<std::ostringstream> getInfo() const;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initPostProcess();
public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );

    BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    int nBlockMatrixGraph() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    void updateParameterValues();

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    //elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    heattransfer_model_ptrtype const& heatTransferModel() const { return M_heatTransferModel; }
    heattransfer_model_ptrtype heatTransferModel() { return M_heatTransferModel; }

    fluid_model_ptrtype const& fluidModel() const { return M_fluidModel; }
    fluid_model_ptrtype fluidModel() { return M_fluidModel; }


    backend_ptrtype const& backend() const { return M_backendMonolithic; }
    BlocksBaseVector<double> const& blockVectorSolutionMonolithic() const { return M_blockVectorSolutionMonolithic; }
    BlocksBaseVector<double> & blockVectorSolutionMonolithic() { return M_blockVectorSolutionMonolithic; }

    size_type startBlockSpaceIndex( std::string const& name )
        {
            auto itFind = M_startBlockSpaceIndex.find( name );
            if ( itFind != M_startBlockSpaceIndex.end() )
                return itFind->second;
            return invalid_size_type_value;
        }
    //___________________________________________________________________________________//

    boost::shared_ptr<TSBase> timeStepBase() { return this->heatTransferModel()->timeStepBase(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->heatTransferModel()->timeStepBase(); }
    void updateTimeStep();

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const;

    void updateLinearPDE( DataUpdateLinear & data ) const;

    void updateNewtonInitialGuess( vector_ptrtype& U ) const;
    void updateJacobian( DataUpdateJacobian & data ) const;
    void updateResidual( DataUpdateResidual & data ) const;

private :
    heattransfer_model_ptrtype M_heatTransferModel;
    fluid_model_ptrtype M_fluidModel;

    //bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    //elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;
    // materials range
    std::map<std::string, elements_reference_wrapper_t<mesh_type> > M_rangeMeshElementsByMaterial;

    // physical parameter
    bool M_useNaturalConvection;
    double M_BoussinesqRefTemperature;
    vector_field_expression<nDim,1,2> M_gravityForce;
    //bool M_modelUseJouleEffect;

    // solver
    //std::string M_solverName;
    //bool M_solverNewtonInitialGuessUseLinearThermoElectric,M_solverNewtonInitialGuessUseLinearHeatTransfer,M_solverNewtonInitialGuessUseLinearElectric;

    // algebraic data/tools
    backend_ptrtype M_backendMonolithic;
    model_algebraic_factory_ptrtype M_algebraicFactoryMonolithic;
    BlocksBaseVector<double> M_blockVectorSolutionMonolithic;
    std::map<std::string,size_type> M_startBlockSpaceIndex;

    // post-process
    export_ptrtype M_exporter;
    std::set<std::string> M_postProcessFieldExportedHeatTransfert, M_postProcessFieldExportedFluid;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_HEATFLUID_HPP
