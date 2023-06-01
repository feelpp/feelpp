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
/**
 * @brief class for conjuguate heat transfer toolbox
 * @ingroup HeatFluid
 *
 * @tparam HeatType type of the heat transfer toolbox
 * @tparam FluidType type of the fluid mechanics toolbox
 *
 * @code {.cpp}
 * using heat_t FeelModels::Heat< Simplex<nDim,1>,
 *                           Lagrange<OrderT, Scalar,Continuous,PointSetFekete> >;
 * using fluid _t = FeelModels::FluidMechanics< Simplex<nDim,1>,
 *                                     Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
 *                                     Lagrange<OrderP, Scalar,Continuous,PointSetFekete> >;
 * using heatfluid_t = FeelModels::HeatFluid<heat_t, fluid_t>;
 * auto heatFluid = std::make_shared<heatfluid_t>("heat-fluid");
 * heatFluid->init();
 * heatFluid->printAndSaveInfo(); * 
 * if (heatFluid->isStationary() )
 * {
 *     heatFluid->solve();
 *     heatFluid->exportResults();
 * }
 * @endcode
 *
 */
template< typename HeatType, typename FluidType>
class HeatFluid : public ModelNumerical,
                  public ModelPhysics<HeatType::convex_type::nDim>
{
    typedef ModelPhysics<HeatType::convex_type::nDim> super_physics_type;
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

    // materials properties
    typedef MaterialsProperties<mesh_type::nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,mesh_type::nOrder> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    HeatFluid( std::string const& prefix,
               std::string const& keyword = "heatfluid",
               worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
               std::string const& subPrefix = "",
               ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;


private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initPostProcess() override;
    void initAlgebraicModel();
    void initAlgebraicFactory();
    void updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models ) override;
public :
    // update for use
    void init( bool buildModelAlgebraicFactory = true );

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    //___________________________________________________________________________________//

    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    void applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp = std::make_shared<RemeshInterpolation>() );

    heat_model_ptrtype const& heatModel() const { return M_heatModel; }
    heat_model_ptrtype heatModel() { return M_heatModel; }

    fluid_model_ptrtype const& fluidModel() const { return M_fluidModel; }
    fluid_model_ptrtype fluidModel() { return M_fluidModel; }

    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//

    std::shared_ptr<TSBase> timeStepBase() { return this->heatModel()->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->heatModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    template <typename SymbolsExprType>
    void updateInitialConditions( SymbolsExprType const& se );

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( this->heatModel()->keyword() ),
                                                  this->fluidModel()->modelFields( this->fluidModel()->keyword() ),
                                                  this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVectorHeat, size_type rowStartInVectorFluid, std::string const& prefix = "" ) const
        {
            return this->modelFields( sol, rowStartInVectorHeat, sol, rowStartInVectorFluid, prefix );
        }
    auto modelFields( vector_ptrtype solHeat, size_type rowStartInVectorHeat, vector_ptrtype solFluid, size_type rowStartInVectorFluid, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( solHeat, rowStartInVectorHeat, this->heatModel()->keyword() ),
                                                  this->fluidModel()->modelFields( solFluid, rowStartInVectorFluid, this->fluidModel()->keyword() ),
                                                  this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }
    auto modelFields( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataHeat,
                      std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataFluid,
                      std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( vectorDataHeat, this->heatModel()->keyword() ),
                                                  this->fluidModel()->modelFields( vectorDataFluid, this->fluidModel()->keyword() ),
                                                  this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }
    template <typename ModelFieldsHeatType,typename ModelFieldsFluidType>
    auto modelFields( ModelFieldsHeatType const& mfieldsHeat, ModelFieldsFluidType const& mfieldsFluid, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( mfieldsHeat, mfieldsFluid, this->template modelFieldsMeshes<mesh_type>( prefix ) );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndexHeat, size_type startBlockSpaceIndexFluid ) const
        {
            return Feel::FeelModels::selectorModelFields( this->heatModel()->trialSelectorModelFields( startBlockSpaceIndexHeat ),
                                                          this->fluidModel()->trialSelectorModelFields( startBlockSpaceIndexFluid ) );
        }


    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seHeat = this->heatModel()->symbolsExprToolbox( mfields );
            auto seFluid = this->fluidModel()->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMeshes = this->template symbolsExprMeshes<mesh_type,false>();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            auto sePhysics = this->symbolsExprPhysics( this->physics() );
            return Feel::vf::symbolsExpr( seHeat,seFluid,seParam,seMeshes,seMat,seFields,sePhysics );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }

    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    // template <typename ModelFieldsType>
    // auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
    //     {
    //         return Feel::FeelModels::modelContext( mfields, this->symbolsExpr( mfields ) );
    //     }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, size_type startBlockSpaceIndexHeat, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
        {
            return this->modelContext( sol, startBlockSpaceIndexHeat, sol, startBlockSpaceIndexFluid, prefix );
        }
    auto modelContext( vector_ptrtype solHeat, size_type startBlockSpaceIndexHeat, vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( solHeat, startBlockSpaceIndexHeat, solFluid, startBlockSpaceIndexFluid, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, this->trialSelectorModelFields( startBlockSpaceIndexHeat, startBlockSpaceIndexFluid ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }
    auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type startBlockSpaceIndexHeat, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
        {
            return this->modelContextNoTrialSymbolsExpr( sol, startBlockSpaceIndexHeat, sol, startBlockSpaceIndexFluid, prefix );
        }
    auto modelContextNoTrialSymbolsExpr( vector_ptrtype solHeat, size_type startBlockSpaceIndexHeat, vector_ptrtype solFluid, size_type startBlockSpaceIndexFluid, std::string const& prefix = "" ) const
        {
            // auto mfields = this->modelFields( solHeat, startBlockSpaceIndexHeat, solFluid, startBlockSpaceIndexFluid, prefix );
            // auto se = this->symbolsExpr( mfields );
            // return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
            return this->modelContextNoTrialSymbolsExpr( { { "solution", std::make_tuple( solHeat, startBlockSpaceIndexHeat ) } },
                                                         { { "solution", std::make_tuple( solFluid, startBlockSpaceIndexFluid ) } },
                                                         prefix );
        }
    auto modelContextNoTrialSymbolsExpr( std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataHeat,
                                         std::map<std::string,std::tuple<vector_ptrtype,size_type> > const& vectorDataFluid,
                                         std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( vectorDataHeat, vectorDataFluid, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
    auto modelContext( vector_ptrtype sol, heat_model_ptrtype const& heatModel, fluid_model_ptrtype const& fluidModel, std::string const& prefix = "" ) const
        {
            return this->modelContext( sol, heatModel->startBlockSpaceIndexVector(), fluidModel->startBlockSpaceIndexVector(), prefix );
        }

    //___________________________________________________________________________________//
    // apply assembly and solver
    //___________________________________________________________________________________//
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


    bool checkResults() const override
        {
            // several calls (not do in on line) to be sure that all check have been run
            bool checkHeatFluid = super_type::checkResults();
            bool checkHeat = this->heatModel()->checkResults();
            bool checkFluid = this->fluidModel()->checkResults();
            return checkHeatFluid && checkHeat && checkFluid;
        }

private :
    void updateLinear_Heat( DataUpdateLinear & data ) const;
    void updateResidual_Heat( DataUpdateResidual & data ) const;
    void updateJacobian_Heat( DataUpdateJacobian & data ) const;

    void updateLinear_Fluid( DataUpdateLinear & data ) const;
    void updateResidual_Fluid( DataUpdateResidual & data ) const;
    void updateJacobian_Fluid( DataUpdateJacobian & data ) const;

    void updateInHousePreconditioner_Fluid( DataUpdateLinear & data ) const;
    void updateInHousePreconditioner_Fluid( DataUpdateJacobian & data ) const;

    void updateTimeStepCurrentResidual();

private :

    heat_model_ptrtype M_heatModel;
    fluid_model_ptrtype M_fluidModel;

    // physical parameter
    bool M_useNaturalConvection;
    double M_BoussinesqRefTemperature;
    vector_field_expression<nDim,1,2> M_gravityForce;
    materialsproperties_ptrtype M_materialsProperties;

    // solver
    std::string M_solverName;
    bool M_useSemiImplicitTimeScheme;

    vector_ptrtype M_timeStepThetaSchemePreviousContrib, M_timeStepThetaSchemePreviousSolution;
    std::map<std::string,double> M_currentParameterValues;

    // post-process
    export_ptrtype M_exporter;
};

template< typename HeatType, typename FluidType>
template <typename SymbolsExprType>
void
HeatFluid<HeatType,FluidType>::updateInitialConditions( SymbolsExprType const& se )
{
    M_heatModel->updateInitialConditions( se );
    M_fluidModel->updateInitialConditions( se );
}


} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_HEATFLUID_HPP
