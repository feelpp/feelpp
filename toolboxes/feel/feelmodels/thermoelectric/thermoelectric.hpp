/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-12-12

  Copyright (C) 2016 Feel++ Consortium

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
   \file thermoelectric.hpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-12-12
 */

#ifndef FEELPP_TOOLBOXES_THERMOELECTRIC_HPP
#define FEELPP_TOOLBOXES_THERMOELECTRIC_HPP 1

#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feelmodels/electric/electric.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename HeatType, typename ElectricType>
class ThermoElectric : public ModelNumerical,
                       public ModelPhysics<HeatType::convex_type::nDim>,
                       public std::enable_shared_from_this< ThermoElectric<HeatType,ElectricType> >
{

public:
    typedef ModelNumerical super_type;
    typedef ThermoElectric<HeatType,ElectricType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef HeatType heat_model_type;
    typedef std::shared_ptr<heat_model_type> heat_model_ptrtype;

    typedef ElectricType electric_model_type;
    typedef std::shared_ptr<electric_model_type> electric_model_ptrtype;

    // mesh
    typedef typename heat_model_type::mesh_type mesh_heat_type;
    typedef typename electric_model_type::mesh_type mesh_electric_type;
    typedef mesh_heat_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef MaterialsProperties<mesh_type::nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,mesh_type::nOrder> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    //___________________________________________________________________________________//
    // constructor
    ThermoElectric( std::string const& prefix,
                    std::string const& keyword = "thermo-electric",
                    worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                    std::string const& subPrefix = "",
                    ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;


private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initPostProcess() override;
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
    //elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    heat_model_ptrtype const& heatModel() const { return M_heatModel; }
    heat_model_ptrtype heatModel() { return M_heatModel; }

    electric_model_ptrtype const& electricModel() const { return M_electricModel; }
    electric_model_ptrtype electricModel() { return M_electricModel; }

    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }


    //___________________________________________________________________________________//

    std::shared_ptr<TSBase> timeStepBase() { return this->heatModel()->timeStepBase(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->heatModel()->timeStepBase(); }
    void startTimeStep();
    void updateTimeStep();

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( prefixvm( prefix, this->heatModel()->keyword() ) ),
                                                  this->electricModel()->modelFields( prefixvm( prefix, this->electricModel()->keyword() ) ) );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVectorHeat, size_type rowStartInVectorElectric, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( sol, rowStartInVectorHeat, prefixvm( prefix,this->heatModel()->keyword() ) ),
                                                  this->electricModel()->modelFields( sol, rowStartInVectorElectric, prefixvm( prefix,this->electricModel()->keyword() ) ) );
        }
    template <typename ModelFieldsHeatType,typename ModelFieldsElectricType>
    auto modelFields( ModelFieldsHeatType const& mfieldsHeat, ModelFieldsElectricType const& mfieldsElectric, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( mfieldsHeat, mfieldsElectric );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndexHeat, size_type startBlockSpaceIndexElectric ) const
        {
            return Feel::FeelModels::selectorModelFields( this->heatModel()->trialSelectorModelFields( startBlockSpaceIndexHeat ),
                                                          this->electricModel()->trialSelectorModelFields( startBlockSpaceIndexElectric ) );
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
    auto modelContext( vector_ptrtype sol, size_type startBlockSpaceIndexHeat, size_type startBlockSpaceIndexElectric, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, startBlockSpaceIndexHeat, startBlockSpaceIndexElectric, prefix );
            auto se = this->symbolsExpr( mfields ).template createTensorContext<mesh_type>();
            auto tse =  this->trialSymbolsExpr( mfields, trialSelectorModelFields( startBlockSpaceIndexHeat, startBlockSpaceIndexElectric ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seHeat = this->heatModel()->symbolsExprToolbox( mfields );
            auto seElectric = this->electricModel()->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            return Feel::vf::symbolsExpr( seHeat,seElectric,seParam,seMat,seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType, typename TrialSelectorModelFieldsType>
    auto trialSymbolsExpr( ModelFieldsType const& mfields, TrialSelectorModelFieldsType const& tsmf ) const
        {
            return mfields.trialSymbolsExpr( tsmf );
        }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinear_Heat( DataUpdateLinear & data ) const;
    void updateResidual_Heat( DataUpdateResidual & data ) const;
    void updateLinear_Electric( DataUpdateLinear & data ) const;

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

    //___________________________________________________________________________________//

    bool checkResults() const override
        {
            // several calls (not do in on line) to be sure that all check have been run
            bool checkThermoElectric = super_type::checkResults();
            bool checkHeat = this->heatModel()->checkResults();
            bool checkElectric = this->electricModel()->checkResults();
            return checkThermoElectric && checkHeat && checkElectric;
        }

private :
    heat_model_ptrtype M_heatModel;
    electric_model_ptrtype M_electricModel;

    // physical parameter
    std::string M_modelName;
    bool M_modelUseJouleEffect;
    materialsproperties_ptrtype M_materialsProperties;

    // solver
    std::string M_solverName;
    bool M_solverNewtonInitialGuessUseLinearThermoElectric,M_solverNewtonInitialGuessUseLinearHeat,M_solverNewtonInitialGuessUseLinearElectric;

    // post-process
    export_ptrtype M_exporter;
};

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_TOOLBOXES_THERMOELECTRIC_HPP
