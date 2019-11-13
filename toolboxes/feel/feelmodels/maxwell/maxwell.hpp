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
   \file maxwell.hpp
   \author Romain Hild <romain.hild@unistra.fr>
   \date 2018-05-03
 */

#ifndef FEELPP_TOOLBOXES_MAXWELL_HPP
#define FEELPP_TOOLBOXES_MAXWELL_HPP 1

#include <boost/mpl/equal.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

#include <feel/feelmodels/maxwell/maxwellpropertiesdescription.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

enum class MaxwellPostProcessFieldExported
{
    MagneticPotential = 0, MagneticField, Pid
};

template< typename ConvexType>//, typename BasisPotentialType>
class Maxwell : public ModelNumerical,
                public MarkerManagementDirichletBC,
                public MarkerManagementNeumannBC,
                public MarkerManagementRobinBC,
                public std::enable_shared_from_this< Maxwell<ConvexType>>//,BasisPotentialType> >
{

public:
    typedef ModelNumerical super_type;
    typedef Maxwell<ConvexType>/*,BasisPotentialType>*/ self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // function space magnetic-potential
    // typedef BasisPotentialType basis_magneticpotential_type;
#if MAXWELL_H1
    using basis_magneticpotential_type = Lagrange<1, Vectorial, Continuous>;
#else
    using basis_magneticpotential_type = Nedelec<0, NedelecKind::NED1>;
#endif
    static const uint16_type nOrderPolyMagneticPotential = basis_magneticpotential_type::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_magneticpotential_type> > space_magneticpotential_type;
    typedef std::shared_ptr<space_magneticpotential_type> space_magneticpotential_ptrtype;
    typedef typename space_magneticpotential_type::element_type element_magneticpotential_type;
    typedef std::shared_ptr<element_magneticpotential_type> element_magneticpotential_ptrtype;
    typedef typename space_magneticpotential_type::element_external_storage_type element_magneticpotential_external_storage_type;

    // function space magnetic-field
    using basis_magneticfield_type = typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
                                                       RaviartThomas<0>,
                                                       Lagrange<0, Scalar, Discontinuous> >::type;
    typedef FunctionSpace<mesh_type, bases<basis_magneticfield_type> > space_magneticfield_type;
    typedef std::shared_ptr<space_magneticfield_type> space_magneticfield_ptrtype;
    typedef typename space_magneticfield_type::element_type element_magneticfield_type;
    typedef std::shared_ptr<element_magneticfield_type> element_magneticfield_ptrtype;

    // mechanical properties desc
    typedef bases<Lagrange<0, Scalar,Discontinuous> > basis_scalar_P0_type;
    typedef FunctionSpace<mesh_type, basis_scalar_P0_type> space_scalar_P0_type;
    typedef std::shared_ptr<space_scalar_P0_type> space_scalar_P0_ptrtype;
    typedef MaxwellPropertiesDescription<space_scalar_P0_type> maxwellproperties_type;
    typedef std::shared_ptr<maxwellproperties_type> maxwellproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_magneticpotential_type,space_magneticfield_type> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;
    // // context for evaluation
    // typedef typename space_magneticpotential_type::Context context_magneticpotential_type;
    // typedef std::shared_ptr<context_magneticpotential_type> context_magneticpotential_ptrtype;

    using map_dirichlet_field = typename mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<3> >,
                                                   map_vector_field<nDim>,
                                                   map_scalar_field<2> >::type;


    //___________________________________________________________________________________//
    // constructor
    Maxwell( std::string const& prefix,
             std::string const& keyword = "maxwell",
             worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
             std::string const& subPrefix = "",
             ModelBaseRepository const& modelRep = ModelBaseRepository() );
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"MaxwellMesh.path"); }
    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initBoundaryConditions();
    template<typename convex = convex_type>
    void initDirichlet(std::enable_if_t<convex::nDim==2>* = nullptr) { this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getScalarFields<nDim>( "magnetic-potential", "Dirichlet" ); }
    template<typename convex = convex_type>
    void initDirichlet(std::enable_if_t<convex::nDim==3>* = nullptr) { this->M_bcDirichlet = this->modelProperties().boundaryConditions().template getVectorFields<nDim>( "magnetic-potential", "Dirichlet" ); }
    void initPostProcess() override;

    template <typename FieldMagneticPotentialType>
    constexpr auto symbolsExprField( FieldMagneticPotentialType const& v ) const
        {
            return Feel::vf::symbolsExpr( symbolExpr("magnetic_A",idv(v) ) );
        }

    template <typename SymbExprType>
    auto symbolsExprFit( SymbExprType const& se ) const { return super_type::symbolsExprFit( se ); }

public :
    void setMesh(mesh_ptrtype const& mesh) { M_mesh = mesh; }
    // update for use
    void init( bool buildModelAlgebraicFactory = true );
    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    int nBlockMatrixGraph() const;
    void initAlgebraicFactory();

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr );
    template <typename TupleFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr );
    measure_points_evaluation_ptrtype& measurePointsEvaluation() { return M_measurePointsEvaluation; }

    void updateParameterValues();

    auto allFields() const
        {
            return hana::make_tuple( std::make_pair( "magnetic-potential",this->fieldMagneticPotentialPtr() ),
                                     std::make_pair( "magnetic-field",this->fieldMagneticFieldPtr() ),
                                     std::make_pair( "magnetic-permeability",M_maxwellProperties->fieldMagneticPermeabilityPtr() )
                                     );
        }

    template <typename FieldMagneticPotentialType>
    /*constexpr*/auto symbolsExpr( FieldMagneticPotentialType const& v ) const
        {
            auto seField = this->symbolsExprField( v );
            auto seFit = this->symbolsExprFit( seField );
            auto seMat = this->symbolsExprMaterial( Feel::vf::symbolsExpr( seField, seFit ) );
            return Feel::vf::symbolsExpr( seField, seFit, seMat );
        }
    auto symbolsExpr() const { return this->symbolsExpr( this->fieldMagneticPotential() ); }

    constexpr auto symbolsExprField() const { return this->symbolsExprField( this->fieldMagneticPotential() ); }

    template <typename SymbExprType>
    auto symbolsExprMaterial( SymbExprType const& se ) const
        {
            typedef decltype(expr(scalar_field_expression<2>{},se)) _expr_type;
            std::vector<std::pair<std::string,_expr_type>> matPropSymbs;
            for ( auto const& [_matname, _expr] : this->maxwellProperties()->magneticPermeabilityByMaterial() )
            {
                matPropSymbs.push_back( std::make_pair( (boost::format("magnetic_%1%_mu")%_matname).str(), expr( _expr.expr(), se ) ) );
            }
            return Feel::vf::symbolsExpr( symbolExpr( matPropSymbs ) );
        }

    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }

    space_magneticpotential_ptrtype const& spaceMagneticPotential() const { return M_XhMagneticPotential; }
    element_magneticpotential_ptrtype const& fieldMagneticPotentialPtr() const { return M_fieldMagneticPotential; }
    element_magneticpotential_type const& fieldMagneticPotential() const { return *M_fieldMagneticPotential; }

    space_magneticfield_ptrtype const& spaceMagneticField() const { return M_XhMagneticField; }
    element_magneticfield_ptrtype const& fieldMagneticFieldPtr() const { return M_fieldMagneticField; }
    element_magneticfield_type const& fieldMagneticField() const { return *M_fieldMagneticField; }

    maxwellproperties_ptrtype const& maxwellProperties() const { return M_maxwellProperties; }

    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    //___________________________________________________________________________________//
    // apply assembly and solver
    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename SymbolsExpr>
    void updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const;
    // void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename SymbolsExpr>
    void updateLinearPDEWeakBC( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr, hana::int_<3> ) const;
    template <typename SymbolsExpr>
    void updateLinearPDEWeakBC( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr, hana::int_<2> ) const;
    // void updateLinearPDEStrongDirichletBC( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const;

    //___________________________________________________________________________________//
    template <typename SymbolsExpr>
    void updateFields( SymbolsExpr const& symbolsExpr )
        {
            this->maxwellProperties()->updateFields( symbolsExpr );
            this->updateMagneticField();
        }
    void updateMagneticField();
    template<typename SymbolsExpr>
    void updateMagneticField( SymbolsExpr const& symbolsExpr )
        {
            auto const& v = this->fieldMagneticPotential();
            for ( auto const& rangeData : this->maxwellProperties()->rangeMeshElementsByMaterial() )
            {
                std::string const& matName = rangeData.first;
                auto const& range = rangeData.second;
                auto const& magneticPermeability = this->maxwellProperties()->magneticPermeability( matName );
                auto muExpr = expr( magneticPermeability.expr(), symbolsExpr );
                M_fieldMagneticField->on(_range=range, _expr=curlv(v)/muExpr);
            }
        }

private :
    bool M_hasBuildFromMesh, M_isUpdatedForUse;

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    space_magneticpotential_ptrtype M_XhMagneticPotential;
    element_magneticpotential_ptrtype M_fieldMagneticPotential;
    space_magneticfield_ptrtype M_XhMagneticField;
    element_magneticfield_ptrtype M_fieldMagneticField;
    // physical parameter
    maxwellproperties_ptrtype M_maxwellProperties;
    // boundary conditions
    map_dirichlet_field M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    map_vector_field<nDim> M_volumicForcesProperties;

    // regularization
    double M_epsilon;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
    // std::map<std::string,std::set<size_type> > M_dofsWithValueImposed;
    // // start dof index fields in matrix (temperature,maxwell-potential,...)
    // std::map<std::string,size_type> M_startBlockIndexFieldsInMatrix;

    // post-process
    export_ptrtype M_exporter;
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;
    // std::set<std::string> M_postProcessFieldExported;
    // std::set<std::string> M_postProcessUserFieldExported;

}; // class Maxwell

template< typename ConvexType>
template <typename SymbolsExpr>
void
Maxwell<ConvexType>::exportResults( double time, SymbolsExpr const& symbolsExpr )
{
    this->log("Maxwell","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );

    auto fields = this->allFields();
    this->executePostProcessExports( M_exporter, time, fields );
    this->executePostProcessMeasures( time, fields, symbolsExpr );
    this->executePostProcessSave( invalid_uint32_type_value, fields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("Maxwell","exportResults", "finish");
}


template< typename ConvexType>
template <typename TupleFieldsType,typename SymbolsExpr>
void
Maxwell<ConvexType>::executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;
    bool hasMeasureNorm = this->executePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasureStatistics = this->executePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, tupleFields, symbolsExpr );
    bool hasMeasurePoint = this->executePostProcessMeasuresPoint( M_measurePointsEvaluation, tupleFields );
    if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
        hasMeasure = true;

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}



} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/maxwell/maxwellassembly.hpp>

#endif // FEELPP_TOOLBOXES_MAXWELL_HPP
