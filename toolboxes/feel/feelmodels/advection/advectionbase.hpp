/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2016-05-04

 Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

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
 \file advectionbase.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-05-04
 */

#ifndef _ADVECTIONBASE_HPP
#define _ADVECTIONBASE_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>
#include <feel/feelmodels/advection/diffusionreactionmodel.hpp>

namespace Feel {
namespace FeelModels {

namespace detail {

template<uint16_type, typename T> struct ChangeBasisOrder;

template<
    uint16_type NewOrder,
    template<uint16_type, template<uint16_type> class, typename, template<class, uint16_type, class> class, uint16_type > class BasisType,
    uint16_type Order,
    template<uint16_type> class PolySetType,
    typename ContinuityType,
    template<class, uint16_type, class> class Pts,
    uint16_type Tag
        >
struct ChangeBasisOrder<NewOrder, BasisType<Order, PolySetType, ContinuityType, Pts, Tag>>
{
    typedef BasisType<NewOrder, PolySetType, ContinuityType, Pts, Tag> type;
};

template<template<uint16_type> class, typename T> struct ChangeBasisPolySet;

template<
    template<uint16_type> class NewPolySetType,
    template<uint16_type, template<uint16_type> class, typename, template<class, uint16_type, class> class, uint16_type > class BasisType,
    uint16_type Order,
    template<uint16_type> class PolySetType,
    typename ContinuityType,
    template<class, uint16_type, class> class Pts,
    uint16_type Tag
    >
struct ChangeBasisPolySet<NewPolySetType, BasisType<Order, PolySetType, ContinuityType, Pts, Tag>>
{
    typedef BasisType<Order, NewPolySetType, ContinuityType, Pts, Tag> type;
};

} // namespace detail

enum class AdvectionStabMethod { NONE=0, GALS, CIP, SUPG, SGS };

template< 
    typename ConvexType, typename BasisAdvectionType, 
    typename PeriodicityType = NoPeriodicity,
    typename BasisDiffusionCoeffType = typename detail::ChangeBasisPolySet<Scalar, BasisAdvectionType>::type,
    typename BasisReactionCoeffType = typename detail::ChangeBasisPolySet<Scalar, BasisAdvectionType>::type
        >
class AdvectionBase : 
    public ModelNumerical,
    public MarkerManagementDirichletBC,
    public MarkerManagementNeumannBC

{
public :
    typedef ModelNumerical super_type;

    typedef AdvectionBase< ConvexType, BasisAdvectionType, PeriodicityType, 
            BasisDiffusionCoeffType, BasisReactionCoeffType > self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    //--------------------------------------------------------------------//
    // Mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //--------------------------------------------------------------------//
    // Space advection
    typedef BasisAdvectionType basis_advection_type;
    static const uint16_type nOrder = basis_advection_type::nOrder;

    typedef PeriodicityType periodicity_type;
    
    typedef FunctionSpace< mesh_type, bases<basis_advection_type>, Periodicity<periodicity_type> > space_advection_type;
    typedef boost::shared_ptr<space_advection_type> space_advection_ptrtype;
    
    typedef typename space_advection_type::element_type element_advection_type;
    typedef boost::shared_ptr<element_advection_type> element_advection_ptrtype;

    typedef typename space_advection_type::value_type value_type;
    typedef typename space_advection_type::periodicity_type periodicity_advection_type;

    static constexpr bool is_scalar = space_advection_type::is_scalar;
    static constexpr bool is_vectorial = space_advection_type::is_vectorial;
    static constexpr bool is_continuous = space_advection_type::is_continuous;

    //--------------------------------------------------------------------//
    // Space advection velocity
    //typedef typename basis_advection_type::vectorial_basis_type basis_vectorial_type;
    //typedef Lagrange<nOrder, Vectorial, basis_advection_continuity_type, basis_advection_pointset_type> basis_vectorial_type;
    typedef Lagrange<nOrder, Vectorial, Continuous, PointSetFekete> basis_vectorial_type;

    typedef bases<basis_vectorial_type> basis_advection_velocity_type;
    typedef FunctionSpace< 
        mesh_type, 
        basis_advection_velocity_type, 
        value_type, 
        periodicity_advection_type > space_advection_velocity_type;
    typedef boost::shared_ptr<space_advection_velocity_type> space_advection_velocity_ptrtype;
    typedef typename space_advection_velocity_type::element_type element_advection_velocity_type;
    typedef boost::shared_ptr<element_advection_velocity_type> element_advection_velocity_ptrtype;

    //--------------------------------------------------------------------//
    // Diffusion-reaction model
    typedef BasisDiffusionCoeffType basis_diffusioncoeff_type;
    typedef BasisReactionCoeffType basis_reactioncoeff_type;
    static const uint16_type nOrderDiffusionCoeff = BasisDiffusionCoeffType::nOrder;
    static const uint16_type nOrderReactionCoeff = BasisReactionCoeffType::nOrder;
    typedef FunctionSpace< mesh_type, bases<basis_diffusioncoeff_type> > space_diffusioncoeff_type;
    typedef boost::shared_ptr<space_diffusioncoeff_type> space_diffusioncoeff_ptrtype;
    typedef FunctionSpace< mesh_type, bases<basis_reactioncoeff_type> > space_reactioncoeff_type;
    typedef boost::shared_ptr<space_reactioncoeff_type> space_reactioncoeff_ptrtype;

    typedef DiffusionReactionModel<space_diffusioncoeff_type, space_reactioncoeff_type> diffusionreaction_model_type;
    typedef boost::shared_ptr<diffusionreaction_model_type> diffusionreaction_model_ptrtype;

    //--------------------------------------------------------------------//
    // Backend
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //--------------------------------------------------------------------//
    // Algebraic tools
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef boost::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    typedef typename model_algebraic_factory_type::graph_type graph_type;
    typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    typedef typename model_algebraic_factory_type::indexsplit_type indexsplit_type;
    typedef typename model_algebraic_factory_type::indexsplit_ptrtype indexsplit_ptrtype;

    //--------------------------------------------------------------------//
    // Time
    typedef Bdf<space_advection_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    //--------------------------------------------------------------------//
    // Exporter
    typedef Exporter<mesh_type, nOrderGeo> exporter_type;
    typedef boost::shared_ptr<exporter_type> exporter_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    AdvectionBase( 
            std::string const& prefix,
            WorldComm const& _worldComm = Environment::worldComm(),
            std::string const& subPrefix = "",
            std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );

    AdvectionBase( self_type const& A ) = default;

    void build();
    void build( mesh_ptrtype const& mesh );
    void build( space_advection_ptrtype const& space );

    //--------------------------------------------------------------------//
    // Initialization
    void init( bool buildModelAlgebraicFactory, model_algebraic_factory_type::appli_ptrtype const& app );
    //void initFromMesh( 
            //mesh_ptrtype const& mesh,
            //bool buildModelAlgebraicFactory, 
            //model_algebraic_factory_type::appli_ptrtype const& app );

    //--------------------------------------------------------------------//
    // Periodicity
    void setPeriodicity( periodicity_type const& p );
    periodicity_type const& periodicity() const { return M_periodicity; }

    //--------------------------------------------------------------------//
    // Model and solver
    std::string const& modelName() const { return M_modelName; }
    void setModelName( std::string const& type );

    bool hasDiffusion() const;
    bool hasReaction() const;

    std::string const& solverName() const { return M_solverName; }
    void setSolverName( std::string const& type );

    //--------------------------------------------------------------------//
    virtual std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"AdvectionMesh.path"); }

    mesh_ptrtype const& mesh() const { return M_mesh; }

    space_advection_ptrtype const& functionSpace() const { return M_Xh; }
    space_advection_velocity_ptrtype const& functionSpaceAdvectionVelocity() const { return M_XhAdvectionVelocity; }
    space_diffusioncoeff_ptrtype const& functionSpaceDiffusionCoeff() const { return this->diffusionReactionModel()->functionSpaceDiffusion(); }
    space_reactioncoeff_ptrtype const& functionSpaceReactionCoeff() const { return this->diffusionReactionModel()->functionSpaceReaction(); }
    
    bool useExtendedDofTable() const;

    element_advection_ptrtype & fieldSolutionPtr() { return M_fieldSolution; }
    element_advection_ptrtype const& fieldSolutionPtr() const { return M_fieldSolution; }
    element_advection_type & fieldSolution() { return *M_fieldSolution; }
    element_advection_type const& fieldSolution() const { return *M_fieldSolution; }

    element_advection_velocity_type const& fieldAdvectionVelocity() const { return *M_fieldAdvectionVelocity; }
    element_advection_velocity_ptrtype const& fieldAdvectionVelocityPtr() const { return M_fieldAdvectionVelocity; }

    //--------------------------------------------------------------------//
    // Algebraic data
    backend_ptrtype const& backend() const { return M_backend; }
    size_type matrixPattern() const;
    virtual int nBlockMatrixGraph() const;
    virtual BlocksBaseGraphCSR buildBlockMatrixGraph() const;
    graph_ptrtype buildMatrixGraph() const;
    virtual int buildBlockVectorSolution();
    void buildVectorSolution();
    //indexsplit_ptrtype buildIndexSplit() const;
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    virtual size_type nLocalDof() const;

    //--------------------------------------------------------------------//
    // Time scheme
    bdf_ptrtype timeStepBDF() { return M_bdf; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    virtual void updateTimeStepBDF();
    void updateTimeStep() { this->updateTimeStepBDF(); }
    void initTimeStep();

    //--------------------------------------------------------------------//
    // Stabilization
    AdvectionStabMethod stabilizationMethod() const { return M_stabMethod; }
    double stabilizationCoefficient() const { return M_stabCoeff; }
    void setStabilizationCoefficient(double coeff) { M_stabCoeff = coeff; }

    //--------------------------------------------------------------------//
    // Algebraic model updates
    // Linear PDE
    void updateLinearPDE( DataUpdateLinear & data ) const;
    virtual void updateLinearPDEAdditional( sparse_matrix_ptrtype & A, vector_ptrtype & F, bool _BuildCstPart ) const {}
    virtual void updateLinearPDEStabilization( sparse_matrix_ptrtype& A, vector_ptrtype& F, bool buildCstPart ) const;
    virtual void updateSourceTermLinearPDE( element_advection_ptrtype& fieldSource, bool buildCstPart ) const {}
    virtual bool hasSourceTerm() const =0;
    virtual void updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const =0;
    virtual void updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const =0;
    
    void updateBCNeumannLinearPDE( vector_ptrtype& F ) const;
    
    //--------------------------------------------------------------------//
    // Advection velocity update
    template<typename ExprT>
    void updateAdvectionVelocity(vf::Expr<ExprT> const& expr);
    //--------------------------------------------------------------------//
    // Source update
    template<typename ExprT>
    void updateSourceAdded(vf::Expr<ExprT> const& expr);
    bool hasSourceAdded() const { return M_hasSourceAdded; }
    //--------------------------------------------------------------------//
    // Diffusion-reaction parameters update
    diffusionreaction_model_ptrtype & diffusionReactionModel() { return M_diffusionReactionModel; }
    diffusionreaction_model_ptrtype const& diffusionReactionModel() const { return M_diffusionReactionModel; }

    void updateDiffusionCoeff(double D)
    {
        this->diffusionReactionModel()->setCstDiffusionCoeff(D);
    }
    void updateReactionCoeff(double R)
    {
        this->diffusionReactionModel()->setCstReactionCoeff(R);
    }

    template<typename ExprT>
    void updateDiffusionCoeff(vf::Expr<ExprT> const& expr)
    {
        this->diffusionReactionModel()->updateDiffusionCoeff(expr);
    }
    template<typename ExprT>
    void updateReactionCoeff(vf::Expr<ExprT> const& expr)
    {
        this->diffusionReactionModel()->updateReactionCoeff(expr);
    }

    //--------------------------------------------------------------------//
    // Solve
    virtual void solve();

    //--------------------------------------------------------------------//
    // Export results
    void initPostProcess();
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    void exportMeasures( double time );

    exporter_ptrtype getExporter() { return M_exporter; }
    exporter_ptrtype const& getExporter() const { return M_exporter; }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
protected:
    void loadParametersFromOptionsVm();

    void createMesh();
    void createFunctionSpaces();
    void createAlgebraicData();
    void createTimeDiscretization();
    void createExporters();
    void createOthers();

    virtual void updateLinearPDETransient( sparse_matrix_ptrtype& A, vector_ptrtype& F, bool buildCstPart ) const;

    virtual std::string geoExportType() const { return "static"; }
    virtual void exportResultsImpl( double time );
    virtual void exportMeasuresImpl( double time );

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    bool M_isBuilt;
    bool M_isUpdatedForUse;
    //--------------------------------------------------------------------//
    // Model and solver
    std::string M_modelName;
    std::string M_solverName;
    //--------------------------------------------------------------------//
    // Mesh
    mesh_ptrtype M_mesh;
    // Periodicity
    periodicity_type M_periodicity;
    // Advection space
    space_advection_ptrtype M_Xh;
    // Time discretization
    bdf_ptrtype M_bdf;
    //--------------------------------------------------------------------//
    // Algebraic data
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
    //--------------------------------------------------------------------//
    // Advection velocity
    space_advection_velocity_ptrtype M_XhAdvectionVelocity;
    element_advection_velocity_ptrtype M_fieldAdvectionVelocity;
    boost::optional<vector_field_expression<nDim,1,2> > M_exprAdvectionVelocity;
    //--------------------------------------------------------------------//
    // Physical parameters (diffusivity and reaction coefficient)
    diffusionreaction_model_ptrtype M_diffusionReactionModel;
    //--------------------------------------------------------------------//
    // Source added
    element_advection_ptrtype M_fieldSource;
    bool M_hasSourceAdded;
    //--------------------------------------------------------------------//
    // Solution
    element_advection_ptrtype M_fieldSolution;
    //--------------------------------------------------------------------//
    // Export
    exporter_ptrtype M_exporter;
    bool M_doExportAll;
    bool M_doExportAdvectionVelocity;
    bool M_doExportDiffusionCoefficient;
    bool M_doExportReactionCoefficient;
    bool M_doExportSourceField;
    //--------------------------------------------------------------------//
    // Stabilization
    static const std::map<std::string, AdvectionStabMethod> AdvectionStabMethodIdMap;
    AdvectionStabMethod M_stabMethod;
    double M_stabCoeff;
    double M_gamma1;
    
};//AdvectionBase

//----------------------------------------------------------------------------//
// Advection velocity update
template< 
    typename ConvexType, typename BasisAdvectionType, 
    typename PeriodicityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ExprT>
void
AdvectionBase<ConvexType, BasisAdvectionType, PeriodicityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateAdvectionVelocity(
        vf::Expr<ExprT> const& v_expr)
{
    M_exprAdvectionVelocity.reset(); // remove symbolic expr
    M_fieldAdvectionVelocity->on(_range=elements(this->mesh()), _expr=v_expr );
}

//----------------------------------------------------------------------------//
// Source update
template< 
    typename ConvexType, typename BasisAdvectionType, 
    typename PeriodicityType,
    typename BasisDiffusionCoeffType,
    typename BasisReactionCoeffType
        >
template<typename ExprT>
void 
AdvectionBase<ConvexType, BasisAdvectionType, PeriodicityType, BasisDiffusionCoeffType, BasisReactionCoeffType>::updateSourceAdded(
        vf::Expr<ExprT> const& f_expr)
{
    if (!M_fieldSource)
    {
        M_fieldSource.reset( new element_advection_type(M_Xh, "SourceAdded") );
    }
    M_fieldSource->on(_range=elements( this->mesh() ), _expr=f_expr );
    M_hasSourceAdded=true;
}

} // namespace FeelModels
} // namespace Feel


#endif
