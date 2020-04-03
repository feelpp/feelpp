/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP 1


#include <feel/feelmodels/coefficientformpdes/coefficientformpdebase.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>


#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
class CoefficientFormPDE : public CoefficientFormPDEBase<ConvexType>,
                           public std::enable_shared_from_this< CoefficientFormPDE<ConvexType,BasisUnknownType> >
{
public:
    typedef CoefficientFormPDEBase<ConvexType> super_type;
    using size_type = typename super_type::size_type;
    typedef CoefficientFormPDE<ConvexType,BasisUnknownType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // basis
    typedef BasisUnknownType basis_unknown_type;
    static const uint16_type nOrderUnknown = basis_unknown_type::nOrder;
    // function space unknown
    typedef FunctionSpace<mesh_type, bases<basis_unknown_type> > space_unknown_type;
    typedef std::shared_ptr<space_unknown_type> space_unknown_ptrtype;
    typedef typename space_unknown_type::element_type element_unknown_type;
    typedef std::shared_ptr<element_unknown_type> element_unknown_ptrtype;
    typedef typename space_unknown_type::element_external_storage_type element_unknown_external_storage_type;
    static constexpr bool unknown_is_scalar = space_unknown_type::is_scalar;
    static constexpr bool unknown_is_vectorial = space_unknown_type::is_vectorial;
    // materials properties
    typedef MaterialsProperties<mesh_type> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;
    // time scheme
    typedef Bdf<space_unknown_type> bdf_unknown_type;
    typedef std::shared_ptr<bdf_unknown_type> bdf_unknown_ptrtype;
    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_unknown_type> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;


    CoefficientFormPDE( typename super_type::super2_type const& genericPDE,
                        std::string const& prefix,
                        std::string const& keyword = "cfpde",
                        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                        std::string const& subPrefix  = "",
                        ModelBaseRepository const& modelRep = ModelBaseRepository() );

    CoefficientFormPDE( std::string const& prefix,
                        std::string const& keyword = "cfpde",
                        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                        std::string const& subPrefix  = "",
                        ModelBaseRepository const& modelRep = ModelBaseRepository() )
        :
        CoefficientFormPDE( typename super_type::super2_type(), prefix, keyword, worldComm, subPrefix, modelRep )
        {}

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"CoefficientFormPDEMesh.path"); }

    //___________________________________________________________________________________//
    // mesh, space, element unknown
    space_unknown_ptrtype const& spaceUnknown() const { return M_Xh; }
    element_unknown_ptrtype const& fieldUnknownPtr() const { return M_fieldUnknown; }
    element_unknown_type const& fieldUnknown() const { return *M_fieldUnknown; }
    //___________________________________________________________________________________//
#if 0
    // boundary condition + body forces
    map_scalar_field<2> const& bcDirichlet() const { return M_bcDirichlet; }
    map_scalar_field<2> const& bcNeumann() const { return M_bcNeumann; }
    map_scalar_fields<2> const& bcRobin() const { return M_bcRobin; }
    map_scalar_field<2> const& bodyForces() const { return M_volumicForcesProperties; }
#endif

    //___________________________________________________________________________________//
    // time step scheme
    std::string const& timeStepping() const { return M_timeStepping; }
    bdf_unknown_ptrtype const& timeStepBdfUnknown() const { return M_bdfUnknown; }
    //std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfUnknown(); }
    std::shared_ptr<TSBase> timeStepBase() const override { return this->timeStepBdfUnknown(); }
    void startTimeStep();
    void updateTimeStep();
    //___________________________________________________________________________________//

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;


    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();

    //___________________________________________________________________________________//
    // execute post-processing
    //___________________________________________________________________________________//

    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );

    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr )
        {
            return this->exportResults( time, this->modelFields(), symbolsExpr, this->exprPostProcessExports( symbolsExpr ) );
        }

    template <typename ModelFieldsType,typename SymbolsExpr>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr );

    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            return this->materialsProperties()->exprPostProcessExports( this->physics(),se );
        }

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    struct FieldTag
    {
        static auto unknown( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldUnknownPtr(), prefix );
        }

    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_t = this->spaceUnknown()->elementPtr( *sol, rowStartInVector + this->startSubBlockSpaceIndex( this->unknownName() ) );
            return this->modelFields( field_t, prefix );
        }

    template <typename TheUnknownFieldType>
    auto modelFields( TheUnknownFieldType const& field_u, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( modelField<FieldCtx::ID|FieldCtx::GRAD|FieldCtx::GRAD_NORMAL>( FieldTag::unknown(this), prefix, this->unknownName(), field_u, this->unknownSymbol(), this->keyword() ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            //auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols heat_T, heat_grad_T(_x,_y,_z), heat_dn_T
            return Feel::vf::symbolsExpr( /*seToolbox,*/ seParam, seMat, seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }


    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues ) override;

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;

    //___________________________________________________________________________________//
    // apply assembly and solver
    //___________________________________________________________________________________//

    void solve();

    template <typename ModelContextType>
    void updateLinearPDE( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType,typename RangeType>
    void updateLinearPDEStabilizationGLS( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const;

    template <typename ModelContextType>
    void updateNewtonInitialGuess( ModelAlgebraic::DataNewtonInitialGuess & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateJacobian( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType,typename RangeType>
    void updateJacobianStabilizationGLS( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const;
    void updateJacobianDofElimination( ModelAlgebraic::DataUpdateJacobian & data ) const override;
    template <typename ModelContextType>
    void updateResidual( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType,typename RangeType>
    void updateResidualStabilizationGLS( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx, std::string const& matName, RangeType const& range ) const;
    void updateResidualDofElimination( ModelAlgebraic::DataUpdateResidual & data ) const override;


private :
    void initFunctionSpaces();
    void initBoundaryConditions();
    void initTimeStep();
    void initInitialConditions();
    void initPostProcess() override;

private :


    space_unknown_ptrtype M_Xh;
    element_unknown_ptrtype M_fieldUnknown;

    // time discretisation
    std::string M_timeStepping;
    bdf_unknown_ptrtype M_bdfUnknown;
    double M_timeStepThetaValue;
    vector_ptrtype M_timeStepThetaSchemePreviousContrib;

    // boundary conditions
    using map_field_dirichlet = typename mpl::if_c<unknown_is_scalar,  map_scalar_field<2>,  map_vector_field<nDim,1,2> >::type;
    map_field_dirichlet/*map_scalar_field<2>*/ M_bcDirichlet;
    map_scalar_field<2> M_bcNeumann;
    map_scalar_fields<2> M_bcRobin;
    MarkerManagementDirichletBC M_bcDirichletMarkerManagement;
    MarkerManagementNeumannBC M_bcNeumannMarkerManagement;
    MarkerManagementRobinBC M_bcRobinMarkerManagement;

    // post-process
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;
};

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("CoefficientFormPDE","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->executePostProcessExports( this->M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr );
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : this->timeStepBase()->iteration(), mfields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("CoefficientFormPDE","exportResults", "finish");
}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelFieldsType, typename SymbolsExpr>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;
    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), this->rangeMeshElements(), symbolsExpr, mfields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), this->rangeMeshElements(), symbolsExpr, mfields );
    bool hasMeasurePoint = this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, mfields );
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

} // namespace Feel
} // namespace FeelModels

#include <feel/feelmodels/coefficientformpdes/coefficientformpdeassembly.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdeassemblystabilizationgls.hpp>

#endif
