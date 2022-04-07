/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDE_HPP 1


#include <feel/feelmodels/coefficientformpdes/coefficientformpdebase.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelpoly/nedelec.hpp>

#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdeboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisUnknownType>
class CoefficientFormPDE : public CoefficientFormPDEBase<ConvexType>
{
public:
    typedef CoefficientFormPDEBase<ConvexType> super_type;
    using Coefficient = typename super_type::Coefficient;

    using size_type = typename super_type::size_type;
    typedef CoefficientFormPDE<ConvexType,BasisUnknownType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nRealDim = convex_type::nRealDim;
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
    // time scheme
    typedef Bdf<space_unknown_type> bdf_unknown_type;
    typedef std::shared_ptr<bdf_unknown_type> bdf_unknown_ptrtype;


    CoefficientFormPDE( typename super_type::super2_type::infos_ptrtype const& infosPDE,
                        std::string const& prefix,
                        std::string const& keyword = "cfpde",
                        worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                        std::string const& subPrefix  = "",
                        ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    //! return true is the unknown is scalar
    bool unknownIsScalar() const override { return unknown_is_scalar; }

    //___________________________________________________________________________________//
    // mesh, space, element unknown
    space_unknown_ptrtype const& spaceUnknown() const { return M_Xh; }
    element_unknown_ptrtype const& fieldUnknownPtr() const { return M_fieldUnknown; }
    element_unknown_type const& fieldUnknown() const { return *M_fieldUnknown; }

    //___________________________________________________________________________________//
    // time step scheme
    bdf_unknown_ptrtype const& timeStepBdfUnknown() const { return M_bdfUnknown; }
    //std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBdfUnknown(); }
    std::shared_ptr<TSBase> timeStepBase() const override { return this->timeStepBdfUnknown(); }
    void startTimeStep() override;
    void updateTimeStep() override;

    //! update initial conditions with symbols expression \se
    template <typename SymbolsExprType>
    void updateInitialConditions( SymbolsExprType const& se );

    //___________________________________________________________________________________//

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();


    template <typename SymbolsExprType>
    bool hasSymbolDependencyInBoundaryConditions( std::set<std::string> const& symbs, SymbolsExprType const& se ) const;

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

    template <typename ModelFieldsType,typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities );

    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            return this->materialsProperties()->exprPostProcessExports( this->mesh(), this->physicsAvailable(), se );
        }

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    struct FieldTag
    {
        static auto unknown( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
        static auto unknown_previous( self_type const* t ) { return ModelFieldTag<self_type,1>( t ); }
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
            return Feel::FeelModels::modelFields( modelField<FieldCtx::FULL>( FieldTag::unknown(this), prefix, this->unknownName(), field_u, this->unknownSymbol(), this->keyword() ),
                                                  modelField<FieldCtx::FULL>( FieldTag::unknown_previous(this), prefix, this->unknownName()+"_previous", this->fieldUnknownPtr(), this->unknownSymbol() + "_previous", this->keyword() ),
                                                  modelField<FieldCtx::FULL>( FieldTag::unknown(this), prefix, this->unknownName()+"_remove_trial", field_u, this->unknownSymbol() + "_rt", this->keyword() )
                                                  );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::unknown(this), this->unknownName(), startBlockSpaceIndex ) );
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr(); // generate symbols heat_T, heat_grad_T(_x,_y,_z), heat_dn_T
            return Feel::vf::symbolsExpr( seToolbox, seParam, seMat, seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            using _expr_first_time_derivative_rhs_type = std::decay_t<decltype(idv(M_bdfUnknown->polyDeriv()))>;
            symbol_expression_t<_expr_first_time_derivative_rhs_type> se_firstTimeDerivative_rhs;
            using _expr_first_time_derivative_lhs_type = std::decay_t<decltype( expr<space_unknown_type::nComponents1,space_unknown_type::nComponents2>( "" )*cst(1.) )>;
            symbol_expression_t<_expr_first_time_derivative_lhs_type> se_firstTimeDerivative_lhs;

            using _expr_first_time_derivative_type = std::decay_t<decltype( _expr_first_time_derivative_lhs_type{} - _expr_first_time_derivative_rhs_type{} )>;
            symbol_expression_t<_expr_first_time_derivative_type> se_firstTimeDerivative;

            if ( M_bdfUnknown )
            {
                SymbolExprComponentSuffix secs( space_unknown_type::nComponents1,space_unknown_type::nComponents2 );

                // first time derivative : rhs expression
                auto exprFirstTimeDerivative_rhs = idv(M_bdfUnknown->polyDeriv());

                // first time derivative : lhs expression
                std::string symbolEvalUnknown = prefixvm( this->keyword(), this->unknownSymbol(), "_" );
                std::string symbolicExprFirstTimeDerivative_lhs = (boost::format("%1%:%1%")%symbolEvalUnknown).str();
                if ( space_unknown_type::nComponents1 == 1 && space_unknown_type::nComponents2 == 1 )
                    symbolicExprFirstTimeDerivative_lhs = (boost::format("%1%:%1%")%symbolEvalUnknown).str();
                else
                {
                    // TODO : move this part in SymbolExprComponentSuffix
                    std::vector<std::string> vecOfCompSuffix( secs.nComp1()*secs.nComp2() );
                    for ( auto const& [_suffix,compArray] : secs )
                    {
                        uint16_type c1 = compArray[0];
                        uint16_type c2 = compArray[1];
                        vecOfCompSuffix[c1*secs.nComp2()+c2] = _suffix;
                    }
                    symbolicExprFirstTimeDerivative_lhs = "{";
                    for ( int k=0;k<vecOfCompSuffix.size();++k )
                    {
                        if ( k>0 )
                            symbolicExprFirstTimeDerivative_lhs += ",";
                        symbolicExprFirstTimeDerivative_lhs += symbolEvalUnknown + vecOfCompSuffix[k];
                    }
                    symbolicExprFirstTimeDerivative_lhs += "}";
                    for ( int k=0;k<vecOfCompSuffix.size();++k )
                        symbolicExprFirstTimeDerivative_lhs += ":" + symbolEvalUnknown + vecOfCompSuffix[k];
                }
                //std::cout << "symbolicExprFirstTimeDerivative_lhs = " << symbolicExprFirstTimeDerivative_lhs << std::endl;
                auto exprFirstTimeDerivative_lhs = expr<space_unknown_type::nComponents1,space_unknown_type::nComponents2>( symbolicExprFirstTimeDerivative_lhs, "", this->worldComm(), this->repository().expr() )*cst(M_bdfUnknown->polyDerivCoefficient(0));

                // first time derivative : rhs-lhs expression
                auto exprFirstTimeDerivative = exprFirstTimeDerivative_lhs - exprFirstTimeDerivative_rhs;

                //<eqname>_d<unknown>_dt_rhs (example:  heat_dT_dt_rhs)
                std::string symbolFirstTimeDerivative_rhs = prefixvm( this->keyword(), (boost::format("d%1%_dt_rhs")%this->unknownSymbol()).str(), "_");
                se_firstTimeDerivative_rhs.add( symbolFirstTimeDerivative_rhs, std::move(exprFirstTimeDerivative_rhs), secs );
                //<eqname>_d<unknown>_dt_lhs (example:  heat_dT_dt_lhs)
                std::string symbolFirstTimeDerivative_lhs = prefixvm( this->keyword(), (boost::format("d%1%_dt_lhs")%this->unknownSymbol()).str(), "_");
                se_firstTimeDerivative_lhs.add( symbolFirstTimeDerivative_lhs, std::move(exprFirstTimeDerivative_lhs), secs );
                //<eqname>_d<unknown>_dt (example:  heat_dT_dt)
                std::string symbolFirstTimeDerivative = prefixvm( this->keyword(), (boost::format("d%1%_dt")%this->unknownSymbol()).str(), "_");
                se_firstTimeDerivative.add( symbolFirstTimeDerivative,std::move(exprFirstTimeDerivative),secs );
            }

            return Feel::vf::symbolsExpr( se_firstTimeDerivative_rhs,se_firstTimeDerivative_lhs,se_firstTimeDerivative );
        }

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
    void updateLinearPDEStabilizationGLS( ModelAlgebraic::DataUpdateLinear & data, ModelContextType const& mctx,
                                          std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                          std::string const& matName, RangeType const& range ) const;

    template <typename ModelContextType>
    void updateNewtonInitialGuess( ModelAlgebraic::DataNewtonInitialGuess & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateJacobian( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType,typename RangeType>
    void updateJacobianStabilizationGLS( ModelAlgebraic::DataUpdateJacobian & data, ModelContextType const& mctx,
                                         std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                         std::string const& matName, RangeType const& range ) const;
    void updateJacobianDofElimination( ModelAlgebraic::DataUpdateJacobian & data ) const override;
    template <typename ModelContextType>
    void updateResidual( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType,typename RangeType>
    void updateResidualStabilizationGLS( ModelAlgebraic::DataUpdateResidual & data, ModelContextType const& mctx,
                                         std::shared_ptr<ModelPhysicCoefficientFormPDE<nDim>> physicCFPDEData,
                                         std::string const& matName, RangeType const& range ) const;
    void updateResidualDofElimination( ModelAlgebraic::DataUpdateResidual & data ) const override;


private :
    void initFunctionSpaces();
    void initBoundaryConditions();
    void initTimeStep();
    void initPostProcess() override;

private :


    space_unknown_ptrtype M_Xh;
    element_unknown_ptrtype M_fieldUnknown;

    // time discretisation
    bdf_unknown_ptrtype M_bdfUnknown;

    // boundary conditions
    using boundary_conditions_type = CoefficientFormPDEBoundaryConditions<nRealDim,unknown_is_scalar?0:1>;
    std::shared_ptr<boundary_conditions_type> M_boundaryConditions;
};

template< typename ConvexType, typename BasisUnknownType>
template <typename SymbolsExprType>
bool
CoefficientFormPDE<ConvexType,BasisUnknownType>::hasSymbolDependencyInBoundaryConditions( std::set<std::string> const& symbs, SymbolsExprType const& se ) const
{
    bool hasDependency = false;
    for ( auto const& [bcname,bcData] : M_boundaryConditions->neumann() )
    {
        auto neumannExpr = bcData->expr();
        if ( neumannExpr.hasSymbolDependency( symbs, se ) )
        {
            hasDependency = true;
            break;
        }
    }
    if ( hasDependency )
        return true;

    for ( auto const& [bcname,bcData] : M_boundaryConditions->robin() )
    {
        auto theExpr1 = bcData->expr1();
        if ( theExpr1.hasSymbolDependency( symbs, se ) )
        {
            hasDependency = true;
            break;
        }
        auto theExpr2 = bcData->expr2();
        if ( theExpr2.hasSymbolDependency( symbs, se ) )
        {
            hasDependency = true;
            break;
        }
    }
    return hasDependency;
}

template< typename ConvexType, typename BasisUnknownType>
template <typename SymbolsExprType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::updateInitialConditions( SymbolsExprType const& se )
{
    if ( !this->doRestart() )
    {
        std::vector<element_unknown_ptrtype> icFields;
        std::map<int, double> icPriorTimes;
        if ( this->isStationary() )
        {
            icFields = { this->fieldUnknownPtr() };
            icPriorTimes = {{0,0}};
        }
        else
        {
            icFields = this->timeStepBdfUnknown()->unknowns();
            icPriorTimes = this->timeStepBdfUnknown()->priorTimes();
        }

        // auto paramValues = this->modelProperties().parameters().toParameterValues();
        // this->modelProperties().initialConditions().setParameterValues( paramValues );

        super_type::updateInitialConditions( this->unknownName(), this->rangeMeshElements(), se, icFields, icPriorTimes );

        if ( !this->isStationary() )
            *this->fieldUnknownPtr() = this->timeStepBdfUnknown()->unknown(0);
    }
}

template< typename ConvexType, typename BasisUnknownType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("CoefficientFormPDE","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->executePostProcessExports( this->M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr, model_measures_quantities_empty_t{} );
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
template <typename ModelFieldsType, typename SymbolsExpr, typename ModelMeasuresQuantitiesType>
void
CoefficientFormPDE<ConvexType,BasisUnknownType>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ModelMeasuresQuantitiesType const& mquantities )
{
    // execute common post process and save measures
    super_type::executePostProcessMeasures( time, this->mesh(), this->rangeMeshElements(), symbolsExpr, mfields, mquantities );
}

} // namespace Feel
} // namespace FeelModels

#include <feel/feelmodels/coefficientformpdes/coefficientformpdeassembly.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdeassemblystabilizationgls.hpp>

#endif
