/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDES_HPP 1

#include <feel/feelmodels/coefficientformpdes/coefficientformpde.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename... BasisUnknownType>
class CoefficientFormPDEs : public ModelNumerical,
                            public ModelGenericPDEs<ConvexType::nDim>,
                            public std::enable_shared_from_this< CoefficientFormPDEs<ConvexType,BasisUnknownType...> >
{
public :
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;

    using self_type = CoefficientFormPDEs<ConvexType,BasisUnknownType...>;
    using self_ptrtype = std::shared_ptr<self_type>;

    using coefficient_form_pde_base_type = CoefficientFormPDEBase<ConvexType>;
    using coefficient_form_pde_base_ptrtype = std::shared_ptr<coefficient_form_pde_base_type>;

    using mesh_type = typename coefficient_form_pde_base_type::mesh_type;
    using mesh_ptrtype = typename coefficient_form_pde_base_type::mesh_ptrtype;
    using convex_type = typename coefficient_form_pde_base_type::convex_type;
    static const uint16_type nDim = coefficient_form_pde_base_type::nDim;
    static const uint16_type nOrderGeo = coefficient_form_pde_base_type::nOrderGeo;

    // materials properties
    typedef MaterialsProperties<mesh_type::nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;


    static constexpr auto tuple_type_unknown_basis = hana::to_tuple(hana::tuple_t<BasisUnknownType...>);

private :
    struct internal_detail
    {
        template <typename FilterBasisUnknownType, int Index,typename ResType, typename... SomeType>
        static constexpr auto applyFilterBasisUnknownImpl( ResType && res, hana::tuple<SomeType...> const& t )
        {
            constexpr int nTupleElt = std::decay_t<decltype(hana::size( t ))>::value;
            if constexpr ( Index < nTupleElt )
            {
                auto const& currentBasis = hana::at( t, hana::int_c<Index> );
                //using current_basis_type = typename traits::template remove_hana_type_t< std::decay_t<decltype( currentBasis ) > >;
                using current_basis_type = typename std::decay_t<decltype( currentBasis ) >::type;
                if constexpr ( FilterBasisUnknownType::template apply<current_basis_type>::value )
                    return applyFilterBasisUnknownImpl<FilterBasisUnknownType,Index+1>( hana::append( std::forward<ResType>( res ), currentBasis ), t );
                else
                    return applyFilterBasisUnknownImpl<FilterBasisUnknownType,Index+1>( std::forward<ResType>( res ), t );
            }
            else
                return std::move( res );
        }
    };

    template <typename FilterBasisUnknownType>
    using tuple_type_unknown_basis_filtered_t = std::decay_t<decltype(internal_detail::template applyFilterBasisUnknownImpl<FilterBasisUnknownType,0>( hana::tuple<>{}, tuple_type_unknown_basis ) )>;

    template <typename FilterBasisUnknownType>
    static constexpr tuple_type_unknown_basis_filtered_t<FilterBasisUnknownType> tuple_type_unknown_basis_filtered = internal_detail::template applyFilterBasisUnknownImpl<FilterBasisUnknownType,0>( hana::tuple<>{}, tuple_type_unknown_basis );

public :

    struct traits
    {
        template <typename TheType>
        using remove_hana_type_t = typename std::decay_t<TheType>::type;

        template <typename TheBasisType>
        using coefficient_form_pde_t = CoefficientFormPDE<ConvexType,remove_hana_type_t<TheBasisType>>;

        template <typename... TheType>
        static constexpr auto
        variant_from_tuple( hana::tuple<TheType...> const& t )
            {
                return std::variant<TheType...>{};
            }

    private :
        struct TransformModelFields
        {
            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(t)>;
                    std::shared_ptr<coefficient_form_pde_type> dummycfpde;
                    return dummycfpde->modelFields( std::string("") );
                }

            struct View
            {
                template <typename T>
                constexpr auto operator()(T const& t) const
                    {
                        using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(t)>;
                        std::shared_ptr<coefficient_form_pde_type> dummycfpde;
                        return dummycfpde->modelFields( vector_ptrtype{}, size_type(0), std::string("") );
                    }
            };

            template <typename... TheType>
            static constexpr auto
            toModelFields(  hana::tuple<TheType...> const& t )
                {
                    return model_fields_t<TheType...>{};
                }
        };

        struct TransformTrialSelectorModelFields
        {
            template <typename T>
            constexpr auto operator()(T const& t) const
                {
                    using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(t)>;
                    std::shared_ptr<coefficient_form_pde_type> dummycfpde;
                    return dummycfpde->trialSelectorModelFields( 0 );
                }

            template <typename... TheType>
            static constexpr auto
            toTrialSelectorModelFields( hana::tuple<TheType...> const& t )
                {
                    return selector_model_fields_t<TheType...>{};
                }

        };


    public :
        template <typename FilterBasisUnknownType>
        using model_fields_type = std::decay_t<decltype( TransformModelFields::toModelFields( hana::transform( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, TransformModelFields{} ) ) )>;

        template <typename FilterBasisUnknownType>
        using model_fields_view_type = std::decay_t<decltype( TransformModelFields::toModelFields( hana::transform( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, typename TransformModelFields::View{} ) ) )>;

        template <typename FilterBasisUnknownType>
        using trial_selector_model_fields_type = std::decay_t<decltype( TransformTrialSelectorModelFields::toTrialSelectorModelFields( hana::transform( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, TransformTrialSelectorModelFields{} ) ) )>;
    };


    struct FilterBasisUnknownAll {
        template<typename T>
        struct apply { static constexpr bool value = true; };
    };
    template<typename TheBasisType>
    struct FilterBasisUnknown {
        template<typename T>
        struct apply { static constexpr bool value = std::is_same_v<T,TheBasisType>; };
    };

    using variant_unknown_basis_type = std::decay_t<decltype(traits::variant_from_tuple(tuple_type_unknown_basis))>;

    CoefficientFormPDEs( std::string const& prefix,
                         std::string const& keyword = "cfpdes",
                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                         std::string const& subPrefix  = "",
                         ModelBaseRepository const& modelRep = ModelBaseRepository() );

    static
    std::string const& unknowBasisTag( variant_unknown_basis_type const& vb );

    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();

    std::shared_ptr<std::ostringstream> getInfo() const override;
    void updateInformationObject( pt::ptree & p ) override;

    //! return all subtoolboxes related to each equation
    std::vector<std::shared_ptr<coefficient_form_pde_base_type>> const& coefficientFormPDEs() const { return M_coefficientFormPDEs; }

    //! return toolbox related to an equation from the name of the equation (empty ptr if not found)
    coefficient_form_pde_base_ptrtype coefficientFormPDE( std::string const& nameEq ) const
        {
            for (auto const& cfpdeBase : M_coefficientFormPDEs )
                if ( cfpdeBase->equationName() == nameEq )
                    return cfpdeBase;
            return std::shared_ptr<coefficient_form_pde_base_type>{};
        }

    //! return toolbox related to an equation from the base object and the fe basis (empty ptr if not found)
    template <typename TheBasisType>
    auto coefficientFormPDE( coefficient_form_pde_base_ptrtype const& cfpdeBase, TheBasisType const& e ) const
        {
            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                return std::shared_ptr<coefficient_form_pde_type>{};
            return std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
        }

    //! return toolbox related to an equation from the name of the equation and the fe basis (empty ptr if not found)
    template <typename TheBasisType>
    auto coefficientFormPDE( std::string const& nameEq, TheBasisType const& e ) const
        {
            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
            auto cfpdeBase = this->coefficientFormPDE( nameEq );
            if ( !cfpdeBase )
                return std::shared_ptr<coefficient_form_pde_type>{};
            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                return std::shared_ptr<coefficient_form_pde_type>{};
            return std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
        }

    //___________________________________________________________________________________//
    // mesh
    mesh_ptrtype const& mesh() const { return M_mesh; }
    void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }

    //___________________________________________________________________________________//
    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//
    // time stepping
    std::shared_ptr<TSBase> timeStepBase() const;
    void startTimeStep();
    void updateTimeStep();

    //___________________________________________________________________________________//
    // post-process
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename ModelFieldsType, typename SymbolsExpr>
    void executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr );

    bool checkResults() const override;

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

private :

    template <typename FilterBasisUnknownType,typename ResType>
    void modelFieldsImpl( std::string const& prefix, ResType && res ) const
        {
            hana::for_each( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, [this,&prefix,&res]( auto const& e )
                            {
                                for (auto const& cfpdeBase : M_coefficientFormPDEs )
                                {
                                    if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                        continue;

                                    using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                    auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                                    res = Feel::FeelModels::modelFields( res, cfpde->modelFields( prefixvm( prefix, cfpde->keyword() ) ) );
                                }
                            });
        }

    template <typename FilterBasisUnknownType,typename ResType>
    void modelFieldsImpl( vector_ptrtype sol, size_type rowStartInVector, std::string const& prefix, ResType && res ) const
        {
            hana::for_each( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, [this,&sol,&rowStartInVector,&prefix,&res]( auto const& e )
                            {
                                for (auto const& cfpdeBase : M_coefficientFormPDEs )
                                {
                                    if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                        continue;

                                    using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                    auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );
                                    res = Feel::FeelModels::modelFields( res, cfpde->modelFields( sol,
                                                                                                  rowStartInVector + this->startSubBlockSpaceIndex( cfpdeBase->physicDefault() ),
                                                                                                  prefixvm( prefix, cfpde->keyword() ) ) );
                                }
                            });
        }

    template <typename FilterBasisUnknownType,typename ResType>
    void trialSelectorModelFieldsImpl( size_type startBlockSpaceIndex, ResType && res ) const
        {
            hana::for_each( tuple_type_unknown_basis_filtered<FilterBasisUnknownType>, [this,&startBlockSpaceIndex,&res]( auto const& e )
                            {
                                for (auto const& cfpdeBase : M_coefficientFormPDEs )
                                {
                                    if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                        continue;

                                    using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                                    auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );

                                    res = Feel::FeelModels::selectorModelFields( res, cfpde->trialSelectorModelFields( startBlockSpaceIndex + this->startSubBlockSpaceIndex( cfpdeBase->physicDefault() ) ) );
                                }
                            });
        }

public :

    template <typename FilterBasisUnknownType = FilterBasisUnknownAll>
    auto modelFields( std::string const& prefix = "" ) const
        {
            typename traits::template model_fields_type<FilterBasisUnknownType> res;
            this->modelFieldsImpl<FilterBasisUnknownType>( prefix, std::move(res) );
            return res;
        }

    template <typename FilterBasisUnknownType = FilterBasisUnknownAll>
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            typename traits::template model_fields_view_type<FilterBasisUnknownType> res;
            this->modelFieldsImpl<FilterBasisUnknownType>( sol, rowStartInVector, prefix, std::move(res) );
            return res;
        }

    template <typename FilterBasisUnknownType = FilterBasisUnknownAll>
    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            typename traits::template trial_selector_model_fields_type<FilterBasisUnknownType> res;
            this->trialSelectorModelFieldsImpl<FilterBasisUnknownType>( startBlockSpaceIndex, std::move(res) );
            return res;
        }

    //___________________________________________________________________________________//
    // symbols expressions
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            //auto seHeat = this->heatModel()->symbolsExprToolbox( mfields );
            //auto seElectric = this->electricModel()->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            return Feel::vf::symbolsExpr( /*seHeat,seElectric,*/seParam,seMat,seFields );
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

    template <typename ModelFieldsType>
    auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelContext( mfields, this->symbolsExpr( mfields ) );
        }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            auto se = this->symbolsExpr( mfields );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }
   auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields );
            auto tse =  this->trialSymbolsExpr( mfields, trialSelectorModelFields( rowStartInVector ) );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ), std::move( tse ) );
        }
   auto modelContextNoTrialSymbolsExpr( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            auto se = this->symbolsExpr( mfields );
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }

    //___________________________________________________________________________________//
    // algebraic data and solver
    backend_ptrtype const& backend() const { return  M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    size_type nLocalDof() const;
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    //int nBlockMatrixGraph() const { return 1; }

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    void solve();

    void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename FilterBasisUnknownType,typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename FilterBasisUnknownType,typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const;

    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    template <typename FilterBasisUnknownType,typename ModelContextType>
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const;

    void updateJacobian( DataUpdateJacobian & data ) const override;
    template <typename FilterBasisUnknownType,typename ModelContextType>
    void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;

    void updateResidual( DataUpdateResidual & data ) const override;
    template <typename FilterBasisUnknownType,typename ModelContextType>
    void updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const;
    template <typename ModelContextType>
    void updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

private :
    void initMesh();
    void initMaterialProperties();
    void initPostProcess() override;

    void updateAutomaticSolverSelection();

    void updateTimeStepCurrentResidual();

    template <typename FilterBasisUnknownType>
    void updateLinearPDE_spec( DataUpdateLinear & data, std::any const& mctxAsAny ) const;
    template <typename FilterBasisUnknownType>
    void updateLinearPDEDofElimination_spec( DataUpdateLinear & data, std::any const& mctxAsAny ) const;
    template <typename FilterBasisUnknownType>
    void updateNewtonInitialGuess_spec( DataNewtonInitialGuess & data, std::any const& mctxAsAny ) const;
    template <typename FilterBasisUnknownType>
    void updateJacobian_spec( DataUpdateJacobian & data, std::any const& mctxAsAny ) const;
    template <typename FilterBasisUnknownType>
    void updateResidual_spec( DataUpdateResidual & data, std::any const& mctxAsAny ) const;
private :

    static const std::vector<std::string> S_unknownBasisTags;

    mesh_ptrtype M_mesh;

    // physical parameters
    materialsproperties_ptrtype M_materialsProperties;

    std::string M_solverName;

    // post-process
    export_ptrtype M_exporter;

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;

    vector_ptrtype M_timeStepThetaSchemePreviousContrib;

    std::map<std::string,double> M_currentParameterValues;

    std::vector<std::shared_ptr<coefficient_form_pde_base_type>> M_coefficientFormPDEs;

};


template< typename ConvexType, typename... BasisUnknownType>
template <typename ModelFieldsType,typename SymbolsExprType, typename ExportsExprType>
void
CoefficientFormPDEs<ConvexType,BasisUnknownType...>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExprType const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    if ( M_coefficientFormPDEs.empty() )
        return;
    this->log("CoefficientFormPDEs","exportResults", "start");
    this->timerTool("PostProcessing").start();

    //std::cout << "holalla \n "<< symbolsExpr.names() << std::endl; 

    hana::for_each( tuple_type_unknown_basis, [this,&time,&mfields,&symbolsExpr]( auto const& e )
                    {
                        for ( auto const& cfpdeBase : M_coefficientFormPDEs )
                        {
                            if ( this->unknowBasisTag( e ) != cfpdeBase->unknownBasis() )
                                continue;

                            using coefficient_form_pde_type = typename self_type::traits::template coefficient_form_pde_t<decltype(e)>;
                            auto cfpde = std::dynamic_pointer_cast<coefficient_form_pde_type>( cfpdeBase );

                            cfpde->exportResults( time,symbolsExpr );
                        }
                    });


    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );
    this->executePostProcessMeasures( time, mfields, symbolsExpr );
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_coefficientFormPDEs.front()->timeStepBase()->iteration(), mfields );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("CoefficientFormPDEs","exportResults", "finish");
}

template< typename ConvexType, typename... BasisUnknownType>
template <typename ModelFieldsType, typename SymbolsExpr>
void
CoefficientFormPDEs<ConvexType,BasisUnknownType...>::executePostProcessMeasures( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;
#if 0
    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
    //bool hasMeasurePoint = this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, mfields );
    if ( hasMeasureNorm || hasMeasureStatistics /*|| hasMeasurePoint*/ )
        hasMeasure = true;
#endif
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

#include <feel/feelmodels/coefficientformpdes/coefficientformpdesassembly.hpp>

#endif

