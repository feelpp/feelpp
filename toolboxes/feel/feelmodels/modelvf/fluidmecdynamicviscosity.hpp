#ifndef FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H
#define FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>
#include <feel/feelmodels/modelvf/exprevaluatefieldoperators.hpp>

#include <feel/feelmodels/fluid/fluidmechanicsmaterialproperties.hpp>

namespace Feel
{
namespace FeelModels
{

// Forward declaration
template<typename ExprEvaluateFieldOperatorsType, typename FiniteElementVelocityType, typename SymbolsExprType, ExprApplyType ExprApplied, ExprOperatorType ExprOp >
class FluidMecDynamicViscosityImpl;

template< typename ExprType>
struct FluidMecDynamicViscosityBase
{
    using expr_type = ExprType;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    using tensor_main_type = typename expr_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    using tensor_base_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                        typename expr_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>::shape,
                                        typename expr_type::value_type>;

    FluidMecDynamicViscosityBase( expr_type const& expr ):
        M_expr( expr )
    {}
    FluidMecDynamicViscosityBase( FluidMecDynamicViscosityBase const& ) = default;
    FluidMecDynamicViscosityBase( FluidMecDynamicViscosityBase && ) = default;

    expr_type const& expr() const { return M_expr; }

    virtual bool dependsOnVelocityField() const = 0;

    virtual size_type dynamicContext() const = 0;

    virtual uint16_type polynomialOrder() const = 0;

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename... TheArgsType>
    std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
    evaluator( tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, const TheArgsType&... theInitArgs/*Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu*/ ) const;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType,typename... TheArgsType>
    std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
    evaluator( std::true_type/**/, tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
               TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, const TheArgsType&... theInitArgs ) const;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
    void updateEvaluator( std::true_type /**/, std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> > & tensorToUpdate, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                          Geo_t const& geom, const TheArgsType&... theUpdateArgs );

private :
    expr_type const& M_expr;
};

/**
 * Newtonian
 */
template< typename ExprType>
class FluidMecDynamicViscosityNewtonian : public FluidMecDynamicViscosityBase<ExprType>
{
    using super_type = FluidMecDynamicViscosityBase<ExprType>;
public :
    using expr_type = typename super_type::expr_type;
    using this_type = FluidMecDynamicViscosityNewtonian<ExprType>;
    using material_property_scalar_expr_type = typename expr_type::template material_property_expr_type<1,1>;


    // special trick to avoid a infinite compilation loop.
    // it's caused by grad<expr_type::nDim>( material_property_scalar_expr_type{} ) because fluidmecdynaviscosity is alrady in se).
    // TODO : diff of symbol expr should be done before (when se is fully built)
    struct DeferToExprOpId
    {
        template<class Ts1>
        using execute = Ts1;
    };
    struct DeferToExprOpGrad
    {
        template<class Ts1>
        using execute = std::decay_t<decltype( grad<expr_type::nDim>( Ts1{} ) )>;
    };
    template< class Prog, class... Ts >
    using runRRR = typename Prog::template execute<Ts...>;
    //static const bool is_same_expr = expr_type::isExprOpID();//std::is_same_v<the_expr_type,the_expr_expand_type>;
    using choice = typename std::conditional< expr_type::isExprOpID(), DeferToExprOpId, DeferToExprOpGrad  >::type;

    using expr_dynamic_viscosity_newtonian_type = runRRR< choice, material_property_scalar_expr_type >;


#if 0
    using expr_grad_material_property_scalar_type = std::decay_t<decltype( grad<expr_type::nDim>( material_property_scalar_expr_type{} ) )>;
    using expr_dynamic_viscosity_newtonian_type = typename mpl::if_c< expr_type::isExprOpID(),
                                                                      material_property_scalar_expr_type, //material_property_scalar_expr_type
                                                                      expr_grad_material_property_scalar_type >::type;
#endif
    template <ExprOperatorType TheExprOp = expr_type::exprOp, std::enable_if_t< TheExprOp == ExprOperatorType::ID, bool> = true>
    FluidMecDynamicViscosityNewtonian( ExprType const& expr ):
        super_type( expr ),
        M_exprDynamicViscosity( expr.template materialPropertyExpr<1,1>("dynamic-viscosity") )
    {}

    template <ExprOperatorType TheExprOp = expr_type::exprOp, std::enable_if_t< TheExprOp == ExprOperatorType::GRAD, bool> = true>
    FluidMecDynamicViscosityNewtonian( ExprType const& expr ):
        super_type( expr ),
        M_exprDynamicViscosity( grad<expr_type::nDim>( expr.template materialPropertyExpr<1,1>("dynamic-viscosity")/*, "", world, dirLibExpr*/ /*TODO*/ ) )
    {}

    FluidMecDynamicViscosityNewtonian( FluidMecDynamicViscosityNewtonian const& ) = default;
    FluidMecDynamicViscosityNewtonian( FluidMecDynamicViscosityNewtonian && ) = default;

    bool dependsOnVelocityField() const override { return false; }

    size_type dynamicContext() const override { return Feel::vf::dynamicContext( M_exprDynamicViscosity ); }

    uint16_type polynomialOrder() const override { return M_exprDynamicViscosity.polynomialOrder(); }

    expr_dynamic_viscosity_newtonian_type const& exprDynamicViscosity() const { return M_exprDynamicViscosity; }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return M_exprDynamicViscosity.hasSymbolDependency( symb, se );
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicViscosity().evaluator( geom ) )
        {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev ):
            super_type( geom,fev ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicViscosity().evaluator( geom ) )
        {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom ):
            super_type( geom ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicViscosity().evaluator( geom ) )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs ):
            super_type( geom ),
            M_expr( expr ),
            M_muExprTensor( std::true_type{}, exprExpanded.exprDynamicViscosity(), ttse, expr.exprDynamicViscosity(), geom, theInitArgs... )
        {}


        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
        {
            M_muExprTensor.update( geom );
            this->updateImpl();
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_muExprTensor.update( std::true_type{}, exprExpanded.exprDynamicViscosity(), ttse, geom, theUpdateArgs... );
            this->updateImpl();
        }

        void updateImpl()
        {
            uint16_type nPoints = this->gmc()->nPoints();
            if ( M_localEval.size() != nPoints )
                M_localEval.setConstant( nPoints, super_type::matrix_shape_type::Zero() );

            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                if constexpr ( expr_type::isExprOpID() )
                    M_localEval[q](0,0) = M_muExprTensor.evalq(0,0,q);
                else
                {
                    for ( uint16_type c1 = 0; c1 < expr_type::nDim; ++c1 )
                        M_localEval[q](0,c1) = M_muExprTensor.evalq(0,c1,q);
                }
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                return value_type(0);
            }
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                this->locMatrixShape().setZero();
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "not allow";
                return value_type(0);
            }
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                CHECK( false) << "not allow";
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
            evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
            {
                return M_localEval[q](c1,c2);
            }
        ret_type
            evalq( uint16_type q ) const override
            {
                return ret_type(M_localEval[q].data());
            }


    private :
        this_type const& M_expr;
        typename expr_dynamic_viscosity_newtonian_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/>  M_muExprTensor;
        array_shape_type M_localEval;
    };


private:
    expr_dynamic_viscosity_newtonian_type M_exprDynamicViscosity;
};


/**
 * Power Law
 */
template< typename ExprType>
class FluidMecDynamicViscosityPowerLaw : public FluidMecDynamicViscosityBase<ExprType>
{
    using super_type = FluidMecDynamicViscosityBase<ExprType>;
public :
    using expr_type = typename super_type::expr_type;
    using this_type = FluidMecDynamicViscosityPowerLaw<ExprType>;
    using material_property_scalar_expr_type = typename expr_type::template material_property_expr_type<1,1>;

    FluidMecDynamicViscosityPowerLaw( ExprType const& expr ):
        super_type( expr ),
        M_kExpr( expr.template materialPropertyExpr<1,1>("consistency-index") ),
        M_nExpr( expr.template materialPropertyExpr<1,1>("power-law-index") ),
        M_muMinExpr( expr.template materialPropertyExpr<1,1>("viscosity-min") ),
        M_muMaxExpr( expr.template materialPropertyExpr<1,1>("viscosity-max") )
    {
        expr.exprEvaluateVelocityOperatorsPtr()->setEnableGrad( true );
    }

    FluidMecDynamicViscosityPowerLaw( FluidMecDynamicViscosityPowerLaw const& ) = default;
    FluidMecDynamicViscosityPowerLaw( FluidMecDynamicViscosityPowerLaw && ) = default;

    bool dependsOnVelocityField() const override { return true; }

    size_type dynamicContext() const override
    {
        return Feel::vf::dynamicContext( M_kExpr ) | Feel::vf::dynamicContext( M_nExpr ) | Feel::vf::dynamicContext( M_muMinExpr ) | Feel::vf::dynamicContext( M_muMaxExpr );
    }

    uint16_type polynomialOrder() const override { return 2+this->expr().exprEvaluateVelocityOperatorsPtr()->exprGrad().polynomialOrder(); }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return (symb == "x") || (symb == "y") || (symb == "z") ||
            M_kExpr.hasSymbolDependency( symb, se ) || M_nExpr.hasSymbolDependency( symb, se ) ||
            M_muMinExpr.hasSymbolDependency( symb, se )  || M_muMaxExpr.hasSymbolDependency( symb, se ) ;
    }

    material_property_scalar_expr_type const& kExpr() const { return M_kExpr; }
    material_property_scalar_expr_type const& nExpr() const { return M_nExpr; }
    material_property_scalar_expr_type const& muMinExpr() const { return M_muMinExpr; }
    material_property_scalar_expr_type const& muMaxExpr() const { return M_muMaxExpr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        using tensor_main_type = typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t>;
        using tensor_expr_evaluate_velocity_opertors_type = typename tensor_main_type::tensor_expr_evaluate_velocity_opertors_type;


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( this->expr().kExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_muMinExprTensor( this->expr().muMinExpr().evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().muMaxExpr().evaluator( geom ) )
        {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev ):
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( this->expr().kExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_muMinExprTensor( this->expr().muMinExpr().evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().muMaxExpr().evaluator( geom ) )
        {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom ):
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( this->expr().kExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_muMinExprTensor( this->expr().muMinExpr().evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().muMaxExpr().evaluator( geom ) )
        {}
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs ):
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( std::true_type{}, exprExpanded.kExpr(), ttse, expr.kExpr(), geom, theInitArgs... ),
            M_nExprTensor( std::true_type{}, exprExpanded.nExpr(), ttse, expr.nExpr(), geom, theInitArgs... ),
            M_muMinExprTensor( std::true_type{}, exprExpanded.muMinExpr(), ttse, expr.muMinExpr(), geom, theInitArgs... ),
            M_muMaxExprTensor( std::true_type{}, exprExpanded.muMaxExpr(), ttse, expr.muMaxExpr(), geom, theInitArgs... )
        {}

        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
        {
            this->setGmc( geom );
            M_kExprTensor.update( geom );
            M_nExprTensor.update( geom );
            M_muMinExprTensor.update( geom );
            M_muMaxExprTensor.update( geom );
            this->updateImpl();
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            this->setGmc( geom );
            M_kExprTensor.update( std::true_type{}, exprExpanded.kExpr(), ttse, geom, theUpdateArgs... );
            M_nExprTensor.update( std::true_type{}, exprExpanded.nExpr(), ttse, geom, theUpdateArgs... );
            M_muMinExprTensor.update( std::true_type{}, exprExpanded.muMinExpr(), ttse, geom, theUpdateArgs... );
            M_muMaxExprTensor.update( std::true_type{}, exprExpanded.muMaxExpr(), ttse, geom, theUpdateArgs... );
            this->updateImpl();
        }


        void updateImpl()
        {
            uint16_type nPoints = this->gmc()->nPoints();
            if ( M_localEval.size() != nPoints )
            {
                M_localEval.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
                M_muIsInInterval.resize( nPoints, false );
                if constexpr ( expr_type::is_applied_as_jacobian )
                    M_localEvalPrecompute.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
            }

            value_type gammapoint2v(0);
            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                const value_type power_k_generic = M_kExprTensor.evalq(0,0,q);
                const value_type power_n_generic = ( M_nExprTensor.evalq(0,0,q) - 1.)/2.;
                const value_type muMin = M_muMinExprTensor.evalq(0,0,q);
                const value_type muMax = M_muMaxExprTensor.evalq(0,0,q);

                //auto const mu_powerlaw = power_k_generic*pow( 2.0*inner(defv,defv) /*+chiInv*/ , cst( power_n_generic ) )/**chiSup*/;
                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );

                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                    gammapoint2v = 2*DxD;
                }
                else
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    const value_type DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                        0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                    gammapoint2v = 2*DxD;
                }

                value_type muEval = power_k_generic*math::pow( gammapoint2v , power_n_generic );

                if ( muEval < muMin )
                {
                    muEval = muMin;
                    M_muIsInInterval[q] = false;
                }
                else if ( muEval > muMax )
                {
                    muEval = muMax;
                    M_muIsInInterval[q] = false;
                }
                else
                    M_muIsInInterval[q] = true;

                M_localEval[q](0,0) = muEval;

                if constexpr ( expr_type::is_applied_as_jacobian )
                    M_localEvalPrecompute[q](0,0) = power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "TODO";
                return value_type(0);
            }
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                matrix_shape_type & locMat = this->locMatrixShape();
                if ( M_muIsInInterval[q] )
                {
                    auto const& gradTrial = this->fecTrial()->grad( j, q );
                    auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                    value_type gammapoint2t;
                    if constexpr ( gmc_type::nDim == 2 )
                    {
                        const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
                        const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                        gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
                    }
                    else
                    {
                        const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
                        const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
                        const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                        const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                        const value_type DxDt01 = du1tdy+du2tdx;
                        const value_type DxDv01 = du1vdy+du2vdx;
                        const value_type DxDt02 = du1tdz+du3tdx;
                        const value_type DxDv02 = du1vdz+du3vdx;
                        const value_type DxDt12 = du2tdz+du3tdy;
                        const value_type DxDv12 = du2vdz+du3vdy;
                        gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                            DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
                    }
                    //const value_type muEval = M_localEval[q](0,0);
                    //const value_type mut = gammapoint2t*power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                    const value_type mut = gammapoint2t*M_localEvalPrecompute[q](0,0);
                    locMat(0,0) = mut;
                }
                else
                {
                    locMat(0,0) = 0.;
                }

                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "TODO";
                return value_type(0);
            }
        }

        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                CHECK( false) << "TODO";
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            return M_localEval[q](c1,c2);
        }

        ret_type
        evalq( uint16_type q ) const override
        {
            return ret_type(M_localEval[q].data());
        }


    private :
        this_type const& M_expr;
        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        typename material_property_scalar_expr_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/> M_kExprTensor, M_nExprTensor, M_muMinExprTensor, M_muMaxExprTensor;
        std::vector<bool> M_muIsInInterval;
        array_shape_type M_localEval, M_localEvalPrecompute;
    };

private :
    material_property_scalar_expr_type  M_kExpr, M_nExpr, M_muMinExpr, M_muMaxExpr;
};


template< typename ExprType>
class FluidMecDynamicViscosityCarreau : public FluidMecDynamicViscosityBase<ExprType>
{
    using super_type = FluidMecDynamicViscosityBase<ExprType>;
public :
    using expr_type = typename super_type::expr_type;
    using this_type = FluidMecDynamicViscosityCarreau<ExprType>;
    using material_property_scalar_expr_type = typename expr_type::template material_property_expr_type<1,1>;

    FluidMecDynamicViscosityCarreau( ExprType const& expr ):
        super_type( expr ),
        M_mu0Expr( expr.template materialPropertyExpr<1,1>("viscosity-zero-shear") ),
        M_muInfExpr( expr.template materialPropertyExpr<1,1>("viscosity-infinite-shear") ),
        M_lambdaExpr( expr.template materialPropertyExpr<1,1>("carreau-law-lambda") ),
        M_nExpr( expr.template materialPropertyExpr<1,1>("carreau-law-n") )
    {
        expr.exprEvaluateVelocityOperatorsPtr()->setEnableGrad( true );
    }

    FluidMecDynamicViscosityCarreau( FluidMecDynamicViscosityCarreau const& ) = default;
    FluidMecDynamicViscosityCarreau( FluidMecDynamicViscosityCarreau && ) = default;

    bool dependsOnVelocityField() const override { return true; }

    size_type dynamicContext() const override
    {
        return Feel::vf::dynamicContext( M_mu0Expr ) | Feel::vf::dynamicContext( M_muInfExpr ) | Feel::vf::dynamicContext( M_lambdaExpr ) | Feel::vf::dynamicContext( M_nExpr );
    }

    uint16_type polynomialOrder() const override { return 2+this->expr().exprEvaluateVelocityOperatorsPtr()->exprGrad().polynomialOrder(); }


    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return (symb == "x") || (symb == "y") || (symb == "z") ||
            M_mu0Expr.hasSymbolDependency( symb, se ) || M_muInfExpr.hasSymbolDependency( symb, se ) ||
            M_lambdaExpr.hasSymbolDependency( symb, se ) || M_nExpr.hasSymbolDependency( symb, se ) ;
    }

    material_property_scalar_expr_type const& mu0Expr() const { return M_mu0Expr; }
    material_property_scalar_expr_type const& muInfExpr() const { return M_muInfExpr; }
    material_property_scalar_expr_type const& lambdaExpr() const { return M_lambdaExpr; }
    material_property_scalar_expr_type const& nExpr() const { return M_nExpr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        using tensor_main_type = typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t>;
        using tensor_expr_evaluate_velocity_opertors_type = typename tensor_main_type::tensor_expr_evaluate_velocity_opertors_type;


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) )
        {}

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev ):
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) )
        {}

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom ):
        super_type( geom ),
        M_expr( expr ),
        M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
        M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
        M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
        M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
        M_nExprTensor( this->expr().nExpr().evaluator( geom ) )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs ):
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( std::true_type{}, exprExpanded.mu0Expr(), ttse, expr.mu0Expr(), geom, theInitArgs... ),
            M_muInfExprTensor( std::true_type{}, exprExpanded.muInfExpr(), ttse, expr.muInfExpr(), geom, theInitArgs... ),
            M_lambdaExprTensor( std::true_type{}, exprExpanded.lambdaExpr(), ttse, expr.lambdaExpr(), geom, theInitArgs... ),
            M_nExprTensor( std::true_type{}, exprExpanded.nExpr(), ttse, expr.nExpr(), geom, theInitArgs... )
        {}

        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
        {
            this->setGmc( geom );
            M_mu0ExprTensor.update( geom );
            M_muInfExprTensor.update( geom );
            M_lambdaExprTensor.update( geom );
            M_nExprTensor.update( geom );
            this->updateImpl();
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            this->setGmc( geom );
            M_mu0ExprTensor.update( std::true_type{}, exprExpanded.mu0Expr(), ttse, geom, theUpdateArgs... );
            M_muInfExprTensor.update( std::true_type{}, exprExpanded.muInfExpr(), ttse, geom, theUpdateArgs... );
            M_lambdaExprTensor.update( std::true_type{}, exprExpanded.lambdaExpr(), ttse, geom, theUpdateArgs... );
            M_nExprTensor.update( std::true_type{}, exprExpanded.nExpr(), ttse, geom, theUpdateArgs... );
            this->updateImpl();
        }

        void updateImpl()
        {
            uint16_type nPoints = this->gmc()->nPoints();
            if ( M_localEval.size() != nPoints )
            {
                M_localEval.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
                if constexpr ( expr_type::is_applied_as_jacobian )
                     M_localEvalPrecompute.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
            }

            value_type DxD(0);
            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreau_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreau_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreau_lambda_pow2_time2 = math::pow(carreau_lambda,2.)*2.0;
                const value_type carreauValPower =  (carreau_n-1)/2.0;
                const value_type carreau_lambda2 = math::pow(carreau_lambda,2);
                const value_type part1_carreauLaw = ( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2;

                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );

                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                }
                else
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                        0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                }
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreau_lambda_pow2_time2*DxD, carreauValPower );
                M_localEval[q](0,0) = muEval;

                if constexpr ( expr_type::is_applied_as_jacobian )
                {
                    const value_type gammapoint2v = 2*DxD;
                    M_localEvalPrecompute[q](0,0) = part1_carreauLaw*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                }
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "TODO";
                return value_type(0);
            }
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                matrix_shape_type & locMat = this->locMatrixShape();
                auto const& gradTrial = this->fecTrial()->grad( j, q );
                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                value_type gammapoint2t;
                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
                }
                else
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
                    const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    const value_type DxDt01 = du1tdy+du2tdx;
                    const value_type DxDv01 = du1vdy+du2vdx;
                    const value_type DxDt02 = du1tdz+du3tdx;
                    const value_type DxDv02 = du1vdz+du3vdx;
                    const value_type DxDt12 = du2tdz+du3tdy;
                    const value_type DxDv12 = du2vdz+du3vdy;
                    gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) + DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
                }
                const value_type mut = gammapoint2t*M_localEvalPrecompute[q](0,0);
                locMat(0,0) = mut;

                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false ) << "not allow";
                return value_type(0);
            }
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                CHECK( false ) << "not allow";
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            return M_localEval[q](c1,c2);
        }

        ret_type
        evalq( uint16_type q ) const override
        {
            return ret_type(M_localEval[q].data());
        }


    private :
        this_type const& M_expr;
        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        typename material_property_scalar_expr_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/> M_mu0ExprTensor, M_muInfExprTensor, M_lambdaExprTensor, M_nExprTensor;
        array_shape_type M_localEval, M_localEvalPrecompute;
    };

private :
    material_property_scalar_expr_type M_mu0Expr, M_muInfExpr, M_lambdaExpr, M_nExpr;
};

template< typename ExprType>
class FluidMecDynamicViscosityCarreauYasuda : public FluidMecDynamicViscosityBase<ExprType>
{
    using super_type = FluidMecDynamicViscosityBase<ExprType>;
public :
    using expr_type = typename super_type::expr_type;
    using this_type = FluidMecDynamicViscosityCarreauYasuda<ExprType>;
    using material_property_scalar_expr_type = typename expr_type::template material_property_expr_type<1,1>;

    FluidMecDynamicViscosityCarreauYasuda( ExprType const& expr ):
        super_type( expr ),
        M_mu0Expr( expr.template materialPropertyExpr<1,1>("viscosity-zero-shear") ),
        M_muInfExpr( expr.template materialPropertyExpr<1,1>("viscosity-infinite-shear") ),
        M_lambdaExpr( expr.template materialPropertyExpr<1,1>("carreau-yasuda-law-lambda") ),
        M_nExpr( expr.template materialPropertyExpr<1,1>("carreau-yasuda-law-n") ),
        M_aExpr( expr.template materialPropertyExpr<1,1>("carreau-yasuda-law-a") )
    {
        expr.exprEvaluateVelocityOperatorsPtr()->setEnableGrad( true );
    }

    FluidMecDynamicViscosityCarreauYasuda( FluidMecDynamicViscosityCarreauYasuda const& ) = default;
    FluidMecDynamicViscosityCarreauYasuda( FluidMecDynamicViscosityCarreauYasuda && ) = default;

    bool dependsOnVelocityField() const override { return true; }

    size_type dynamicContext() const override
    {
        return Feel::vf::dynamicContext( M_mu0Expr ) | Feel::vf::dynamicContext( M_muInfExpr ) | Feel::vf::dynamicContext( M_lambdaExpr ) | Feel::vf::dynamicContext( M_nExpr ) | Feel::vf::dynamicContext( M_aExpr );
    }

    uint16_type polynomialOrder() const override { return 2+this->expr().exprEvaluateVelocityOperatorsPtr()->exprGrad().polynomialOrder(); }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return (symb == "x") || (symb == "y") || (symb == "z") ||
            M_mu0Expr.hasSymbolDependency( symb, se ) || M_muInfExpr.hasSymbolDependency( symb, se ) ||
            M_lambdaExpr.hasSymbolDependency( symb, se ) || M_nExpr.hasSymbolDependency( symb, se ) || M_aExpr.hasSymbolDependency( symb, se );
    }

    material_property_scalar_expr_type const& mu0Expr() const { return M_mu0Expr; }
    material_property_scalar_expr_type const& muInfExpr() const { return M_muInfExpr; }
    material_property_scalar_expr_type const& lambdaExpr() const { return M_lambdaExpr; }
    material_property_scalar_expr_type const& nExpr() const { return M_nExpr; }
    material_property_scalar_expr_type const& aExpr() const { return M_aExpr; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        using tensor_main_type = typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t>;
        using tensor_expr_evaluate_velocity_opertors_type = typename tensor_main_type::tensor_expr_evaluate_velocity_opertors_type;


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_aExprTensor( this->expr().aExpr().evaluator( geom ) )
        {}

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev ):
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_aExprTensor( this->expr().aExpr().evaluator( geom ) )
        {}

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom ):
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_aExprTensor( this->expr().aExpr().evaluator( geom ) )
        {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs ):
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( std::true_type{}, exprExpanded.mu0Expr(), ttse, expr.mu0Expr(), geom, theInitArgs... ),
            M_muInfExprTensor( std::true_type{}, exprExpanded.muInfExpr(), ttse, expr.muInfExpr(), geom, theInitArgs... ),
            M_lambdaExprTensor( std::true_type{}, exprExpanded.lambdaExpr(), ttse, expr.lambdaExpr(), geom, theInitArgs... ),
            M_nExprTensor( std::true_type{}, exprExpanded.nExpr(), ttse, expr.nExpr(), geom, theInitArgs... ),
            M_aExprTensor( std::true_type{}, exprExpanded.aExpr(), ttse, expr.aExpr(), geom, theInitArgs... )
        {}

        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
        {
            this->setGmc( geom );
            M_mu0ExprTensor.update( geom );
            M_muInfExprTensor.update( geom );
            M_lambdaExprTensor.update( geom );
            M_nExprTensor.update( geom );
            M_aExprTensor.update( geom );
            this->updateImpl();
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            this->setGmc( geom );
            M_mu0ExprTensor.update( std::true_type{}, exprExpanded.mu0Expr(), ttse, geom, theUpdateArgs... );
            M_muInfExprTensor.update( std::true_type{}, exprExpanded.muInfExpr(), ttse, geom, theUpdateArgs... );
            M_lambdaExprTensor.update( std::true_type{}, exprExpanded.lambdaExpr(), ttse, geom, theUpdateArgs... );
            M_nExprTensor.update( std::true_type{}, exprExpanded.nExpr(), ttse, geom, theUpdateArgs... );
            M_aExprTensor.update( std::true_type{}, exprExpanded.aExpr(), ttse, geom, theUpdateArgs... );
            this->updateImpl();
        }

        void updateImpl()
        {
            uint16_type nPoints = this->gmc()->nPoints();
            if ( M_localEval.size() != nPoints )
            {
                M_localEval.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
                if constexpr ( expr_type::is_applied_as_jacobian )
                     M_localEvalPrecompute.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
            }

            value_type gammapoint2v(0);
            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_a = M_aExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda_pow_a = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type carreauYasudaValPower = carreauYasuda_a/2.;
                const value_type carreauYasudaValPower2 = (carreauYasuda_n-1)/carreauYasuda_a;
                const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);

                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );

                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                    gammapoint2v = 2*DxD;
                }
                else
                {
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    const value_type DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                        0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                    gammapoint2v = 2*DxD;
                }
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreauYasuda_lambda_pow_a*math::pow( gammapoint2v, carreauYasudaValPower) , carreauYasudaValPower2 );
                M_localEval[q](0,0) = muEval;

                if constexpr ( expr_type::is_applied_as_jacobian )
                {
                    const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
                    const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
                    M_localEvalPrecompute[q](0,0) = part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                }
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "TODO";
                return value_type(0);
            }
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                matrix_shape_type & locMat = this->locMatrixShape();

                auto const& gradTrial = this->fecTrial()->grad( j, q );
                auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                value_type gammapoint2t;
                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
                }
                else
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
                    const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    const value_type DxDt01 = du1tdy+du2tdx;
                    const value_type DxDv01 = du1vdy+du2vdx;
                    const value_type DxDt02 = du1tdz+du3tdx;
                    const value_type DxDv02 = du1vdz+du3vdx;
                    const value_type DxDt12 = du2tdz+du3tdy;
                    const value_type DxDv12 = du2vdz+du3vdy;
                    gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                        DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
                }
                const value_type mut = gammapoint2t*M_localEvalPrecompute[q](0,0);
                locMat(0,0) = mut;

                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false ) << "not allow";
                return value_type(0);
            }
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                CHECK( false ) << "not allow";
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            return M_localEval[q](c1,c2);
        }
        ret_type
        evalq( uint16_type q ) const override
        {
            return ret_type(M_localEval[q].data());
        }


    private :
        this_type const& M_expr;
        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        typename material_property_scalar_expr_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/> M_mu0ExprTensor, M_muInfExprTensor, M_lambdaExprTensor, M_nExprTensor, M_aExprTensor;
        array_shape_type M_localEval, M_localEvalPrecompute;
    };

private :
    material_property_scalar_expr_type M_mu0Expr, M_muInfExpr, M_lambdaExpr, M_nExpr, M_aExpr;
};

template< typename ExprType >
class FluidMecDynamicViscosityMultifluid : public FluidMecDynamicViscosityBase<ExprType>
{
    using super_type = FluidMecDynamicViscosityBase<ExprType>;
public :
    using expr_type = typename super_type::expr_type;
    using this_type = FluidMecDynamicViscosityMultifluid<ExprType>;

    using material_property_scalar_expr_type = typename expr_type::template material_property_expr_type<1,1>;

    using fluidmec_dynamicviscosity_impl_type = FluidMecDynamicViscosityImpl<
        typename expr_type::expr_evaluate_velocity_opertors_type,
        typename expr_type::finite_element_velocity_type,
        typename expr_type::symbols_expr_type,
        expr_type::expr_apply_type,
        expr_type::expr_operator_type
            >;
    using fluidmec_dynamicviscosity_impl_ptrtype = std::shared_ptr<fluidmec_dynamicviscosity_impl_type>;

    FluidMecDynamicViscosityMultifluid( ExprType const& ex ):
        super_type( ex ),
        M_fluidMecOuterDynamicViscosityImpl(),
        M_fluidMecInnerDynamicViscosityImpls(),
        M_exprHeavisides()
    {
        auto const & matProps = ex.materialProperties();
        auto const multifluidDynamicViscosityLawPtr = matProps.template law<DynamicViscosityLaw>( "dynamic-viscosity" );
        auto const & multifluidLaw = multifluidDynamicViscosityLawPtr->template law<DynamicViscosityLaw::MultifluidLaw>();
        // outer dynamic viscosity
        auto const & outerSubMatProps = ex.materialProperties().subMaterialProperties( multifluidLaw.outerFluid );
        auto const outerDynamicViscosityLawPtr = outerSubMatProps.template law<DynamicViscosityLaw>( "dynamic-viscosity" );
        M_fluidMecOuterDynamicViscosityImpl.reset( new fluidmec_dynamicviscosity_impl_type( ex.exprEvaluateVelocityOperatorsPtr(), *outerDynamicViscosityLawPtr, outerSubMatProps, ex.polynomialOrder(), ex.symbolsExpr() ) );

        // inner dynamic viscosities and heaviside functions
        for( std::string const & innerFluid: multifluidLaw.innerFluids )
        {
            auto const & innerSubMatProps = ex.materialProperties().subMaterialProperties( innerFluid );
            auto const innerDynamicViscosityLawPtr = innerSubMatProps.template law<DynamicViscosityLaw>( "dynamic-viscosity" );
            M_fluidMecInnerDynamicViscosityImpls.emplace( innerFluid, new fluidmec_dynamicviscosity_impl_type( ex.exprEvaluateVelocityOperatorsPtr(), *innerDynamicViscosityLawPtr, innerSubMatProps, ex.polynomialOrder(), ex.symbolsExpr() ) );

            M_exprHeavisides.emplace( innerFluid, expr( innerSubMatProps.property( "heaviside" ).template expr<1,1>(), ex.symbolsExpr() ) );
        }
    }

    FluidMecDynamicViscosityMultifluid( FluidMecDynamicViscosityMultifluid const& ) = default;
    FluidMecDynamicViscosityMultifluid( FluidMecDynamicViscosityMultifluid && ) = default;

    fluidmec_dynamicviscosity_impl_type const& exprOuterDynamicViscosity() const { return *M_fluidMecOuterDynamicViscosityImpl; }
    fluidmec_dynamicviscosity_impl_type const& exprInnerDynamicViscosity( std::string const& innerMat ) const { return *(M_fluidMecInnerDynamicViscosityImpls.at( innerMat )); }
    material_property_scalar_expr_type const& exprHeaviside( std::string const& innerMat ) const { return M_exprHeavisides.at( innerMat ); }
    std::map<std::string, fluidmec_dynamicviscosity_impl_ptrtype> const& exprInnerDynamicViscosities() const { return M_fluidMecInnerDynamicViscosityImpls; }
    std::map<std::string, material_property_scalar_expr_type> const& exprHeavisides() const { return M_exprHeavisides; }

    bool dependsOnVelocityField() const override 
    { 
        return M_fluidMecOuterDynamicViscosityImpl->dependsOnVelocityField() || std::any_of( M_fluidMecInnerDynamicViscosityImpls.begin(), M_fluidMecInnerDynamicViscosityImpls.end(), []( auto const & innerDynamicViscosityImplPair ) { return innerDynamicViscosityImplPair.second->dependsOnVelocityField(); } ); 
    }

    size_type dynamicContext() const override 
    { 
        size_type context = M_fluidMecOuterDynamicViscosityImpl->dynamicContext();
        for( auto const& [name, innerDynamicViscosityImpl]: M_fluidMecInnerDynamicViscosityImpls )
            context |= innerDynamicViscosityImpl->dynamicContext();
        for( auto const& [name, exprHeaviside]: M_exprHeavisides )
            context |= Feel::vf::dynamicContext( exprHeaviside );
        return context;
    }

    uint16_type polynomialOrder() const override 
    { 
        auto polynomialOrderPtrCmp = []( auto const& lhs, auto const& rhs ) 
        { 
            return lhs.second->polynomialOrder() < rhs.second->polynomialOrder(); 
        };
        auto polynomialOrderCmp = []( auto const& lhs, auto const& rhs ) 
        { 
            return lhs.second.polynomialOrder() < rhs.second.polynomialOrder(); 
        };

        return std::max( 
                M_fluidMecOuterDynamicViscosityImpl->polynomialOrder(), 
                std::max_element( M_fluidMecInnerDynamicViscosityImpls.begin(), M_fluidMecInnerDynamicViscosityImpls.end(), polynomialOrderPtrCmp )->second->polynomialOrder() )
            + std::max_element( M_exprHeavisides.begin(), M_exprHeavisides.end(), polynomialOrderCmp )->second.polynomialOrder();
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        bool hasDep = M_fluidMecOuterDynamicViscosityImpl->hasSymbolDependency( symb, se );
        for( auto const& [name, innerDynamicViscosityImpl]: M_fluidMecInnerDynamicViscosityImpls )
            hasDep |= innerDynamicViscosityImpl->hasSymbolDependency( symb, se );
        for( auto const& [name, exprHeaviside]: M_exprHeavisides )
            hasDep |= exprHeaviside.hasSymbolDependency( symb, se );
        return hasDep;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>
    {
        using super_type = typename this_type::super_type::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t>;
    public :
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;

        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        using mu_expr_tensor_type = typename fluidmec_dynamicviscosity_impl_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using heaviside_expr_tensor_type = typename material_property_scalar_expr_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/>;

        struct inner_expr_tensors_type
        {
            mu_expr_tensor_type muExprTensor;
            heaviside_expr_tensor_type heavisideExprTensor;
        };

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_muOuterExprTensor( this->expr().exprOuterDynamicViscosity(), geom )
        {
            for( auto const & [innerMatName, innerDynamicViscosity]: this->expr().exprInnerDynamicViscosities() )
            {
                M_innerExprTensors.emplace( innerMatName, inner_expr_tensors_type{ 
                        mu_expr_tensor_type( *innerDynamicViscosity, geom ),
                        heaviside_expr_tensor_type( this->expr().exprHeaviside( innerMatName ).evaluator( geom ) )
                        } );
            }
        }
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev ):
            super_type( geom,fev ),
            M_expr( expr ),
            M_muOuterExprTensor( this->expr().exprOuterDynamicViscosity(), geom )
        {
            for( auto const & [innerMatName, innerDynamicViscosity]: this->expr().exprInnerDynamicViscosities() )
            {
                M_innerExprTensors.emplace( innerMatName, inner_expr_tensors_type{ 
                        mu_expr_tensor_type( *innerDynamicViscosity, geom ),
                        heaviside_expr_tensor_type( this->expr().exprHeaviside( innerMatName ).evaluator( geom ) )
                        } );
            }
        }
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom ):
            super_type( geom ),
            M_expr( expr ),
            M_muOuterExprTensor( this->expr().exprOuterDynamicViscosity(), geom )
        {
            for( auto const & [innerMatName, innerDynamicViscosity]: this->expr().exprInnerDynamicViscosities() )
            {
                M_innerExprTensors.emplace( innerMatName, inner_expr_tensors_type{ 
                        mu_expr_tensor_type( *innerDynamicViscosity, geom ),
                        heaviside_expr_tensor_type( this->expr().exprHeaviside( innerMatName ).evaluator( geom ) )
                        } );
            }
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs ):
            super_type( geom ),
            M_expr( expr ),
            M_muOuterExprTensor( std::true_type{}, exprExpanded.exprOuterDynamicViscosity(), ttse, expr.exprOuterDynamicViscosity(), geom, theInitArgs... )
        {
            for( auto const & [innerMatName, innerDynamicViscosity]: this->expr().exprInnerDynamicViscosities() )
            {
                M_innerExprTensors.emplace( innerMatName, inner_expr_tensors_type{ 
                        mu_expr_tensor_type( std::true_type{}, exprExpanded.exprInnerDynamicViscosity( innerMatName ), ttse, expr.exprInnerDynamicViscosity( innerMatName ), geom, theInitArgs... ),
                        heaviside_expr_tensor_type( std::true_type{}, exprExpanded.exprHeaviside( innerMatName ), ttse, expr.exprHeaviside( innerMatName ), geom, theInitArgs... )
                        } );
            }
        }


        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
        {
            M_muOuterExprTensor.update( geom );
            for( auto & [innerMatName, innerExprTensors]: M_innerExprTensors )
            {
                innerExprTensors.heavisideExprTensor.update( geom );
                innerExprTensors.muExprTensor.update( geom );
                //TODO: optimization->do not compute inner tensor if heaviside is zero
            }

            this->updateImpl();
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            M_muOuterExprTensor.update( std::true_type{}, exprExpanded.exprOuterDynamicViscosity(), ttse, geom, theUpdateArgs... );
            for( auto & [innerMatName, innerExprTensors]: M_innerExprTensors )
            {
                innerExprTensors.heavisideExprTensor.update( std::true_type{}, exprExpanded.exprHeaviside(innerMatName), ttse, geom, theUpdateArgs... );
                innerExprTensors.muExprTensor.update( std::true_type{}, exprExpanded.exprInnerDynamicViscosity(innerMatName), ttse, geom, theUpdateArgs... );
                //TODO: optimization->do not compute inner tensor if heaviside is zero
            }

            this->updateImpl();
        }

        void updateImpl()
        {
            uint16_type nPoints = this->gmc()->nPoints();
            if ( M_localEval.size() != nPoints )
                M_localEval.setConstant( nPoints, super_type::matrix_shape_type::Zero() );

            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                if constexpr ( expr_type::isExprOpID() )
                {
                    value_type muOuter = M_muOuterExprTensor.evalq(0,0,q);
                    M_localEval[q](0,0) = muOuter;
                    for( auto const& [innerMatName, innerExprTensors]: M_innerExprTensors )
                    {
                        value_type H = innerExprTensors.heavisideExprTensor.evalq(0,0,q);
                        value_type muInner = innerExprTensors.muExprTensor.evalq(0,0,q);
                        M_localEval[q](0,0) += ( muInner - muOuter ) * (1.0 - H);
                    }
                }
                else
                {
                    CHECK( false ) << "todo: grad(H) not implemented";
                    //for ( uint16_type c1 = 0; c1 < expr_type::nDim; ++c1 )
                        //M_localEval[q](0,c1) = M_muInnerExprTensor.evalq(0,c1,q) * M_heavisideExprTensor.evalq(0,c1,q);
                }
            }
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                return value_type(0);
            }
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                this->locMatrixShape().setZero();
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "not allow";
                return value_type(0);
            }
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( expr_type::is_applied_as_eval )
                return this->evalq( q );
            else
            {
                CHECK( false) << "not allow";
                return ret_type(this->locMatrixShape().data());
            }
        }

        value_type
            evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
            {
                return M_localEval[q](c1,c2);
            }
        ret_type
            evalq( uint16_type q ) const override
            {
                return ret_type(M_localEval[q].data());
            }


    private :
        this_type const& M_expr;
        mu_expr_tensor_type M_muOuterExprTensor;
        std::map<std::string, inner_expr_tensors_type> M_innerExprTensors;
        array_shape_type M_localEval;
    };


private:
    fluidmec_dynamicviscosity_impl_ptrtype M_fluidMecOuterDynamicViscosityImpl;
    std::map<std::string, fluidmec_dynamicviscosity_impl_ptrtype> M_fluidMecInnerDynamicViscosityImpls;
    std::map<std::string, material_property_scalar_expr_type> M_exprHeavisides;

};


template< typename ExprType>
template <typename TheSymbolExprType>
bool
FluidMecDynamicViscosityBase<ExprType>::hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
{
    auto const& dynamicViscosityLaw = this->expr().dynamicViscosityLaw();
    if ( dynamicViscosityLaw.isNewtonianLaw() )
        return static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosityLaw.isPowerLaw() )
        return static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosityLaw.isCarreauLaw() )
        return static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosityLaw.isCarreauYasudaLaw() )
        return static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosityLaw.isMultifluidLaw() )
        return static_cast<FluidMecDynamicViscosityMultifluid<expr_type> const&>(*this).hasSymbolDependency( symb, se );

    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosityLaw.lawName() <<"\n";
    return false;
}



template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename... TheArgsType>
std::shared_ptr<typename FluidMecDynamicViscosityBase<ExprType>::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
FluidMecDynamicViscosityBase<ExprType>::evaluator( tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, const TheArgsType&... theInitArgs ) const
{
    auto const& dynamicViscosityLaw = this->expr().dynamicViscosityLaw();
    if ( dynamicViscosityLaw.isNewtonianLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosityLaw.isPowerLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosityLaw.isCarreauLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosityLaw.isCarreauYasudaLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosityLaw.isMultifluidLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityMultifluid<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityMultifluid<expr_type> const&>(*this), tensorExprMain, theInitArgs... );

    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosityLaw.lawName() <<"\n";
    return std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >{};
}

template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType,typename... TheArgsType>
std::shared_ptr<typename FluidMecDynamicViscosityBase<ExprType>::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
FluidMecDynamicViscosityBase<ExprType>::evaluator( std::true_type/**/, tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                                                   TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, const TheArgsType&... theInitArgs ) const
{
    using expr_expanded_type = typename TheExprExpandedType::expr_type;

    auto const& dynamicViscosityLaw = this->expr().dynamicViscosityLaw();
    if ( dynamicViscosityLaw.isNewtonianLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{}, static_cast<FluidMecDynamicViscosityNewtonian<expr_expanded_type> const&>(exprExpanded), ttse, static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this), tensorExprMain, theInitArgs... );

    else if ( dynamicViscosityLaw.isPowerLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{}, static_cast<FluidMecDynamicViscosityPowerLaw<expr_expanded_type> const&>(exprExpanded), ttse, static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this), tensorExprMain, theInitArgs... );

    else if ( dynamicViscosityLaw.isCarreauLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{}, static_cast<FluidMecDynamicViscosityCarreau<expr_expanded_type> const&>(exprExpanded), ttse, static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this), tensorExprMain, theInitArgs... );

    else if ( dynamicViscosityLaw.isCarreauYasudaLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{}, static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_expanded_type> const&>(exprExpanded), ttse, static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosityLaw.isMultifluidLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityMultifluid<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{}, static_cast<FluidMecDynamicViscosityMultifluid<expr_expanded_type> const&>(exprExpanded), ttse, static_cast<FluidMecDynamicViscosityMultifluid<expr_type> const&>(*this), tensorExprMain, theInitArgs... );

    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosityLaw.lawName() <<"\n";
    return std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >{};
}


template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
void
FluidMecDynamicViscosityBase<ExprType>::updateEvaluator( std::true_type /**/, std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> > & tensorToUpdate, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                                         Geo_t const& geom, const TheArgsType&... theUpdateArgs )
{
    using expr_expanded_type = typename TheExprExpandedType::expr_type;
    auto const& dynamicViscosityLaw = this->expr().dynamicViscosityLaw();
    if ( dynamicViscosityLaw.isNewtonianLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{}, static_cast<FluidMecDynamicViscosityNewtonian<expr_expanded_type> const&>(exprExpanded), ttse, geom, theUpdateArgs... );

    else if ( dynamicViscosityLaw.isPowerLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{}, static_cast<FluidMecDynamicViscosityPowerLaw<expr_expanded_type> const&>(exprExpanded), ttse, geom, theUpdateArgs... );

    else if ( dynamicViscosityLaw.isCarreauLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{}, static_cast<FluidMecDynamicViscosityCarreau<expr_expanded_type> const&>(exprExpanded), ttse, geom, theUpdateArgs... );

    else if ( dynamicViscosityLaw.isCarreauYasudaLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{}, static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_expanded_type> const&>(exprExpanded), ttse, geom, theUpdateArgs... );

    else if ( dynamicViscosityLaw.isMultifluidLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityMultifluid<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{}, static_cast<FluidMecDynamicViscosityMultifluid<expr_expanded_type> const&>(exprExpanded), ttse, geom, theUpdateArgs... );

    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosityLaw.lawName() <<"\n";
}


template<typename ExprEvaluateFieldOperatorsType, typename FiniteElementVelocityType, typename SymbolsExprType, ExprApplyType ExprApplied, ExprOperatorType ExprOp = ExprOperatorType::ID >
class FluidMecDynamicViscosityImpl// : public Feel::vf::ExprDynamicBase
{
public:

    typedef FluidMecDynamicViscosityImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,SymbolsExprType,ExprApplied,ExprOp> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_velocity|vm::DYNAMIC;

    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorsType;
    using finite_element_velocity_type = FiniteElementVelocityType;
    using dynamic_viscosity_law_type = DynamicViscosityLaw;
    using symbols_expr_type = SymbolsExprType;
    static constexpr ExprApplyType expr_apply_type = ExprApplied;
    static constexpr ExprOperatorType expr_operator_type = ExprOp;

    static constexpr bool is_applied_as_eval = (ExprApplied == ExprApplyType::EVAL);
    static constexpr bool is_applied_as_jacobian = (ExprApplied == ExprApplyType::JACOBIAN);

    static constexpr ExprOperatorType exprOp = ExprOp;

    using value_type = typename expr_evaluate_velocity_opertors_type::value_type;
    static constexpr uint16_type nDim = expr_evaluate_velocity_opertors_type::nDim;

    template <int M,int N>
    using material_property_expr_type = std::decay_t<decltype( expr( ModelExpression{}.template expr<M,N>(),symbols_expr_type{} ) )>;

    using expr_dynamic_viscosity_base_type = FluidMecDynamicViscosityBase<this_type>;

    static const bool is_terminal = true;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = is_applied_as_jacobian && std::is_same_v<Func,FiniteElementVelocityType>;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    static constexpr bool isExprOpID() { return ExprOp == ExprOperatorType::ID; }
    static constexpr bool isExprOpGRAD() { return ExprOp == ExprOperatorType::GRAD; }


    FluidMecDynamicViscosityImpl( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperators,
                                  dynamic_viscosity_law_type const& dynamicViscosityLaw,
                                  MaterialProperties const& matProps,
                                  uint16_type polyOrder,
                                  SymbolsExprType const& se ):
        M_exprEvaluateVelocityOperators( exprEvaluateVelocityOperators ),
        M_dynamicViscosityLaw( dynamicViscosityLaw ),
        M_matProps( matProps ),
        M_polynomialOrder( polyOrder ),
        M_se( se )
    {
        this->initExprDynamicViscosityBase();
    }

    FluidMecDynamicViscosityImpl( FluidMecDynamicViscosityImpl const & op ):
        M_exprEvaluateVelocityOperators( op.M_exprEvaluateVelocityOperators ),
        M_dynamicViscosityLaw( op.M_dynamicViscosityLaw ),
        M_matProps( op.M_matProps ),
        M_polynomialOrder( op.M_polynomialOrder ),
        M_se( op.M_se )
    {
        this->initExprDynamicViscosityBase();
    }

    FluidMecDynamicViscosityImpl( FluidMecDynamicViscosityImpl && op ):
        M_exprEvaluateVelocityOperators( std::move( op.M_exprEvaluateVelocityOperators ) ),
        M_dynamicViscosityLaw( std::move( op.M_dynamicViscosityLaw ) ),
        M_matProps( std::move( op.M_matProps ) ),
        M_polynomialOrder( std::move( op.M_polynomialOrder ) ),
        M_se( std::move( op.M_se ) )
    {
        this->initExprDynamicViscosityBase();
    }

    ~FluidMecDynamicViscosityImpl() {}

    // allow to use the current SymbolsExpr in a concat operation
    // TODO : improve management SymbolsExpr with tensor ctx with concat
    auto const& symbolsExpressionWithoutTensorContext() const
    {
        if constexpr ( is_symbols_expression_v<SymbolsExprType> )
            return M_se;
        else
            return M_se.symbolsExpression();
    }

    symbols_expr_type const& symbolsExpr() const { return M_se; }

    bool dependsOnVelocityField() const { return M_exprDynamicViscosityBase->dependsOnVelocityField(); }

    size_type dynamicContext() const
    {
        return M_exprDynamicViscosityBase->dynamicContext();
    }

    //! polynomial order
    uint16_type polynomialOrder() const
    {
        if ( M_polynomialOrder != invalid_uint16_type_value )
            return M_polynomialOrder;

#if 0 //Thibaut:WHY ?
        uint16_type orderGradVelocity = M_exprEvaluateVelocityOperators->exprGrad().polynomialOrder();
        uint16_type res = 2*(orderGradVelocity+1); // default value for non newtonian
        if ( this->dynamicViscosityLaw().isNewtonianLaw() )
        {
            res = M_exprDynamicViscosityBase->polynomialOrder();
        }
        return res;
#endif
        if( M_exprDynamicViscosityBase )
            return M_exprDynamicViscosityBase->polynomialOrder();
        else
            return M_polynomialOrder;
    }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        // TODO mv this information in SymbolsExpr class
        typedef typename SymbolsExprType::symbols_expr_type::tuple_type symbols_expression_tuple_type;
        static const int nSymbolsExpr = std::decay_t<decltype(hana::size(symbols_expression_tuple_type{}))>::value;
        if constexpr ( nSymbolsExpr > 0 )
        {
            return this->exprDynamicViscosityBase()->hasSymbolDependency( symb, Feel::vf::symbolsExpr( this->symbolsExpressionWithoutTensorContext(), se ) );
        }
        else
            return this->exprDynamicViscosityBase()->hasSymbolDependency( symb, se );
    }

    template <typename OtherSymbolsExprType>
    auto applySymbolsExpr( OtherSymbolsExprType const& se ) const
    {
        auto newse = Feel::vf::symbolsExpr( this->symbolsExpressionWithoutTensorContext(), se );
        using new_se_type = std::decay_t<decltype(newse)>;
        using new_this_type = FluidMecDynamicViscosityImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,new_se_type,ExprApplied,ExprOp>;
        return new_this_type( M_exprEvaluateVelocityOperators,M_dynamicViscosityLaw,M_matProps,M_polynomialOrder,newse );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        CHECK( false ) << "not implemented";
        return *this;
    }


    std::shared_ptr<expr_dynamic_viscosity_base_type> exprDynamicViscosityBase() const { return M_exprDynamicViscosityBase; }
    std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperatorsPtr() const { return M_exprEvaluateVelocityOperators; }

    dynamic_viscosity_law_type const& dynamicViscosityLaw() const { return M_dynamicViscosityLaw; }
    MaterialProperties const& materialProperties() const { return M_matProps; }

    template <int M,int N>
    material_property_expr_type<M,N> materialPropertyExpr( std::string const& prop ) const { return expr( this->materialProperties().property( prop ).template expr<M,N>(), M_se ); }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using tensor_expr_evaluate_velocity_opertors_type = typename expr_evaluate_velocity_opertors_type::template tensor<Geo_t,Basis_i_t,Basis_j_t>;

#if 1
        using shape = typename mpl::if_c< ExprOp == ExprOperatorType::ID,
                                          Shape<tensor_expr_evaluate_velocity_opertors_type::shape_grad_type::nDim, Scalar, false, false>,
                                          Shape<tensor_expr_evaluate_velocity_opertors_type::shape_grad_type::nDim, Vectorial, true, false> >::type;
#else
        using shape = Shape<tensor_expr_evaluate_velocity_opertors_type::shape_grad_type::nDim, Scalar, false, false>;
#endif
        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t,
                           shape,//Shape<tensor_expr_evaluate_velocity_opertors_type::shape_grad_type::nDim, Scalar, false, false>,
                           typename expr_evaluate_velocity_opertors_type::value_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        //typedef typename this_type::value_type value_type;
        using value_type = typename tensorbase_type::value_type;
        //using key_type = typename tensorbase_type::key_type;
        using gmc_type = typename tensorbase_type::gmc_type;
        using gm_type = typename tensorbase_type::gm_type;
        //using shape = typename tensorbase_type::shape_type;
        using matrix_shape_type = typename tensorbase_type::matrix_shape_type;
        using ret_type = typename tensorbase_type::ret_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            M_expr( expr )
        {
            this->initTensor( expr, true, geom, fev, feu );
        }

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev ):
            M_expr( expr )
        {
            this->initTensor( expr, true, geom, fev );
        }

        tensor( this_type const& expr, Geo_t const& geom ):
            M_expr( expr )
        {
            this->initTensor( expr, true, geom );
        }

        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators,  Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ):
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom, fev, feu );
        }

        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators, Geo_t const& geom, Basis_i_t const& fev ):
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom, fev );
        }

        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators, Geo_t const& geom ):
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs ):
            M_expr( expr )
        {
            this->initTensor( std::true_type{}, true, exprExpanded, ttse, expr, geom, theInitArgs... );
        }


        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperatorsPtr() const { return M_tensorExprEvaluateVelocityOperators; }
        tensor_expr_evaluate_velocity_opertors_type const& tensorExprEvaluateVelocityOperators() const { return *M_tensorExprEvaluateVelocityOperators; }

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom, bool upEvaluateVelocityOperators = true )
        {
            if ( upEvaluateVelocityOperators )
                M_tensorExprEvaluateVelocityOperators->update( geom );
            M_tensorbase->update( geom );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            this->update( std::true_type{}, true, exprExpanded, ttse, geom, theUpdateArgs... );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, bool upEvaluateVelocityOperators, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
        {
            if ( upEvaluateVelocityOperators )
                M_tensorExprEvaluateVelocityOperators->update( geom, theUpdateArgs... );
            M_expr.exprDynamicViscosityBase()->template updateEvaluator<Geo_t, Basis_i_t, Basis_j_t>( std::true_type{}, M_tensorbase, *(exprExpanded.exprDynamicViscosityBase()), ttse,  geom, theUpdateArgs...);
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i,c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return M_tensorbase->evalq( q );
        }

        template<typename... TheArgsType>
        void initTensor( this_type const& expr, bool initEvaluateVelocityOperators, const TheArgsType&... theInitArgs )
        {
            if ( initEvaluateVelocityOperators )
                M_tensorExprEvaluateVelocityOperators = std::make_shared<tensor_expr_evaluate_velocity_opertors_type>( *(expr.exprEvaluateVelocityOperatorsPtr()), theInitArgs... );
            M_tensorbase = expr.exprDynamicViscosityBase()->template evaluator<Geo_t, Basis_i_t, Basis_j_t>( *this, theInitArgs... );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void initTensor( std::true_type /**/, bool initEvaluateVelocityOperators, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                         this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
        {
            if ( initEvaluateVelocityOperators )
                M_tensorExprEvaluateVelocityOperators = std::make_shared<tensor_expr_evaluate_velocity_opertors_type>( *(expr.exprEvaluateVelocityOperatorsPtr()), geom, theInitArgs... );
            M_tensorbase = expr.exprDynamicViscosityBase()->template evaluator<Geo_t, Basis_i_t, Basis_j_t>( std::true_type{}, *this, *(exprExpanded.exprDynamicViscosityBase()), ttse, geom, theInitArgs... );
        }

    private:
        this_type const& M_expr;
        tensorbase_ptrtype M_tensorbase;
        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
    };
private :
    void initExprDynamicViscosityBase()
    {
        auto const& dynamicViscosityLaw = this->dynamicViscosityLaw();
        if ( dynamicViscosityLaw.isNewtonianLaw() )
            M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityNewtonian<this_type>( *this ) );
        else if ( dynamicViscosityLaw.isPowerLaw() )
            M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityPowerLaw<this_type>( *this ) );
        else if ( dynamicViscosityLaw.isCarreauLaw() )
            M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityCarreau<this_type>( *this ) );
        else if ( dynamicViscosityLaw.isCarreauYasudaLaw() )
            M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityCarreauYasuda<this_type>( *this ) );
        else if( dynamicViscosityLaw.isMultifluidLaw() )
            M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityMultifluid<this_type>( *this ) );
        else
            CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosityLaw.lawName() <<"\n";
    }

private:
    std::shared_ptr<expr_dynamic_viscosity_base_type> M_exprDynamicViscosityBase;
    std::shared_ptr<expr_evaluate_velocity_opertors_type> M_exprEvaluateVelocityOperators;
    dynamic_viscosity_law_type const& M_dynamicViscosityLaw;
    MaterialProperties const& M_matProps;
    uint16_type M_polynomialOrder;
    SymbolsExprType M_se;
};

template<class ExprGradVelocityType, typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecViscosity( Expr<ExprGradVelocityType> const& grad_u,
                   DynamicViscosityLaw const& dynamicViscosityLaw,
                   MaterialProperties const& matProps,
                   SymbolsExprType const& se = symbols_expression_empty_t{},
                   uint16_type polyOrder = invalid_uint16_type_value )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorGradFromExpr<Expr<ExprGradVelocityType>>;
    typedef FluidMecDynamicViscosityImpl<expr_evaluate_velocity_opertors_type, std::nullptr_t, SymbolsExprType, ExprApplyType::EVAL> fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( grad_u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators, dynamicViscosityLaw, matProps, polyOrder, se ) );
}

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H
