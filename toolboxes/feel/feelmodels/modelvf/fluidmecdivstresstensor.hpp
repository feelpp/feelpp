#ifndef FEELPP_MODELS_VF_FLUIDMEC_DIVSTRESSTENSOR_H
#define FEELPP_MODELS_VF_FLUIDMEC_DIVSTRESSTENSOR_H 1

#include <feel/feelmodels/modelvf/fluidmecdynamicviscosity.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename ExprEvaluateFieldOperatorsType, typename FiniteElementVelocityType, typename ModelPhysicFluidType, typename SymbolsExprType, ExprApplyType ExprApplied>
class FluidMecDivStressTensorImpl // : public Feel::vf::ExprDynamicBase
{
public:

    typedef FluidMecDivStressTensorImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ModelPhysicFluidType,SymbolsExprType,ExprApplied> this_type;


    using model_physic_fluid_type = ModelPhysicFluidType;
    using symbols_expr_type = SymbolsExprType;
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorsType;
    using value_type = typename expr_evaluate_velocity_opertors_type::value_type;

    static const size_type context_velocity = expr_evaluate_velocity_opertors_type::laplacian_is_zero? vm::JACOBIAN|vm::KB|vm::GRAD : vm::JACOBIAN|vm::KB|vm::GRAD|vm::LAPLACIAN|vm::SECOND_DERIVATIVE;
    static const size_type context = context_velocity|vm::DYNAMIC;
    static const bool is_terminal = true;

    static constexpr bool is_applied_as_eval = ExprApplied == ExprApplyType::EVAL;
    static constexpr bool is_applied_as_jacobian = ExprApplied == ExprApplyType::JACOBIAN;
    static constexpr bool is_applied_as_linear_trial = ExprApplied == ExprApplyType::LINEAR_TRIAL;
    static constexpr bool is_applied_as_linear_test = ExprApplied == ExprApplyType::LINEAR_TEST;
    static constexpr bool is_applied_as_jacobian_or_linear_trial = is_applied_as_jacobian || is_applied_as_linear_trial;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = is_applied_as_linear_test && std::is_same_v<Func,FiniteElementVelocityType>;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = is_applied_as_jacobian_or_linear_trial && std::is_same_v<Func,FiniteElementVelocityType>;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    template <int M,int N>
    using material_property_expr_type = std::decay_t<decltype( expr( ModelExpression{}.template expr<M,N>(),symbols_expr_type{} ) )>;

    using material_property_scalar_expr_type = material_property_expr_type<1,1>;

    using expr_dynamic_viscosity_impl_type = FluidMecDynamicViscosityImpl<expr_evaluate_velocity_opertors_type,FiniteElementVelocityType,ModelPhysicFluidType,SymbolsExprType,
                                                                          is_applied_as_jacobian? ExprApplyType::JACOBIAN : ExprApplyType::EVAL,
                                                                          //mpl::int_< this_type::specific_expr_type::value == ExprApplyType::EVAL? ExprApplyType::EVAL : ExprApplyType::JACOBIAN>
                                                                          ExprOperatorType::ID>;
    using expr_dynamic_viscosity_type = Expr<expr_dynamic_viscosity_impl_type>;

    using expr_grad_dynamic_viscosity_impl_type = FluidMecDynamicViscosityImpl<expr_evaluate_velocity_opertors_type,FiniteElementVelocityType,ModelPhysicFluidType,SymbolsExprType,
                                                                               is_applied_as_jacobian? ExprApplyType::JACOBIAN : ExprApplyType::EVAL,
                                                                               //mpl::int_< this_type::specific_expr_type::value == ExprApplyType::EVAL? ExprApplyType::EVAL : ExprApplyType::JACOBIAN>,
                                                                               ExprOperatorType::GRAD>;
    using expr_grad_dynamic_viscosity_type = Expr<expr_grad_dynamic_viscosity_impl_type>;

    using expr_grad_material_property_scalar_type = std::decay_t<decltype( grad<expr_evaluate_velocity_opertors_type::field_clean_type::nDim>( material_property_scalar_expr_type{} ) )>;
    // using expr_grad_turbulent_dynamic_viscosity_type = expr_grad_material_property_scalar_type;


    FluidMecDivStressTensorImpl( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvalVelocityOperators,
                                 model_physic_fluid_type const& physicFluid,
                                 MaterialProperties const& matProps,
                                 WorldComm const& world, std::string const& dirLibExpr,
                                 uint16_type polyOrder,
                                 SymbolsExprType const& se,
                                 bool checkViscosityDependencyOnCoordinates )
        :
        M_exprEvaluateVelocityOperators( exprEvalVelocityOperators ),
        M_physicFluid( physicFluid ),
        M_matProps( matProps ),
        M_polynomialOrder( polyOrder ),
        //M_withPressureTerm( withPressure ),
        M_se( se )
    {
        M_exprEvaluateVelocityOperators->setEnableLaplacian( true );
        M_exprDynamicViscosity.emplace( expr_dynamic_viscosity_impl_type( M_exprEvaluateVelocityOperators,physicFluid,matProps,invalid_uint16_type_value,se ) );

        bool dynamicViscosityDependsOnCoordinatesInSpace = M_exprDynamicViscosity->template hasSymbolDependencyOnCoordinatesInSpace<expr_evaluate_velocity_opertors_type::field_clean_type::nRealDim>() && checkViscosityDependencyOnCoordinates;
        if ( dynamicViscosityDependsOnCoordinatesInSpace )
        {
            M_exprEvaluateVelocityOperators->setEnableGrad( true );
            M_exprGradDynamicViscosity.emplace( expr_grad_dynamic_viscosity_impl_type( M_exprEvaluateVelocityOperators,physicFluid,matProps,invalid_uint16_type_value,se ) );
        }

        if ( this->turbulence().isEnabled() )
        {
            M_exprEvaluateVelocityOperators->setEnableGrad( true );
            M_exprTurbulentDynamicVisocity.emplace( this->materialPropertyExpr<1,1>("turbulent-dynamic-viscosity") );
            M_exprGradTurbulentDynamicVisocity.emplace( grad<expr_evaluate_velocity_opertors_type::field_clean_type::nDim>( *M_exprTurbulentDynamicVisocity, "", world, dirLibExpr ) );
            if ( this->turbulence().hasTurbulentKineticEnergy() )
            {
                M_exprTurbulentKineticEnergy.emplace( this->materialPropertyExpr<1,1>("turbulent-kinetic-energy") );
                M_exprDensity.emplace( this->materialPropertyExpr<1,1>("density") );
                M_exprGradTurbulentKineticEnergy.emplace( grad<expr_evaluate_velocity_opertors_type::field_clean_type::nDim>( *M_exprTurbulentKineticEnergy, "", world, dirLibExpr ) );
            }
        }
    }

    FluidMecDivStressTensorImpl( FluidMecDivStressTensorImpl const & op ) = default;
    FluidMecDivStressTensorImpl( FluidMecDivStressTensorImpl && op ) = default;

    size_type dynamicContext() const
        {
            size_type res = M_exprDynamicViscosity->dynamicContext();
            if ( M_exprGradDynamicViscosity )
                res = res | Feel::vf::dynamicContext( *M_exprGradDynamicViscosity );
            if ( this->turbulence().isEnabled() )
            {
                res= res | Feel::vf::dynamicContext( *M_exprTurbulentDynamicVisocity ) | Feel::vf::dynamicContext( *M_exprGradTurbulentDynamicVisocity );
                if (  M_exprGradTurbulentKineticEnergy )
                    res= res | Feel::vf::dynamicContext( *M_exprGradTurbulentKineticEnergy ) |  Feel::vf::dynamicContext( *M_exprDensity );
            }
            return res;
        }

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            if ( M_polynomialOrder != invalid_uint16_type_value )
                return M_polynomialOrder;

            uint16_type orderGradVelocity =  M_exprEvaluateVelocityOperators->isEnabledGrad()?  M_exprEvaluateVelocityOperators->exprGrad().polynomialOrder() : 0;
            uint16_type orderLaplacianVelocity = expr_evaluate_velocity_opertors_type::laplacian_is_zero? 0 : M_exprEvaluateVelocityOperators->exprLaplacian().polynomialOrder();

            uint16_type res = 0;
            if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                 res = M_exprDynamicViscosity->polynomialOrder() + orderLaplacianVelocity;

            if ( M_exprGradDynamicViscosity )
                res = std::max( res, (uint16_type) (M_exprGradDynamicViscosity->polynomialOrder() + orderGradVelocity ) );

            if ( this->turbulence().isEnabled() )
            {
                res = std::max( res, (uint16_type) (M_exprTurbulentDynamicVisocity->polynomialOrder() + orderLaplacianVelocity) );
                res = std::max( res, (uint16_type) (M_exprGradTurbulentDynamicVisocity->polynomialOrder() + orderGradVelocity) );
                // TODO turbument kinematic energy
            }

            return res;
        }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    template <typename OtherSymbolsExprType>
    auto applySymbolsExpr( OtherSymbolsExprType const& se ) const
        {
            auto newse = Feel::vf::symbolsExpr( M_se, se );
            using new_se_type = std::decay_t<decltype(newse)>;
            using new_this_type = FluidMecDivStressTensorImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ModelPhysicFluidType,new_se_type,ExprApplied>;
            return new_this_type( M_exprEvaluateVelocityOperators,M_physicFluid,M_matProps,M_polynomialOrder,newse );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            CHECK( false ) << "not implemented";
            return *this;
        }

    bool isZero() const
        {
            if ( expr_evaluate_velocity_opertors_type::laplacian_is_zero /*&& !M_withPressureTerm*/ && !turbulence().isEnabled() && !M_exprGradDynamicViscosity )
                return true;
            return false;
        }

    MaterialProperties const& materialProperties() const { return M_matProps; }

    template <int M,int N>
    material_property_expr_type<M,N> materialPropertyExpr( std::string const& prop ) const { return expr( this->materialProperties().property( prop ).template expr<M,N>(), M_se ); }

    std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperatorsPtr() const { return M_exprEvaluateVelocityOperators; }
    expr_dynamic_viscosity_type const& exprDynamicVisocity() const { return *M_exprDynamicViscosity; }
    expr_grad_dynamic_viscosity_type const& exprGradDynamicViscosity() const { return *M_exprGradDynamicViscosity; }
    material_property_scalar_expr_type const& exprTurbulentDynamicVisocity() const { return *M_exprTurbulentDynamicVisocity; }
    expr_grad_material_property_scalar_type const& exprGradTurbulentDynamicVisocity() const { return *M_exprGradTurbulentDynamicVisocity; }
    expr_grad_material_property_scalar_type const& exprGradTurbulentKineticEnergy() const { return *M_exprGradTurbulentKineticEnergy; }
    material_property_scalar_expr_type const& exprDensity() const { return *M_exprDensity; }
    //bool withPressureTerm() const { return M_withPressureTerm; }

    bool hasExprGradDynamicViscosity() const { return M_exprGradDynamicViscosity? true : false; }

    auto const& turbulence() const { return M_physicFluid.turbulence(); }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                      Shape<expr_evaluate_velocity_opertors_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape_grad_type::nDim,Vectorial, false, false>,
                                      typename this_type::value_type>
    {
        using super_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                      Shape<expr_evaluate_velocity_opertors_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape_grad_type::nDim,Vectorial, false, false>,
                                      typename this_type::value_type>;
        typedef typename this_type::value_type value_type;

        using key_type = typename super_type::key_type;
        using gmc_type = typename super_type::gmc_type;
        using gm_type = typename super_type::gm_type;

        using shape = typename super_type::shape_type;
        using matrix_shape_type = typename super_type::matrix_shape_type;
        using array_shape_type = typename super_type::new_array_shape_type;
        using ret_type = typename super_type::ret_type;

        using tensor_expr_material_property_scalar_type = typename this_type::material_property_scalar_expr_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_evaluate_velocity_opertors_type = typename expr_evaluate_velocity_opertors_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_dynamic_viscosity_type = typename expr_dynamic_viscosity_impl_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_grad_dynamic_viscosity_type = typename expr_grad_dynamic_viscosity_impl_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_grad_material_property_scalar_type = typename expr_grad_material_property_scalar_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_grad_turbulent_dynamic_viscosity_type = tensor_expr_grad_material_property_scalar_type;

        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_hasExprGradDynamicViscosity( expr.hasExprGradDynamicViscosity() ),
            M_useTurbulenceModel( expr.turbulence().isEnabled() ),
            M_dynamicViscosityDependsOnVelocityField( false )
            {
                this->initTensor( expr, geom, fev, feu );
            }

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_hasExprGradDynamicViscosity( expr.hasExprGradDynamicViscosity() ),
            M_useTurbulenceModel( expr.turbulence().isEnabled() ),
            M_dynamicViscosityDependsOnVelocityField( false )
            {
                this->initTensor( expr, geom, fev );
            }
        tensor( this_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_hasExprGradDynamicViscosity( expr.hasExprGradDynamicViscosity() ),
            M_useTurbulenceModel( expr.turbulence().isEnabled() ),
            M_dynamicViscosityDependsOnVelocityField( false )
            {
                this->initTensor( expr, geom );
            }
        template<typename IM>
        void init( IM const& im ) {}

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) override
        {
            this->update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ ) override
        {
            this->update(geom);
        }
        void update( Geo_t const& geom ) override
        {
            this->setGmc( geom );
            M_tensorExprEvaluateVelocityOperators->update( geom );
            if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                M_tensorDynamicViscosity->update( geom, false );
            if ( M_tensorGradDynamicViscosity )
                M_tensorGradDynamicViscosity->update( geom, false );
            if ( M_useTurbulenceModel )
            {
                if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                    M_tensorExprTurbulentDynamicVisocity->update( geom );
                M_tensorExprGradTurbulentDynamicVisocity->update( geom );
                if ( M_tensorGradTurbulentKineticEnergy )
                {
                    M_tensorGradTurbulentKineticEnergy->update( geom );
                    M_tensorExprDensity->update( geom );
                }
            }
            this->updateImpl();
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            this->setGmc( geom );
            M_tensorExprEvaluateVelocityOperators->update( geom, face );
            if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                  M_tensorDynamicViscosity->update( geom, face, false );
            if ( M_tensorGradDynamicViscosity )
                M_tensorGradDynamicViscosity->update( geom, face, false );
            if ( M_useTurbulenceModel )
            {
                if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                    M_tensorExprTurbulentDynamicVisocity->update( geom, face );
                M_tensorExprGradTurbulentDynamicVisocity->update( geom, face );
                if ( M_tensorGradTurbulentKineticEnergy )
                {
                    M_tensorGradTurbulentKineticEnergy->update( geom, face );
                    M_tensorExprDensity->update( geom, face );
                }
            }
            this->updateImpl();
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( this_type::is_applied_as_eval )
                return this->evalq( q );
            else if constexpr ( this_type::is_applied_as_jacobian_or_linear_trial )
                return this->evalImpl_iq_or_jq( j, q );
            else
                return this->evalImpl_iq_or_jq( i, q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( this_type::is_applied_as_eval )
                return this->evalq( c1,c2,q );
            else
            {
                return this->evalijq(i,j,q)(c1,c2);
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( this_type::is_applied_as_eval )
                  return this->evalq( c1,c2,q );
            else if constexpr ( this_type::is_applied_as_linear_test )
                  return this->evaliq( i,q )(c1,c2);
            CHECK( false ) << "not allow";
            return value_type(0);// this->evalq( c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( this_type::is_applied_as_eval )
                  return this->evalq( q );
            else if constexpr ( this_type::is_applied_as_linear_test )
                  return this->evalImpl_iq_or_jq( i, q );
            CHECK( false ) << "not allow";
            return this->evalq( q );
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
        template<typename... TheArgsType>
        void initTensor( this_type const& expr, const TheArgsType&... theInitArgs )
            {
                M_tensorExprEvaluateVelocityOperators = std::make_shared<tensor_expr_evaluate_velocity_opertors_type>( *(expr.exprEvaluateVelocityOperatorsPtr()), theInitArgs... );
                if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                    M_tensorDynamicViscosity.emplace( expr.exprDynamicVisocity().expression(), M_tensorExprEvaluateVelocityOperators, theInitArgs... );

                if ( M_hasExprGradDynamicViscosity )
                    M_tensorGradDynamicViscosity.emplace( expr.exprGradDynamicViscosity().expression(), M_tensorExprEvaluateVelocityOperators, theInitArgs... );

                M_dynamicViscosityDependsOnVelocityField = expr.exprDynamicVisocity().expression().dependsOnVelocityField();
                CHECK( !M_dynamicViscosityDependsOnVelocityField ) << "TODO : non newtonian case";

                if ( M_useTurbulenceModel )
                {
                    if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                        M_tensorExprTurbulentDynamicVisocity.emplace( expr.exprTurbulentDynamicVisocity(), theInitArgs... );
                    M_tensorExprGradTurbulentDynamicVisocity.emplace( expr.exprGradTurbulentDynamicVisocity(), theInitArgs... );
                    if ( expr.turbulence().hasTurbulentKineticEnergy() )
                    {
                        M_tensorGradTurbulentKineticEnergy.emplace( expr.exprGradTurbulentKineticEnergy(), theInitArgs... );
                        M_tensorExprDensity.emplace( expr.exprDensity(), theInitArgs... );
                    }
                }
            }
        void updateImpl()
            {
                if constexpr ( !this_type::is_applied_as_eval )
                      return;

                uint16_type nPoints = this->gmc()->nPoints();
                if ( M_localEval.size() != nPoints )
                    M_localEval.resize( nPoints );

                bool hasTurbulentKineticEnergy = ( M_useTurbulenceModel && M_tensorGradTurbulentKineticEnergy )? true : false;

                for ( uint16_type q=0;q< nPoints;++q )
                {
                    matrix_shape_type & locMat = M_localEval[q];

                    if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                    {
                        value_type muEval = M_tensorDynamicViscosity->evalq(0,0,q);
                        if ( M_useTurbulenceModel )
                            muEval += M_tensorExprTurbulentDynamicVisocity->evalq(0,0,q);
                        auto const& laplacianVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalLaplacian( q );
                        locMat = muEval*laplacianVelocityEval;
                    }
                    else
                        locMat.setZero();

                    if ( M_hasExprGradDynamicViscosity || M_useTurbulenceModel )
                    {
                        auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        if constexpr ( gmc_type::nDim == 2 )
                        {
                            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                            value_type gradMuTurb_0 = 0, gradMuTurb_1 = 0;
                            if ( M_hasExprGradDynamicViscosity )
                            {
                                gradMuTurb_0 += M_tensorGradDynamicViscosity->evalq(0,0,q);
                                gradMuTurb_1 += M_tensorGradDynamicViscosity->evalq(0,1,q);
                            }
                            if ( M_useTurbulenceModel )
                            {
                                gradMuTurb_0 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,0,q);
                                gradMuTurb_1 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                            }
                            locMat(0,0) += 2*du1vdx*gradMuTurb_0 + ( du2vdx + du1vdy )*gradMuTurb_1;
                            locMat(1,0) += ( du2vdx + du1vdy )*gradMuTurb_0 + 2*du2vdy*gradMuTurb_1;
                        }
                        else
                        {
                            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                            value_type gradMuTurb_0 = 0, gradMuTurb_1 = 0, gradMuTurb_2 = 0;
                            if ( M_hasExprGradDynamicViscosity )
                            {
                                gradMuTurb_0 += M_tensorGradDynamicViscosity->evalq(0,0,q);
                                gradMuTurb_1 += M_tensorGradDynamicViscosity->evalq(0,1,q);
                                gradMuTurb_2 += M_tensorGradDynamicViscosity->evalq(0,2,q);
                            }
                            if ( M_useTurbulenceModel )
                            {
                                gradMuTurb_0 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,0,q);
                                gradMuTurb_1 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                                gradMuTurb_2 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,2,q);
                            }
                            locMat(0,0) += 2*du1vdx*gradMuTurb_0 + ( du2vdx + du1vdy )*gradMuTurb_1 + ( du3vdx + du1vdz )*gradMuTurb_2;
                            locMat(1,0) += ( du2vdx + du1vdy )*gradMuTurb_0 + 2*du2vdy*gradMuTurb_1 + ( du2vdz + du3vdy )*gradMuTurb_2;
                            locMat(2,0) += ( du3vdx + du1vdz )*gradMuTurb_0 + ( du2vdz + du3vdy )*gradMuTurb_1 + 2*du3vdz*gradMuTurb_2;
                        }

                        if ( hasTurbulentKineticEnergy )
                        {
                            const value_type rho = M_tensorExprDensity->evalq(0,0,q);
                            const value_type gradTKE_0 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,0,q);
                            const value_type gradTKE_1 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                            const value_type rhoWithFactor = (2./3.)*rho;
                            locMat(0,0) -= rhoWithFactor*gradTKE_0;
                            locMat(1,0) -= rhoWithFactor*gradTKE_1;
                            if constexpr ( gmc_type::nDim == 3 )
                            {
                                const value_type gradTKE_2 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,2,q);
                                locMat(2,0) -= rhoWithFactor*gradTKE_2;
                            }
                        }
                    }
                }

            }

#if 0
        template <bool IsTrial = this_type::is_applied_as_jacobian_or_linear_trial, std::enable_if_t< IsTrial, bool> = true>
        auto const& fecTrialOrTest() const { return this->fecTrial(); }
        template <bool IsTest = this_type::is_applied_as_linear_test, std::enable_if_t< IsTest, bool> = true>
        auto const& fecTrialOrTest() const { return this->fecTest(); }
#else
        template <class TT = this_type, std::enable_if_t< TT::is_applied_as_jacobian_or_linear_trial, bool> = true>
        auto const& fecTrialOrTest() const { return this->fecTrial(); }
        template <class TT = this_type, std::enable_if_t< TT::is_applied_as_linear_test, bool> = true>
        auto const& fecTrialOrTest() const { return this->fecTest(); }
#endif

        ret_type
        evalImpl_iq_or_jq( uint16_type i_or_j, uint16_type q ) const
        {
            if constexpr ( ( this_type::is_applied_as_jacobian_or_linear_trial && std::is_same_v<typename super_type::basis_fec_trial_type::reference_element_type,FiniteElementVelocityType> ) ||
                                ( this_type::is_applied_as_linear_test && std::is_same_v<typename super_type::basis_fec_test_type::reference_element_type,FiniteElementVelocityType> ) ) // NEED TO UNDERSTAND
            {
                matrix_shape_type & locMat = this->locMatrixShape();

                if constexpr ( !expr_evaluate_velocity_opertors_type::laplacian_is_zero )
                {
                    value_type muEval = M_tensorDynamicViscosity->evalq(0,0,q);
                    if ( M_useTurbulenceModel )
                        muEval += M_tensorExprTurbulentDynamicVisocity->evalq(0,0,q);
                    auto const& laplacianTrial = this->fecTrialOrTest()->laplacian( i_or_j, q );
                    for ( uint16_type d=0;d<gmc_type::nDim;++d )
                        locMat(d,0) = muEval*laplacianTrial(d,0);
                }
                else
                    locMat.setZero();

                if ( M_hasExprGradDynamicViscosity || M_useTurbulenceModel )
                {
                    auto const& gradTrial = this->fecTrialOrTest()->grad( i_or_j, q );
                    if constexpr ( gmc_type::nDim == 2 )
                    {
                        const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
                        const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
                        value_type gradMuTurb_0 = 0, gradMuTurb_1 = 0;
                        if ( M_hasExprGradDynamicViscosity )
                        {
                            gradMuTurb_0 += M_tensorGradDynamicViscosity->evalq(0,0,q);
                            gradMuTurb_1 += M_tensorGradDynamicViscosity->evalq(0,1,q);
                        }
                        if ( M_useTurbulenceModel )
                        {
                            gradMuTurb_0 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,0,q);
                            gradMuTurb_1 += M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                        }
                        locMat(0,0) += 2*du1tdx*gradMuTurb_0 + ( du2tdx + du1tdy )*gradMuTurb_1;
                        locMat(1,0) += ( du2tdx + du1tdy )*gradMuTurb_0 + 2*du2tdy*gradMuTurb_1;
                    }
                    else
                    {
                        const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
                        const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
                        const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
                        value_type gradMuTurb_0 = 0, gradMuTurb_1 = 0, gradMuTurb_2 = 0;
                        if ( M_hasExprGradDynamicViscosity )
                        {
                            gradMuTurb_0 += M_tensorGradDynamicViscosity->evalq(0,0,q);
                            gradMuTurb_1 += M_tensorGradDynamicViscosity->evalq(0,1,q);
                            gradMuTurb_2 += M_tensorGradDynamicViscosity->evalq(0,1,q);
                        }
                        if ( M_useTurbulenceModel )
                        {
                            gradMuTurb_0 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,0,q);
                            gradMuTurb_1 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                            gradMuTurb_2 = M_tensorExprGradTurbulentDynamicVisocity->evalq(0,1,q);
                        }
                        locMat(0,0) += 2*du1tdx*gradMuTurb_0 + ( du2tdx + du1tdy )*gradMuTurb_1 + ( du3tdx + du1tdz )*gradMuTurb_2;
                        locMat(1,0) += ( du2tdx + du1tdy )*gradMuTurb_0 + 2*du2tdy*gradMuTurb_1 + ( du2tdz + du3tdy )*gradMuTurb_2;
                        locMat(2,0) += ( du3tdx + du1tdz )*gradMuTurb_0 + ( du2tdz + du3tdy )*gradMuTurb_1 + 2*du3tdz*gradMuTurb_2;
                    }
                }
                return ret_type(this->locMatrixShape().data());
            }
            CHECK( false ) << "not allow";
            return ret_type(this->locMatrixShape().data());
        }

    private:
        this_type const& M_expr;

        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        std::optional<tensor_expr_dynamic_viscosity_type> M_tensorDynamicViscosity;
        std::optional<tensor_expr_grad_dynamic_viscosity_type> M_tensorGradDynamicViscosity;
        std::optional<tensor_expr_material_property_scalar_type> M_tensorExprTurbulentDynamicVisocity, M_tensorExprDensity;
        std::optional<tensor_expr_grad_turbulent_dynamic_viscosity_type> M_tensorExprGradTurbulentDynamicVisocity, M_tensorGradTurbulentKineticEnergy;

        array_shape_type M_localEval;

        bool M_hasExprGradDynamicViscosity;
        bool M_useTurbulenceModel;
        bool M_dynamicViscosityDependsOnVelocityField;

    };
private :
    std::shared_ptr<expr_evaluate_velocity_opertors_type> M_exprEvaluateVelocityOperators;

    model_physic_fluid_type const& M_physicFluid;
    MaterialProperties const& M_matProps;
    uint16_type M_polynomialOrder;
    //bool M_withPressureTerm;
    SymbolsExprType M_se;

    std::optional<expr_dynamic_viscosity_type> M_exprDynamicViscosity;
    std::optional<expr_grad_dynamic_viscosity_type> M_exprGradDynamicViscosity;
    std::optional<material_property_scalar_expr_type> M_exprTurbulentDynamicVisocity, M_exprTurbulentKineticEnergy, M_exprDensity;
    std::optional<expr_grad_material_property_scalar_type> M_exprGradTurbulentDynamicVisocity, M_exprGradTurbulentKineticEnergy;
};

template<class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecDivViscousStressTensor( VelocityFieldType const& u,
                                ModelPhysicFluidType const& physicFluid,
                                MaterialProperties const& matProps,
                                WorldComm const& world, std::string const& dirLibExpr,
                                SymbolsExprType const& se = symbols_expression_empty_t{},
                                bool checkViscosityDependencyOnCoordinates = true,
                                uint16_type polyOrder = invalid_uint16_type_value
                                )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperators<VelocityFieldType>;
    typedef FluidMecDivStressTensorImpl<expr_evaluate_velocity_opertors_type,std::nullptr_t,ModelPhysicFluidType,SymbolsExprType,ExprApplyType::EVAL > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,physicFluid,matProps,world,dirLibExpr,polyOrder,se,checkViscosityDependencyOnCoordinates ) );
}

template<ExprApplyType ExprApplied,class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t>
inline
auto
fluidMecDivViscousStressTensorTestTrialApplied( VelocityFieldType const& u,
                                                ModelPhysicFluidType const& physicFluid,
                                                MaterialProperties const& matProps,
                                                WorldComm const& world, std::string const& dirLibExpr,
                                                SymbolsExprType const& se = symbols_expression_empty_t{},
                                                bool checkViscosityDependencyOnCoordinates = true,
                                                uint16_type polyOrder = invalid_uint16_type_value
                                                )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperators<VelocityFieldType>;
    typedef FluidMecDivStressTensorImpl<expr_evaluate_velocity_opertors_type,typename unwrap_ptr_t<VelocityFieldType>::functionspace_type::reference_element_type,ModelPhysicFluidType,SymbolsExprType,ExprApplied > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,physicFluid,matProps,world,dirLibExpr,polyOrder,se,checkViscosityDependencyOnCoordinates ) );
}

template<class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecDivViscousStressTensorJacobian( VelocityFieldType const& u,
                                        ModelPhysicFluidType const& physicFluid,
                                        MaterialProperties const& matProps,
                                        WorldComm const& world, std::string const& dirLibExpr,
                                        SymbolsExprType const& se = symbols_expression_empty_t{},
                                        bool checkViscosityDependencyOnCoordinates = true,
                                        uint16_type polyOrder = invalid_uint16_type_value
                                        )
{
    return fluidMecDivViscousStressTensorTestTrialApplied<ExprApplyType::JACOBIAN>(u,physicFluid,matProps,world,dirLibExpr,se,checkViscosityDependencyOnCoordinates,polyOrder);
}

template<class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecDivViscousStressTensorLinearTrial( VelocityFieldType const& u,
                                           ModelPhysicFluidType const& physicFluid,
                                           MaterialProperties const& matProps,
                                           WorldComm const& world, std::string const& dirLibExpr,
                                           SymbolsExprType const& se = symbols_expression_empty_t{},
                                           bool checkViscosityDependencyOnCoordinates = true,
                                           uint16_type polyOrder = invalid_uint16_type_value
                                           )
{
    return fluidMecDivViscousStressTensorTestTrialApplied<ExprApplyType::LINEAR_TRIAL>(u,physicFluid,matProps,world,dirLibExpr,se,checkViscosityDependencyOnCoordinates,polyOrder);
}

template<class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecDivViscousStressTensorLinearTest( VelocityFieldType const& u,
                                          ModelPhysicFluidType const& physicFluid,
                                          MaterialProperties const& matProps,
                                          WorldComm const& world, std::string const& dirLibExpr,
                                          SymbolsExprType const& se = symbols_expression_empty_t{},
                                          bool checkViscosityDependencyOnCoordinates = true,
                                          uint16_type polyOrder = invalid_uint16_type_value
                                          )
{
    return fluidMecDivViscousStressTensorTestTrialApplied<ExprApplyType::LINEAR_TEST>(u,physicFluid,matProps,world,dirLibExpr,se,checkViscosityDependencyOnCoordinates,polyOrder);
}

} // namespace FeelModels
} // namespace Feel
#endif /* FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H */
