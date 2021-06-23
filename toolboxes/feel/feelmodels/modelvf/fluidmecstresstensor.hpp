/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 2013-11-19

  Copyright (C) 2012 Universite Joseph Fourier (Grenoble I)

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
   \file fluidmecstresstensor.hpp
   \author Vincent Chabannes
   \date 2013-11-19
 */
#ifndef FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H
#define FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H 1

#include <feel/feelmodels/modelvf/fluidmecdynamicviscosity.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename ExprEvaluateFieldOperatorsType, typename FiniteElementVelocityType, typename ExprIdPressureType, typename ModelPhysicFluidType, typename SymbolsExprType, typename SpecificExprType>
class FluidMecStressTensorImpl : public Feel::vf::ExprDynamicBase
{
public:

    typedef FluidMecStressTensorImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ExprIdPressureType,ModelPhysicFluidType,SymbolsExprType,SpecificExprType> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_pressure = vm::JACOBIAN;
    static const size_type context = context_velocity|vm::DYNAMIC;

    using expr_id_pressure_type = ExprIdPressureType;
    using model_physic_fluid_type = ModelPhysicFluidType;
    using symbols_expr_type = SymbolsExprType;
    using specific_expr_type = SpecificExprType;
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorsType;
    using value_type = typename expr_evaluate_velocity_opertors_type::value_type;

    static const bool is_terminal = true;


    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<ExprApplyType::JACOBIAN> >,
                                            typename mpl::if_<boost::is_same<Func,FiniteElementVelocityType>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type,
                                            mpl::bool_<false> >::type::value;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    template <int M,int N>
    using material_property_expr_type = std::decay_t<decltype( expr( ModelExpression{}.template expr<M,N>(),symbols_expr_type{} ) )>;

    using material_property_scalar_expr_type = material_property_expr_type<1,1>;

    using expr_dynamic_viscosity_type = FluidMecDynamicViscosityImpl<expr_evaluate_velocity_opertors_type,FiniteElementVelocityType,ModelPhysicFluidType,SymbolsExprType,
                                                                     mpl::int_< this_type::specific_expr_type::value == ExprApplyType::EVAL? ExprApplyType::EVAL : ExprApplyType::JACOBIAN> >;

    FluidMecStressTensorImpl( std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvalVelocityOperators,
                              expr_id_pressure_type const& exprIdPressure,
                              model_physic_fluid_type const& physicFluid,
                              MaterialProperties const& matProps,
                              uint16_type polyOrder,
                              bool withPressure,
                              SymbolsExprType const& se)
        :
        M_exprEvaluateVelocityOperators( exprEvalVelocityOperators ),
        M_exprIdPressure( exprIdPressure ),
        M_physicFluid( physicFluid ),
        M_matProps( matProps ),
        M_polynomialOrder( polyOrder ),
        M_withPressureTerm( withPressure ),
        M_se( se )
    {
        M_exprEvaluateVelocityOperators->setEnableGrad( true );
        M_exprDynamicVisocity.emplace( M_exprEvaluateVelocityOperators,physicFluid,matProps,invalid_uint16_type_value,se );
        if ( this->turbulence().isEnabled() )
        {
            M_exprTurbulentDynamicVisocity.emplace( this->materialPropertyExpr<1,1>("turbulent-dynamic-viscosity") );
            if ( this->turbulence().hasTurbulentKineticEnergy() )
            {
                M_exprTurbulentKineticEnergy.emplace( this->materialPropertyExpr<1,1>("turbulent-kinetic-energy") );
                M_exprDensity.emplace( this->materialPropertyExpr<1,1>("density") );
            }
        }
    }

    FluidMecStressTensorImpl( FluidMecStressTensorImpl const & op ) = default;
    FluidMecStressTensorImpl( FluidMecStressTensorImpl && op ) = default;

    //~FluidMecStressTensorImpl() {}

    size_type dynamicContext() const
        {
            size_type res = M_exprDynamicVisocity->dynamicContext();
            if ( this->turbulence().isEnabled() )
            {
                res= res | Feel::vf::dynamicContext( *M_exprTurbulentDynamicVisocity );
                if ( M_exprTurbulentKineticEnergy )
                    res= res | Feel::vf::dynamicContext( *M_exprTurbulentKineticEnergy ) |  Feel::vf::dynamicContext( *M_exprDensity );
            }
            return res;
        }

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            if ( M_polynomialOrder != invalid_uint16_type_value )
                return M_polynomialOrder;

            uint16_type orderGradVelocity = M_exprEvaluateVelocityOperators->exprGrad().polynomialOrder();

            uint16_type res = M_exprDynamicVisocity->polynomialOrder() + orderGradVelocity; // TODO Jacobian case

            if ( this->turbulence().isEnabled() )
            {
                res = std::max( res, (uint16_type) (M_exprTurbulentDynamicVisocity->polynomialOrder() + orderGradVelocity) );
                //if ( M_exprTurbulentKineticEnergy )
                //TODO
            }

            if ( this->withPressureTerm() )
                res = std::max( res, M_exprIdPressure.polynomialOrder() );

            return res;
        }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }


    template <typename OtherSymbolsExprType>
    auto applySymbolsExpr( OtherSymbolsExprType const& se ) const
        {
            auto newse = Feel::vf::symbolsExpr( M_se, se );
            using new_se_type = std::decay_t<decltype(newse)>;
            using new_this_type = FluidMecStressTensorImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ExprIdPressureType,ModelPhysicFluidType,new_se_type,SpecificExprType>;
            return new_this_type( M_exprEvaluateVelocityOperators,M_exprIdPressure,M_physicFluid,M_matProps,M_polynomialOrder,M_withPressureTerm,newse );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            CHECK( false ) << "not implemented";
            return *this;
        }

    std::shared_ptr<expr_evaluate_velocity_opertors_type> exprEvaluateVelocityOperatorsPtr() const { return M_exprEvaluateVelocityOperators; }

    expr_id_pressure_type const& expressionIdPressure() const { return M_exprIdPressure; }
    expr_dynamic_viscosity_type const& exprDynamicVisocity() const { return *M_exprDynamicVisocity; }
    material_property_scalar_expr_type const& exprTurbulentDynamicVisocity() const { return *M_exprTurbulentDynamicVisocity; }
    material_property_scalar_expr_type const& exprTurbulentKineticEnergy() const { return *M_exprTurbulentKineticEnergy; }
    material_property_scalar_expr_type const& exprDensity() const { return *M_exprDensity; }


    auto const& turbulence() const { return M_physicFluid.turbulence(); }
    MaterialProperties const& materialProperties() const { return M_matProps; }
    bool withPressureTerm() const { return M_withPressureTerm; }


    template <int M,int N>
    material_property_expr_type<M,N> materialPropertyExpr( std::string const& prop ) const { return expr( this->materialProperties().property( prop ).template expr<M,N>(), M_se ); }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                      Shape<expr_evaluate_velocity_opertors_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape_grad_type::nDim,Tensor2, false, false>,
                                      typename this_type::value_type>
    {
        using super_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                      Shape<expr_evaluate_velocity_opertors_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape_grad_type::nDim,Tensor2, false, false>,
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
        using tensor_expr_dynamic_viscosity_type = typename expr_dynamic_viscosity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_id_pressure_type = typename expr_id_pressure_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_dynamicViscosityDependsOnVelocityField( false )
        {
            this->initTensor( expr, geom, fev, feu );
        }

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_dynamicViscosityDependsOnVelocityField( false )
        {
            this->initTensor( expr, geom, fev );
        }
        tensor( this_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
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
            M_tensorDynamicViscosity->update( geom, false );
            if ( M_tensorExprTurbulentDynamicVisocity )
                M_tensorExprTurbulentDynamicVisocity->update( geom );
            if ( M_tensorExprTurbulentKineticEnergy)
            {
                M_tensorExprTurbulentKineticEnergy->update( geom );
                M_tensorExprDensity->update( geom );
            }
            if ( M_tensorExprIdPressure )
                M_tensorExprIdPressure->update( geom );
            this->updateImpl();
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            this->setGmc( geom );
            M_tensorExprEvaluateVelocityOperators->update( geom, face );
            M_tensorDynamicViscosity->update( geom, face, false );
            if ( M_tensorExprTurbulentDynamicVisocity )
                M_tensorExprTurbulentDynamicVisocity->update( geom, face );
            if ( M_tensorExprTurbulentKineticEnergy)
            {
                M_tensorExprTurbulentKineticEnergy->update( geom, face );
                M_tensorExprDensity->update( geom, face );
            }
            if ( M_tensorExprIdPressure )
                M_tensorExprIdPressure->update( geom, face );
            this->updateImpl();
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            if constexpr ( this_type::specific_expr_type::value == ExprApplyType::EVAL )
                return this->evalq( q );
            else
            {
                value_type muEval = M_tensorDynamicViscosity->evalq(0,0,q);
                if ( M_tensorExprTurbulentDynamicVisocity )
                    muEval += M_tensorExprTurbulentDynamicVisocity->evalq(0,0,q);

                auto const& gradTrial = this->fecTrial()->grad( j, q );
                matrix_shape_type & locMat = this->locMatrixShape();
                if constexpr ( gmc_type::nDim == 2 )
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
                    locMat(0,0) = 2*muEval*du1tdx;
                    locMat(1,0) = muEval*( du2tdx + du1tdy );
                    locMat(1,1) = 2*muEval*du2tdy;
                    if ( M_dynamicViscosityDependsOnVelocityField )
                    {
                        auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                        value_type mut = M_tensorDynamicViscosity->evalijq(i,j,q)(0,0);
                        locMat(0,0) += 2*mut*du1vdx;
                        locMat(1,0) += mut*( du2vdx + du1vdy );
                        locMat(1,1) += 2*mut*du2vdy;
                    }
                    locMat(0,1) = locMat(1,0);
                }
                else
                {
                    const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
                    const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
                    const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
                    locMat(0,0) = 2*muEval*du1tdx;
                    locMat(1,0) = muEval*( du2tdx + du1tdy );
                    locMat(2,0) = muEval*( du3tdx + du1tdz );
                    locMat(1,1) = 2*muEval*du2tdy;
                    locMat(2,1) = muEval*( du2tdz + du3tdy );
                    locMat(2,2) = 2*muEval*du3tdz;
                    if ( M_dynamicViscosityDependsOnVelocityField )
                    {
                        auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                        const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                        value_type mut = M_tensorDynamicViscosity->evalijq(i,j,q)(0,0);
                        locMat(0,0) += 2*mut*du1vdx;
                        locMat(1,0) += mut*( du2vdx + du1vdy );
                        locMat(2,0) += mut*( du3vdx + du1vdz );
                        locMat(1,1) += 2*mut*du2vdy;
                        locMat(2,1) += mut*( du2vdz + du3vdy );
                        locMat(2,2) += 2*mut*du3vdz;
                    }
                    locMat(0,1) = locMat(1,0);
                    locMat(0,2) = locMat(2,0);
                    locMat(1,2) = locMat(2,1);
                }
                return ret_type(this->locMatrixShape().data());
            }
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( this_type::specific_expr_type::value == ExprApplyType::EVAL )
                return this->evalq( c1,c2,q );
            else
            {
                CHECK( false) << "TODO";
                return value_type(0);
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if constexpr ( this_type::specific_expr_type::value != ExprApplyType::EVAL )
                 CHECK( false ) << "not allow";
            return this->evalq( c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const override
        {
            if constexpr ( this_type::specific_expr_type::value != ExprApplyType::EVAL )
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
                M_tensorDynamicViscosity.emplace( expr.exprDynamicVisocity(), M_tensorExprEvaluateVelocityOperators, theInitArgs... );
                M_dynamicViscosityDependsOnVelocityField = expr.exprDynamicVisocity().dependsOnVelocityField();
                if ( expr.turbulence().isEnabled() )
                {
                    M_tensorExprTurbulentDynamicVisocity.emplace( expr.exprTurbulentDynamicVisocity(), theInitArgs... );
                    if ( expr.turbulence().hasTurbulentKineticEnergy() )
                    {
                        M_tensorExprTurbulentKineticEnergy.emplace( expr.exprTurbulentKineticEnergy(), theInitArgs... );
                        M_tensorExprDensity.emplace( expr.exprDensity(), theInitArgs... );
                    }
                }
                if ( expr.withPressureTerm() )
                    M_tensorExprIdPressure.emplace( expr.expressionIdPressure(),theInitArgs... );
            }

        void updateImpl()
            {
                if constexpr ( this_type::specific_expr_type::value != ExprApplyType::EVAL )
                                 return;

                uint16_type nPoints = this->gmc()->nPoints();
                if ( M_localEval.size() != nPoints )
                    M_localEval.resize( nPoints );

                bool withPressure = M_expr.withPressureTerm();
                bool hasTurbulentDynamicViscosity = ( M_tensorExprTurbulentDynamicVisocity )? true : false;
                bool hasTurbulentKineticEnergy = ( M_tensorExprTurbulentKineticEnergy)? true : false;

                for ( uint16_type q=0;q< nPoints;++q )
                {
                    value_type muEval = M_tensorDynamicViscosity->evalq(0,0,q);
                    if ( hasTurbulentDynamicViscosity )
                        muEval += M_tensorExprTurbulentDynamicVisocity->evalq(0,0,q);
                    matrix_shape_type & locMat = M_localEval[q];
                    auto const& gradVelocityEval = M_tensorExprEvaluateVelocityOperators->localEvalGrad( q );

                    if constexpr ( gmc_type::nDim == 2 )
                    {
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                        locMat(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                        locMat(1,0) = muEval*( du2vdx + du1vdy );
                        locMat(0,1) = locMat(1,0);
                        locMat(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    }
                    else
                    {
                        const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                        const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                        const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                        locMat(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                        locMat(1,0) = muEval*( du2vdx + du1vdy );
                        locMat(2,0) = muEval*( du3vdx + du1vdz );
                        locMat(0,1) = locMat(1,0);
                        locMat(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                        locMat(1,2) = muEval*( du2vdz + du3vdy );
                        locMat(2,1) = locMat(1,2);
                        locMat(0,2) = locMat(2,0);
                        locMat(2,2) = /*-pres +*/ 2*muEval*du3vdz;
                    }

                    if ( hasTurbulentKineticEnergy )
                    {
                        const value_type tke = M_tensorExprTurbulentKineticEnergy->evalq(0,0,q);
                        const value_type rho = M_tensorExprDensity->evalq(0,0,q);
                        const value_type tkeContrib = (2./3.)*rho*tke;
                        locMat(0,0) -= tkeContrib;
                        locMat(1,1) -= tkeContrib;
                        if constexpr ( gmc_type::nDim == 3 )
                        {
                             locMat(2,2) -= tkeContrib;
                        }
                    }

                    if ( withPressure )
                    {
                        const value_type pres = M_tensorExprIdPressure->evalq( 0,0,q );
                        locMat(0,0) -= pres;
                        locMat(1,1) -= pres;
                        if constexpr ( gmc_type::nDim == 3 )
                        {
                             locMat(2,2) -= pres;
                        }
                    }

                }
            }


    private:
        this_type const& M_expr;

        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> M_tensorExprEvaluateVelocityOperators;
        std::optional<tensor_expr_dynamic_viscosity_type> M_tensorDynamicViscosity;
        std::optional<tensor_expr_id_pressure_type> M_tensorExprIdPressure;
        std::optional<tensor_expr_material_property_scalar_type> M_tensorExprTurbulentDynamicVisocity, M_tensorExprTurbulentKineticEnergy, M_tensorExprDensity;

        array_shape_type M_localEvalStrainTensor;
        array_shape_type M_localEval;

        bool M_dynamicViscosityDependsOnVelocityField;
    };

private:
    std::shared_ptr<expr_evaluate_velocity_opertors_type> M_exprEvaluateVelocityOperators;
    expr_id_pressure_type M_exprIdPressure;


    model_physic_fluid_type const& M_physicFluid;
    MaterialProperties const& M_matProps;
    uint16_type M_polynomialOrder;
    bool M_withPressureTerm;
    SymbolsExprType M_se;

    std::optional<expr_dynamic_viscosity_type> M_exprDynamicVisocity;
    std::optional<material_property_scalar_expr_type> M_exprTurbulentDynamicVisocity, M_exprTurbulentKineticEnergy, M_exprDensity;
};


template<class ExprGradVelocityType,class ExprIdPressureType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecStressTensor( Expr<ExprGradVelocityType> const& grad_u, Expr<ExprIdPressureType> const& p,
                      ModelPhysicFluidType const& physicFluid,
                      MaterialProperties const& matProps,
                      bool withPressure,
                      SymbolsExprType const& se = symbols_expression_empty_t{},
                      uint16_type polyOrder = invalid_uint16_type_value )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorGradFromExpr<Expr<ExprGradVelocityType>>;
    typedef FluidMecStressTensorImpl<expr_evaluate_velocity_opertors_type,std::nullptr_t,Expr<ExprIdPressureType>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<ExprApplyType::EVAL> > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( grad_u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,p,physicFluid,matProps,polyOrder,withPressure,se ) );
}

template<class VelocityFieldType,class ExprIdPressureType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecStressTensor( VelocityFieldType const& u, Expr<ExprIdPressureType> const& p,
                      ModelPhysicFluidType const& physicFluid,
                      MaterialProperties const& matProps,
                      bool withPressure,
                      SymbolsExprType const& se = symbols_expression_empty_t{},
                      uint16_type polyOrder = invalid_uint16_type_value,
                      typename std::enable_if_t<is_functionspace_element_v<unwrap_ptr_t<VelocityFieldType>> >* = nullptr
                      )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperators<VelocityFieldType>;
    typedef FluidMecStressTensorImpl<expr_evaluate_velocity_opertors_type,std::nullptr_t,Expr<ExprIdPressureType>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<ExprApplyType::EVAL> > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,p,physicFluid,matProps,polyOrder,withPressure,se ) );
}

template<class ExprGradVelocityType,class ElementVelocityType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecViscousStressTensorJacobian( Expr<ExprGradVelocityType> const& grad_u, ElementVelocityType const& /*u*/,
                                     ModelPhysicFluidType const& physicFluid,
                                     MaterialProperties const& matProps,
                                     SymbolsExprType const& se = symbols_expression_empty_t{},
                                     uint16_type polyOrder = invalid_uint16_type_value )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorGradFromExpr<Expr<ExprGradVelocityType>>;
    typedef FluidMecStressTensorImpl<expr_evaluate_velocity_opertors_type,typename unwrap_ptr_t<ElementVelocityType>::functionspace_type::reference_element_type,Expr<Cst<double>>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<ExprApplyType::JACOBIAN> > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( grad_u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,/*p*/cst(0.),physicFluid,matProps,polyOrder,false,se ) );
}

template<class VelocityFieldType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecViscousStressTensorJacobian( VelocityFieldType const& u,
                                     ModelPhysicFluidType const& physicFluid,
                                     MaterialProperties const& matProps,
                                     SymbolsExprType const& se = symbols_expression_empty_t{},
                                     uint16_type polyOrder = invalid_uint16_type_value )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperators<VelocityFieldType>;
    typedef FluidMecStressTensorImpl<expr_evaluate_velocity_opertors_type,typename unwrap_ptr_t<VelocityFieldType>::functionspace_type::reference_element_type,Expr<Cst<double>>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<ExprApplyType::JACOBIAN> > fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,/*p*/cst(0.),physicFluid,matProps,polyOrder,false,se ) );
}

} // namespace FeelModels
} // namespace Feel
#endif /* FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H */
