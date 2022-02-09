/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_MODELS_VF_SHOCKCAPTURING_EXPR_H
#define FEELPP_MODELS_VF_SHOCKCAPTURING_EXPR_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{


template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ShapeType, typename ExprType>
struct tensorShockCapturingGaleaoCarmo : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type >
{
    typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type> super_type;
    //typedef tensorStabGLSParameterEigenValue<Geo_t,Basis_i_t,Basis_j_t,ShapeType,ExprType,HasConvectionExpr,HasDiffusionExpr,HasReactionExpr,VariantType> self_type;
    typedef ExprType expr_type;

    typedef typename ExprType::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;

    typedef typename expr_type::expression_residual_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_residual_type;
    typedef typename expr_type::expression_residual_derivative_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_residual_derivative_type;
    typedef typename expr_type::expression_tau_parameter_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_tau_parameter_type;
    typedef typename expr_type::expression_convection_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_convection_type;
    typedef typename expr_type::expression_gradv_u_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_gradv_u_type;

    tensorShockCapturingGaleaoCarmo( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( geom, fev, feu ),
        M_expr( expr ),
        M_tensorResidual( tensor_residual_type( expr.exprResidual(), geom, fev, feu ) ),
        M_tensorTauParameter( tensor_tau_parameter_type( expr.exprTauParameter(), geom, fev, feu ) ),
        M_tensorConvection( tensor_convection_type( expr.exprConvection(), geom, fev, feu ) ),
        M_tensorGradv_u( tensor_gradv_u_type( expr.exprGradv_u(), geom, fev, feu ) )
        {
            if constexpr ( ExprType::is_applied_in_jacobian )
                M_tensorResidualDerivative.emplace( tensor_residual_derivative_type( expr.exprResidualDerivative(), geom, fev, feu ) );
        }

    tensorShockCapturingGaleaoCarmo( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( geom, fev ),
        M_expr( expr ),
        M_tensorResidual( tensor_residual_type( expr.exprResidual(), geom, fev ) ),
        M_tensorTauParameter( tensor_tau_parameter_type( expr.exprTauParameter(), geom, fev ) ),
        M_tensorConvection( tensor_convection_type( expr.exprConvection(), geom, fev ) ),
        M_tensorGradv_u( tensor_gradv_u_type( expr.exprGradv_u(), geom, fev ) )
        {}

    tensorShockCapturingGaleaoCarmo( expr_type const& expr, Geo_t const& geom )
        :
        super_type( geom ),
        M_expr( expr ),
        M_tensorResidual( tensor_residual_type( expr.exprResidual(), geom ) ),
        M_tensorTauParameter( tensor_tau_parameter_type( expr.exprTauParameter(), geom ) ),
        M_tensorConvection( tensor_convection_type( expr.exprConvection(), geom ) ),
        M_tensorGradv_u( tensor_gradv_u_type( expr.exprGradv_u(), geom ) )
        {}

    void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) override
        {
            if constexpr ( !ExprType::is_applied_in_jacobian )
                this->update( geom );

            this->setGmc( geom );
            M_tensorResidual.update( geom );
            M_tensorTauParameter.update( geom );
            M_tensorConvection.update( geom );
            M_tensorGradv_u.update( geom );
            M_tensorResidualDerivative->update( geom, fev, feu );

            uint16_type nPoints = this->gmc()->nPoints();

            M_localMatrixIsZero.resize( nPoints );
            std::fill( M_localMatrixIsZero.begin(), M_localMatrixIsZero.end(), false );

            M_localAssemblyJacobian_partB1a.resize( nPoints );
            M_localAssemblyJacobian_partB2.resize( nPoints );
            M_localAssemblyJacobian_partB1b.resize( nPoints );
            M_localAssemblyZh.resize( nPoints );
            M_localAssemblyJacobianZhPrime_part1.resize( nPoints );
            M_localAssemblyJacobianZhPrime_part2.resize( nPoints );
            M_localAssemblyJacobianZhPrime_part3.resize( nPoints );
            M_localAssemblyJacobianTauscPrime_part1.resize( nPoints );
            M_localAssemblyJacobianTermWithoutTausc.resize( nPoints );

            typename super_type::matrix_shape_type zh;
            for ( uint16_type q = 0; q < nPoints; ++q )
            {
                auto evalGradv_u = M_tensorGradv_u.evalq( q );
                value_type norm_gradv_u_square = evalGradv_u.squaredNorm();
                if ( norm_gradv_u_square < 1e-10 )
                {
                    M_localMatrixIsZero[q] = true;
                    continue;
                }

                value_type residualValue = M_tensorResidual.evalq(0,0,q);
                zh = residualValue*(evalGradv_u.transpose()/norm_gradv_u_square);

                value_type norm_zh = zh.norm();

                if ( norm_zh < 1e-14 )
                {
                    M_localMatrixIsZero[q] = true;
                    continue;
                }

                value_type norm_convection_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_convection_square += math::pow(M_tensorConvection.evalq( c,0,q),2);
                value_type norm_convection = math::sqrt( norm_convection_square );

                value_type zeta = 1;
                if ( true )
                {
                    //value_type tau_sc = M_tensorTauParameter.evalq(0,0,q)*std::max( 0., norm_convection/norm_zh -1.0);
                }
                else
                {
                    value_type beta_dot_gradv_u = 0;
                    for ( uint16_type c=0;c<tensor_convection_type::shape::M;++c )
                        beta_dot_gradv_u += M_tensorConvection.evalq( c,0,q)*evalGradv_u(0,c);

                    if ( residualValue > 1e-10 )
                        zeta = std::max( 1., beta_dot_gradv_u/residualValue );
                }

                //value_type tau_sc = M_tensorTauParameter.evalq(0,0,q)*std::max( 0., norm_convection/norm_zh - zeta);
                value_type tau_sc = 0;
                value_type tauScaling = norm_convection/norm_zh - zeta;
                if ( tauScaling > 0 )
                    tau_sc = M_tensorTauParameter.evalq(0,0,q)*tauScaling;
                else
                {
                    M_localMatrixIsZero[q] = true;
                    continue;
                }

                M_localAssemblyJacobian_partB1a[q] = (tau_sc*2*residualValue/norm_gradv_u_square)*evalGradv_u.transpose();
                M_localAssemblyJacobian_partB1b[q] = tau_sc*math::pow(residualValue,2)/norm_gradv_u_square;
                M_localAssemblyJacobian_partB2[q] = (tau_sc*math::pow(residualValue,2)/math::pow(norm_gradv_u_square,2))*evalGradv_u.transpose();

                M_localAssemblyZh[q] = zh;

                M_localAssemblyJacobianZhPrime_part1[q] = (1./norm_gradv_u_square)*evalGradv_u.transpose();
                M_localAssemblyJacobianZhPrime_part2[q] = residualValue/norm_gradv_u_square;
                M_localAssemblyJacobianZhPrime_part3[q] = (residualValue/math::pow(norm_gradv_u_square,2))*evalGradv_u.transpose();

                M_localAssemblyJacobianTauscPrime_part1[q] = M_tensorTauParameter.evalq(0,0,q)*norm_convection/math::pow(norm_zh,3);
                M_localAssemblyJacobianTermWithoutTausc[q] = residualValue*zh;

            }

        }

    void update( Geo_t const& geom ) override
        {
            this->setGmc( geom );
            M_tensorResidual.update( geom );
            M_tensorTauParameter.update( geom );
            M_tensorConvection.update( geom );
            M_tensorGradv_u.update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );

            typename super_type::matrix_shape_type zh;
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                auto evalGradv_u = M_tensorGradv_u.evalq( q );
                value_type norm_gradv_u_square = evalGradv_u.squaredNorm();
                if ( norm_gradv_u_square < 1e-10 )
                {
                    M_localMatrixInGeoContext[q].setZero();
                    continue;
                }

                value_type residualValue = M_tensorResidual.evalq(0,0,q);
                zh = residualValue*(evalGradv_u.transpose()/norm_gradv_u_square);

                value_type norm_zh = zh.norm();

                if ( norm_zh < 1e-14 )
                {
                    M_localMatrixInGeoContext[q].setZero();
                    continue;
                }

                value_type norm_convection_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_convection_square += math::pow(M_tensorConvection.evalq( c,0,q),2);
                value_type norm_convection = math::sqrt( norm_convection_square );

                value_type zeta = 1;
                if ( true )
                {
                    //value_type tau_sc = M_tensorTauParameter.evalq(0,0,q)*std::max( 0., norm_convection/norm_zh -1.0);
                }
                else
                {
                    value_type beta_dot_gradv_u = 0;
                    for ( uint16_type c=0;c<tensor_convection_type::shape::M;++c )
                        beta_dot_gradv_u += M_tensorConvection.evalq( c,0,q)*evalGradv_u(0,c);

                    if ( residualValue > 1e-10 )
                        zeta = std::max( 1., beta_dot_gradv_u/residualValue );
                }

                value_type tauScaling = norm_convection/norm_zh - zeta;
                if ( tauScaling > 0 )
                {
                    value_type tau_sc = M_tensorTauParameter.evalq(0,0,q)*tauScaling;
                    M_localMatrixInGeoContext[q] = tau_sc*residualValue*zh;
                }
                else
                    M_localMatrixInGeoContext[q].setZero();
            }

        }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            if ( M_localMatrixIsZero[q] )
                return value_type(0);

            auto const& gradTrial = Feel::vf::detail::convertEigenMatrixTensor( this->fecTrial()->grad( j, q ) );

            auto gradv_u = M_tensorGradv_u.evalq( q );
            value_type norm_gradv_u_square_prime = 0;
            for (uint16_type c=0;c<tensor_gradv_u_type::shape::N;++c)
                norm_gradv_u_square_prime += gradv_u(0,c)*gradTrial(0,c/*,0*/);
            norm_gradv_u_square_prime *= 2;

            value_type residualDerivativeValue = M_tensorResidualDerivative->evalijq(i,j,c1,c2,q);

            value_type res =
                M_localAssemblyJacobian_partB1a[q](c1,c2)*residualDerivativeValue
                + M_localAssemblyJacobian_partB1b[q]*gradTrial(c2,c1/*,0*/) /* transpose grad*/
                - M_localAssemblyJacobian_partB2[q](c1,c2)*norm_gradv_u_square_prime;

            if ( true )
            {
                auto zhprime =
                    residualDerivativeValue*M_localAssemblyJacobianZhPrime_part1[q]
                    + M_localAssemblyJacobianZhPrime_part2[q]*gradTrial.transpose()
                    - norm_gradv_u_square_prime*M_localAssemblyJacobianZhPrime_part3[q];

                auto const& zh = M_localAssemblyZh[q];

                // the derivative of this term has been simplified but still exact
                // warning not include the division by norm_zh, do with M_localAssemblyJacobianTauscPrime_part1
                value_type partial_norm_zh_prime = zh.adjoint()*zhprime;

                value_type tau_sc_prime = -M_localAssemblyJacobianTauscPrime_part1[q]*partial_norm_zh_prime;

                res += tau_sc_prime*M_localAssemblyJacobianTermWithoutTausc[q](c1,c2);
            }
             return res;
        }

    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            return this->evalq( c1,c2,q );
        }

    value_type
    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
        {
            return M_localMatrixInGeoContext[q](c1,c2);
        }
#if 0
    ret_type
    evalq( uint16_type q ) const
        {
            return ret_type(M_localMatrixInGeoContext[q].data());
        }
#endif

private :
    expr_type const& M_expr;
    tensor_residual_type M_tensorResidual;
    std::optional<tensor_residual_derivative_type> M_tensorResidualDerivative;
    tensor_tau_parameter_type M_tensorTauParameter;
    tensor_convection_type M_tensorConvection;
    tensor_gradv_u_type M_tensorGradv_u;
    std::vector<bool> M_localMatrixIsZero;
    std::vector<matrix_shape_type> M_localMatrixInGeoContext;

    std::vector<matrix_shape_type> M_localAssemblyJacobian_partB1a, M_localAssemblyJacobian_partB2;
    std::vector<value_type>  M_localAssemblyJacobian_partB1b;

    std::vector<matrix_shape_type> M_localAssemblyZh, M_localAssemblyJacobianZhPrime_part1, M_localAssemblyJacobianZhPrime_part3;
    std::vector<value_type> M_localAssemblyJacobianZhPrime_part2;

    std::vector<value_type> M_localAssemblyJacobianTauscPrime_part1;
    std::vector<matrix_shape_type> M_localAssemblyJacobianTermWithoutTausc;
};


template<typename ExprResidualType, typename ExprResidualDerivativeType, typename ExprTauParameterType, typename ExprConvectionType,typename ElementType, ExprApplyType TheExprType>
class ShockCapturingExpr : public Feel::vf::ExprDynamicBase
{
public:

    typedef ShockCapturingExpr<ExprResidualType,ExprResidualDerivativeType,ExprTauParameterType,ExprConvectionType,ElementType,TheExprType> this_type;
    typedef ExprResidualType expression_residual_type;
    typedef ExprResidualDerivativeType expression_residual_derivative_type;
    typedef ExprTauParameterType expression_tau_parameter_type;
    typedef ExprConvectionType expression_convection_type;

    static constexpr bool is_applied_in_jacobian = (TheExprType == ExprApplyType::JACOBIAN);
    using expression_gradv_u_type = std::decay_t<decltype(Feel::vf::gradv(ElementType{}))>;
    typedef ElementType element_type;

    static const size_type context_base = vm::DYNAMIC|expression_residual_type::context|expression_tau_parameter_type::context|expression_convection_type::context|expression_gradv_u_type::context;
    static const size_type context = mpl::if_c<is_applied_in_jacobian,
                                               mpl::integral_c<size_type,context_base|expression_residual_derivative_type::context>,
                                               mpl::integral_c<size_type,context_base> >::type::value;
    static const bool is_terminal = false; // TODO optimize
    typedef double value_type;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = is_applied_in_jacobian &&  std::is_same_v<Func,typename unwrap_ptr_t<element_type>::functionspace_type::reference_element_type>;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    ShockCapturingExpr( expression_residual_type const& exprResidual,
                        expression_tau_parameter_type const& exprTauParameter,
                        expression_convection_type const& exprConvection,
                        element_type const& u )
        :
        M_exprResidual( exprResidual ),
        M_exprTauParameter( exprTauParameter ),
        M_exprConvection( exprConvection ),
        M_u( u ),
        M_exprGradv_u( Feel::vf::gradv(M_u) )
        {}

    ShockCapturingExpr( expression_residual_type const& exprResidual,
                        expression_residual_derivative_type const& exprResidualDerivative,
                        expression_tau_parameter_type const& exprTauParameter,
                        expression_convection_type const& exprConvection,
                        element_type const& u )
        :
        M_exprResidual( exprResidual ),
        M_exprResidualDerivative( exprResidualDerivative ),
        M_exprTauParameter( exprTauParameter ),
        M_exprConvection( exprConvection ),
        M_u( u ),
        M_exprGradv_u( Feel::vf::gradv(M_u) )
        {}

    expression_residual_type const& exprResidual() const { return M_exprResidual; }
    expression_residual_derivative_type const& exprResidualDerivative() const { CHECK( this->hasResidualDerivative() ) << "no residual derivative"; return *M_exprResidualDerivative; }
    expression_tau_parameter_type const& exprTauParameter() const { return M_exprTauParameter; }
    expression_convection_type const& exprConvection() const { return M_exprConvection; }
    element_type const& u() const { return M_u; }
    expression_gradv_u_type const& exprGradv_u() const { return M_exprGradv_u; }

    bool hasResidualDerivative() const { return M_exprResidualDerivative.has_value(); }

    size_type dynamicContext() const
        {
            size_type res = 0;
            res = res | Feel::vf::dynamicContext( this->exprResidual() );
            res = res | Feel::vf::dynamicContext( this->exprTauParameter() );
            res = res | Feel::vf::dynamicContext( this->exprConvection() );
            res = res | Feel::vf::dynamicContext( this->exprGradv_u() );
            if ( this->hasResidualDerivative() )
                res = res | Feel::vf::dynamicContext( this->exprResidualDerivative() );
            return res;
        }

    //! polynomial order
    uint16_type polynomialOrder() const
    {
        return this->exprResidual().polynomialOrder()+1; // TODO
    }

    //! expression is polynomial?
    bool isPolynomial() const { return false; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef tensor<Geo_t,Basis_i_t,Basis_j_t> self_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,Feel::vf::detail::gmc<0> >, \
                                  mpl::identity<Feel::vf::detail::gmc<0> >, \
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type key_type; \
        //typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename gmc_type::value_type value_type;
        //typedef  value_type evaluate_type;
        typedef Shape<gmc_type::NDim, Vectorial/*Scalar*/, false,false> shape;
        struct is_zero { static const bool value = false; };

        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t,shape,value_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename tensorbase_type::ret_type;
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            {
                M_tensorbase.reset( new tensorShockCapturingGaleaoCarmo<Geo_t, Basis_i_t, Basis_j_t, shape, this_type>(expr,geom,fev,feu) );
            }
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            {
                M_tensorbase.reset( new tensorShockCapturingGaleaoCarmo<Geo_t, Basis_i_t, Basis_j_t, shape, this_type>(expr,geom,fev) );
            }
        tensor( this_type const& expr, Geo_t const& geom )
            {
                M_tensorbase.reset( new tensorShockCapturingGaleaoCarmo<Geo_t, Basis_i_t, Basis_j_t, shape, this_type>(expr,geom) );
            }
        tensor( tensor const& t ) = default;

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensorbase->update( geom, fev, feu );
            //M_tensorbase->update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            //M_tensorbase->update( geom, fev );
            M_tensorbase->update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_tensorbase->update( geom );
        }
#if 0
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }
#endif
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
#if 0
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }
#endif
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
#if 0
        ret_type
        evalq( uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }
#endif
    private:
        tensorbase_ptrtype M_tensorbase;
    };

private :
    expression_residual_type M_exprResidual;
    std::optional<expression_residual_derivative_type> M_exprResidualDerivative;
    expression_tau_parameter_type M_exprTauParameter;
    expression_convection_type M_exprConvection;
    element_type const& M_u;
    expression_gradv_u_type M_exprGradv_u;
};

template<typename ExprResidualType, typename ExprTauParameterType, typename ExprConvectionType,typename ElementType>
inline
Expr< ShockCapturingExpr<ExprResidualType,ExprResidualType,ExprTauParameterType,ExprConvectionType, ElementType, ExprApplyType::EVAL> >
shockCapturingExpr( ExprResidualType const& exprResidual, ExprTauParameterType const& exprTau, ExprConvectionType const& exprConvection, ElementType const& u )
{
    typedef ShockCapturingExpr<ExprResidualType,ExprResidualType,ExprTauParameterType,ExprConvectionType, ElementType, ExprApplyType::EVAL > expr_sc_t;
    return Expr< expr_sc_t >( expr_sc_t(exprResidual,exprTau,exprConvection,u) );
}

template<typename ExprResidualType, typename ExprResidualDerivativeType, typename ExprTauParameterType, typename ExprConvectionType,typename ElementType>
inline
Expr< ShockCapturingExpr<ExprResidualType,ExprResidualDerivativeType,ExprTauParameterType,ExprConvectionType, ElementType, ExprApplyType::JACOBIAN> >
shockCapturingJacobianExpr( ExprResidualType const& exprResidual, ExprResidualDerivativeType const& exprResidualDerivative, ExprTauParameterType const& exprTau, ExprConvectionType const& exprConvection, ElementType const& u )
{
    typedef ShockCapturingExpr<ExprResidualType,ExprResidualDerivativeType,ExprTauParameterType,ExprConvectionType, ElementType, ExprApplyType::JACOBIAN > expr_sc_t;
    return Expr< expr_sc_t >( expr_sc_t(exprResidual,exprResidualDerivative,exprTau,exprConvection,u) );
}


} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_VF_STABILIZATION_GLS_PARMATER_EXPR_H
