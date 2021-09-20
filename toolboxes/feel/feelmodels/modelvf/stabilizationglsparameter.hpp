/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_MODELS_VF_STABILIZATION_GLS_PARMATER_EXPR_H
#define FEELPP_MODELS_VF_STABILIZATION_GLS_PARMATER_EXPR_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{

namespace detail
{
template <typename ExprTensorDiffusion>
typename ExprTensorDiffusion::value_type
tensorStabGLSParameter_evaluateDiffusion( ExprTensorDiffusion const& tensorDiffusion, uint16_type q,
                                               typename std::enable_if_t< ExprTensorDiffusion::shape::is_scalar >* = nullptr  )
{
    return tensorDiffusion.evalq(0,0,q);
}

template <typename ExprTensorDiffusion>
typename ExprTensorDiffusion::value_type
tensorStabGLSParameter_evaluateDiffusion( ExprTensorDiffusion const& tensorDiffusion, uint16_type q,
                                               typename std::enable_if_t< ExprTensorDiffusion::shape::is_tensor2 >* = nullptr  )
{
    // take the min of diag (TODO biblio, maybe min of eigen value or a norm)
    typename ExprTensorDiffusion::value_type res = tensorDiffusion.evalq(0,0,q);
    for ( int c1 = 1; c1 < ExprTensorDiffusion::shape::M ;++c1 )
        return res = std::min( res, tensorDiffusion.evalq(c1,c1,q) );
    return res;
}

} // namespace detail


template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ShapeType, typename ExprType, bool HasConvectionExpr, bool HasDiffusionExpr,bool HasReactionExpr,int VariantType>
struct tensorStabGLSParameterEigenValue : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type >
{
    typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type> super_type;
    typedef tensorStabGLSParameterEigenValue<Geo_t,Basis_i_t,Basis_j_t,ShapeType,ExprType,HasConvectionExpr,HasDiffusionExpr,HasReactionExpr,VariantType> self_type;
    typedef ExprType expr_type;
    static const bool has_convection_expr = HasConvectionExpr;
    static const bool has_diffusion_expr = HasDiffusionExpr;
    static const bool has_reaction_expr = HasReactionExpr;


    typedef typename ExprType::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;

    typedef typename expr_type::expression_convection_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_convection_type;
    typedef typename expr_type::expression_diffusion_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_diffusion_type;
    typedef typename expr_type::expression_reaction_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_reaction_type;

    tensorStabGLSParameterEigenValue( expr_type const& expr, Geo_t const& geom )
        :
        super_type( geom ),
        M_expr( expr )
        {
            if constexpr ( HasConvectionExpr )
                M_tensorConvection.emplace( tensor_convection_type( expr.exprConvection(), geom ) );
            if constexpr ( HasDiffusionExpr )
                M_tensorDiffusion.emplace( tensor_diffusion_type( expr.exprDiffusion(), geom ) );
            if constexpr ( HasReactionExpr )
                M_tensorReaction.emplace( tensor_reaction_type( expr.exprReaction(), geom ) );
        }

    void update( Geo_t const& geom )
        {
            this->update( geom, mpl::bool_<has_convection_expr>(), mpl::bool_<has_diffusion_expr>() );
        }
    void update( Geo_t const& geom, uint16_type face ) { CHECK( false ) << "TODO"; }

    value_type
    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_localMatrixInGeoContext[q](0,0);
        }
    ret_type
    evalq( uint16_type q ) const
        {
            return ret_type(M_localMatrixInGeoContext[q].data());
        }
private :
    void update( Geo_t const& geom, mpl::true_, mpl::true_ )
        {
            this->setGmc( geom );
            M_tensorConvection->update( geom );
            M_tensorDiffusion->update( geom );
            if constexpr ( has_reaction_expr )
                 M_tensorReaction->update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            value_type hSize2 = math::pow(hSize,2);
            value_type mK = (VariantType==0)? M_expr.stabGLSParameter().mK( eltId ) : value_type(1./3.);

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                value_type norm_u_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_u_square += math::pow(M_tensorConvection->evalq( c,0,q),2);

                value_type norm_u = math::sqrt( norm_u_square );
                //value_type evalDiffusion = M_tensorDiffusion.evalq(0,0,q);
                value_type evalDiffusion = Feel::FeelModels::detail::tensorStabGLSParameter_evaluateDiffusion( *M_tensorDiffusion, q );
                value_type Re = mK*norm_u*hSize/(2*evalDiffusion);
                if constexpr ( has_reaction_expr )
                {
                    value_type evalReaction = M_tensorReaction->evalq( 0,0,q );
                    value_type alpha_s = 2*evalDiffusion/(mK*hSize2*evalReaction);
                    value_type partReaction = (alpha_s < 1.)? evalReaction*hSize2 : 2*evalDiffusion/mK;
                    value_type partConvection = ( Re < 1. )? 2*evalDiffusion/mK : 2*hSize*norm_u;
                    M_localMatrixInGeoContext[q](0,0) = hSize2/(partReaction+partConvection);
                }
                else
                {
                    if ( Re < 1. )
                        M_localMatrixInGeoContext[q](0,0) = hSize2*mK/(4*evalDiffusion);
                    else
                        M_localMatrixInGeoContext[q](0,0) = hSize/(2*norm_u);
                }
            }
        }
    void update( Geo_t const& geom, mpl::false_, mpl::true_ )
        {
            this->setGmc( geom );
            M_tensorDiffusion->update( geom );
            if constexpr ( has_reaction_expr )
                 M_tensorReaction->update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            value_type mK = (VariantType==0)? M_expr.stabGLSParameter().mK( eltId ) : value_type(1./3.);
            value_type hSize2 = math::pow(hSize,2);
            value_type mK_times_hSize2 = mK*hSize2;
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                //value_type evalDiffusion = M_tensorDiffusion.evalq(0,0,q);
                value_type evalDiffusion = Feel::FeelModels::detail::tensorStabGLSParameter_evaluateDiffusion( *M_tensorDiffusion, q );

                if constexpr ( has_reaction_expr )
                {
                    value_type evalReaction = M_tensorReaction->evalq( 0,0,q );
                    value_type alpha_s = 2*evalDiffusion/(mK_times_hSize2*evalReaction);
                    if ( alpha_s < 1. )
                        M_localMatrixInGeoContext[q](0,0) = hSize2/(evalReaction*hSize2 + (2*evalDiffusion)/mK);
                    else
                        M_localMatrixInGeoContext[q](0,0) = mK_times_hSize2/(4*evalDiffusion);
                }
                else
                {
                    M_localMatrixInGeoContext[q](0,0) = mK_times_hSize2/(4*evalDiffusion);
                }
            }
        }
    void update( Geo_t const& geom, mpl::true_, mpl::false_ )
        {
            this->setGmc( geom );
            M_tensorConvection->update( geom );
            if constexpr ( has_reaction_expr )
                 M_tensorReaction->update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            value_type hSize2 = math::pow(hSize,2);
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                value_type norm_u_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_u_square += math::pow(M_tensorConvection->evalq( c,0,q),2);
                value_type norm_u = math::sqrt( norm_u_square );
                if constexpr ( has_reaction_expr )
                {
                    value_type evalReaction = M_tensorReaction->evalq( 0,0,q );
                    value_type partReactionConvection = evalReaction*hSize2 + 2*hSize*norm_u;
                    if ( math::abs( partReactionConvection ) > 1e-10 )
                        M_localMatrixInGeoContext[q](0,0) = hSize2/partReactionConvection;
                    else
                        M_localMatrixInGeoContext[q](0,0) = 0.;
                }
                else
                {
                    if ( math::abs( norm_u ) > 1e-10 )
                        M_localMatrixInGeoContext[q](0,0) = hSize/(2*norm_u);
                    else
                        M_localMatrixInGeoContext[q](0,0) = 0.;
                }
            }
        }

    void update( Geo_t const& geom, mpl::false_, mpl::false_ )
        {}

private :
    expr_type const& M_expr;
    std::optional<tensor_convection_type> M_tensorConvection;
    std::optional<tensor_diffusion_type> M_tensorDiffusion;
    std::optional<tensor_reaction_type> M_tensorReaction;
    std::vector<matrix_shape_type> M_localMatrixInGeoContext;
};
#if 0
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ShapeType, typename ExprType, bool HasConvectionExpr, bool HasDiffusionExpr,bool HasReactionExpr>
struct tensorStabGLSParameterDoublyAsymptoticApproximation : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type >
{
    typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,ShapeType,typename ExprType::value_type> super_type;
    typedef tensorStabGLSParameterDoublyAsymptoticApproximation<Geo_t,Basis_i_t,Basis_j_t,ShapeType,ExprType,HasConvectionExpr,HasDiffusionExpr,HasReactionExpr> self_type;
    typedef ExprType expr_type;
    static const bool has_convection_expr = HasConvectionExpr;
    static const bool has_diffusion_expr = HasDiffusionExpr;
    static const bool has_reaction_expr = HasReactionExpr;

    typedef typename ExprType::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;
    typedef typename expr_type::expression_convection_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_convection_type;
    typedef typename expr_type::expression_diffusion_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_diffusion_type;
    typedef typename expr_type::expression_reaction_type::template tensor<Geo_t,Basis_i_t,Basis_j_t> tensor_reaction_type;

    tensorStabGLSParameterDoublyAsymptoticApproximation( expr_type const& expr, Geo_t const& geom )
        :
        super_type( geom ),
        M_expr( expr ),
        M_tensorConvection( expr.exprConvection(),geom ),
        M_tensorDiffusion( expr.exprDiffusion(),geom )
        {}

    void update( Geo_t const& geom )
        {
            this->update( geom, mpl::bool_<has_convection_expr>(), mpl::bool_<has_diffusion_expr>() );
        }
    void update( Geo_t const& geom, uint16_type face ) { CHECK( false ) << "TODO"; }

    value_type
    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_localMatrixInGeoContext[q](0,0);
        }
    ret_type
    evalq( uint16_type q ) const
        {
            return ret_type(M_localMatrixInGeoContext[q].data());
        }
private :
    void update( Geo_t const& geom, mpl::true_, mpl::true_ )
        {
            this->setGmc( geom );
            M_tensorConvection.update( geom );
            M_tensorDiffusion.update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            value_type mK = (1./3.);
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                value_type norm_u_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_u_square += math::pow(M_tensorConvection.evalq( c,0,q),2);
                value_type norm_u = math::sqrt( norm_u_square );
                //value_type evalDiffusion = M_tensorDiffusion.evalq(0,0,q);
                value_type evalDiffusion = Feel::FeelModels::detail::tensorStabGLSParameter_evaluateDiffusion( M_tensorDiffusion, q );
#if 0
                value_type xi_Re = std::min(1.,(1./3.)*hSize*norm_u/(2*evalDiffusion));
                M_localMatrixInGeoContext[q](0,0) = hSize*xi_Re/( std::max( 1e-10, 2*norm_u) );
#else
                value_type Re = mK*hSize*norm_u/(2*evalDiffusion);
                if ( Re < 1 )
                    M_localMatrixInGeoContext[q](0,0) = math::pow(hSize,2)*mK/(4*evalDiffusion);
                else
                    M_localMatrixInGeoContext[q](0,0) = hSize/(2*norm_u);
#endif
            }
        }
    void update( Geo_t const& geom, mpl::false_, mpl::true_ )
        {
            this->setGmc( geom );
            M_tensorDiffusion.update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            value_type mK = (1./3.);
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                //value_type evalDiffusion = M_tensorDiffusion.evalq(0,0,q);
                value_type evalDiffusion = Feel::FeelModels::detail::tensorStabGLSParameter_evaluateDiffusion( M_tensorDiffusion, q );
                M_localMatrixInGeoContext[q](0,0) = math::pow(hSize,2)*mK/(4*evalDiffusion);
            }
        }
    void update( Geo_t const& geom, mpl::true_, mpl::false_ )
        {
            this->setGmc( geom );
            M_tensorConvection.update( geom );
            M_localMatrixInGeoContext.resize( this->gmc()->nPoints() );
            size_type eltId = this->gmc()->id();
            value_type hSize = M_expr.stabGLSParameter().hSize( eltId );
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                value_type norm_u_square = 0.;
                for (uint16_type c=0;c<tensor_convection_type::shape::M;++c)
                    norm_u_square += math::pow(M_tensorConvection.evalq( c,0,q),2);
                value_type norm_u = math::sqrt( norm_u_square );
                if ( math::abs( norm_u ) > 1e-10 )
                    M_localMatrixInGeoContext[q](0,0) = hSize/(2*norm_u);
                else
                    M_localMatrixInGeoContext[q](0,0) = 0.;
            }
        }
    void update( Geo_t const& geom, mpl::false_, mpl::false_ )
        {}

private :
    expr_type const& M_expr;
    tensor_convection_type M_tensorConvection;
    tensor_diffusion_type M_tensorDiffusion;
    std::vector<matrix_shape_type> M_localMatrixInGeoContext;
};

#endif


template<typename StabilizationGLSParameterType, typename ExprConvectionType, typename ExprDiffusionType, typename ExprReactionType, int QuadOrder>
class StabilizationGLSParameterExpr : public Feel::vf::ExprDynamicBase
{
public:

    typedef StabilizationGLSParameterExpr<StabilizationGLSParameterType,ExprConvectionType,ExprDiffusionType,ExprReactionType,QuadOrder> this_type;
    typedef ExprConvectionType expression_convection_type;
    typedef ExprDiffusionType expression_diffusion_type;
    typedef ExprReactionType expression_reaction_type;
    typedef StabilizationGLSParameterType stabilization_glsparameter_type;

    static const size_type context = vm::DYNAMIC;//ExprConvectionType::context|ExprDiffusionType::context|expression_reaction_type::context;
    static const bool is_terminal = true;
    typedef double value_type;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    StabilizationGLSParameterExpr( stabilization_glsparameter_type const& stabGLSParameter )
        :
        M_stabGLSParameter( stabGLSParameter )
        {}

    StabilizationGLSParameterExpr( stabilization_glsparameter_type const& stabGLSParameter,
                                   expression_convection_type const& exprConvection,
                                   expression_diffusion_type const& exprDiffusion,
                                   expression_reaction_type const& exprReaction,
                                   bool hasConvection, bool hasDiffusion, bool hasReaction )
        :
        M_stabGLSParameter( stabGLSParameter )
        {
            if ( hasConvection )
                this->setConvection( exprConvection );
            if ( hasDiffusion )
                this->setDiffusion( exprDiffusion );
            if ( hasReaction )
                this->setReaction( exprReaction );
        }

    StabilizationGLSParameterExpr( StabilizationGLSParameterExpr const& ) = default;
    StabilizationGLSParameterExpr( StabilizationGLSParameterExpr && ) = default;

    stabilization_glsparameter_type const& stabGLSParameter() const { return M_stabGLSParameter; }

    expression_convection_type const& exprConvection() const { CHECK( this->hasConvection() ) << "no convection expr"; return *M_exprConvection; }
    expression_diffusion_type const& exprDiffusion() const { CHECK( this->hasDiffusion() ) << "no convection expr"; return *M_exprDiffusion; }
    expression_reaction_type const& exprReaction() const { CHECK( this->hasReaction() ) << "no convection expr"; return *M_exprReaction; }

    void setConvection( expression_convection_type const& exprConvection ) { M_exprConvection.emplace( exprConvection ); }
    void setDiffusion( expression_diffusion_type const& exprDiffusion ) { M_exprDiffusion.emplace( exprDiffusion ); }
    void setReaction( expression_reaction_type const& exprReaction ) { M_exprReaction.emplace( exprReaction ); }

    bool hasConvection() const { return M_exprConvection.has_value(); }
    bool hasDiffusion() const { return M_exprDiffusion.has_value(); }
    bool hasReaction() const { return M_exprReaction.has_value(); }

    size_type dynamicContext() const
        {
            size_type res = 0;
            if ( this->hasConvection() )
                res = res | Feel::vf::dynamicContext( this->exprConvection() );
            if ( this->hasDiffusion() )
                res = res | Feel::vf::dynamicContext( this->exprDiffusion() );
            if ( this->hasReaction() )
                res = res | Feel::vf::dynamicContext( this->exprReaction() );
            return res;
        }

    //! polynomial order
    uint16_type polynomialOrder() const
    {
        if ( QuadOrder>=0 )
            return QuadOrder;
        uint16_type polyOrderConvection = this->hasConvection() ? this->exprConvection().polynomialOrder() : uint16_type(0);
        uint16_type polyOrderDiffusion = this->hasDiffusion()? this->exprDiffusion().polynomialOrder() : uint16_type(0);
        uint16_type polyOrderReaction = this->hasReaction()? this->exprReaction().polynomialOrder() : uint16_type(0);
        return std::max( polyOrderConvection, std::max( polyOrderDiffusion,polyOrderReaction ) );
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
        typedef Shape<gmc_type::NDim, Scalar, false,false> shape;
        struct is_zero { static const bool value = false; };

        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t,shape,value_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename tensorbase_type::ret_type;
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            self_type( expr,geom )
            {}
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            self_type( expr,geom )
            {}
        tensor( this_type const& expr, Geo_t const& geom )
            {
                if ( expr.hasReaction() )
                {
                    if ( expr.hasConvection() && expr.hasDiffusion() )
                        this->initTensorBase<true,true,true>( expr, geom );
                    else if ( expr.hasConvection() )
                        this->initTensorBase<true,false,true>( expr, geom );
                    else if ( expr.hasDiffusion() )
                        this->initTensorBase<false,true,true>( expr, geom );
                    else
                        this->initTensorBase<false,false,true>( expr, geom );
                }
                else
                {
                    if ( expr.hasConvection() && expr.hasDiffusion() )
                        this->initTensorBase<true,true,false>( expr, geom );
                    else if ( expr.hasConvection() )
                        this->initTensorBase<true,false,false>( expr, geom );
                    else if ( expr.hasDiffusion() )
                        this->initTensorBase<false,true,false>( expr, geom );
                    else
                        this->initTensorBase<false,false,false>( expr, geom );
                }
            }
        tensor( tensor const& t ) = default;

        template<bool HasConvec,bool HasDiff,bool HasReaction>
        void
        initTensorBase( this_type const& expr, Geo_t const& geom )
            {
                if ( expr.stabGLSParameter().method() == "eigenvalue" )
                    M_tensorbase.reset( new tensorStabGLSParameterEigenValue<Geo_t, Basis_i_t, Basis_j_t, shape, this_type,HasConvec,HasDiff,HasReaction,0>(expr,geom) );
                else if ( expr.stabGLSParameter().method() == "doubly-asymptotic-approximation" )
                    M_tensorbase.reset( new tensorStabGLSParameterEigenValue<Geo_t, Basis_i_t, Basis_j_t, shape, this_type,HasConvec,HasDiff,HasReaction,1>(expr,geom) );
                else
                    CHECK( false ) << "invalid method " << expr.stabGLSParameter().method();
            }

        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom )
        {
            M_tensorbase->update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensorbase->update( geom, face );
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return ret_type(M_tensorbase->evalq( q ).data());
        }

    private:
        tensorbase_ptrtype M_tensorbase;
    };

private :
    stabilization_glsparameter_type const& M_stabGLSParameter;
    std::optional<expression_convection_type> M_exprConvection;
    std::optional<expression_diffusion_type> M_exprDiffusion;
    std::optional<expression_reaction_type> M_exprReaction;
};

template<typename ExprConvectionType, typename ExprDiffusionType, typename ExprReactionType, int QuadOrder = -1, typename StabilizationGLSParameterType >
inline
    Expr< StabilizationGLSParameterExpr<StabilizationGLSParameterType,ExprConvectionType,ExprDiffusionType,ExprReactionType,QuadOrder > >
    stabilizationGLSParameterExpr( StabilizationGLSParameterType const& stabilizationGLSParameter )
{
    typedef StabilizationGLSParameterExpr<StabilizationGLSParameterType,ExprConvectionType,ExprDiffusionType,ExprReactionType,QuadOrder> expr_stabgls_parameter_t;
    return Expr< expr_stabgls_parameter_t >( expr_stabgls_parameter_t( stabilizationGLSParameter ) );
}

template<int QuadOrder = -1, typename StabilizationGLSParameterType, typename ExprConvectionType, typename ExprDiffusionType >
inline
    auto //Expr< StabilizationGLSParameterExpr<StabilizationGLSParameterType,ExprConvectionType,ExprDiffusionType,QuadOrder > >
stabilizationGLSParameterExpr( StabilizationGLSParameterType const& stabilizationGLSParameter,
                               ExprConvectionType const& exprConvection,
                               ExprDiffusionType const& exprDiffusion,
                               bool hasConvection = true,
                               bool hasDiffusion = true)
{
    bool hasReaction = false;
    auto exprReaction = cst(0.);
    using ExprReactionType = std::decay_t<decltype(exprReaction)>;
    typedef StabilizationGLSParameterExpr<StabilizationGLSParameterType,ExprConvectionType,ExprDiffusionType,ExprReactionType,QuadOrder> expr_stabgls_parameter_t;
    return Expr< expr_stabgls_parameter_t >(  expr_stabgls_parameter_t( stabilizationGLSParameter,exprConvection,exprDiffusion,exprReaction,
                                                                        hasConvection,hasDiffusion,hasReaction ) );
}

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_VF_STABILIZATION_GLS_PARMATER_EXPR_H
