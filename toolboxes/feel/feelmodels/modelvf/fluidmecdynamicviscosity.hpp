#ifndef FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H
#define FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{
enum class ExprOperatorType { ID=0,GRAD };

template<typename FieldType>
class ExprEvaluateFieldOperators
{
public :
    using this_type = ExprEvaluateFieldOperators<FieldType>;
    using field_type = FieldType;
    using field_clean_type = unwrap_ptr_t<field_type>;
    using expr_grad_type = std::decay_t<decltype(gradv(field_type{}))>;
    using expr_laplacian_type = std::decay_t<decltype(laplacianv(field_type{}))>;
    using value_type = typename field_clean_type::value_type;

    static const bool laplacian_is_zero = field_clean_type::fe_type::nOrder < 2;

    ExprEvaluateFieldOperators( field_type const& field )
        :
        M_field( field  ),
        M_enableGrad( false ), M_enableLaplacian( false )
        {}

    ExprEvaluateFieldOperators( ExprEvaluateFieldOperators const& ) = default;
    ExprEvaluateFieldOperators( ExprEvaluateFieldOperators && ) = default;

    bool isEnabledGrad() const { return M_enableGrad; }
    void setEnableGrad( bool b )
        {
            M_enableGrad = b;
            if ( b && !M_exprGrad )
                M_exprGrad.emplace( gradv( M_field ) );
        }
    bool isEnabledLaplacian() const { return M_enableLaplacian; }
    void setEnableLaplacian( bool b )
        {
            M_enableLaplacian = b;
            if ( b && !M_exprLaplacian )
                M_exprLaplacian.emplace( laplacianv( M_field ) );
        }

    expr_grad_type const& exprGrad() const { return *M_exprGrad; }
    expr_laplacian_type const& exprLaplacian() const { return *M_exprLaplacian; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using tensor_expr_grad_type = typename expr_grad_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_laplacian_type = typename expr_laplacian_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        using shape_grad_type = typename tensor_expr_grad_type::shape;
        using tensor_base_grad_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t, shape_grad_type, typename this_type::value_type>;
        using matrix_shape_grad_type = typename tensor_base_grad_type::matrix_shape_type;
        using array_shape_grad_type = typename tensor_base_grad_type::new_array_shape_type;

        using shape_laplacian_type = typename tensor_expr_laplacian_type::shape;
        using tensor_base_laplacian_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t, shape_laplacian_type, typename this_type::value_type>;
        using matrix_shape_laplacian_type = typename tensor_base_laplacian_type::matrix_shape_type;
        using array_shape_laplacian_type = typename tensor_base_laplacian_type::new_array_shape_type;

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr )
            {
                this->initTensor( expr, geom, fev, feu );
            }
        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr )
            {
                this->initTensor( expr, geom, fev );
            }
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr )
            {
                this->initTensor( expr, geom );
            }

        array_shape_grad_type const& localEvalGrad() const { return M_localEvalGrad; }
        matrix_shape_grad_type const& localEvalGrad(  uint16_type q ) const { return M_localEvalGrad[q]; }

        array_shape_laplacian_type const& localEvalLaplacian() const { return M_localEvalLaplacian; }
        matrix_shape_laplacian_type const& localEvalLaplacian(  uint16_type q ) const { return M_localEvalLaplacian[q]; }

        void update( Geo_t const& geom )
            {
                bool hasGrad = false, hasLaplacian = false;
                if ( M_tensorExprGrad )
                {
                    M_tensorExprGrad->update( geom );
                    hasGrad = true;
                }
                if constexpr ( !laplacian_is_zero )
                {
                    if ( M_tensorExprLaplacian )
                    {
                        M_tensorExprLaplacian->update( geom );
                        hasLaplacian = true;
                    }
                }
                this->updateImpl( geom, hasGrad, hasLaplacian );
            }
        void update( Geo_t const& geom, uint16_type face )
            {
                bool hasGrad = false, hasLaplacian = false;
                if ( M_tensorExprGrad )
                {
                    M_tensorExprGrad->update( geom, face );
                    hasGrad = true;
                }
                if constexpr ( !laplacian_is_zero )
                {
                    if ( M_tensorExprLaplacian )
                    {
                        M_tensorExprLaplacian->update( geom, face );
                        hasLaplacian = true;
                    }
                }
                this->updateImpl( geom, hasGrad, hasLaplacian );
            }
    private :
        template<typename... TheArgsType>
        void initTensor( this_type const& expr, const TheArgsType&... theInitArgs )
            {
                if ( expr.isEnabledGrad() )
                    M_tensorExprGrad.emplace( expr.exprGrad(), theInitArgs... );
                if constexpr ( !laplacian_is_zero )
                {
                    if ( expr.isEnabledLaplacian() )
                        M_tensorExprLaplacian.emplace( expr.exprLaplacian(), theInitArgs... );
                }
            }
        void updateImpl(  Geo_t const& geom, bool hasGrad, bool hasLaplacian )
            {
                if ( !hasGrad && !hasLaplacian )
                    return;

                auto gmc = fusion::at_key<typename tensor_base_grad_type::key_type>( geom );

                uint16_type nPoints = gmc->nPoints();
                if ( hasGrad && M_localEvalGrad.size() != nPoints )
                {
                    M_localEvalGrad.resize( nPoints );
                    //M_localEvalGrad.setConstant( nPoints, this->M_zeroLocTensor2 );
                }
                if constexpr ( !laplacian_is_zero )
                {
                    if ( hasLaplacian && M_localEvalLaplacian.size() != nPoints )
                        M_localEvalLaplacian.resize( nPoints );
                }

                for ( uint16_type q=0;q< nPoints;++q )
                {
                    if ( hasGrad )
                    {
                        if constexpr ( false /*expr_grad_type::is_terminal*/ )
                                     {
                                         M_localEvalGrad[q] = M_tensorExprGrad->evalq( q ); //not compile, need to investigate
                                     }
                        else
                        {
                            matrix_shape_grad_type& locData = M_localEvalGrad[q];
                            for (uint16_type c1=0;c1<shape_grad_type::M;++c1 )
                                for (uint16_type c2=0;c2<shape_grad_type::N;++c2 )
                                    locData(c1,c2) = M_tensorExprGrad->evalq( c1,c2,q );
                        }
                    }

                    if constexpr ( !laplacian_is_zero )
                    {
                        if ( hasLaplacian )
                        {
                            matrix_shape_laplacian_type& locData = M_localEvalLaplacian[q];
                            for (uint16_type c1=0;c1<shape_laplacian_type::M;++c1 )
                                for (uint16_type c2=0;c2<shape_laplacian_type::N;++c2 )
                                    locData(c1,c2) = M_tensorExprLaplacian->evalq( c1,c2,q );
                        }
                    }
                }
            }



    private :
        this_type const& M_expr;
        std::optional<tensor_expr_grad_type> M_tensorExprGrad;
        std::optional<tensor_expr_laplacian_type> M_tensorExprLaplacian;
        array_shape_grad_type M_localEvalGrad;
        array_shape_laplacian_type M_localEvalLaplacian;
    };
private :
    field_type const& M_field;
    bool M_enableGrad, M_enableLaplacian;
    std::optional<expr_grad_type> M_exprGrad;
    std::optional<expr_laplacian_type> M_exprLaplacian;
};

// only grad currently
template<typename ExprGradType>
class ExprEvaluateFieldOperatorGradFromExpr
{
public :
    using this_type = ExprEvaluateFieldOperatorGradFromExpr<ExprGradType>;
    using expr_grad_type = ExprGradType;
    using value_type = typename ExprGradType::value_type;
    ExprEvaluateFieldOperatorGradFromExpr( expr_grad_type const& exprGrad )
        :
        M_enableGrad( false ),
        M_exprGrad( exprGrad )
        {}

    ExprEvaluateFieldOperatorGradFromExpr( ExprEvaluateFieldOperatorGradFromExpr const& ) = default;
    ExprEvaluateFieldOperatorGradFromExpr( ExprEvaluateFieldOperatorGradFromExpr && ) = default;

    bool isEnabledGrad() const { return M_enableGrad; }
    void setEnableGrad( bool b ) { M_enableGrad = b; }

    expr_grad_type const& exprGrad() const { return M_exprGrad; }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using tensor_expr_grad_type = typename expr_grad_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        using shape_grad_type = Shape<tensor_expr_grad_type::shape::nDim,Tensor2, false, false>;
        using tensor_base_grad_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                                 shape_grad_type, typename this_type::value_type>;

        using matrix_shape_grad_type = typename tensor_base_grad_type::matrix_shape_type;
        using array_shape_grad_type = typename tensor_base_grad_type::new_array_shape_type;


        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr ),
            M_tensorExprGrad( expr.exprGrad(), geom, fev, feu )
            {}
        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr ),
            M_tensorExprGrad( expr.exprGrad(), geom, fev )
            {}
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr ),
            M_tensorExprGrad( expr.exprGrad(), geom )
            {}

        array_shape_grad_type const& localEvalGrad() const { return M_localEvalGrad; }
        matrix_shape_grad_type const& localEvalGrad(  uint16_type q ) const { return M_localEvalGrad[q]; }

        void update( Geo_t const& geom )
            {
                if ( !M_expr.isEnabledGrad() )
                    return;
                M_tensorExprGrad.update( geom );
                this->updateImpl( geom );
            }
        void update( Geo_t const& geom, uint16_type face )
            {
                if ( !M_expr.isEnabledGrad() )
                    return;
                M_tensorExprGrad.update( geom, face );
                this->updateImpl( geom );
            }
    private :
        void updateImpl(  Geo_t const& geom )
            {
                auto gmc = fusion::at_key<typename tensor_base_grad_type::key_type>( geom );

                uint16_type nPoints = gmc->nPoints();
                if ( M_localEvalGrad.size() != nPoints )
                {
                    M_localEvalGrad.resize( nPoints );
                    //M_localEvalGrad.setConstant( nPoints, this->M_zeroLocTensor2 );
                }

                for ( uint16_type q=0;q< nPoints;++q )
                {
                    if constexpr ( false /*expr_grad_type::is_terminal*/ )
                                 {
                                     M_localEvalGrad[q] = M_tensorExprGrad.evalq( q ); //not compile, need to investigate
                                 }
                    else
                    {
                        matrix_shape_grad_type& locData = M_localEvalGrad[q];
                        for (uint16_type c1=0;c1<shape_grad_type::M;++c1 )
                            for (uint16_type c2=0;c2<shape_grad_type::N;++c2 )
                                locData(c1,c2) = M_tensorExprGrad.evalq( c1,c2,q );
                    }
                }
            }



    private :
        this_type const& M_expr;
        tensor_expr_grad_type M_tensorExprGrad;
        array_shape_grad_type M_localEvalGrad;
    };
private :
    bool M_enableGrad;
    expr_grad_type M_exprGrad;
};


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

    FluidMecDynamicViscosityBase( expr_type const& expr )
        :
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
    FluidMecDynamicViscosityNewtonian( ExprType const& expr )
        :
        super_type( expr ),
        M_exprDynamicVisocsity( expr.template materialPropertyExpr<1,1>("dynamic-viscosity") )
        {}

    template <ExprOperatorType TheExprOp = expr_type::exprOp, std::enable_if_t< TheExprOp == ExprOperatorType::GRAD, bool> = true>
    FluidMecDynamicViscosityNewtonian( ExprType const& expr )
        :
        super_type( expr ),
        M_exprDynamicVisocsity( grad<expr_type::nDim>( expr.template materialPropertyExpr<1,1>("dynamic-viscosity")/*, "", world, dirLibExpr*/ /*TODO*/ ) )
        {}

    FluidMecDynamicViscosityNewtonian( FluidMecDynamicViscosityNewtonian const& ) = default;
    FluidMecDynamicViscosityNewtonian( FluidMecDynamicViscosityNewtonian && ) = default;

    bool dependsOnVelocityField() const override { return false; }

    size_type dynamicContext() const override { return Feel::vf::dynamicContext( M_exprDynamicVisocsity ); }

    uint16_type polynomialOrder() const override { return M_exprDynamicVisocsity.polynomialOrder(); }

    expr_dynamic_viscosity_newtonian_type const& exprDynamicVisocsity() const { return M_exprDynamicVisocsity; }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_exprDynamicVisocsity.hasSymbolDependency( symb, se );
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

        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicVisocsity().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicVisocsity().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_muExprTensor( this->expr().exprDynamicVisocsity().evaluator( geom ) )
            {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            super_type( geom ),
            M_expr( expr ),
            M_muExprTensor( std::true_type{}, exprExpanded.exprDynamicVisocsity(), ttse, expr.exprDynamicVisocsity(), geom, theInitArgs... )
            {}


        this_type const& expr() const { return M_expr; }

        void update( Geo_t const& geom ) override
            {
                M_muExprTensor.update( geom );
                this->updateImpl();
            }

        void update( Geo_t const& geom, uint16_type face ) override
            {
                M_muExprTensor.update( geom, face );
                this->updateImpl();
            }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_muExprTensor.update( std::true_type{}, exprExpanded.exprDynamicVisocsity(), ttse, geom, theUpdateArgs... );
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
    expr_dynamic_viscosity_newtonian_type M_exprDynamicVisocsity;
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

    FluidMecDynamicViscosityPowerLaw( ExprType const& expr )
        :
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


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( this->expr().kExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_muMinExprTensor( this->expr().muMinExpr().evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().muMaxExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_kExprTensor( this->expr().kExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_muMinExprTensor( this->expr().muMinExpr().evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().muMaxExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom )
            :
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
                Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
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

        void update( Geo_t const& geom, uint16_type face ) override
            {
                this->setGmc( geom );
                M_kExprTensor.update( geom, face );
                M_nExprTensor.update( geom, face );
                M_muMinExprTensor.update( geom, face );
                M_muMaxExprTensor.update( geom, face );
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

    FluidMecDynamicViscosityCarreau( ExprType const& expr )
        :
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


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom )
            :
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
                Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
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

        void update( Geo_t const& geom, uint16_type face ) override
            {
                this->setGmc( geom );
                M_mu0ExprTensor.update( geom, face );
                M_muInfExprTensor.update( geom, face );
                M_lambdaExprTensor.update( geom, face );
                M_nExprTensor.update( geom, face );
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

    FluidMecDynamicViscosityCarreauYasuda( ExprType const& expr )
        :
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


        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_aExprTensor( this->expr().aExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprEvaluateVelocityOperators( tensorExprMain.tensorExprEvaluateVelocityOperatorsPtr() ),
            M_mu0ExprTensor( this->expr().mu0Expr().evaluator( geom ) ),
            M_muInfExprTensor( this->expr().muInfExpr().evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().lambdaExpr().evaluator( geom ) ),
            M_nExprTensor( this->expr().nExpr().evaluator( geom ) ),
            M_aExprTensor( this->expr().aExpr().evaluator( geom ) )
            {}
        tensor( this_type const& expr, typename this_type::super_type::template tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, Geo_t const& geom )
            :
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
                Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
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

        void update( Geo_t const& geom, uint16_type face ) override
            {
                this->setGmc( geom );
                M_mu0ExprTensor.update( geom, face );
                M_muInfExprTensor.update( geom, face );
                M_lambdaExprTensor.update( geom, face );
                M_nExprTensor.update( geom, face );
                M_aExprTensor.update( geom, face );
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


template< typename ExprType>
template <typename TheSymbolExprType>
bool
FluidMecDynamicViscosityBase<ExprType>::hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
{
    auto const& dynamicViscosity = this->expr().dynamicViscosity();
    if ( dynamicViscosity.isNewtonianLaw() )
        return static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosity.isPowerLaw() )
        return static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosity.isCarreauLaw() )
        return static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else if ( dynamicViscosity.isCarreauYasudaLaw() )
        return static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this).hasSymbolDependency( symb, se );
    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
    return false;
}



template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename... TheArgsType>
std::shared_ptr<typename FluidMecDynamicViscosityBase<ExprType>::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
FluidMecDynamicViscosityBase<ExprType>::evaluator( tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain, const TheArgsType&... theInitArgs ) const
{
    auto const& dynamicViscosity = this->expr().dynamicViscosity();
    if ( dynamicViscosity.isNewtonianLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isPowerLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isCarreauLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isCarreauYasudaLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this), tensorExprMain, theInitArgs... );
    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
    return std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >{};
}

template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType,typename... TheArgsType>
std::shared_ptr<typename FluidMecDynamicViscosityBase<ExprType>::template tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >
FluidMecDynamicViscosityBase<ExprType>::evaluator( std::true_type/**/, tensor_main_type<Geo_t,Basis_i_t,Basis_j_t> const& tensorExprMain,
                                                   TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, const TheArgsType&... theInitArgs ) const
{
    using expr_expanded_type = typename TheExprExpandedType::expr_type;

    auto const& dynamicViscosity = this->expr().dynamicViscosity();
    if ( dynamicViscosity.isNewtonianLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{},
                                                                                                                                     static_cast<FluidMecDynamicViscosityNewtonian<expr_expanded_type> const&>(exprExpanded) ,
                                                                                                                                     ttse,
                                                                                                                                     static_cast<FluidMecDynamicViscosityNewtonian<expr_type> const&>(*this),
                                                                                                                                     tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isPowerLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{},
                                                                                                                                    static_cast<FluidMecDynamicViscosityPowerLaw<expr_expanded_type> const&>(exprExpanded) ,
                                                                                                                                    ttse,
                                                                                                                                    static_cast<FluidMecDynamicViscosityPowerLaw<expr_type> const&>(*this),
                                                                                                                                    tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isCarreauLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{},
                                                                                                                                   static_cast<FluidMecDynamicViscosityCarreau<expr_expanded_type> const&>(exprExpanded) ,
                                                                                                                                   ttse,
                                                                                                                                   static_cast<FluidMecDynamicViscosityCarreau<expr_type> const&>(*this),
                                                                                                                                   tensorExprMain, theInitArgs... );
    else if ( dynamicViscosity.isCarreauYasudaLaw() )
        return std::make_shared< typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( std::true_type{},
                                                                                                                                         static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_expanded_type> const&>(exprExpanded) ,
                                                                                                                                         ttse,
                                                                                                                                         static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_type> const&>(*this),
                                                                                                                                         tensorExprMain, theInitArgs... );
    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
    return std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> >{};
}


template< typename ExprType>
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
void
FluidMecDynamicViscosityBase<ExprType>::updateEvaluator( std::true_type /**/, std::shared_ptr<tensor_base_type<Geo_t,Basis_i_t,Basis_j_t> > & tensorToUpdate, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                                                         Geo_t const& geom, const TheArgsType&... theUpdateArgs )
{
    using expr_expanded_type = typename TheExprExpandedType::expr_type;
    auto const& dynamicViscosity = this->expr().dynamicViscosity();
    if ( dynamicViscosity.isNewtonianLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityNewtonian<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{},
                                                                                                                                                               static_cast<FluidMecDynamicViscosityNewtonian<expr_expanded_type> const&>(exprExpanded),
                                                                                                                                                               ttse, geom, theUpdateArgs... );
    else if ( dynamicViscosity.isPowerLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityPowerLaw<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{},
                                                                                                                                                               static_cast<FluidMecDynamicViscosityPowerLaw<expr_expanded_type> const&>(exprExpanded),
                                                                                                                                                               ttse, geom, theUpdateArgs... );
    else if ( dynamicViscosity.isCarreauLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityCarreau<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{},
                                                                                                                                                             static_cast<FluidMecDynamicViscosityCarreau<expr_expanded_type> const&>(exprExpanded),
                                                                                                                                                             ttse, geom, theUpdateArgs... );

    else if ( dynamicViscosity.isCarreauYasudaLaw() )
        std::static_pointer_cast<typename FluidMecDynamicViscosityCarreauYasuda<expr_type>::template tensor<Geo_t,Basis_i_t,Basis_j_t>>( tensorToUpdate )->update( std::true_type{},
                                                                                                                                                                   static_cast<FluidMecDynamicViscosityCarreauYasuda<expr_expanded_type> const&>(exprExpanded),
                                                                                                                                                                   ttse, geom, theUpdateArgs... );
    else
        CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
}


template<typename ExprEvaluateFieldOperatorsType, typename FiniteElementVelocityType, typename ModelPhysicFluidType, typename SymbolsExprType, ExprApplyType ExprApplied, ExprOperatorType ExprOp = ExprOperatorType::ID>
class FluidMecDynamicViscosityImpl// : public Feel::vf::ExprDynamicBase
{
public:

    typedef FluidMecDynamicViscosityImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ModelPhysicFluidType,SymbolsExprType,ExprApplied,ExprOp> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_velocity|vm::DYNAMIC;

    using model_physic_fluid_type = ModelPhysicFluidType;
    using symbols_expr_type = SymbolsExprType;

    static constexpr bool is_applied_as_eval = (ExprApplied == ExprApplyType::EVAL);
    static constexpr bool is_applied_as_jacobian = (ExprApplied == ExprApplyType::JACOBIAN);

    static constexpr ExprOperatorType exprOp = ExprOp;

    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorsType;
    using value_type = typename expr_evaluate_velocity_opertors_type::value_type;

    static constexpr uint16_type nDim = model_physic_fluid_type::nDim;

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
                                  model_physic_fluid_type const& physicFluid,
                                  MaterialProperties const& matProps,
                                  uint16_type polyOrder,
                                  SymbolsExprType const& se )
        :
        M_exprEvaluateVelocityOperators( exprEvaluateVelocityOperators ),
        M_physicFluid( physicFluid ),
        M_matProps( matProps ),
        M_polynomialOrder( polyOrder ),
        M_se( se )
        {
            this->initExprDynamicViscosityBase();
        }

    FluidMecDynamicViscosityImpl( FluidMecDynamicViscosityImpl const & op )
        :
        M_exprEvaluateVelocityOperators( op.M_exprEvaluateVelocityOperators ),
        M_physicFluid( op.M_physicFluid ),
        M_matProps( op.M_matProps ),
        M_polynomialOrder( op.M_polynomialOrder ),
        M_se( op.M_se )
        {
            this->initExprDynamicViscosityBase();
        }
    FluidMecDynamicViscosityImpl( FluidMecDynamicViscosityImpl && op )
        :
        M_exprEvaluateVelocityOperators( std::move( op.M_exprEvaluateVelocityOperators ) ),
        M_physicFluid( std::move( op.M_physicFluid ) ),
        M_matProps( std::move( op.M_matProps ) ),
        M_polynomialOrder( std::move( op.M_polynomialOrder ) ),
        M_se( std::move( op.M_se ) )
        {
            this->initExprDynamicViscosityBase();
        }

    ~FluidMecDynamicViscosityImpl()
    {}


    // allow to use the current SymbolsExpr in a concat operation
    // TODO : improve management SymbolsExpr with tensor ctx with concat
    auto const& symbolsExpressionWithoutTensorContext() const
        {
            if constexpr ( is_symbols_expression_v<SymbolsExprType> )
                return M_se;
            else
                return M_se.symbolsExpression();
        }

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
            uint16_type orderGradVelocity = M_exprEvaluateVelocityOperators->exprGrad().polynomialOrder();

            uint16_type res = 2*(orderGradVelocity+1); // default value for non newtonian
            if ( this->dynamicViscosity().isNewtonianLaw() )
            {
                res = M_exprDynamicViscosityBase->polynomialOrder();
            }
            return res;
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
            using new_this_type = FluidMecDynamicViscosityImpl<ExprEvaluateFieldOperatorsType,FiniteElementVelocityType,ModelPhysicFluidType,new_se_type,ExprApplied,ExprOp>;
            return new_this_type( M_exprEvaluateVelocityOperators,M_physicFluid,M_matProps,M_polynomialOrder,newse );
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

    auto const& dynamicViscosity() const { return M_physicFluid.dynamicViscosity(); }
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

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr )
        {
            this->initTensor( expr, true, geom, fev, feu );
        }
        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr )
        {
            this->initTensor( expr, true, geom, fev );
        }
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_expr( expr )
        {
            this->initTensor( expr, true, geom );
        }
        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators,  Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom, fev, feu );
        }
        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom, fev );
        }
        tensor( this_type const& expr, std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperators, Geo_t const& geom )
            :
            M_expr( expr )
        {
            CHECK( tensorExprEvaluateVelocityOperators ) << "tensorExprEvaluateVelocityOperators not init";
            M_tensorExprEvaluateVelocityOperators = tensorExprEvaluateVelocityOperators;
            this->initTensor( expr, false, geom );
        }

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_expr( expr )
            {
                this->initTensor( std::true_type{}, true, exprExpanded, ttse, expr, geom, theInitArgs... );
            }


        std::shared_ptr<tensor_expr_evaluate_velocity_opertors_type> tensorExprEvaluateVelocityOperatorsPtr() const { return M_tensorExprEvaluateVelocityOperators; }
        tensor_expr_evaluate_velocity_opertors_type const& tensorExprEvaluateVelocityOperators() const { return *M_tensorExprEvaluateVelocityOperators; }

        template<typename IM>
        void init( IM const& im ) {}

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
        void update( Geo_t const& geom, uint16_type face, bool upEvaluateVelocityOperators = true  )
        {
            if ( upEvaluateVelocityOperators )
                M_tensorExprEvaluateVelocityOperators->update( geom, face );
            M_tensorbase->update( geom, face );
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
            auto const& dynamicViscosity = this->dynamicViscosity();
            if ( dynamicViscosity.isNewtonianLaw() )
                M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityNewtonian<this_type>( *this ) );
            else if ( dynamicViscosity.isPowerLaw() )
                M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityPowerLaw<this_type>( *this ) );
            else if ( dynamicViscosity.isCarreauLaw() )
                M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityCarreau<this_type>( *this ) );
            else if ( dynamicViscosity.isCarreauYasudaLaw() )
                M_exprDynamicViscosityBase.reset( new FluidMecDynamicViscosityCarreauYasuda<this_type>( *this ) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
        }

private:
    std::shared_ptr<expr_dynamic_viscosity_base_type> M_exprDynamicViscosityBase;
    std::shared_ptr<expr_evaluate_velocity_opertors_type> M_exprEvaluateVelocityOperators;
    model_physic_fluid_type const& M_physicFluid;
    MaterialProperties const& M_matProps;
    uint16_type M_polynomialOrder;
    SymbolsExprType M_se;
};
/// \endcond

template<class ExprGradVelocityType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecViscosity( Expr<ExprGradVelocityType> const& grad_u,
                   ModelPhysicFluidType const& physicFluid,
                   MaterialProperties const& matProps,
                   SymbolsExprType const& se = symbols_expression_empty_t{},
                   uint16_type polyOrder = invalid_uint16_type_value )
{
    using expr_evaluate_velocity_opertors_type = ExprEvaluateFieldOperatorGradFromExpr<Expr<ExprGradVelocityType>>;
    typedef FluidMecDynamicViscosityImpl<expr_evaluate_velocity_opertors_type,std::nullptr_t,ModelPhysicFluidType,SymbolsExprType,ExprApplyType::EVAL> fmstresstensor_t;
    auto exprEvaluateVelocityOperators = std::make_shared<expr_evaluate_velocity_opertors_type>( grad_u );
    return Expr< fmstresstensor_t >(  fmstresstensor_t( exprEvaluateVelocityOperators,physicFluid,matProps,polyOrder,se ) );
}

} // namespace FeelModels
} // namespace Feel

#endif // FEELPP_MODELS_VF_FLUIDMEC_DYNAMIC_VISCOSITY_H
