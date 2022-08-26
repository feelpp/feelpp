#ifndef FEELPP_MODELS_VF_EXPREVALUATEFIELDOPERATORS_H
#define FEELPP_MODELS_VF_EXPREVALUATEFIELDOPERATORS_H

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
    static constexpr uint16_type nDim = field_clean_type::nDim;

    using expr_id_type = std::decay_t<decltype(idv(field_type{}))>;
    using expr_grad_type = std::decay_t<decltype(gradv(field_type{}))>;
    using expr_laplacian_type = std::decay_t<decltype(laplacianv(field_type{}))>;
    using value_type = typename field_clean_type::value_type;

    static const bool laplacian_is_zero = field_clean_type::fe_type::nOrder < 2;

    ExprEvaluateFieldOperators( field_type const& field )
        :
        M_field( field  ),
        M_enableId( false ),M_enableGrad( false ), M_enableLaplacian( false )
        {}

    ExprEvaluateFieldOperators( ExprEvaluateFieldOperators const& ) = default;
    ExprEvaluateFieldOperators( ExprEvaluateFieldOperators && ) = default;

    bool isEnabledId() const { return M_enableId; }
    void setEnableId( bool b )
        {
            M_enableId = b;
            if ( b && !M_exprId )
                M_exprId.emplace( idv( M_field ) );
        }

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

    expr_id_type const& exprId() const { return *M_exprId; }
    expr_grad_type const& exprGrad() const { return *M_exprGrad; }
    expr_laplacian_type const& exprLaplacian() const { return *M_exprLaplacian; }

    size_type dynamicContext() const
        {
            size_type res = 0;
            if ( M_exprId )
                res = res | Feel::vf::dynamicContext( *M_exprId );
            if ( M_exprGrad )
                res = res | Feel::vf::dynamicContext( *M_exprGrad );
            if ( M_exprLaplacian )
                res = res | Feel::vf::dynamicContext( *M_exprLaplacian );
            return res;
        }

    static uint16_type polynomialOrderId() { return expr_id_type::expression_type::element_type::functionspace_type::basis_type::nOrder; }
    static uint16_type polynomialOrderGrad() { return std::max( expr_id_type::expression_type::element_type::functionspace_type::basis_type::nOrder-1,0); }
    static uint16_type polynomialOrderLaplacian() { return std::max( expr_id_type::expression_type::element_type::functionspace_type::basis_type::nOrder-2,0); }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        using tensor_expr_id_type = typename expr_id_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_grad_type = typename expr_grad_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using tensor_expr_laplacian_type = typename expr_laplacian_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;

        using shape_id_type = typename tensor_expr_id_type::shape;
        using tensor_base_id_type = tensorBase<Geo_t,Basis_i_t,Basis_j_t, shape_id_type, typename this_type::value_type>;
        using matrix_shape_id_type = typename tensor_base_id_type::matrix_shape_type;
        using array_shape_id_type = typename tensor_base_id_type::new_array_shape_type;


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

        array_shape_id_type const& localEvalId() const { return M_localEvalId; }
        matrix_shape_id_type const& localEvalId( uint16_type q ) const { return M_localEvalId[q]; }

        array_shape_grad_type const& localEvalGrad() const { return M_localEvalGrad; }
        matrix_shape_grad_type const& localEvalGrad( uint16_type q ) const { return M_localEvalGrad[q]; }

        array_shape_laplacian_type const& localEvalLaplacian() const { return M_localEvalLaplacian; }
        matrix_shape_laplacian_type const& localEvalLaplacian( uint16_type q ) const { return M_localEvalLaplacian[q]; }

        void update( Geo_t const& geom )
            {
                bool hasId = false, hasGrad = false, hasLaplacian = false;
                if ( M_tensorExprId )
                {
                    M_tensorExprId->update( geom );
                    hasId = true;
                }
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
                this->updateImpl( geom, hasId, hasGrad, hasLaplacian );
            }
    private :
        template<typename... TheArgsType>
        void initTensor( this_type const& expr, const TheArgsType&... theInitArgs )
            {
                if ( expr.isEnabledId() )
                    M_tensorExprId.emplace( expr.exprId(), theInitArgs... );
                if ( expr.isEnabledGrad() )
                    M_tensorExprGrad.emplace( expr.exprGrad(), theInitArgs... );
                if constexpr ( !laplacian_is_zero )
                {
                    if ( expr.isEnabledLaplacian() )
                        M_tensorExprLaplacian.emplace( expr.exprLaplacian(), theInitArgs... );
                }
            }
        void updateImpl(  Geo_t const& geom, bool hasId, bool hasGrad, bool hasLaplacian )
            {
                if ( !hasId && !hasGrad && !hasLaplacian )
                    return;

                auto gmc = fusion::at_key<typename tensor_base_grad_type::key_type>( geom );

                uint16_type nPoints = gmc->nPoints();
                if ( hasId && M_localEvalId.size() != nPoints )
                {
                    M_localEvalId.resize( nPoints );
                }
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
                    if ( hasId )
                    {
                            matrix_shape_id_type& locData = M_localEvalId[q];
                            for (uint16_type c1=0;c1<shape_id_type::M;++c1 )
                                for (uint16_type c2=0;c2<shape_id_type::N;++c2 )
                                    locData(c1,c2) = M_tensorExprId->evalq( c1,c2,q );
                    }
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
        std::optional<tensor_expr_id_type> M_tensorExprId;
        std::optional<tensor_expr_grad_type> M_tensorExprGrad;
        std::optional<tensor_expr_laplacian_type> M_tensorExprLaplacian;
        array_shape_id_type M_localEvalId;
        array_shape_grad_type M_localEvalGrad;
        array_shape_laplacian_type M_localEvalLaplacian;
    };
private :
    field_type const& M_field;
    bool M_enableId, M_enableGrad, M_enableLaplacian;
    std::optional<expr_id_type> M_exprId;
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
    static constexpr uint16_type nDim = expr_grad_type::functionspace_type::nDim;

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


} // namespace FeelModels
} // namespace Feel
#endif
