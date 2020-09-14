#ifndef _CAUCHY_GREEN_INVARIANTS_EXPR_HPP
#define _CAUCHY_GREEN_INVARIANTS_EXPR_HPP 1

namespace Feel {
namespace FeelModels {

template<typename CauchyGreenTensorExprType>
class CauchyGreenInvariant1Expr
{
public:
    typedef CauchyGreenInvariant1Expr<CauchyGreenTensorExprType> this_type;
    typedef CauchyGreenTensorExprType expr_cauchygreentensor_type;
    typedef typename expr_cauchygreentensor_type::value_type value_type;
    typedef value_type evaluate_type;

    typedef Eigen::Matrix<value_type, 1, 1> matrix_shape_type;

    static const size_type context_cauchygreentensor = expr_cauchygreentensor_type::context;
    static const size_type context = expr_cauchygreentensor_type::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = expr_cauchygreentensor_type::template HasTestFunction<Func>::result;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = expr_cauchygreentensor_type::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = expr_cauchygreentensor_type::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = expr_cauchygreentensor_type::template has_trial_basis<Func>;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    //--------------------------------------------------------------------//
    // Constructors
    CauchyGreenInvariant1Expr( expr_cauchygreentensor_type const& ex )
        : M_cauchyGreenTensorExpr( ex )
    {}

    CauchyGreenInvariant1Expr( CauchyGreenInvariant1Expr const& c )
        : M_cauchyGreenTensorExpr( c.M_cauchyGreenTensorExpr )
    {}

    //! polynomial order
    uint16_type polynomialOrder() const { return M_cauchyGreenTensorExpr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_cauchyGreenTensorExpr.isPolynomial(); }

    //--------------------------------------------------------------------//
    // Accessors
    expr_cauchygreentensor_type const& cauchyGreenTensorExpr() const { return M_cauchyGreenTensorExpr; }

    //--------------------------------------------------------------------//
    // Expr tensor
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expr_cauchygreentensor_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        struct is_zero { static const bool value = tensor_expr_type::is_zero::value; };

        typedef typename tensor_expr_type::shape tensor_expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<tensor_expr_shape::M>,mpl::int_<tensor_expr_shape::N> >::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, ( mpl::int_<tensor_expr_shape::M>, mpl::int_<tensor_expr_shape::N> ) );
        typedef Shape<tensor_expr_shape::nDim, Tensor2, false, false> shape_cauchygreentensor;
        typedef Shape<tensor_expr_shape::nDim, Scalar, false,false> shape;

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_tensorExpr( expr.cauchyGreenTensorExpr(), geom, fev, feu )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensorExpr( expr.cauchyGreenTensorExpr(), geom, fev )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensorExpr( expr.cauchyGreenTensorExpr(), geom )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensorExpr.init( im );
        }

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_tensorExpr.update( geom );
            computeCauchyGreenInvariant1( mpl::int_<shape_cauchygreentensor::N>() );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensorExpr.update( geom, face );
            computeCauchyGreenInvariant1( mpl::int_<shape_cauchygreentensor::N>() );
        }

        matrix_shape_type const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return this->evalq( q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1,c2,q );
        }
        matrix_shape_type const&
        evaliq( uint16_type i, uint16_type q ) const
        {
            return this->evalq( q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_locRes[q](0,0);
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_locRes[q];
        }

    private:
        void computeCauchyGreenInvariant1( mpl::int_<1> )
        {
        }
        void computeCauchyGreenInvariant1( mpl::int_<2> )
        {
            /*
             * if A = [[a,b],[c,d]] then Cof(A) = [[d,-c],[-b,a]]
             * => Tr(Cof(A)) = a + d
             * => I1 = sqrt(Tr(Cof(A))) = sqrt(a+d)
             */
            for( uint16_type q = 0; q < M_locRes.size(); ++q )
            {
                double a = M_tensorExpr.evalq( 0, 0, q );
                double d = M_tensorExpr.evalq( 1, 1, q );
                M_locRes[q](0,0) = math::sqrt( a + d );
            }
        }
        void computeCauchyGreenInvariant1( mpl::int_<3> )
        {
            /*
             * if A = [[a,b,c],[d,e,f],[g,h,i]] then Tr(Cof(A)) = a*e + a*i - b*d - c*g + e*i - f*h
             */
            for( uint16_type q = 0; q < M_locRes.size(); ++q )
            {
                double a = M_tensorExpr.evalq( 0, 0, q );
                double b = M_tensorExpr.evalq( 0, 1, q );
                double c = M_tensorExpr.evalq( 0, 2, q );
                double d = M_tensorExpr.evalq( 1, 0, q );
                double e = M_tensorExpr.evalq( 1, 1, q );
                double f = M_tensorExpr.evalq( 1, 2, q );
                double g = M_tensorExpr.evalq( 2, 0, q );
                double h = M_tensorExpr.evalq( 2, 1, q );
                double i = M_tensorExpr.evalq( 2, 2, q );
                M_locRes[q](0,0) = math::sqrt( a*e+a*i-b*d-c*g+e*i-f*h );
            }
        }

    private:
        tensor_expr_type M_tensorExpr;
        std::vector<matrix_shape_type> M_locRes;
    };


private:
    expr_cauchygreentensor_type M_cauchyGreenTensorExpr;
};

template<typename CauchyGreenTensorExprType>
inline Expr< CauchyGreenInvariant1Expr<CauchyGreenTensorExprType> >
cauchyGreenInvariant1Expr( CauchyGreenTensorExprType const& tensorexpr )
{
    typedef CauchyGreenInvariant1Expr< CauchyGreenTensorExprType > expr_cauchygreeninvariant1_type;
    return Expr<expr_cauchygreeninvariant1_type>( expr_cauchygreeninvariant1_type( tensorexpr ) );
}

}
}

#endif
