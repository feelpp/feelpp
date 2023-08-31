#ifndef _LEVELSET_HEAVISIDE_EXPR_HPP
#define _LEVELSET_HEAVISIDE_EXPR_HPP 1

namespace Feel {
namespace FeelModels {

template<typename LSExprType, typename EpsExprType, int IMOrder>
class LevelsetHeavisideExpr
{
public:
    typedef LevelsetHeavisideExpr<LSExprType, EpsExprType, IMOrder> this_type;
    typedef LSExprType expr_levelsetphi_type;
    typedef EpsExprType expr_epsilon_type;
    typedef typename expr_levelsetphi_type::value_type value_type;
    typedef value_type evaluate_type;

    typedef Eigen::Matrix<value_type, 1, 1> matrix_shape_type;

    static const size_type context_levelsetphi = expr_levelsetphi_type::context;
    static const size_type context = expr_levelsetphi_type::context;
    static const bool is_terminal = false;

    //static inline const uint16_type imorder = expr_levelsetphi_type::imorder;
    //static inline const uint16_type imorderDefault = 2*expr_levelsetphi_type::imorder;
    //static inline const uint16_type imorder = (IMOrder>=0) ? IMOrder: imorderDefault;
    //static const bool imIsPoly = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = expr_levelsetphi_type::template HasTestFunction<Func>::result;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = expr_levelsetphi_type::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = expr_levelsetphi_type::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = expr_levelsetphi_type::template has_trial_basis<Func>;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    //--------------------------------------------------------------------//
    // Constructors
    LevelsetHeavisideExpr( expr_levelsetphi_type const& ex, expr_epsilon_type const& eps )
        : M_levelsetPhiExpr( ex )
        , M_thicknessHeaviside( eps )
    {}

    LevelsetHeavisideExpr( LevelsetHeavisideExpr const& c )
        : M_levelsetPhiExpr( c.M_levelsetPhiExpr )
        , M_thicknessHeaviside( c.M_thicknessHeaviside )
    {}

    //--------------------------------------------------------------------//
    //! polynomial order
    uint16_type polynomialOrder() const { return (IMOrder>=0) ? IMOrder : 2*M_levelsetPhiExpr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return false; }

    //--------------------------------------------------------------------//
    // Accessors
    expr_levelsetphi_type const& levelsetPhiExpr() const { return M_levelsetPhiExpr; }
    expr_epsilon_type const& thickness() const { return M_thicknessHeaviside; }

    //--------------------------------------------------------------------//
    // Expr tensor
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expr_levelsetphi_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename expr_epsilon_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_epsilon_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape tensor_expr_shape;
        BOOST_MPL_ASSERT_MSG( ( tensor_expr_shape::M == 1 && tensor_expr_shape::N == 1 ), 
                INVALID_TENSOR_SHOULD_BE_RANK_0, ( mpl::int_<tensor_expr_shape::M>, mpl::int_<tensor_expr_shape::N> ) );
        typedef Shape<tensor_expr_shape::nDim, Scalar, false,false> shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_phiExpr( expr.levelsetPhiExpr(), geom, fev, feu )
            , M_epsilon( expr.thickness(), geom, fev, feu )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_phiExpr( expr.levelsetPhiExpr(), geom, fev )
            , M_epsilon( expr.thickness(), geom, fev )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_phiExpr( expr.levelsetPhiExpr(), geom )
            , M_epsilon( expr.thickness(), geom )
            , M_locRes( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_phiExpr.init( im );
            M_epsilon.init( im );
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
            M_phiExpr.update( geom );
            M_epsilon.update( geom );
            computeHeaviside();
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_phiExpr.update( geom, face );
            M_epsilon.update( geom, face );
            computeHeaviside();
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
        void computeHeaviside()
        {
            for( uint16_type q = 0; q < M_locRes.size(); ++q )
            {
                double phi = M_phiExpr.evalq( 0, 0, q );
                double eps = M_epsilon.evalq( 0, 0, q );
                if( phi < -eps )
                    M_locRes[q](0,0) = 0.;
                else if( phi > eps )
                    M_locRes[q](0,0) = 1.;
                else
                    M_locRes[q](0,0) = 0.5 * ( 1 + phi/eps + 1/M_PI*math::sin(M_PI*phi/eps) );
            }
        }

    private:
        tensor_expr_type M_phiExpr;
        tensor_expr_epsilon_type M_epsilon;
        std::vector<matrix_shape_type> M_locRes;
    };


private:
    expr_levelsetphi_type M_levelsetPhiExpr;
    expr_epsilon_type M_thicknessHeaviside;
};

template<int IMOrder = -1, typename LSExprT, typename EpsExprT>
inline
Expr< LevelsetHeavisideExpr<LSExprT, EpsExprT, IMOrder> >
levelsetHeaviside( LSExprT phi, EpsExprT eps )
{
    typedef LevelsetHeavisideExpr<LSExprT, EpsExprT, IMOrder> lsheaviside_t;
    return Expr< lsheaviside_t >(  lsheaviside_t( phi, eps ) );
}

}
}

#endif
