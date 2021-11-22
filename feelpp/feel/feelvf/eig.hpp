/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2019-12-07

  Copyright (C) 2019 Feel++ Consortium

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
#ifndef __FEELPP_VF_EIG_H
#define __FEELPP_VF_EIG_H 1

#include <feel/feelvf/expr.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Eig
 * \brief Eig of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprT>
class Eig : public ExprDynamicBase
{
  public:
    using super = ExprDynamicBase;
    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template <typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef Eig<ExprT> this_type;

    using evaluate_type = Eigen::Matrix<value_type,
                                        (expression_type::evaluate_type::SizeAtCompileTime == Eigen::Dynamic)? Eigen::Dynamic : expression_type::evaluate_type::RowsAtCompileTime,
                                        1>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Eig( expression_type const& __expr )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr )
    {
    }
    Eig( Eig const& te )
        : super( te ), M_expr( te.M_expr )
    {
    }
    ~Eig() = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //! polynomial order
    uint16_type polynomialOrder() const { return 2 * M_expr.polynomialOrder(); }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return false; }

    expression_type const& expression() const
    {
        return M_expr;
    }

    evaluate_type
    evaluate( bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            auto evalExpr = M_expr.evaluate( p, worldcomm );
            Eigen::SelfAdjointEigenSolver<typename expression_type::evaluate_type> es( evalExpr );
            return es.eigenvalues();
        }

    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptrvf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>>::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, (mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>));
        typedef Shape<expr_shape::nDim, Vectorial, false, false> shape;

        typedef Eigen::Matrix<value_type, expr_shape::M, expr_shape::N> matrix_type;
        typedef Eigen::Matrix<value_type, shape::M, 1> vector_type;
        typedef std::vector<matrix_type> eig_matrix_type;
        typedef std::vector<vector_type> eig_vector_type;

        template <class Args>
        struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            : M_tensor_expr( expr.expression(), geom, fev, feu ),
              M_eig( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() ),
              M_eigv( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev ),
              M_eig( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() ),
              M_eigv( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom ),
              M_eig( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() ),
              M_eigv( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }
        template <typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
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
            M_tensor_expr.update( geom );
            computeEig();
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
            computeEig();
        }

        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }

        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_eigv[q]( c1, c2 );
        }

        Eigen::Map<const vector_type>
        evalq( uint16_type q ) const
        {
            return Eigen::Map<const vector_type>(M_eigv[q].data());
        }

      private:
        void  computeEig()
            {
                for ( int q = 0; q < M_eig.size(); ++q )
                {
                    M_eig[q](0,0) = M_tensor_expr.evalq( 0, 0, q );
                    if constexpr ( expr_shape::N > 1 )
                    {
                        M_eig[q](1,1) = M_tensor_expr.evalq( 1, 1, q );
                        M_eig[q](1,0) = M_tensor_expr.evalq( 1, 0, q );
                        if constexpr ( expr_shape::N > 2 )
                        {
                            M_eig[q](2,2) = M_tensor_expr.evalq( 2, 2, q );
                            M_eig[q](2,0) = M_tensor_expr.evalq( 2, 0, q );
                            M_eig[q](2,1) = M_tensor_expr.evalq( 2, 1, q );
                        }
                    }
                    Eigen::SelfAdjointEigenSolver<matrix_type> es(M_eig[q]);
                    M_eigv[q] = es.eigenvalues();
                }
            }

      private:
        tensor_expr_type M_tensor_expr;
        eig_matrix_type M_eig;
        eig_vector_type M_eigv;
    };

  private:
    mutable expression_type M_expr;
};
/// \endcond

/**
 * \brief eig of the expression tensor (expression is supposed to be self adjoint)
 */
template <typename ExprT>
inline Expr<Eig<ExprT>>
eig( ExprT v, std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    typedef Eig<ExprT> eig_t;
    return Expr<eig_t>( eig_t( v ) );
}

} // namespace vf
} // namespace Feel
#endif /* __Inv_H */
