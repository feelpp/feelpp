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
#ifndef __FEELPP_VF_VONMISES_H
#define __FEELPP_VF_VONMISES_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class VonMises
 * \brief VonMises of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprT>
class VonMises : public ExprDynamicBase
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
    typedef value_type evaluate_type;
    typedef VonMises<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit VonMises( expression_type const& __expr )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr )
    {
    }
    VonMises( VonMises const& te )
        : super( te ), M_expr( te.M_expr )
    {
    }
    ~VonMises() = default;

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

    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptrvf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( std::is_same_v<mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>> ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, (mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>));
        typedef Shape<expr_shape::nDim, Scalar, false, false> shape;

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
            : M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom )
        {
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
                if constexpr ( expr_shape::N == 1 )
                    return M_tensor_expr.evalq( 0, 0, q );
                else if constexpr ( expr_shape::N == 2 )
                    return std::sqrt( std::pow( M_tensor_expr.evalq( 0, 0, q ), 2 ) -
                                      M_tensor_expr.evalq( 0, 0, q )*M_tensor_expr.evalq( 1, 1, q )+
                                      std::pow( M_tensor_expr.evalq( 1, 1, q ), 2 ) +
                                      3 * std::pow( M_tensor_expr.evalq( 0, 1, q ), 2 ) 
                                      );
                else
                    return std::sqrt( ( std::pow( M_tensor_expr.evalq( 0, 0, q )-M_tensor_expr.evalq( 1, 1, q ), 2 ) +
                                        std::pow( M_tensor_expr.evalq( 1, 1, q )-M_tensor_expr.evalq( 2, 2, q ), 2 ) +
                                        std::pow( M_tensor_expr.evalq( 2, 2, q )-M_tensor_expr.evalq( 0, 0, q ), 2 ) +
                                        6*(std::pow( M_tensor_expr.evalq( 0, 1, q ),2)+
                                           std::pow( M_tensor_expr.evalq( 0, 2, q ),2)+
                                           std::pow( M_tensor_expr.evalq( 1, 2, q ),2) ) )/2
                                      );
        }
      private:
        tensor_expr_type M_tensor_expr;
    };

  private:
    mutable expression_type M_expr;
};
/// \endcond

/**
 * \brief compute the VonMises yield criterion
 */
template <typename ExprT>
inline Expr<VonMises<ExprT>>
vonmises( ExprT v, std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr )
{
    typedef VonMises<ExprT> vonmises_t;
    return Expr<vonmises_t>( vonmises_t( v ) );
}

} // namespace vf
} // namespace Feel
#endif /* __Inv_H */
