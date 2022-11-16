/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-04-26

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
   \file det.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-04-26
 */
#ifndef __FEELPP_VF_Det_H
#define __FEELPP_VF_Det_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Det
 * \brief det of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template <typename ExprT>
class Det : public ExprDynamicBase
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
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;
    typedef Det<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Det( expression_type const& __expr )
        : super( Feel::vf::dynamicContext( __expr ) ),
          M_expr( __expr )
    {
    }
    Det( Det const& te )
        : super( te ),
          M_expr( te.M_expr )
    {
    }
    ~Det() = default;

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
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            auto eval = M_expr.evaluate(p,worldcomm);
            CHECK( eval.rows() == eval.cols() ) << "only square matrix";
            return evaluate_type::Constant( eval.determinant() );
       }
    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptr<vf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>>::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, (mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N>));
        typedef Shape<expr_shape::nDim, Scalar, false, false> shape;

        typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> vector_type;

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
              M_det( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            : M_tensor_expr( expr.expression(), geom, fev ),
              M_det( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            : M_tensor_expr( expr.expression(), geom ),
              M_det( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
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
            computeDet( mpl::int_<expr_shape::nDim>() );
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
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q ) const
        {
            return M_det( q );
        }

      private:
        void
            computeDet( mpl::int_<1> )
        {
            for ( int q = 0; q < M_det.rows(); ++q )
            {
                M_det( q ) = M_tensor_expr.evalq( 0, 0, q );
            }
        }
        void
            computeDet( mpl::int_<2> )
        {
            for ( int q = 0; q < M_det.rows(); ++q )
            {
                double a = M_tensor_expr.evalq( 0, 0, q );
                double b = M_tensor_expr.evalq( 0, 1, q );
                double c = M_tensor_expr.evalq( 1, 0, q );
                double d = M_tensor_expr.evalq( 1, 1, q );
                M_det( q ) = a * d - b * c;
            }
        }
        void
            computeDet( mpl::int_<3> )
        {

            for ( int q = 0; q < M_det.rows(); ++q )
            {
                double a = M_tensor_expr.evalq( 0, 0, q );
                double b = M_tensor_expr.evalq( 0, 1, q );
                double c = M_tensor_expr.evalq( 0, 2, q );
                double d = M_tensor_expr.evalq( 1, 0, q );
                double e = M_tensor_expr.evalq( 1, 1, q );
                double f = M_tensor_expr.evalq( 1, 2, q );
                double g = M_tensor_expr.evalq( 2, 0, q );
                double h = M_tensor_expr.evalq( 2, 1, q );
                double i = M_tensor_expr.evalq( 2, 2, q );
                M_det( q ) = a * e * i + b * f * g + c * d * h - ( c * e * g + b * d * i + a * f * h );
            }
        }

      private:
        tensor_expr_type M_tensor_expr;
        vector_type M_det;
    };

  private:
    mutable expression_type M_expr;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */
template <typename ExprT>
inline Expr<Det<ExprT>>
det( ExprT v )
{
    typedef Det<ExprT> det_t;
    return Expr<det_t>( det_t( v ) );
}

} // namespace vf
} // namespace Feel
#endif /* __Det_H */
