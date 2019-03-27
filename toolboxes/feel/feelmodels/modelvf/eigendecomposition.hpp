/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2019-03-22

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
#ifndef FEELPP_VF_EigenDecomposition_HPP
#define FEELPP_VF_EigenDecomposition_HPP 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class EigenDecomposition
 * \brief Modify a tensor by taking absolute value of eigen value
 *
 */
template<typename ExprT>
class EigenDecomposition
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;


    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    typedef EigenDecomposition<ExprT> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit EigenDecomposition( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    EigenDecomposition( EigenDecomposition const & te )
        :
        M_expr( te.M_expr )
    {}
    ~EigenDecomposition()
    {}

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
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    //template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptr<vf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape expr_shape;
        BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>,mpl::int_<expr_shape::N> >::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, ( mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N> ) );
        typedef expr_shape shape;

        typedef Eigen::Matrix<value_type,shape::M,shape::N> eval_matrix_type;
        typedef Eigen::Matrix< eval_matrix_type, Eigen::Dynamic,1> vector_type;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu ),
            M_eval( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev ),
            M_eval( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom ),
            M_eval( vf::detail::ExtractGm<Geo_t>::get( geom )->nPoints() )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
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
            M_tensor_expr.update( geom );
            this->computeEigenDecomposition();
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
            this->computeEigenDecomposition();
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
            return M_eval(q)(c1,c2);
        }

    private:
        void
        computeEigenDecomposition()
            {
                for( int q = 0; q < M_eval.rows(); ++q )
                {
                    eval_matrix_type A;
                    for ( int i=0;i<shape::M;++i )
                        for ( int j=0;j<shape::N;++j )
                            A(i,j) = M_tensor_expr.evalq( i, j, q );
                    Eigen::EigenSolver<eval_matrix_type> es(A);
                    eval_matrix_type D = es.pseudoEigenvalueMatrix();
                    eval_matrix_type V = es.pseudoEigenvectors();
                    M_eval(q) = V * D.cwiseAbs() * V.inverse();
                }
            }

    private:
        tensor_expr_type M_tensor_expr;
        vector_type M_eval;
    };

private:
    mutable expression_type  M_expr;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */
template<typename ExprT>
inline
Expr< EigenDecomposition<ExprT> >
eigenDecomposition( ExprT v )
{
    typedef EigenDecomposition<ExprT> ed_t;
    return Expr< ed_t >( ed_t( v ) );
}

}
}
#endif /* __Det_H */
