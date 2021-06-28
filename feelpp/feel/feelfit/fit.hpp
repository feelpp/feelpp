/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Huber <vincent.huber@cemosis.Fr>
Date: 01-04-2016

Copyright (C) 2007-2008 Universite Joseph Fourier (Grenoble I)

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

#include <feel/feelfit/interpolator.hpp>

#ifndef FEELPP_FIT_HPP
#define FEELPP_FIT_HPP 1

namespace Feel
{
namespace vf
{
template<typename ExprT, int InterpOperator>
class Fit
{
public:
    typedef ExprT expression_type;

    // list of mathematical objet that has to be precomputed
    static const size_type context = expression_type::context;
    // idv(T) is terminal
    // idv(T)+idv(T) is not
    // That is for optimisation
    // if not know: false
    static const bool is_terminal = false;//expression_type::is_terminal;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = expression_type::template HasTestFunction<Func>::result;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = expression_type::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = expression_type::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = expression_type::template has_trial_basis<Func>;
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;

    typedef typename expression_type::value_type value_type;
    using evaluate_type = Eigen::Matrix<value_type,1,1>;

    typedef Fit<ExprT,InterpOperator> this_type;

    Fit( expression_type const & expr,
         std::shared_ptr<Interpolator> const& interpolator )
        :
        M_expr( expr ),
        M_interpolator( interpolator )
        {}
    Fit(Fit const & te)
        :
        M_expr( te.M_expr ),
        M_interpolator( te.M_interpolator )
        {}


    //~Fit(){}
    expression_type const& expression() const
    {
          return M_expr;
    }

    Interpolator const& interpolator() const { return *M_interpolator; }
    std::shared_ptr<Interpolator> const& interpolatorPtr() const { return M_interpolator; }

    //! dynamic context
    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_expr ); }

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    evaluate_type
    evaluate( bool parallel, worldcomm_ptr_t const& worldcomm ) const
        {
            if constexpr ( InterpOperator == 0 )
                return evaluate_type::Constant( this->interpolator()( M_expr.evaluate( parallel,worldcomm )(0,0) ) );
            else
                return evaluate_type::Constant( this->interpolator().diff( M_expr.evaluate( parallel,worldcomm )(0,0) ) );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            return Fit<std::decay_t<decltype(M_expr.applySymbolsExpr( se ))>,InterpOperator>( M_expr.applySymbolsExpr( se ), this->interpolatorPtr() );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            return M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        CHECK( false ) << "TODO";
        return *this;
    }


    // geo_t : transformation geométrique
    // basis_i_t : fonctions tests
    // basis_j_t : fonctions trial
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        // type of the expression given in parameter
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        // shape = scalar, vectial, tensorial ...
        typedef typename tensor_expr_type::shape expr_shape;

        //BOOST_MPL_ASSERT_MSG( ( boost::is_same<mpl::int_<expr_shape::M>,mpl::int_<expr_shape::N> >::value ), INVALID_TENSOR_SHOULD_BE_RANK_2_OR_0, ( mpl::int_<expr_shape::M>, mpl::int_<expr_shape::N> ) );
        /// SCALAR : what will be returned by the expression
        typedef Shape<expr_shape::nDim,Scalar,false,false> shape;

        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vector_type;

        // ???
        template <class Args> struct sig
        {
            typedef value_type type;
        };

        // is the expression null ?
        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        // u = TRIAL
        // v = TEST
        //
        // bilinear form
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_exprFit( expr ),
            M_tensor_expr( expr.expression(), geom, fev, feu)
            {}

        // linear form
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_exprFit( expr ),
            M_tensor_expr( expr.expression(), geom, fev )
            {}

        // evaluation
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_exprFit( expr ),
            M_tensor_expr( expr.expression(), geom )
            {}

        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_exprFit( expr ),
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
            {}


        // IM =
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }

        // precompute the part in the bilinear form
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update(geom);
        }
        // precompute the part in the linear form
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        // precompute the part in the evaluate part
        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }
        // precompute the part in the evaluate part for a face
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensor_expr.update( geom, face );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
            }

        // i : local index of basis test function
        // j : local index of basis trial function
        // c1 : id x, y  z of the first component (= 0 in scalar case)
        // c2 : id x, y  z of the second component (= 0 in scalar and vectorial case)
        // q : id of the current point in the geomap (basically, quadrature point or interpolation point)

        //local assembling - bilinear form a(i,j)
        // value_type
        // evalij( uint16_type i, uint16_type j ) const
        // {
        //     CHECK( false ) << "not implemented";
        //     return 1.; //M_tensor_expr.evalij( i, j );
        // }

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
            return this->evalq( c1,c2,q, mpl::int_<InterpOperator>() );
        }

    private :
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> /**/ ) const
            {
                return M_exprFit.interpolator()( M_tensor_expr.evalq(c1, c2, q) );
            }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> /**/ ) const
            {
                return M_exprFit.interpolator().diff( M_tensor_expr.evalq(c1, c2, q) );
            }

    private :
        this_type const& M_exprFit;
        tensor_expr_type M_tensor_expr;
    };
    /// end of tensor
private:
    expression_type M_expr;
    std::shared_ptr<Interpolator> M_interpolator;
};

/**
 * \brief Fit
 **/

template<typename ExprT>
inline
Expr< Fit<ExprT,0> >
fit( ExprT const& v,
     std::string const& dataFile = soption("fit.datafile"),
     std::string const& abscissa = "", std::string const& ordinate = "",
     std::string const& type = soption("fit.type"),
     WorldComm const& worldComm = Environment::worldComm() )
{
    auto itFindType = InterpolationTypeMap.find( type );
    CHECK( itFindType != InterpolationTypeMap.end() ) << "invalid interpolator type " << type;
    InterpolationType interpolatorEnumType = itFindType->second;
    LOG(INFO) << "Fit "<< dataFile << " with " << type;
    typedef Fit<ExprT,0> fit_t;
    return Expr< fit_t >(  fit_t( v, Interpolator::New( interpolatorEnumType, dataFile, abscissa, ordinate, worldComm ) ) );
}

template<typename ExprT>
inline
Expr< Fit<ExprT,0> >
fit( ExprT const& v,
     std::shared_ptr<Interpolator> const& interpolator )
{
    typedef Fit<ExprT,0> fit_t;
    return Expr< fit_t >(  fit_t( v, interpolator ) );
}


template<typename ExprT>
inline
Expr< Fit<ExprT,1> >
fitDiff( ExprT const& v,
         std::string const& dataFile = soption("fit.datafile"),
         std::string const& abscissa = "", std::string const& ordinate = "",
         std::string const& type = soption("fit.type"),
         WorldComm const& worldComm = Environment::worldComm() )
{
    auto itFindType = InterpolationTypeMap.find( type );
    CHECK( itFindType != InterpolationTypeMap.end() ) << "invalid interpolator type " << type;
    InterpolationType interpolatorEnumType = itFindType->second;
    LOG(INFO) << "Fit Diff"<< dataFile << " with " << type;
    typedef Fit<ExprT,1> fit_t;
    return Expr< fit_t >(  fit_t( v, Interpolator::New( interpolatorEnumType, dataFile, abscissa, ordinate, worldComm ) ) );
}

template<typename ExprT>
inline
Expr< Fit<ExprT,1> >
fitDiff( ExprT const& v,
         std::shared_ptr<Interpolator> const& interpolator )
{
    typedef Fit<ExprT,1> fit_t;
    return Expr< fit_t >(  fit_t( v, interpolator ) );
}


}
}

#endif //FEELPP_FIT_HPP
