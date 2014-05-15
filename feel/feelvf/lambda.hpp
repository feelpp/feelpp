/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file lambda.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_LAMBDA_HPP
#define FEELPP_VF_LAMBDA_HPP 1

namespace Feel
{
namespace vf
{
enum { NONE             = 0x00, // Notice we are using bits as flags here.
       FIRST            = 0x01,
       SECOND           = 0x02,
       THIRD            = 0x04,
       EXCEPTION        = 0x08,
       RETHROW          = 0x10};

class LambdaExprBase {};

template<int I = FIRST>
class LambdaExpr : public LambdaExprBase
{
public:
    static const int kind = I;

    static const size_type context = 0;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = false;

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

    typedef double value_type;
    typedef value_type evaluate_type;

    template<typename TheExpr1, typename TheExpr2 = boost::none_t, typename TheExpr3 = boost::none_t>
    struct Lambda
    {
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<kind>,mpl::int_<FIRST>>,
                                  mpl::identity<TheExpr1>,
                                  typename mpl::if_<mpl::equal_to<mpl::int_<kind>,mpl::int_<SECOND>>,
                                                    mpl::identity<TheExpr2>,
                                                    mpl::identity<TheExpr3>>::type>::type::type _type;
        typedef typename _type::expression_type type;
    };

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e )
        {
            return e.expression();
        }
    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    LambdaImpl( ExprT1 const& e1, ExprT2 const& e2, mpl::int_<FIRST>)
        {
            return e1.expression();
        }
    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    LambdaImpl( ExprT1 const& e1, ExprT2 const& e2, mpl::int_<SECOND>)
        {
            return e2.expression();
        }

    template<typename ExprT1, typename ExprT2,typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    LambdaImpl( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3, mpl::int_<FIRST>)
        {
            return e1.expression();
        }
    template<typename ExprT1, typename ExprT2,typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    LambdaImpl( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3, mpl::int_<SECOND>)
        {
            return e2.expression();
        }
    template<typename ExprT1, typename ExprT2,typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    LambdaImpl( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3, mpl::int_<THIRD>)
        {
            return e3.expression();
        }

    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2)
        {
            return LambdaImpl( e1, e2, mpl::int_<kind>() );
        }
    template<typename ExprT1, typename ExprT2, typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3 )
        {
            return LambdaImpl(e1,e2,e3,mpl::int_<kind>());
        }

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e ) const
        {
            return e.expression();
        }

    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2 ) const
        {
            return LambdaImpl(e1,e2,mpl::int_<kind>());
        }

    template<typename ExprT1, typename ExprT2, typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3 ) const
        {
            return LambdaImpl(e1,e2,e3,mpl::int_<kind>());
        }


    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,boost::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef LambdaExpr<I> expression_type;
        typedef typename expression_type::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<0> >, mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;


        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
        }
        tensor( expression_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& /*fev*/ )
        {
        }
        tensor( expression_type const& expr, Geo_t const& /*geom*/ )
        {
        }
        template<typename IM>
        void init( IM const& /*im*/ )
        {
        }
        void update( Geo_t const&, Basis_i_t const& , Basis_j_t const&  )
        {
        }
        void update( Geo_t const& , Basis_i_t const&  )
        {
        }
        void update( Geo_t const& )
        {
        }
        void update( Geo_t const&, uint16_type )
        {
        }
        template<typename CTX>
        void updateContext( CTX const& ctx )
        {
        }


        value_type
        evalij( uint16_type /*i*/, uint16_type /*j*/ ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return 0;
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/,
                 mpl::int_<PatternContext> ) const
        {
            return 0;
        }

        value_type
        evaliq( uint16_type /*i*/, uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/  ) const
        {
            return 0;
        }
        value_type
        evalq( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type /*q*/ ) const
        {
            return 0;
        }
    };

};

using LambdaExpr1 = LambdaExpr<FIRST>;
using LambdaExpr2 = LambdaExpr<SECOND>;
using LambdaExpr3 = LambdaExpr<THIRD>;

} // vf
} // Feel
#endif /* FEELPP_VF_LAMBDA_HPP */
