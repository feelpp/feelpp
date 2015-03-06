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
   \file placeholder.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_PLACEHOLDER_HPP
#define FEELPP_VF_PLACEHOLDER_HPP 1

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

template<bool B, typename Elt> struct GetElementType2 {};
template<typename Elt> struct GetElementType2<true,Elt> { typedef typename Elt::expression_type type; };
template<typename Elt> struct GetElementType2<false,Elt> { typedef Elt type; };



class PlaceHolderBase {};

template<int I = FIRST>
class PlaceHolder : public PlaceHolderBase
{
public:
    static const int kind = I;

    static const size_type context = 0;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = false;

    typedef PlaceHolder<I> this_type;
    typedef this_type functionspace_type;
    typedef this_type reference_element_type;
    typedef this_type* fe_ptrtype;
    typedef this_type fe_type;
    typedef this_type basis_type;
    typedef this_type mortar_fe_type;
    typedef this_type geoelement_type;
    typedef this_type gm_type;
    typedef this_type polyset_type;

    static const uint16_type nDim = 0;
    static const uint16_type nOrder = 0;
    static const uint16_type rank = 0;
    static const uint16_type nComponents = 0;
    static const uint16_type nComponents1 = 0;
    static const uint16_type nComponents2 = 0;

    struct PreCompute {};
    template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g = context_v>
    struct Context {};
    static const bool is_mortar = false;

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
        static const bool is_expression = boost::is_base_of<ExprBase,_type>::value;
        typedef boost::is_base_of<ExprBase,_type> is_expression_type;
        typedef typename GetElementType2<is_expression,_type>::type type;

        template<typename ExprT>
        static type impl( mpl::bool_<true>, ExprT const& e  )
            {
                std::cout << "return expression\n";
                return e.expression();
            }
        template<typename ExprT>
        static type impl( mpl::bool_<false>, ExprT const& e  )
            {
                std::cout << "return element\n";
                return e;
            }

        template<typename ExprT>
        static type expr( ExprT const& e  )
            {
                return impl( is_expression_type(), e );
            }

        /**
         * two args
         */
        template<typename ExprT1, typename ExprT2>
        static type impl( mpl::bool_<true>, mpl::int_<FIRST>, ExprT1 const& e1, ExprT2 const& e2  )
            {
                return e1.expression();
            }
        template<typename ExprT1, typename ExprT2>
        static type impl( mpl::bool_<true>, mpl::int_<SECOND>, ExprT1 const& e1, ExprT2 const& e2  )
            {
                return e2.expression();
            }

        template<typename ExprT1, typename ExprT2>
        static type impl( mpl::bool_<false>, mpl::int_<FIRST>, ExprT1 const& e1, ExprT2 const& e2  )
            {
                return e1;
            }
        template<typename ExprT1, typename ExprT2>
        static type impl( mpl::bool_<false>, mpl::int_<SECOND>, ExprT1 const& e1, ExprT2 const& e2  )
            {
                return e2;
            }

        template<typename ExprT1, typename ExprT2>
        static type expr( ExprT1 const& e1, ExprT2 const& e2  )
            {
                return impl( is_expression_type(), mpl::int_<kind>(), e1, e2 );
            }
        /**
         * three args
         */
        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<true>, mpl::int_<FIRST>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e1.expression();
            }
        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<true>, mpl::int_<SECOND>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e2.expression();
            }
        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<true>, mpl::int_<THIRD>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e3.expression();
            }

        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<false>, mpl::int_<FIRST>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e1;
            }
        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<false>, mpl::int_<SECOND>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e2;
            }
        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type impl( mpl::bool_<false>, mpl::int_<THIRD>, ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return e3;
            }

        template<typename ExprT1, typename ExprT2, typename ExprT3>
        static type expr( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3  )
            {
                return impl( is_expression_type(), mpl::int_<kind>(), e1, e2, e3 );
            }
    };

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e )
        {
            return  Lambda<ExprT>::expr( e );
        }

    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2)
        {
            return  Lambda<ExprT1, ExprT2>::expr( e1, e2 );
        }
    template<typename ExprT1, typename ExprT2, typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3 )
        {
            return Lambda<ExprT1, ExprT2, ExprT3>::expr( e1, e2, e3 );
            //return LambdaImpl(e1,e2,e3,mpl::int_<kind>());
        }

    template<typename ExprT>
    typename Lambda<ExprT>::type
    operator()( ExprT const& e ) const
        {
            return  Lambda<ExprT>::expr( e );
            //return typename Lambda<ExprT>::type ( e );
        }

    template<typename ExprT1, typename ExprT2>
    typename Lambda<ExprT1,ExprT2>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2 ) const
        {
            return  Lambda<ExprT1, ExprT2>::expr( e1, e2 );
            //return LambdaImpl(e1,e2,mpl::int_<kind>());
        }

    template<typename ExprT1, typename ExprT2, typename ExprT3>
    typename Lambda<ExprT1,ExprT2,ExprT3>::type
    operator()( ExprT1 const& e1, ExprT2 const& e2, ExprT3 const& e3 ) const
        {
            return  Lambda<ExprT1, ExprT2, ExprT3>::expr( e1, e2, e3 );
            //return LambdaImpl(e1,e2,e3,mpl::int_<kind>());
        }


    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > >,fusion::pair<vf::detail::gmc<1>,boost::shared_ptr<vf::detail::gmc<1> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef PlaceHolder<I> expression_type;
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

using PlaceHolder1 = PlaceHolder<FIRST>;
using PlaceHolder2 = PlaceHolder<SECOND>;
using PlaceHolder3 = PlaceHolder<THIRD>;


} // vf
} // Feel
#endif /* FEELPP_VF_PLACEHOLDER_HPP */
