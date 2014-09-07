/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-07

  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file twovalued.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-07
 */
#ifndef __TwoValued_H
#define __TwoValued_H 1


#include <boost/fusion/container/generation/make_map.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/*!
  \class SumvExpr
  \brief Sumvpose expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT>
class SumvExpr
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

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


    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    typedef SumvExpr<ExprT> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit SumvExpr( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~SumvExpr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
        //BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,Geo_t ));

        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
               vf::detail::gmc<1>,
               vf::detail::gmc<0> >::type gmc1;

        typedef typename fusion::result_of::value_at_key<Geo_t,vf::detail::gmc<0> >::type left_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,vf::detail::gmc<0> >::type::element_type left_gmc_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmc1 >::type right_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmc1 >::type::element_type right_gmc_type;

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;

        typedef typename expression_type::template tensor<map_left_gmc_type, Basis_i_t, Basis_j_t> left_tensor_expr_type;
        typedef typename expression_type::template tensor<map_right_gmc_type, Basis_i_t, Basis_j_t> right_tensor_expr_type;
        typedef typename left_tensor_expr_type::value_type value_type;

        typedef typename left_tensor_expr_type::shape shape;

        struct is_zero
        {
            static const bool value = false;
        };

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map, fev, feu ),
            M_tensor_expr_right( expr.expression(), M_right_map, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map, fev ),
            M_tensor_expr_right( expr.expression(), M_right_map, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map ),
            M_tensor_expr_right( expr.expression(), M_right_map )
        {
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr_left.init( im );
            M_tensor_expr_right.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map );
            M_tensor_expr_right.update( M_right_map );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map );
            M_tensor_expr_right.update( M_right_map );
        }
        void update( Geo_t const& geom )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map );
            M_tensor_expr_right.update( M_right_map );


        }
        void update( Geo_t const& geom, uint16_type face )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map, face );
            M_tensor_expr_right.update( M_right_map, face  );


        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr_left.evalij( i, j ) + M_tensor_expr_right.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensor_expr_left.evalq( c1, c2, q ) + M_tensor_expr_right.evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_tensor_expr_left.evalq( c1, c2, q, mpl::int_<PatternContext>() ) +
                   M_tensor_expr_right.evalq( c1, c2, q, mpl::int_<PatternContext>() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            //DVLOG(2) << "sumv_left= " << M_tensor_expr_left.evalq( q, c1, c2 ) << "\n"
            //<< "sumv_right= " << M_tensor_expr_right.evalq( q, c1, c2 ) << "\n";
            Feel::detail::ignore_unused_variable_warning( i );
            value_type resl = M_tensor_expr_left.evalq( c1, c2, q );
            value_type resr = M_tensor_expr_right.evalq( c1, c2, q );
            value_type res = resl + resr;
            //DVLOG(2) << "resl( " << i << "," << c1 << "," << c2 << "," << q << ")=" << resl << "\n";
            //DVLOG(2) << "resr( " << i << "," << c1 << "," << c2 << "," << q << ")=" << resr << "\n";
            //DVLOG(2) << "res( " << i << "," << c1 << "," << c2 << "," << q << ")=" << res << "\n";
            return res;
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            //DVLOG(2) << "sumv_left= " << M_tensor_expr_left.evalq( q, c1, c2 ) << "\n"
            //<< "sumv_right= " << M_tensor_expr_right.evalq( q, c1, c2 ) << "\n";
            return M_tensor_expr_left.evalq( c1, c2, q )+M_tensor_expr_right.evalq( c1, c2, q );
        }


        left_gmc_ptrtype M_gmc_left;
        right_gmc_ptrtype M_gmc_right;
        map_left_gmc_type M_left_map;
        map_right_gmc_type M_right_map;
        left_tensor_expr_type M_tensor_expr_left;
        right_tensor_expr_type M_tensor_expr_right;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //@}

protected:

private:

    expression_type  M_expr;
};



/*!
  \class SumExpr
  \brief Sumpose expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT, int Side>
class SumExpr
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = false;
    };


    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    typedef SumExpr<ExprT,Side> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit SumExpr( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~SumExpr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
        //BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,Geo_t ));

        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
               vf::detail::gmc<1>,
               vf::detail::gmc<0> >::type gmc1;

        typedef typename mpl::if_<typename fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        template<typename T>
        struct ttt
        {
            typedef typename T::geometric_mapping_context_type type;
        };
        struct yyy
        {
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type >::type::element_type type;
        };
        typedef typename fusion::result_of::value_at_key<Basis_i_t,key_type >::type::element_type e_type;
        typedef typename mpl::if_<boost::is_same<e_type,key_type >,
                mpl::identity<yyy >,
                mpl::identity<ttt<e_type> > >::type::type::type gmc_type;

        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;


        typedef fusion::map<fusion::pair<key_type, gmc_ptrtype> > map_gmc_type;

        typedef typename expression_type::template tensor<map_gmc_type, Basis_i_t, Basis_j_t> tensor_expr_type;

        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = !fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >::value;
            //static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& fev, Basis_j_t const& /*feu*/ )
            :
            M_gmc( fusion::at_key<key_type >( fev )->gmContext() ),
            M_map( fusion::make_map<key_type >( M_gmc ) ),
            M_tensor_expr( expr.expression(), M_map, fev )
        {
            DVLOG(2) << "expr SumExpr<" << Side << "> is_zero " << is_zero::value << "\n";
        }

        tensor( this_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& fev )
            :
            M_gmc( fusion::at_key<key_type >( fev )->gmContext() ),
            M_map( fusion::make_map<key_type >( M_gmc ) ),
            M_tensor_expr( expr.expression(), M_map, fev )
        {
            DVLOG(2) << "expr SumExpr is_zero " << is_zero::value << "\n";
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& /*feu*/ )
        {
            update( geom, fev, fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >() );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            update( geom, fev, fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >() );
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& fev, mpl::true_ )
        {
            M_gmc = fusion::at_key<vf::detail::gmc<Side> >( fev )->gmContext();
            M_map = fusion::make_map<vf::detail::gmc<Side> >( M_gmc );
            M_tensor_expr.update(  M_map, fev );
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& /*fev*/, mpl::false_ )
        {

        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( j );
            return evaliq( i, c1, c2, q, fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >() );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( j );
            return evaliq( i, c1, c2, q, fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >() );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evaliq( i, c1, c2, q, fusion::result_of::has_key<Basis_i_t,vf::detail::gmc<Side> >() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::true_ ) const
        {
            return M_tensor_expr.evaliq( i, c1, c2, q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::false_ ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            Feel::detail::ignore_unused_variable_warning( q );
            return value_type( 0 );
        }

        gmc_ptrtype M_gmc;
        map_gmc_type M_map;
        tensor_expr_type M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //@}

protected:

private:

    expression_type  M_expr;
};

/*!
  \class SumTExpr
  \brief SumTpose expression

  @author Christophe Prud'homme
  @see
*/
template<typename ExprT, int Side>
class SumTExpr
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    typedef SumTExpr<ExprT,Side> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit SumTExpr( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~SumTExpr()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
        //BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,Geo_t ));

        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
               vf::detail::gmc<1>,
               vf::detail::gmc<0> >::type gmc1;

        typedef typename mpl::if_<typename fusion::result_of::has_key<Basis_j_t,vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;

        template<typename T>
        struct ttt
        {
            typedef typename T::geometric_mapping_context_type type;
        };
        struct yyy
        {
            typedef typename fusion::result_of::value_at_key<Geo_t,key_type >::type::element_type type;
        };
        typedef typename fusion::result_of::value_at_key<Basis_j_t,key_type >::type::element_type e_type;
        typedef typename mpl::if_<boost::is_same<e_type,key_type >,
                mpl::identity<yyy >,
                mpl::identity<ttt<e_type> > >::type::type::type gmc_type;

        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

        typedef fusion::map<fusion::pair<key_type, gmc_ptrtype> > map_gmc_type;

        typedef typename expression_type::template tensor<map_gmc_type, Basis_i_t, Basis_j_t> tensor_expr_type;

        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = !fusion::result_of::has_key<Basis_j_t,vf::detail::gmc<Side> >::value;
            //static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& /*geom*/, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc( fusion::at_key<key_type>( feu )->gmContext() ),
            M_map( fusion::make_map<key_type>( M_gmc ) ),
            M_tensor_expr( expr.expression(), M_map, fev, feu )
        {
            DVLOG(2) << "expr SumTExpr<" << Side << "> is_zero " << is_zero::value << "\n";
        }

        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            update( geom, fev, feu, fusion::result_of::has_key<Basis_j_t,vf::detail::gmc<Side> >() );
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& fev, Basis_j_t const& feu, mpl::true_ )
        {
            M_gmc = fusion::at_key<key_type>( feu )->gmContext();
            M_map = fusion::make_map<key_type>( M_gmc );
            M_tensor_expr.update(  M_map, fev, feu );
        }
        void update( Geo_t const& /*geom*/, Basis_i_t const& fev, Basis_j_t const& feu, mpl::false_ )
        {
            Feel::detail::ignore_unused_variable_warning( fev );
            Feel::detail::ignore_unused_variable_warning( feu );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return M_tensor_expr.evalij( i, j );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalijq( i, j, c1,c2, q, fusion::result_of::has_key<Basis_j_t,vf::detail::gmc<Side> >() );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1,c2, q, mpl::int_<PatternContext>(), fusion::result_of::has_key<Basis_j_t,vf::detail::gmc<Side> >() );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::true_ ) const
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::false_ ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            Feel::detail::ignore_unused_variable_warning( q );
            return value_type( 0 );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext>, mpl::true_ ) const
        {
            return M_tensor_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext>, mpl::false_ ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            Feel::detail::ignore_unused_variable_warning( q );
            return value_type( 0 );
        }

    private:

        gmc_ptrtype M_gmc;
        map_gmc_type M_map;
        tensor_expr_type M_tensor_expr;
    };

    //@}

    /** @name Accessors
     */
    //@{

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //@}

protected:

private:

    expression_type  M_expr;
};
/// \endcond


/// \cond detail
/**
 * \class FaceExprV
 * \brief FaceExprV expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprT, double func( double,double )>
class FaceExprV
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = ExprT::imorder;
    static const bool imIsPoly = ExprT::imIsPoly;

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

    /** @name Typedefs
     */
    //@{

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    typedef FaceExprV<ExprT, func> this_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit FaceExprV( expression_type const & __expr )
        :
        M_expr( __expr )
    {}
    ~FaceExprV()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename Geo_t, typename Basis_i_t = fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<vf::detail::gmc<0> > > >, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
        //BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,Geo_t ));

        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
               vf::detail::gmc<1>,
               vf::detail::gmc<0> >::type gmc1;

        typedef typename fusion::result_of::value_at_key<Geo_t,vf::detail::gmc<0> >::type left_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,vf::detail::gmc<0> >::type::element_type left_gmc_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmc1 >::type right_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,gmc1 >::type::element_type right_gmc_type;

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;

        typedef typename expression_type::template tensor<map_left_gmc_type, Basis_i_t, Basis_j_t> left_tensor_expr_type;
        typedef typename expression_type::template tensor<map_right_gmc_type, Basis_i_t, Basis_j_t> right_tensor_expr_type;
        typedef typename left_tensor_expr_type::value_type value_type;

        typedef typename left_tensor_expr_type::shape shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map, fev, feu ),
            M_tensor_expr_right( expr.expression(), M_right_map, fev, feu )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map, fev ),
            M_tensor_expr_right( expr.expression(), M_right_map, fev )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_gmc_left( fusion::at_key<vf::detail::gmc<0> >( geom ) ),
            M_gmc_right( fusion::at_key<gmc1 >( geom ) ),
            M_left_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_left ) ),
            M_right_map( fusion::make_map<vf::detail::gmc<0> >( M_gmc_right ) ),
            M_tensor_expr_left( expr.expression(), M_left_map ),
            M_tensor_expr_right( expr.expression(), M_right_map )
        {
        }
        template<typename IM>
        void init( IM const& im )
        {
            M_tensor_expr_left.init( im );
            M_tensor_expr_right.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map );
            M_tensor_expr_right.update( M_right_map );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map );
            M_tensor_expr_right.update( M_right_map );
        }
        void update( Geo_t const& geom )
        {
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            typedef mpl::int_<fusion::result_of::template size<Geo_t>::type::value> map_size;
            FEELPP_ASSERT( map_size::value == 2 )( map_size::value ).error( "invalid map size (should be 2)" );

            M_gmc_left = fusion::at_key<vf::detail::gmc<0> >( geom );
            M_gmc_right =  fusion::at_key<gmc1 >( geom );
            FEELPP_ASSERT( M_gmc_left != M_gmc_right )( M_gmc_left->id() )( M_gmc_right->id() ).error( "same geomap, something is wrong" );

            M_left_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_left );
            M_right_map = fusion::make_map<vf::detail::gmc<0> >( M_gmc_right );
            M_tensor_expr_left.update(  M_left_map, face );
            M_tensor_expr_right.update( M_right_map, face );


        }

        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return func( M_tensor_expr_left.evalij( i, j ), M_tensor_expr_right.evalij( i, j ) );
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return func( M_tensor_expr_left.evalijq( i, j, c1, c2, q ),
                         M_tensor_expr_right.evalijq( i, j, c1, c2, q ) );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return func( M_tensor_expr_left.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() ),
                         M_tensor_expr_right.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() ) );
        }


        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return func( M_tensor_expr_left.evaliq( i, c1, c2, q ),
                         M_tensor_expr_right.evaliq( i, c1, c2, q ) );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return func( M_tensor_expr_left.evalq( c1, c2, q ),
                         M_tensor_expr_right.evalq( c1, c2, q ) );
        }


        left_gmc_ptrtype M_gmc_left;
        right_gmc_ptrtype M_gmc_right;
        map_left_gmc_type M_left_map;
        map_right_gmc_type M_right_map;
        left_tensor_expr_type M_tensor_expr_left;
        right_tensor_expr_type M_tensor_expr_right;
    };

    //@}

    /** @name Accessors
     */
    //@{

    bool isSymetric() const
    {
        return M_expr.isSymetric();
    }

    expression_type const& expression() const
    {
        return M_expr;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    //@}

protected:

private:

    expression_type  M_expr;
};

namespace detail
{
inline double max( double a, double b )
{
    return std::max( a, b );
}
inline double min( double a, double b )
{
    return std::min( a, b );
}

inline double left( double a, double /*b*/ )
{
    return a;
}
inline double right( double /*a*/, double b )
{
    return b;
}
}

/// \endcond

template<typename ExprT>
inline
auto
sumv( ExprT const& v ) -> decltype( expr( SumvExpr<ExprT>( v ) ) )
{
    typedef SumvExpr<ExprT> sumv_t;
    return expr(  sumv_t( v ) );
}

template<typename ExprT>
inline
auto
leftface( ExprT const& v ) -> decltype( expr(  SumExpr<ExprT,0>( v ) ) )
{
    typedef SumExpr<ExprT,0> sumT_t;
    return expr(  sumT_t( v ) );
}

template<typename ExprT>
inline
auto
rightface( ExprT const& v ) -> decltype( expr( SumExpr<ExprT,1>( v ) ) )
{
    typedef SumExpr<ExprT,1> sumT_t;
    return expr(  sumT_t( v ) );
}


/**
 *
 */
template<typename ExprT>
inline
Expr< SumTExpr<ExprT, 0> >
leftfacet( ExprT const& v )
{
    typedef SumTExpr<ExprT,0> sumT_t;
    return expr(  sumT_t( v ) );
}

/**
 *
 */
template<typename ExprT>
inline
Expr< SumTExpr<ExprT, 1> >
rightfacet( ExprT const& v )
{
    typedef SumTExpr<ExprT,1> sumT_t;
    return expr(  sumT_t( v ) );
}

#if 0
/**
 *
 */
template<typename U>
auto
jump( U const& u ) -> decltype( leftface( u*N() ) + rightface( u*N() ) )
{
    return ( leftface( u*N() ) + rightface( u *N() ) );
}
/**
 *
 */
template<typename U>
auto
jumpt( U const& u ) -> decltype( leftfacet( u*N() ) + rightfacet( u*N() ) )
{
    return ( leftfacet( u*N() ) + rightfacet( u*N() ) );
}
/**
 *
 */
template<typename U>
auto
average( U const& u ) -> decltype( 0.5*( leftface( u )+rightface( u ) ) )
{
    return 0.5*( leftface( u )+rightface( u ) );
}
/**
 *
 */
template<typename U>
auto
averaget( U const&u ) -> decltype( 0.5*( leftfacet( u )+rightfacet( u ) ) )
{
    return 0.5*( leftfacet( u )+rightfacet( u ) );
}
#else
#define jump( u )  (::Feel::vf::leftface((u)*::Feel::vf::N())+::Feel::vf::rightface((u)*::Feel::vf::N()))
#define average( u ) (.5*(::Feel::vf::leftface((u))+::Feel::vf::rightface((u))))

#define jumpt( u )  (::Feel::vf::leftfacet((u)*::Feel::vf::N())+::Feel::vf::rightfacet((u)*::Feel::vf::N()))
#define averaget( u ) (.5*(::Feel::vf::leftfacet((u))+::Feel::vf::rightfacet((u))))

#endif

template<typename ExprT>
inline
Expr< FaceExprV<ExprT,vf::detail::max > >
maxface( ExprT const& v )
{
    typedef FaceExprV<ExprT,vf::detail::max > maxface_t;
    return Expr< maxface_t >(  maxface_t( v ) );
}
template<typename ExprT>
inline
Expr< FaceExprV<ExprT,vf::detail::min > >
minface( ExprT const& v )
{
    typedef FaceExprV<ExprT,vf::detail::min > minface_t;
    return Expr< minface_t >(  minface_t( v ) );
}


template<typename ExprT>
inline
Expr< FaceExprV<ExprT,vf::detail::left> >
leftfacev( ExprT const& v )
{
    typedef FaceExprV<ExprT,vf::detail::left> leftface_t;
    return Expr< leftface_t >(  leftface_t( v ) );
}
template<typename ExprT>
inline
Expr< FaceExprV<ExprT,vf::detail::right> >
rightfacev( ExprT const& v )
{
    typedef FaceExprV<ExprT,vf::detail::right> rightface_t;
    return Expr< rightface_t >(  rightface_t( v ) );
}

#if 0
/**
 *
 */
template<typename U>
auto
jumpv( U const& u ) -> decltype( ( leftfacev( u*N() )+rightfacev( u*N() ) ) )
{
    return leftfacev( u*N() )+rightfacev( u*N() );
}
/**
 *
 */
template<typename U>
auto
averagev( U const& u ) -> decltype( .5*( leftfacev( u )+rightfacev( u ) ) )
{
    return .5*( leftfacev( u )+rightfacev( u ) );
}
#else
#define jumpv( u )  (::Feel::vf::leftfacev((u)*::Feel::vf::N())+::Feel::vf::rightfacev((u)*::Feel::vf::N()))
/**
 *
 */
#define averagev( u ) (.5*(::Feel::vf::leftfacev((u))+::Feel::vf::rightfacev((u))))

#endif
} // vf
} // Feel
#endif /* __TwoValued_H */
