/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-20

  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
   \file function.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-20
 */
#ifndef __Function_H
#define __Function_H 1

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class Function
 * \brief allow runtime function in expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename Func>
class Function
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;
    static const uint16_type imorder = Func::imorder;
    static const bool imIsPoly = Func::imIsPoly;

    template<typename Funct>
    struct HasTestFunction
    {
        static const bool result = false;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static const bool result = false;
    };

    typedef Func expression_type;
    typedef Function<Func> this_type;
    typedef typename expression_type::value_type value_type;
    typedef value_type evaluate_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Function( expression_type const & fun )
        :
        M_fun( fun )
    {}
    Function( Function const & fun )
        :
        M_fun( fun.M_fun )
    {}
    ~Function()
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

    const expression_type& fun() const
    {
        return M_fun;
    }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::value_type value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<Func::rank>,mpl::int_<0> >,
                mpl::identity<Shape<gmc_type::nDim, Scalar, false, false> >,
                typename mpl::if_<mpl::equal_to<mpl::int_<Func::rank>,mpl::int_<1> >,
                mpl::identity<Shape<gmc_type::nDim, Vectorial, false, false> >,
                mpl::identity<Shape<gmc_type::nDim, Tensor2, false, false> > >::type >::type::type shape;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        template<typename IM>
        void init( IM const& im )
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
            M_gmc =  fusion::at_key<key_type>( geom ).get();
        }

        void update( Geo_t const& geom, uint16_type /*face*/ )
        {
            M_gmc =  fusion::at_key<key_type>( geom ).get();
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun( c1, c2, M_gmc->xReal( q ), M_gmc->unitNormal( q ) );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_fun( c1, c2, M_gmc->xReal( q ), M_gmc->unitNormal( q ) );
        }



        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun( c1, c2, M_gmc->xReal( q ), M_gmc->unitNormal( q ) );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun( c1, c2, M_gmc->xReal( q ), M_gmc->unitNormal( q ) );
        }

        Func M_fun;
        gmc_ptrtype M_gmc;
    };

private:
    mutable expression_type  M_fun;
};
/// \endcond

/**
 * \brief functor enabling function
 *
 */
template<typename Func>
inline
Expr< Function<Func> >
idf( Func f )
{
    typedef Function<Func> func_t;
    return Expr< func_t >(  func_t( f ) );
}

} // vf
} // Feel
#endif /* __Function_H */
