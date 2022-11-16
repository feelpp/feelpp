/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2016-02-10

  Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_FUNCTION2_HPP
#define FEELPP_FUNCTION2_HPP 1

namespace Feel
{
namespace vf
{
/// \cond detail
namespace details
{
/**
 * \class Function2
 * \brief Helper class to allow functors in expressions
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename Func>
class Function2
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;

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
    template<typename Funct>
    static const bool has_test_basis = false;
    template<typename Funct>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    typedef Func expression_type;
    typedef Function2<Func> this_type;
    typedef typename expression_type::value_type value_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Function2( expression_type const & fun )
        :
        M_fun( fun )
    {}
    Function2( Function2 const & fun )
        :
        M_fun( fun.M_fun )
    {}
    ~Function2()
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
    uint16_type polynomialOrder() const { return M_fun.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_fun.isPolynomial(); }

    const expression_type& fun() const
    {
        return M_fun;
    }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::value_type value_type;

        using key_type = key_t<Geo_t>;
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
            tensor( expr, geom )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            tensor( expr, geom )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
            M_fun.init( M_gmc );
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
            M_fun.update( M_gmc );
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun.evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_fun.evalq( c1, c2, q );
        }



        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun.evalq( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_fun.evalq( c1, c2, q );
        }

        Func M_fun;
        gmc_ptrtype M_gmc;
    };

private:
    mutable expression_type  M_fun;
};
} // detail
/// \endcond

/**
 * \brief functor enabling function
 *
 */
template<typename Func>
inline
Expr< Feel::vf::details::Function2<Func> >
idf2( Func f )
{
    typedef Feel::vf::details::Function2<Func> func_t;
    return Expr< func_t >(  func_t( f ) );
}

} // vf
} // Feel
#endif
