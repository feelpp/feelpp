/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-03-05

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file symm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-03-05
 */
/**
   \file symm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-02-10
 */
#ifndef __Symm_H
#define __Symm_H 1

namespace Life
{
namespace vf
{
/// \cond detail
/**
 * \class Symm
 * \brief symm of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprT>
class Symm
{
public:

    static const size_type context = ExprT::context|vm::SYMM;
    static const bool is_symmetric = true;
    static const int pattern = vm::SYMM;

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
    typedef Symm<ExprT> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Symm( expression_type const & __expr )
        :
        _M_expr( __expr )
    {}
    Symm( Symm const & te )
        :
        _M_expr( te._M_expr )
    {}
    ~Symm()
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

    expression_type const& expression() const { return _M_expr; }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;

        template <class Args> struct sig { typedef value_type type; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q, mpl::bool_<PatternContext == pattern>() );
        }

    private:

        /*
         * evalijq
         */
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true> ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<false> ) const
        {
            return value_type(0.0);
        }

    private:
        tensor_expr_type _M_tensor_expr;
    };

private:
    mutable expression_type  _M_expr;
};
/// \endcond

/**
 * \brief symmetric expression
 */
template<typename ExprT>
inline
Expr< Symm<ExprT> >
symm( ExprT v )
{
    typedef Symm<ExprT> symm_t;
    return Expr< symm_t >(  symm_t( v ) );
}

/// \cond detail
//
// unsymm
//
/**
 * \class Unsymm
 * \brief unsymm of a matrix
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename ExprT>
class Unsymm
{
public:

    static const size_type context = ExprT::context|vm::UNSYMM;
    static const bool is_unsymmetric = true;
    static const int pattern = vm::UNSYMM;

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
    typedef Unsymm<ExprT> this_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit Unsymm( expression_type const & __expr )
        :
        _M_expr( __expr )
    {}
    Unsymm( Unsymm const & te )
        :
        _M_expr( te._M_expr )
    {}
    ~Unsymm()
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

    expression_type const& expression() const { return _M_expr; }

    //@}

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;

        typedef typename tensor_expr_type::shape shape;


        template <class Args> struct sig { typedef value_type type; };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            _M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            _M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            _M_tensor_expr( expr.expression(), geom )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            _M_tensor_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            _M_tensor_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            _M_tensor_expr.update( geom );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalij( IndexI const& i, IndexJ const& j ) const
        {
            return _M_tensor_expr.evalij( i, j );
        }

        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }

        template<typename IndexI, typename IndexJ, int PatternContext>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q, mpl::bool_<PatternContext == pattern>() );
        }

    private:

        /*
         * evalijq
         */
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true> ) const
        {
            return _M_tensor_expr.evalijq( i, j, c1, c2, q );
        }
        template<typename IndexI, typename IndexJ>
        value_type
        evalijq( IndexI const& i, IndexJ const& j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<false> ) const
        {
            return value_type(0.0);
        }

    private:
        tensor_expr_type _M_tensor_expr;
    };

private:
    mutable expression_type  _M_expr;
};
/// \endcond

/**
 * \brief unsymmetric expression
 */
template<typename ExprT>
inline
Expr< Unsymm<ExprT> >
unsymm( ExprT v )
{
    typedef Unsymm<ExprT> unsymm_t;
    return Expr< unsymm_t >(  unsymm_t( v ) );
}

}
}
#endif /* __Unsymm_H */

