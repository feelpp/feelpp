/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-10-15

  Copyright (C) 2012 Université de Strasbourg

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
   \file ginac.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-10-15
 */
#ifndef __Ginac_H
#define __Ginac_H 1

#include <ginac/ginac.h>

namespace GiNaC
{
matrix grad( ex const& f, std::vector<symbol> const& l );
matrix grad( matrix const& f, std::vector<symbol> const& l );
matrix div( matrix const& f, std::vector<symbol> const& l );
ex laplacian( ex const& f, std::vector<symbol> const& l );
matrix laplacian( matrix const& f, std::vector<symbol> const& l );

} // GiNaC

namespace Feel
{
using  GiNaC::matrix;
using  GiNaC::symbol;
using  GiNaC::lst;
using  GiNaC::ex;
using  GiNaC::parser;

template<int Dim> inline std::vector<symbol> symbols() { return {symbol("x")}; }
template<> inline std::vector<symbol> symbols<1>() { return {symbol("x")}; }
template<> inline std::vector<symbol> symbols<2>() { return {symbol("x"),symbol("y") };}
template<> inline std::vector<symbol> symbols<3>() { return {symbol("x"),symbol("y"),symbol("z") };}

ex parse( std::string const& str, std::vector<symbol> const& syms );

namespace vf
{

/// \cond detail
/**
 * \class Ginac
 * \brief allow runtime ginac in expression
 *
 * @author Christophe Prud'homme
 * @see
 */
template<int Order = 2>
class GinacEx
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;
    static const uint16_type imorder = Order;
    static const bool imIsPoly = false;

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

    typedef GiNaC::ex expression_type;
    typedef GinacEx<Order> this_type;
    typedef double value_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit GinacEx( expression_type const & fun, std::vector<GiNaC::symbol> const& syms )
    :
    M_fun( fun ),
    M_syms( syms),
    M_cfun()
    {
        GiNaC::lst exprs(fun);
        GiNaC::lst syml;
        std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
        GiNaC::compile_ex(exprs, syml, M_cfun);

    }

   GinacEx( GinacEx const & fun )
        :
        M_fun( fun.M_fun ),
        M_syms( fun.M_syms),
        M_cfun()
    {
        GiNaC::lst exprs(M_fun);
        GiNaC::lst syml;
        std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );


        GiNaC::compile_ex(exprs, syml, M_cfun);
    }

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

    const GiNaC::FUNCP_CUBA& fun() const
    {
        return M_cfun;
    }

    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        // change 0 into rank
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<0>,mpl::int_<0> >,
                mpl::identity<Shape<gmc_type::nDim, Scalar, false, false> >,
                typename mpl::if_<mpl::equal_to<mpl::int_<0>,mpl::int_<1> >,
                mpl::identity<Shape<gmc_type::nDim, Vectorial, false, false> >,
                mpl::identity<Shape<gmc_type::nDim, Tensor2, false, false> > >::type >::type::type shape;

        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vec_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
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

            int no = 1;
            int ni = M_nsyms;///gmc_type::nDim;
            Eigen::VectorXd xi( M_nsyms );
            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    xi[k]=M_gmc->xReal( q )[k];
                for( int k = gmc_type::nDim; k < xi.size(); ++k )
                    xi[k] = 0;
                M_fun(&ni,xi.data(),&no,&M_y[q]);
            }

        }

        void update( Geo_t const& geom, uint16_type /*face*/ )
        {
            M_gmc =  fusion::at_key<key_type>( geom ).get();

            int no = 1;
            int ni = M_nsyms;//gmc_type::nDim;
            Eigen::VectorXd xi( M_nsyms );
            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    xi[k]=M_gmc->xReal( q )[k];
                for( int k = gmc_type::nDim; k < xi.size(); ++k )
                    xi[k] = 0;
                M_fun(&ni,xi.data(),&no,&M_y[q]);
            }
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }



        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q];
        }

        GiNaC::FUNCP_CUBA M_fun;
        gmc_ptrtype M_gmc;


        vec_type M_y;
        int M_nsyms;
    };

private:
    mutable expression_type  M_fun;
    std::vector<GiNaC::symbol> M_syms;
    GiNaC::FUNCP_CUBA M_cfun;
};

inline
Expr< GinacEx<2> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym )
{
    return Expr< GinacEx<2> >(  GinacEx<2>( f, lsym ) );
}

/**
 * \brief functor enabling ginac
 *
 */
template<int Order>
inline
Expr< GinacEx<Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym )
{
    return Expr< GinacEx<Order> >(  GinacEx<Order>( f, lsym ) );
}

template<int M=1, int N=1, int Order = 2>
class GinacMatrix
{
public:


    /** @name Typedefs
     */
    //@{


    static const size_type context = vm::POINT;
    static const bool is_terminal = false;
    static const uint16_type imorder = Order;
    static const bool imIsPoly = false;

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

    typedef GiNaC::ex expression_type;
    typedef GinacMatrix<M,N,Order> this_type;
    typedef double value_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    explicit GinacMatrix( GiNaC::matrix const & fun, std::vector<GiNaC::symbol> const& syms )
    :
        M_fun( fun.evalm() ),
        M_syms( syms),
        M_cfun()
        {
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
            GiNaC::compile_ex(exprs, syml, M_cfun);
        }
    explicit GinacMatrix( GiNaC::ex const & fun, std::vector<GiNaC::symbol> const& syms )
        :
        M_fun(fun.evalm()),
        M_syms( syms),
        M_cfun()
        {
            GiNaC::lst exprs;
            for( int i = 0; i < M_fun.nops(); ++i ) exprs.append( M_fun.op(i) );

            GiNaC::lst syml;
            std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
            GiNaC::compile_ex(exprs, syml, M_cfun);
        }

    GinacMatrix( GinacMatrix const & fun )
        :
        M_fun( fun.M_fun ),
        M_syms( fun.M_syms),
        M_cfun()
    {
        GiNaC::lst exprs;
        for( int i = 0; i < fun.M_fun.nops(); ++i ) exprs.append( fun.M_fun.op(i) );

        GiNaC::lst syml;
        std::for_each( M_syms.begin(),M_syms.end(), [&]( GiNaC::symbol const& s ) { syml.append(s); } );
        GiNaC::compile_ex(exprs, syml, M_cfun);
    }


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

    const GiNaC::FUNCP_CUBA& fun() const
    {
        return M_cfun;
    }

    std::vector<GiNaC::symbol> const& syms() const { return M_syms; }
    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        //typedef typename expression_type::value_type value_type;
        typedef double value_type;

        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t,vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

        typedef typename mn_to_shape<gmc_type::nDim,M,N>::type shape;
        typedef Eigen::Matrix<value_type,shape::M,shape::N,Eigen::RowMajor> vec_type;
        typedef std::vector<vec_type> loc_type;

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
        {}

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& /*fev*/ )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
        {}

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_fun( expr.fun() ),
            M_gmc( fusion::at_key<key_type>( geom ).get() ),
            M_y( M_gmc->nPoints() ),
            M_nsyms( expr.syms().size() )
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

            int no = M*N;
            int ni = M_nsyms;//gmc_type::nDim;
            Eigen::VectorXd xi( M_nsyms );
            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    xi[k]=M_gmc->xReal( q )[k];
                for( int k = gmc_type::nDim; k < xi.size(); ++k )
                    xi[k] = 0;
                M_fun(&ni,xi.data(),&no,M_y[q].data());
            }

        }

        void update( Geo_t const& geom, uint16_type /*face*/ )
        {
            M_gmc =  fusion::at_key<key_type>( geom ).get();

            int no = M*N;
            int ni = M_nsyms;//gmc_type::nDim;
            Eigen::VectorXd xi( M_nsyms );
            for(int q = 0; q < M_gmc->nPoints();++q )
            {
                for(int k = 0;k < gmc_type::nDim;++k )
                    xi[k]=M_gmc->xReal( q )[k];
                for( int k = gmc_type::nDim; k < xi.size(); ++k )
                    xi[k] = 0;
                M_fun(&ni,xi.data(),&no,M_y[q].data());
            }
        }


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q](c1,c2);
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q](c1,c2);
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_y[q](c1,c2);
        }

        GiNaC::FUNCP_CUBA M_fun;
        gmc_ptrtype M_gmc;
        loc_type M_y;
        int M_nsyms;
    };

private:
    mutable expression_type  M_fun;
    std::vector<GiNaC::symbol> M_syms;
    GiNaC::FUNCP_CUBA M_cfun;
}; // GinacMatrix
/// \endcond

inline
Expr< GinacMatrix<1,1,2> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym )
{
    return Expr< GinacMatrix<1,1,2> >(  GinacMatrix<1,1,2>( f, lsym ) );
}

/**
 * \brief functor enabling ginac
 *
 */
template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::matrix const& f, std::vector<GiNaC::symbol> const& lsym )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym ) );
}

template<int M, int N, int Order>
inline
Expr< GinacMatrix<M,N,Order> >
expr( GiNaC::ex const& f, std::vector<GiNaC::symbol> const& lsym )
{
    return Expr< GinacMatrix<M,N,Order> >(  GinacMatrix<M,N,Order>( f, lsym ) );
}


} // vf
} // Feel
#endif /* __Ginac_H */

