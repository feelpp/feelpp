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

#include "interpolator.hpp"

#ifndef FEELPP_FITDIFF_H
#define FEELPP_FITDIFF_H 1

namespace Feel
{
namespace vf
{
template<typename ExprT>
    class FitDiff
    {
    public:
    typedef ExprT expression_type;
    
    // list of mathematical objet that has to be precomputed
    static const size_type context = expression_type::context;
    // idv(T) is terminal
    // idv(T)+idv(T) is not
    // That is for optimisation
    // if not know: false
    static const bool is_terminal = expression_type::is_terminal;
    // im: integration method
    // l'ordre polynomial de la représentaiton
    static const uint16_type imorder = expression_type::imorder;
    // ...
    static const bool imIsPoly = expression_type::imIsPoly;
    
    using test_basis = typename expression_type::test_basis;
    using trial_basis = typename expression_type::trial_basis;
    
    typedef typename expression_type::value_type value_type;
    // scalar, vectorial (to be checked)
    typedef typename expression_type::evaluate_type evaluate_type;

    typedef FitDiff<ExprT> this_type;

    explicit FitDiff( expression_type const & __expr, std::string aDataFile = soption("fit.datafile"), int T = 3 ) : M_expr(__expr), dataFile(aDataFile), type(T)
    {
    }
    FitDiff(FitDiff const & te) : M_expr(te.M_expr), dataFile(te.dataFile), type(te.type)
        {
        }
    //~FitDiff(){}
    expression_type const& expression() const
    {
          return M_expr;
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
            M_tensor_expr( expr.expression(), geom, fev, feu)
        {
            M_interpolator = Interpolator::New(static_cast<InterpolationType>(expr.type), expr.dataFile);
        }

        // linear form
        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {
            M_interpolator = Interpolator::New(static_cast<InterpolationType>(expr.type), expr.dataFile);
        }

        // evaluation
        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
            M_interpolator = Interpolator::New(static_cast<InterpolationType>(expr.type), expr.dataFile);
        }

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

        // i : local index of basis test function
        // j : local index of basis trial function
        // c1 : id x, y  z of the first component (= 0 in scalar case)
        // c2 : id x, y  z of the second component (= 0 in scalar and vectorial case)
        // q : id of the current point in the geomap (basically, quadrature point or interpolation point)
        
        //local assembling - bilinear form a(i,j) 
        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 1.; //M_tensor_expr.evalij( i, j );
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
            // here : call the interpolator
            //return M_interpolator(M_tensor_expr.evalq(c1, c2, q));
            return M_interpolator->diff(M_tensor_expr.evalq(c1, c2, q));
        }

    private:
        tensor_expr_type M_tensor_expr;
        std::unique_ptr<Interpolator> M_interpolator;
    };
    /// end of tensor
        private:
            mutable expression_type  M_expr;
            std::string dataFile;
            int type;
    };

/**
 * \brief FitDiff
 **/
template<typename ExprT>
inline
Expr< FitDiff<ExprT> >
fitDiff( ExprT v,
        std::string dataFile = soption("fit.datafile"), 
        int intType = ioption("fit.kind") )
{
    typedef FitDiff<ExprT> fit_t;
    return Expr< fit_t >(  fit_t( v, dataFile, intType  ) );
}

}
}

#endif //FEELPP_FITDIFF_H
