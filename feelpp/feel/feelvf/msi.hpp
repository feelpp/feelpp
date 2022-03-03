/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 May 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_VF_MSI_HPP
#define FEELPP_VF_MSI_HPP 1

#include <feel/feelconfig.h>


#if defined( FEELPP_HAS_FFTW )
#include <boost/multi_array.hpp>
#include <feel/feeldiscr/multiscaleimage.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/shape.hpp>
#include <unsupported/Eigen/Splines>

namespace Eigen {
  
// lets do some explicit instantiations and thus
// force the compilation of all spline functions...
template class Spline<double, 2, Dynamic>;
template class Spline<double, 3, Dynamic>;

template class Spline<double, 2, 2>;
template class Spline<double, 2, 3>;
}

namespace Feel
{
namespace vf
{
/// \cond DETAIL
namespace detail
{
/**
 * @brief Return multiscale image accessor from coarse to fine grid
 * \code
 * // coarse level with respect to fine level
 * int L = 4;
 * auto l = form1(_test=Xh);
 * l=integrate( _range=elements(coarse_mesh), _expr=msi(f)*id(v), _quad=im_msi(L) );
 * \endcode
 * @author Christophe Prud'homme
 */
template<typename T, int _Options = 0>
class MSI
{
public:


    /** @name Typedefs
     */
    //@{
    static const size_type context = vm::POINT;

    static const int Options = _Options;

    static const bool is_terminal = true;

    typedef Feel::MultiScaleImage<T,Options> msi_type;
    using needs_gradient_t = typename msi_type::needs_gradient_t;
    using do_compute_gradient_t = typename msi_type::do_compute_gradient_t;
    using no_compute_gradient_t = typename msi_type::no_compute_gradient_t;

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

    template<typename Func>
    static const bool has_test_basis = false;
    template<typename Func>
    static const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;


    typedef MSI<T,Options> this_type;
    typedef T value_type;

    using image_type = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MSI( image_type const& m, int L )
        :
        M_coarse2fine( m, L )
    {

    }

    MSI( MSI const &  ) = default;

    ~MSI() = default;

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
    uint16_type polynomialOrder() const { return 0; }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    //blitz::Array<value_type,2> ones() const { return M_values; }

    value_type coarse2fine( ublas::vector<T> const& real, ublas::vector<T> const& ref   ) const
    {
        return M_coarse2fine( real, ref );
    }
    value_type coarse2fine( int c, ublas::vector<T> const& real, ublas::vector<T> const& ref   ) const
        {
            return M_coarse2fine( c, real, ref );
        }


private:
    
    msi_type M_coarse2fine;
    
public:
    //@}
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef typename expression_type::value_type value_type;
        typedef value_type return_value_type;
        using key_type = key_t<Geo_t>;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type* gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef typename mpl::if_<needs_gradient_t,
                                  mpl::identity<Shape<gmc_type::nDim, Vectorial, true, false>>,
                                  mpl::identity<Shape<gmc_type::nDim, Scalar, false, false>> >::type::type shape;

        template <class Args> struct sig
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
            :
            M_expr( expr ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
            //std::cout << "tensor::ones = " << M_expr.ones() << "\n";
        }

        tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const& )
            :
            M_expr( expr ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        tensor( this_type const& expr, Geo_t const&  geom )
            :
            M_expr( expr ),
            M_gmc( fusion::at_key<key_type>( geom ).get() )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const&, Basis_j_t const& )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& )
        {
            update( geom );
            
        }

        void update( Geo_t const& geom)
        {
            M_gmc = fusion::at_key<key_type>( geom ).get();
        }


        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, q );
        }

        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( j );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( q );
            return eval( c1, c2, q );
        }
    private:
        
        value_type
        eval( int c1, int c2, uint16_type q ) const
        {
            return eval( c1, c2, q, needs_gradient_t() );
        }
        /**
         * compute the gradient and return the component \c c at point \c q
         */
        value_type
        eval( int c1, int c2, uint16_type q, do_compute_gradient_t ) const
            {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            
            return M_expr.coarse2fine( c2, M_gmc->xReal(q), M_gmc->xRef(q) );
            }
        /**
         * compute the value of the image of point \c q
         */
        value_type
        eval( int c1, int c2, uint16_type q, no_compute_gradient_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
            
                return M_expr.coarse2fine( M_gmc->xReal(q), M_gmc->xRef(q) );
            }
        this_type M_expr;
        gmc_ptrtype M_gmc;
        
    };

};
} // detail
/// \endcond

/**
 *
 */
template<typename T, int Options = 0>
inline
Expr<vf::detail::MSI<T, Options> >
msi( Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> const& f, int level )
{
    return Expr<vf::detail::MSI<T,Options> >( vf::detail::MSI<T,Options>(f,level ));
}

} // vf
} // Feel

#endif // FEELPP_HAS_FFTW

#endif /* __MSI_H */
