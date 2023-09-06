/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Universite Joseph Fourier (Grenoble I)

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
   \file operators2.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#if !defined( __FEELPP_VF_OPERATORS2_HPP )
#define __FEELPP_VF_OPERATORS2_HPP 1

# include <boost/preprocessor/stringize.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
template <class Element1, class Element2>
class OpMass
{
public:

    static const size_type context = vm::MASS|vm::JACOBIAN;
    static const bool is_terminal = false;
    typedef Element1 test_element_type;
    typedef Element2 trial_element_type;
    typedef OpMass<test_element_type, trial_element_type> this_type;
    typedef this_type self_type;

    typedef typename test_element_type::polyset_type return_value_type;
    typedef strongest_numeric_type<typename test_element_type::value_type,
                                   typename trial_element_type::value_type> value_type;
    typedef typename test_element_type::functionspace_type test_functionspace_type; 
    typedef typename test_functionspace_type::reference_element_type test_fe_t;
    typedef typename trial_element_type::functionspace_type trial_functionspace_type; 
    typedef typename trial_functionspace_type::reference_element_type trial_fe_t; 
    
    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = true;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = true;
    };

    template<typename Func>
    static const bool has_test_basis = true;
    template<typename Func>
    static const bool has_trial_basis = true;
    using test_basis = test_fe_t;
    using trial_basis = trial_fe_t;
    //using value_type = typename test_functionspace_type::value_type;

    typedef ublas::matrix<value_type> matrix_type;

    OpMass ( test_element_type const& v,
             trial_element_type const& u )
        :
        M_v ( v ),
        M_u ( u ),
        M_exact_mass( M_v.functionSpace()->basis()->coeff() )
    {
        DVLOG(2) << "[" BOOST_PP_STRINGIZE( OpMass ) "] default constructorn";

        M_exact_mass = ublas::prod( test_element_type::polyset_type::toMatrix( M_v.functionSpace()->basis()->coeff() ),
                                    ublas::trans( trial_element_type::polyset_type::toMatrix( M_v.functionSpace()->basis()->coeff() ) ) );
        //Feel::cout << "m=" << M_exact_mass << std::endl;

#if 0
        ublas::scalar_vector<double> one (M_exact_mass.size2(),1);
        ublas::vector<double> w (M_exact_mass.size1());
        ublas::vector<double> wo (M_exact_mass.size1());
        wo = one;
        w = ublas::prod(M_exact_mass,wo);
        Feel::cout << "meas = " << ublas::inner_prod( w, one ) << std::endl;
#endif
    }
    OpMass( OpMass const& op )
        :
        M_v ( op.M_v ),
        M_u ( op.M_u ),
        M_exact_mass( op.M_exact_mass )
        //M_quad_mass() TO BE USED IF QUADRATURE IS NEEDED (transformation order >= 2)
    {
        DVLOG(2) << "[" BOOST_PP_STRINGIZE( OpMass ) "] copy constructorn";

    }

    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return Element1::functionspace_type::basis_type::nOrder+Element2::functionspace_type::basis_type::nOrder; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    test_element_type const& testFunction() const
    {
        return M_v;
    }
    trial_element_type const& trialFunction() const
    {
        return M_u;
    }

    value_type exactMass( uint16_type i, uint16_type j ) const
    {
        return M_exact_mass( i, j );
    }
    matrix_type const& exactMass() const
    {
        return M_exact_mass;
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef this_type expression_type;
        typedef Basis_i_t test_basis_context_type;
        typedef Basis_j_t trial_basis_context_type;

        using key_type = key_t<Geo_t>;
        using gmc_ptrtype = gmc_ptr_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        typedef Shape<gmc_type::nDim, Scalar, false, false> shape;
        using value_type = typename test_element_type::value_type;
        //typedef typename test_basis_context_type::value_type value_type;
        //typedef typename test_basis_context_type::polyset_type return_value_type;
        template<typename Func>
        struct HasTestFunction
        {
            static const bool result = true;
        };
        template<typename Func>
        struct HasTrialFunction
        {
            static const bool result = true;
        };

        template<typename Func>
        static const bool has_test_basis = true;
        template<typename Func>
        static const bool has_trial_basis = true;

        //static inline const uint16_type nComponents = return_value_type::nComponents;
        using test_basis = test_fe_t;
        using trial_basis = trial_fe_t;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            typedef value_type type;
        };

        struct is_zero
        {
            static const bool value = false;
        };
        
        tensor( this_type const& expr,
                Geo_t const& /*geom*/,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_fev( fev ),
            M_feu( feu ),
            M_mat( expr.exactMass() )
        {}

        void update( Geo_t const& geom, Basis_i_t const& , Basis_j_t const&  )
            {
                M_gmc = fusion::at_key<key_type>( geom ).get();
            }
        void update( Geo_t const& , Basis_i_t const&  )
            {
            }
        void update( Geo_t const& )
            {
            }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
            {
            }
        value_type
        operator()( uint16_type i, uint16_type j ) const
        {
            return M_mat( i, j )*M_gmc->J( 0 );
        }
        value_type
        evalij( uint16_type i, uint16_type j ) const
            {
                return M_mat( i, j )*M_gmc->J( 0 );
            }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
            {
                return M_mat( i, j )*M_gmc->J( 0 );
            }
        value_type
        evalijq( uint16_type i, uint16_type j, int c1, int c2, uint16_type q ) const
            {
                return M_mat( i, j )*M_gmc->J( 0 );
            }


        value_type
        operator()( uint16_type i, uint16_type j, int q ) const
        {
            return M_mat( i, j );
            //return 0;//M_expr.quadratureMass( q, i, j );
        }

        test_basis_context_type const& M_fev;
        trial_basis_context_type const& M_feu;
        //this_type const& M_expr;
        matrix_type const& M_mat;
        gmc_ptrtype M_gmc;
    };

protected:
    OpMass () {}

    test_element_type const& M_v;
    trial_element_type const& M_u;
    ublas::matrix<value_type> M_exact_mass;
    
};
/// \endcond
/**
 * \brief mass term
 */
template <class Element1, class Element2>
inline Expr< OpMass< Element1, Element2> >
mass( Element1 const& el1, Element2 const& el2 )
{
    typedef OpMass< Element1, Element2> expr_t;
    return Expr< expr_t >(  expr_t( el1, el2 ) );
}

} // vf
} // feel

#endif /* __FEELPP_VF_OPERATORS2_HPP */
