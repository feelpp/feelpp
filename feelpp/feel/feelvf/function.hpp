/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2007-2012 Universite Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2007-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
    static inline const bool is_terminal = false;

    template<typename Funct>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Funct>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Funct>
    static inline const bool has_test_basis = false;
    template<typename Funct>
    static inline const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    typedef Func expression_type;
    typedef Function<Func> this_type;
    typedef typename expression_type::value_type value_type;

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
            static inline const bool value = false;
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


        value_type
        evalij( uint16_type i, uint16_type j ) const
        {
            return 0;
        }


        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evaluate( c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return this->evaluate( c1, c2, q );
        }



        value_type
        evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evaluate( c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evaluate( c1, c2, q );
        }

    private:
        value_type
        evaluate( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if constexpr ( expression_type::context & vm::NORMAL )
                return M_fun( c1, c2, M_gmc->xReal( q ),
                              M_gmc->unitNormal( q ) );
            else
            {
                // The legacy functor signature always accepts a normal, even
                // when its context declares that the argument is unused.
                using normal_type =
                    eigen_vector_type<gmc_type::nDim, value_type>;
                static const normal_type unusedNormal = normal_type::Zero();
                return M_fun( c1, c2, M_gmc->xReal( q ), unusedNormal );
            }
        }

    public:
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
