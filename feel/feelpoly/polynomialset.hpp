/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Gilles Steiner <gilles.steiner@epfl.ch>
       Date: 2005-10-11

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Université Joseph Fourier Grenoble 1

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
   \file polynomialset.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-11
 */
#ifndef __PolynomialSet_H
#define __PolynomialSet_H 1

#include <vector>
#include <boost/multi_array.hpp>
#include <boost/multi_array/extent_gen.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/mpl/min_max.hpp>
#include <Eigen/Core>

#include <feel/feelpoly/context.hpp>

#include <feel/feelalg/svd.hpp>
#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/legendre.hpp>
#include <feel/feelpoly/boundadapted.hpp>
#include <feel/feelpoly/tensorisedboundadapted.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feelpoly/quadmapped.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>

namespace Feel
{

/**
 * \class PolynomialSet
 *  \brief a Set of polynomials
 *
 * This class represents a set of polynomials \f$ {p_i}_{i=1...N} \f$
 * defined in a certain basis given by the template argument. The
 * coefficients of the polynomials in the basis are represented by a
 * matrix whose line represents the polynomials and columns the basis
 * index \f$ C_{i,j} = \mathcal{R}( p_i)_j \f$ where \f$\mathcal{R}\f$
 * is the mapping between the polynomial and its coefficients.
 * We have that the polynomial set is represented as follows:
 * \f[p_i = \sum_j=1^N \mathcal{R}( p_i )_j \phi_j\f]
 *
 *  \ingroup Polynomial
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename Poly, template<uint16_type> class PolySetType = Scalar >
class PolynomialSet
    :
        public boost::enable_shared_from_this<PolynomialSet<Poly, PolySetType > >
{
    typedef boost::enable_shared_from_this<PolynomialSet<Poly, PolySetType > > super_enable_shared_from_this;
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Poly::nDim;
    static const uint16_type nRealDim = Poly::nRealDim;
    static const uint16_type nOrder = Poly::nOrder;

    //@}

    /** @name Typedefs
     */
    //@{
    typedef PolynomialSet<Poly, PolySetType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;
    static const bool is_product = Poly::is_product;


    typedef PolySetType<nRealDim> polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static const uint16_type rank = polyset_type::rank;

    typedef PolynomialSet<Poly, Scalar> component_type;
    typedef Polynomial<Poly, PolySetType> polynomial_type;
    //typedef Polynomial<Poly, PolySetType, ublas::vector_range<ublas::vector<value_type> > > polynomial_type;
    typedef polynomial_type polynomial_view_type;

    typedef typename Poly::convex_type convex_type;
    typedef typename basis_type::matrix_type matrix_type;
    typedef typename basis_type::points_type points_type;


    typedef PolynomialSet<Poly, Vectorial> gradient_polynomialset_type;

    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, value_type>::value ) );
    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, typename points_type::value_type>::value ) );

    //@}

    /** @name Constructors, destructor
     */
    //@{

    PolynomialSet()
        :
        M_basis(),
        M_coeff(),
        M_fname( "pset" )
    {
        //std::cout << "[PolynomialSet::default] dim = " << nDim << " order = " << nOrder << "\n";
    }
    PolynomialSet( Poly const& p )
        :
        M_basis( p.basis() ),
        M_coeff( p.coeff() ),
        M_fname( p.familyName() )
    {}
    /**
     */
    //template<typename AE>
    //PolynomialSet( Poly const& p, ublas::matrix_expression<AE> const& c )
    PolynomialSet( Poly const& p, matrix_type const& c, bool __as_is = false )
        :
        M_basis( p.basis() ),
        M_coeff( p.coeff() ),
        M_fname( p.familyName() )
    {
        setCoefficient( c, __as_is );
        //FEELPP_ASSERT( c.size2() == p.coeff().size1() )( c.size2() )( p.coeff().size1() ).error( "invalid dimension\n" );
        //std::cout << "[PolynomialSet] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet] c = " << c << "\n";
        //std::cout << "[PolynomialSet] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet] coeff = " << M_coeff << "\n";
    }

    /**
     */
    //template<typename AE>
    //PolynomialSet( Poly const& p, ublas::matrix_expression<AE> const& c )
    PolynomialSet( matrix_type const& c, bool __as_is = false )
        :
        M_basis(),
        M_coeff( c ),
        M_fname( "pset" )
    {
        setCoefficient( c, __as_is );
        //FEELPP_ASSERT( c.size2() == p.coeff().size1() )( c.size2() )( p.coeff().size1() ).error( "invalid dimension\n" );
        //std::cout << "[PolynomialSet] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet] c = " << c << "\n";
        //std::cout << "[PolynomialSet] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet] coeff = " << M_coeff << "\n";
    }


    PolynomialSet( PolynomialSet const & p )
        :
        M_basis( p.M_basis ),
        M_coeff( p.M_coeff ),
        M_fname( p.M_fname )
    {
        //std::cout << "[PolynomialSet::copy] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet::copy] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet::copy] coeff = " << M_coeff << "\n";
    }

    virtual ~PolynomialSet()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& pset )
    {
        if ( this != &pset )
        {
            M_basis = pset.M_basis;
            M_coeff = pset.M_coeff;
        }

        return *this;
    }

    /**
     * \brief extract the i-th component of a vectorial polynomial set
     *
     * \return the i-th component of the polynomial set
     */
    component_type operator[]( uint16_type i ) const
    {
        BOOST_STATIC_ASSERT( is_vectorial );
        FEELPP_ASSERT( i < nComponents )( i )( nComponents ).error ( "invalid component index" );
        const int nrows = M_coeff.size1()/nComponents;
        const int ncols = M_coeff.size2();
        return component_type( Poly(), ublas::project( M_coeff,
                               ublas::slice( nrows*i+i, nComponents, nrows/nComponents ),
                               ublas::slice( 0, 1, ncols ) ), true );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the degree of the polynomials in the set
     */
    uint16_type degree() const
    {
        return nOrder;
    }

    /**
     * \return coefficient of the polynomials of the polynomial set in
     * the basis associated with Poly
     */
    matrix_type const&  coeff() const
    {
        return  M_coeff;
    }

    /**
     * \return the basis associated with the polynomial set
     */
    basis_type const& basis() const
    {
        return M_basis;
    }

    /**
     * \return true if the polynomial set is scalar, false otherwise
     */
    static bool isScalar()
    {
        return is_scalar;
    }

    /**
     * \return true if the polynomial set is vectorial, false otherwise
     */
    static bool isVectorial()
    {
        return is_vectorial;
    }

    /**
     *
     */
    static uint16_type numberOfComponents()
    {
        return nComponents;
    }

    /**
     * \return the polynomial dimension
     */
    size_type polynomialDimension() const
    {
        return M_coeff.size1()/nComponents;
    }

    /**
     * \return the polynomial dimension per component
     */
    size_type polynomialDimensionPerComponent() const
    {
        return M_coeff.size2();
    }

    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
        return ublas::norm_frobenius( M_coeff ) < Feel::type_traits<value_type>::epsilon();
    }

    /**
     * the \c familyName() identifies the finite element
     * \return the family name of a finite element
     */
    virtual std::string familyName() const
    {
        return M_fname;
    }

    /**
     * the name of a finite element is a string composed by:
     *
     * -# a prefix which identifies the family (eg lagrange)
     * -# the dimension
     * -# the order of the finite element
     *
     * \param sep separator between family name, dimension and order (by default it is '.')
     *
     * \return the name of the finite element
     */
    std::string name( std::string sep = "." ) const
    {
        std::ostringstream os;
        os << this->familyName() << sep << nDim << sep << nOrder;
        return os.str();
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the coefficient of the polynomial set in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {
        if ( is_scalar )
        {
            if ( !__as_is )
                M_coeff = ublas::prod( __c, M_coeff );

            else
                M_coeff = __c;
        }

        else
        {
            if ( !__as_is )
            {
                M_coeff = ublas::prod( __c, polyset_type::toMatrix( M_coeff ) );
                M_coeff = polyset_type::toType( M_coeff );
            }

            else
                M_coeff = __c;
        }
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Extract the polynomials whose indices are listed in \p list_p
     *
     * \param list_p list of indices of polynomials to extract
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly, PolySetType> polynomials( std::vector<uint16_type> const& list_p  ) const
    {
        size_type dim_p = this->polynomialDimension();
        size_type new_dim_p = nComponents*list_p.size();
        matrix_type coeff( nComponents*nComponents*list_p.size(), M_coeff.size2() );
        int j = 0;
        BOOST_FOREACH( uint16_type i, list_p )
        {
            for ( int c = 0; c < nComponents; ++c )
            {
                ublas::project( coeff,
                                ublas::range( c*new_dim_p+nComponents*j, c*dim_p+nComponents*j+nComponents ),
                                ublas::range( 0, M_coeff.size2() ) ) =
                                    ublas::project( M_coeff,
                                                    ublas::range( c*dim_p+nComponents*i, c*dim_p+nComponents*i+nComponents ),
                                                    ublas::range( 0, M_coeff.size2() ) );
            }

            ++j;
        }
        return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials up to dimension \p dim
     *
     * \param dim_p polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly, PolySetType> polynomialsUpToDimension( int dim_p  ) const
    {
        matrix_type coeff( nComponents*nComponents*dim_p, M_coeff.size2() );

        for ( int c = 0; c < nComponents; ++c )
        {
            size_type nc = c*this->polynomialDimension();
            ublas::project( coeff,
                            ublas::range( c*nComponents*dim_p, ( c+1 )*nComponents*dim_p ),
                            ublas::range( 0, M_coeff.size2() ) ) =
                                ublas::project( M_coeff,
                                                ublas::range( nc, nc+nComponents*dim_p ),
                                                ublas::range( 0, M_coeff.size2() ) );
        }

        return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the polynomials between dimension \p dim_bot up to dimension \p dim_top
     *
     * \param dim_bot polynomial dimension
     * \param dim_top polynomial dimension
     * \return the polynomial set extracted
     */
    PolynomialSet<Poly, PolySetType>
    polynomialsRange( uint16_type dim_bot, uint16_type dim_top  ) const
    {
        uint16_type dim_p = dim_top-dim_bot;
        matrix_type coeff( nComponents*nComponents*dim_p,
                           M_coeff.size2() );

        for ( int c = 0; c < nComponents; ++c )
        {
            size_type nc = c*this->polynomialDimension();
            ublas::project( coeff,
                            ublas::range( c*nComponents*dim_bot, ( c+1 )*nComponents*dim_top ),
                            ublas::range( 0, M_coeff.size2() ) ) =
                                ublas::project( M_coeff,
                                                ublas::range( nc+dim_bot, nc+nComponents*dim_top ),
                                                ublas::range( 0, M_coeff.size2() ) );
        }

        return PolynomialSet<Poly, PolySetType>( Poly(), coeff, true );
    }

    /**
     * Extract the \p i -th polynomial
     *
     * \param i index of the polynomial to extract
     * \return the polynomial extracted
     */
    Polynomial<Poly, PolySetType> polynomial( uint16_type i  ) const
    {
        size_type dim_p = this->polynomialDimension();
        matrix_type coeff( nComponents, M_coeff.size2() );
        //for ( int c = 0; c < nComponents; ++c )
        {
            ublas::project( coeff,
                            ublas::range( 0, nComponents ),
                            ublas::range( 0, M_coeff.size2() ) ) =
                                ublas::project( M_coeff,
                                                //ublas::range( c*dim_p+nComponents*i, c*dim_p+nComponents*i+nComponents ),
                                                ublas::range( nComponents*i, nComponents*( i+1 ) ),
                                                ublas::range( 0, M_coeff.size2() ) );
        }
        return Polynomial<Poly, PolySetType> ( Poly(), coeff, true );
    }


    /**
     * evaluate the i-th polynomial at node __pt
     *
     * \warning this function is not efficient at all, the preferred
     * method is to evaluate at a set of points
     */
    template<typename AE>
    ublas::vector<value_type> evaluate( uint16_type i, ublas::vector_expression<AE> const& __pt ) const
    {
        return ublas::row( ublas::prod( M_coeff, M_basis( __pt ) ),  i );
    }

    /**
     * evaluate all polynomials at node __pt
     * \return a column matrix
     */
    template<typename AE>
    ublas::matrix<value_type> evaluate( ublas::vector_expression<AE> const& __pt ) const
    {
        return ublas::prod( M_coeff, M_basis( __pt ) );
    }


    /**
     * evaluate all polynomials of the set at a set of nodes
     *
     * Constructs \f$A_{i,j} = p_i(x_j) = \sum_{k=1}^N \mathcal{R}(p_i)_k \phi_k(x_j)\f$
     *
     * \arg __pts a column oriented matrix contained the node
     * coordinates (in the columns).
     *
     * \return the matrix \f$A\f$
     */
    template<typename AE>
    matrix_type evaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        matrix_type m ( M_basis.evaluate( __pts ) );
        FEELPP_ASSERT( M_coeff.size2() == m.size1() )( M_coeff.size2() )( m.size1() ).error( "invalid size" );
        return ublas::prod( M_coeff, m );
    }

    /**
     * Derivate with respect to the \f$\ell\f$ coordinates at the
     * nodes where the polynomials basis have been constructed.
     *
     * We construct the matrix \f{eqnarray*} A_{i,j} &= \frac{\partial
     * p_i(x_j)}{\partial x_\ell}\\ &= \sum_{k=1}^N \mathcal{R}(p_i)_k
     * \frac{\partial \phi_k(x_j)}{\partial x_\ell} \f}
     *
     * \arg l the derivation index \f$\ell\f$
     * \arg __pts a column oriented matrix contained the node
     * coordinates (in the columns).
     *
     * \return the matrix \f$A\f$
     */

    /**
     * \brief derivatives of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
        return M_basis.d( i );
    }

    matrix_type d( uint16_type i, uint16_type j ) const
    {
        return ublas::prod( M_basis.d( i ), M_basis.d( j ) );
    }

    /**
     * \brief Derivate with respect to the l-th direction.
     *
     * \return the polynomial set associated with the derivation in the l-direction
     */
    self_type derivate( uint16_type l ) const
    {
        return self_type( Poly(), ublas::prod(  M_coeff, M_basis.d( l ) ), true );
    }

    template<typename AE>
    ublas::vector<matrix_type> derivate( ublas::matrix_expression<AE> const& pts ) const
    {
        ublas::vector<matrix_type> der( M_basis.derivate( pts ) );
        ublas::vector<matrix_type> res( nDim );

        for ( uint16_type i = 0; i < nDim; ++i )
        {
            res[i].resize( M_coeff.size1(), pts().size2() );
            ublas::axpy_prod( M_coeff, der[i], res[i] );
        }

        return res;
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, ublas::matrix_expression<AE> const& pts ) const
    {
        ublas::vector<matrix_type> der( M_basis.derivate( pts ) );
        matrix_type res( M_coeff.size1(), pts().size2() );
        ublas::axpy_prod( M_coeff, der[i], res );
        return res;
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, uint16_type j, ublas::matrix_expression<AE> const& pts ) const
    {
        //std::cout << "[derivate2] M_coeff = " << M_coeff << "\n";
        matrix_type eval( M_basis.evaluate( pts ) );
        //matrix_type res( M_coeff.size1(), pts().size2() );
        //ublas::axpy_prod( M_coeff, der[i], res );
        matrix_type p1 = ublas::prod( M_coeff, M_basis.d( i ) );
        matrix_type p2 = ublas::prod( p1, M_basis.d( j ) );
        return ublas::prod( p2, eval );
    }
    /**
     * Gradient of the polynomial set
     *
     * Computes the gradient of the polynomial set.
     */
    gradient_polynomialset_type
    gradient() const
    {
        return gradient( mpl::int_<polyset_type::rank>() );
    }

    gradient_polynomialset_type
    gradient( mpl::int_<0> ) const
    {
        const int n1 = M_coeff.size1();
        const int n2 = M_coeff.size2();
        ublas::matrix<value_type> c ( nDim*nDim*n1, n2 );
        c.clear();

        for ( int i = 0; i <nDim; ++i )
        {
            ublas::project( c,
                            ublas::slice( nDim*n1*i+i, nDim, n1 ),
                            ublas::slice( 0, 1, n2 ) )  = ublas::prod( M_coeff, M_basis.d( i ) );
        }

        return gradient_polynomialset_type( c, true );
    }

    gradient_polynomialset_type
    gradient( mpl::int_<1> ) const
    {
        // we deal with a VECTORIAL polynomial set: nComponents = nDim*nDim
        // number of basis functions times the numb er of components(was vertorial function)
        const int n1 = M_coeff.size1();
        const int n2 = M_coeff.size2();
        ublas::matrix<value_type> c ( nComponents*nComponents*n1, n2 );
        c.clear();
#if 0
        std::cout << "M_coeff = " << M_coeff << "\n"
                  << "c = " << c << "\n";
#endif

        for ( int i = 0; i <nDim; ++i )
        {
#if 0
            std::cout << "i=" << i << "\n"
                      << " prod = "
                      << ublas::prod( M_coeff, M_basis.d( i ) ) << "\n"
                      << " slice = "
                      << ublas::project( c,
                                         ublas::slice( nDim*n1*i+i, nDim, n1 ),
                                         ublas::slice( 0, 1, n2 ) ) << "\n";
            // c[i][c1][c2][q]
#endif
            ublas::project( c,
                            ublas::slice( nDim*n1*i+i, nDim, n1 ),
                            ublas::slice( 0, 1, n2 ) )  = ublas::prod( M_coeff, M_basis.d( i ) );
        }

        //std::cout << "coeff = " << c << "\n";
        return gradient_polynomialset_type( c, true );
    }




    /**
     * \return the number of degrees of freedom
     */
    uint16_type nbDof() const
    {
        return M_coeff.size1();
    }


    /**
     * insert the polynomial set \p p at the end of the set
     */
    void insert( PolynomialSet<Poly,PolySetType> const& p, bool erase = false )
    {
        FEELPP_ASSERT( p.coeff().size2() == coeff().size2() )( p.coeff().size2() )( coeff().size2() ).warn( "invalid polynomial set" );

#if 0
        std::cout << "coeff = " << M_coeff << "\n"
                  << "p     = " << p.coeff() << "\n";
#endif

        if ( erase )
        {
            M_coeff = p.M_coeff;
            //std::cout << "insert after erase ="  << M_coeff << "\n";
            return;
        }

        matrix_type oldcoeff = M_coeff;
        M_coeff.resize( M_coeff.size1() + p.coeff().size1(), p.coeff().size2(), false );


        if ( oldcoeff.size1() > 0 )
        {
            ublas::project( M_coeff,
                            ublas::range( 0, oldcoeff.size1() ),
                            ublas::range( 0, oldcoeff.size2() ) ) = oldcoeff;
        }

        ublas::project( M_coeff,
                        ublas::range( oldcoeff.size1(), oldcoeff.size1()+p.coeff().size1() ),
                        ublas::range( 0, oldcoeff.size2() ) ) = p.coeff();
        //std::cout << "insert ="  << M_coeff << "\n";

    }


    //@}


    /******************************************************/
    /**
     * \class PreCompute
     */
    class PreCompute
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef PolynomialSet<Poly, PolySetType > reference_element_type;
        typedef boost::shared_ptr<reference_element_type> reference_element_ptrtype;

        typedef typename reference_element_type::value_type value_type;

        static const uint16_type nDim = reference_element_type::nDim;
        static const uint16_type nComponents1 = reference_element_type::nComponents1;
        static const uint16_type nComponents2 = reference_element_type::nComponents2;
        //static const uint16_type nComponents3 = reference_element_type::nComponents3;
        static const uint16_type nComponents = reference_element_type::nComponents;

        typedef typename reference_element_type::points_type matrix_node_t_type;
        typedef ublas::matrix<value_type> matrix_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nDim> g_type;
        typedef Eigen::Matrix<value_type,nDim,nDim> h_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> l_type;
        typedef boost::multi_array<id_type,2> functionvalue_type;
        typedef boost::multi_array<g_type,2> grad_type;
        typedef boost::multi_array<h_type,2> hessian_type;
        typedef boost::multi_array<l_type,2> laplacian_type;

        PreCompute() {}

#if 0
        /**
         * precompute values of basis functions given by __ref_ele at
         * points \c __pts
         */
        PreCompute( matrix_node_t_type const& __pts )
            :
            M_ref_ele( new reference_element_type() ),
            M_nodes( __pts ),
            M_phi(),
            M_grad(),
            M_hessian(),
            M_laplacian()
        {
            init( M_ref_ele, __pts, mpl::int_<rank>() );
        }
#endif

        /**
         * precompute values of basis functions given by __ref_ele at
         * points \c __pts
         */
        PreCompute( reference_element_ptrtype const& __ref_ele,
                    matrix_node_t_type const& __pts )
            :
            M_ref_ele( __ref_ele ),
            M_nodes( __pts ),
            M_phi(),
            M_grad(),
            M_hessian(),
            M_laplacian()
        {
            init( M_ref_ele, __pts, mpl::int_<rank>() );
        }

        /** copy constructor (deep copy) */
        PreCompute( PreCompute const& __pc )
            :
            M_ref_ele( __pc.M_ref_ele ),
            M_nodes( __pc.M_nodes ),
            M_phi( __pc.M_phi ),
            M_grad( __pc.M_grad ),
            M_hessian( __pc.M_hessian ),
            M_laplacian( __pc.M_laplacian )
        {}

        /** */
        ~PreCompute()
        {}

        /** copy operator (deep copy) */
        PreCompute& operator=( PreCompute const& __pc )
        {
            if ( this != &__pc )
            {
                M_ref_ele = __pc.M_ref_ele;
                M_nodes = __pc.M_nodes;
                M_phi = __pc.M_phi;
                M_grad = __pc.M_grad;
                M_hessian = __pc.M_hessian;
                M_laplacian = __pc.M_laplacian;
            }

            return *this;
        }

        void update( matrix_node_t_type const& __pts )
        {
            M_nodes = __pts;
            init( M_ref_ele, __pts, mpl::int_<rank>() );
        }

        //! \return the finite element
        reference_element_ptrtype fePtr() { return M_ref_ele; }

        /**
           \return the dimension of the space where the nodes are defined
        */
        uint16_type dim() const
        {
            return reference_element_type::nDim;
        }

        /**
           \return the number of nodes in the reference element where
           the preconputation is done
        */
        uint16_type nComputedNodes() const
        {
            return M_nodes.size2();
        }

        /**
           \return the number of nodes in the reference element where
           the preconputation is done
        */
        uint16_type nPoints() const
        {
            return M_nodes.size2();
        }

        /**
           \return the nodes at which the basis functions, gradient
           and hessian have been computed in the reference element
        */
        matrix_node_t_type const& nodes() const
        {
            return M_nodes;
        }
        matrix_node_t_type const& points() const
        {
            return M_nodes;
        }

        /**
           \return the nodes at which the basis functions, gradient
           and hessian have been computed in the reference element
        */
        ublas::matrix_column<matrix_node_t_type const> node( uint16_type __i ) const
        {
            return ublas::column( M_nodes, __i );
        }

        /**
         *  Return the matrix evaluation the basis functions (rows) at
         *  a set of points (columns). The matrix is column oriented,
         *  so performance wise it is better to iterate over the columns
         *
         */
        functionvalue_type const& phi() const
        {
            return M_phi;
        }

        functionvalue_type const*  phiPtr() const
            {
                return &M_phi;
            }
        /**
         * Returns the value of the q-th node of the i-th basis
         * functions.
         */
        value_type phi( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            CHECK( q < nComputedNodes() ) << "invalid node index q=" << q << " nComputedNodes="<< nComputedNodes();
            CHECK( i < M_ref_ele->nbDof() ) << "invalid dof index i=" << i << " ndof =" << M_ref_ele->nbDof();

            return M_phi[i][q]( c1,c2 );
        }
        id_type const& phi( uint16_type i, uint16_type q ) const
        {
            CHECK( q < nComputedNodes() ) << "invalid node index q=" << q << " nComputedNodes="<< nComputedNodes();
            CHECK( i < M_ref_ele->nbDof() ) << "invalid dof index i=" << i << " ndof =" << M_ref_ele->nbDof();

            return M_phi[i][q];
        }
        id_type& phi( uint16_type i, uint16_type q )
        {
            CHECK( q < nComputedNodes() ) << "invalid node index q=" << q << " nComputedNodes="<< nComputedNodes();
            CHECK( i < M_ref_ele->nbDof() ) << "invalid dof index i=" << i << " ndof =" << M_ref_ele->nbDof();

            return M_phi[i][q];
        }

        grad_type const& grad() const
        {
            return M_grad;
        }
        grad_type const* gradPtr() const
            {
                return &M_grad;
            }

        value_type grad( size_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_grad[i][q]( c1,c2 );
        }

        hessian_type const& hessian() const
        {
            return M_hessian;
        }

        value_type hessian( size_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_hessian[i][q]( c1,c2 );
        }



        laplacian_type const& laplacian() const
            {
                return M_laplacian;
            }

        value_type laplacian( size_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
            {
                return M_laplacian[i][q]( c1,c2 );
            }

        void print()
            {
                typedef typename grad_type::index index;
                const index I = M_ref_ele->nbDof();
                const index Q = M_nodes.size2();

                LOG(INFO) << "Precompute nodes : " << M_nodes;
                for ( index i = 0; i < I; ++i )
                {
                    for ( index q = 0; q < Q; ++q )
                        LOG(INFO) << "Precompute phi[" << i << ", " << q << "]= " << M_phi[i][q](0,0);
                }
            }
    private:

        void
        init( reference_element_ptrtype const& M_ref_ele,
              matrix_node_t_type const& __pts,
              mpl::int_<0> )
        {
            M_phi.resize( boost::extents[M_ref_ele->nbDof()][__pts.size2()] );
            M_grad.resize( boost::extents[M_ref_ele->nbDof()][__pts.size2()] );
            M_hessian.resize( boost::extents[M_ref_ele->nbDof()][__pts.size2()] );
            M_laplacian.resize( boost::extents[M_ref_ele->nbDof()][__pts.size2()] );

            matrix_type phiv = M_ref_ele->evaluate( __pts );
            ublas::vector<matrix_type> __grad( M_ref_ele->derivate( __pts ) );
            matrix_type __hessian( M_ref_ele->gradient().gradient().evaluate( __pts ) );

            typedef typename grad_type::index index;
            const index I = M_ref_ele->nbDof();
            const index Q = __pts.size2();

            for ( index i = 0; i < I; ++i )
            {
                for ( index q = 0; q < Q; ++q )
                    M_phi[i][q]( 0,0 ) = phiv( i, q );

                for ( index q = 0; q < Q; ++q )
                    for ( index j = 0; j < nDim; ++j )
                        M_grad[i][q]( 0,j ) = __grad[j]( i, q );

                for ( index q = 0; q < Q; ++q )
                    for ( index j = 0; j < nDim; ++j )
                        for ( index k = j; k < nDim; ++k )
                        {
                            value_type t = __hessian( nDim*nDim*I*( nDim*k+j )+nDim*nDim*i+nDim*j+k, q );
                            M_hessian[i][q]( j,k ) = t;
                            M_hessian[i][q]( k,j ) = t;
                        }
                for ( index q = 0; q < Q; ++q )
                    for ( index j = 0; j < nDim; ++j )
                    {
                        value_type t = __hessian( nDim*nDim*I*( nDim*j+j )+nDim*nDim*i+nDim*j+j, q );
                        M_laplacian[i][q]( 0,0 ) += t;
                    }
            }

        }
        void
        init( reference_element_ptrtype const& M_ref_ele,
              matrix_node_t_type const& __pts,
              mpl::int_<1> )
        {
            typedef typename grad_type::index index;
#if 0
            std::cout << "family = " << M_ref_ele->familyName() << std::endl;
            std::cout << "is_product = " << reference_element_type::is_product << std::endl;
            std::cout << "nbDof = " << M_ref_ele->nbDof() << std::endl;
#endif
            const index I = ( reference_element_type::is_product?M_ref_ele->nbDof()/nRealDim/nRealDim:M_ref_ele->nbDof()/nRealDim );
            //std::cout << "I = " << I << std::endl;
            //std::cout << "nbDof = " << M_ref_ele->nbDof() << std::endl;

            const index Q = __pts.size2();
            //std::cout << "Q = " << I << std::endl;

            const int ncdof= ( reference_element_type::is_product?nRealDim:1 );
            int nldof= ( reference_element_type::is_product?I*nRealDim:I );
            //std::cout << "ncdof = " << ncdof << ", nldof = " << nldof << "\n";
            M_phi.resize( boost::extents[nldof][__pts.size2()] );
            M_grad.resize( boost::extents[nldof][__pts.size2()] );

            matrix_type phiv = M_ref_ele->evaluate( __pts );
            ublas::vector<matrix_type> __grad( M_ref_ele->derivate( __pts ) );

            for ( index i = 0; i < I; ++i )
            {

                for ( index c1 = 0; c1 < ncdof; ++c1 )
                {
                    for ( index q = 0; q < Q; ++q )
                        for ( index j = 0; j < nRealDim; ++j )
                            M_phi[I*c1+i][q]( j,0 ) = phiv( nldof*c1+nRealDim*i+j, q );

                    //M_phi[I*c1+i][q](j,0) = phiv( nDim*I*c1+nDim*i+j, q );

                    for ( index q = 0; q < Q; ++q )
                        for ( index j = 0; j < nRealDim; ++j )
                            for ( index l = 0; l < nDim; ++l )
                            {
                                //M_grad[I*c1+i][j](l,q) = __grad[l]( nDim*I*c1+nDim*i+j, q );
                                M_grad[I*c1+i][q]( j,l ) = __grad[l]( nldof*c1+nRealDim*i+j, q );
                                //M_grad[I*c1+i][j][nRealDim-1][q] = __grad[l]( nldof*c1+nRealDim*i+j, q );
                                //std::cout << "grad(" << i << "," << c1 << "," << j << "," << l << "," << q << ")=" <<  M_grad[I*c1+i][j](l,q) << "\n";
                            }
#if 0
                    for ( index q = 0; q < Q; ++q )
                        for ( index j = 0; j < nDim; ++j )
                        {
                            value_type t = __hessian( nDim*nDim*I*( nDim*j+j )+nDim*nDim*i+nDim*j+j, q );
                            M_laplacian[i][q]( 0,0 ) += t;
                        }
#endif //

                }
            }
        }
    private:

        reference_element_ptrtype M_ref_ele;
        matrix_node_t_type M_nodes;
        functionvalue_type M_phi;
        grad_type M_grad;
        hessian_type M_hessian;
        laplacian_type M_laplacian;
    }; /** class PreCompute **/

    typedef PreCompute precompute_type;
    typedef boost::shared_ptr<precompute_type> precompute_ptrtype;

    precompute_ptrtype
    preCompute( self_ptrtype p, points_type const& P )
    {
        return precompute_ptrtype( new PreCompute( p, P ) );
    }

    precompute_ptrtype
    preCompute( points_type const& P )
    {
        return precompute_ptrtype( new PreCompute( super_enable_shared_from_this::shared_from_this(), P ) );
    }

    typedef std::vector<std::map<typename convex_type::permutation_type, precompute_ptrtype> > faces_precompute_type;

    std::vector<std::map<typename convex_type::permutation_type, precompute_ptrtype> >
    preComputeOnFaces( self_ptrtype p, points_type const& P )
    {
#if 0
        QuadMapped<pointset_type> qm;
        typedef typename QuadMapped<pointset_type>::permutation_type permutation_type;
        typename QuadMapped<pointset_type>::permutation_points_type ppts( qm( P ) );
#endif
        typedef typename convex_type::permutation_type permutation_type;
        std::vector<std::map<permutation_type, precompute_ptrtype> > geopc( convex_type::numTopologicalFaces );

        for ( uint16_type __f = 0; __f < convex_type::numTopologicalFaces; ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                geopc[__f][__p] = precompute_ptrtype(  new precompute_type( p, P ) );
            }
        }

        return geopc;
    }


    template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g = context_v>
    class Context
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        static const int context = context_v;
        // reference space dimension
        static const uint16_type PDim = ElementType::nDim;
        // real space dimension
        static const uint16_type NDim = ElementType::nRealDim;
        static const uint16_type nDof = Basis_t::nLocalDof;
        static const bool is_product = Basis_t::is_product;
        static const bool is_scalar = Basis_t::is_scalar;
        static const bool is_vectorial = Basis_t::is_vectorial;
        static const bool is_tensor2 = Basis_t::is_tensor2;
        static const uint16_type nComponents = Basis_t::nComponents;
        static const uint16_type nComponents1 = Basis_t::nComponents1;
        static const uint16_type nComponents2 = Basis_t::nComponents2;

        static const bool is_hdiv_conforming = Feel::is_hdiv_conforming<Basis_t>::value;
        static const bool is_hcurl_conforming = Feel::is_hcurl_conforming<Basis_t>::value;

        static const bool do_opt= ( nOrder<=1 ) && ( Geo_t::nOrder==1 ) && ( convex_type::is_simplex );

        typedef typename Basis_t::polyset_type polyset_type;
        static const uint16_type rank = polyset_type::rank;

        typedef Basis_t reference_element_type;
        typedef boost::shared_ptr<Basis_t> reference_element_ptrtype;

        typedef typename reference_element_type::value_type value_type;

        typedef ElementType geometric_element_type;
        typedef typename Geo_t::template Context<context_g, ElementType> geometric_mapping_context_type;
        typedef boost::shared_ptr<geometric_mapping_context_type> geometric_mapping_context_ptrtype;

        typedef typename node<value_type>::type node_type;

        typedef ublas::matrix<value_type> matrix_type;
        typedef points_type matrix_node_t_type;
#if 0
        typedef value_type dn_type;
        typedef value_type grad_type;
        typedef value_type div_type;


        typedef typename matrix_type::value_type phi_type;
        typedef typename matrix_type::value_type dphi_type;
        typedef typename matrix_type::value_type id_type;
#else
        typedef Eigen::Matrix<value_type,nComponents1,1> id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nDim> ref_grad_type;
        typedef Eigen::Matrix<value_type,nComponents1,NDim> grad_type;
        typedef Eigen::Matrix<value_type,NDim,NDim> hess_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> laplacian_type;
        typedef Eigen::Matrix<value_type,1,1> div_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> dn_type;
        typedef Eigen::Matrix<value_type,3,1> curl_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> dx_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> dy_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> dz_type;

#if 0
        typedef boost::multi_array<id_type,2> functionvalue_type;
        typedef boost::multi_array<g_type,2> grad_type;
        typedef boost::multi_array<h_type,2> hessian_type;
        typedef boost::multi_array<n_type,2> dn_type;
        typedef boost::multi_array<d_type,2> div_type;
#endif
#endif

        typedef geometric_mapping_context_type gmc_type;
        typedef Eigen::Matrix<value_type,Eigen::Dynamic, Eigen::Dynamic> matrix_eigen_type;
        typedef Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::PDim> matrix_eigen_NP_type;
        typedef Eigen::Matrix<value_type,gmc_type::PDim,gmc_type::NDim> matrix_eigen_PN_type;
        typedef Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::NDim> matrix_eigen_NN_type;
        typedef Eigen::Matrix<value_type,nComponents1,NDim> matrix_eigen_grad_type;
        typedef typename Eigen::Map<const Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::PDim,Eigen::RowMajor> > matrix_eigen_ublas_NP_type;
        typedef typename Eigen::Map<const Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > matrix_eigen_ublas_type;


        template<uint16_type TheRank = polyset_type::rank+2>
        struct Index
        {
            static const uint16_type rank = TheRank;
            typedef boost::array<size_type,rank> index_type;
            typedef boost::detail::multi_array::extent_gen<rank> extents_type;
            template<uint16_type> friend class Rank;
            Index()
            {
                //std::cout << "Index::rank = " << rank << "\n";
                init( mpl::int_<rank>() );
            }
            Index( Index const& i )
                :
                M_index( i.M_index ),
                M_extents( i.M_extents ),
                M_comp( 0 )
            {
                //std::cout << "Index::rank = " << rank << "\n";
            }

            // rank + 1
            Index( Index<rank-1> const& __index )
                :
                M_extents( __index.extents()[RankUp<polyset_type>::type::nComponentsLast] )
            {
                std::copy( __index.beginIndex(), __index.endIndex(), M_index.begin() );
                //std::cout << "[rankup] index<" << rank << ">: ";
                //std::for_each( M_index.begin(), M_index.end(), std::cout << lambda::_1 << " " );
                //std::cout << "\n";
            }
#if 0
            // rank reduction wrt the last extent
            Index( Index<rank+1> const& index )
                :
                M_extents()
            {
#if 0
                // copy only up to rank
                std::copy( index.M_extents.begin(),
                           boost::prior( index.M_extents.end() ),
                           M_extents.begin() );
#endif
            }
#endif
            ~Index() {}

            Index const& operator=( Index const& i )
            {
                if ( this != &i )
                {
                    M_index = i.M_index;
                    M_extents = i.M_extents;
                }

                return *this;
            }

            typename index_type::iterator beginIndex()
            {
                return M_index.begin();
            }
            typename index_type::const_iterator beginIndex() const
            {
                return M_index.begin();
            }
            typename index_type::iterator endIndex()
            {
                return M_index.end();
            }
            typename index_type::const_iterator endIndex() const
            {
                return M_index.end();
            }

            template<typename Tuple>
            void setIndex( Tuple const& tu )
            {
                setIndex( tu, mpl::int_<rank>() );
            }

            void setIndex( uint16_type c, size_type i )
            {
                M_index[c] = i;
            }


            size_type index() const
            {
                return index( mpl::int_<rank>() );
            }

            size_type div() const
            {
                return nComponents*nDof*M_index[1] + nComponents*M_index[0] + M_index[1];
            }

            uint16_type component() const
            {
                return M_comp;
            }

            Index<rank+1> rankUp() const
            {
                return Index<rank+1>( M_extents[RankUp<polyset_type>::type::nComponentsLast] );
            }

            operator size_type() const
            {
                return index( mpl::int_<rank>() );
            }

            extents_type extents() const
            {
                return M_extents;
            }

        private:
            void init( mpl::int_<2> )
            {
                M_extents = boost::extents[nDof][nComponents];
                M_index[0] = size_type( -1 );
            }
            void init( mpl::int_<3> )
            {
                M_extents = boost::extents[nDof][nComponents][nComponents];
                M_index[0] = size_type( -1 );
                M_index[1] = size_type( -1 );
                M_index[2] = size_type( -1 );
            }
            void init( mpl::int_<4> )
            {
                M_extents = boost::extents[nDof][nComponents][nComponents][nComponents];
                M_index[0] = size_type( -1 );
                M_index[1] = size_type( -1 );
                M_index[2] = size_type( -1 );
            }
            size_type index( mpl::int_<2> ) const
            {
                return M_index[0];
            }
            size_type index( mpl::int_<3> ) const
            {
                return nDof*M_index[1] + M_index[0];// + M_index[2];
            }
            size_type index( mpl::int_<4> ) const
            {
                const size_type d = M_extents.ranges_[1].size();
                const size_type d2 = d*d;
                const size_type n = M_extents.ranges_[0].size();
#if 0
                return ( d2*n*( d*M_index[1]+M_index[2] ) +
                         d2*M_index[0] +
                         d*M_index[1] + M_index[2] );
#else
                return ( d2*n*M_index[1]*M_index[2] +
                         d2*M_index[0]*M_index[2] +
                         d*M_index[1]*M_index[2] +
                         d2*n*M_index[1]+
                         d2*M_index[0] +
                         d*M_index[1] +
                         M_index[2] );
#endif
            }
            template<typename Tuple>
            void setIndex( Tuple const& tu, mpl::int_<2> )
            {
                M_index[ 0 ] = boost::get<0>( tu );
                M_comp = 0;
            }
            template<typename Tuple>
            void setIndex( Tuple const& tu, mpl::int_<3> )
            {
                M_index[ 0 ] = boost::get<0>( tu );
                M_index[ 1 ] = boost::get<1>( tu );
                M_index[ 2 ] = boost::get<2>( tu );
                M_comp = M_index[ 1 ];
            }
            template<typename Tuple>
            void setIndex( Tuple const& tu, mpl::int_<4> )
            {
                M_index[ 0 ] = boost::get<0>( tu );
                M_index[ 1 ] = boost::get<1>( tu );
                M_index[ 2 ] = boost::get<2>( tu );
            }

        private:
            index_type M_index;

            extents_type M_extents;

            uint16_type M_comp;
        };

        Context( reference_element_ptrtype const& __RefEle,
                 geometric_mapping_context_ptrtype const& __gmc,
                 precompute_ptrtype const& __pc )
            :
            M_pc( __pc ),
            M_npoints( __pc->nPoints() ),

            M_ipt( 0 ),
            M_ref_ele( __RefEle ),

            M_gmc( __gmc ),
            M_phi( __pc->phi() ),
            M_gradphi( __pc->gradPtr() ),
            M_dn(),
            M_grad(),
            M_dx(),
            M_dy(),
            M_dz()
        {
            LOG(INFO) << " Polynomial derivatives optimized for P1: " << do_opt;
            if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
            {
                const int ntdof = nDof*nComponents1;
                if ( do_opt )
                    M_grad.resize( boost::extents[ntdof][1] );
                else
                    M_grad.resize( boost::extents[ntdof][M_npoints] );

                if ( vm::has_first_derivative_normal<context>::value )
                {
                    M_dn.resize( boost::extents[ntdof][M_npoints] );
                }

                if ( vm::has_div<context>::value )
                {
                    if ( do_opt )
                        M_div.resize( boost::extents[ntdof][1] );
                    else
                        M_div.resize( boost::extents[ntdof][M_npoints] );
                }

                if ( vm::has_curl<context>::value )
                {
                    if ( do_opt )
                        M_curl.resize( boost::extents[ntdof][1] );
                    else
                        M_curl.resize( boost::extents[ntdof][M_npoints] );       
                }

                if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value  )
                {
                    M_hessian.resize( boost::extents[ntdof][M_npoints] );
                }
                if ( vm::has_laplacian<context>::value )
                {
                    M_laplacian.resize( boost::extents[ntdof][M_npoints] );
                }
            }

            update( __gmc );
        }
        Context( Context const& c )
            :
            M_pc( c.M_pc ),
            M_npoints( c.M_npoints ),

            M_ipt( c.M_ipt ),
            M_ref_ele( c.M_ref_ele ),

            M_gmc( c.M_gmc ),
            M_phi( c.M_phi ),
            M_gradphi( c.M_gradphi ),
            M_dn( c.M_dn ),
            M_grad( c.M_grad ),
            M_dx( c.M_dx ),
            M_dy( c.M_dy ),
            M_dz( c.M_dz )
            {

            }
        virtual ~Context() {}
        /**
         * if isTransformationEquivalent is set to true in basis then no
         * transformation is required
         */
        void transformationEquivalence( geometric_mapping_context_ptrtype const& __gmc,
                                        mpl::bool_<true> )
        {
            // M_phi = phi;
            // if ( vm::has_grad<context>::value || vm::has_first_derivative<context>::value  )
            //     M_gradphi = M_pc->grad();
            // if ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value  )
            //     M_hessphi = M_pc->hessian();
        }

        /**
         * if isTransformationEquivalent is set to false in basis then a
         * transformation is required to ensure that the basis function are
         * equivalent on the reference and real elements.
         * We deleguate the transformation to the basis
         */
        void transformationEquivalence( geometric_mapping_context_ptrtype const& __gmc,
                                        mpl::bool_<false> )
        {
#if 0
            //M_ref_ele->transform( __gmc, phi, M_phi, M_gradphi, M_hessphi );
            M_ref_ele->transform( __gmc, M_pc.get(), M_phi,
                                   M_gradphi, ( vm::has_div<context>::value || vm::has_curl<context>::value || vm::has_grad<context>::value || vm::has_first_derivative<context>::value ),
                                   M_hessphi, ( vm::has_hessian<context>::value || vm::has_second_derivative<context>::value  )
                                 );
#endif
        }

        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc,
                                        precompute_ptrtype const& __pc );

        void update( geometric_mapping_context_ptrtype const& __gmc )
        {
            //M_phi = M_pc->get()->phi();
            //M_gradphi = M_pc->get()->grad();
            transformationEquivalence( __gmc, mpl::bool_<Basis_t::isTransformationEquivalent>() );
#if 0

            for ( int i = 0; i < M_gradphi.num_elements(); ++i )
                std::cout << "M_gradphi[" << i << "]=" << *( M_gradphi.data()+i ) << "\n";

#endif
            update( __gmc, mpl::int_<rank>(), mpl::bool_<do_opt>() );
            //update( __gmc, mpl::int_<rank>(), mpl::bool_<false>() );
        }
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, mpl::bool_<true> );
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<0>, mpl::bool_<false> );
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, mpl::bool_<true> );
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<1>, mpl::bool_<false> );
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<2>, mpl::bool_<true> );
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, mpl::int_<2>, mpl::bool_<false> );

        /**
         * \return the number of points at which the basis functions
         * have been evaluated
         */
        uint16_type nPoints() const
        {
            return M_npoints;
        }

        /**
         *  \return the geometric mapping context
         */
        geometric_mapping_context_ptrtype const& gmContext() const
        {
            return M_gmc;
        }

        //! \return the element id
        size_type eId() const
        {
            return M_gmc->id();
        }

        //! \return the points in the reference element
        matrix_node_t_type const& xRefs() const
        {
            return M_gmc->xRefs();
        }

        //! \return the precomputation of the basis function in the reference element
        precompute_ptrtype const& pc() const
        {
            return M_pc.get();
        }

        boost::multi_array<value_type,4> const& id() const
        {
            return M_phi;
        }

        value_type const& id( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q  ) const
        {
            return id( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type const& id( uint32_type i,
                              uint16_type /*c1*/,
                              uint16_type /*c2*/,
                              uint32_type q,
                              mpl::int_<0> ) const
        {
            return M_phi[i][q]( 0,0 );
        }

        value_type const& id( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q,
                              mpl::int_<1> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c2 );
            return M_phi[i][q]( c1,0 );
        }

        id_type const& id( uint32_type i, uint32_type q ) const
        {
            return M_phi[i][q];
        }

        value_type const& id( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q,
                              mpl::int_<2> ) const
        {
            return M_phi[i][q]( c1,c2 );
        }

        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return d( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0]( c1,c2 );
            }
        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false>  ) const
            {
                return M_grad[i][q]( c1,c2 );
            }

        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return dx( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][0]( c1,0 );
            }
        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][q]( c1,0 );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return dy( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][0]( c1,1 );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][q]( c1,1 );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return dz( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][0]( c1,2 );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                BOOST_MPL_ASSERT_MSG( nDim >= 1, INVALID_DIM, ( mpl::int_<nDim>, mpl::int_<1> ) );
                return M_grad[i][q]( c1,2 );
            }

        /**
         * \return the matrix containing the value of the normal
         * derivative of the basis function \p i at the node \p q
         */
        value_type const& dn( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q  ) const
        {
            return M_dn[i][q]( c1,c2 );
        }

        dn_type const& dn( uint16_type i, uint32_type q ) const
        {
            return M_dn[i][q];
        }
        grad_type const& grad( uint16_type i, uint32_type q ) const
            {
                return grad( i, q, mpl::bool_<do_opt>() );
            }
        grad_type const& grad( uint16_type i, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0];
            }
        grad_type const& grad( uint16_type i, uint32_type q, mpl::bool_<false> ) const
            {
                return M_grad[i][q];
            }

        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return grad(i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0]( c1,c2 );
            }
        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                return M_grad[i][q]( c1,c2 );
            }
        dx_type  dx( uint16_type i, uint32_type q ) const
            {
                return dx( i, q, mpl::bool_<do_opt>() );
            }
        dx_type  dx( uint16_type i, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0].col(0);
            }
        dx_type  dx( uint16_type i, uint32_type q, mpl::bool_<false> ) const
            {
                return M_grad[i][q].col(0);
            }
        dy_type  dy( uint16_type i, uint32_type q ) const
            {
                return dy( i, q, mpl::bool_<do_opt>() );
            }
        dy_type  dy( uint16_type i, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0].col(1);
            }
        dy_type  dy( uint16_type i, uint32_type q, mpl::bool_<false> ) const
            {
                return M_grad[i][q].col(1);
            }
        dz_type  dz( uint16_type i, uint32_type q ) const
            {
                return dz( i, q, mpl::bool_<do_opt>() );
            }
        dz_type  dz( uint16_type i, uint32_type q, mpl::bool_<true> ) const
            {
                return M_grad[i][0].col(1);
            }
        dz_type  dz( uint16_type i, uint32_type q, mpl::bool_<false> ) const
            {
                return M_grad[i][q].col(1);
            }

        /**
         * divergence of the basis function at the q-th point.
         *
         * \param q index of the points to evaluate the divergence
         * \param i index containing current function and component indices
         * \return divergence of the basis function at the q-th point
         */

        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return div( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type div( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            /**
             * divergence of a scalar function is undefined.
             */
            CHECK( 0 ) << "divergence of a scalar function is undefined.";
            return 0;
        }

        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<1> ) const
            {
                return div( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][0]( 0,0 );
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][q]( 0,0 );
            }

        /**
         * curl of the basis function at the q-th point.
         *
         * \param q index of the points to evaluate the curl
         * \param i index containing current function and component indices
         * \return curl of the basis function at the q-th point
         */

        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curl( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type curl( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curl operator, field must be vectorial" );
            return 0;
        }

        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<1> ) const
            {
                return curl( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( c1 );
            }
        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( c1 );
            }


        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curlx( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type curlx( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curlx operator, field must be vectorial" );
            return 0;
        }

        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<1> ) const
            {
                return curlx( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 0 );
            }
        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 0 );
            }

        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curly( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type curly( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curly operator, field must be vectorial" );
            return 0;
        }

        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<1> ) const
            {
                return curly( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 1 );
            }
        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 1 );
            }

        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curlz( i, c1, c2, q, mpl::int_<rank>() );
        }

        value_type curlz( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curlz operator, field must be vectorial" );
            return 0;
        }

        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<1> ) const
            {
                return curlz( i, c1, c2, q, mpl::bool_<do_opt>() );
            }
        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 2 );
            }
        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::bool_<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 2 );
            }


        value_type hess( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            return hess( i, c1, c2, q, mpl::int_<rank>() );
        }
        value_type hess( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, mpl::int_<0> ) const
        {
            return M_hessian[i][q]( c1,c2 );
        }
        value_type laplacian( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                return M_laplacian[i][q](c1,c2);
            }
        //    private:
        Context() {}

        void print()
            {
                const uint16_type Q = M_npoints;
                const uint16_type I = nDof;
                for ( uint16_type i = 0; i < I; ++i )
                {
                    for ( uint16_type q = 0; q < Q; ++q )
                    {
                        LOG(INFO) << "M_phi[" << i << "," << q << "]=" << M_phi[i][q];
                    }
                }

            }
    private:

        boost::optional<precompute_ptrtype> M_pc;
        uint16_type M_npoints;

        uint16_type M_ipt;

        reference_element_ptrtype M_ref_ele;

        geometric_mapping_context_ptrtype M_gmc;

        boost::multi_array<id_type,2> M_phi;
        boost::multi_array<ref_grad_type,2> const* M_gradphi;
        boost::multi_array<hess_type,2> M_hessphi;
        boost::multi_array<dn_type,2> M_dn;
        boost::multi_array<grad_type,2> M_grad;
        boost::multi_array<dx_type,2> M_dx;
        boost::multi_array<dy_type,2> M_dy;
        boost::multi_array<dz_type,2> M_dz;
        boost::multi_array<div_type,2> M_div;
        boost::multi_array<curl_type,2> M_curl;
        boost::multi_array<hess_type,2> M_hessian;
        boost::multi_array<id_type,2> M_laplacian;
    };

    template<size_type context_v, typename BasisType, typename GeoType, typename ElementType>
    boost::shared_ptr<Context<context_v,BasisType, GeoType, ElementType> >
    context( boost::shared_ptr<BasisType> b, boost::shared_ptr<GeoType> gm, precompute_ptrtype& pc )
    {
        return boost::shared_ptr<Context<context_v,BasisType, GeoType, ElementType> >(
                   new Context<context_v, BasisType, GeoType, ElementType>( context_v,
                           b,
                           gm,
                           pc ) );
    }

    template<int contextv, int contextg, typename BasisType, typename GeoType,typename ElementType>
    boost::shared_ptr<Context<contextv,BasisType, GeoType, ElementType> >
    ctx( boost::shared_ptr<BasisType> const& b,
         boost::shared_ptr<typename GeoType::template Context<contextg, ElementType> > const& gm,
         precompute_ptrtype pc, ElementType& e )
    {
        typedef Context<contextv,BasisType, GeoType, ElementType> ctx_type;
        return boost::shared_ptr<ctx_type>( new ctx_type( b, gm, pc ) );

    }
    template<int contextv, typename BasisType, typename GeoType,typename ElementType>
    boost::shared_ptr<Context<contextv,BasisType, GeoType, ElementType> >
    ctx( boost::shared_ptr<BasisType> const& b,
         boost::shared_ptr<typename GeoType::template Context<contextv, ElementType> > const& gm,
         precompute_ptrtype pc, ElementType& e )
    {
        typedef Context<contextv,BasisType, GeoType, ElementType> ctx_type;
        return boost::shared_ptr<ctx_type>( new ctx_type( b, gm, pc ) );

    }

protected:

private:

private:

    basis_type M_basis;
    matrix_type M_coeff;
    std::string M_fname;
};

template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_scalar;
template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_vectorial;
template<typename Poly,template<uint16_type> class PolySetType> const bool PolynomialSet<Poly,PolySetType>::is_tensor2;
template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents;
template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents1;
template<typename Poly,template<uint16_type> class PolySetType> const uint16_type PolynomialSet<Poly,PolySetType>::nComponents2;

} // Feel

#include <feel/feelpoly/polynomialset_impl.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>

#endif /* __PolynomialSet_H */
