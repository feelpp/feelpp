// SPDX-FileCopyrightText: 2005-2006 EPFL
// SPDX-FileCopyrightText: 2011 Universite de Grenoble
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
   \file polynomial.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __Polynomial_H
#define __Polynomial_H 1

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on


#include <type_traits>
#include <variant>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/polynomialbasis.hpp>
#include <feel/feelpoly/polynomialcoefficients.hpp>
#include <feel/feelpoly/policy.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

template<typename, template<uint16_type> class PolySetType, int OrderSpec> class PolynomialSet;

/**
 * \class Polynomial
 * \brief polynomial class
 *
 * The polynomial is expressed in the basis from \p Poly. The
 * coefficients of the polynomial in this basis are represented by a
 * matrix whose lines are the polymomial components coefficients (1 if
 * \code is_scalar == true \endcode, \p nDim if \code is_vectorial ==
 * true \endcode and columns are the basis
 *
 * Evaluating the polynomial at a set of points(or just one point) is
 * then simply a matrix-matrix product.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<typename Poly,
         template<uint16_type> class PolySetType = Scalar,
         typename Container = typename Poly::basis_type::matrix_type,
         int OrderSpec = orderSpecification<Poly>()>
class Polynomial
{
public:

    /** @name Constants
     */
    //@{

    static constexpr uint16_type nDim = Poly::nDim;
    /** Authoritative signed order specification. */
    static constexpr int order_spec = OrderSpec;
    static constexpr bool is_order_static = ( OrderSpec != Dynamic );
    static constexpr bool is_order_dynamic = !is_order_static;
    static constexpr uint16_type nOrder = is_order_static ? static_cast<uint16_type>( OrderSpec ) : 1;

    //@}


    /** @name Typedefs
     */
    //@{

    typedef Polynomial<Poly, PolySetType, Container, OrderSpec> self_type;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;



    typedef PolySetType<nDim> polyset_type;
    static constexpr bool is_tensor2 = polyset_type::is_tensor2;
    static constexpr bool is_vectorial = polyset_type::is_vectorial;
    static constexpr bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;

    /** Scalar polynomial obtained by extracting a component. */
    using component_type = Polynomial<Poly, Scalar, typename Poly::basis_type::matrix_type, OrderSpec>;
    typedef Polynomial<Poly, Scalar, typename Poly::basis_type::matrix_type, OrderSpec> scalar_component_type;

    typedef typename basis_type::points_type points_type;
    typedef typename basis_type::matrix_type matrix_type;
    typedef Container container_type;
    /** Canonical coefficient owner keeps the expansion axis independent of component rows. */
    using CoefficientStorage = ::Feel::detail::PolynomialCoefficientStorage<value_type,
        ::Feel::detail::polynomialCoefficientExtent<Poly, OrderSpec>()>;

    typedef typename node<value_type>::type node_type;

    static_assert( std::is_same_v<typename matrix_type::value_type, value_type> );
    static_assert( std::is_same_v<typename matrix_type::value_type, typename points_type::value_type> );

    //@}

    /** @name Constructors, destructor
     */
    //@

    /**
     * default constructor
     */
    Polynomial()
        :
        M_basis(),
        M_coeff( M_basis.coeff() )
    {
        setOrderValue( nOrder );
    }

    /** @brief Construct a polynomial with a runtime expansion degree. */
    Polynomial( uint16_type order )
        requires( is_order_dynamic )
        :
        M_basis(),
        M_coeff( M_basis.coeff() )
    {
        setOrderValue( order );
    }


    /**
     * constructor giving only the underlying basis
     * \param __poly polynomial whose we take the basis
     */
    template<class Source>
        requires requires( Source const& p ) { p.basis(); p.coeff(); }
    Polynomial( Source const& __poly )
        :
        M_basis( __poly.basis() ),
        M_coeff( ::Feel::detail::polynomialCoefficientStorageFrom<CoefficientStorage>( __poly ) )
    {
        setOrderValue( ::Feel::detail::polynomialOrder( __poly ) );
    }

    /** @brief Construct from source-relative or modal Eigen/legacy coefficients. */
    template<class Source, typename Coefficients>
        requires requires( Source const& p ) { p.basis(); p.coeff(); } &&
                 ( EigenMatrix<Coefficients> || UBlasMatrix<Coefficients> )
    Polynomial( Source const& source, Coefficients const& coefficients, bool asIs = false )
        : M_basis( source.basis() ),
          M_coeff( ::Feel::detail::polynomialCoefficientStorageFrom<CoefficientStorage, polyset_type, true>(
              source, coefficients, asIs ) )
    {
        setOrderValue( ::Feel::detail::polynomialOrder( source ) );
    }

    /** @brief Construct coefficients in the default modal expansion. */
    template<typename Coefficients>
        requires ( EigenMatrix<Coefficients> || UBlasMatrix<Coefficients> )
    Polynomial( Coefficients const& coefficients, bool asIs = false )
        : M_basis(), M_coeff( M_basis.coeff() )
    {
        setOrderValue( nOrder );
        setCoefficient( coefficients, asIs );
    }

    /** @brief Copy coefficients and share immutable differentiation operators. */
    Polynomial( Polynomial const & p )
        : Polynomial( p, CoefficientStorage( p.M_coeff ) )
    {}

    /** @brief Move coefficients, basis and order policy together. */
    Polynomial( Polynomial&& ) = default;

    /** @brief Move-assign a polynomial without retaining its former order. */
    self_type& operator=( Polynomial&& ) = default;

    ~Polynomial()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& __p )
    {
        if ( this != &__p )
        {
            M_order = __p.M_order;
            M_basis = __p.M_basis;
            M_coeff = __p.M_coeff;
        }

        return *this;
    }

    /** @brief Subtract coefficients in the same expansion, rejecting different degrees. */
    self_type const& operator-=( self_type const& __p )
    {
        if ( order() != __p.order() )
            throw std::invalid_argument("cannot subtract polynomials with different expansion orders");
        M_coeff.assign( M_coeff.matrix() - __p.M_coeff.matrix() );
        return *this;
    }

    /**
     * \brief extract the i-th component of a vectorial polynomial
     *
     * \return the i-th component of the polynomial
     */
    component_type operator[]( int i ) const
    {
        if ( i < 0 || i >= nComponents )
            throw std::out_of_range( "polynomial component index is outside the component range" );
        return component_type( *this, M_coeff.matrix().row( i ), true );
    }

    /** @brief Extract a tensor entry as a scalar polynomial with the same degree. */
    scalar_component_type operator()( uint16_type i, uint16_type j ) const requires( is_tensor2 )
    {
        if ( i >= nComponents1 || j >= nComponents2 )
            throw std::out_of_range("tensor component index is outside the component range");
        return (*this)[i*nComponents2+j];
    }

    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    operator()( node_type const& __x ) const
    {
        return evaluate( __x );
    }

    /**
     * Evaluate polynomial at points \p __pts
     *
     * \param __pts the matrix of coordinates of the points
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type operator()( points_type const& __pts ) const
    {
        return evaluate( __pts );
    }

    //@}

    /** @name Accessors
     */
    //@{


    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
        return M_coeff.matrix().norm() < Feel::type_traits<value_type>::epsilon();

    }
    /**
     * \return the dof
     */
    matrix_type coeff() const
    {
        return M_coeff.template legacy<matrix_type>();
    }

    /**
     * \return the dof
     */
    auto coefficients() const&
    {
        return M_coeff.matrix();
    }

    /** @brief Reject a coefficient view of a temporary polynomial. */
    auto coefficients() const&& = delete;

    /**
     * \return the basis in which the polynomial is expressed
     */
    basis_type const& basis() const
    {
        return M_basis;
    }

    /** @brief Name of the underlying modal expansion family. */
    std::string familyName() const { return M_basis.familyName(); }

    /** @brief Effective expansion degree, retained under differentiation. */
    uint16_type degree() const
    {
        return orderValue();
    }

    /** @brief Effective degree selected by the static or runtime order policy. */
    uint16_type order() const
    {
        return orderValue();
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /** @brief Change runtime degree; coefficients must match the requested expansion. */
    void setOrder( uint16_type order )
        requires( is_order_dynamic )
    {
        setOrderValue( order );
    }

    /**
     * set the coefficient of the polynomial in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {
        setCoefficient( eigenMap( __c ), __as_is );
    }


    /** @brief Compose component coefficients or directly own modal Eigen coefficients. */
    template<typename Derived>
    void setCoefficient( Eigen::MatrixBase<Derived> const& coefficients, bool asIs = false )
    {
        if ( asIs ) M_coeff.assign( coefficients );
        else
        {
            decltype(auto) packed = polyset_type::toMatrix( coefficients.derived() );
            M_coeff.template set<polyset_type>( packed, false );
        }
    }

    /** @brief Materialize an unmapped legacy expression at its input boundary. */
    template<typename Derived>
    void setCoefficient( ublas::matrix_expression<Derived> const& expression, bool asIs = false )
    {
        matrix_type owned( expression );
        setCoefficient( owned, asIs );
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    evaluate( node_type const& __x ) const
    {
        points_type points( __x.size(), 1 );
        ublas::column( points, 0 ) = __x;
        return evaluate( points );
    }

    /** @brief Evaluate this polynomial at a batch of points. */
    matrix_type evaluate( points_type const& points ) const
    {
        return M_order.withValues( M_basis, points, [&]( auto const& modal )
        { return M_coeff.template multiplyAs<matrix_type>( modal ); } );
    }

    /** @brief Evaluate C Phi in canonical Eigen storage. */
    template<typename Points>
    [[nodiscard]] auto values( Points const& points ) const
    {
        return M_order.withValues( M_basis, points, [&]( auto const& modal )
        { return M_coeff.multiply( modal ); } );
    }

    /** @brief Own requested reference values and derivatives through one modal recurrence. */
    template <Evaluation Request = Evaluation::All, typename Points>
    [[nodiscard]] auto tabulate( Points const& points ) const
    {
        return ::Feel::tabulate<Request>( *this, points );
    }

    /** @brief Reuse caller-owned reference tables and modal recurrence storage. */
    template <typename Points, typename Tables, typename Workspace>
    void tabulateInto( Points const& points, Tables& output, Workspace& workspace ) const
    {
        ::Feel::tabulateInto( *this, points, output, workspace );
    }

    /** @brief Evaluate one reference gradient direction directly in the legacy result. */
    template<typename AE>
    matrix_type derivate( uint16_type i, ublas::matrix_expression<AE> const& pts ) const
    {
        if ( i >= nDim ) throw std::out_of_range( "invalid polynomial derivative direction" );
        matrix_type result;
        M_order.forEachGradient( M_basis, pts(), [&]( int direction, auto const& modal )
        {
            if ( direction == i ) result = M_coeff.template multiplyAs<matrix_type>( modal );
        } );
        return result;
    }

    /** @brief Evaluate a mixed second reference derivative using native modal values. */
    template<typename AE>
    matrix_type derivate( uint16_type i, uint16_type j, ublas::matrix_expression<AE> const& pts ) const
    {
        if ( i >= nDim || j >= nDim ) throw std::out_of_range( "invalid polynomial Hessian direction" );

        // Retain one coefficient-operator path for compatibility and low-order tensors.
        auto withOperators = [&]()
        {
            auto first = M_coeff.multiply( eigenMap( d( i ) ) );
            auto second = ( first * eigenMap( d( j ) ) ).eval();
            return M_order.withValues( M_basis, pts(), [&]( auto const& modal )
            { return ::Feel::detail::materializePolynomialMatrix<matrix_type>( second * modal ); } );
        };

        // Share backend selection with polynomial sets at the same effective degree.
        if constexpr ( order_storage_type::template hasDirectHessians<AE> )
        {
            if ( M_order.useDirectHessians() )
            {
                matrix_type result;
                M_order.forEachHessian( M_basis, pts(), [&]( int first, int second, auto const& modal )
                                        {
                    if ( first == i && second == j ) result = M_coeff.template multiplyAs<matrix_type>( modal ); } );
                return result;
            }
        }

        return withOperators();
    }


    /**
     * \brief differentiation matrix of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
        return M_order.d( M_basis, i );
    }

    /** @brief Return the mixed second differentiation operator. */
    matrix_type d( uint16_type i, uint16_type j ) const
    {
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( eigenMap( d( i ) ) * eigenMap( d( j ) ) );
    }

    /**
     * derivative of the polynomial in the l-th direction
     *
     * \param l direction of derivation
     * \return the derivative of the polynomial
     * \note Retains the expansion state and shared immutable operators.
     */
    self_type derivative( uint16_type l ) const
    {
        return self_type( *this, M_coeff.rightProduct( eigenMap( d( l ) ) ) );
    }

    /** @brief Gradient or Jacobian, with one result for this polynomial. */
    auto gradient() const { return toSet().gradient().polynomial(0); }

    /**
     * transform a \p Polynomial into a PolynomialSet with one polynomial
     *
     * \return a \p PolynomialSet
     */
    PolynomialSet<Poly, PolySetType, OrderSpec> toSet( [[maybe_unused]] bool asis = false ) const
    {
        PolynomialSet<Poly, PolySetType, OrderSpec> result( *this, M_coeff.matrix(), true );
        if constexpr ( is_order_dynamic )
            result.setOrder( orderValue() );
        return result;
    }
    //@}


protected:

private:

    /** @brief Retain basis/order state while adopting independently owned coefficients. */
    Polynomial( self_type const& source, CoefficientStorage&& coefficients )
        : M_order( source.M_order ), M_basis( source.M_basis ), M_coeff( std::move( coefficients ) ) {}

    using order_storage_type = ::Feel::detail::PolynomialBasisState<basis_type, OrderSpec>;

    /** @brief Effective expansion degree from the shared order policy. */
    [[nodiscard]] uint16_type orderValue() const
    {
        if constexpr ( is_order_static )
            return nOrder;
        else
            return M_order.order();
    }

    /** @brief Set runtime degree or validate a static degree. */
    void setOrderValue( uint16_type order )
    {
        M_order.setOrder( order );
    }

    [[no_unique_address]] order_storage_type M_order{};
    basis_type M_basis;
    /** Sole Eigen owner; compatibility access materializes an independent legacy matrix. */
    CoefficientStorage M_coeff;
};
/** @brief Subtract polynomials while retaining their common basis and order policy. */
template<typename Poly, template<uint16_type> class PolySetType, typename Container, int OrderSpec>
Polynomial<Poly, PolySetType, Container, OrderSpec> operator-( Polynomial<Poly, PolySetType, Container, OrderSpec> const& p1,Polynomial<Poly, PolySetType, Container, OrderSpec> const& p2 )
{
    if ( p1.order() != p2.order() )
        throw std::invalid_argument("cannot subtract polynomials with different expansion orders");
    typename Polynomial<Poly, PolySetType, Container, OrderSpec>::matrix_type c = p1.coeff()-p2.coeff();
    //std::cout << "operator- c=" << c << "\n";
    Polynomial<Poly, PolySetType, Container, OrderSpec> result( p1, c, true );
    if constexpr ( Polynomial<Poly, PolySetType, Container, OrderSpec>::is_order_dynamic )
        result.setOrder( p1.order() );
    return result;
}

}
#endif /* __Polynomial_H */
