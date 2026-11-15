/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Gilles Steiner <gilles.steiner@epfl.ch>
 *
 * SPDX-FileCopyrightText: 2006-2011 Université Joseph Fourier Grenoble 1
 * SPDX-FileCopyrightText: 2011-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_POLYNOMIALSET_HPP
#define FEELPP_POLYNOMIALSET_HPP 1

#include <array>
#include <stdexcept>
#include <vector>
#include <feel/feelcore/boostmultiarray.hpp>

#include <boost/optional.hpp>
#include <type_traits>
#include <variant>
#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>

#include <feel/feelpoly/context.hpp>

#include <feel/feelalg/svd.hpp>

#include <feel/feelpoly/boundadapted.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feelpoly/hcurlpolynomialset.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/legendre.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/piolatransform.hpp>
#include <feel/feelpoly/polynomial.hpp>
#include <feel/feelpoly/polynomialbasis.hpp>
#include <feel/feelpoly/polynomialcoefficients.hpp>
#include <feel/feelpoly/precomputetabulation.hpp>
#include <feel/feelpoly/quadmapped.hpp>
#include <feel/feelpoly/tensorisedboundadapted.hpp>
#include <feel/feelpoly/traits.hpp>

namespace Feel
{
template<bool B>
using bool_t = std::bool_constant<B>;
template<int N>
using int_t = std::integral_constant<int, N>;

using optimizable_t = bool_t<true>;
using default_t = bool_t<false>;
template<int r>
using rank_t = int_t<r>;
using scalar_t = rank_t<0>;
using vectorial_t = rank_t<1>;
using matricial_t = rank_t<2>;

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
template<typename Poly, template<uint16_type> class PolySetType = Scalar, int OrderSpec = orderSpecification<Poly>()>
class PolynomialSet
    :
        public std::enable_shared_from_this<PolynomialSet<Poly, PolySetType, OrderSpec > >,
        public PolySetType<Poly::nRealDim>
{
    typedef std::enable_shared_from_this<PolynomialSet<Poly, PolySetType, OrderSpec > > super_enable_shared_from_this;
    template<typename, template<uint16_type> class, int>
    friend class PolynomialSet;
public:

    /** @name Constants
     */
    //@{

    static constexpr uint16_type nDim = Poly::nDim;
    static constexpr uint16_type nRealDim = Poly::nRealDim;
    /** Authoritative signed order specification. */
    static constexpr int order_spec = OrderSpec;
    static constexpr bool is_order_static = ( OrderSpec != Dynamic );
    static constexpr bool is_order_dynamic = !is_order_static;
    //! Static order value (or placeholder 1 when dynamic for template instantiation)
    static constexpr uint16_type nOrder = is_order_static ? static_cast<uint16_type>( OrderSpec ) : 1;

    //@}

    /** @name Typedefs
     */
    //@{
    typedef PolynomialSet<Poly, PolySetType, OrderSpec> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;
    static inline const bool is_product = Poly::is_product;


    typedef PolySetType<nRealDim> polyset_type;
    static constexpr bool is_tensor2symm = polyset_type::is_tensor2 && is_symm_v<polyset_type>;
    static constexpr bool is_tensor2 = polyset_type::is_tensor2;
    static constexpr bool is_vectorial = polyset_type::is_vectorial;
    static constexpr bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static const uint16_type rank = polyset_type::rank;

    typedef PolynomialSet<Poly, Scalar, OrderSpec> component_type;
    typedef Polynomial<Poly, PolySetType, typename basis_type::matrix_type, OrderSpec> polynomial_type;
    //typedef Polynomial<Poly, PolySetType, ublas::vector_range<ublas::vector<value_type> > > polynomial_type;
    typedef polynomial_type polynomial_view_type;

    typedef typename Poly::convex_type convex_type;
    typedef typename basis_type::matrix_type matrix_type;
    /** Canonical Eigen coefficient owner, independent of the number of functions. */
    using CoefficientStorage = ::Feel::detail::PolynomialCoefficientStorage<value_type,
        ::Feel::detail::polynomialCoefficientExtent<Poly, OrderSpec>()>;
    typedef typename basis_type::points_type points_type;

    using gradient_polynomialset_type = if_t<is_scalar_polynomial<polyset_type>::value,
                                             PolynomialSet<Poly, Vectorial, OrderSpec>,
                                             if_t<is_vector_polynomial<polyset_type>::value,
                                                  PolynomialSet<Poly, Tensor2, OrderSpec>,
                                                  PolynomialSet<Poly, Tensor3, OrderSpec>>>;

    static_assert( std::is_same_v<typename matrix_type::value_type, value_type> );
    static_assert( std::is_same_v<typename matrix_type::value_type, typename points_type::value_type> );

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct an empty set at the default expansion degree. */
    PolynomialSet()
        :
        M_basis(),
        M_coeff(),
        M_fname( "pset" )
    {
        setOrderValue( nOrder );
        //std::cout << "[PolynomialSet::default] dim = " << nDim << " order = " << nOrder << "\n";
    }
    /** @brief Construct an empty runtime set at the supplied degree. */
    PolynomialSet( uint16_type order )
        requires( is_order_dynamic )
        :
        M_basis(),
        M_coeff(),
        M_fname( "pset" )
    {
        setOrderValue( order );
    }
    /** @brief Copy expansion, coefficients, family and effective degree from a source. */
    template<class Source>
        requires requires( Source const& p ) { p.basis(); p.coeff(); }
    PolynomialSet( Source const& p )
        :
        M_basis( p.basis() ),
        M_coeff( ::Feel::detail::polynomialCoefficientStorageFrom<CoefficientStorage>( p ) ),
        M_fname( p.familyName() )
    {
        setOrderValue( ::Feel::detail::polynomialOrder( p ) );
    }
    /** @brief Construct from source-relative coefficients, or modal coefficients when as-is. */
    template<class Source, typename Coefficients>
        requires requires( Source const& p ) { p.basis(); p.coeff(); } &&
                 ( EigenMatrix<Coefficients> || UBlasMatrix<Coefficients> )
    PolynomialSet( Source const& p, Coefficients const& c, bool asIs = false )
        :
        M_basis( p.basis() ),
        M_coeff( ::Feel::detail::polynomialCoefficientStorageFrom<CoefficientStorage, polyset_type>( p, c, asIs ) ),
        M_fname( p.familyName() )
    {
        setOrderValue( ::Feel::detail::polynomialOrder( p ) );
        //std::cout << "[PolynomialSet] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet] c = " << c << "\n";
        //std::cout << "[PolynomialSet] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet] coeff = " << M_coeff << "\n";
    }

    /** @brief Construct from coefficients in the default modal expansion. */
    template<typename Coefficients>
        requires ( EigenMatrix<Coefficients> || UBlasMatrix<Coefficients> )
    PolynomialSet( Coefficients const& c, bool __as_is = false )
        :
        M_basis(),
        M_coeff( c ),
        M_fname( "pset" )
    {
        setOrderValue( nOrder );
        setCoefficient( c, __as_is );
        //std::cout << "[PolynomialSet] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet] c = " << c << "\n";
        //std::cout << "[PolynomialSet] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet] coeff = " << M_coeff << "\n";
    }


    /** @brief Copy all value state, sharing immutable runtime derivative operators. */
    PolynomialSet( PolynomialSet const & p )
        : PolynomialSet( p, CoefficientStorage( p.M_coeff ) )
    {
        //std::cout << "[PolynomialSet::copy] dim = " << nDim << " order = " << nOrder << "\n";
        //std::cout << "[PolynomialSet::copy] p.coeff = " << p.coeff() << "\n";
        //std::cout << "[PolynomialSet::copy] coeff = " << M_coeff << "\n";
    }

    /** @brief Destroy a polynomial set through its public interface. */
    virtual ~PolynomialSet()
    {}

    /** @brief Move the basis, coefficients, family and order policy together. */
    PolynomialSet( PolynomialSet&& ) = default;

    /** @brief Move-assign all polynomial-set value state. */
    self_type& operator=( PolynomialSet&& ) = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& pset )
    {
        if ( this != &pset )
        {
            M_order = pset.M_order;
            M_basis = pset.M_basis;
            M_coeff = pset.M_coeff;
            M_fname = pset.M_fname;
        }

        return *this;
    }

    /**
     * \brief extract the i-th component of a vectorial polynomial set
     *
     * \param i component index in `[0,nComponents)`
     * \return the i-th component of the polynomial set
     * \throws std::out_of_range if \p i is not a valid component
     */
    component_type operator[]( uint16_type i ) const
    {
        static_assert( is_vectorial, "component extraction requires a vector polynomial set" );
        return extractComponent( i );
    }

    /** @brief Extract one tensor entry for every function, preserving the expansion. */
    component_type operator()( uint16_type i, uint16_type j ) const
        {
        static_assert( is_tensor2, "tensor entry extraction requires a rank-two polynomial set" );
        if ( i >= nComponents1 || j >= nComponents2 )
            throw std::out_of_range( "tensor component index is outside the component range" );
        return extractComponent( nComponents2 * i + j );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the semantic degree of the polynomials in the set
     *
     * Static types use their compile-time degree without runtime storage.
     */
    [[nodiscard]] uint16_type degree() const
    {
        return M_order.order();
    }

    /**
     * \return the semantic polynomial order
     *
     * Static types use their compile-time degree without runtime storage.
     */
    [[nodiscard]] uint16_type order() const
    {
        return M_order.order();
    }

    /** @brief Change runtime expansion degree; coefficients must match the new extent. */
    void setOrder( uint16_type order )
        requires is_dynamic_order<OrderSpec>
    {
        setOrderValue( order );
    }

    /** @brief Preserve expansion and coefficients while adopting a runtime order policy. */
    [[nodiscard]] PolynomialSet<Poly, PolySetType, Dynamic> toDynamic() const
        requires is_static_order<OrderSpec>
    {
        PolynomialSet<Poly, PolySetType, Dynamic> dynamic_set;
        dynamic_set.M_basis = M_basis;
        dynamic_set.M_coeff.assign( M_coeff.matrix() );
        dynamic_set.M_fname = M_fname;
        dynamic_set.setOrder( order() );
        return dynamic_set;
    }

    /** @brief Convert to static degree N when the effective degree agrees. */
    template<int N>
    [[nodiscard]] boost::optional<PolynomialSet<Poly, PolySetType, N>> toStatic() const
        requires is_dynamic_order<OrderSpec>
    {
        if ( order() != static_cast<uint16_type>( N ) )
            return boost::none;
        PolynomialSet<Poly, PolySetType, N> static_set;
        static_set.M_basis = M_basis;
        static_set.M_coeff.assign( M_coeff.matrix() );
        static_set.M_fname = M_fname;
        return static_set;
    }

    /**
     * \return an owning legacy snapshot of coefficients in the expansion basis.
     * @note Prefer coefficients() for a borrowed native Eigen view.
     */
    matrix_type coeff() const
    {
        return M_coeff.template legacy<matrix_type>();
    }

    /** @brief Borrow native coefficients C(component/function, mode) until this owner changes. */
    [[nodiscard]] auto coefficients() const& { return M_coeff.matrix(); }
    /** @brief Reject a coefficient view of a temporary polynomial set. */
    auto coefficients() const&& = delete;

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
        return M_coeff.rows()/nComponents;
    }

    /**
     * \return the polynomial dimension per component
     */
    size_type polynomialDimensionPerComponent() const
    {
        return M_coeff.cols();
    }

    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
        return M_coeff.matrix().norm() < Feel::type_traits<value_type>::epsilon();
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
        os << this->familyName() << sep << nDim << sep << orderValue();
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
        setCoefficient( eigenMap( __c ), __as_is );
    }

    /** @brief Compose or directly assign native Eigen coefficients with checked dimensions. */
    template<typename Derived>
    void setCoefficient( Eigen::MatrixBase<Derived> const& coefficients, bool asIs = false )
    { M_coeff.template set<polyset_type>( coefficients, asIs ); }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Extract the polynomials whose indices are listed in \p list_p
     *
     * \param indices list of indices of polynomials to extract
     * \return the polynomial set extracted
     */
    self_type polynomials( std::vector<int> const& indices ) const
    {
        auto const dimension = polynomialDimension();
        auto const blockSize = nComponents * indices.size();
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> selected( nComponents * blockSize, M_coeff.cols() );
        for ( std::size_t j = 0; j < indices.size(); ++j )
        {
            int const index = indices[j];
            if ( index < 0 || static_cast<size_type>( index ) >= dimension / nComponents )
                throw std::out_of_range( "polynomial block index is outside the set" );
            for ( int component = 0; component < nComponents; ++component )
                selected.middleRows( component * blockSize + nComponents * j, nComponents ) =
                    M_coeff.matrix().middleRows( component * dimension + nComponents * index, nComponents );
        }
        return self_type( *this, selected, true );
    }

    /** @brief Extract the first count polynomial blocks using the shared range implementation. */
    self_type polynomialsUpToDimension( int count ) const
    {
        if ( count < 0 || static_cast<size_type>( count ) > polynomialDimension() / nComponents )
            throw std::out_of_range( "polynomial block dimension is outside the set" );
        return polynomialsRange( 0, count );
    }

    /** @brief Extract a half-open range of polynomial blocks, preserving component order. */
    self_type polynomialsRange( int begin, int end ) const
    {
        if ( begin < 0 || begin > end || static_cast<size_type>( end ) > polynomialDimension() / nComponents )
            throw std::out_of_range( "polynomial block range is outside the set" );
        std::vector<int> indices( end - begin );
        std::iota( indices.begin(), indices.end(), begin );
        return polynomials( indices );
    }

    /** @brief Extract one polynomial as an owning value in the same expansion. */
    polynomial_type polynomial( uint16_type i ) const
    {
        if ( i >= polynomialDimension() )
            throw std::out_of_range( "polynomial index is outside the set" );
        return polynomial_type( *this, M_coeff.matrix().middleRows( nComponents * i, nComponents ), true );
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
        matrix_type values = polynomial(i).evaluate( __pt );
        return ublas::column( values, 0 );
    }

    /**
     * evaluate all polynomials at node __pt
     * \return a column matrix
     */
    template<typename AE>
    ublas::matrix<value_type> evaluate( ublas::vector_expression<AE> const& __pt ) const
    {
        points_type points( __pt().size(), 1 );
        ublas::column( points, 0 ) = __pt;
        return evaluate( points );
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
     * \throws std::invalid_argument if the basis and coefficient dimensions
     * do not match
     */
    template<typename AE>
    matrix_type evaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        return M_order.withValues( M_basis, __pts(), [&]( auto const& modal )
        { return M_coeff.template multiplyAs<matrix_type>( modal ); } );
    }

    /** @brief Evaluate C Phi in Eigen storage; legacy points are borrowed by the modal boundary. */
    template <typename Points>
    [[nodiscard]] auto values( Points const& points ) const
    {
        return M_order.withValues( M_basis, points, [&]( auto const& modal )
                                   { return M_coeff.multiply( modal ); } );
    }

    /** @brief Evaluate all reference gradient directions in fixed-dimensional standard storage. */
    template <typename Points>
    [[nodiscard]] auto gradients( Points const& points ) const
    {
        using Matrix = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
        std::array<Matrix, nDim> result;
        M_order.forEachGradient( M_basis, points, [&]( int direction, auto const& modal )
                                 { result[direction] = M_coeff.multiply( modal ); } );
        return result;
    }

    /** @brief Evaluate ordered second reference derivatives in direction-major storage. */
    template <typename Points>
    [[nodiscard]] auto hessians( Points const& points ) const
    {
        using Matrix = Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>;
        std::array<Matrix, nDim * nDim> result;
        forEachHessian( points, [&]( int i, int j, auto const& expression )
                        { result[i * nDim + j].noalias() = expression; } );
        return result;
    }

    /** @brief Own requested reference tables from the common modal-to-nodal tabulator. */
    template <Evaluation Request = Evaluation::All, typename Points>
    [[nodiscard]] auto tabulate( Points const& points ) const
    {
        return ::Feel::tabulate<Request>( *this, points );
    }

    /** @brief Reuse caller-owned requested tables and modal recurrence storage. */
    template <typename Points, typename Tables, typename Workspace>
    void tabulateInto( Points const& points, Tables& output, Workspace& workspace ) const
    {
        ::Feel::tabulateInto( *this, points, output, workspace );
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
        return basisD( i );
    }

    /** @brief Return the mixed second differentiation operator at the effective degree. */
    matrix_type d( uint16_type i, uint16_type j ) const
    {
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( eigenMap( basisD( i ) ) * eigenMap( basisD( j ) ) );
    }

    /**
     * \brief Derivate with respect to the l-th direction.
     *
     * \return the polynomial set associated with the derivation in the direction
     * \note Copies retain the effective order and shared immutable operators.
     */
    self_type derivate( uint16_type direction ) const
    {
        return self_type( *this, M_coeff.rightProduct( eigenMap( basisD( direction ) ) ) );
    }

    /** @brief Mathematical derivative in a reference coordinate direction. */
    self_type derivative( uint16_type direction ) const { return derivate( direction ); }

    /** @brief Evaluate every reference gradient direction in owning legacy matrices. */
    template<typename AE>
    ublas::vector<matrix_type> derivate( ublas::matrix_expression<AE> const& pts ) const
    {
        ublas::vector<matrix_type> result( nDim );
        M_order.forEachGradient( M_basis, pts(), [&]( int direction, auto const& modal )
        { result[direction] = M_coeff.template multiplyAs<matrix_type>( modal ); } );
        return result;
    }

    /** @brief Evaluate one reference gradient direction directly into its final matrix. */
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

    /** @brief Evaluate all ordered reference Hessians through the shared destination path. */
    template<typename AE>
    ublas::matrix<matrix_type> derivate2( ublas::matrix_expression<AE> const& pts ) const
        {
        ublas::matrix<matrix_type> result( nDim, nDim );
        forEachHessian( pts(), [&]( int i, int j, auto const& expression )
        { result( i, j ) = ::Feel::detail::materializePolynomialMatrix<matrix_type>( expression ); } );
        return result;
    }
    /** @brief Evaluate an ordered second reference derivative with checked directions. */
    template<typename AE>
    matrix_type derivate( uint16_type i, uint16_type j, ublas::matrix_expression<AE> const& pts ) const
    {
        if ( i >= nDim || j >= nDim ) throw std::out_of_range( "invalid polynomial Hessian direction" );
        auto native = hessians( pts );
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( native[i * nDim + j] );
    }
    /**
     * Gradient of the polynomial set
     *
     * Computes the gradient of the polynomial set.
     */
    gradient_polynomialset_type
    gradient() const
    {
        return gradient( int_t<polyset_type::rank>() );
    }

    gradient_polynomialset_type
    gradient( rank_t<0> ) const
    {
        return gradientCoefficients();
    }

    /** @brief Form the Jacobian of each vector polynomial. */
    gradient_polynomialset_type gradient( rank_t<1> ) const
    {
        return gradientCoefficients();
    }

private:
    /** @brief Extract one interleaved component using the common coefficient indexing rule. */
    component_type extractComponent( uint16_type index ) const
    {
        if ( index >= nComponents ) throw std::out_of_range( "polynomial component index is outside the component range" );
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> selected( M_coeff.rows() / nComponents, M_coeff.cols() );
        for ( Eigen::Index row = 0; row < selected.rows(); ++row )
            selected.row( row ) = M_coeff.matrix().row( row * nComponents + index );
        return component_type( *this, selected, true );
    }

    /** @brief Assemble gradients in function/component/direction order. */
    gradient_polynomialset_type gradientCoefficients() const
    {
        static_assert( nDim == nRealDim, "reference gradient requires equal topological and real dimensions" );
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> coefficients( M_coeff.rows() * nDim, M_coeff.cols() );
        for ( int direction = 0; direction < nDim; ++direction )
        {
            auto differentiated = M_coeff.multiply( eigenMap( basisD( direction ) ) );
            for ( Eigen::Index row = 0; row < M_coeff.rows(); ++row )
                coefficients.row( row * nDim + direction ) = differentiated.row( row );
        }
        return gradient_polynomialset_type( *this, coefficients, true );
    }
public:



    /**
     * \return the number of degrees of freedom
     */
    uint16_type nbDof() const
    {
        return M_coeff.rows();
    }


    /**
     * insert the polynomial set \p p at the end of the set
     * \param p polynomial set to append or assign
     * \param erase replace the current coefficients when true
     * \throws std::invalid_argument if the coefficient dimensions differ
     */
    void insert( self_type const& p, bool erase = false )
    {
        if ( M_coeff.cols() != p.M_coeff.cols() )
            throw std::invalid_argument( "polynomial insertion requires matching expansion dimensions" );
        if ( erase ) { M_coeff = p.M_coeff; return; }
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> combined( M_coeff.rows() + p.M_coeff.rows(), M_coeff.cols() );
        combined.topRows( M_coeff.rows() ) = M_coeff.matrix();
        combined.bottomRows( p.M_coeff.rows() ) = p.M_coeff.matrix();
        M_coeff.assign( combined );
    }

    //!
    //! compute the first moment value or first-moment of the polynomial set.
    //! @note the values are not divided by the measure of the domain.
    //! @return a column matrix containing the first moment value of the basis functions
    //!
    matrix_type firstMoments()
        {
        points_type points( nDim, 1 );
        eigenMap( points ).setZero();
        auto modal = basisEvaluate( points );
        return ::Feel::detail::materializePolynomialMatrix<matrix_type>( M_coeff.matrix().col( 0 ) / modal( 0, 0 ) );
    }

    //@}


    /******************************************************/

    /** @brief Requested reference tabulation adapted to existing FE/geometry layouts.
     * Native column-major tables are canonical owners. Compatibility tensors
     * use function/component/direction indexing without re-evaluating the basis.
     * Updates reuse storage; borrowed native views expire on update or destruction.
     */
    class PreCompute
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef PolynomialSet<Poly, PolySetType, OrderSpec > reference_element_type;
        typedef std::shared_ptr<reference_element_type> reference_element_ptrtype;

        typedef typename reference_element_type::value_type value_type;
        using polyset_type = typename reference_element_type::polyset_type;

        static const uint16_type nDim = reference_element_type::nDim;
        static const uint16_type nComponents1 = reference_element_type::nComponents1;
        static const uint16_type nComponents2 = reference_element_type::nComponents2;
        //static const uint16_type nComponents3 = reference_element_type::nComponents3;
        static const uint16_type nComponents = reference_element_type::nComponents;

        typedef typename reference_element_type::points_type matrix_node_t_type;
        typedef ublas::matrix<value_type> matrix_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>> id_type;

        static constexpr int g_comp2=(rank==2)?nComponents2:(rank==1)?nDim:nDim;
        static constexpr int g_comp3=(rank==2)?nDim:(rank==1)?1:1;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,g_comp2,g_comp3>> g_type;

        static constexpr int h_comp1=(rank==2)?nComponents1:(rank==1)?nComponents1:nDim;
        static constexpr int h_comp2=(rank==2)?nComponents2:(rank==1)?nDim:nDim;
        static constexpr int h_comp3=(rank==2)?nDim:(rank==1)?nDim:1;

        /** Existing compatibility Hessian tensor type; full matrix-valued Hessians use hessianEigen(). */
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<h_comp1,h_comp2,h_comp3>> h_type;
        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> mean_value_type;
        typedef boost::multi_array<id_type,2> functionvalue_type;
        typedef boost::multi_array<g_type,2> grad_type;
        typedef boost::multi_array<h_type,2> hessian_type;
        using derivative_matrix_type = eigen_matrix_xx_col_type<value_type>;
        using gradient_component_matrices_type =
            std::array<derivative_matrix_type, nDim>;
        using hessian_component_matrices_type =
            std::array<derivative_matrix_type, nDim * nDim>;

        /** @brief Construct empty reference precomputation state. */
        PreCompute() = default;

        /**
         * precompute values of basis functions given by __ref_ele at
         * points \c __pts
         */
        PreCompute( reference_element_ptrtype const& __ref_ele,
                    matrix_node_t_type const& __pts, Evaluation request = Evaluation::All )
            : M_ref_ele( __ref_ele ),
              M_nodes( __pts ),
              M_mean( __ref_ele->nbDof() ),
              M_reference( RuntimeOrder{ ::Feel::detail::polynomialOrder( *__ref_ele ) } )
        {
            // First moments depend on the reference element, independently of sample points.
            auto fm = __ref_ele->firstMoments();
            for( int i = 0; i < __ref_ele->nbDof(); ++i )
            {
                M_mean[i] = fm(i,0);
            }

            // Populate requested native tables and existing FE compatibility tensors.
            initialize( request );
        }

        /** copy constructor (deep copy) */
        PreCompute( PreCompute const& __pc ) = default;

        /** */
        ~PreCompute() = default;

        /** copy operator (deep copy) */
        PreCompute& operator=( PreCompute const& __pc ) = default;

        /** @brief Reuse owners and recurrence storage with the current request. */
        void update( matrix_node_t_type const& __pts )
        {
            update( __pts, M_reference.request() );
        }

        /** @brief Change evaluation points and requested groups; unrequested access throws. */
        void update( matrix_node_t_type const& points, Evaluation request )
        {
            if ( points.size1() != nDim || !::Feel::detail::validEvaluation( request ) )
                throw std::invalid_argument( "invalid precompute coordinates or request" );
            M_nodes = points;
            initialize( request );
        }

        /** @brief Active reference-table request. */
        [[nodiscard]] Evaluation requestedTables() const noexcept { return M_reference.request(); }

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

        //!
        //! @return the first moment values of the basis functions
        //!
        mean_value_type const& firstMoments() const
            {
                return M_mean;
            }
        //!
        //! @return the first moment value of the i-th basis function
        //!
        value_type  firstMoment( int i ) const
            {
                return M_mean[i];
            }
        /**
         *  Return the matrix evaluation the basis functions (rows) at
         *  a set of points (columns). The matrix is column oriented,
         *  so performance wise it is better to iterate over the columns
         *
         */
        functionvalue_type const& phi() const
        {
            M_reference.require( Evaluation::Values );
            return M_phi;
        }

        /** @brief Borrow the requested legacy value collection until update or destruction. */
        functionvalue_type const* phiPtr() const
        {
            M_reference.require( Evaluation::Values );
            return &M_phi;
        }

        /** @brief Borrow canonical values with interleaved component rows and sample columns. */
        [[nodiscard]] eigen_matrix_xx_col_type<value_type> const& phiEigen() const
        {
            M_reference.require( Evaluation::Values );
            return M_reference.tables().values().matrixStorage();
        }
        /**
         * Returns the value of the q-th node of the i-th basis
         * functions.
         */
        value_type phi( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return phi( i, c1, c2, q, rank_t<rank>() );
        }
        /** @brief Read one requested scalar reference value. */
        value_type phi( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, rank_t<0> ) const
        {
            M_reference.require( Evaluation::Values );
            return *M_phi[i][q].data();
        }
        /** @brief Read one requested vector-component reference value. */
        value_type phi( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, rank_t<1> ) const
        {
            M_reference.require( Evaluation::Values );
            return *( M_phi[i][q].data() + c1 );
        }
        /** @brief Read one requested matrix-component reference value. */
        value_type phi( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, rank_t<2> ) const
        {
            M_reference.require( Evaluation::Values );
            return M_phi[i][q]( c1, c2 );
        }
        /** @brief Borrow one requested legacy value until update or destruction. */
        id_type const& phi( uint16_type i, uint16_type q ) const
        {
            M_reference.require( Evaluation::Values );
            DCHECK( q < nComputedNodes() ) << "invalid node index q=" << q << " nComputedNodes=" << nComputedNodes();
            DCHECK( i < M_ref_ele->nbDof() ) << "invalid dof index i=" << i << " ndof =" << M_ref_ele->nbDof();

            return M_phi[i][q];
        }
        /** @brief Borrow mutable legacy value storage using the existing compatibility API. */
        id_type& phi( uint16_type i, uint16_type q )
        {
            M_reference.require( Evaluation::Values );
            DCHECK( q < nComputedNodes() ) << "invalid node index q=" << q << " nComputedNodes=" << nComputedNodes();
            DCHECK( i < M_ref_ele->nbDof() ) << "invalid dof index i=" << i << " ndof =" << M_ref_ele->nbDof();

            return M_phi[i][q];
        }

        /** @brief Borrow the requested legacy reference-gradient collection. */
        grad_type const& grad() const
        {
            M_reference.require( Evaluation::Gradients );
            return M_grad;
        }
        /** @brief Borrow the requested legacy reference-gradient collection by pointer. */
        grad_type const* gradPtr() const
        {
            M_reference.require( Evaluation::Gradients );
            return &M_grad;
        }

        /**
         * @brief Return one reference-gradient component for all basis
         * functions and evaluation points.
         * @param direction reference-coordinate direction
         * @return column-major matrix with shape `(nDof,nPoints)`
         * @throws std::out_of_range if @p direction is invalid
         */
        [[nodiscard]] derivative_matrix_type const& gradientEigen( uint16_type direction ) const
        {
            if ( direction >= nDim )
                throw std::out_of_range(
                    "PolynomialSet gradient direction is outside the reference dimension" );
            M_reference.require( Evaluation::Gradients );
            return M_reference.tables().gradient( direction ).matrixStorage();
        }

        /**
         * @brief Return one reference-Hessian component for all basis
         * functions and evaluation points.
         * @param firstDirection first reference-coordinate direction
         * @param secondDirection second reference-coordinate direction
         * @return column-major matrix with shape `(nDof,nPoints)`
         * @throws std::out_of_range if either direction is invalid
         */
        [[nodiscard]] derivative_matrix_type const& hessianEigen(
            uint16_type firstDirection, uint16_type secondDirection ) const
        {
            if ( firstDirection >= nDim || secondDirection >= nDim )
                throw std::out_of_range(
                    "PolynomialSet Hessian direction is outside the reference dimension" );
            M_reference.require( Evaluation::Hessians );
            return M_reference.tables().hessian( firstDirection, secondDirection ).matrixStorage();
        }

        /** @brief Read one requested scalar/vector reference-gradient entry. */
        value_type grad( size_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            M_reference.require( Evaluation::Gradients );
            return M_grad[i][q]( c1, c2, 0 );
        }

        /** @brief Read one requested matrix-component reference-gradient entry. */
        value_type grad( size_type i, uint16_type c1, uint16_type c2, uint16_type c3, uint16_type q ) const
        {
            M_reference.require( Evaluation::Gradients );
            return M_grad[i][q]( c1, c2, c3 );
        }

        /** @brief Borrow the requested legacy reference-Hessian collection. */
        hessian_type const& hessian() const
        {
            M_reference.require( Evaluation::Hessians );
            return M_hessian;
        }
        /** @brief Borrow one requested legacy Hessian tensor. */
        h_type const& hessian( uint16_type i, uint16_type q ) const
        {
            M_reference.require( Evaluation::Hessians );
            return M_hessian[i][q];
        }
        /** @brief Read one requested scalar reference-Hessian entry. */
        value_type hessian( size_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            M_reference.require( Evaluation::Hessians );
            return M_hessian[i][q]( c1, c2, 0 );
        }
        /** @brief Borrow one requested legacy Hessian using a dof-sized index. */
        h_type const& hessian( size_type i, uint16_type q ) const
        {
            M_reference.require( Evaluation::Hessians );
            return M_hessian[i][q];
        }
        /** @brief Read one requested vector-component reference-Hessian entry. */
        value_type hessian( size_type i, uint16_type c1, uint16_type c2, uint16_type c3, uint16_type q ) const
        {
            M_reference.require( Evaluation::Hessians );
            return M_hessian[i][q]( c1, c2, c3 );
        }

        /** @brief Reference Hessian of one matrix-valued function component. */
        value_type hessian( size_type function, uint16_type c1, uint16_type c2,
                            uint16_type first, uint16_type second, uint16_type point ) const
            requires( rank == 2 )
        {
            M_reference.require( Evaluation::Hessians );
            return hessianEigen( first, second )( function * nComponents + c1 * nComponents2 + c2, point );
        }

        /** @brief Log requested legacy values using their function/sample shape. */
        void print()
        {
            M_reference.require( Evaluation::Values );
            typedef typename grad_type::index index;
            const index I = M_phi.shape()[0];
            const index Q = M_nodes.size2();

            LOG( INFO ) << "Precompute nodes : " << M_nodes;
            for ( index i = 0; i < I; ++i )
            {
                for ( index q = 0; q < Q; ++q )
                    LOG( INFO ) << "Precompute phi[" << i << ", " << q << "]= " << M_phi[i][q]( 0, 0 );
            }
        }

      private:
        /** @brief Compute selected native tables once and fill compatibility tensors. */
        void initialize( Evaluation request )
        {
            // Preserve the scalar, vector and matrix compatibility layout contracts.
            static_assert( rank <= 2, "legacy precompute tensors support scalar, vector and matrix fields" );
            static_assert( nComponents == nComponents1 * nComponents2 );

            if ( !M_ref_ele || M_nodes.size1() != nDim )
            {
                throw std::invalid_argument( "invalid reference precompute source or coordinate dimension" );
            }

            // Evaluate the requested mathematical tables once through the common tabulator.
            M_reference.compute( *M_ref_ele, M_nodes, request );
            auto const functions = M_ref_ele->nbDof() / nComponents;
            auto const points = M_nodes.size2();

            // Keep compatibility allocations stable when the request and shape are unchanged.
            auto resize = [&]( auto& array, Evaluation group )
            {
                auto const rows = includes( request, group ) ? functions : 0;
                auto const cols = includes( request, group ) ? points : 0;

                if ( array.shape()[0] != rows || array.shape()[1] != cols )
                {
                    array.resize( boost::extents[rows][cols] );
                }
            };

            resize( M_phi, Evaluation::Values );
            resize( M_grad, Evaluation::Gradients );
            resize( M_hessian, Evaluation::Hessians );

            // Keep the legacy rank-three placeholder defined for matrix-valued fields.
            if constexpr ( rank == 2 )
            {
                // The retained rank-three compatibility type has no second
                // derivative axis for a matrix field. Complete derivatives live
                // in canonical hessianEigen(i,j) tables and the six-index accessor.
                for ( std::size_t entry = 0; entry < M_hessian.num_elements(); ++entry )
                {
                    M_hessian.data()[entry].setZero();
                }
            }

            // Dispatch once, then pack all field ranks with the same component-row rule.
            auto packSelected = [&]( auto selection )
            {
                constexpr auto selected = decltype( selection )::value;

                for ( std::size_t function = 0; function < functions; ++function )
                {
                    for ( std::size_t point = 0; point < points; ++point )
                    {
                        for ( int component = 0; component < nComponents; ++component )
                        {
                            int const c1 = component / nComponents2;
                            int const c2 = component % nComponents2;
                            auto const row = function * nComponents + component;

                            // Pack values into the existing function/component representation.
                            if constexpr ( includes( selected, Evaluation::Values ) )
                            {
                                M_phi[function][point]( c1, c2 ) = M_reference.tables().values().matrix()( row, point );
                            }

                            // Place the reference direction on the legacy rank-dependent axis.
                            if constexpr ( includes( selected, Evaluation::Gradients ) )
                            {
                                for ( int direction = 0; direction < nDim; ++direction )
                                {
                                    auto const value = M_reference.tables().gradient( direction ).matrix()( row, point );

                                    if constexpr ( rank == 0 )
                                    {
                                        M_grad[function][point]( 0, direction, 0 ) = value;
                                    }
                                    else if constexpr ( rank == 1 )
                                    {
                                        M_grad[function][point]( c1, direction, 0 ) = value;
                                    }
                                    else
                                    {
                                        M_grad[function][point]( c1, c2, direction ) = value;
                                    }
                                }
                            }

                            // Scalar/vector tensors retain both ordered derivative axes.
                            if constexpr ( includes( selected, Evaluation::Hessians ) && rank < 2 )
                            {
                                for ( int first = 0; first < nDim; ++first )
                                {
                                    for ( int second = 0; second < nDim; ++second )
                                    {
                                        auto const value = M_reference.tables().hessian( first, second ).matrix()( row, point );

                                        if constexpr ( rank == 0 )
                                        {
                                            M_hessian[function][point]( first, second, 0 ) = value;
                                        }
                                        else
                                        {
                                            M_hessian[function][point]( c1, first, second ) = value;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } };

            ::Feel::withEvaluationRequest( request, packSelected );
        }

      /** Shared immutable source of modal coefficients and effective degree. */
      reference_element_ptrtype M_ref_ele;

      /** Owned reference coordinates, independent of caller lifetime. */
      matrix_node_t_type M_nodes;

      /** Existing first-moment contract, independent of the point batch. */
      mean_value_type M_mean;

      /** Existing function/component tensors consumed by FE contexts. */
      functionvalue_type M_phi;

      /** Existing function/component/direction gradient tensors. */
      grad_type M_grad;

      /** Existing Hessian tensors, including both matrix-component axes. */
      hessian_type M_hessian;

      /** Sole native owners and reusable modal recurrence state. */
      ::Feel::detail::PrecomputeTabulationState<reference_element_type> M_reference;
    }; /** class PreCompute **/

    typedef PreCompute precompute_type;
    typedef std::shared_ptr<precompute_type> precompute_ptrtype;

    precompute_ptrtype
    preCompute( self_ptrtype p, points_type const& P )
    {
        return preCompute( p, P, Evaluation::All );
    }

    /** @brief Construct only the requested reference tables through the existing factory. */
    precompute_ptrtype preCompute( self_ptrtype p, points_type const& points, Evaluation request )
    {
        return std::make_shared<PreCompute>( p, points, request );
    }

    precompute_ptrtype
    preCompute( points_type const& P )
    {
        return std::make_shared<PreCompute>( super_enable_shared_from_this::shared_from_this(), P );
    }

    /** @brief Select reference derivative groups while retaining shared source ownership. */
    precompute_ptrtype preCompute( points_type const& points, Evaluation request )
    {
        return preCompute( super_enable_shared_from_this::shared_from_this(), points, request );
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
                geopc[__f][__p] = precompute_ptrtype(  new precompute_type( p, P ) );
            }
        }

        return geopc;
    }

    /** @brief Precompute caller-supplied reference face points for every orientation.
     * A generator accepting (face,permutation) supplies oriented coordinates;
     * the existing one-argument form supplies permutation-independent coordinates.
     */
    template <typename FacePointGenT>
    std::vector<std::map<typename convex_type::permutation_type, precompute_ptrtype> >
    preComputeOnFaces( self_ptrtype p, FacePointGenT const& P )
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
                if constexpr ( requires { P( __f, __p ); } )
                    geopc[__f][__p] = std::make_shared<precompute_type>( p, P( __f, __p ) );
                else
                    geopc[__f][__p] = std::make_shared<precompute_type>( p, P( __f ) );
            }
        }

        return geopc;
    }

    std::vector<std::map<typename convex_type::permutation_type, precompute_ptrtype> >
    preComputeOnFaces( self_ptrtype p, std::vector<std::map<uint16_type,points_type > > const& P )
    {
        typedef typename convex_type::permutation_type permutation_type;
        std::vector<std::map<permutation_type, precompute_ptrtype> > geopc( convex_type::numTopologicalFaces );
        CHECK( P.size() == convex_type::numTopologicalFaces ) << "invalid number of face : " << P.size() << " vs " << convex_type::numTopologicalFaces;
        for ( uint16_type __f = 0; __f < convex_type::numTopologicalFaces; ++__f )
        {
            auto const& ptsOnFace = P[__f];
            for ( permutation_type __p( permutation_type::IDENTITY );
                    __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                auto itFindPerm = ptsOnFace.find( __p.value() );
                CHECK( itFindPerm != ptsOnFace.end() ) << "permutation not find";
                geopc[__f][__p] = precompute_ptrtype(  new precompute_type( p, itFindPerm->second ) );
            }
        }

        return geopc;
    }

    /** @brief Precompute oriented edge samples with the established result-map type.
     * A two-argument generator supplies each orientation explicitly. For the
     * legacy ordered one-argument point batch, reversal uses the same column
     * permutation as QuadMapped. The legacy undefined key aliases identity.
     */
    template<typename PAtEdges>
    std::vector<std::map<typename convex_type::permutation_type, precompute_ptrtype>>
    preComputeOnEdges( self_ptrtype gm, PAtEdges const& p )
    {
        using Permutation = typename convex_type::permutation_type;
        using EdgePermutation = typename convex_type::edge_permutation_type;
        std::vector<std::map<Permutation, precompute_ptrtype>> result( convex_type::numEdges );

        // Build each oriented reference batch using the caller's established map type.
        for ( uint16_type edge = 0; edge < convex_type::numEdges; ++edge )
        {
            for ( EdgePermutation orientation{ EdgePermutation::IDENTITY };
                  orientation < EdgePermutation{ EdgePermutation::N_PERMUTATIONS }; ++orientation )
            {
                auto const key = Permutation{ orientation.value() };

                // Explicit generators own their orientation; legacy batches reverse columns.
                if constexpr ( requires { p( edge, orientation ); } )
                    result[edge][key] = std::make_shared<precompute_type>( gm, p( edge, orientation ) );
                else
                {
                    points_type points( p( edge ) );
                    if ( orientation.value() != EdgePermutation::IDENTITY ) eigenMap( points ).rowwise().reverseInPlace();
                    result[edge][key] = std::make_shared<precompute_type>( gm, points );
                }
            }

            // Preserve the historical undefined-key lookup as an identity alias.
            result[edge][Permutation{ Permutation::NO_PERMUTATION }] = result[edge].at( Permutation{ EdgePermutation::IDENTITY } );
        }
        return result;
    }

    template<size_type context_v, typename Basis_t, typename Geo_t, typename ElementType, size_type context_g = context_v, int SubEntityCoDim = 0>
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
        static constexpr bool is_product = Basis_t::is_product;
        static constexpr bool is_scalar = Basis_t::is_scalar;
        static constexpr bool is_vectorial = Basis_t::is_vectorial;
        static constexpr bool is_tensor2 = Basis_t::is_tensor2;
        static constexpr bool is_tensor2symm = Basis_t::is_tensor2symm;
        static const uint16_type nComponents = Basis_t::nComponents;
        static const uint16_type nComponents1 = Basis_t::nComponents1;
        static const uint16_type nComponents2 = Basis_t::nComponents2;

        static inline const bool is_hdiv_conforming = Feel::is_hdiv_conforming<Basis_t>::value;
        static inline const bool is_hcurl_conforming = Feel::is_hcurl_conforming<Basis_t>::value;

        static const bool do_optimization_p1 =
            !Basis_t::is_order_dynamic && !Geo_t::is_order_dynamic &&
            ( nOrder <= 1 ) && ( Geo_t::nOrder == 1 ) && ( convex_type::is_simplex );
        using do_optimization_p1_t = bool_t<do_optimization_p1>;
        using optimization_p1_t = bool_t<true>;
        using no_optimization_p1_t = bool_t<false>;


        static const bool do_optimization_p2 =
            !Basis_t::is_order_dynamic && !Geo_t::is_order_dynamic &&
            ( nOrder == 2 ) && ( Geo_t::nOrder == 1 ) && ( convex_type::is_simplex );
        using do_optimization_p2_t = bool_t<do_optimization_p2>;
        using optimization_p2_t = bool_t<true>;
        using no_optimization_p2_t = bool_t<false>;

        static const bool do_pt_q1 =
            !Basis_t::is_order_dynamic && !Geo_t::is_order_dynamic &&
            ( nOrder <= 2 ) && ( Geo_t::nOrder == 1 ) && ( convex_type::is_hypercube );


        typedef typename Basis_t::polyset_type polyset_type;
        static const uint16_type rank = polyset_type::rank;


        typedef Basis_t reference_element_type;
        typedef std::shared_ptr<Basis_t> reference_element_ptrtype;

        typedef typename reference_element_type::value_type value_type;

        typedef ElementType geometric_element_type;
        typedef typename Geo_t::template Context<ElementType,SubEntityCoDim> geometric_mapping_context_type;
        typedef std::shared_ptr<geometric_mapping_context_type> geometric_mapping_context_ptrtype;

        static inline const bool second_derivative_require_grad = !geometric_mapping_context_type::is_linear;//  Geo_t::nOrder > 1 || !convex_type::is_simplex;

        typedef typename node<value_type>::type node_type;

        typedef ublas::matrix<value_type> matrix_type;
        typedef points_type matrix_node_t_type;

#if 0
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nDim> ref_grad_type;
        typedef Eigen::Matrix<value_type,nComponents1,NDim> grad_type;
        typedef Eigen::Matrix<value_type,NDim,NDim,1> hess_type;
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> laplacian_type;
        typedef Eigen::Matrix<value_type,nComponents2,1> div_type;
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> dn_type;
        typedef Eigen::Matrix<value_type,3,1> curl_type;
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> dx_type;
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> dy_type;
        typedef Eigen::Matrix<value_type,nComponents1,nComponents2> dz_type;
#else
        using id_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>>;
        using ref_grad_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2,NDim>>;
        static constexpr int g_comp2=(rank==2)?nComponents2:(rank==1)?NDim:NDim;
        static constexpr int g_comp3=(rank==2)?NDim:(rank==1)?1:1;
        using grad_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,g_comp2,g_comp3>>;
        using symm_grad_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,NDim>>;
        static constexpr int h_comp1=(rank==2)?nComponents1:(rank==1)?nComponents1:NDim;
        static constexpr int h_comp2=(rank==2)?nComponents2:(rank==1)?NDim:NDim;
        static constexpr int h_comp3=(rank==2)?NDim:(rank==1)?NDim:1;
        using hess_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<h_comp1,h_comp2,h_comp3>>;
        using laplacian_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>>;
        using div_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents2,1>>;
        using dn_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>>;
        using curl_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<(nComponents1==3)?nComponents1:1,1>>;
        using dx_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents2,1>>;
        using dy_type = dx_type;
        using dz_type = dx_type;
        using normal_component_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<(rank==2)? nComponents1:1,1>>;
#endif
        typedef geometric_mapping_context_type gmc_type;
        typedef Eigen::Matrix<value_type,Eigen::Dynamic, Eigen::Dynamic> matrix_eigen_type;
        typedef Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::PDim> matrix_eigen_NP_type;
        typedef Eigen::Matrix<value_type,gmc_type::PDim,gmc_type::NDim> matrix_eigen_PN_type;
        typedef Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::NDim> matrix_eigen_NN_type;
        typedef Eigen::Matrix<value_type,nComponents1,NDim> matrix_eigen_grad_type;
        typedef typename Eigen::Map<const Eigen::Matrix<value_type,gmc_type::NDim,gmc_type::PDim,Eigen::ColMajor> > matrix_eigen_ublas_NP_type;
        typedef typename Eigen::Map<const Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > matrix_eigen_ublas_type;
        using quadrature_batch_table_type =
            Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using quadrature_batch_row_type = typename quadrature_batch_table_type::ConstRowXpr;
        using quadrature_batch_component_table_type =
            Eigen::Map<quadrature_batch_table_type const, Eigen::Unaligned,
                       Eigen::OuterStride<>>;

        Context( reference_element_ptrtype const& __RefEle,
                 geometric_mapping_context_ptrtype const& __gmc,
                 precompute_ptrtype const& __pc )
            : M_pc( __pc ),
              M_npoints( 0 /*__pc->nPoints()*/ ),

              M_ipt( 0 ),
              M_ref_ele( __RefEle ),

              M_gmc( __gmc ),
              M_dofs( nDofs() ),
              M_phi( __pc->phi() ),
              M_gradphi( includes( __pc->requestedTables(), Evaluation::Gradients ) ? __pc->gradPtr() : nullptr ),
              M_normal_component(),
              M_trace(),
              M_dn(),
              M_grad(),
              M_symm_grad(),
              M_dx(),
              M_dy(),
              M_dz(),
              M_id_batch(),
              M_grad_batch(),
              M_div_batch(),
              M_curl_batch()
        {
            std::iota( M_dofs.begin(), M_dofs.end(), 0 );

            update( __gmc, __pc );
        }
        Context( Context const& c )
            :
            M_pc( c.M_pc ),
            M_npoints( c.M_npoints ),

            M_ipt( c.M_ipt ),
            M_ref_ele( c.M_ref_ele ),

            M_gmc( c.M_gmc ),

            M_dofs( c.M_dofs ),
            M_phi( c.M_phi ),
            M_gradphi( c.M_gradphi ),
            M_normal_component( c.M_normal_component ),
            M_trace( c.M_trace ),
            M_dn( c.M_dn ),
            M_grad( c.M_grad ),
            M_symm_grad( c.M_symm_grad ),
            M_dx( c.M_dx ),
            M_dy( c.M_dy ),
            M_dz( c.M_dz ),
            M_id_batch( c.M_id_batch ),
            M_grad_batch( c.M_grad_batch ),
            M_div_batch( c.M_div_batch ),
            M_curl_batch( c.M_curl_batch ),
            M_id_batch_valid( c.M_id_batch_valid ),
            M_grad_batch_valid( c.M_grad_batch_valid ),
            M_div_batch_valid( c.M_div_batch_valid ),
            M_curl_batch_valid( c.M_curl_batch_valid )
            {

            }
        virtual ~Context() {}
        /**
         * if isTransformationEquivalent is set to true in basis then no
         * transformation is required
         */
        void transformationEquivalence( geometric_mapping_context_ptrtype const& __gmc,
                                        bool_t<true> )
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
                                        bool_t<false> )
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
            validateReferenceGradients( __gmc->dynamicContext() );
            this->invalidateQuadratureBatchTables();
            //M_phi = M_pc->get()->phi();
            //M_gradphi = M_pc->get()->grad();
            transformationEquivalence( __gmc, bool_t<Basis_t::isTransformationEquivalent>() );
#if 0

            for ( int i = 0; i < M_gradphi.num_elements(); ++i )
                std::cout << "M_gradphi[" << i << "]=" << *( M_gradphi.data()+i ) << "\n";

#endif
            update( __gmc, rank_t<rank>() );
            //update( __gmc, int_t<rank>(), bool_t<false>() );
        }

        //
        // update scalar valued polynomials
        //
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, rank_t<0> );
        //
        // update vector valued polynomials
        //
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, rank_t<1> );
        //
        // update matrix valued polynomials
        //
        FEELPP_DONT_INLINE void update( geometric_mapping_context_ptrtype const& __gmc, rank_t<2>  );

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

        /**
         *  set the geometric mapping context
         */
        void setGmContext( geometric_mapping_context_ptrtype const& gmc )
        {
            M_gmc = gmc;
        }

        //! \return the element id
        size_type eId() const
        {
            return M_gmc->id();
        }

        //! @return the dynamic context associated
        size_type dynamicContext() const { return M_gmc->dynamicContext()/* M_dynamic_context*/; }

        [[nodiscard]] uint16_type basisOrder() const
        {
            if constexpr ( requires { Basis_t::is_order_dynamic; } )
            {
                if constexpr ( Basis_t::is_order_dynamic )
                {
                    if ( M_ref_ele )
                    {
                        if constexpr ( requires( reference_element_type const& fe ) { fe.polynomialDegree(); } )
                            return static_cast<uint16_type>( M_ref_ele->polynomialDegree() );
                        else if constexpr ( requires( reference_element_type const& fe ) { fe.order(); } )
                            return static_cast<uint16_type>( M_ref_ele->order() );
                    }
                }
            }
            return static_cast<uint16_type>( nOrder );
        }

        /**
         * @brief Report whether the active geometric mapping is affine.
         *
         * Static mappings preserve the compile-time fast path. Dynamic-order
         * mappings query their runtime strategy, preventing a higher-order
         * geometry from inheriting the placeholder compile-time order one.
         */
        [[nodiscard]] bool geometryIsLinear() const
        {
            if constexpr ( Geo_t::is_order_dynamic )
            {
                return M_gmc && M_gmc->geometricMapping()->isLinear();
            }
            else
            {
                return Geo_t::nOrder == 1;
            }
        }

        /**
         * @brief Report whether first derivatives may be cached at one point.
         */
        [[nodiscard]] bool canUseP1Optimization() const
        {
            if constexpr ( do_optimization_p1 )
            {
                return true;
            }
            else if constexpr ( convex_type::is_simplex )
            {
                return this->geometryIsLinear() &&
                       this->basisOrder() <= 1;
            }
            return false;
        }

        /**
         * @brief Report whether second derivatives may be cached at one point.
         */
        [[nodiscard]] bool canUseP2Optimization() const
        {
            if constexpr ( do_optimization_p2 )
            {
                return true;
            }
            else if constexpr ( convex_type::is_simplex )
            {
                return this->geometryIsLinear() &&
                       this->basisOrder() == 2;
            }
            return false;
        }

        [[nodiscard]] uint16_type nPointsFirstDerivative() const
        {
            return this->canUseP1Optimization() ? 1 : M_npoints;
        }

        [[nodiscard]] uint16_type nPointsSecondDerivative() const
        {
            return this->canUseP2Optimization() ? 1 : M_npoints;
        }


        /**
         * @return the number of basis functions
         *
         * Uses the reference element's semantic local cardinality when the
         * object is available and the compile-time value during setup.
         */
        uint16_type nDofs() const
        {
            const uint16_type localDof = M_ref_ele ? M_ref_ele->localDof() : nDof;
            constexpr uint16_type factor = ( reference_element_type::is_product ? nComponents1 * nComponents2 : 1 );
            return static_cast<uint16_type>( localDof * factor );
        }

        /**
         * @return the degrees of freedom
         */
        std::vector<uint16_type> const& dofs() const { return M_dofs; }

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
            return id( i, c1, c2, q, int_t<rank>() );
        }

        value_type const& id( uint32_type i,
                              uint16_type /*c1*/,
                              uint16_type /*c2*/,
                              uint32_type q,
                              rank_t<0> ) const
        {
            return *M_phi[i][q].data();
        }

        value_type const& id( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q,
                              rank_t<1> ) const
        {
            Feel::detail::ignore_unused_variable_warning( c2 );
            return *(M_phi[i][q].data()+c1);
        }

        id_type const& id( uint32_type i, uint32_type q ) const
        {
            return M_phi[i][q];
        }

        value_type const& id( uint32_type i,
                              uint16_type c1,
                              uint16_type c2,
                              uint32_type q,
                              rank_t<2> ) const
        {
            return M_phi[i][q]( c1,c2 );
        }
        value_type const& normalComponent( uint32_type i,
                                           uint16_type c1,
                                           uint16_type c2,
                                           uint32_type q ) const
            {
                return M_normal_component[i][q]( c1,c2 );
            }
        normal_component_type const& normalComponent( uint32_type i,
                                                      uint32_type q ) const
            {
                return M_normal_component[i][q];
            }
        value_type const& trace( uint32_type i,
                                 uint16_type c1,
                                 uint16_type c2,
                                 uint32_type q ) const
            {
                return M_trace[i][q];
            }
        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return d( i, c1, c2, q, optimization_p1_t() );
                return d( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0]( c1,c2,0 );
            }
        value_type d( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false>  ) const
            {
                return M_grad[i][q]( c1,c2,0 );
            }

        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dx( i, c1, c2, q, optimization_p1_t() );
                return dx( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][0]( c1,0,0 );
            }
        value_type dx( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][q]( c1,0,0 );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dy( i, c1, c2, q, optimization_p1_t() );
                return dy( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][0]( c1,1,0 );
            }
        value_type dy( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][q]( c1,1,0 );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dz( i, c1, c2, q, optimization_p1_t() );
                return dz( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][0]( c1,2,0 );
            }
        value_type dz( uint32_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                static_assert( nDim >= 1, "INVALID_DIM" );
                return M_grad[i][q]( c1,2,0 );
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
                if ( this->canUseP1Optimization() )
                    return grad( i, q, optimization_p1_t() );
                return grad( i, q, no_optimization_p1_t() );
            }
        grad_type const& grad( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0];
            }
        grad_type const& grad( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_grad[i][q];
            }

        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return grad( i, c1, c2, q, optimization_p1_t() );
                return grad( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0]( c1,c2,0 );
            }
        value_type grad( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                return M_grad[i][q]( c1,c2,0 );
            }


        symm_grad_type const& symmetricGradient( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return symmetricGradient( i, q, optimization_p1_t() );
                return symmetricGradient( i, q, no_optimization_p1_t() );
            }
        symm_grad_type const& symmetricGradient( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_symm_grad[i][0];
            }
        symm_grad_type const& symmetricGradient( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_symm_grad[i][q];
            }

        value_type symmetricGradient( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return symmetricGradient( i, c1, c2, q, optimization_p1_t() );
                return symmetricGradient( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type symmetricGradient( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                return M_symm_grad[i][0]( c1,c2 );
            }
        value_type symmetricGradient( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                return M_symm_grad[i][q]( c1,c2 );
            }
        
        dx_type  dx( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dx( i, q, optimization_p1_t() );
                return dx( i, q, no_optimization_p1_t() );
            }
        dx_type  dx( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0].chip(0,1);
            }
        dx_type  dx( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_grad[i][q].chip(0,1);
            }
        dy_type  dy( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dy( i, q, optimization_p1_t() );
                return dy( i, q, no_optimization_p1_t() );
            }
        dy_type  dy( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0].chip(1,1);
            }
        dy_type  dy( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_grad[i][q].chip(1,1);
            }
        dz_type  dz( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return dz( i, q, optimization_p1_t() );
                return dz( i, q, no_optimization_p1_t() );
            }
        dz_type  dz( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_grad[i][0].chip(2,1);
            }
        dz_type  dz( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_grad[i][q].chip(2,1);
            }

        /**
         * divergence of the basis function at the q-th point.
         *
         * \param q index of the points to evaluate the divergence
         * \param i index containing current function and component indices
         * \param c1 first component index
         * \param c2 second component index
         * \return divergence of the basis function at the q-th point
         */

        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
             return div( i, c1, c2, q, rank_t<rank>() );
        }

        value_type div( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, rank_t<0> ) const
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

        div_type const& div( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return div( i, q, optimization_p1_t() );
                return div( i, q, no_optimization_p1_t() );
            }
        div_type const& div( uint16_type i, uint32_type q, bool_t<true> ) const
            {
                return M_div[i][0];
            }
        div_type const& div( uint16_type i, uint32_type q, bool_t<false> ) const
            {
                return M_div[i][q];
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<1> ) const
            {
                if ( this->canUseP1Optimization() )
                    return div( i, c1, c2, q, optimization_p1_t() );
                return div( i, c1, c2, q, no_optimization_p1_t() );
            }

        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][0]( 0,0 );
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][q]( 0,0 );
            }

        /**
         * divergence of a matrix field
         */
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<2> ) const
            {
                if ( this->canUseP1Optimization() )
                    return div( i, c1, c2, q, rank_t<2>(), optimization_p1_t() );
                return div( i, c1, c2, q, rank_t<2>(), no_optimization_p1_t() );
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<2>, optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][0]( c1,0 );
            }
        value_type div( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<2>, no_optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_div[i][q]( c1,0 );
            }
        /**
         * curl of the basis function at the q-th point.
         *
         * \param q index of the points to evaluate the curl
         * \param i index containing current function and component indices
         * \return curl of the basis function at the q-th point
         */
        curl_type const& curl( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return curl( i, q, optimization_p1_t() );
                return curl( i, q, no_optimization_p1_t() );
            }
        curl_type const& curl( uint16_type i, uint32_type q, optimization_p1_t ) const
            {
                return M_curl[i][0];
            }
        curl_type const& curl( uint16_type i, uint32_type q, no_optimization_p1_t ) const
            {
                return M_curl[i][q];
            }
        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curl( i, c1, c2, q, int_t<rank>() );
        }

        value_type curl( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, rank_t<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curl operator, polynomial must be vectorial" );
            return 0;
        }

        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<1> ) const
            {
                if ( this->canUseP1Optimization() )
                    return curl( i, c1, c2, q, optimization_p1_t() );
                return curl( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( c1 );
            }
        value_type curl( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, no_optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( c1 );
            }

        value_type curlx( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP1Optimization() )
                    return curlx( i, q, optimization_p1_t() );
                return curlx( i, q, no_optimization_p1_t() );
            }
        value_type curlx( uint16_type i, uint32_type q, optimization_p1_t ) const
            {
                return M_curl[i][0](0);
            }
        value_type curlx( uint16_type i, uint32_type q, no_optimization_p1_t ) const
            {
                return M_curl[i][q](0);
            }

        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curlx( i, c1, c2, q, int_t<rank>() );
        }

        value_type curlx( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, rank_t<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curlx operator, field must be vectorial" );
            return 0;
        }

        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<1> ) const
            {
                if ( this->canUseP1Optimization() )
                    return curlx( i, c1, c2, q, optimization_p1_t() );
                return curlx( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<true> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 0 );
            }
        value_type curlx( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, bool_t<false> ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 0 );
            }

        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curly( i, c1, c2, q, int_t<rank>() );
        }

        value_type curly( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, rank_t<0> ) const
        {
            Feel::detail::ignore_unused_variable_warning( i );
            Feel::detail::ignore_unused_variable_warning( q );
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            throw std::logic_error( "invalid use of curly operator, field must be vectorial" );
            return 0;
        }

        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<1> ) const
            {
                if ( this->canUseP1Optimization() )
                    return curly( i, c1, c2, q, optimization_p1_t() );
                return curly( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 1 );
            }
        value_type curly( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, no_optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 1 );
            }

        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
        {
            Feel::detail::ignore_unused_variable_warning( c1 );
            Feel::detail::ignore_unused_variable_warning( c2 );
            return curlz( i, c1, c2, q, int_t<rank>() );
        }

        value_type curlz( uint16_type i,  uint16_type c1, uint16_type c2, uint32_type q, rank_t<0> ) const
            {
                Feel::detail::ignore_unused_variable_warning( i );
                Feel::detail::ignore_unused_variable_warning( q );
                Feel::detail::ignore_unused_variable_warning( c1 );
                Feel::detail::ignore_unused_variable_warning( c2 );
                throw std::logic_error( "invalid use of curlz operator, field must be vectorial" );
                return 0;
            }

        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<1> ) const
            {
                if ( this->canUseP1Optimization() )
                    return curlz( i, c1, c2, q, optimization_p1_t() );
                return curlz( i, c1, c2, q, no_optimization_p1_t() );
            }
        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][0]( 2 );
            }
        value_type curlz( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, no_optimization_p1_t ) const
            {
                Feel::detail::ignore_unused_variable_warning( c2 );
                return M_curl[i][q]( 2 );
            }

        hess_type const& hess( uint16_type i, uint32_type q ) const
            {
                return M_hessian[i][q];
            }
        value_type hess( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP2Optimization() )
                    return hess( i, c1, c2, q, int_t<rank>(), optimization_p2_t() );
                return hess( i, c1, c2, q, int_t<rank>(), no_optimization_p2_t() );
            }
        value_type hess( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<0>, optimization_p2_t ) const
            {
                return M_hessian[i][0]( c1,c2,0 );
            }
        value_type hess( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, rank_t<0>, no_optimization_p2_t ) const
            {
                return M_hessian[i][q]( c1,c2,0 );
            }
        laplacian_type const& laplacian( uint16_type i, uint32_type q ) const
            {
                if ( this->canUseP2Optimization() )
                    return laplacian( i, q, optimization_p2_t() );
                return laplacian( i, q, no_optimization_p2_t() );
            }
        laplacian_type const& laplacian( uint16_type i, uint32_type q, optimization_p2_t ) const
            {
                return M_laplacian[i][0];
            }
        laplacian_type const& laplacian( uint16_type i, uint32_type q, no_optimization_p2_t ) const
            {
                return M_laplacian[i][q];
            }
        value_type laplacian( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q ) const
            {
                if ( this->canUseP2Optimization() )
                    return laplacian( i, c1, c2, q, optimization_p2_t() );
                return laplacian( i, c1, c2, q, no_optimization_p2_t() );
            }
        value_type laplacian( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, optimization_p2_t ) const
            {
                return M_laplacian[i][0](c1,c2);
            }
        value_type laplacian( uint16_type i, uint16_type c1, uint16_type c2, uint32_type q, no_optimization_p2_t ) const
            {
                return M_laplacian[i][q](c1,c2);
            }

        /**
         * @brief Values of one basis-function component at all quadrature points.
         *
         * The returned row is a non-owning view of lazy component-major Eigen
         * storage. It remains valid until the next context update.
         *
         * @param i local basis-function index
         * @param c1 first value component
         * @param c2 second value component
         * @return row expression indexed by quadrature point
         */
        quadrature_batch_row_type idBatch( uint16_type i, uint16_type c1 = 0,
                                           uint16_type c2 = 0 ) const
        {
            this->ensureIdBatchTable();
            CHECK_LT( i, this->nDofs() );
            CHECK_LT( c1, nComponents1 );
            CHECK_LT( c2, nComponents2 );
            auto const component = static_cast<Eigen::Index>( c1 + nComponents1*c2 );
            return M_id_batch.row( static_cast<Eigen::Index>( i )*
                                   nComponents1*nComponents2 + component );
        }

        /**
         * @brief Complete component-major basis-value table.
         *
         * Rows are ordered by local basis function and value component;
         * columns are quadrature points. The table is constructed lazily and
         * remains valid until the next finite-element context update.
         *
         * @return read-only Eigen table of physical basis values
         */
        quadrature_batch_table_type const& valueBatchTable() const
        {
            this->ensureIdBatchTable();
            return M_id_batch;
        }

        /**
         * @brief One physical value component for every local DoF and point.
         *
         * The returned strided Eigen view has shape `(nDof,nPoint)` and avoids
         * copying interleaved component rows before a local tensor product.
         * It remains valid until the next context update.
         *
         * @param c1 first value component
         * @param c2 second value component
         * @return non-owning DoF-by-point Eigen view
         */
        quadrature_batch_component_table_type
        valueComponentBatchTable( uint16_type c1 = 0,
                                  uint16_type c2 = 0 ) const
        {
            this->ensureIdBatchTable();
            CHECK_LT( c1, nComponents1 );
            CHECK_LT( c2, nComponents2 );
            auto const component = static_cast<Eigen::Index>(
                c1 + nComponents1*c2 );
            return this->componentBatchTable(
                M_id_batch, nComponents1*nComponents2, component );
        }

        /**
         * @brief Physical gradient component at all quadrature points.
         *
         * Affine compact derivative storage is broadcast into the batch table,
         * so callers use the same row API for static/runtime and affine/curved
         * mappings.
         *
         * @param i local basis-function index
         * @param c1 first field component
         * @param c2 second field or derivative component
         * @param c3 third tensor component
         * @return row expression indexed by quadrature point
         */
        quadrature_batch_row_type gradientBatch( uint16_type i, uint16_type c1,
                                                 uint16_type c2,
                                                 uint16_type c3 = 0 ) const
        {
            this->ensureGradientBatchTable();
            CHECK_LT( i, this->nDofs() );
            CHECK_LT( c1, nComponents1 );
            CHECK_LT( c2, g_comp2 );
            CHECK_LT( c3, g_comp3 );
            auto const component = static_cast<Eigen::Index>(
                ( c1*g_comp2 + c2 )*g_comp3 + c3 );
            return M_grad_batch.row( static_cast<Eigen::Index>( i )*
                                     nComponents1*g_comp2*g_comp3 + component );
        }

        /**
         * @brief One physical gradient component for every DoF and point.
         * @param c1 field component
         * @param c2 physical derivative direction
         * @param c3 third tensor component
         * @return non-owning DoF-by-point Eigen view
         */
        quadrature_batch_component_table_type
        gradientComponentBatchTable( uint16_type c1, uint16_type c2,
                                     uint16_type c3 = 0 ) const
        {
            this->ensureGradientBatchTable();
            CHECK_LT( c1, nComponents1 );
            CHECK_LT( c2, g_comp2 );
            CHECK_LT( c3, g_comp3 );
            auto const component = static_cast<Eigen::Index>(
                ( c1*g_comp2 + c2 )*g_comp3 + c3 );
            return this->componentBatchTable(
                M_grad_batch, nComponents1*g_comp2*g_comp3, component );
        }

        /**
         * @brief Divergence component at all quadrature points.
         * @param i local basis-function index
         * @param component divergence component
         * @return row expression indexed by quadrature point
         */
        quadrature_batch_row_type divergenceBatch( uint16_type i,
                                                   uint16_type component = 0 ) const
        {
            this->ensureDivergenceBatchTable();
            CHECK_LT( i, this->nDofs() );
            CHECK_LT( component, nComponents2 );
            return M_div_batch.row( static_cast<Eigen::Index>( i )*nComponents2 +
                                    component );
        }

        /**
         * @brief One divergence component for every local DoF and point.
         * @param component divergence component
         * @return non-owning DoF-by-point Eigen view
         */
        quadrature_batch_component_table_type
        divergenceComponentBatchTable( uint16_type component = 0 ) const
        {
            this->ensureDivergenceBatchTable();
            CHECK_LT( component, nComponents2 );
            return this->componentBatchTable(
                M_div_batch, nComponents2, component );
        }

        /**
         * @brief Curl component at all quadrature points.
         * @param i local basis-function index
         * @param component curl component
         * @return row expression indexed by quadrature point
         */
        quadrature_batch_row_type curlBatch( uint16_type i,
                                             uint16_type component = 0 ) const
        {
            constexpr uint16_type curlComponents = ( nComponents1 == 3 ) ? nComponents1 : 1;
            this->ensureCurlBatchTable();
            CHECK_LT( i, this->nDofs() );
            CHECK_LT( component, curlComponents );
            return M_curl_batch.row( static_cast<Eigen::Index>( i )*curlComponents +
                                     component );
        }

        /**
         * @brief One curl component for every local DoF and point.
         * @param component physical curl component
         * @return non-owning DoF-by-point Eigen view
         */
        quadrature_batch_component_table_type
        curlComponentBatchTable( uint16_type component = 0 ) const
        {
            constexpr uint16_type curlComponents =
                ( nComponents1 == 3 ) ? nComponents1 : 1;
            this->ensureCurlBatchTable();
            CHECK_LT( component, curlComponents );
            return this->componentBatchTable(
                M_curl_batch, curlComponents, component );
        }

        void resizeAndSet( rank_t<0> );
        void resizeAndSet( rank_t<1> );
        void resizeAndSet( rank_t<2> );

        Context() = default;

        void print()
            {
                const uint16_type Q = M_npoints;
                const uint16_type I = this->nDofs();
                for ( uint16_type i = 0; i < I; ++i )
                {
                    for ( uint16_type q = 0; q < Q; ++q )
                    {
                        LOG(INFO) << "M_phi[" << i << "," << q << "]=" << M_phi[i][q];
                    }
                }

            }
    private :
        /** @brief Build a strided DoF-by-point view of one component table. */
        quadrature_batch_component_table_type
        componentBatchTable( quadrature_batch_table_type const& table,
                             Eigen::Index componentCount,
                             Eigen::Index component ) const
        {
            auto const points = static_cast<Eigen::Index>( M_npoints );
            return quadrature_batch_component_table_type(
                table.data() + component*points,
                static_cast<Eigen::Index>( this->nDofs() ), points,
                Eigen::OuterStride<>( componentCount*points ) );
        }

        FEELPP_STRONG_INLINE void updateGrad( geometric_mapping_context_type* thegmc, rank_t<0> );
        FEELPP_STRONG_INLINE void updateFirstDerivativeNormal( geometric_mapping_context_type* thegmc, rank_t<0> );
        FEELPP_STRONG_INLINE void updateHessian( geometric_mapping_context_type* thegmc, rank_t<0> );
        FEELPP_STRONG_INLINE void updateLaplacian( geometric_mapping_context_type* thegmc, rank_t<0> );

        FEELPP_STRONG_INLINE void updateNormalComponent( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateGrad( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateSymm( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateDiv( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateCurl( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateFirstDerivativeNormal( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateHessian( geometric_mapping_context_type* thegmc, rank_t<1> );
        FEELPP_STRONG_INLINE void updateLaplacian( geometric_mapping_context_type* thegmc, rank_t<1> );

        FEELPP_STRONG_INLINE void updateNormalComponent( geometric_mapping_context_type* thegmc, rank_t<2> );
        FEELPP_STRONG_INLINE void updateTrace( geometric_mapping_context_type* thegmc, rank_t<2> );
        FEELPP_STRONG_INLINE void updateGrad( geometric_mapping_context_type* thegmc, rank_t<2> );
        FEELPP_STRONG_INLINE void updateDiv( geometric_mapping_context_type* thegmc, rank_t<2> );

        /** @brief Invalidate lazy quadrature-component tables after an update. */
        void invalidateQuadratureBatchTables() noexcept
        {
            M_id_batch_valid = false;
            M_grad_batch_valid = false;
            M_div_batch_valid = false;
            M_curl_batch_valid = false;
        }

        /** @brief Build the lazy component-major value table. */
        void ensureIdBatchTable() const
        {
            if ( M_id_batch_valid )
                return;
            auto const rows = static_cast<Eigen::Index>( this->nDofs() )*
                              nComponents1*nComponents2;
            M_id_batch.resize( rows, M_npoints );
            for ( uint16_type i = 0; i < this->nDofs(); ++i )
                for ( uint16_type q = 0; q < M_npoints; ++q )
                    for ( uint16_type c2 = 0; c2 < nComponents2; ++c2 )
                        for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                            M_id_batch( static_cast<Eigen::Index>( i )*
                                            nComponents1*nComponents2 +
                                        c1 + nComponents1*c2,
                                        q ) = M_phi[i][q]( c1, c2 );
            M_id_batch_valid = true;
        }

        /** @brief Build the lazy component-major physical-gradient table. */
        void ensureGradientBatchTable() const
        {
            if ( M_grad_batch_valid )
                return;
            CHECK_GT( M_grad.num_elements(), 0 )
                << "gradient batch requested without a gradient context";
            auto const rows = static_cast<Eigen::Index>( this->nDofs() )*
                              nComponents1*g_comp2*g_comp3;
            M_grad_batch.resize( rows, M_npoints );
            for ( uint16_type i = 0; i < this->nDofs(); ++i )
                for ( uint16_type q = 0; q < M_npoints; ++q )
                {
                    auto const storedQ = this->canUseP1Optimization() ? 0 : q;
                    for ( uint16_type c1 = 0; c1 < nComponents1; ++c1 )
                        for ( uint16_type c2 = 0; c2 < g_comp2; ++c2 )
                            for ( uint16_type c3 = 0; c3 < g_comp3; ++c3 )
                                M_grad_batch(
                                    static_cast<Eigen::Index>( i )*
                                            nComponents1*g_comp2*g_comp3 +
                                        ( c1*g_comp2 + c2 )*g_comp3 + c3,
                                    q ) = M_grad[i][storedQ]( c1, c2, c3 );
                }
            M_grad_batch_valid = true;
        }

        /** @brief Build the lazy component-major divergence table. */
        void ensureDivergenceBatchTable() const
        {
            if ( M_div_batch_valid )
                return;
            CHECK_GT( M_div.num_elements(), 0 )
                << "divergence batch requested without a divergence context";
            M_div_batch.resize( static_cast<Eigen::Index>( this->nDofs() )*
                                    nComponents2,
                                M_npoints );
            for ( uint16_type i = 0; i < this->nDofs(); ++i )
                for ( uint16_type q = 0; q < M_npoints; ++q )
                {
                    auto const storedQ = this->canUseP1Optimization() ? 0 : q;
                    for ( uint16_type c = 0; c < nComponents2; ++c )
                        M_div_batch( static_cast<Eigen::Index>( i )*nComponents2 + c,
                                     q ) = M_div[i][storedQ]( c, 0 );
                }
            M_div_batch_valid = true;
        }

        /** @brief Build the lazy component-major curl table. */
        void ensureCurlBatchTable() const
        {
            if ( M_curl_batch_valid )
                return;
            CHECK_GT( M_curl.num_elements(), 0 )
                << "curl batch requested without a curl context";
            constexpr uint16_type curlComponents = ( nComponents1 == 3 ) ? nComponents1 : 1;
            M_curl_batch.resize( static_cast<Eigen::Index>( this->nDofs() )*
                                     curlComponents,
                                 M_npoints );
            for ( uint16_type i = 0; i < this->nDofs(); ++i )
                for ( uint16_type q = 0; q < M_npoints; ++q )
                {
                    auto const storedQ = this->canUseP1Optimization() ? 0 : q;
                    for ( uint16_type c = 0; c < curlComponents; ++c )
                        M_curl_batch( static_cast<Eigen::Index>( i )*curlComponents + c,
                                      q ) = M_curl[i][storedQ]( c, 0 );
                }
            M_curl_batch_valid = true;
        }
    private:
      /** @brief Whether the active evaluator needs reference first derivatives.
       * Curved-map Hessians also require gradients for the existing chain rule.
       */
      [[nodiscard]] bool needsReferenceGradients( size_type dynamicFlags ) const
      {
          // First-order physical operators consume the reference gradients directly.
          bool needed = vm::has_grad_v<context> || vm::has_first_derivative_v<context> ||
                        vm::has_div_v<context> || vm::has_curl_v<context> ||
                        vm::has_symm_v<context> || vm::has_first_derivative_normal_v<context>;

          // Curved-map second derivatives also need gradients for the chain rule.
          if constexpr ( second_derivative_require_grad )
          {
              needed = needed || vm::has_hessian_v<context> ||
                       vm::has_second_derivative_v<context> || vm::has_laplacian_v<context>;
          }

          // Runtime evaluator flags supplement the compile-time operator selection.
          if constexpr ( vm::has_dynamic_basis_function_v<context> )
          {
              auto const flags = dynamicFlags;
              needed = needed || hasGRAD( flags ) || hasFIRST_DERIVATIVE( flags ) ||
                       hasDIV( flags ) || hasCURL( flags ) || hasSYMM( flags ) ||
                       hasFIRST_DERIVATIVE_NORMAL( flags );

              if constexpr ( second_derivative_require_grad )
              {
                  needed = needed || hasHESSIAN( flags ) ||
                           hasSECOND_DERIVATIVE( flags ) || hasLAPLACIAN( flags );
              }
          }

          return needed;
      }

      /** @brief Reject an insufficient request before a physical derivative reads its tables. */
      void validateReferenceGradients( size_type dynamicFlags ) const
      {
          if ( needsReferenceGradients( dynamicFlags ) &&
               ( !M_gradphi || M_gradphi->shape()[0] < this->nDofs() ) )
          {
              throw std::invalid_argument( "reference gradients were not requested for this FE context" );
          }
      }

        boost::optional<precompute_ptrtype> M_pc;
        uint16_type M_npoints;

        uint16_type M_ipt;

        reference_element_ptrtype M_ref_ele;

        geometric_mapping_context_ptrtype M_gmc;
        std::vector<uint16_type> M_dofs;
#if 0
        boost::multi_array<id_type,2> M_phi;
        boost::multi_array<ref_grad_type,2> const* M_gradphi;
#else
        typename precompute_type::functionvalue_type M_phi;
        typename precompute_type::grad_type const* M_gradphi;
#endif
        boost::multi_array<hess_type,2> M_hessphi;
        boost::multi_array<normal_component_type,2> M_normal_component;
        boost::multi_array<value_type,2> M_trace;
        boost::multi_array<dn_type,2> M_dn;
        boost::multi_array<grad_type,2> M_grad;
        boost::multi_array<symm_grad_type,2> M_symm_grad;
        boost::multi_array<dx_type,2> M_dx;
        boost::multi_array<dy_type,2> M_dy;
        boost::multi_array<dz_type,2> M_dz;
        boost::multi_array<div_type,2> M_div;
        boost::multi_array<curl_type,2> M_curl;
        boost::multi_array<hess_type,2> M_hessian;
        boost::multi_array<id_type,2> M_laplacian;
        mutable quadrature_batch_table_type M_id_batch;
        mutable quadrature_batch_table_type M_grad_batch;
        mutable quadrature_batch_table_type M_div_batch;
        mutable quadrature_batch_table_type M_curl_batch;
        mutable bool M_id_batch_valid = false;
        mutable bool M_grad_batch_valid = false;
        mutable bool M_div_batch_valid = false;
        mutable bool M_curl_batch_valid = false;
    };

    template<size_type context_v, typename BasisType, typename GeoType, typename ElementType>
    std::shared_ptr<Context<context_v,BasisType, GeoType, ElementType> >
    context( std::shared_ptr<BasisType> b, std::shared_ptr<GeoType> gm, precompute_ptrtype& pc )
    {
        return std::shared_ptr<Context<context_v,BasisType, GeoType, ElementType> >(
                   new Context<context_v, BasisType, GeoType, ElementType>( context_v,
                           b,
                           gm,
                           pc ) );
    }

    template<int contextv, int contextg, typename BasisType, typename GeoType,typename ElementType>
    std::shared_ptr<Context<contextv,BasisType, GeoType, ElementType> >
    ctx( std::shared_ptr<BasisType> const& b,
         std::shared_ptr<typename GeoType::template Context<contextg, ElementType> > const& gm,
         precompute_ptrtype pc, ElementType& e )
    {
        typedef Context<contextv,BasisType, GeoType, ElementType> ctx_type;
        return std::shared_ptr<ctx_type>( new ctx_type( b, gm, pc ) );

    }
    template<int contextv, typename BasisType, typename GeoType,typename ElementType>
    std::shared_ptr<Context<contextv,BasisType, GeoType, ElementType> >
    ctx( std::shared_ptr<BasisType> const& b,
         std::shared_ptr<typename GeoType::template Context<contextv, ElementType> > const& gm,
         precompute_ptrtype pc, ElementType& e )
    {
        typedef Context<contextv,BasisType, GeoType, ElementType> ctx_type;
        return std::shared_ptr<ctx_type>( new ctx_type( b, gm, pc ) );

    }

public:
    /** @brief Evaluate modal values through the common order policy. */
    template<typename AE>
    matrix_type basisEvaluate( ublas::matrix_expression<AE> const& pts ) const
    {
        return M_order.evaluate( M_basis, pts );
    }

protected:
    /** @brief Consume each C partial_i partial_j Phi expression while its operands are alive.
     * @param points Column-oriented sample coordinates.
     * @param operation Callable consuming direction indices and an ephemeral product expression.
     * @note Direct modal Hessians avoid coefficient operators where available.
     * Native and legacy results share this path and differ only in their destination.
     */
    template<typename Points, typename Operation>
    void forEachHessian( Points const& points, Operation&& operation ) const
    {
        // Keep the existing coefficient-operator implementation in one callable.
        auto withOperators = [&]
        {
            M_order.withValues( M_basis, points, [&]( auto const& modal )
            {
                for ( int i = 0; i < nDim; ++i )
                {
                    auto differentiated = M_coeff.multiply( eigenMap( basisD( i ) ) );
                    for ( int j = 0; j < nDim; ++j )
                        operation( i, j, differentiated * eigenMap( basisD( j ) ) * modal );
                }
            } );
        };

        // The common order policy selects the measured sampled-Hessian backend.
        if constexpr ( order_storage_type::template hasDirectHessians<Points> )
        {
            if ( M_order.useDirectHessians() )
            {
                M_order.forEachHessian( M_basis, points, [&]( int i, int j, auto const& modal )
                                        { operation( i, j, M_coeff.matrix() * modal ); } );
                return;
            }
        }

        withOperators();
    }

    /** @brief Evaluate modal derivatives at the effective order. */
    template<typename AE>
    ublas::vector<matrix_type> basisDerivate( ublas::matrix_expression<AE> const& pts ) const
    {
        return M_order.derivate( M_basis, pts );
    }

    /** @brief Return a stable differentiation operator. */
    matrix_type const& basisD( uint16_type i ) const
    {
        return M_order.d( M_basis, i );
    }

    /** Shared evaluation policy; matching static bases need no state. */
    using order_storage_type = ::Feel::detail::PolynomialBasisState<basis_type, OrderSpec>;

    /** @brief Effective expansion degree from the shared order policy. */
    [[nodiscard]] uint16_type orderValue() const
    {
        return M_order.order();
    }

    /** @brief Set runtime degree or validate a static degree. */
    void setOrderValue( uint16_type order )
    {
        M_order.setOrder( order );
    }


private:
    /** @brief Retain basis/order state while adopting independently owned coefficients. */
    PolynomialSet( self_type const& source, CoefficientStorage&& coefficients )
        : M_order( source.M_order ), M_basis( source.M_basis ),
          M_coeff( std::move( coefficients ) ), M_fname( source.M_fname ) {}

    /** Order policy and owned differentiation operators. */
    [[no_unique_address]] order_storage_type M_order{};
    basis_type M_basis;
    /** Sole owning coefficient representation; legacy access returns an owning copy. */
    CoefficientStorage M_coeff;
    std::string M_fname;
};

namespace detail
{
/**
 * @brief Materialize a polynomial subset with static or runtime-order storage.
 *
 * Family builders use this as their only storage-policy seam.  Static sets
 * retain their compile-time extent while runtime sets receive the effective
 * semantic order explicitly; the mathematical coefficient path is otherwise
 * identical.
 */
template<typename PolynomialSetType, typename SourcePolynomialSet>
[[nodiscard]] PolynomialSetType
polynomialSetFrom( SourcePolynomialSet const& source, uint16_type order )
{
    PolynomialSetType result = [&]()
    {
        if constexpr ( std::constructible_from<PolynomialSetType,
                                                SourcePolynomialSet const&> )
            return PolynomialSetType( source );
        else
            return PolynomialSetType( source.coeff(), true );
    }();
    if constexpr ( requires( PolynomialSetType& p ) { p.setOrder( order ); } )
        result.setOrder( order );
    return result;
}
} // namespace detail

template<typename Poly, template<uint16_type> class PolySetType = Scalar>
using DynamicPolynomialSet = PolynomialSet<Poly, PolySetType, Dynamic>;

template<typename Poly,template<uint16_type> class PolySetType, int OrderSpec>
const uint16_type PolynomialSet<Poly,PolySetType,OrderSpec>::nComponents;
template<typename Poly,template<uint16_type> class PolySetType, int OrderSpec>
const uint16_type PolynomialSet<Poly,PolySetType,OrderSpec>::nComponents1;
template<typename Poly,template<uint16_type> class PolySetType, int OrderSpec>
const uint16_type PolynomialSet<Poly,PolySetType,OrderSpec>::nComponents2;

} // Feel

#include <feel/feelpoly/polynomialset_impl.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>

#endif /* __PolynomialSet_H */
