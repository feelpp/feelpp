/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2009 Universite Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2009-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_ORTHONORMALPOLYNOMIALSET_HPP
#define FEELPP_ORTHONORMALPOLYNOMIALSET_HPP 1

#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/fefamilyshell.hpp>

#include <boost/mp11/utility.hpp>

#include <feel/feelpoly/referencedoflayout.hpp>

namespace Feel
{
/// \cond DETAIL
namespace detail
{
/**
 * \internal
 * \class OrthonormalPolynomialSet
 * \brief a set of orthonormal polynomials over a convex
 *
 * On the simplicies we use the Dubiner basis
 *
 */
template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType = Scalar,
         typename T = double,
         uint16_type TheTAG = 0,
         template<int,int,int> class Convex = Simplex>
class OrthonormalPolynomialSet
{};

template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG>
class OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Simplex>
    :
public PolynomialSet<Dubiner<Dim, RealDim, (Order >= 0 ? Order : 1), Normalized<true>, T, StorageUBlas>, PolySetType, Order >,
public ::Feel::detail::ReferenceDofLayout<OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Simplex>, true>
{
    // For Dynamic order, use Order=1 as compile-time placeholder; actual order is runtime
    static constexpr int CompileTimeOrder = (Order >= 0 ? Order : 1);
    typedef PolynomialSet<Dubiner<Dim, RealDim, CompileTimeOrder, Normalized<true>, T, StorageUBlas>, PolySetType, Order > super;
public:

    static const uint16_type TAG = TheTAG;
    static const uint16_type nDim = Dim;
    /** Authoritative signed order specification. */
    static constexpr int order_spec = Order;
    static const int nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static inline const bool isTransformationEquivalent = true;

    //! True if order is determined at runtime (Order == Dynamic)
    static constexpr bool is_order_dynamic = (Order == Dynamic);
    //! True if order is determined at compile-time
    static constexpr bool is_order_static = !is_order_dynamic;

    typedef OrthonormalPolynomialSet<Dim, Order,RealDim, PolySetType, T, TheTAG, Simplex> self_type;
    /** Scalar component family retaining the same order policy. */
    using component_basis_type = OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Simplex>;

    typedef typename super::polyset_type polyset_type;
    static inline const bool is_tensor2 = polyset_type::is_tensor2;
    static inline const bool is_tensor2symm = polyset_type::is_tensor2 && is_symm_v<polyset_type>;
    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const bool is_continuous = false;
    static inline const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static inline const bool is_product = true;
    static inline const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Dubiner<Dim, RealDim, CompileTimeOrder, Normalized<true>, T, StorageUBlas> basis_type;
    typedef Simplex<Dim, CompileTimeOrder, /*RealDim*/Dim> convex_type;
    /** @brief Rebind the reference simplex degree. */
    template<int O>
    struct Convex
    {
        typedef Simplex<Dim, O, /*RealDim*/Dim> type;
    };
    /** Compatibility alias for the reference-cell rebind protocol. */
    template<int O> using convex = Convex<O>;
    typedef Reference<convex_type, nDim, CompileTimeOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex (compile-time, use dofPerVertex() for dynamic)
    static const uint16_type nDofPerVertex = nDim == 0 ? convex_type::numPoints : 0;
    //!< Number of degrees  of freedom per edge (compile-time, use dofPerEdge() for dynamic)
    static const uint16_type nDofPerEdge = nDim == 1 ? convex_type::numPoints : 0;
    //!< Number of degrees  of freedom per face (compile-time, use dofPerFace() for dynamic)
    static const uint16_type nDofPerFace = nDim == 2 ? convex_type::numPoints : 0;

    //!< Number of degrees  of freedom per volume (compile-time, use dofPerVolume() for dynamic)
    static const uint16_type nDofPerVolume = nDim == 3 ? convex_type::numPoints : 0;

    //!< Compile-time local DOF count (for static order) - use localDof() for dynamic
    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    /**
     * @brief Get polynomial order (semantic runtime accessor)
     *
     * The signed order policy selects static or runtime behavior independently
     * of the degree used to instantiate the legacy modal kernel.
     */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return static_cast<uint16_type>( super::order() );
    }

    /**
     * @brief Get local DOF count (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type localDof() const
    {
        return ::Feel::detail::checkedLocalIndexValue(
            ::Feel::detail::simplexTotal( nDim, this->order() ) );
    }


    typedef typename matrix_node<value_type>::type points_type;

    /**
     * local interpolant is undefined
     */
    using  local_interpolant_type = std::monostate;
    using  local_interpolants_type = std::monostate;

    /** @brief Legacy lower-degree approximation space with a constant floor. */
    struct SSpace
    {
        /** Degree of the lower approximation, never below zero. */
        static constexpr int TheOrder = Order == Dynamic ? Dynamic : (Order > 1 ? Order-1 : 0);
        /** Simplex modal implementation at the lower approximation degree. */
        using type = OrthonormalPolynomialSet<Dim, TheOrder, RealDim, PolySetType, T, TheTAG, Simplex>;

    };
    template<int OtherOrder>
    struct ChangeOrder
    {
        typedef OrthonormalPolynomialSet<Dim, OtherOrder, RealDim, PolySetType, T, TheTAG,Simplex> type;
    };

    /**
     * @brief Default constructor for static order
     *
     * For Dynamic order, use the RuntimeOrder constructor instead.
     */
    OrthonormalPolynomialSet()
        :
        super( basis_type() )
    {
        if constexpr ( is_order_static )
        {
            const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
                static_cast<size_type>( nComponents ) * convex_type::polyDims( Order ) );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }
        else
        {
            // For Dynamic with default constructor, use order 1 as default
            const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
                static_cast<size_type>( nComponents ) * this->localDof() );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }

        initSymmetricMapping();
    }

    /**
     * @brief Constructor with runtime order specification
     *
     * Use this constructor when Order == Dynamic to specify the polynomial order at runtime.
     *
     * @param ro The runtime order specification
     */
    explicit OrthonormalPolynomialSet( RuntimeOrder ro )
        :
        super( basis_type() )
    {
        ro = validateOrderAndCardinality( ro );
        // Validate static degree or store the effective runtime degree.
        this->setOrderValue( ro.value );

        const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
            static_cast<size_type>( nComponents ) * this->localDof() );
        ublas::matrix<value_type> m( n, n );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
        mMap.setIdentity();
        this->setCoefficient( polyset_type::toType( m ), true );

        initSymmetricMapping();
    }

private:
    [[nodiscard]] static RuntimeOrder
    validateOrderAndCardinality( RuntimeOrder order )
    {
        order = ::Feel::detail::validatedNodalOrder<
            ::Feel::detail::SimplexCardinality>( nDim, order, nComponents );
        return order;
    }

    /**
     * @brief Initialize symmetric index mapping for tensor2symm case
     */
    void initSymmetricMapping()
    {
        if constexpr ( is_tensor2symm )
        {
            const uint16_type localDof = this->localDof();
            M_unsymm2symm.resize( nComponents * localDof );
            for ( uint16_type l = 0; l < localDof; ++l )
            {
                for ( int c1 = 0; c1 < nComponents1; ++c1 )
                {
                    for ( int c2 = c1 + 1; c2 < nComponents2; ++c2 )
                    {
                        const int k = Feel::detail::symmetricIndex( c1, c2, nComponents1 );
                        M_unsymm2symm[localDof * ( nComponents1 * c1 + c2 ) + l] = localDof * k + l;
                        M_unsymm2symm[localDof * ( nComponents1 * c2 + c1 ) + l] = localDof * k + l;
                    }
                    const int k = Feel::detail::symmetricIndex( c1, c1, nComponents1 );
                    M_unsymm2symm[localDof * ( nComponents1 * c1 + c1 ) + l] = localDof * k + l;
                }
            }
        }
    }

public:
    OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar,T, TheTAG, Simplex > toScalar() const
    {
        return makeOrderAware<OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar,T, TheTAG, Simplex >>(RuntimeOrder{this->order()});
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const override
    {
        return "dubiner";
    }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const override
        {
            if ( !is_tensor2symm )
                return i;
            DCHECK( M_unsymm2symm.size() > i ) << "invalid size of unsymm2symm container";
            return M_unsymm2symm[i];
        }

    points_type points() const
    {
        return points_type();
    }
    points_type points( int f ) const
    {
        return points_type();
    }

private :
    std::vector<uint16_type> M_unsymm2symm;
};

template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG>
const uint16_type OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType,T, TheTAG, Simplex>::nLocalDof;


template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG>
class OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Hypercube>
    :
public PolynomialSet<Legendre<Dim, RealDim, (Order >= 0 ? Order : 1), Normalized<true>, T>, PolySetType, Order >,
public ::Feel::detail::ReferenceDofLayout<OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Hypercube>, true>
{
    // For Dynamic order, use Order=1 as compile-time placeholder; actual order is runtime
    static constexpr int CompileTimeOrder = (Order >= 0 ? Order : 1);
    typedef PolynomialSet<Legendre<Dim, RealDim, CompileTimeOrder, Normalized<true>, T>, PolySetType, Order > super;
public:

    static const uint16_type TAG = TheTAG;
    static const uint16_type nDim = Dim;
    /** Authoritative signed order specification. */
    static constexpr int order_spec = Order;
    static const int nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static inline const bool isTransformationEquivalent = true;

    //! True if order is determined at runtime (Order == Dynamic)
    static constexpr bool is_order_dynamic = (Order == Dynamic);
    //! True if order is determined at compile-time
    static constexpr bool is_order_static = !is_order_dynamic;

    typedef OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Hypercube> self_type;
    /** Scalar component family retaining the same order policy. */
    using component_basis_type = OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Hypercube>;

    typedef typename super::polyset_type polyset_type;
    static inline const bool is_tensor2 = polyset_type::is_tensor2;
    static inline const bool is_tensor2symm = polyset_type::is_tensor2 && is_symm_v<polyset_type>;
    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const bool is_continuous = false;
    static inline const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static inline const bool is_product = true;
    static inline const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    typedef typename super::component_type component_type;
    typedef T value_type;
    typedef Legendre<Dim, RealDim, CompileTimeOrder, Normalized<true>, T> basis_type;
    typedef Hypercube<Dim, CompileTimeOrder, /*RealDim*/Dim> convex_type;
    /** @brief Rebind the reference hypercube degree. */
    template<int O>
    struct Convex
    {
        typedef Hypercube<Dim, O, nDim/*RealDim*/> type;
    };
    /** Compatibility alias for the reference-cell rebind protocol. */
    template<int O> using convex = Convex<O>;
    typedef Reference<convex_type, nDim, CompileTimeOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex (compile-time, use dofPerVertex() for dynamic)
    static const uint16_type nDofPerVertex = nDim == 0 ? convex_type::numPoints : 0;
    //!< Number of degrees  of freedom per edge (compile-time, use dofPerEdge() for dynamic)
    static const uint16_type nDofPerEdge = nDim == 1 ? convex_type::numPoints : 0;
    //!< Number of degrees  of freedom per face (compile-time, use dofPerFace() for dynamic)
    static const uint16_type nDofPerFace = nDim == 2 ? convex_type::numPoints : 0;

    //!< Number of degrees  of freedom per volume (compile-time, use dofPerVolume() for dynamic)
    static const uint16_type nDofPerVolume = nDim == 3 ? convex_type::numPoints : 0;

    //!< Compile-time local DOF count (for static order) - use localDof() for dynamic
    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    typedef typename matrix_node<value_type>::type points_type;

    /**
     * local interpolant is undefined
     */
    using  local_interpolant_type = std::monostate;
    using  local_interpolants_type = std::monostate;

    /** @brief Legacy lower-degree approximation space with a constant floor. */
    struct SSpace
    {
        /** Degree of the lower approximation, never below zero. */
        static constexpr int TheOrder = Order == Dynamic ? Dynamic : (Order > 1 ? Order-1 : 0);
        /** Hypercube modal implementation at the lower approximation degree. */
        using type = OrthonormalPolynomialSet<Dim, TheOrder, RealDim, PolySetType, T, TheTAG, Hypercube>;

    };
    template<int OtherOrder>
    struct ChangeOrder
    {
        typedef OrthonormalPolynomialSet<Dim, OtherOrder, RealDim, PolySetType, T, TheTAG, Hypercube> type;
    };

    /**
     * @brief Get polynomial order (semantic runtime accessor)
     *
     * The signed order policy selects static or runtime behavior independently
     * of the degree used to instantiate the legacy modal kernel.
     */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return static_cast<uint16_type>( super::order() );
    }

    /**
     * @brief Get local DOF count (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type localDof() const
    {
        return ::Feel::detail::checkedLocalIndexValue(
            ::Feel::detail::hypercubeTotal( nDim, this->order() ) );
    }


    /**
     * @brief Default constructor for static order
     *
     * For Dynamic order, use the RuntimeOrder constructor instead.
     */
    OrthonormalPolynomialSet()
        :
        super( basis_type() )
    {
        if constexpr ( is_order_static )
        {
            const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
                static_cast<size_type>( nComponents ) * convex_type::polyDims( Order ) );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }
        else
        {
            // For Dynamic with default constructor, use order 1 as default
            const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
                static_cast<size_type>( nComponents ) * this->localDof() );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }

        initSymmetricMapping();
    }

    /**
     * @brief Constructor with runtime order specification
     *
     * Use this constructor when Order == Dynamic to specify the polynomial order at runtime.
     *
     * @param ro The runtime order specification
     */
    explicit OrthonormalPolynomialSet( RuntimeOrder ro )
        :
        super( basis_type() )
    {
        ro = validateOrderAndCardinality( ro );
        // Validate static degree or store the effective runtime degree.
        this->setOrderValue( ro.value );

        const uint16_type n = ::Feel::detail::checkedLocalIndexValue(
            static_cast<size_type>( nComponents ) * this->localDof() );
        ublas::matrix<value_type> m( n, n );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
        mMap.setIdentity();
        this->setCoefficient( polyset_type::toType( m ), true );

        initSymmetricMapping();
    }

private:
    [[nodiscard]] static RuntimeOrder
    validateOrderAndCardinality( RuntimeOrder order )
    {
        order = ::Feel::detail::validatedNodalOrder<
            ::Feel::detail::HypercubeCardinality>( nDim, order, nComponents );
        return order;
    }

    /**
     * @brief Initialize symmetric index mapping for tensor2symm case
     */
    void initSymmetricMapping()
    {
        if constexpr ( is_tensor2symm )
        {
            const uint16_type localDof = this->localDof();
            M_unsymm2symm.resize( nComponents * localDof );
            for ( uint16_type l = 0; l < localDof; ++l )
            {
                for ( int c1 = 0; c1 < nComponents1; ++c1 )
                {
                    for ( int c2 = c1 + 1; c2 < nComponents2; ++c2 )
                    {
                        const int k = Feel::detail::symmetricIndex( c1, c2, nComponents1 );
                        M_unsymm2symm[localDof * ( nComponents1 * c1 + c2 ) + l] = localDof * k + l;
                        M_unsymm2symm[localDof * ( nComponents1 * c2 + c1 ) + l] = localDof * k + l;
                    }
                    const int k = Feel::detail::symmetricIndex( c1, c1, nComponents1 );
                    M_unsymm2symm[localDof * ( nComponents1 * c1 + c1 ) + l] = localDof * k + l;
                }
            }
        }
    }

public:
    OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Hypercube > toScalar() const
    {
        return makeOrderAware<OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Hypercube >>(RuntimeOrder{this->order()});
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const override
    {
        return "legendre";
    }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const override
        {
            if ( !is_tensor2symm )
                return i;
            DCHECK( M_unsymm2symm.size() > i ) << "invalid size of unsymm2symm container";
            return M_unsymm2symm[i];
        }

    points_type points() const
    {
        return points_type();
    }
    points_type points( int f ) const
    {
        return points_type();
    }

private:
    std::vector<uint16_type> M_unsymm2symm;
};

template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG>
const uint16_type OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType,T, TheTAG, Hypercube>::nLocalDof;
} // detail
/// \endcond

/** @brief Public descriptor selecting an orthonormal modal polynomial space.
 * @tparam Order Static expansion degree or `Dynamic`.
 * @tparam PolySetType Scalar, vector or tensor value-shape policy.
 * @tparam TheTAG Compatibility tag distinguishing basis instances.
 * @note The inherited signed `order_spec` is authoritative for dynamic order;
 * the inherited `nOrder` retains the common descriptor's legacy static value.
 */
template<int Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         uint16_type TheTAG=0 >
class OrthonormalPolynomialSet : public FiniteElementFamilyDescriptor<Order, TheTAG>
{
public:
    /** @brief Select the concrete modal space for a reference cell.
     * @tparam N Reference dimension.
     * @tparam RealDim Physical-space dimension.
     * @tparam T Scalar type.
     * @tparam Convex Simplex or hypercube cell type.
     */
    template<uint16_type N,
             uint16_type RealDim,
             typename T = double,
             typename Convex = Simplex<N> >
    struct Apply
    {
        /** Concrete implementation preserving the descriptor's order policy. */
        using result_type = boost::mp11::mp_if_c<Convex::is_simplex,
            Feel::detail::OrthonormalPolynomialSet<N, Order, RealDim, PolySetType, T, TheTAG, Simplex>,
            Feel::detail::OrthonormalPolynomialSet<N, Order, RealDim, PolySetType, T, TheTAG, Hypercube>>;
        /** Compatibility alias for the selected implementation. */
        typedef result_type type;
    };

    /** Concrete-element application retained for the established descriptor protocol. */
    template<uint16_type N, uint16_type RealDim, typename T = double, typename Convex = Simplex<N>>
    using apply = Apply<N, RealDim, T, Convex>;

    /** @brief Rebind the tag without changing expansion order or value shape.
     * @tparam TheNewTAG Replacement compatibility tag.
     */
    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        /** Descriptor retaining the existing order policy and the new tag. */
        typedef OrthonormalPolynomialSet<Order,PolySetType,TheNewTAG> type;
    };

    /** Scalar component descriptor preserving order specification and tag. */
    typedef OrthonormalPolynomialSet<Order,Scalar,TheTAG> component_basis_type;

    /** Modal bases do not define a nodal local-interpolant buffer. */
    using  local_interpolant_type = std::monostate;
    /** Modal bases do not define a collection of nodal interpolant buffers. */
    using  local_interpolants_type = std::monostate;
};

} // Feel
#endif /* __OrthonormalPolynomialSet_H */
