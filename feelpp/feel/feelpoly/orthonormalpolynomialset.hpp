/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
#ifndef FEELPP_ORTHONORMALPOLYNOMIALSET_HPP
#define FEELPP_ORTHONORMALPOLYNOMIALSET_HPP 1

#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/concepts.hpp>

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
public PolynomialSet<Dubiner<Dim, RealDim, (Order >= 0 ? Order : 1), Normalized<true>, T, StorageUBlas>, PolySetType >
{
    // For Dynamic order, use Order=1 as compile-time placeholder; actual order is runtime
    static constexpr int CompileTimeOrder = (Order >= 0 ? Order : 1);
    typedef PolynomialSet<Dubiner<Dim, RealDim, CompileTimeOrder, Normalized<true>, T, StorageUBlas>, PolySetType > super;
public:

    static const uint16_type TAG = TheTAG;
    static const uint16_type nDim = Dim;
    static const int nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static inline const bool isTransformationEquivalent = true;

    //! True if order is determined at runtime (Order == Dynamic)
    static constexpr bool is_order_dynamic = (Order == Dynamic);
    //! True if order is determined at compile-time
    static constexpr bool is_order_static = !is_order_dynamic;

    typedef OrthonormalPolynomialSet<Dim, Order,RealDim, PolySetType, T, TheTAG, Simplex> self_type;
    typedef self_type component_basis_type;

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
    template<int O>
    struct convex
    {
        typedef Simplex<Dim, O, /*RealDim*/Dim> type;
    };
    typedef Reference<convex_type, nDim, CompileTimeOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex (compile-time, use runtimeDofPerVertex() for dynamic)
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge (compile-time, use runtimeDofPerEdge() for dynamic)
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face (compile-time, use runtimeDofPerFace() for dynamic)
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume (compile-time, use runtimeDofPerVolume() for dynamic)
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    //!< Compile-time local DOF count (for static order) - use runtimeLocalDof() for dynamic
    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    /**
     * @brief Get polynomial order (semantic runtime accessor)
     *
     * Returns the runtime semantic order, including when this type is used as
     * a low-order compile-time placeholder in dynamic FE producer code paths.
     */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return static_cast<uint16_type>( super::order() );
    }

    /**
     * @brief Get runtime polynomial order
     * @deprecated Use order() instead - unified interface handles both static and dynamic cases
     */
    [[nodiscard]] uint16_type runtimeOrder() const noexcept
    {
        return order();
    }

    /**
     * @brief Get local DOF count (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return static_cast<uint16_type>( ::Feel::detail::simplexTotal( nDim, this->order() ) );
    }

    /**
     * @brief Get runtime local DOF count
     * @deprecated Use localDof() instead - unified interface handles both static and dynamic cases
     */
    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
    {
        return localDof();
    }

    /**
     * @brief Get DOF per vertex (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerVertex() const noexcept
    {
        return ::Feel::detail::simplexPerVertex( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerVertex() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerVertex() const noexcept
    {
        return dofPerVertex();
    }

    /**
     * @brief Get DOF per edge (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return ::Feel::detail::simplexPerEdge( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerEdge() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerEdge() const noexcept
    {
        return dofPerEdge();
    }

    /**
     * @brief Get DOF per face (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        return ::Feel::detail::simplexPerFace( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerFace() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerFace() const noexcept
    {
        return dofPerFace();
    }

    /**
     * @brief Get DOF per volume (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return ::Feel::detail::simplexPerVolume( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerVolume() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
    {
        return dofPerVolume();
    }
    typedef typename matrix_node<value_type>::type points_type;

    /**
     * local interpolant is undefined
     */
    using  local_interpolant_type = std::monostate;
    using  local_interpolants_type = std::monostate;

    struct SSpace
    {
        static constexpr uint16_type TheOrder = (Order > 1)?Order-1:0;
        typedef typename mpl::if_<mpl::less_equal<mpl::int_<Order>, mpl::int_<1> >,
                                  mpl::identity<OrthonormalPolynomialSet<Dim, 0, RealDim, PolySetType, T, TheTAG,Simplex> >,
                                  mpl::identity<OrthonormalPolynomialSet<Dim, TheOrder, RealDim, PolySetType, T, TheTAG,Simplex> > >::type::type type;

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
            const uint16_type n = static_cast<uint16_type>( nComponents * convex_type::polyDims( Order ) );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }
        else
        {
            // For Dynamic with default constructor, use order 1 as default
            const uint16_type n = static_cast<uint16_type>( nComponents * this->localDof() );
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
        // Set runtime order in base PolynomialSet class for isUsingDynamicOrder() check
        this->set_order_value( ro.value );

        const uint16_type n = static_cast<uint16_type>( nComponents * this->localDof() );
        ublas::matrix<value_type> m( n, n );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
        mMap.setIdentity();
        this->setCoefficient( polyset_type::toType( m ), true );

        initSymmetricMapping();
    }

private:
    /**
     * @brief Initialize symmetric index mapping for tensor2symm case
     */
    void initSymmetricMapping()
    {
        if constexpr ( is_tensor2symm )
        {
            const uint16_type localDof = runtimeLocalDof();
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

    /**
     * @brief Evaluate the underlying basis at given points (Simplex specialization)
     *
     * For static order types, uses compile-time evaluation. For dynamic order
     * scenarios (when semantic runtime order differs from compile-time order), uses
     * runtime evaluation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Basis evaluation matrix (nBasis x nPoints)
     */
    template<typename AE>
    typename super::matrix_type basisEvaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        // Check if runtime order differs from compile-time order
        if ( runtime_order != CompileTimeOrder )
        {
            // Explicit low-order runtime dispatch to preserve fast kernels for P0/P1/P2.
            return basisEvaluateRuntimeLowOrderDispatch( __pts, runtime_order );
        }
        else
        {
            // Use static path for matching order
            return this->basis()( __pts );
        }
    }

    /**
     * @brief Evaluate the polynomial set at given points
     *
     * For dynamic order, uses the runtime order for basis evaluation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Evaluation matrix (nLocalDof x nPoints)
     */
    template<typename AE>
    typename super::matrix_type evaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        if ( runtime_order != CompileTimeOrder )
        {
            // Use runtime order when it differs from compile-time order.
            typename super::matrix_type m( basisEvaluateRuntimeLowOrderDispatch( __pts, runtime_order ) );
            return ublas::prod( this->coeff(), m );
        }
        else
        {
            // Compile-time/static path.
            return super::evaluate( __pts );
        }
    }

    /**
     * @brief Derivate the polynomial set at given points
     *
     * For dynamic order, uses the runtime order for basis derivation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Vector of derivation matrices (one per dimension)
     */
    template<typename AE>
    ublas::vector<typename super::matrix_type> derivate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        if ( runtime_order != CompileTimeOrder )
        {
            // Use runtime order when it differs from compile-time order.
            ublas::vector<typename super::matrix_type> der(
                basisDerivateRuntimeLowOrderDispatch( __pts, runtime_order ) );
            ublas::vector<typename super::matrix_type> res( nDim );

            for ( uint16_type i = 0; i < nDim; ++i )
            {
                res[i].resize( this->coeff().size1(), __pts().size2() );
                ublas::axpy_prod( this->coeff(), der[i], res[i] );
            }

            return res;
        }
        else
        {
            // Compile-time/static path.
            return super::derivate( __pts );
        }
    }

private:
    template<typename AE>
    static typename super::matrix_type
    basisEvaluateRuntimeLowOrderDispatch( ublas::matrix_expression<AE> const& __pts, uint16_type runtimeOrder )
    {
        using basis_o0_type = Dubiner<Dim, RealDim, 0, Normalized<true>, T, StorageUBlas>;
        using basis_o1_type = Dubiner<Dim, RealDim, 1, Normalized<true>, T, StorageUBlas>;
        using basis_o2_type = Dubiner<Dim, RealDim, 2, Normalized<true>, T, StorageUBlas>;
        switch ( runtimeOrder )
        {
        case 0:
            return basis_o0_type::evaluate( __pts );
        case 1:
            return basis_o1_type::evaluate( __pts );
        case 2:
            return basis_o2_type::evaluate( __pts );
        default:
            return basis_type::evaluate( __pts, runtimeOrder );
        }
    }

    template<typename AE>
    static ublas::vector<typename super::matrix_type>
    basisDerivateRuntimeLowOrderDispatch( ublas::matrix_expression<AE> const& __pts, uint16_type runtimeOrder )
    {
        using basis_o0_type = Dubiner<Dim, RealDim, 0, Normalized<true>, T, StorageUBlas>;
        using basis_o1_type = Dubiner<Dim, RealDim, 1, Normalized<true>, T, StorageUBlas>;
        using basis_o2_type = Dubiner<Dim, RealDim, 2, Normalized<true>, T, StorageUBlas>;
        switch ( runtimeOrder )
        {
        case 0:
            return basis_o0_type::derivate( __pts );
        case 1:
            return basis_o1_type::derivate( __pts );
        case 2:
            return basis_o2_type::derivate( __pts );
        default:
            return basis_type::derivate( __pts, runtimeOrder );
        }
    }

public:
    OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar,T, TheTAG, Simplex > toScalar() const
    {
        return OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar,T, TheTAG, Simplex >();
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const override
    {
        return "dubiner";
    }

    //! \return the component of a local dof
    uint16_type component( uint16_type localDofId ) const
        {
            uint16_type comp = localDofId/nLocalDof;
            DCHECK( comp < nComponents ) << "invalid localDofId " << localDofId;
            return comp;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const
        {
            uint16_type ldofParent = localDofId % nLocalDof;
            return ldofParent;
        }

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const
        {
            return 1;
        }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const
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
public PolynomialSet<Legendre<Dim, RealDim, (Order >= 0 ? Order : 1), Normalized<true>, T>, PolySetType >
{
    // For Dynamic order, use Order=1 as compile-time placeholder; actual order is runtime
    static constexpr int CompileTimeOrder = (Order >= 0 ? Order : 1);
    typedef PolynomialSet<Legendre<Dim, RealDim, CompileTimeOrder, Normalized<true>, T>, PolySetType > super;
public:

    static const uint16_type TAG = TheTAG;
    static const uint16_type nDim = Dim;
    static const int nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static inline const bool isTransformationEquivalent = true;

    //! True if order is determined at runtime (Order == Dynamic)
    static constexpr bool is_order_dynamic = (Order == Dynamic);
    //! True if order is determined at compile-time
    static constexpr bool is_order_static = !is_order_dynamic;

    typedef OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Hypercube> self_type;
    typedef self_type component_basis_type;

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
    template<int O>
    struct convex
    {
        typedef Hypercube<Dim, O, nDim/*RealDim*/> type;
    };
    typedef Reference<convex_type, nDim, CompileTimeOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex (compile-time, use runtimeDofPerVertex() for dynamic)
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge (compile-time, use runtimeDofPerEdge() for dynamic)
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face (compile-time, use runtimeDofPerFace() for dynamic)
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume (compile-time, use runtimeDofPerVolume() for dynamic)
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    //!< Compile-time local DOF count (for static order) - use runtimeLocalDof() for dynamic
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

    struct SSpace
    {
        static constexpr uint16_type TheOrder = (Order > 1)?Order-1:0;
        typedef typename mpl::if_<mpl::less_equal<mpl::int_<Order>, mpl::int_<1> >,
                                  mpl::identity<OrthonormalPolynomialSet<Dim, 0, RealDim, PolySetType, T, TheTAG, Hypercube> >,
                                  mpl::identity<OrthonormalPolynomialSet<Dim, TheOrder, RealDim, PolySetType, T, TheTAG, Hypercube> > >::type::type type;

    };
    template<int OtherOrder>
    struct ChangeOrder
    {
        typedef OrthonormalPolynomialSet<Dim, OtherOrder, RealDim, PolySetType, T, TheTAG, Hypercube> type;
    };

    /**
     * @brief Get polynomial order (semantic runtime accessor)
     *
     * Returns the runtime semantic order, including when this type is used as
     * a low-order compile-time placeholder in dynamic FE producer code paths.
     */
    [[nodiscard]] uint16_type order() const noexcept
    {
        return static_cast<uint16_type>( super::order() );
    }

    /**
     * @brief Get runtime polynomial order
     * @deprecated Use order() instead - unified interface handles both static and dynamic cases
     */
    [[nodiscard]] uint16_type runtimeOrder() const noexcept
    {
        return order();
    }

    /**
     * @brief Get local DOF count (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return hypercubePolyDims( nDim, this->order() );
    }

    /**
     * @brief Get runtime local DOF count
     * @deprecated Use localDof() instead - unified interface handles both static and dynamic cases
     */
    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
    {
        return localDof();
    }

    /**
     * @brief Get DOF per vertex (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerVertex() const noexcept
    {
        return ::Feel::detail::hypercubePerVertex( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerVertex() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerVertex() const noexcept
    {
        return dofPerVertex();
    }

    /**
     * @brief Get DOF per edge (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return ::Feel::detail::hypercubePerEdge( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerEdge() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerEdge() const noexcept
    {
        return dofPerEdge();
    }

    /**
     * @brief Get DOF per face (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        return ::Feel::detail::hypercubePerFace( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerFace() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerFace() const noexcept
    {
        return dofPerFace();
    }

    /**
     * @brief Get DOF per volume (semantic runtime accessor)
     */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return ::Feel::detail::hypercubePerVolume( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerVolume() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
    {
        return dofPerVolume();
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
            const uint16_type n = static_cast<uint16_type>( nComponents * convex_type::polyDims( Order ) );
            ublas::matrix<value_type> m( n, n );
            Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
            mMap.setIdentity();
            this->setCoefficient( polyset_type::toType( m ), true );
        }
        else
        {
            // For Dynamic with default constructor, use order 1 as default
            const uint16_type n = static_cast<uint16_type>( nComponents * this->localDof() );
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
        // Set runtime order in base PolynomialSet class for isUsingDynamicOrder() check
        this->set_order_value( ro.value );

        const uint16_type n = static_cast<uint16_type>( nComponents * this->localDof() );
        ublas::matrix<value_type> m( n, n );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mMap( m.data().begin(), m.size1(), m.size2() );
        mMap.setIdentity();
        this->setCoefficient( polyset_type::toType( m ), true );

        initSymmetricMapping();
    }

private:
    /**
     * @brief Compute hypercube polynomial dimensions at runtime
     * @param dim Spatial dimension
     * @param order Polynomial order
     * @return (order+1)^dim
     */
    static uint16_type hypercubePolyDims( uint16_type dim, uint16_type order )
    {
        uint16_type result = 1;
        for ( uint16_type d = 0; d < dim; ++d )
            result *= ( order + 1 );
        return result;
    }

    /**
     * @brief Initialize symmetric index mapping for tensor2symm case
     */
    void initSymmetricMapping()
    {
        if constexpr ( is_tensor2symm )
        {
            const uint16_type localDof = runtimeLocalDof();
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

    /**
     * @brief Evaluate the underlying basis at given points (Hypercube specialization)
     *
     * For static order types, uses compile-time evaluation. For dynamic order
     * scenarios (when semantic runtime order differs from compile-time order), uses
     * runtime evaluation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Basis evaluation matrix (nBasis x nPoints)
     */
    template<typename AE>
    typename super::matrix_type basisEvaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        // Check if runtime order differs from compile-time order
        if ( runtime_order != CompileTimeOrder )
        {
            // Explicit low-order runtime dispatch to preserve fast kernels for Q0/Q1/Q2.
            return basisEvaluateRuntimeLowOrderDispatch( __pts, runtime_order );
        }
        else
        {
            // Use static path for matching order
            return this->basis()( __pts );
        }
    }

    /**
     * @brief Evaluate the polynomial set at given points
     *
     * For dynamic order, uses the runtime order for basis evaluation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Evaluation matrix (nLocalDof x nPoints)
     */
    template<typename AE>
    typename super::matrix_type evaluate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        if ( runtime_order != CompileTimeOrder )
        {
            // Use runtime order when it differs from compile-time order.
            typename super::matrix_type m( basisEvaluateRuntimeLowOrderDispatch( __pts, runtime_order ) );
            return ublas::prod( this->coeff(), m );
        }
        else
        {
            // Compile-time/static path.
            return super::evaluate( __pts );
        }
    }

    /**
     * @brief Derivate the polynomial set at given points
     *
     * For dynamic order, uses the runtime order for basis derivation.
     *
     * @param __pts Points to evaluate at (nDim x nPoints matrix)
     * @return Vector of derivation matrices (one per dimension)
     */
    template<typename AE>
    ublas::vector<typename super::matrix_type> derivate( ublas::matrix_expression<AE> const& __pts ) const
    {
        const auto runtime_order = this->order();
        if ( runtime_order != CompileTimeOrder )
        {
            // Use runtime order when it differs from compile-time order.
            ublas::vector<typename super::matrix_type> der(
                basisDerivateRuntimeLowOrderDispatch( __pts, runtime_order ) );
            ublas::vector<typename super::matrix_type> res( nDim );

            for ( uint16_type i = 0; i < nDim; ++i )
            {
                res[i].resize( this->coeff().size1(), __pts().size2() );
                ublas::axpy_prod( this->coeff(), der[i], res[i] );
            }

            return res;
        }
        else
        {
            // Compile-time/static path.
            return super::derivate( __pts );
        }
    }

private:
    template<typename AE>
    static typename super::matrix_type
    basisEvaluateRuntimeLowOrderDispatch( ublas::matrix_expression<AE> const& __pts, uint16_type runtimeOrder )
    {
        using basis_o0_type = Legendre<Dim, RealDim, 0, Normalized<true>, T>;
        using basis_o1_type = Legendre<Dim, RealDim, 1, Normalized<true>, T>;
        using basis_o2_type = Legendre<Dim, RealDim, 2, Normalized<true>, T>;
        switch ( runtimeOrder )
        {
        case 0:
            return basis_o0_type::evaluate( __pts );
        case 1:
            return basis_o1_type::evaluate( __pts );
        case 2:
            return basis_o2_type::evaluate( __pts );
        default:
            return basis_type::evaluate( __pts, runtimeOrder );
        }
    }

    template<typename AE>
    static ublas::vector<typename super::matrix_type>
    basisDerivateRuntimeLowOrderDispatch( ublas::matrix_expression<AE> const& __pts, uint16_type runtimeOrder )
    {
        using basis_o0_type = Legendre<Dim, RealDim, 0, Normalized<true>, T>;
        using basis_o1_type = Legendre<Dim, RealDim, 1, Normalized<true>, T>;
        using basis_o2_type = Legendre<Dim, RealDim, 2, Normalized<true>, T>;
        switch ( runtimeOrder )
        {
        case 0:
            return basis_o0_type::derivate( __pts );
        case 1:
            return basis_o1_type::derivate( __pts );
        case 2:
            return basis_o2_type::derivate( __pts );
        default:
            return basis_type::derivate( __pts, runtimeOrder );
        }
    }

public:
    OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Hypercube > toScalar() const
    {
        return OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar, T, TheTAG, Hypercube >();
    }

    /**
     * \return the family name of the polynomial set
     */
    std::string familyName() const override
    {
        return "legendre";
    }

    //! \return the component of a local dof
    uint16_type component( uint16_type localDofId ) const
        {
            uint16_type comp = localDofId/nLocalDof;
            DCHECK( comp < nComponents ) << "invalid localDofId " << localDofId;
            return comp;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const
        {
            uint16_type ldofParent = localDofId % nLocalDof;
            return ldofParent;
        }

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const
        {
            return 1;
        }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const
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
/// \encond

template<int Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         uint16_type TheTAG=0 >
class OrthonormalPolynomialSet
{
public:
    template<uint16_type N,
             uint16_type RealDim,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                                  mpl::identity<Feel::detail::OrthonormalPolynomialSet<N,Order,RealDim,PolySetType,T,TheTAG,Simplex> >,
                                  mpl::identity<Feel::detail::OrthonormalPolynomialSet<N,Order,RealDim,PolySetType,T,TheTAG,Hypercube> > >::type::type result_type;
    typedef result_type type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef OrthonormalPolynomialSet<Order,PolySetType,TheNewTAG> type;
    };

    typedef OrthonormalPolynomialSet<Order,Scalar,TheTAG> component_basis_type;

    static const uint16_type nOrder =  Order;
    static const uint16_type TAG = TheTAG;
    using  local_interpolant_type = std::monostate;
    using  local_interpolants_type = std::monostate;
};

} // Feel
#endif /* __OrthonormalPolynomialSet_H */
