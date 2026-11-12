/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2008-2012 Université de Grenoble 1 (Joseph Fourier)
 * SPDX-FileCopyrightText: 2005-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file lagrange.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-18
 */
#ifndef __lagrange_H
#define __lagrange_H 1

#include <boost/ptr_container/ptr_vector.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/fekete.hpp>

#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/pointfunctionals.hpp>
#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/fefamilyshell.hpp>
#include <feel/feelpoly/isp0continuous.hpp>
#include <feel/feelpoly/order.hpp>

#include <concepts>




namespace Feel
{

namespace fem
{

/// \cond detail
namespace details
{
/** @brief Point-set families that accept a runtime degree and expose entity coordinates. */
template<template<class, int, class> class PointSetT, typename ConvexType, typename ValueType>
concept RuntimePointSet = requires( RuntimeOrder ro,
                                    PointSetT<ConvexType, Dynamic, ValueType> dyn_pset )
{
    { PointSetT<ConvexType, Dynamic, ValueType>( ro ) };
    { dyn_pset.order() } -> std::convertible_to<uint16_type>;
    { dyn_pset.nPoints() } -> std::convertible_to<uint32_type>;
    dyn_pset.points();
    dyn_pset.pointsBySubEntity( uint16_type{}, uint16_type{}, uint16_type{} );
};

/**
 * @brief Nodal dual owning interpolation points with the primal order policy.
 * @tparam Basis Primal polynomial space supplying the effective expansion degree.
 * @tparam PointSetType Interpolation-node family used to construct the functionals.
 */
template<typename Basis, template<class, int, class> class PointSetType>
class LagrangeDual
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    inline static const uint16_type nDim = super::nDim;
    inline static const uint16_type nOrder= super::nOrder;

    typedef typename super::primal_space_type primal_space_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    typedef typename primal_space_type::convex_type convex_type;
    typedef typename primal_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;

    /** Requested interpolation-node policy with the primal order specification. */
    using requested_pointset_type=PointSetType<convex_type,super::order_spec,value_type>;
    /** Owned point set, retaining the established equispaced fallback when necessary. */
    using pointset_type=mp::mp_if_c<super::is_order_dynamic &&
        !RuntimePointSet<PointSetType,convex_type,value_type>,
        PointSetEquiSpaced<convex_type,Dynamic,value_type>,requested_pointset_type>;

    template< template<class, int, class> class TestPointSetType >
    static constexpr bool is_pointset_v = std::is_base_of_v<TestPointSetType<convex_type, super::order_spec, value_type>,pointset_type >;
    inline static const uint16_type numPoints = reference_convex_type::numPoints;
    inline static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    inline static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    inline static const uint16_type nbPtsPerFace = reference_convex_type::nbPtsPerFace;
    inline static const uint16_type nbPtsPerVolume = reference_convex_type::nbPtsPerVolume;

    static const uint16_type nVertices = reference_convex_type::numVertices;
    static const uint16_type nFaces = reference_convex_type::numFaces;
    static const uint16_type nGeometricFaces = reference_convex_type::numFaces;
    static const uint16_type nEdges = reference_convex_type::numEdges;
    static const uint16_type nNormals = reference_convex_type::numNormals;


    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = nbPtsPerVertex;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume = nbPtsPerVolume;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static const uint16_type nLocalDof = numPoints;

    static constexpr uint16_type nFacesInConvex =
        ( nDim == 1 ) ? nVertices : ( nDim == 2 ) ? nEdges : nFaces;

    /** @brief Copy the owned primal, interpolation nodes and functional representations. */
    LagrangeDual( LagrangeDual const& d )
        :
        super( d ),
        M_convex_ref(),
        M_eid( d.M_eid ),
        M_pts( d.M_pts ),
        M_points_face( d.M_points_face ),
        M_points_by_subentity( d.M_points_by_subentity ),
        M_fset( d.M_fset ),
        M_pset( d.M_pset )
        {}

    /** @brief Move the complete dual while retaining point-set/order consistency. */
    LagrangeDual( LagrangeDual && d ) = default;

    /** @brief Construct the owned point set at the validated effective degree. */
    LagrangeDual( primal_space_type const& primal, RuntimeOrder order )
        : LagrangeDual(primal,makePointSet(primal,order)) {}

    /** @brief Construct the dual using the primal's effective expansion degree. */
    LagrangeDual( primal_space_type const& primal )
        : LagrangeDual(primal,RuntimeOrder{primal.order()}) {}

    /** @brief Own interpolation nodes and build the Ciarlet functionals once. */
    LagrangeDual( primal_space_type const& primal, pointset_type pts )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( pts.points() ),
        M_points_face( nFacesInConvex ),
        M_fset( primal ),
        M_pset( std::move(pts) )
    {
        if(M_pts.size2()!=primal.localDof())
            throw std::invalid_argument("interpolation point count disagrees with the primal degree");
        cacheSubEntityPoints( M_pset );

        if ( primal.order() > 0 )
        {
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                    e < M_convex_ref.entityRange( nDim-1 ).end();
                    ++e )
            {
                M_points_face[e] = M_points_by_subentity[nDim-1][e];
                DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
            }
        }

        M_fset.setFunctionalSet( functional::makePointEvaluationFunctionals( primal, M_pts ) );
    }

    ~LagrangeDual() = default;
    LagrangeDual& operator=( LagrangeDual const& ) = default;

    /** @brief Interpolation point set owned by this dual, including its runtime order. */
    pointset_type const& pointSet() const noexcept { return M_pset; }

    points_type const& points() const
    {
        return M_pts;
    }

    points_type const& points( uint16_type f ) const
    {
        return M_points_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_points_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_points_face[f], __i );
    }

#if 0
    std::vector<point_type> points( int topodim ) const
        {
            std::vector<point_type> pts( ;
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                  e < M_convex_ref.entityRange( nDim-1 ).end();
                  ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
            }
            return M_pset.pointsBySubEntity( topodim, edge, 1 );
        }
#endif
     points_type points( int topodim, int entity ) const
            {
                if ( topodim == nDim )
                    return M_pts;
                return M_points_by_subentity.at( topodim ).at( entity );
            }


    points_type edgePoints(int edge) const
        {
            if constexpr ( nDim == 1 )
                return M_pts;
            return M_points_by_subentity.at( 1 ).at( edge );
        }


    points_type vertexPoints(int vertex) const
        {
            return M_points_by_subentity.at( 0 ).at( vertex );
        }

    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }
private:

    /**
     * set the pointset at face \c f using points \c n
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_points_face[f].resize( n.size1(), n.size2(), false );
        M_points_face[f] = n;
    }

private:
    template<typename PointSetT>
    uint16_type pointSetPtsPerEntity( PointSetT const& pset,
                                      uint16_type entityDim ) const
    {
        using cardinality_type = std::conditional_t<
            convex_type::is_simplex,
            ::Feel::detail::SimplexCardinality,
            ::Feel::detail::HypercubeCardinality>;
        const uint16_type effectiveOrder = [&]()
        {
            if constexpr ( requires { pset.order(); } )
                return static_cast<uint16_type>( pset.order() );
            else
                return nOrder;
        }();
        return cardinality_type::onEntity(
            nDim, entityDim, effectiveOrder );
    }

    template<typename PointSetT>
    uint16_type subEntityPointCount( PointSetT const& pset, uint16_type topodim ) const
    {
        const auto ptsPerVertex = pointSetPtsPerEntity( pset, 0 );
        const auto ptsPerEdge = pointSetPtsPerEntity( pset, 1 );
        const auto ptsPerFace = pointSetPtsPerEntity( pset, 2 );
        const auto ptsPerVolume = pointSetPtsPerEntity( pset, 3 );

        switch ( topodim )
        {
        case 0:
            return ptsPerVertex;
        case 1:
            return static_cast<uint16_type>( 2 * ptsPerVertex + ptsPerEdge );
        case 2:
            if constexpr ( convex_type::is_simplex )
                return static_cast<uint16_type>( 3 * ptsPerVertex + 3 * ptsPerEdge + ptsPerFace );
            else
                return static_cast<uint16_type>( 4 * ptsPerVertex + 4 * ptsPerEdge + ptsPerFace );
        case 3:
            if constexpr ( convex_type::is_simplex )
                return static_cast<uint16_type>( 4 * ptsPerVertex + 6 * ptsPerEdge + 4 * ptsPerFace + ptsPerVolume );
            else
                return static_cast<uint16_type>( 8 * ptsPerVertex + 12 * ptsPerEdge + 6 * ptsPerFace + ptsPerVolume );
        default:
            return 0;
        }
    }

    /** @brief Validate primal/tag agreement and select static or runtime construction. */
    static pointset_type makePointSet(primal_space_type const& primal,RuntimeOrder order)
    {
        (void)validatedOrder<primal_space_type>(order);
        if(primal.order()!=order.value)
            throw std::invalid_argument("point-set order disagrees with the primal degree");
        if constexpr(super::is_order_dynamic) return pointset_type{order};
        else return pointset_type{};
    }

    template<typename PointSetT>
    void cacheSubEntityPoints( PointSetT const& pset )
    {
        M_points_by_subentity.clear();
        M_points_by_subentity.resize( nDim + 1 );
        for ( uint16_type d = 0; d < nDim; ++d )
        {
            auto const range = M_convex_ref.entityRange( d );
            M_points_by_subentity[d].resize( range.end() );
            auto const nSubEntityPoints = subEntityPointCount( pset, d );
            for ( uint16_type e = range.begin(); e < range.end(); ++e )
            {
                if ( nSubEntityPoints == 0 )
                    M_points_by_subentity[d][e].resize( nDim, 0, false );
                else
                    M_points_by_subentity[d][e] = pset.pointsBySubEntity( d, e, 1 );
            }
        }
    }

    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_points_face;
    std::vector<std::vector<points_type> > M_points_by_subentity;
    FunctionalSet<primal_space_type> M_fset;
    /** Interpolation nodes with the same static/runtime policy as the primal. */
    pointset_type M_pset;

};
}// details
/// \endcond detail

    class LagrangePolynomialSet {};
/**
 * \class Lagrange
 * \brief Lagrange polynomial set
 *
 * The \p Lagrange polynomial set is parametrized by
 *
 * -# dimension of the geometrical space
 * -# order of the Lagrange polynomials
 * -# the numerical type
 * -# the geometry it applies to (convexes such as simplices or product of simplices)
 *
 * Supports both static order (compile-time) and dynamic order (runtime).
 * For dynamic order, use O = Dynamic and provide RuntimeOrder at construction.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type N,
         uint16_type RealDim,
         int O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType = Continuous,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         template<class, int, class> class Pts = PointSetFekete,
         uint16_type TheTAG = 0 >
class Lagrange
    :
    public LagrangePolynomialSet,
    public FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, O, RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts >
{
    //! Effective order for base class (use placeholder 1 for dynamic)
    /** Authoritative signed polynomial order. */
    static constexpr int order_spec = O;
    typedef FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, O, RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts > super;
public:

    BOOST_STATIC_ASSERT( ( boost::is_same<PolySetType<N>, Scalar<N> >::value ||
                           boost::is_same<PolySetType<N>, Vectorial<N> >::value ||
                           boost::is_same<PolySetType<N>, Tensor2<N> >::value ||
                           boost::is_same<PolySetType<N>, Tensor2Symm<N> >::value ) );

    /** @name Typedefs
     */
    //@{

    //! @brief True if Order is known at compile time
    static constexpr bool is_order_static = (O >= 0);
    //! @brief True if Order is determined at runtime
    static constexpr bool is_order_dynamic = !is_order_static;
    //! @brief Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = O;

    inline static const uint16_type nDim = N;
    inline static const uint16_type nRealDim = RealDim;
    //! Static order (or placeholder 1 for dynamic case)
    inline static const uint16_type nOrder = is_order_static ? static_cast<uint16_type>(O) : 1;
    static constexpr bool isTransformationEquivalent = true;
    static constexpr bool isContinuous = ContinuityType::is_continuous;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef ContinuityType continuity_type;
    inline static const uint16_type TAG = TheTAG;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    inline static constexpr bool is_symm_v  = Feel::is_symm_v<polyset_type>;
    using is_symm  = Feel::is_symm<polyset_type>;
    static constexpr bool is_tensor2 = polyset_type::is_tensor2;
    static constexpr bool is_tensor2symm = is_tensor2 && is_symm_v;
    static constexpr bool is_vectorial = polyset_type::is_vectorial;
    static constexpr bool is_scalar = polyset_type::is_scalar;
    inline static const uint16_type nComponents = polyset_type::nComponents;
    inline static const uint16_type nComponents1 = polyset_type::nComponents1;
    inline static const uint16_type nComponents2 = polyset_type::nComponents2;

    static inline const bool is_product = true;
    static constexpr int Nm2 = (N>2)?N-2:0;
    static constexpr int Nm1 = (N>0)?N-1:0;

    typedef Lagrange<N, RealDim, O, PolySetType, ContinuityType, T, Convex,  Pts, TheTAG> this_type;
    typedef Lagrange<N, RealDim, O, Scalar, continuity_type, T, Convex,  Pts, TheTAG> component_basis_type;


    using face_basis_type = if_t<nDim <= 1,
                                 boost::none_t,
                                 Lagrange<Nm1, RealDim, O, Scalar, continuity_type, T, Convex, Pts, TheTAG>>;

    typedef std::shared_ptr<face_basis_type> face_basis_ptrtype;
    using edge_basis_type = if_t<nDim <= 2,
                                 boost::none_t,
                                 Lagrange<Nm2, RealDim, O, Scalar, continuity_type, T, Convex, Pts, TheTAG>>;

    typedef std::shared_ptr<edge_basis_type> edge_basis_ptrtype;

    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;
    typedef typename convex_type::edge_type edge_type;

    //! @name Static DOF counts (compile-time, use runtime methods for dynamic order)
    //! @{
    inline static const uint16_type numPoints = reference_convex_type::numPoints;
    inline static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    inline static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    inline static const uint16_type nbPtsPerFace = reference_convex_type::nbPtsPerFace;
    inline static const uint16_type nbPtsPerVolume = reference_convex_type::nbPtsPerVolume;
    inline static const uint16_type nLocalDof = dual_space_type::nLocalDof;
    inline static const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    inline static const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    inline static const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    inline static const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    inline static const uint16_type nLocalFaceDof = ( face_type::numVertices * nDofPerVertex +
                                               face_type::numEdges * nDofPerEdge +
                                               face_type::numFaces * nDofPerFace );
    inline static const uint16_type nLocalEdgeDof = ( edge_type::numVertices * nDofPerVertex +
                                                      edge_type::numEdges * nDofPerEdge);
    inline static const uint16_type nLocalVertexDof = nDofPerVertex;
    //! @}

    //! @name DOF accessors (unified interface using C++20 requires clauses)
    //! @{

    /**
     * @brief Get polynomial order (unified interface)
     *
     * For static order: constexpr, returns compile-time O
     * For dynamic order: runtime semantic value provided by the primal space.
     */
    [[nodiscard]] constexpr uint16_type order() const noexcept
        requires is_static_order<O>
    {
        return static_cast<uint16_type>( O );
    }

    [[nodiscard]] uint16_type order() const noexcept
        requires is_dynamic_order<O>
    {
        return static_cast<uint16_type>( super::order() );
    }

    /**
     * @brief Get total number of local DOFs (unified interface)
     */
    [[nodiscard]] constexpr uint16_type localDof() const noexcept
        requires is_static_order<O>
    {
        return nLocalDof;
    }

    [[nodiscard]] uint16_type localDof() const
        requires is_dynamic_order<O>
    {
        const auto runtime_order = this->order();
        if constexpr ( convex_type::is_simplex )
            return ::Feel::detail::checkedLocalIndexValue(
                ::Feel::detail::simplexTotal( nDim, runtime_order ) );
        else
            return ::Feel::detail::checkedLocalIndexValue(
                ::Feel::detail::hypercubeTotal( nDim, runtime_order ) );
    }

    /**
     * @brief Get DOFs per vertex (unified interface)
     */
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept
        requires is_static_order<O>
    {
        return nDofPerVertex;
    }

    [[nodiscard]] uint16_type dofPerVertex() const noexcept
        requires is_dynamic_order<O>
    {
        if constexpr ( convex_type::is_simplex )
            return ::Feel::detail::simplexPerVertex( nDim, this->order() );
        else
            return ::Feel::detail::hypercubePerVertex( nDim, this->order() );
    }

    /**
     * @brief Get DOFs per edge (unified interface)
     */
    [[nodiscard]] constexpr uint16_type dofPerEdge() const noexcept
        requires is_static_order<O>
    {
        return nDofPerEdge;
    }

    [[nodiscard]] uint16_type dofPerEdge() const noexcept
        requires is_dynamic_order<O>
    {
        if constexpr ( convex_type::is_simplex )
            return ::Feel::detail::simplexPerEdge( nDim, this->order() );
        else
            return ::Feel::detail::hypercubePerEdge( nDim, this->order() );
    }

    /**
     * @brief Get DOFs per face (unified interface)
     */
    [[nodiscard]] constexpr uint16_type dofPerFace() const noexcept
        requires is_static_order<O>
    {
        return nDofPerFace;
    }

    [[nodiscard]] uint16_type dofPerFace() const noexcept
        requires is_dynamic_order<O>
    {
        if constexpr ( convex_type::is_simplex )
            return ::Feel::detail::simplexPerFace( nDim, this->order() );
        else
            return ::Feel::detail::hypercubePerFace( nDim, this->order() );
    }

    /**
     * @brief Get DOFs per volume (unified interface)
     */
    [[nodiscard]] constexpr uint16_type dofPerVolume() const noexcept
        requires is_static_order<O>
    {
        return nDofPerVolume;
    }

    [[nodiscard]] uint16_type dofPerVolume() const noexcept
        requires is_dynamic_order<O>
    {
        if constexpr ( convex_type::is_simplex )
            return ::Feel::detail::simplexPerVolume( nDim, this->order() );
        else
            return ::Feel::detail::hypercubePerVolume( nDim, this->order() );
    }

    //! @}
    template<int subN>
    struct SubSpace
    {
        typedef Lagrange<N-1, RealDim, O, PolySetType, continuity_type, T, Convex,  Pts, TheTAG> type;
    };

    struct SSpace
    {
        //! For dynamic order, keep dynamic; otherwise compute reduced order
        static constexpr int TheOrder = is_order_dynamic ? Dynamic : ((O > 1) ? O-1 : 0);
        using type = if_t<is_order_dynamic,
                          Lagrange<N, RealDim, Dynamic, PolySetType, continuity_type, T, Convex, Pts, TheTAG>,
                          if_t<O <= 1,
                               Lagrange<N, RealDim, 0, PolySetType, Discontinuous, T, Convex, Pts, TheTAG>,
                               Lagrange<N, RealDim, TheOrder, PolySetType, continuity_type, T, Convex, Pts, TheTAG>>>;
    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef Lagrange<NewDim, RealDim, O, PolySetType, continuity_type, T, Convex,  Pts, TheTAG> type;
    };

    static constexpr bool isLagrangeP0Continuous = isP0Continuous<this_type>::result;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * @brief Default constructor
     *
     * For static order, uses the compile-time order.
     * For dynamic order, defaults to order 1 (use RuntimeOrder constructor instead).
     */
    Lagrange()
        :
        super( buildDualSpace( RuntimeOrder{O == Dynamic ? 1 : O} ) ),
        M_refconvex()
    {
        initSymmetricMapping();
    }

    /**
     * @brief Constructor with runtime order
     *
     * For dynamic order types, this constructor specifies the actual polynomial order.
     * For static order types, it validates that the supplied value agrees with
     * the compile-time order.
     *
     * @param order Runtime order specification
     */
    explicit Lagrange( RuntimeOrder order )
        :
        super( buildDualSpace( order ) ),
        M_refconvex()
    {
        initSymmetricMapping();
    }

private:
    /**
     * @brief Helper to build the dual space with appropriate order
     *
     * For dynamic order, creates dual space with runtime-order-aware constructor.
     * For static order, uses the default constructor.
     */
    static dual_space_type buildDualSpace( RuntimeOrder order )
    {
        order = validatedOrder<this_type>( order );
        if constexpr ( is_order_dynamic )
        {
            if constexpr ( convex_type::is_simplex )
                order = ::Feel::detail::validatedNodalOrder<
                    ::Feel::detail::SimplexCardinality>( nDim, order,
                                                         nComponents );
            else
                order = ::Feel::detail::validatedNodalOrder<
                    ::Feel::detail::HypercubeCardinality>( nDim, order,
                                                           nComponents );

            // Use runtime order for both primal and dual spaces
            primal_space_type primal( order );
            return dual_space_type( primal, order );
        }
        else
        {
            // Static order was validated above; use compile-time construction.
            return dual_space_type( primal_space_type() );
        }
    }

public:

    ~Lagrange() override {}

private:
    /**
     * @brief Initialize the symmetric tensor mapping (called by constructors)
     */
    void initSymmetricMapping()
    {
        if ( is_tensor2symm )
        {
            const auto localDof = this->localDof();
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

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the reference convex associated with the lagrange polynomials
     */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const override
    {
        return "lagrange";
    }

    /** @brief Return the effective cardinality of one field component. */
    uint16_type localDofPerComponent() const override
    {
        return localDof();
    }

    typename super::DofAttachment dofAttachment( uint16_type localDofId ) const override
    {
        const uint16_type parentLocalDofId = this->dofParent( localDofId );
        const uint16_type nV = static_cast<uint16_type>( reference_convex_type::numVertices * dofPerVertex() );
        const uint16_type nE = static_cast<uint16_type>( reference_convex_type::numEdges * dofPerEdge() );
        const uint16_type nF = static_cast<uint16_type>( reference_convex_type::numFaces * dofPerFace() );

        if ( parentLocalDofId < nV && dofPerVertex() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 0,
                .entityId = static_cast<uint16_type>( parentLocalDofId / dofPerVertex() ),
                .ordinal = static_cast<uint16_type>( parentLocalDofId % dofPerVertex() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterVertex = static_cast<uint16_type>( parentLocalDofId - nV );
        if ( parentAfterVertex < nE && dofPerEdge() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 1,
                .entityId = static_cast<uint16_type>( parentAfterVertex / dofPerEdge() ),
                .ordinal = static_cast<uint16_type>( parentAfterVertex % dofPerEdge() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterEdge = static_cast<uint16_type>( parentAfterVertex - nE );
        if ( parentAfterEdge < nF && dofPerFace() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 2,
                .entityId = static_cast<uint16_type>( parentAfterEdge / dofPerFace() ),
                .ordinal = static_cast<uint16_type>( parentAfterEdge % dofPerFace() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type nDofPerVolume = dofPerVolume();
        if ( nDofPerVolume > 0 )
        {
            const uint16_type parentAfterFace = static_cast<uint16_type>( parentAfterEdge - nF );
            return typename super::DofAttachment{
                .entityDim = 3,
                .entityId = 0,
                .ordinal = static_cast<uint16_type>( parentAfterFace % nDofPerVolume ),
                .kind = this->dofType( localDofId ) };
        }

        return typename super::DofAttachment{
            .entityDim = -1,
            .entityId = super::DofAttachment::invalid_id,
            .ordinal = super::DofAttachment::invalid_id,
            .kind = this->dofType( localDofId ) };
    }


    //! \return the component of a local dof
    uint16_type component( uint16_type localDofId ) const override
        {
            const uint16_type localDofPerComp = localDof();
            uint16_type comp = localDofId/localDofPerComp;
            DCHECK( comp < nComponents ) << "invalid localDofId " << localDofId;
            return comp;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const override
        {
            const uint16_type localDofPerComp = localDof();
            uint16_type ldofParent = localDofId % localDofPerComp;
            return ldofParent;
        }

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const override
        {
            return 1;
        }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const override
        {
            if ( !is_tensor2symm )
                return i;
            DCHECK( M_unsymm2symm.size() > i ) << "invalid size of unsymm2symm container";
            return M_unsymm2symm[i];
        }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //!
    //! build interpolant object on \p n finite elements
    //!
    using  local_interpolant_type = Eigen::VectorXd;
    local_interpolant_type
    localInterpolant( int n = 1 ) const
        {
            // Use the effective local DOF count for both order modes.
            return local_interpolant_type::Zero( n*nComponents*localDof() );
        }

    //!
    //! build \p p  interpolants object on \p n finite elements
    //!
    using  local_interpolants_type = Eigen::MatrixXd;
    local_interpolants_type
    localInterpolants( int p, int n = 1 ) const
        {
            // Use the effective local DOF count for both order modes.
            return local_interpolants_type::Zero( n*nComponents*localDof(), p );
        }

    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            // Use the effective local DOF count for both order modes.
            const int rtLocalDof = localDof();
            for( int q = 0; q < rtLocalDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( (c2+nComponents2*c1)*rtLocalDof+q ) = expr.evalq( c1, c2, q );
                            Ihloc( (c1+nComponents2*c2)*rtLocalDof+q ) = Ihloc( (c2+nComponents2*c1)*rtLocalDof+q );
                        }
                        // diagonal
                        Ihloc( (c1+nComponents2*c1)*rtLocalDof+q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( (c2+nComponents2*c1)*rtLocalDof+q ) = expr.evalq( c1, c2, q );
                    }
                }
        }
    local_interpolant_type
    faceLocalInterpolant() const
        {
            const int rtLocalFaceDof =
                face_type::numVertices * dofPerVertex() +
                face_type::numEdges * dofPerEdge() +
                face_type::numFaces * dofPerFace();
            return local_interpolant_type::Zero( nComponents * rtLocalFaceDof, 1 );
        }
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalFaceDof =
                face_type::numVertices * dofPerVertex() +
                face_type::numEdges * dofPerEdge() +
                face_type::numFaces * dofPerFace();
            for( int q = 0; q < rtLocalFaceDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalFaceDof + q ) = expr.evalq( c1, c2, q );
                            Ihloc( ( c1 + nComponents2 * c2 ) * rtLocalFaceDof + q ) = Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalFaceDof + q );
                        }
                        Ihloc( ( c1 + nComponents2 * c1 ) * rtLocalFaceDof + q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalFaceDof + q ) = expr.evalq( c1, c2, q );
                    }
                }
        }

    local_interpolant_type
    edgeLocalInterpolant() const
        {
            const int rtLocalEdgeDof =
                edge_type::numVertices * dofPerVertex() +
                edge_type::numEdges * dofPerEdge();
            return local_interpolant_type::Zero( nComponents * rtLocalEdgeDof, 1 );
        }
    template<typename ExprType>
    void
    edgeInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalEdgeDof =
                edge_type::numVertices * dofPerVertex() +
                edge_type::numEdges * dofPerEdge();
            for( int q = 0; q < rtLocalEdgeDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalEdgeDof + q ) = expr.evalq( c1, c2, q );
                            Ihloc( ( c1 + nComponents2 * c2 ) * rtLocalEdgeDof + q ) = Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalEdgeDof + q );
                        }
                        Ihloc( ( c1 + nComponents2 * c1 ) * rtLocalEdgeDof + q ) = expr.evalq( c1, c1, q );

                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalEdgeDof + q ) = expr.evalq( c1, c2, q );
                    }
                }
        }
    local_interpolant_type
    vertexLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nComponents * dofPerVertex(), 1 );
        }
    template<typename ExprType>
    void
    vertexInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalVertexDof = dofPerVertex();
            for( int q = 0; q < rtLocalVertexDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalVertexDof + q ) = expr.evalq( c1, c2, q );
                            Ihloc( ( c1 + nComponents2 * c2 ) * rtLocalVertexDof + q ) = Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalVertexDof + q );
                        }
                        Ihloc( ( c1 + nComponents2 * c1 ) * rtLocalVertexDof + q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( ( c2 + nComponents2 * c1 ) * rtLocalVertexDof + q ) = expr.evalq( c1, c2, q );
                    }

                }
        }
    template<typename ExprType>
    void
    evaluateBasisFunction/*interpolateBasisFunction*/( ExprType&& expr, local_interpolants_type & Ihloc ) const
    {
        using shape = typename std::decay_t<ExprType>::shape;
        //for ( int cc1 = 0; cc1 < nComponents1; ++cc1 )
        using expr_basis_t = typename std::decay_t<ExprType>::expr_type::test_basis;

        const int nPoints = expr.nPoints();
        const int ncomp1 = ( expr_basis_t::is_product ? expr_basis_t::nComponents1 : 1 );
        CHECK( ncomp1 > 0 );
        CHECK( Ihloc.rows() % ncomp1 == 0 ) << Ihloc.rows() << " % " << ncomp1;
        const int exprLocalDofPerComp = Ihloc.rows() / ncomp1;
        for( int q = 0; q < nPoints; ++q )
        {
            for( int i = 0; i < exprLocalDofPerComp; ++i )
            {
                for ( uint16_type c = 0; c < ncomp1; ++c )
                {
                    const uint16_type I = static_cast<uint16_type>( exprLocalDofPerComp * c + i );
                    for( int c1 = 0; c1 < shape::M; ++c1 )
                    {
                        if ( is_symm_v )
                        {
                            for( int c2 = 0; c2 < c1; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*nPoints + q;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                                int ldof2 = (c1+nComponents2*c2)*nPoints + q;
                                Ihloc( I, ldof2) = Ihloc( I, ldof);
                            }
                            int ldof = (c1+nComponents2*c1)*nPoints + q;
                            Ihloc( I, ldof) = expr.evaliq( I, c1, c1, q );
                        }
                        else
                        {
                            for( int c2 = 0; c2 < shape::N; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*nPoints + q;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                            }
                        }
                    }
                }
            }
        }
    }

    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc, std::vector<uint16_type> const& mapExprPointToDofPoint ) const
    {
        using shape = typename std::decay_t<ExprType>::shape;
        //for ( int cc1 = 0; cc1 < nComponents1; ++cc1 )
        using expr_basis_t = typename std::decay_t<ExprType>::expr_type::test_basis;

        const int nPoints = expr.nPoints();
        const int ncomp1 = ( expr_basis_t::is_product ? expr_basis_t::nComponents1 : 1 );
        CHECK( ncomp1 > 0 );
        CHECK( Ihloc.rows() % ncomp1 == 0 ) << Ihloc.rows() << " % " << ncomp1;
        const int exprLocalDofPerComp = Ihloc.rows() / ncomp1;
        const int localDofPerComp = localDof();
        for( int q = 0; q < nPoints; ++q )
        {
            uint16_type q2 = mapExprPointToDofPoint[q];
            for( int i = 0; i < exprLocalDofPerComp; ++i )
            {
                for ( uint16_type c = 0; c < ncomp1; ++c )
                {
                    const uint16_type I = static_cast<uint16_type>( exprLocalDofPerComp * c + i );
                    for( int c1 = 0; c1 < shape::M; ++c1 )
                    {
                        if ( is_symm_v )
                        {
                            for( int c2 = 0; c2 < c1; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*localDofPerComp + q2;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                                int ldof2 = (c1+nComponents2*c2)*localDofPerComp + q2;
                                Ihloc( I, ldof2) = Ihloc( I, ldof);
                            }
                            int ldof = (c1+nComponents2*c1)*localDofPerComp + q2;
                            Ihloc( I, ldof) = expr.evaliq( I, c1, c1, q );
                        }
                        else
                        {
                            for( int c2 = 0; c2 < shape::N; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*localDofPerComp + q2;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                            }
                        }
                    }
                }
            }
        }
    }

    template<typename ExprType,
             std::enable_if_t< std::decay_t<ExprType>::gmc_type::subEntityCoDim == 0 ,bool> = true >
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc ) const
        {
            std::vector<uint16_type> mapExprPointToDofPoint( localDof() );
            std::iota( mapExprPointToDofPoint.begin(), mapExprPointToDofPoint.end(), 0 );
            this->interpolateBasisFunction( expr, Ihloc, mapExprPointToDofPoint );
        }

    template<typename ExprType,
             std::enable_if_t< std::decay_t<ExprType>::gmc_type::subEntityCoDim == 1 ,bool> = true >
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc ) const
        {
            int nPoints = expr.nPoints();
            const uint16_type rtDofPerVertex = dofPerVertex();
            const uint16_type rtDofPerEdge = dofPerEdge();
            const uint16_type rtDofPerFace = dofPerFace();
            const uint16_type rtLocalFacetDof =
                static_cast<uint16_type>( face_type::numVertices * rtDofPerVertex +
                                          face_type::numEdges * rtDofPerEdge +
                                          face_type::numFaces * rtDofPerFace );
            CHECK( rtLocalFacetDof == nPoints ) << rtLocalFacetDof << " vs "<< nPoints;
            auto gmc = expr.geom();
            uint16_type facetIdInElt = gmc->faceId();
            std::vector<uint16_type> mapExprPointToDofPoint( rtLocalFacetDof );
            for ( int lfd = 0;lfd < rtLocalFacetDof;++lfd )
            {
                if ( reference_convex_type::nDim == 2 )
                {
                    const uint16_type vertexBlock = static_cast<uint16_type>( face_type::numVertices * rtDofPerVertex );
                    if ( lfd < vertexBlock )
                    {
                        const uint16_type facetVertex = static_cast<uint16_type>( lfd / rtDofPerVertex );
                        const uint16_type vertexOrdinal = static_cast<uint16_type>( lfd % rtDofPerVertex );
                        const uint16_type elementVertex = reference_convex_type::e2p( facetIdInElt, facetVertex );
                        mapExprPointToDofPoint[lfd] = static_cast<uint16_type>( elementVertex * rtDofPerVertex + vertexOrdinal );
                    }
                    else if ( lfd < vertexBlock + face_type::numEdges * rtDofPerEdge )
                    {
                        const uint16_type edgeOrdinal = static_cast<uint16_type>( ( lfd - vertexBlock ) % rtDofPerEdge );
                        mapExprPointToDofPoint[lfd] =
                            static_cast<uint16_type>( reference_convex_type::numVertices * rtDofPerVertex +
                                                      facetIdInElt * rtDofPerEdge +
                                                      edgeOrdinal );
                    }
                    else
                        CHECK( false ) << "TODO";
                }
                else
                    CHECK( false ) << "TODO";
            }

            this->interpolateBasisFunction( expr, Ihloc, mapExprPointToDofPoint );
        }

    //@}

private:

    reference_convex_type M_refconvex;
    face_basis_ptrtype M_bdylag;
    std::vector<uint16_type> M_unsymm2symm;
};
} // namespace fem

/**
 * @class Lagrange
 * @brief Lagrange finite element factory
 *
 * Supports both static order (default, backward compatible) and dynamic order.
 * For dynamic order, use Lagrange<Dynamic, PolySetType>.
 *
 * @tparam Order Polynomial order (can be Dynamic for runtime order)
 * @tparam PolySetType Field type (Scalar, Vectorial, Tensor2)
 * @tparam ContinuityType Continuous or Discontinuous
 * @tparam Pts Point set type for interpolation
 * @tparam TheTAG Tag for multiple instances
 */
template<int Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         typename ContinuityType = Continuous,
         template<class, int, class> class Pts = PointSetFekete,
         uint16_type TheTAG = 0>
class Lagrange : public FiniteElementFamilyDescriptor<Order, TheTAG>
{
public:
    /**
     * @brief Apply metafunction to get concrete finite element type
     *
     * For static order, returns fem::Lagrange<N, RealDim, Order, ...>
     * For dynamic order, requires runtime order at construction
     */
    template<uint16_type N,
             uint16_type RealDim,
             typename T = double,
             typename Convex = Simplex<N>>
    struct Apply
    {
        // Preserve the Order template parameter (including Dynamic = -1)
        // For dynamic order, the actual polynomial degree is set at construction via RuntimeOrder
        using result_type = if_t<Convex::is_simplex,
                                 fem::Lagrange<N, RealDim, Order, PolySetType, ContinuityType, T, Simplex, Pts, TheTAG>,
                                 fem::Lagrange<N, RealDim, Order, PolySetType, ContinuityType, T, Hypercube, Pts, TheTAG>>;
        typedef result_type type;
    };

    /** Concrete-element application retained for the established descriptor protocol. */
    template<uint16_type N, uint16_type RealDim, typename T = double, typename Convex = Simplex<N>>
    using apply = Apply<N, RealDim, T, Convex>;

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef Lagrange<Order, PolySetType, ContinuityType, Pts, TheNewTAG> type;
    };

    typedef Lagrange<Order, Scalar, ContinuityType, Pts, TheTAG> component_basis_type;
};



template<typename P>
using is_lagrange_polynomialset = std::is_base_of<fem::LagrangePolynomialSet,P>;
template<typename P>
constexpr bool is_lagrange_polynomialset_v = boost::is_base_of<fem::LagrangePolynomialSet,P>::value;

template<typename P>
constexpr bool is_lagrange_polynomialset_P0d_v = is_lagrange_polynomialset_v<P> && (P::nOrder == 0) && !P::continuity_type::is_continuous;

} // namespace Feel
#endif /* __lagrange_H */
