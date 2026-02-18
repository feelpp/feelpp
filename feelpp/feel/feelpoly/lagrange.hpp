/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008-2012 Université de Grenoble 1 (Joseph Fourier)

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
#include <feel/feelpoly/functionals.hpp>
#include <feel/feelpoly/fe.hpp>
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
template<template<class, int, class> class PointSetT, typename ConvexType, typename ValueType>
concept RuntimePointSet = requires( RuntimeOrder ro,
                                    PointSetT<ConvexType, Dynamic, ValueType> dyn_pset )
{
    { PointSetT<ConvexType, Dynamic, ValueType>( ro ) };
    { dyn_pset.runtimeNumPoints() } -> std::convertible_to<uint32_type>;
    dyn_pset.points();
    dyn_pset.pointsBySubEntity( uint16_type{}, uint16_type{}, uint16_type{} );
};

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

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    template< template<class, int, class> class TestPointSetType >
    static constexpr bool is_pointset_v = std::is_base_of_v<TestPointSetType<convex_type, nOrder, value_type>,pointset_type >;
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

    LagrangeDual( LagrangeDual const& d )
        :
        super( d ),
        M_convex_ref(),
        M_eid( d.M_eid ),
        M_pts( d.M_pts ),
        M_points_face( d.M_points_face ),
        M_fset( d.M_fset )
        {}

    LagrangeDual( LagrangeDual && d ) = default;

    /**
     * @brief Constructor for dynamic order support
     *
     * Uses PointSetEquiSpaced with runtime order to create proper interpolation points.
     * This constructor enables dynamic polynomial order selection at runtime.
     *
     * @param primal The primal space (OrthonormalPolynomialSet with matching runtime order)
     * @param order Runtime order specification
     */
    LagrangeDual( primal_space_type const& primal, RuntimeOrder order )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts(),  // Will be set from dynamic pointset
        M_points_face( nFacesInConvex ),
        M_fset( primal )
    {
        DVLOG(2) << "Lagrange finite element (dynamic order): \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << order.value << " (runtime)\n";

        auto initFromPointSet = [&]( auto const& dyn_pset )
        {
            const uint32_type runtime_numPoints = dyn_pset.runtimeNumPoints();
            DVLOG(2) << " o- numPoints = " << runtime_numPoints << " (runtime)\n";

            M_pts.resize( nDim, runtime_numPoints );
            M_pts = dyn_pset.points();

            if ( order.value > 0 )
            {
                for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                        e < M_convex_ref.entityRange( nDim-1 ).end();
                        ++e )
                {
                    M_points_face[e] = dyn_pset.pointsBySubEntity( nDim-1, e, 1 );
                    DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
                }
            }
        };

        using dynamic_pointset_type = PointSetType<convex_type, Dynamic, value_type>;
        constexpr bool has_runtime_pointset =
            RuntimePointSet<PointSetType, convex_type, value_type>;

        if constexpr ( has_runtime_pointset )
        {
            dynamic_pointset_type dyn_pset( order );
            initFromPointSet( dyn_pset );
        }
        else
        {
            PointSetEquiSpaced<convex_type, Dynamic, value_type> dyn_pset( order );
            initFromPointSet( dyn_pset );
        }

        setFset( primal, M_pts, bool_c<primal_space_type::is_scalar>{} );
    }

    LagrangeDual( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_points_face( nFacesInConvex ),
        M_fset( primal )
    {
        DVLOG(2) << "Lagrange finite element: \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << nOrder << "\n";
        DVLOG(2) << " o- numPoints      = " << numPoints << "\n";
        DVLOG(2) << " o- nbPtsPerVertex = " << nbPtsPerVertex << "\n";
        DVLOG(2) << " o- nbPtsPerEdge   = " << nbPtsPerEdge << "\n";
        DVLOG(2) << " o- nbPtsPerFace   = " << nbPtsPerFace << "\n";
        DVLOG(2) << " o- nbPtsPerVolume = " << nbPtsPerVolume << "\n";

        M_pts = M_pset.points();

        if constexpr ( nOrder > 0 )
        {
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                    e < M_convex_ref.entityRange( nDim-1 ).end();
                    ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
                DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
            }
        }

        setFset( primal, M_pts, bool_c<primal_space_type::is_scalar>{} );
    }

    LagrangeDual( primal_space_type const& primal, pointset_type const& pts )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( pts.points() ),
        M_points_face( nFacesInConvex ),
        M_fset( primal ),
        M_pset( pts )
    {
        DVLOG(2) << "Lagrange finite element: \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << nOrder << "\n";
        DVLOG(2) << " o- numPoints      = " << numPoints << "\n";
        DVLOG(2) << " o- nbPtsPerVertex = " << nbPtsPerVertex << "\n";
        DVLOG(2) << " o- nbPtsPerEdge   = " << nbPtsPerEdge << "\n";
        DVLOG(2) << " o- nbPtsPerFace   = " << nbPtsPerFace << "\n";
        DVLOG(2) << " o- nbPtsPerVolume = " << nbPtsPerVolume << "\n";

        if constexpr ( nOrder > 0 )
        {
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                    e < M_convex_ref.entityRange( nDim-1 ).end();
                    ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
                DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
            }
        }

        setFset( primal, M_pts, bool_c<primal_space_type::is_scalar>{} );
    }

    ~LagrangeDual() = default;
    LagrangeDual& operator=( LagrangeDual const& ) = default;

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
                return M_pset.pointsBySubEntity( topodim, entity, 1 );
            }


    points_type edgePoints(int edge) const
        {
            return M_pset.pointsBySubEntity( 1, edge, 1 );
        }


    points_type vertexPoints(int vertex) const
        {
            return M_pset.pointsBySubEntity( 0, vertex, 1 );
        }

    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }
private:

    void setFset( primal_space_type const& primal, points_type const& __pts, bool_c<true> )
    {
        M_fset.setFunctionalSet( functional::PointsEvaluation<primal_space_type>( primal,
                                  __pts ) );
    }

    void setFset( primal_space_type const& primal, points_type const& __pts, bool_c<false> )
    {
        M_fset.setFunctionalSet( functional::ComponentsPointsEvaluation<primal_space_type>( primal,
                                  __pts ) );
    }

    /**
     * set the pointset at face \c f using points \c n
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_points_face[f].resize( n.size1(), n.size2(), false );
        M_points_face[f] = n;
    }

private:
    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_points_face;
    FunctionalSet<primal_space_type> M_fset;
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
    public FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, (O >= 0 ? O : 1), RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts >
{
    //! Effective order for base class (use placeholder 1 for dynamic)
    static constexpr int EffectiveOrder = (O >= 0) ? O : 1;
    typedef FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, EffectiveOrder, RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts > super;
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
     * @deprecated Use order() instead
     */
    [[nodiscard]] uint16_type runtimeOrder() const noexcept
    {
        return order();
    }

    /**
     * @brief Get total number of local DOFs (unified interface)
     */
    [[nodiscard]] constexpr uint16_type localDof() const noexcept
        requires is_static_order<O>
    {
        return nLocalDof;
    }

    [[nodiscard]] uint16_type localDof() const noexcept
        requires is_dynamic_order<O>
    {
        const auto runtime_order = this->order();
        return static_cast<uint16_type>(
            ::Feel::detail::simplexTotal( nDim, runtime_order ) );
    }

    /**
     * @deprecated Use localDof() instead
     */
    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
    {
        return localDof();
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
        return ::Feel::detail::simplexPerVolume( nDim, this->order() );
    }

    /**
     * @deprecated Use dofPerVolume() instead
     */
    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
    {
        return dofPerVolume();
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
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {
        initSymmetricMapping();
    }

    /**
     * @brief Constructor with runtime order
     *
     * For dynamic order types, this constructor specifies the actual polynomial order.
     * For static order types, the RuntimeOrder is ignored (compile-time order is used).
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
        if constexpr ( is_order_dynamic )
        {
            // Use runtime order for both primal and dual spaces
            primal_space_type primal( order );
            return dual_space_type( primal, order );
        }
        else
        {
            // Static order: ignore RuntimeOrder and use compile-time order
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
            const auto localDof = runtimeLocalDof();
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

    uint16_type localDofPerComponent() const override
    {
        return runtimeLocalDof();
    }

    uint16_type localDofId( uint16_type parentLocalDofId, uint16_type component = 0 ) const override
    {
        FEELPP_ASSERT( component < nComponents )
            ( component )( parentLocalDofId )( nComponents ).error( "invalid component index" );
        const uint16_type localDofPerComp = runtimeLocalDof();
        FEELPP_ASSERT( parentLocalDofId < localDofPerComp )
            ( parentLocalDofId )( localDofPerComp ).error( "invalid parent local dof index" );
        return static_cast<uint16_type>( localDofPerComp * component + parentLocalDofId );
    }

    typename super::DofAttachment dofAttachment( uint16_type localDofId ) const override
    {
        const uint16_type parentLocalDofId = this->dofParent( localDofId );
        const uint16_type nV = static_cast<uint16_type>( reference_convex_type::numVertices * runtimeDofPerVertex() );
        const uint16_type nE = static_cast<uint16_type>( reference_convex_type::numEdges * runtimeDofPerEdge() );
        const uint16_type nF = static_cast<uint16_type>( reference_convex_type::numFaces * runtimeDofPerFace() );

        if ( parentLocalDofId < nV && runtimeDofPerVertex() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 0,
                .entityId = static_cast<uint16_type>( parentLocalDofId / runtimeDofPerVertex() ),
                .ordinal = static_cast<uint16_type>( parentLocalDofId % runtimeDofPerVertex() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterVertex = static_cast<uint16_type>( parentLocalDofId - nV );
        if ( parentAfterVertex < nE && runtimeDofPerEdge() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 1,
                .entityId = static_cast<uint16_type>( parentAfterVertex / runtimeDofPerEdge() ),
                .ordinal = static_cast<uint16_type>( parentAfterVertex % runtimeDofPerEdge() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterEdge = static_cast<uint16_type>( parentAfterVertex - nE );
        if ( parentAfterEdge < nF && runtimeDofPerFace() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 2,
                .entityId = static_cast<uint16_type>( parentAfterEdge / runtimeDofPerFace() ),
                .ordinal = static_cast<uint16_type>( parentAfterEdge % runtimeDofPerFace() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type nDofPerVolume = runtimeDofPerVolume();
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
            const uint16_type localDofPerComp = runtimeLocalDof();
            uint16_type comp = localDofId/localDofPerComp;
            DCHECK( comp < nComponents ) << "invalid localDofId " << localDofId;
            return comp;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const override
        {
            const uint16_type localDofPerComp = runtimeLocalDof();
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
            // Use runtime local DOF count for dynamic order support
            return local_interpolant_type::Zero( n*nComponents*runtimeLocalDof() );
        }

    //!
    //! build \p p  interpolants object on \p n finite elements
    //!
    using  local_interpolants_type = Eigen::MatrixXd;
    local_interpolants_type
    localInterpolants( int p, int n = 1 ) const
        {
            // Use runtime local DOF count for dynamic order support
            return local_interpolants_type::Zero( n*nComponents*runtimeLocalDof(), p );
        }

    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            // Use runtime local DOF count for dynamic order support
            const int rtLocalDof = runtimeLocalDof();
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
                face_type::numVertices * runtimeDofPerVertex() +
                face_type::numEdges * runtimeDofPerEdge() +
                face_type::numFaces * runtimeDofPerFace();
            return local_interpolant_type::Zero( nComponents * rtLocalFaceDof, 1 );
        }
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalFaceDof =
                face_type::numVertices * runtimeDofPerVertex() +
                face_type::numEdges * runtimeDofPerEdge() +
                face_type::numFaces * runtimeDofPerFace();
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
                edge_type::numVertices * runtimeDofPerVertex() +
                edge_type::numEdges * runtimeDofPerEdge();
            return local_interpolant_type::Zero( nComponents * rtLocalEdgeDof, 1 );
        }
    template<typename ExprType>
    void
    edgeInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalEdgeDof =
                edge_type::numVertices * runtimeDofPerVertex() +
                edge_type::numEdges * runtimeDofPerEdge();
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
            return local_interpolant_type::Zero( nComponents * runtimeDofPerVertex(), 1 );
        }
    template<typename ExprType>
    void
    vertexInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            static_assert( nComponents1 == ExprType::shape::M, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            static_assert( nComponents2 == ExprType::shape::N, "INCOMPATIBLE_NUMBER_OF_COMPONENTS" );
            const int rtLocalVertexDof = runtimeDofPerVertex();
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
        const int localDofPerComp = runtimeLocalDof();
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
            std::vector<uint16_type> mapExprPointToDofPoint( runtimeLocalDof() );
            std::iota( mapExprPointToDofPoint.begin(), mapExprPointToDofPoint.end(), 0 );
            this->interpolateBasisFunction( expr, Ihloc, mapExprPointToDofPoint );
        }

    template<typename ExprType,
             std::enable_if_t< std::decay_t<ExprType>::gmc_type::subEntityCoDim == 1 ,bool> = true >
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc ) const
        {
            int nPoints = expr.nPoints();
            CHECK( nLocalFaceDof == nPoints ) << nLocalFaceDof << " vs "<< nPoints;
            auto gmc = expr.geom();
            uint16_type faceIdInElt = gmc->faceId();
            std::vector<uint16_type> mapExprPointToDofPoint( nLocalFaceDof );
            for ( int lfd = 0;lfd < nLocalFaceDof;++lfd )
            {
                if ( reference_convex_type::nDim == 2 )
                {
                    if ( lfd < (face_type::numVertices * nDofPerVertex) )
                        mapExprPointToDofPoint[lfd] = reference_convex_type::e2p(faceIdInElt,lfd);
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
class Lagrange
{
public:
    //! @name Order type detection
    //! @{
    static constexpr bool is_order_static = (Order != Dynamic);
    static constexpr bool is_order_dynamic = !is_order_static;
    //! @}

    //! Static order value (0 if dynamic, use runtime order in that case)
    inline static const uint16_type nOrder = is_order_static ? static_cast<uint16_type>(Order) : 0;
    inline static const uint16_type TAG = TheTAG;

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
    struct apply
    {
        // Preserve the Order template parameter (including Dynamic = -1)
        // For dynamic order, the actual polynomial degree is set at construction via RuntimeOrder
        using result_type = if_t<Convex::is_simplex,
                                 fem::Lagrange<N, RealDim, Order, PolySetType, ContinuityType, T, Simplex, Pts, TheTAG>,
                                 fem::Lagrange<N, RealDim, Order, PolySetType, ContinuityType, T, Hypercube, Pts, TheTAG>>;
        typedef result_type type;
    };

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
