/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-06

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2011-2021 Feel++ Consortium

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
#ifndef FEELPP_FE_HPP
#define FEELPP_FE_HPP 1

#include <cstdint>
#include <limits>

#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/polynomialset.hpp>

namespace Feel
{
template<typename Poly, template<uint16_type> class PolySetType, int OrderSpec> class PolynomialSet;
namespace detail
{
template<uint16_type Dim,
         int Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG,
         template<int,int,int> class Convex>
class OrthonormalPolynomialSet;
}
/**
 * \class FiniteElement
 * \brief Finite element following Ciarlet framework
 *
 *  \ingroup Polynomial
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename P,
         template<class Pr,  template<class,int,class> class Pt> class PDual,
         template<class,int,class> class Pts>
class FiniteElement :
    public mpl::if_<mpl::bool_<P::is_scalar>,
                    mpl::identity<PolynomialSet<P, Scalar, P::nOrder> >,
                    typename mpl::if_<mpl::bool_<P::is_vectorial>,
                                      mpl::identity<PolynomialSet<P, Vectorial, P::nOrder> >,
                                      typename mpl::if_<mpl::bool_<P::is_tensor2 && is_symm_v<typename P::polyset_type>>,
                                                        mpl::identity<PolynomialSet<P, Tensor2Symm, P::nOrder>> ,
                                                        mpl::identity<PolynomialSet<P, Tensor2, P::nOrder>>
                                                        >::type
                                      >::type
                    >::type::type
{
    using super = typename mpl::if_<mpl::bool_<P::is_scalar>,
                                    mpl::identity<PolynomialSet<P, Scalar, P::nOrder> >,
                                    typename mpl::if_<mpl::bool_<P::is_vectorial>,
                                                      mpl::identity<PolynomialSet<P, Vectorial, P::nOrder> >,
                                                      typename mpl::if_<mpl::bool_<P::is_tensor2 && is_symm_v<typename P::polyset_type>>,
                                                                        mpl::identity<PolynomialSet<P, Tensor2Symm, P::nOrder>>,
                                                                        mpl::identity<PolynomialSet<P, Tensor2, P::nOrder>>
                                                                        >::type
                                                      >::type
                                    >::type::type;

public:

    /** @name Typedefs
     */
    //@{

    typedef FiniteElement<P, PDual, Pts> self_type;

    using value_type = typename P::value_type;

    typedef P primal_space_type;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename primal_space_type::polyset_type polyset_type;

    static inline const bool is_modal = false;

    typedef PDual<P, Pts> dual_space_type;

    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;
    typedef typename super::self_type polynomialset_type;

    typedef typename super::polynomial_type polynomial_type;
    typedef typename super::polynomial_view_type polynomial_view_type;

    //!< Total number of degrees of freedom (equal to refEle::nDof)
    static inline const uint16_type nLocalDof = dual_space_type::nLocalDof;
    //!< Number of degrees of freedom per vertex
    static inline const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    //!< Number of degrees  of freedom per edge
    static inline const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    //!< Number of degrees  of freedom per face
    static inline const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    //!< Number of degrees  of freedom per volume
    static inline const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    //! Compile-time order placeholder (for dynamic order, this keeps compatibility placeholder semantics)
    static constexpr int nOrder = super::nOrder;

    static constexpr uint16_type nDof = nLocalDof;
    static constexpr uint16_type nNodes = nDof;
    static constexpr uint16_type nDofGrad = super::nDim*nDof;
    static constexpr uint16_type nDofHess = super::nDim*super::nDim*nDof;
    static constexpr bool islinear_simplex = P::convex_type::is_simplex && (super::nOrder == 1);
    static constexpr bool islinear_hypercube = P::convex_type::is_hypercube && (super::nDim==1)  && (super::nOrder == 1);
    static constexpr bool islinear  = islinear_hypercube || islinear_simplex;
    static constexpr fem::transformation_type trans = islinear?fem::LINEAR:fem::NONLINEAR;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    FiniteElement( dual_space_type const& pdual )
        :
        super( pdual.primalSpace() ),
        M_dual( pdual ),
        M_primal( M_dual.primalSpace() )
    {
        DVLOG(2) << "============================================================\n";
        DVLOG(2) << "New FE \n";
        ublas::matrix<value_type> A( M_dual( M_primal ) );
        //std::cout << "[FiniteElement] A = " << A << "\n";

        ublas::matrix<value_type> D = ublas::identity_matrix<value_type>( A.size1(), A.size2() );
        LU<ublas::matrix<value_type> > lu( A );
        ublas::matrix<value_type> C = lu.solve( D );
        //std::cout << "[FiniteElement] D = " << D << "\n";
        //std::cout << "[FiniteElement] C = " << C << "\n";
        DVLOG(2) << "is singular : " << lu.isNonsingular() << "\n"
                      << "det(A) =  " << lu.det() << "\n";
#if 0

        if ( !lu.isNonsingular() )
        {
            std::cout << "A=" << A << "\n"
                      << "D=" << D << "\n"
                      << "C=" << C << "\n";
        }

#endif
        FEELPP_ASSERT( lu.isNonsingular() )( A )( D )( C ).error( "vandermonde matrix is singular" );

        this->setCoefficient( ublas::trans( C ) );

        //M_pset = polynomialset_type( M_primal, C );

        //std::cout << "coeff = " << M_pset.coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        DVLOG(2) << "============================================================\n";
    }
    FiniteElement( FiniteElement const & fe )
        :
        super( fe ),
        M_dual( fe.M_dual ),
        M_primal( fe.M_primal )
    {}

    ~FiniteElement() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& fe )
    {
        if ( this != &fe )
        {
            super::operator=( fe );
            M_primal = fe.M_primal;
            M_dual = fe.M_dual;
        }

        return *this;
    }

    template<typename AE>
    value_type operator()( uint16_type i, ublas::vector_expression<AE> const& pt ) const
    {
        return this->evaluate( i, pt );
    }

    template<typename AE>
    value_type operator()( ublas::vector_expression<AE> const& pt ) const
    {
        matrix_type m( pt().size(), 1 );
        ublas::column( m, 0 ) = pt;
        ublas::vector<value_type> r( ublas::column( this->evaluate( m ), 0 ) );
        return ublas::inner_prod( r, ublas::scalar_vector<value_type>( r.size(), 1.0 ) );
    }


    matrix_type operator()( points_type const& pts ) const
    {
        return this->evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    //!
    //! @return semantic order of the finite element
    //!
    [[nodiscard]] uint16_type order() const
    {
        if constexpr ( requires( primal_space_type const& p ) { p.order(); } )
            return static_cast<uint16_type>( M_primal.order() );
        else if constexpr ( requires( super const& s ) { s.order(); } )
            return static_cast<uint16_type>( static_cast<super const&>( *this ).order() );
        else
            return static_cast<uint16_type>( nOrder );
    }

    /**
     * @deprecated Use order() instead
     */
    [[nodiscard]] uint16_type runtimeOrder() const
    {
        return this->order();
    }

    //! return true if finite element is linear, false otherwise
    static constexpr bool isLinear() { return islinear; }
    
    /**
     * \return the domain shape of the finite element
     */
    void domainShape() const {}

    /**
     * \return the number of points associated with FE
     */
    uint16_type nbPoints() const
    {
        return points().size2();
    }

    /**
     * \return the polynomial set defining the finite element
     */
    // Que devient cette fonction ??
    // polynomialset_type const& functionShape() const { return M_pset; }

    /**
     * \return the dual basis of the finite element
     */
    primal_space_type const& primal() const
    {
        return M_primal;
    }

    /**
     * \return the dual basis of the finite element
     */
    dual_space_type const& dual() const
    {
        return M_dual;
    }

    /**
     * \return points associated with the lagrange finite element
     */
    points_type const& points() const
    {
        return M_dual.points();
    }

    /**
     * get the points associated with the finite element on a face \c
     * f if any
     *
     * \arg f face index
     *
     * \return points associated with a face of the lagrange finite
     * element
     */
    points_type const& points( uint16_type f ) const
    {
        return M_dual.points( f );
    }

    points_type edgePoints( uint16_type e ) const
        {
            return M_dual.edgePoints( e );
        }

    points_type vertexPoints( uint16_type v ) const
        {
            return M_dual.vertexPoints( v );
        }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const override = 0;

    struct DofAttachment
    {
        static constexpr uint16_type invalid_id = std::numeric_limits<uint16_type>::max();
        int8_t entityDim = -1;
        uint16_type entityId = invalid_id;
        uint16_type ordinal = invalid_id;
        uint16_type kind = 0;

        [[nodiscard]] bool isValid() const noexcept { return entityDim >= 0; }
    };

    struct LocalDofLayout
    {
        uint16_type localDofId = 0;
        uint16_type parentLocalDofId = 0;
        uint16_type component = 0;
        DofAttachment attachment;
    };

    /**
     * \return number of local dofs per component.
     *
     * This is FE-owned so consumers do not assume local dof layout.
     */
    virtual uint16_type localDofPerComponent() const
    {
        return nLocalDof;
    }

    /**
     * \return number of local dofs owned by the FE layout.
     *
     * When \p perComponent is false this is the flattened local cardinality
     * exposed to local-to-global tables. When it is true this is the parent
     * local cardinality for one component.
     */
    virtual uint16_type localDofCount( bool perComponent = false ) const
    {
        const uint16_type localDofPerComp = this->localDofPerComponent();
        if ( perComponent )
            return localDofPerComp;

        if constexpr ( P::is_product )
            return static_cast<uint16_type>( super::nComponents * localDofPerComp );
        else
            return localDofPerComp;
    }

    /**
     * \return local dof id from (parent local dof id, component)
     *
     * Layout ownership stays in FE implementation.
     */
    virtual uint16_type localDofId( uint16_type parentLocalDofId, uint16_type component = 0 ) const
    {
        if constexpr ( P::is_product )
        {
            const uint16_type localDofPerComp = this->localDofPerComponent();
            FEELPP_ASSERT( parentLocalDofId < localDofPerComp )
                ( component )( parentLocalDofId )( localDofPerComp ).error( "invalid parent local dof index" );
            const auto localDof = static_cast<uint16_type>( component * localDofPerComp + parentLocalDofId );
            FEELPP_ASSERT( this->dofParent( localDof ) == parentLocalDofId )
                ( component )( parentLocalDofId )( localDof ).error( "invalid component/local dof layout" );
            return localDof;
        }
        else
        {
            FEELPP_ASSERT( component == 0 )
                ( component )( parentLocalDofId ).error( "non-product FE supports component 0 only" );
            return parentLocalDofId;
        }
    }

    /**
     * \return entity attachment metadata for a local dof.
     *
     * Default implementation returns unknown entity attachment and keeps
     * dof type as functional kind.
     */
    virtual DofAttachment dofAttachment( uint16_type localDofId ) const
    {
        return DofAttachment{ .entityDim = -1,
                              .entityId = DofAttachment::invalid_id,
                              .ordinal = DofAttachment::invalid_id,
                              .kind = this->dofType( localDofId ) };
    }

    /**
     * \return full FE-owned layout descriptor for a local dof.
     */
    LocalDofLayout localDofLayout( uint16_type localDofId ) const
    {
        return LocalDofLayout{ .localDofId = localDofId,
                               .parentLocalDofId = this->dofParent( localDofId ),
                               .component = this->component( localDofId ),
                               .attachment = this->dofAttachment( localDofId ) };
    }

    /**
     * \return number of local dofs attached to a reference entity.
     *
     * The default implementation derives the count from FE-owned local dof
     * layouts. FE families may override when a cheaper/runtime-aware path is
     * available.
     */
    virtual uint16_type localDofCountOnEntity( uint16_type topologicalDim,
                                               uint16_type localEntity,
                                               bool perComponent = false ) const
    {
        uint16_type count = 0;
        const uint16_type nLocalDof = this->localDofCount();
        for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
        {
            auto const layout = this->localDofLayout( localDof );
            auto const& attachment = layout.attachment;
            if ( !attachment.isValid() )
                continue;
            if ( perComponent && layout.component != 0 )
                continue;
            if ( attachment.entityDim == static_cast<int8_t>( topologicalDim ) &&
                 attachment.entityId == localEntity )
                ++count;
        }
        return count;
    }

    /**
     * \return number of local dofs attached to a reference codim-1 facet.
     */
    virtual uint16_type localDofCountOnFacet( uint16_type localFacet = 0,
                                              bool perComponent = false ) const
    {
        if constexpr ( super::nDim == 1 )
            return this->localDofCountOnEntity( 0, localFacet, perComponent );
        else
        {
            using convex_type = typename super::convex_type;
            using face_type = typename convex_type::topological_face_type;

            auto isInFacetClosure = [localFacet]( DofAttachment const& attachment ) -> bool
            {
                if ( !attachment.isValid() )
                    return false;

                if ( attachment.entityDim == 0 )
                {
                    for ( uint16_type localVertex = 0; localVertex < face_type::numVertices; ++localVertex )
                        if ( static_cast<uint16_type>( convex_type::f2p( localFacet, localVertex ) ) == attachment.entityId )
                            return true;
                    return false;
                }

                if ( attachment.entityDim == 1 )
                {
                    if constexpr ( super::nDim == 2 )
                        return attachment.entityId == localFacet;
                    else
                    {
                        for ( uint16_type localEdge = 0; localEdge < face_type::numEdges; ++localEdge )
                            if ( static_cast<uint16_type>( convex_type::f2e( localFacet, localEdge ) ) == attachment.entityId )
                                return true;
                        return false;
                    }
                }

                if ( attachment.entityDim == 2 )
                {
                    if constexpr ( super::nDim == 3 )
                        return attachment.entityId == localFacet;
                    else
                        return false;
                }

                return false;
            };

            uint16_type count = 0;
            const uint16_type nLocalDof = this->localDofCount();
            for ( uint16_type localDof = 0; localDof < nLocalDof; ++localDof )
            {
                auto const layout = this->localDofLayout( localDof );
                if ( perComponent && layout.component != 0 )
                    continue;
                if ( isInFacetClosure( layout.attachment ) )
                    ++count;
            }
            return count;
        }
    }

    /**
     * \return true when a local dof has a representative interpolation point.
     *
     * Point-evaluation finite elements keep the historical behavior: the
     * representative point is the FE point indexed by dofParent(localDofId).
     * Moment/modal finite elements should override this method when no such
     * geometric representative exists.
     */
    virtual bool dofHasRepresentativePoint( uint16_type localDofId ) const
    {
        (void)localDofId;
        return true;
    }

    /**
     * \return point index in points() used to represent localDofId.
     *
     * This method is only meaningful when dofHasRepresentativePoint(localDofId)
     * returns true.
     */
    virtual uint16_type dofRepresentativePointIndex( uint16_type localDofId ) const
    {
        return this->dofParent( localDofId );
    }

    /**
     * \return finite-element owned functional kind for a local dof.
     *
     * The default preserves legacy DofAttachment::kind values. FE families with
     * moment functionals can override this without changing their public
     * DofAttachment compatibility value.
     */
    virtual uint16_type dofFunctionalKind( uint16_type localDofId ) const
    {
        return this->dofType( localDofId );
    }

    //! \return the component of a local dof
    virtual uint16_type component( uint16_type localDofId ) const = 0;

    //! \return a parent local dof id for each component (for example, the first component)
    virtual uint16_type dofParent( uint16_type localDofId ) const = 0;

    //! \return the type of a local dof
    virtual uint16_type dofType( uint16_type localDofId ) const = 0;

    //! give an unsymmetric dof index i, provide the symmetric one
    virtual uint16_type unsymmToSymm( uint16_type i ) const { return i; }

    //@}


private:

    dual_space_type M_dual;
    primal_space_type const& M_primal;

};
}  // Feel
#endif /* __FiniteElement_H */
