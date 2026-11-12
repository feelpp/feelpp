/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2005-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FE_HPP
#define FEELPP_FE_HPP 1

#include <cstdint>
#include <limits>

#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/referencedoflayout.hpp>
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
namespace detail
{
/** @brief Select field rank while retaining the primal order policy. */
template<class P>
using FiniteElementPolynomialSet = mp::mp_if_c<P::is_scalar,
    PolynomialSet<P, Scalar, orderSpecification<P>()>,
    mp::mp_if_c<P::is_vectorial, PolynomialSet<P, Vectorial, orderSpecification<P>()>,
    mp::mp_if_c<P::is_tensor2 && is_symm_v<typename P::polyset_type>,
               PolynomialSet<P, Tensor2Symm, orderSpecification<P>()>,
               PolynomialSet<P, Tensor2, orderSpecification<P>()>>>>;
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
class FiniteElement : public ::Feel::detail::FiniteElementPolynomialSet<P>,
                      public ::Feel::detail::ReferenceDofLayout<FiniteElement<P, PDual, Pts>>
{
    /** Polynomial-set selection shared by all field ranks. */
    using super = ::Feel::detail::FiniteElementPolynomialSet<P>;

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
    static constexpr bool islinear_simplex = super::is_order_static && P::convex_type::is_simplex && (super::nOrder == 1);
    static constexpr bool islinear_hypercube = super::is_order_static && P::convex_type::is_hypercube && (super::nDim==1)  && (super::nOrder == 1);
    static constexpr bool islinear  = islinear_hypercube || islinear_simplex;
    static constexpr fem::transformation_type trans = islinear?fem::LINEAR:fem::NONLINEAR;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct the Ciarlet basis by inverting the owned dual-on-primal matrix. */
    FiniteElement( dual_space_type const& pdual )
        :
        super( pdual.primalSpace() ),
        M_dual( pdual )
    {
        DVLOG(2) << "============================================================\n";
        DVLOG(2) << "New FE \n";
        ublas::matrix<value_type> A( M_dual( M_dual.primalSpace() ) );
        if ( A.size1() != A.size2() )
            throw std::invalid_argument("invalid Ciarlet triple: dual-on-primal matrix must be square");
        //std::cout << "[FiniteElement] A = " << A << "\n";

        ublas::matrix<value_type> D = ublas::identity_matrix<value_type>( A.size1(), A.size2() );
        LU<ublas::matrix<value_type> > lu( A );
        if ( !lu.isNonsingular() )
            throw std::runtime_error("finite-element Vandermonde matrix is singular");
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

        this->setCoefficient( ublas::trans( C ) );

        //std::cout << "coeff = " << M_pset.coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        //std::cout << "d_x = " << M_pset.derivate(0).coeff() << "\n";
        DVLOG(2) << "============================================================\n";
    }
    FiniteElement( FiniteElement const& ) = default;
    FiniteElement( FiniteElement&& ) = default;

    ~FiniteElement() override = default;

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& ) = default;
    self_type& operator=( self_type&& ) = default;

    /** @brief Evaluate one basis function as a scalar or a vector of components. */
    template<typename AE>
    auto operator()( uint16_type i, ublas::vector_expression<AE> const& pt ) const
    {
        auto value = this->evaluate( i, pt );
        if constexpr ( super::is_scalar ) return value(0);
        else return value;
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
            return static_cast<uint16_type>( M_dual.primalSpace().order() );
        else if constexpr ( requires( super const& s ) { s.order(); } )
            return static_cast<uint16_type>( static_cast<super const&>( *this ).order() );
        else
            return static_cast<uint16_type>( nOrder );
    }

    /**
     * \return polynomial degree of the primal reference space.
     *
     * This is intentionally distinct from the public family order: for
     * shifted-order families such as RT, BDM, and Nedelec second kind the two
     * values need not coincide.
     */
    [[nodiscard]] uint16_type polynomialDegree() const
    {
        if constexpr ( requires( primal_space_type const& p ) { p.order(); } )
            return static_cast<uint16_type>( M_dual.primalSpace().order() );
        else
            return static_cast<uint16_type>( primal_space_type::nOrder );
    }

    //! return true if finite element is linear, false otherwise
    static constexpr bool isLinear() requires( super::is_order_static ) { return islinear; }

    /** @brief Test linearity using the effective expansion degree. */
    [[nodiscard]] bool isLinear() const requires( super::is_order_dynamic )
    {
        return polynomialDegree() == 1 && (P::convex_type::is_simplex || super::nDim == 1);
    }
    
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
        return M_dual.primalSpace();
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

    /** Shared entity descriptor retained for source compatibility. */
    using DofAttachment = typename ::Feel::detail::ReferenceDofLayout<self_type>::DofAttachment;
    /** Shared flattened local-layout descriptor. */
    using LocalDofLayout = typename ::Feel::detail::ReferenceDofLayout<self_type>::LocalDofLayout;

    //@}


private:

    dual_space_type M_dual;

};
}  // Feel
#endif /* __FiniteElement_H */
