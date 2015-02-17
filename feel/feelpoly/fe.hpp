/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-06

  Copyright (C) 2005,2006 EPFL

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
   \file fe.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __FiniteElement_H
#define __FiniteElement_H 1

#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/polynomialset.hpp>

namespace Feel
{
template<typename Poly, template<uint16_type> class PolySetType > class PolynomialSet;
namespace detail
{
template<uint16_type Dim,
         uint16_type Order,
         uint16_type RealDim,
         template<uint16_type> class PolySetType,
         typename T,
         uint16_type TheTAG,
         template<uint16_type,uint16_type,uint16_type> class Convex>
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
         template<class Pr,  template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
class FiniteElement :
    public mpl::if_<mpl::bool_<P::is_scalar>,
    mpl::identity<PolynomialSet<P, Scalar> >,
    mpl::identity<PolynomialSet<P, Vectorial> > >::type::type
{
    typedef typename mpl::if_<mpl::bool_<P::is_scalar>,
            mpl::identity<PolynomialSet<P, Scalar> >,
            mpl::identity<PolynomialSet<P, Vectorial> > >::type::type super;

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

    static const bool is_modal = false;

    typedef PDual<P, Pts> dual_space_type;

    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;
    typedef typename super::self_type polynomialset_type;

    typedef typename super::polynomial_type polynomial_type;
    typedef typename super::polynomial_view_type polynomial_view_type;

    //!< Total number of degrees of freedom (equal to refEle::nDof)
    static const uint16_type nLocalDof = dual_space_type::nLocalDof;
    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;
    static const bool islinear_simplex = P::convex_type::is_simplex && (super::nOrder == 1);
    static const bool islinear_hypercube = P::convex_type::is_hypercube && (super::nDim==1)  && (super::nOrder == 1);
    static const bool islinear  = islinear_hypercube || islinear_simplex;
    static const fem::transformation_type trans = ( fem::transformation_type )mpl::if_<mpl::bool_<islinear>,
                                                                                       mpl::int_<fem::LINEAR>,
                                                                                       mpl::int_<fem::NONLINEAR> >::type::value;
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

    ~FiniteElement()
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

    /**
     * \return the family name of the finite element
     */
    virtual std::string familyName() const = 0;

    //@}


private:

    dual_space_type M_dual;
    primal_space_type const& M_primal;

};

template<typename P,
         template<class Pr,  template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
const uint16_type FiniteElement<P,PDual,Pts>::nLocalDof;

template<typename P,
         template<class Pr,  template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
const uint16_type FiniteElement<P,PDual,Pts>::nDofPerVertex;

template<typename P,
         template<class Pr, template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
const uint16_type FiniteElement<P,PDual,Pts>::nDofPerEdge;

template<typename P,
         template<class Pr,  template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
const uint16_type FiniteElement<P,PDual,Pts>::nDofPerFace;

template<typename P,
         template<class Pr,  template<class,uint16_type,class> class Pt> class PDual,
         template<class,uint16_type,class> class Pts>
const uint16_type FiniteElement<P,PDual,Pts>::nDofPerVolume;
}
#endif /* __FiniteElement_H */
