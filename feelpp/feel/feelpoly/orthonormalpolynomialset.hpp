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
         template<uint16_type,int,uint16_type> class Convex = Simplex>
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
public PolynomialSet<Dubiner<Dim, RealDim, Order, Normalized<true>, T, StorageUBlas>, PolySetType >
{
    typedef PolynomialSet<Dubiner<Dim, RealDim, Order, Normalized<true>, T, StorageUBlas>, PolySetType > super;
public:

    static const uint16_type TAG = TheTAG;
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static const bool isTransformationEquivalent = true;
    typedef OrthonormalPolynomialSet<Dim, Order,RealDim, PolySetType, T, TheTAG, Simplex> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_tensor2symm = polyset_type::is_tensor2 && is_symm_v<polyset_type>;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static const bool is_product = true;
    static const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    typedef typename super::component_type component_type;

    typedef T value_type;
    typedef Dubiner<Dim, RealDim, Order, Normalized<true>, T, StorageUBlas> basis_type;
    typedef Simplex<Dim, Order, /*RealDim*/Dim> convex_type;
    template<int O>
    struct convex
    {
        typedef Simplex<Dim, O, /*RealDim*/Dim> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;
    typedef typename matrix_node<value_type>::type points_type;

    /**
     * local interpolant is undefined
     */
    typedef boost::none_t local_interpolant_type;

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

    OrthonormalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );

#if 0
        if ( !( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                       m ) < 1e-10 ) )
            std::cout << "m1=" << m << "\n"
                      << "m2=" << polyset_type::toMatrix( polyset_type::toType( m ) ) << "\n"
                      << ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) - m ) << "\n";

        FEELPP_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                              m ) < 1e-10 )( m ).warn ( "invalid transformation" );
#endif
        this->setCoefficient( polyset_type::toType( m ), true );

        if ( is_tensor2symm )
        {
            M_unsymm2symm.resize( nComponents*nLocalDof );
            for ( uint16_type l = 0; l < nLocalDof; ++l )
            {
                for ( int c1 = 0; c1 < nComponents1; ++c1)
                {
                    for ( int c2 = c1+1; c2 < nComponents2; ++c2 )
                    {
                        const int k = Feel::detail::symmetricIndex(c1,c2,nComponents1);
                        M_unsymm2symm[ nLocalDof*(nComponents1*c1+c2) + l ] = nLocalDof*k+l;
                        M_unsymm2symm[ nLocalDof*(nComponents1*c2+c1) + l ] = nLocalDof*k+l;
                    }
                    const int k = Feel::detail::symmetricIndex(c1,c1,nComponents1);
                    M_unsymm2symm[ nLocalDof*(nComponents1*c1+c1) + l ] = nLocalDof*k+l;
                }
            }
        }

    }

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
public PolynomialSet<Legendre<Dim, RealDim, Order, Normalized<true>, T>, PolySetType >
{
    typedef PolynomialSet<Legendre<Dim, RealDim, Order, Normalized<true>, T>, PolySetType > super;
public:

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RealDim;
    static const bool isTransformationEquivalent = true;

    typedef OrthonormalPolynomialSet<Dim, Order, RealDim, PolySetType, T, TheTAG, Hypercube> self_type;
    typedef self_type component_basis_type;

    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_tensor2symm = polyset_type::is_tensor2 && is_symm_v<polyset_type>;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const bool is_continuous = false;
    static const bool is_modal = true;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const bool is_product = true;
    static const bool isContinuous = false;
    typedef Discontinuous continuity_type;

    typedef typename super::component_type component_type;
    typedef T value_type;
    typedef Legendre<Dim, RealDim, Order, Normalized<true>, T> basis_type;
    typedef Hypercube<Dim, Order, /*RealDim*/Dim> convex_type;
    typedef typename matrix_node<value_type>::type points_type;

    /**
     * local interpolant is undefined
     */
    typedef boost::none_t local_interpolant_type;

    template<int O>
    struct convex
    {
        typedef Hypercube<Dim, O, nDim/*RealDim*/> type;
    };
    typedef Reference<convex_type, nDim, nOrder, nDim/*nRealDim*/, value_type> reference_convex_type;

    typedef typename super::polynomial_type polynomial_type;

    //!< Number of degrees of freedom per vertex
    static const uint16_type nDofPerVertex = convex_type::nbPtsPerVertex;
    //!< Number of degrees  of freedom per edge
    static const uint16_type nDofPerEdge = convex_type::nbPtsPerEdge;
    //!< Number of degrees  of freedom per face
    static const uint16_type nDofPerFace = convex_type::nbPtsPerFace;

    //!< Number of degrees  of freedom per volume
    static const uint16_type nDofPerVolume = convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = convex_type::numPoints;

    static const uint16_type nDof = nLocalDof;
    static const uint16_type nNodes = nDof;
    static const uint16_type nDofGrad = super::nDim*nDof;
    static const uint16_type nDofHess = super::nDim*super::nDim*nDof;

    OrthonormalPolynomialSet()
        :
        super( basis_type() )

    {
        ublas::matrix<value_type> m( ublas::identity_matrix<value_type>( nComponents*convex_type::polyDims( nOrder ) ) );
#if 0
        if ( is_tensor2 )
            std::cout << "[orthonormalpolynomialset] m = " << m << "\n";

        FEELPP_ASSERT( ublas::norm_frobenius( polyset_type::toMatrix( polyset_type::toType( m ) ) -
                                              m ) < 1e-10 )( m ).warn ( "invalid transformation" );
#endif
        this->setCoefficient( polyset_type::toType( m ), true );
    }

    OrthonormalPolynomialSet<Dim, Order, RealDim,Scalar,T, TheTAG, Hypercube > toScalar() const
    {
        return OrthonormalPolynomialSet<Dim, Order, RealDim, Scalar,T, TheTAG, Hypercube >();
    }
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

    points_type points() const
    {
        return points_type();
    }
    points_type points( int f ) const
    {
        return points_type();
    }
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
};

} // Feel
#endif /* __OrthonormalPolynomialSet_H */
