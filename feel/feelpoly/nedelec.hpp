/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-01-14

  Copyright (C) 2011 Université de Grenoble 1 (Joseph Fourier)

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
   \file nedelec.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-11-25
 */
#ifndef FEELPP_NEDELEC_HPP
#define FEELPP_NEDELEC_HPP 1

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

#include <boost/type_traits.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>


#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/operations.hpp>
#include <feel/feelpoly/functionals.hpp>
#include <feel/feelpoly/functionals2.hpp>
#include <feel/feelpoly/quadpoint.hpp>
#include <feel/feelpoly/fe.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelpoly/hcurlpolynomialset.hpp>

namespace Feel
{
/**
 * - NED1 H(curl) Nedelec finite element of the first kind
 * - NED2 H(curl) Nedelec finite element of the second kind (full polynomial space)
 */
enum class NedelecKind { NED1, NED2 };

namespace detail
{
template<typename P>
struct times_rotx
{
    typedef typename P::value_type value_type;
    typedef typename P::points_type points_type;
    times_rotx ( P const& p, int c  )
        :
        M_p ( p ),
        M_c( c )
    {
        //VLOG(1) << "component : " << c << std::endl;
    }
    typename ublas::vector<value_type> operator() ( points_type const& __pts ) const
    {
#if 0
        VLOG(1) << "times_rotx(pts) : " << __pts << std::endl;
        VLOG(1) << "times_rotx(pts) : " << M_p.evaluate( __pts ) << std::endl;
        VLOG(1) << "times_rotx(coeff) : " << M_p.coefficients() << std::endl;
#endif

        // __pts[c] * p( __pts )
        if ( M_c == 0 )
            return ublas::element_prod( ublas::row( __pts, 1 ),
                                        ublas::row( M_p.evaluate( __pts ), 0 ) );

        // c == 1
        return ublas::element_prod( -ublas::row( __pts, 0 ),
                                    ublas::row( M_p.evaluate( __pts ), 0 ) );
    }
    //P const& M_p;
    P M_p;
    int M_c;
};

#if 0
template< class T >
struct extract_all_poly_indices
{
    T start;
    extract_all_poly_indices( T comp, T N )
        :
        start( comp*N )
    { }

    T operator()()
    {
        return start++;
    }
};
#endif

}// detail

template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED2,
         typename T = double,
         uint16_type TheTAG = 0 >
class NedelecPolynomialSet{};

template<uint16_type O,
         typename T,
         uint16_type TheTAG>
class NedelecPolynomialSet<2,O,NedelecKind::NED1,T,TheTAG>
    :
    public Feel::detail::OrthonormalPolynomialSet<2, O+1, 2, Vectorial, T, TheTAG, Simplex>
{
    static const int N = 2;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED1,T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nComponents = super::nComponents;
    static const bool is_product = false;
    NedelecPolynomialSet()
        :
        super()
    {
        VLOG(1) << "[Nedelec1stKindset] nOrder = " << nOrder << "\n";
        VLOG(1) << "[Nedelec1stKindset] O = " << O << "\n";
        uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
        uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
        uint16_type dim_Pkm1 = ( nOrder==1 )?0:convex_type::polyDims( nOrder-2 );
#if 1
        VLOG(1) << "[Nedelec1stKindset] dim_Pkp1 = " << dim_Pkp1 << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pk   = " << dim_Pk << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pkm1 = " << dim_Pkm1 << "\n";
#endif
        // (P_k)^d
        Pkp1_v_type Pkp1_v;
        vectorial_polynomialset_type Pk_v( Pkp1_v.polynomialsUpToDimension( dim_Pk ) );
#if 1
        LOG(INFO) << "[Nedelec1stKindset] Pk_v =" << Pk_v.coeff() << "\n";
        VLOG(1) << "[Nedelec1stKindset] Pk_v =" << Pk_v.coeff() << "\n";
#endif
        // P_k
        Pkp1_s_type Pkp1;
        scalar_polynomialset_type Pk ( Pkp1.polynomialsUpToDimension( dim_Pk ) );
#if 1
        LOG(INFO) << "[Nedelec1stKindset] Pk =" << Pk.coeff() << "\n";
        VLOG(1) << "[Nedelec1stKindset] Pk(0) =" << Pk.polynomial( 0 ).coefficients() << "\n";
#endif

        // curl(x) P_k \ P_{k-1}
        IMGeneral<convex_type::nDim, 2*nOrder,value_type> im;
        //VLOG(1) << "[Nedelec1stKindPset] im.points() = " << im.points() << std::endl;
        ublas::matrix<value_type> xPkc( nComponents*( dim_Pk-dim_Pkm1 ),Pk.coeff().size2() );

        //VLOG(1) << "[Nedelec1stKindPset] before xPkc = " << xPkc << "\n";
        for ( int l = dim_Pkm1, i = 0; l < dim_Pk; ++l, ++i )
        {
            for ( int j = 0; j < convex_type::nDim; ++j )
            {
                Feel::detail::times_rotx<scalar_polynomial_type> xp( Pk.polynomial( l ), j );
                ublas::row( xPkc,i*nComponents+j )=
                    ublas::row( Feel::project( Pkp1,
                                               xp,
                                               im ).coeff(), 0 );
            }
        }

        //VLOG(1) << "[Nedelec1stKindPset] after xPkc = " << xPkc << "\n";
        vectorial_polynomialset_type xPk( typename super::basis_type(), xPkc, true );
        LOG(INFO) << "[Nedelec1stKindset] xPk =" << xPk.coeff() << "\n";
        //VLOG(1) << "[Nedelec1stKindPset] here 1\n";
        // (P_k)^d + x P_k
        LOG(INFO) << "[Nedelec1stKindPset] Nedelec1stKind Poly coeff = " << unite( Pk_v, xPk ).coeff() << "\n";
        this->setCoefficient( unite( Pk_v, xPk ).coeff(), true );
        //VLOG(1) << "[Nedelec1stKindPset] here 2\n";
    }


};

template<uint16_type O,
         typename T,
         uint16_type TheTAG >
class NedelecPolynomialSet<3,O,NedelecKind::NED1,T,TheTAG>
    :
        public Feel::detail::OrthonormalPolynomialSet<3, O+1, 3, Vectorial, T, TheTAG, Simplex>
{
    static const int N = 3;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED1, T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nComponents = super::nComponents;
    static const bool is_product = false;
    NedelecPolynomialSet()
        :
        super()
    {
        VLOG(1) << "[Nedelec1stKindset] nDim = " << nDim << "\n";
        VLOG(1) << "[Nedelec1stKindset] nOrder = " << nOrder << "\n";
        VLOG(1) << "[Nedelec1stKindset] O = " << O << "\n";
        uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
        uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
        uint16_type dim_Pkm1 = ( nOrder==1 )?0:convex_type::polyDims( nOrder-2 );
#if 1
        VLOG(1) << "[Nedelec1stKindset] dim_Pkp1 = " << dim_Pkp1 << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pk   = " << dim_Pk << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pkm1 = " << dim_Pkm1 << "\n";
#endif
        // (P_k)^d
        Pkp1_v_type Pkp1_v;
        VLOG(1) << "[Nedelec1stKindset] Pkp1_v =" << Pkp1_v.coeff() << "\n";
        vectorial_polynomialset_type Pk_v( Pkp1_v.polynomialsUpToDimension( dim_Pk ) );
        vectorial_polynomialset_type Pke_v( Pkp1_v.polynomialsUpToDimension( dim_Pk ) );
#if 1
        VLOG(1) << "[Nedelec1stKindset] Pk_v =" << Pk_v.coeff() << "\n";
#endif
        // P_k
        Pkp1_s_type Pkp1;
        VLOG(1) << "[Nedelec1stKindset] Pkp1 =" << Pkp1.coeff() << "\n";
        scalar_polynomialset_type Pk ( Pkp1.polynomialsUpToDimension( dim_Pk ) );
#if 1
        VLOG(1) << "[Nedelec1stKindset] Pk =" << Pk.coeff() << "\n";
        VLOG(1) << "[Nedelec1stKindset] Pk(0) =" << Pk.polynomial( 0 ).coefficients() << "\n";
#endif

        // curl(x) P_k \ P_{k-1}
        IMGeneral<convex_type::nDim, 2*nOrder,value_type> im;
        VLOG(1) << "[Nedelec1stKindPset] im.points() = " << im.points() << std::endl;
        //ublas::matrix<value_type> xPkcV( nComponents*( dim_Pk-dim_Pkm1 ),Pk.coeff().size2() );
        ublas::matrix<value_type> xPkcV( nComponents*nComponents*dim_Pk,Pk.coeff().size2() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> xPkc( xPkcV.data().begin(),
                                                                                                  xPkcV.size1(), xPkcV.size2() );
        xPkc.setZero();

        // evaluate P_{k+1} at quad points
        auto e_Pkp1 = Pkp1.evaluate( im.points() );
        VLOG(1) << "e_Pkp1 = " << e_Pkp1;
        // evaluate Pke at quad points
        auto e_Pkv = Pk_v.evaluate( im.points() );
        VLOG(1) << "e_Pkv = " << e_Pkv;
        auto const& X = im.points();
        auto const& W = im.weights();
        VLOG(1) << "[Nedelec1stKindPset] before xPkc = " << xPkc << "\n";
        for ( int c = 0; c < convex_type::nDim; ++c )
        {

            for ( int l = c*convex_type::nDim+dim_Pkm1; l < c*convex_type::nDim+dim_Pk; ++l )
                //for ( int l = 0, i = 0; l < dim_Pk; ++l, ++i )
            {
                VLOG(1) << "polynomial(c="<<c << ",l=" << l << ")= " << l;
                for ( int j = 0; j < convex_type::nDim; ++j )
                {
                    for( int q  = 0; q < im.nPoints(); ++q )
                    {
                        // compute cross(e_Pk(l,q),x)
                        auto b = W(q)*( X((j+2)%3,q)*e_Pkv(l+(j+1)%3,q)-
                                        X((j+1)%3,q)*e_Pkv(l+(j+2)%3,q) );
                        for( int p=0;p < xPkc.cols(); ++p )
                            xPkc(l+j,p) += b * e_Pkp1(p,q);
                    }
                }
            }
        }
        VLOG(1) << "[Nedelec1stKindPset] after xPkc = " << xPkc << "\n";
        VLOG(1) << "[Nedelec1stKindPset] after xPkcV = " << xPkcV << "\n";
        vectorial_polynomialset_type xPk( typename super::basis_type(), xPkcV, true );
        //VLOG(1) << "[Nedelec1stKindPset] here 1\n";
        // (P_k)^d + x P_k
        VLOG(1) << "[Nedelec1stKindPset] Nedelec1stKind Poly coeff = " << unite( Pk_v, xPk ).coeff() << "\n";
        this->setCoefficient( unite( Pk_v, xPk ).coeff(), true );
        //VLOG(1) << "[Nedelec1stKindPset] here 2\n";
    }


};

template<uint16_type N,
         uint16_type O,
         typename T,
         uint16_type TheTAG >
class NedelecPolynomialSet<N,O,NedelecKind::NED2,T,TheTAG>
    :
        public Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex>
{
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED1,T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nComponents = super::nComponents;
    static const bool is_product = false;
    NedelecPolynomialSet()
        :
        super()
    {
        VLOG(1) << "[Nedelec1stKindset] nOrder = " << nOrder << "\n";
        VLOG(1) << "[Nedelec1stKindset] O = " << O << "\n";
        uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
        uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
        uint16_type dim_Pkm1 = ( nOrder==1 )?0:convex_type::polyDims( nOrder-2 );
#if 1
        VLOG(1) << "[Nedelec1stKindset] dim_Pkp1 = " << dim_Pkp1 << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pk   = " << dim_Pk << "\n";
        VLOG(1) << "[Nedelec1stKindset] dim_Pkm1 = " << dim_Pkm1 << "\n";
#endif
        // (P_k)^d
        Pkp1_v_type Pkp1_v;
        //vectorial_polynomialset_type Pk_v( Pkp1_v );

        //VLOG(1) << "[Nedelec1stKindPset] after xPkc = " << xPkc << "\n";
        vectorial_polynomialset_type Pkv( typename super::basis_type(), Pkp1_v.coeff(), true );
        this->setCoefficient(  Pkv.coeff(), true );
        //VLOG(1) << "[Nedelec1stKindPset] here 2\n";
    }
};

namespace fem
{

namespace detail
{

template<typename Basis,
         template<class, uint16_type, class> class PointSetType>
class NedelecDualFirstKind
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder= super::nOrder;
    static const NedelecKind kind = NedelecKind::NED1;

    typedef typename super::primal_space_type primal_space_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    typedef typename primal_space_type::template convex<nOrder+1>::type convex_type;
    typedef Reference<convex_type, nDim, nOrder+1, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename convex_type::topological_face_type face_type;

    typedef typename primal_space_type::Pkp1_v_type Pkp1_v_type;
    typedef typename primal_space_type::vectorial_polynomialset_type vectorial_polynomialset_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    static const uint16_type nbPtsPerVertex = 0;
    static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace2D = reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerFace3D = 0;
    static const uint16_type nbPtsPerFace = ((nDim==2)?nbPtsPerFace2D:nbPtsPerFace3D)*nDim;
    static const uint16_type nbPtsPerVolume = 0;
    static const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+reference_convex_type::numEdges*nbPtsPerEdge );

    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = 0;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume = nbPtsPerVolume;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static const uint16_type nLocalDof = numPoints;

    static const uint16_type nFacesInConvex = mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<1> >,
                                                        mpl::int_<reference_convex_type::numVertices>,
                                                        typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                                                          mpl::int_<reference_convex_type::numEdges>,
                                                                          mpl::int_<reference_convex_type::numGeometricFaces> >::type >::type::value;

    NedelecDualFirstKind( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_pts_per_face( nFacesInConvex ),
        M_fset( primal )
    {
        VLOG(1) << "Nedelec finite element(dual): \n";
        VLOG(1) << " o- dim   = " << nDim << "\n";
        VLOG(1) << " o- order = " << nOrder << "\n";
        VLOG(1) << " o- kind = " << static_cast<int>(kind) << "\n";
        VLOG(1) << " o- numPoints      = " << numPoints << "\n";
        VLOG(1) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        VLOG(1) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        VLOG(1) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        VLOG(1) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        VLOG(1) << " o- nLocalDof      = " << nLocalDof << "\n";

        if ( nDim == 2 )
            CHECK( nLocalDof == nOrder*(nOrder+2) ) 
                << "Invalid finite element dimension. Expected : " 
                << nOrder*(nOrder+2) << " Actual: " << nLocalDof;

        LOG(INFO) << "Nedelec finite element(dual): \n";
        LOG(INFO) << " o- dim   = " << nDim << "\n";
        LOG(INFO) << " o- order = " << nOrder << "\n";
        LOG(INFO) << " o- kind = " << static_cast<int>(kind) << "\n";
        LOG(INFO) << " o- numPoints      = " << numPoints << "\n";
        LOG(INFO) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        LOG(INFO) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        LOG(INFO) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        LOG(INFO) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        LOG(INFO) << " o- nLocalDof      = " << nLocalDof << "\n";

        size_type nbDofPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                          mpl::int_<face_type::numEdges*nbPtsPerEdge>,
                                          mpl::int_<face_type::numTopologicalFaces*nbPtsPerEdge + nbPtsPerFace>
                                          >::type::value;

        size_type nbEdgesPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                             mpl::int_<face_type::numEdges>,
                                             mpl::int_<face_type::numTopologicalFaces>
                                             >::type::value;

        for(auto& m : M_pts_per_face)
            m.resize(nDim,nbDofPerFace);

        //VLOG(1) << "[Nedelec1stKind Dual] done 1\n";
        // compute  \f$ \ell_e( U ) = (U * n[e]) (edge_pts(e)) \f$
        typedef Functional<primal_space_type> functional_type;
        std::vector<functional_type> fset;

        // loop on each entity forming the convex of topological
        // dimension nDim-1 ( the faces)
        int d = (nDim == 2 )?nDim-1:1;
        for ( int p = 0, e = M_convex_ref.entityRange( d ).begin();
                e < M_convex_ref.entityRange( d ).end();
                ++e )
        {
            // Dof coord on current edge
            points_type Gt ( M_convex_ref.makePoints( d, e ) );

            for( int f = M_convex_ref.entityRange( nDim-1 ).begin(); f < M_convex_ref.entityRange( nDim-1 ).end(); ++f)
                {
                    // M_pts_per_face[f] : set of dofs on a face (face_type::numTopologicalFaces*nbPtsPerFace)
                    //for(int k=0; k<face_type::numTopologicalFaces; ++k)
                    for(int k=0; k<nbEdgesPerFace; ++k)
                        {
                            int curEdge;
                            if( nDim <= 2 )
                                curEdge = convex_type::f2e(f,f);
                            else
                                curEdge = convex_type::f2e(f,k);

                            if( curEdge == e )
                                {
                                    ublas::subrange( M_pts_per_face[f], 0, nDim, nbPtsPerEdge*k, nbPtsPerEdge*(k+1) ) = Gt;
                                }
                        }
                }

            // M_pts : set of dofs on an element
            if ( Gt.size2() )
            {
                VLOG(1) << "Gt = " << Gt << "\n";
                //VLOG(1) << "p = " << p << "\n";
                ublas::subrange( M_pts, 0, nDim, p, p+Gt.size2() ) = Gt;
                //for ( size_type j = 0; j < Gt.size2(); ++j )
                //M_eid[d].push_back( p+j );
                p+=Gt.size2();
            }

            typedef Feel::functional::DirectionalComponentPointsEvaluation<primal_space_type> dcpe_type;
            VLOG(1) << "tangent " << e << ":" << M_convex_ref.tangent( e ) << "\n";
            node_type dir= M_convex_ref.tangent(e);

            //dcpe_type __dcpe( primal, dir, M_pts_per_face[e] );
            dcpe_type __dcpe( primal, dir, Gt );
            std::copy( __dcpe.begin(), __dcpe.end(), std::back_inserter( fset ) );

        }

        // add integrals of tangential component
        if ( nOrder-1 > 0 )
        {
        }
        //VLOG(1) << "[Nedelec1stKind Dual] done 2" << std::endl;
        // add integral of moments
        if ( nOrder-2 >= 0 )
        {
            // we need more equations : add interior moment
            // indeed the space is orthogonal to Pk-1
            uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
            uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
            uint16_type dim_Pm1 = convex_type::polyDims( nOrder-2 );

            Pkp1_v_type Pkp1;

            vectorial_polynomialset_type Pkm1 ( Pkp1.polynomialsUpToDimension( dim_Pm1 ) );

            VLOG(1) << "[nedelec 1st kind] Dim Pkm1 = " << dim_Pm1 << "\n";
            VLOG(1) << "[nedelec 1st kind] Pkm1 = " << Pkm1.coeff() << "\n";
            VLOG(1) << "[nedelec 1st kind] Primal = " << primal.coeff() << "\n";
            for ( int i = 0; i < Pkm1.polynomialDimension(); ++i )
            {
                typedef functional::IntegralMoment<primal_space_type, vectorial_polynomialset_type> fim_type;
                //typedef functional::IntegralMoment<Pkp1_v_type, vectorial_polynomialset_type> fim_type;
                VLOG(1) << "P(" << i << ")=" << Pkm1.polynomial( i ).coeff() << "\n";
                fset.push_back( fim_type( primal, Pkm1.polynomial( i ) ) );
            }
        }

        VLOG(1) << "[Nedelec1stKind Dual] done 3, n fset = " << fset.size() << std::endl;
        M_fset.setFunctionalSet( fset );
        VLOG(1) << "[Nedelec1stKind DUAL matrix] mat = " << M_fset.rep() << "\n";
        VLOG(1) << "[Nedelec1stKind Dual] done 4\n";

    }

    /**
     * \return the point set
     */
    //pointset_type const& pointSet() const { return M_pset; }

    points_type const& points() const
    {
        return M_pts;
    }


    matrix_type operator()( primal_space_type const& pset ) const
    {
        //VLOG(1) << "Nedelec1stKind matrix = " << M_fset( pset ) << std::endl;
        return M_fset( pset );
    }

    points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_pts_per_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_pts_per_face[f], __i );
    }

private:
    /**
     * set the pointset at face \c f using points \c n
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_pts_per_face[f].resize( n.size1(), n.size2(), false );
        M_pts_per_face[f] = n;
    }

private:
    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    FunctionalSet<primal_space_type> M_fset;


};


template<typename Basis,
         template<class, uint16_type, class> class PointSetType>
class NedelecDualSecondKind
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder= super::nOrder;
    static const NedelecKind kind = NedelecKind::NED2;

    typedef typename super::primal_space_type primal_space_type;
    //typedef typename super::primal_space_type::template ChangeOrder<nOrder+1>::type primal_space_type;

    typedef typename primal_space_type::template ChangeOrder<nOrder+1>::type Pkp1_v_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    // typedef typename primal_space_type::template convex<nDim+nOrder-1>::type convex_type;
    // typedef Reference<convex_type, nDim, nDim+nOrder-1, nDim, value_type> reference_convex_type;

    typedef typename primal_space_type::template convex<nOrder+2>::type convex_type;
    typedef Reference<convex_type, nDim, nOrder+2, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename convex_type::topological_face_type face_type;

    uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
    uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
    uint16_type dim_Pkm1 = ( nOrder==1 )?0:convex_type::polyDims( nOrder-2 );

    typedef PolynomialSet<typename Basis::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename Basis::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    static const uint16_type nbPtsPerVertex = 0;
    static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace2D = 0;
    static const uint16_type nbPtsPerFace3D = 0;
    static const uint16_type nbPtsPerFace = (nDim==2)?nbPtsPerFace2D:nbPtsPerFace3D;
    static const uint16_type nbPtsPerVolume = 0;
    static const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+reference_convex_type::numEdges*nbPtsPerEdge );

    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = 0;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume = 0;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    //static const uint16_type nLocalDof = numPoints;
    static const uint16_type nLocalDof = reference_convex_type::numEdges*nDofPerEdge;

    static const uint16_type nFacesInConvex = mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<1> >,
                                                        mpl::int_<reference_convex_type::numVertices>,
                                                        typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                                                          mpl::int_<reference_convex_type::numEdges>,
                                                                          mpl::int_<reference_convex_type::numGeometricFaces> >::type >::type::value;

    NedelecDualSecondKind( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_pts_per_face( nFacesInConvex ),
        M_fset( primal )
    {
        VLOG(1) << "Nedelec finite element(dual): \n";
        VLOG(1) << " o- dim   = " << nDim << "\n";
        VLOG(1) << " o- order = " << nOrder << "\n";
        VLOG(1) << " o- kind = " << static_cast<int>(kind) << "\n";
        VLOG(1) << " o- numPoints      = " << numPoints << "\n";
        VLOG(1) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        VLOG(1) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        VLOG(1) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        VLOG(1) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        VLOG(1) << " o- nLocalDof      = " << nLocalDof << "\n";

        LOG(INFO) << "Nedelec finite element(dual): \n";
        LOG(INFO) << " o- dim   = " << nDim << "\n";
        LOG(INFO) << " o- order = " << nOrder << "\n";
        LOG(INFO) << " o- kind = " << static_cast<int>(kind) << "\n";
        LOG(INFO) << " o- numPoints      = " << numPoints << "\n";
        LOG(INFO) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        LOG(INFO) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        LOG(INFO) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        LOG(INFO) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        LOG(INFO) << " o- nLocalDof      = " << nLocalDof << "\n";

        size_type nbDofPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                          mpl::int_<face_type::numEdges*nbPtsPerEdge>,
                                          mpl::int_<face_type::numTopologicalFaces*nbPtsPerEdge + nbPtsPerFace>
                                          >::type::value;

        size_type nbEdgesPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                                             mpl::int_<face_type::numEdges>,
                                             mpl::int_<face_type::numTopologicalFaces>
                                             >::type::value;

        for(auto& m : M_pts_per_face)
            m.resize(nDim,nbDofPerFace);

        //VLOG(1) << "[Nedelec2ndKind Dual] done 1\n";
        // compute  \f$ \ell_e( U ) = (U * t[e]) (edge_pts(e)) \f$

        typedef Functional<primal_space_type> functional_type;
        std::vector<functional_type> fset;

        // loop on each entity forming the convex of topological
        // dimension 1 ( the edges)
        int d=1;
        for ( int p = 0, e = M_convex_ref.entityRange( d ).begin();
                e < M_convex_ref.entityRange( d ).end();
                ++e )
        {
            // Dof coord on current edge
            points_type Gt ( M_convex_ref.makePoints( d, e ) ); //center of the edge

            for( int f = M_convex_ref.entityRange( nDim-1 ).begin(); f < M_convex_ref.entityRange( nDim-1 ).end(); ++f)
                {
                    // M_pts_per_face[f] : set of dofs on a face (face_type::numTopologicalFaces*nbPtsPerFace)
                    //for(int k=0; k<face_type::numTopologicalFaces; ++k)
                    for(int k=0; k<nbEdgesPerFace; ++k)
                        {
                            int curEdge;
                            if( nDim <= 2 )
                                curEdge = convex_type::f2e(f,f);
                            else
                                curEdge = convex_type::f2e(f,k);

                            if( curEdge == e )
                                {
                                    ublas::subrange( M_pts_per_face[f], 0, nDim, nbPtsPerEdge*k, nbPtsPerEdge*(k+1) ) = Gt;
                                }
                        }
                }

            // M_pts : set of dofs on an element
            if ( Gt.size2() )
            {
                VLOG(1) << "Gt = " << Gt << "\n";
                ublas::subrange( M_pts, 0, nDim, p, p+Gt.size2() ) = Gt;
                p+=Gt.size2();
            }

            //typedef Feel::functional::DirectionalComponentPointsEvaluation<primal_space_type> dcpe_type;
            typedef Feel::functional::SecondDirectionalComponentPointEvaluation<primal_space_type> dcpe2_type;
            //std::cout << "tangent " << e << ":" << M_convex_ref.tangent( e ) << "\n";
            node_type dir= M_convex_ref.tangent(e);

            // Find extrem points of e (current edge)
            points_type Gt_ext(2,nDim);
            ublas::row(Gt_ext,0) = M_convex_ref.edgeVertex( e, 0 );
            ublas::row(Gt_ext,1) = M_convex_ref.edgeVertex( e, 1 );

            //1 : Find basis functions q0 et q1 of P_k(e)
            std::vector< ublas::vector<value_type> > q(2);

            ublas::matrix<value_type> Ax(2,2);
            ublas::matrix<value_type> Ay(2,2);
            ublas::column(Ax,0) = ublas::column(Gt_ext,0); // When q depends only on x coord
            ublas::column(Ay,0) = ublas::column(Gt_ext,1); // When q depends only on y coord

            ublas::matrix<value_type> solx(2,2);
            ublas::matrix<value_type> soly(2,2);

            for(int i=0; i<2; i++)
                {
                    Ax(i,Ax.size2() - 1) = 1;
                    Ay(i,Ay.size2() - 1) = 1;
                }

            // Solve the systems
            LU<ublas::matrix<value_type> > lu_x( Ax );
            LU<ublas::matrix<value_type> > lu_y( Ay );
            ublas::matrix<value_type> B = ublas::identity_matrix<value_type>( 2,2 );
            solx = lu_x.solve( B );
            soly = lu_y.solve( B );

            if(Gt_ext(0,0) == Gt_ext(1,0)) // x is constant on the edge
                {
                    for(int i=0; i<q.size(); ++i)
                        {
                            q[i].resize(3);
                            q[i](0) = 0.0;
                            ublas::project( q[i], ublas::range(1,q[i].size()) ) = ublas::column(soly,i);
                        }
                }
            else if(Gt_ext(0,1) == Gt_ext(1,1)) // y is constant on the edge
                {
                    for(int i=0; i<q.size(); ++i)
                        {
                            q[i].resize(3);
                            q[i](0) = solx(0,i);
                            q[i](1) = 0.0;
                            q[i](2) = solx(1,i);
                        }
                }
            else
                {
                    for(int i=0; i<q.size(); ++i)
                        {
                            q[i].resize(3);
                            q[i](0) = solx(0,i);
                            q[i](1) = soly(0,i);
                            q[i](2) = 0.5*solx(1,i) + 0.5*soly(1,i);
                        }
                }
            LOG(INFO) << "q0 = " << q[0] << std::endl;
            LOG(INFO) << "q1 = " << q[1] << std::endl;

            // Compute and add associated functionnals ( \int (u.t) q[i], i=0,1 )
            // dcpe2_type __dcpe2_0( primal, dir, ublas::trans(Gt_ext), q[0]);
            // dcpe2_type __dcpe2_1( primal, dir, ublas::trans(Gt_ext), q[1]);
            dcpe2_type __dcpe2_0( primal, dir, Gt, q[0]);
            dcpe2_type __dcpe2_1( primal, dir, Gt, q[1]);

            std::copy( __dcpe2_0.begin(), __dcpe2_0.end(), std::back_inserter( fset ) );
            std::copy( __dcpe2_1.begin(), __dcpe2_1.end(), std::back_inserter( fset ) );
        }

#if 0
        // add integrals of tangential component
        if ( nOrder-1 > 0 )
        {
        }
        //VLOG(1) << "[Nedelec2ndKind Dual] done 2" << std::endl;
        // add integral of moments
        if ( nOrder-2 > 0 )
        {
            // we need more equations : add interior moment
            // indeed the space is orthogonal to Pk-1
            uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
            uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
            uint16_type dim_Pm1 = convex_type::polyDims( nOrder-2 );

            Pkp1_v_type Pkp1;

            //TODO : Nedelec 2nd kind uses Raviart-Thomas functionspace instead of Pkm1 (Dkm1)
            vectorial_polynomialset_type Pkm1 ( Pkp1.polynomialsUpToDimension( dim_Pm1 ) );

            VLOG(1) << "Pkm1 = " << Pkm1.coeff() << "\n";
            VLOG(1) << "Primal = " << primal.coeff() << "\n";
            std::cout << "[loop to fill fset] Pkm1.polynomialDimension() = " << Pkm1.polynomialDimension() << std::endl;
            for ( int i = 0; i < Pkm1.polynomialDimension(); ++i )
            {
                typedef functional::IntegralMoment<primal_space_type, vectorial_polynomialset_type> fim_type;
                //typedef functional::IntegralMoment<Pkp1_v_type, vectorial_polynomialset_type> fim_type;
                VLOG(1) << "P(" << i << ")=" << Pkm1.polynomial( i ).coeff() << "\n";
                fset.push_back( fim_type( primal, Pkm1.polynomial( i ) ) );
            }
        }
#endif

        VLOG(1) << "[Nedelec2ndKind Dual] done 3, n fset = " << fset.size() << std::endl;
        LOG(INFO) << "[Nedelec2ndKind Dual] done 3, n fset = " << fset.size() << std::endl;
        M_fset.setFunctionalSet( fset );
        VLOG(1) << "[Nedelec2ndKind DUAL matrix] mat = " << M_fset.rep() << "\n";
        LOG(INFO) << "[Nedelec2ndKind DUAL matrix] mat = " << M_fset.rep() << "\n";
        VLOG(1) << "[Nedelec2ndKind Dual] done 4\n";

    }

    /**
     * \return the point set
     */
    //pointset_type const& pointSet() const { return M_pset; }

    points_type const& points() const
    {
        LOG(INFO) << "ned2 points =" << M_pts << std::endl;
        return M_pts;
    }


    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }

    points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_pts_per_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_pts_per_face[f], __i );
    }

private:
    /**
     * set the pointset at face \c f using points \c n
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_pts_per_face[f].resize( n.size1(), n.size2(), false );
        M_pts_per_face[f] = n;
    }

private:
    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    FunctionalSet<primal_space_type> M_fset;


};

}// detail

template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED2,
         typename T = double,
         uint16_type TheTAG = 0 >
struct NedelecBase
{
    typedef typename mpl::if_<mpl::bool_<(Kind == NedelecKind::NED2)>,
                              FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex>,
                                            fem::detail::NedelecDualSecondKind, PointSetEquiSpaced >,
                              // FiniteElement<NedelecPolynomialSet<N, O, NedelecKind::NED2, T>,
                              //               fem::detail::NedelecDualSecondKind,
                              //               PointSetEquiSpaced >,
                              FiniteElement<NedelecPolynomialSet<N, O, NedelecKind::NED1, T>,
                                            fem::detail::NedelecDualFirstKind, PointSetEquiSpaced > >::type type;
};

/**
 * \class Nedelec
 * \brief Nedelec Finite Element
 *
 * \f$ H(curl)\f$  conforming element
 * -\tparam  N topological dimension
 * -\tparam  O polynomial order
 * -\tparam  Kind Nedelec finite element : 1st or 2nd kind
 * -\tparam  T numerical type
 * -\tparam  TheTAG identifier for the finite element
 *
 * @author Christophe Prud'homme
 */
template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED2,
         typename T = double,
         uint16_type TheTAG = 0 >
class Nedelec
    :
    public HCurlPolynomialSet,
    public NedelecBase<N,O,Kind,T,TheTAG>::type,
    public boost::enable_shared_from_this<Nedelec<N,O,Kind,T,TheTAG> >
{

public:
    typedef typename NedelecBase<N,O,Kind,T,TheTAG>::type super;

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    //static const bool isTransformationEquivalent = false;
    static const bool isTransformationEquivalent = true;
    static const bool isContinuous = true;
    static const NedelecKind kind = Kind;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef Continuous continuity_type;
    static const uint16_type TAG = TheTAG;

    //static const polynomial_transformation_type transformation = POLYNOMIAL_CONTEXT_NEEDS_1ST_PIOLA_TRANSFORMATION;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const bool is_product = false;


    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;

    static const uint16_type nOrder =  dual_space_type::nOrder;
    static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = dual_space_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = dual_space_type::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = dual_space_type::nbPtsPerVolume;
    static const uint16_type numPoints = dual_space_type::numPoints;

    static const uint16_type nLocalDof = dual_space_type::nLocalDof;
    static const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    static const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    static const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    static const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    static const uint16_type nLocalFaceDof = ( face_type::numVertices * nDofPerVertex +
                                               face_type::numEdges * nDofPerEdge +
                                               face_type::numFaces * nDofPerFace );

    template<int subN>
    struct SubSpace
    {
        typedef Nedelec<N-1, O, kind, T, TheTAG> type;
    };

    struct SSpace
    {
        typedef Nedelec<N, O, kind, T, TheTAG> type;

    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef Nedelec<NewDim, O, kind, T,  TheTAG> type;
    };

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Nedelec()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {
#if 1
        VLOG(1) << "[N] nPtsPerEdge = " << nbPtsPerEdge << "\n";
        VLOG(1) << "[N] nPtsPerFace = " << nbPtsPerFace << "\n";
        VLOG(1) << "[N] numPoints = " << numPoints << "\n";

        VLOG(1) << "[N] nDof = " << super::nDof << "\n";

        VLOG(1) << "[N] coeff : " << this->coeff() << "\n";
        VLOG(1) << "[N] pts : " << this->points() << "\n";
        VLOG(1) << "[N] eval at pts : " << this->evaluate( this->points() ) << "\n";
        VLOG(1) << "[N] is_product : " << is_product << "\n";
#endif
    }

    Nedelec( Nedelec const & cr )
        :
        super( cr ),
        M_refconvex()
    {}

    ~Nedelec()
    {}

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

    std::string familyName() const
    {
        return "nedelec";
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
    template<typename ExprType>
    void
    getEdgeTangent(ExprType& expr, int edgeId, ublas::vector<value_type>& t) const
    {
        auto g = expr.geom();
        auto const& K = g->K(0);

        ublas::axpy_prod( K,
                          g->geometricMapping()->referenceConvex().tangent( edgeId ),
                          t,
                          true );
        t *= g->element().hEdge( edgeId )/ublas::norm_2(t);
        //std::cout << "[getEdgeTangent] t = " << t << std::endl;
    }

    typedef Eigen::MatrixXd local_interpolant_type;
    local_interpolant_type
    localInterpolant() const
        {
            return local_interpolant_type::Zero( nLocalDof, 1 );
        }

    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();
            auto g = expr.geom();
            ublas::vector<value_type> t( nDim );

            // edge dof
            for( int e = 0; e < convex_type::numEdges; ++e )
            {
                getEdgeTangent(expr, e, t);

                for ( int l = 0; l < nDofPerEdge; ++l )
                {
                    int q = e*nDofPerEdge+l;
                    for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                        Ihloc(q) += expr.evalq( c1, 0, q )*t(c1);
                }
            }
            // internal dof: internal moment of the expression against
            // polynomials of degree k-1
            if ( nOrder-2 >= 0 )
            {
                int first_dof = convex_type::numEdges*nDofPerEdge+convex_type::numGeometricFaces*nDofPerFace;
                static const uint16_type Or  = (nOrder-2>0)?nOrder-2:0;
                Feel::detail::OrthonormalPolynomialSet<nDim,Or,nDim,Vectorial,value_type,0,Simplex> p;
                IM<nDim,2*(nOrder+1),value_type,Simplex> im;
                auto Pkm1_at_quad_pts = p.evaluate( im.points() );
#if 0
                // we have to evaluate the expression at the quadrature points
                // 1- retrieve the geometric mapping
                // 2- build precompute type for geomap
                // 3- evaluate expression at quadrature points
#endif
                // 4- evaluate functionals
                const uint16_type n_internal_dof = (nDim==2)?nDofPerFace:nDofPerVolume;
                for ( int l = 0; l < n_internal_dof; ++l )
                {
                    int dof = first_dof+l;
                    
                    for( int q=0;q < im.nPoints(); ++q  )
                    {
                        value_type scal = 0.;
                        for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                            scal += expr.evalq( c1, 0, q )*Pkm1_at_quad_pts(nComponents*l+c1,q);
                        Ihloc(dof) += im.weight(q)*scal;
                    }
                }
            }
        }
    local_interpolant_type
    faceLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nLocalFaceDof, 1 );
        }
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();
            ublas::vector<value_type> t( nDim );
            auto g = expr.geom();

            for( int e = 0; e < face_type::numEdges; ++e )
            {
                int edgeid_in_element;
                if( nDim <= 2 )
                    edgeid_in_element = g->element().fToE( g->faceId(), g->faceId());
                else
                    edgeid_in_element = g->element().fToE( g->faceId(), e);
                //std::cout << "face id:  " << g->faceId() << " einf :" << e << " edge id in element : " << edgeid_in_element << std::endl;
                getEdgeTangent( expr, edgeid_in_element, t );

                for ( int l = 0; l < nDofPerEdge; ++l )
                {
                    int q = e*nDofPerEdge+l;
                    for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                        {
                            Ihloc(q) += expr.evalq( c1, 0, q )*t(c1);
                        }
                }
            }

        }

    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        typedef typename ExprType::tensor_expr_type::expression_type::fe_type fe_expr_type;
        Ihloc.setZero();
        //auto g = expr.geom();
        ublas::vector<value_type> t( nDim );

        for( int e = 0; e < convex_type::numEdges; ++e )
        {
            getEdgeTangent(expr, e, t);

            for ( int l = 0; l < nDofPerEdge; ++l )
            {
                int q = e*nDofPerEdge+l;
                for( int i = 0; i < fe_expr_type::nLocalDof; ++i )
                    for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    {
                        int ldof = c1*fe_expr_type::nLocalDof + i;
                        int ldof2 = (fe_expr_type::is_product)? ldof : i;
                        Ihloc(ldof,q) = expr.evaliq( ldof2, c1, 0, q )*t(c1);
                    }
            }
        }

    }


    //@}



protected:
    reference_convex_type M_refconvex;
private:

};
template<uint16_type N,
         uint16_type O,
         NedelecKind Kind,
         typename T,
         uint16_type TheTAG >
const uint16_type Nedelec<N,O,Kind,T,TheTAG>::nDim;
template<uint16_type N,
         uint16_type O,
         NedelecKind Kind,
         typename T,
         uint16_type TheTAG >
const uint16_type Nedelec<N,O,Kind,T,TheTAG>::nOrder;
template<uint16_type N,
         uint16_type O,
         NedelecKind Kind,
         typename T,
         uint16_type TheTAG >
const uint16_type Nedelec<N,O,Kind,T,TheTAG>::nLocalDof;

} // fem
template<uint16_type Order,
         NedelecKind Kind=NedelecKind::NED1,
         uint16_type TheTAG=0>
class Nedelec
{
public:
    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex=Simplex<N>>
    struct apply
    {
        typedef fem::Nedelec<N,Order,Kind,T,TheTAG> result_type;
        typedef result_type type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef Nedelec<Order,Kind,TheNewTAG> type;
    };

    //typedef Lagrange<Order,Scalar> component_basis_type;
    typedef Lagrange<Order+1,Scalar> component_basis_type;

    static const uint16_type nOrder =  Order;
    static const NedelecKind kind =  Kind;
    static const uint16_type TAG = TheTAG;
};
template<uint16_type Order,
         NedelecKind Kind,
         uint16_type TheTAG>
    const uint16_type Nedelec<Order,Kind,TheTAG>::nOrder;

} // Feel
#endif /* FEELPP_NEDELEC_HPP */
