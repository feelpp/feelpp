/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Universite Joseph Fourier

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
   \file pointsetquadrature.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
#ifndef __Quadpoint_H
#define __Quadpoint_H 1

#include <feel/feelmesh/pointset.hpp>
#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/geomap.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

enum IntegrationFaceEnum
{
    ALL_FACES = -1,
    FACE_0 = 0,
    FACE_1 = 1,
    FACE_2,
    FACE_3,
    FACE_4,
    FACE_5,
    FACE_6,
    FACE_7,
    FACE_8
};

/**
 * \class PointSetQuadrature
 * \brief Quadrature point set base class
 *
 *
 * \author Gilles Steiner <gilles.steiner@epfl.ch>
 * \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<class Convex, uint16_type Integration_Degree, typename T>
class PointSetQuadrature : public PointSet<Convex,T>
{
public :
    static const bool is_face_im = false;
    typedef T value_type;
    typedef PointSet<Convex,value_type> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vector_type;
    typedef ublas::vector<value_type> weights_type;

    typedef PointSetQuadrature<Convex, Integration_Degree, T> self_type;

    typedef self_type parent_quadrature_type;
    static const uint16_type I_deg = Integration_Degree;

    PointSetQuadrature(): super(), _M_w(), _M_prod(), _M_exprq() {}

    PointSetQuadrature( const PointSetQuadrature& Qp )
        : super( Qp ),
          _M_w( Qp.weights() ),
          _M_w_face( Qp.allfweights() ),
          _M_n_face( Qp.allfpoints() ),
          _M_prod( Qp.nPoints() ),
          _M_exprq( Qp.nPoints() )
    {}

    PointSetQuadrature( uint32_type Npoints )
        : super( Npoints ), _M_w( Npoints ), _M_prod( Npoints ), _M_exprq( Npoints )
    {}

    PointSetQuadrature( weights_type Wts )
        :
        super( Wts.size() ),
        _M_w( Wts ),
        _M_prod( Wts.size() ),
        _M_exprq( Wts.size() )
    {}


    virtual ~PointSetQuadrature()
    {}

    virtual bool isFaceIm() const
    {
        return is_face_im;
    }

    self_type& operator=( self_type const& q )
    {
        if ( this != &q )
        {
            super::operator=( q );
            _M_w = q._M_w;
            _M_w_face = q._M_w_face;
            _M_n_face = q._M_n_face;
            _M_prod = q._M_prod;
            _M_exprq = q._M_exprq;
        }

        return *this;
    }

    /** Face Quadrature Discussion **/


    std::vector<weights_type> const& allfweights() const
    {
        return _M_w_face;
    }

    std::vector<nodes_type> const& allfpoints() const
    {
        return _M_n_face;
    }


    /**
     * \return all quadrature weights of face f
     */
    weights_type const& weights( uint16_type __f ) const
    {
        if ( __f == uint16_type( -1 ) )
            return _M_w;

        return _M_w_face[__f];
    }

    /**
     * \return all quadrature coordinates of the q-th node of face f
     */
    nodes_type const& fpoints( uint16_type __f ) const
    {
        if ( __f == uint16_type( -1 ) )
            return this->points();

        return _M_n_face[__f];
    }


    /**
     * \return quadrature weight of the q-th node of face f
     */
    value_type  weight( uint16_type __f,  uint32_type q ) const
    {
        return _M_w_face[__f][q];
    }
    /**
     * \return quadrature coordinates of the q-th node of face f
     */
    ublas::matrix_column<nodes_type const>  fpoint( uint16_type __f, uint32_type __q ) const
    {
        return ublas::column( _M_n_face[__f], __q );
    }


    //  quadrature_data_type data() const { return boost::make_tuple( this->_M_n, this->_M_w ); }

    //size_type nFaces() const { return _M_n_face.size(); }
    //size_type nPointsOnFace( int __face = 0 ) const { return _M_n_face[__face].size2(); }

    weights_type const& weights() const
    {
        return _M_w;
    }
    value_type const& weight( int q ) const
    {
        return _M_w[q];
    }

    size_type nFaces() const
    {
        return _M_n_face.size();
    }
    size_type nPointsOnFace( int __face = 0 ) const
    {
        return _M_n_face[__face].size2();
    }

    /**
     * Integrate the function \p f over the convex associated with the
     * quadrature points
     *
     * This function handles the scalar and vectorial case
     */
    template<typename Expression>
    value_type integrate( Expression f ) const
    {
        //BOOST_STATIC_ASSERT( ( boost::is_same<value_type,qd_real>::value ) );
        //BOOST_MPL_ASSERT_MSG( ( boost::is_same<value_type,qd_real>::value ), INVALID_TYPE, (value_type) );
        //node_type res( this->weights()[0]*f( this->point(0) ) );
        value_type res( 0.0 );

        for ( uint32_type k = 0; k < this->nPoints(); ++k )
        {

            value_type fk = f( k );
            res += this->weights()[k]*fk;
        }

        return res;
    }

    /**
     * Integrate the function \p f over the convex associated with the
     * quadrature points
     *
     * This function handles the scalar and vectorial case
     */
    template<typename Expression>
    value_type integrateAtPoints( Expression const& f ) const
    {
        //node_type res( this->weights()[0]*f( this->point(0) ) );
        value_type res( 0.0 );

        for ( uint32_type k = 0; k < this->nPoints(); ++k )
            res += this->weights()[k]*f( this->point( k ) );

        return res;
    }

    /**
     * Integrate \c f on a face of the domain defined by the shape of
     * the element
     *
     * This function handles the scalar and vectorial case
     */
    template<typename Expression>
    value_type integrate( IntegrationFaceEnum __face,
                          Expression const& f ) const
    {
        //std::cout << "integrating using " << nNodes << "\n";
        FEELPP_ASSERT(  __face == ALL_FACES || ( __face >= 0 && __face < nFaces() ) )
        ( __face )( nFaces() ).error( "invalid face index (can be ALL_FACE or FACE_0 <= f < nFaces()" );
        value_type res( 0.0 );

        if ( __face != ALL_FACES )
        {
            // integrate on face f
            for ( uint32_type k = 0; k < this->nPointsOnFace( __face ); ++k )
                res += this->weight( __face, k )*f( k );
        }

        else
        {
            // integrate on all faces
            for ( uint16_type i = 0; i < this->nFaces(); ++i )
            {
                value_type __res_face( 0.0 );

                for ( uint32_type k = 0; k < this->nPointsOnFace( i ); ++k )
                    __res_face += this->weight( i, k )*f( k );

                //std::cout << "res on face " << i << " = " << __res_face << "\n";
                res += __res_face;
            }
        }

        return res;
    }
    /**
     * Integrate \c f on a face of the domain defined by the shape of
     * the element
     *
     * This function handles the scalar and vectorial case
     */
    template<typename Expression>
    value_type integrateAtPoints( IntegrationFaceEnum __face,
                                  Expression const& f ) const
    {
        Debug( 5048 ) << "[PointSetQuadrature] face " << int( __face )<< " integration\n";
        //std::cout << "integrating using " << nNodes << "\n";
        FEELPP_ASSERT(  int( __face ) == int( ALL_FACES ) ||
                        ( int( __face ) >= 0 &&
                          int( __face ) < int( nFaces() ) ) )
        ( __face )( nFaces() ).error( "invalid face index (can be ALL_FACE or FACE_0 <= f < nFaces()" );
        value_type res( 0.0 );

        if ( __face != ALL_FACES )
        {
            // integrate on face f
            for ( uint32_type k = 0; k < this->nPointsOnFace( __face ); ++k )
                res += this->weight( __face, k )*f( this->fpoint( __face, k ) );
        }

        else
        {
            // integrate on all faces
            for ( uint16_type i = 0; i < this->nFaces(); ++i )
            {
                value_type __res_face( 0.0 );

                for ( uint32_type k = 0; k < this->nPointsOnFace( i ); ++k )
                    __res_face += this->weight( i, k )*f( this->fpoint( i, k ) );

                //std::cout << "res on face " << i << " = " << __res_face << "\n";
                res += __res_face;
            }
        }

        return res;
    }

    template<typename IndexTest, typename IndexTrial, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           IndexTrial  const& indj,
                           uint16_type c1,
                           uint16_type c2 ) const
    {
        for ( uint16_type q = 0; q < this->nPoints(); ++q )
        {
            _M_exprq[q] = expr.evalijq( indi, indj, c1, c2, q );
        }

        return _M_prod.dot( _M_exprq );
    }
    template<typename IndexTest, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           uint16_type c1,
                           uint16_type c2 ) const
    {
        for ( uint16_type q = 0; q < this->nPoints(); ++q )
        {
            _M_exprq[q] = expr.evaliq( indi, c1, c2, q );
        }

        return _M_prod.dot( _M_exprq );
    }

    template<typename ExprT>
    value_type operator()( ExprT const& expr,
                           uint16_type c1,
                           uint16_type c2 ) const
    {
        for ( uint16_type q = 0; q < this->nPoints(); ++q )
        {
            _M_exprq[q] = expr.evalq( c1, c2, q );
        }

        return _M_prod.dot( _M_exprq );
    }
    template<typename IndexTest, typename IndexTrial, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           IndexTrial  const& indj,
                           uint16_type c1,
                           uint16_type c2,
                           std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad ) const
    {
        value_type res = value_type( 0 );

        for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
        {
            auto qReal = indexLocalToQuad[q].get<0>();
            const value_type val_expr = expr.evalijq( indi, indj, c1, c2, q );
            res += _M_prod[qReal]*val_expr;
        }

        return res;
    }
    template<typename IndexTest, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           uint16_type c1,
                           uint16_type c2,
                           std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad ) const
    {
        value_type res = value_type( 0 );

        for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
        {
            auto qReal = indexLocalToQuad[q].get<0>();
            const value_type val_expr = expr.evaliq( indi, c1, c2, q );
            res += _M_prod[qReal]*val_expr;
        }

        return res;
    }

    template<typename GMC>
    void update( GMC const& gmc )
    {
        if ( this->isFaceIm() )
            for ( uint16_type q = 0; q < this->nPoints(); ++q )
            {
                _M_prod[q] = _M_w( q )*gmc.J( q )*gmc.normalNorm( q );
            }

        else
            for ( uint16_type q = 0; q < this->nPoints(); ++q )
            {
                _M_prod[q] = _M_w( q )*gmc.J( q );
            }
    }

    template<typename GMC>
    void update( GMC const& gmc,
                 std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad )
    {

        if ( this->isFaceIm() )
            for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
            {
                auto qReal = indexLocalToQuad[q].get<0>();
                _M_prod[qReal] = _M_w( qReal )*gmc.J( q )*gmc.normalNorm( q );
            }

        else
            for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
            {
                auto qReal = indexLocalToQuad[q].get<0>();
                _M_prod[qReal] = _M_w( qReal )*gmc.J( q );
            }
    }


    class Face
        :
    public PointSetQuadrature<typename Convex::topological_face_type,
        Integration_Degree,
        T>
    {
        typedef PointSetQuadrature<typename Convex::topological_face_type,
                Integration_Degree,
                T> super;
    public:
        typedef super parent_quadrature_type;
        static const bool is_face_im = true;
        Face()
            :
            super(),
            _M_f( invalid_uint16_type_value )
        {
        }
        Face( self_type const& quad_elt, uint16_type f = 0 )
            :
            super(),
            _M_f( f )
        {
            this->setPoints( quad_elt.fpoints( f ) );
            this->setWeights( quad_elt.weights( f ) );

        }
        Face( Face const& f )
            :
            super( f ),
            _M_f( f._M_f )

        {}
        Face& operator=( Face const& f )
        {
            if ( this != &f )
            {
                super::operator=( f );
                _M_f = f._M_f;
            }

            return *this;
        }
        bool isFaceIm() const
        {
            return is_face_im;
        }
        uint16_type face() const
        {
            return _M_f;
        }
    private:
        uint16_type _M_f;
    };

    typedef Face face_quadrature_type;
    face_quadrature_type face( uint16_type f ) const
    {
        return Face( *this, f );
    }

    FEELPP_DEFINE_VISITABLE();

protected:

    void setWeights( weights_type const& w )
    {
        _M_prod.resize( w.size() );
        _M_exprq.resize( w.size() );
        _M_w = w;
    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& ref_convex,
                            boost::shared_ptr<GM> const& __gm,
                            boost::shared_ptr<IM> const& __qr_face )
    {
        constructQROnFace( ref_convex, __gm, __qr_face, mpl::bool_<( Convex::nDim > 1 )>() );
    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& /*ref_convex*/,
                            boost::shared_ptr<GM> const& /*__gm*/,
                            boost::shared_ptr<IM> const& /*__qr_face*/,
                            mpl::bool_<false> )
    {
        BOOST_STATIC_ASSERT( Elem::nDim == 1 );
        _M_n_face.resize( Elem::numTopologicalFaces );
        _M_w_face.resize( Elem::numTopologicalFaces );

        _M_n_face[0].resize( Elem::nDim, 1 );
        _M_w_face[0].resize( 1 );
        _M_n_face[0]( 0, 0 ) = -1;
        _M_w_face[0]( 0 ) = 1;

        _M_n_face[1].resize( Elem::nDim, 1 );
        _M_w_face[1].resize( 1 );
        _M_n_face[1]( 0, 0 ) = 1;
        _M_w_face[1]( 0 ) = 1;

    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& ref_convex,
                            boost::shared_ptr<GM> const& __gm,
                            boost::shared_ptr<IM> const& __qr_face,
                            mpl::bool_<true> );

protected:

    weights_type _M_w;

    std::vector<weights_type> _M_w_face;
    std::vector<nodes_type> _M_n_face;

    vector_type _M_prod;
    mutable vector_type _M_exprq;

};

template<class Convex, uint16_type Integration_Degree, typename T>
template<typename Elem, typename GM, typename IM>
void
PointSetQuadrature<Convex,Integration_Degree,T>::constructQROnFace( Elem const& ref_convex,
        boost::shared_ptr<GM> const& __gm,
        boost::shared_ptr<IM> const& __qr_face,
        mpl::bool_<true>  )
{
    _M_n_face.resize( Elem::numTopologicalFaces );
    _M_w_face.resize( Elem::numTopologicalFaces );

    // FIXME: we don't handle the case where faces are not of the
    // same type like for prism or pyramids. In that case geopc
    // should be defined per face
    typedef typename GM::face_gm_type::precompute_type face_pc_type;
    typedef typename GM::face_gm_type::precompute_ptrtype face_pc_ptrtype;
    face_pc_ptrtype __geopc( new face_pc_type( __gm->boundaryMap(),__qr_face->points() ) );

    for ( uint16_type __f = 0; __f < Elem::numTopologicalFaces; ++__f )
    {
        typedef typename Elem::topological_face_type  element_type;
        element_type ref_convex_face = ref_convex.topologicalFace( __f );
        Debug( 5048 ) << "[quadpt] ref_convex_face "  << __f << "=" << ref_convex_face.points() << "\n";
        //toPython( ref_convex_face );

        typename GM::template Context<vm::JACOBIAN|vm::POINT|vm::KB,element_type> __c( __gm->boundaryMap(),
                ref_convex_face,
                __geopc );
        __c.update( ref_convex_face );
        Debug( 5048 ) << "[quadpt] ref_convex_face "  << __f << " xref" << __c.xRefs() << "\n";
        Debug( 5048 ) << "[quadpt] ref_convex_face "  << __f << " xreal" << __c.xReal() << "\n";

        value_type __len = 0.0;
        _M_n_face[__f].resize( Elem::nDim, __qr_face->nPoints() );
        _M_w_face[__f].resize( __qr_face->nPoints() );
        Debug( 5048 ) << "[PointSetQuadrature::constructQROnFace] npoints on face "
                      << __f << " : "
                      << __qr_face->nPoints() << "\n";
        // transform the quad nodes on the boundary _reference_
        // element to the face of the reference element face

        for ( uint16_type __ip = 0; __ip < __qr_face->nPoints(); ++__ip )
        {
            ublas::column( _M_n_face[__f], __ip ) = __c.xReal( __ip );

            // w = w_face * ||B*n|| * J
            _M_w_face[ __f]( __ip ) = __qr_face->weight( __ip )*__c.J( __ip );

            __len += _M_w_face[ __f]( __ip );
            Debug( 5048 ) << "face " << __f << " ip = " << __ip << "       J =" << __c.J( __ip ) << "\n";
            Debug( 5048 ) << "face " << __f << " ip = " << __ip << "  weight =" << __qr_face->weight( __ip ) << "\n";
            Debug( 5048 ) << "face " << __f << " ip = " << __ip << "  weight =" << _M_w_face[ __f]( __ip ) << "\n";
            //            Debug( 5048 ) << "face " << __f << " ip = " << __ip << "  x  ref =" << __c.xRef( __ip ) << "\n";
            //            Debug( 5048 ) << "face " << __f << " ip = " << __ip << "  x real =" << __c.xReal( __ip ) << "\n";
        }

        Debug( 5048 ) << "length = " << __len << "\n";
    }
}


} // Feel


namespace std
{
template<class Convex, Feel::uint16_type Integration_Degree, typename T>
std::ostream& operator<<( std::ostream& os,
                          Feel::PointSetQuadrature<Convex,Integration_Degree,T> const& qps )
{
    os << "quadrature point set:\n"
       << "number of points: " << qps.nPoints() << "\n"
       << "points :  " << qps.points() << "\n"
       << "weights :  " << qps.weights() << "\n";

    for ( Feel::uint16_type f = 0; f < qps.nFaces(); ++f )
    {
        os << " o Face " << f << "\n";
        os << qps.fpoints( f ) << "\n";
        os << qps.weights( f ) << "\n";
    }

    return os;
}

}

#endif /* _QuadPoint_H */
