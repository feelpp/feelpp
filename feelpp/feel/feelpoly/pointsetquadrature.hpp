/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Universite Joseph Fourier
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
#ifndef FEELPP_POINTSETQUADRATURE_HPP
#define FEELPP_POINTSETQUADRATURE_HPP 1

#include <feel/feelmesh/pointset.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/imfactory.hpp>

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
 * @brief Quadrature point set base class
 *
 *
 * @author Gilles Steiner <gilles.steiner@epfl.ch>
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<class Convex, typename T, typename IndexT = uint32_type>
class PointSetQuadrature : public PointSet<Convex,T>
{
public :

    static inline const bool is_face_im = false;
    typedef T value_type;
    using index_type = IndexT;
    typedef PointSet<Convex,value_type> super;
    typedef typename super::return_type return_type;
    typedef typename super::node_type node_type;
    typedef typename super::nodes_type nodes_type;
    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> vector_type;
    typedef ublas::vector<value_type> weights_type;
    static inline const uint16_type nDim = Convex::nDim;
    typedef PointSetQuadrature<Convex, T> self_type;

    typedef self_type parent_quadrature_type;
    using quad_type = IMBase<value_type>;

    PointSetQuadrature(): super(), M_quad(), M_w(), M_prod(), M_exprq() {}

    PointSetQuadrature( const PointSetQuadrature& Qp ) = default;

    explicit PointSetQuadrature( uint16_type order )
        : super( order ),
          M_order( order ),
          M_name( (boost::format("im(%1%,%2%,%3%)")%nDim %order%Convex::type() ).str() ),
          M_quad(),
          M_w(), M_w_sum(0), M_prod(), M_exprq()
        {
            if ( nDim > 0 )
            {
                std::string imNameMaxOrder =  (boost::format("im(%1%,%2%)")%nDim %Convex::type() ).str();
                auto itFindMaxOrder = IMMaxOrderFactory<value_type>::instance().find( imNameMaxOrder );
                CHECK( itFindMaxOrder != IMMaxOrderFactory<value_type>::instance().end() ) << "im type not found with " << imNameMaxOrder;
                uint16_type maxOrder = itFindMaxOrder->second.maxOrder();
                if ( M_order > maxOrder )
                {
                    LOG(INFO) << "quadrature order " << M_order << " is too big, change to max defined which is " << maxOrder;
                    M_order = maxOrder;
                    M_name = (boost::format("im(%1%,%2%,%3%)")%nDim %M_order%Convex::type() ).str();
                }
                DLOG(INFO) << "Quad name: " << M_name << std::endl;
                M_quad = *IMFactory<value_type>::instance().createObject( M_name );
                M_w.resize(M_quad.numberOfPoints());
                M_prod.resize( M_quad.numberOfPoints() );
                M_exprq.resize( M_quad.numberOfPoints() );
                create( M_order );
            }
        }

    PointSetQuadrature( weights_type Wts )
        :
        super( Wts.size() ),
        M_w( Wts ),
        M_w_sum( Wts.sum() ),
        M_prod( Wts.size() ),
        M_exprq( Wts.size() )
    {}

    
    ~PointSetQuadrature() override = default;

    virtual bool isFaceIm() const noexcept
    {
        return is_face_im;
    }
    constexpr uint16_type order() const noexcept { return M_order; }

    std::string const& name() const noexcept { return M_name; }

    /**
     * build a quadrature rule that integrates up to polynomial \c order included.
     */
    virtual void create( uint16_type order ) 
        {
            M_order = order;
            M_name = (boost::format("im(%1%,%2%,%3%)")%nDim %order%Convex::type() ).str();
            M_quad = *IMFactory<T>::instance().createObject( M_name );
            this->M_npoints = M_quad.numberOfPoints();
            this->M_points.resize( nDim, M_quad.numberOfPoints() );
            this->M_w.resize( M_quad.numberOfPoints() );
            for ( size_type i=0; i< M_quad.numberOfPoints(); i++ )
            {
                
                for ( int j = 0; j < nDim; ++j )
                {
                    this->M_points( j, i ) = M_quad.q[( nDim+1 )*i+j];
                }
                
                this->M_w( i ) = M_quad.q[( nDim+1 )*i+nDim];
            }
            M_w_sum = 0;
            for ( size_type i=0; i< M_quad.numberOfPoints(); i++ )
            {
                M_w_sum += this->M_w( i );
            }
        }

#if 0
    self_type& operator=( self_type const& q ) = default;
#else
    self_type& operator=( self_type const& q )
        {
            if ( this == &q )
                return *this;
            super::operator=( q );
            M_order = q.M_order;
            M_name = q.M_name;

            if ( nDim > 0 )
                M_quad = *IMFactory<T>::instance().createObject( M_name );
            else
                M_quad = q.M_quad;

#if 0
            this->M_npoints = M_quad.numberOfPoints();
            this->M_points.resize( nDim, M_quad.numberOfPoints() );
            this->M_w.resize( M_quad.numberOfPoints() );
            this->M_points = q.M_points;
            this->M_w = q.M_w;
            M_prod.resize( M_quad.numberOfPoints() );
            M_exprq.resize( M_quad.numberOfPoints() );
            M_w_sum = q.M_w_sum;
#endif
            M_w = q.M_w;
            M_w_sum = q.M_w_sum;
            M_w_face = q.M_w_face;
            M_n_face = q.M_n_face;
            M_w_edge = q.M_w_edge;
            M_n_edge = q.M_n_edge;
            M_prod = q.M_prod;
            M_exprq = q.M_exprq;

            return *this;
        }
#endif
    /** Face Quadrature Discussion **/

    std::vector<std::map<uint16_type,weights_type> > const& allfweights() const
    {
        return M_w_face;
    }

    std::vector<std::map<uint16_type,nodes_type > > const& allfpoints() const
    {
        return M_n_face;
    }

    /**
     * \return all quadrature weights of face f
     */
    weights_type const& weights( uint16_type __f, uint16_type __p = 1 ) const
    {
        if ( __f == uint16_type( -1 ) )
            return M_w;

        DCHECK( M_w_face[__f].find(__p) != M_w_face[__f].end() ) << "invalid permutation :" << __p << "\n";
        return M_w_face[__f].find(__p)->second;
    }

    /**
     * \return all quadrature coordinates of the q-th node of face f
     */
    nodes_type const& fpoints( uint16_type __f, uint16_type __p = 1 ) const
    {
        if ( __f == uint16_type( -1 ) )
            return this->points();

        DCHECK( M_n_face[__f].find(__p) != M_n_face[__f].end() ) << "invalid permutation :" << __p << "\n";
        return M_n_face[__f].find(__p)->second;
    }
    

    /**
     * \return quadrature weight of the q-th node of face f
     */
    value_type  weight( uint16_type __f,  uint32_type q, uint16_type __p = 1 ) const
    {
        DCHECK( M_w_face[__f].find(__p) != M_w_face[__f].end() ) << "invalid permutation :" << __p << "\n";
        return M_w_face[__f].find(__p)->second[q];
    }
    /**
     * \return quadrature coordinates of the q-th node of face f
     */
    ublas::matrix_column<nodes_type const>  fpoint( uint16_type __f, uint32_type __q, uint16_type __p = 1 ) const
    {
        DCHECK( M_n_face[__f].find(__p) != M_n_face[__f].end() ) << "invalid permutation :" << __p << "\n";
        return ublas::column( M_n_face[__f].find(__p)->second, __q );
    }


    //  quadrature_data_type data() const { return boost::make_tuple( this->M_n, this->M_w ); }

    weights_type const& weights() const
    {
        return M_w;
    }
    value_type weightsSum() const
    {
        return M_w_sum;
    }
    value_type const& weight( int q ) const
    {
        return M_w[q];
    }

    size_type nFaces() const
    {
        return M_n_face.size();
    }
    size_type nPointsOnFace( uint16_type __face = 0, uint16_type __p = 1 ) const
    {
        DCHECK( M_n_face[__face].find(__p) != M_n_face[__face].end() ) << "invalid permutation :" << __p << "\n";
        return M_n_face[__face].find(__p)->second.size2();
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
        DVLOG(2) << "[PointSetQuadrature] face " << int( __face )<< " integration\n";
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
            M_exprq[q] = expr.evalijq( indi, indj, c1, c2, q );
        }

        return M_prod.dot( M_exprq );
    }
    template<typename IndexTest, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           uint16_type c1,
                           uint16_type c2 ) const
    {
        for ( uint16_type q = 0; q < this->nPoints(); ++q )
        {
            M_exprq[q] = expr.evaliq( indi, c1, c2, q );
        }

        return M_prod.dot( M_exprq );
    }

    template<typename ExprT>
    value_type operator()( ExprT const& expr,
                           uint16_type c1,
                           uint16_type c2 ) const
    {
        for ( uint16_type q = 0; q < this->nPoints(); ++q )
        {
            M_exprq[q] = expr.evalq( c1, c2, q );
        }

        return M_prod.dot( M_exprq );
    }
    template<typename IndexTest, typename IndexTrial, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           IndexTrial  const& indj,
                           uint16_type c1,
                           uint16_type c2,
                           std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad ) const
    {
        value_type res = value_type( 0 );

        for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
        {
            auto qReal = indexLocalToQuad[q].template get<0>();
            const value_type val_expr = expr.evalijq( indi, indj, c1, c2, q );
            res += M_prod[qReal]*val_expr;
        }

        return res;
    }
    template<typename IndexTest, typename ExprType>
    value_type operator()( ExprType const& expr,
                           IndexTest  const& indi,
                           uint16_type c1,
                           uint16_type c2,
                           std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad ) const
    {
        value_type res = value_type( 0 );

        for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
        {
            auto qReal = indexLocalToQuad[q].template get<0>();
            const value_type val_expr = expr.evaliq( indi, c1, c2, q );
            res += M_prod[qReal]*val_expr;
        }

        return res;
    }

    template<typename GMC>
    void update( GMC&& gmc )
    {
        using gmc_t = std::decay_t<GMC>;
        // if ( this->isFaceIm() )
        if constexpr ( gmc_t::is_on_face )
        {
            if constexpr ( gmc_t::PDim == gmc_t::NDim-1 )
            {
                for ( uint16_type q = 0; q < this->nPoints(); ++q )
                {
                    M_prod[q] = M_w( q )*gmc.J( q )*gmc.normalNorm( q );
                }
            }
            else
            {
                for ( uint16_type q = 0; q < this->nPoints(); ++q )
                {
                    M_prod[q] = M_w( q ) * gmc.J( q );
                }
            }
        }
        else
        {
            for ( uint16_type q = 0; q < this->nPoints(); ++q )
            {
                M_prod[q] = M_w( q )*gmc.J( q );
            }
        }
    }

    template<typename GMC>
    void update( GMC const& gmc,
                 std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad )
    {

        if ( this->isFaceIm() )
            for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
            {
                auto qReal = indexLocalToQuad[q].template get<0>();
                M_prod[qReal] = M_w( qReal )*gmc.J( q )*gmc.normalNorm( q );
            }

        else
            for ( uint16_type q = 0; q < indexLocalToQuad.size(); ++q )
            {
                auto qReal = indexLocalToQuad[q].template get<0>();
                M_prod[qReal] = M_w( qReal )*gmc.J( q );
            }
    }


    class Face
        :
    public PointSetQuadrature<typename Convex::topological_face_type,T>
    {
        using super = PointSetQuadrature<typename Convex::topological_face_type,T>;
    public:
        using parent_quadrature_type = super;
        static inline const bool is_face_im = true;

        Face()
            :
            super(),
            M_f( invalid_uint16_type_value )
        {
        }
        Face( self_type const& quad_elt, uint16_type f = 0 )
            :
            super( quad_elt.order() ),
            M_f( f )
        {
            VLOG(2) << "Quadrature face: " << quad_elt;
            this->setPoints( quad_elt.fpoints( f ) );
            this->setWeights( quad_elt.weights( f ) );
            
        }
        Face( Face const& f )
            :
            super( f ),
            M_f( f.M_f )

        {}
        Face& operator=( Face const& f )
        {
            if ( this != &f )
            {
                super::operator=( f );
                M_f = f.M_f;
            }

            return *this;
        }
        bool isFaceIm() const noexcept override
        {
            return is_face_im;
        }
        uint16_type face() const noexcept
        {
            return M_f;
        }
    private:
        uint16_type M_f;
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
        M_prod.resize( w.size() );
        M_exprq.resize( w.size() );
        M_w = w;
    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& ref_convex,
                            std::shared_ptr<GM> const& __gm,
                            std::shared_ptr<IM> const& __qr_face )
    {
        constructQROnFace( ref_convex, __gm, __qr_face, mpl::bool_<( Convex::nDim > 1 )>() );
    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& /*ref_convex*/,
                            std::shared_ptr<GM> const& /*__gm*/,
                            std::shared_ptr<IM> const& /*__qr_face*/,
                            mpl::bool_<false> )
    {

        typedef typename Elem::permutation_type permutation_type;
        DCHECK(permutation_type::N_PERMUTATIONS==2) << " number of permutation must be equal to 2 here\n";
        //BOOST_STATIC_ASSERT( Elem::nDim == 1 );

        M_n_face.resize( Elem::numTopologicalFaces );
        M_w_face.resize( Elem::numTopologicalFaces );

        M_n_face[0][permutation_type::IDENTITY].resize( Elem::nDim, 1 );
        M_w_face[0][permutation_type::IDENTITY].resize( 1 );
        M_n_face[0][permutation_type::IDENTITY]( 0, 0 ) = -1;
        M_w_face[0][permutation_type::IDENTITY]( 0 ) = 1;

        M_n_face[1][permutation_type::IDENTITY].resize( Elem::nDim, 1 );
        M_w_face[1][permutation_type::IDENTITY].resize( 1 );
        M_n_face[1][permutation_type::IDENTITY]( 0, 0 ) = 1;
        M_w_face[1][permutation_type::IDENTITY]( 0 ) = 1;

    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnFace( Elem const& ref_convex,
                            std::shared_ptr<GM> const& __gm,
                            std::shared_ptr<IM> const& __qr_face,
                            mpl::bool_<true> );

    template<typename Elem, typename GM, typename IM>
    void constructQROnEdge( Elem const& ref_convex,
                            std::shared_ptr<GM> const& __gm,
                            std::shared_ptr<IM> const& __qr_edge )
    {
        constructQROnEdge( ref_convex, __gm, __qr_edge, mpl::bool_<( Convex::nDim == 3 )>() );
    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnEdge( Elem const& /*ref_convex*/,
                            std::shared_ptr<GM> const& /*__gm*/,
                            std::shared_ptr<IM> const& /*__qr_edge*/,
                            mpl::bool_<false> )
    {

        typedef typename Elem::permutation_type permutation_type;
        DCHECK(permutation_type::N_PERMUTATIONS==2) << " number of permutation must be equal to 2 here\n";
        BOOST_STATIC_ASSERT( Elem::nDim == 1 );

        M_n_edge.resize( Elem::numEdges );
        M_w_edge.resize( Elem::numEdges );

        M_n_edge[0][permutation_type::IDENTITY].resize( Elem::nDim, 1 );
        M_w_edge[0][permutation_type::IDENTITY].resize( 1 );
        M_n_edge[0][permutation_type::IDENTITY]( 0, 0 ) = -1;
        M_w_edge[0][permutation_type::IDENTITY]( 0 ) = 1;

        M_n_edge[1][permutation_type::IDENTITY].resize( Elem::nDim, 1 );
        M_w_edge[1][permutation_type::IDENTITY].resize( 1 );
        M_n_edge[1][permutation_type::IDENTITY]( 0, 0 ) = 1;
        M_w_edge[1][permutation_type::IDENTITY]( 0 ) = 1;

    }

    template<typename Elem, typename GM, typename IM>
    void constructQROnEdge( Elem const& ref_convex,
                            std::shared_ptr<GM> const& __gm,
                            std::shared_ptr<IM> const& __qr_edge,
                            mpl::bool_<true> );

protected:

    uint16_type M_order;

    std::string M_name;

    //std::unique_ptr<quad_type> M_quad;
    quad_type M_quad;

    weights_type M_w;
    value_type M_w_sum;
    
    std::vector<std::map<uint16_type,weights_type> > M_w_face;
    std::vector<std::map<uint16_type,nodes_type > > M_n_face;

    std::vector<std::map<uint16_type,weights_type> > M_w_edge;
    std::vector<std::map<uint16_type,nodes_type > > M_n_edge;

    vector_type M_prod;
    mutable vector_type M_exprq;


};

template<class Convex, typename T, typename IndexT>
template<typename Elem, typename GM, typename IM>
void
PointSetQuadrature<Convex,T,IndexT>::constructQROnFace( Elem const& ref_convex,
        std::shared_ptr<GM> const& __gm,
        std::shared_ptr<IM> const& __qr_face,
        mpl::bool_<true>  )
{
    M_n_face.resize( Elem::numTopologicalFaces );
    M_w_face.resize( Elem::numTopologicalFaces );

    typedef typename Elem::permutation_type permutation_type;
    // FIXME: we don't handle the case where faces are not of the
    // same type like for prism or pyramids. In that case geopc
    // should be defined per face
    typedef typename GM::face_gm_type::precompute_type face_pc_type;
    typedef typename GM::face_gm_type::precompute_ptrtype face_pc_ptrtype;
    face_pc_ptrtype __geopc( new face_pc_type( __gm->boundaryMap(),__qr_face->points() ) );
    
    for ( uint16_type __f = 0; __f < Elem::numTopologicalFaces; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                typedef typename Elem::topological_face_type element_type;
                // get ref convex of face with a permutation given
                element_type ref_convex_face = ref_convex.topologicalFace( __f, __p.value() );
                DVLOG(2) << "[quadpt] ref_convex_face "  << __f << "=" << ref_convex_face.points() << "\n";
                //toPython( ref_convex_face );

                auto __c = __gm->boundaryMap()->template context<vm::JACOBIAN|vm::POINT|vm::KB>( ref_convex_face,
                                                                                                 __geopc );
                __c->template update<vm::JACOBIAN|vm::POINT|vm::KB>( ref_convex_face );
                DVLOG(2) << "[quadpt] ref_convex_face "  << __f << " xref" << __c->xRefs() << "\n";
                DVLOG(2) << "[quadpt] ref_convex_face "  << __f << " xreal" << __c->xReal() << "\n";

                value_type __len = 0.0;
                M_n_face[__f][__p.value()].resize( Elem::nDim, __qr_face->nPoints() );
                M_w_face[__f][__p.value()].resize( __qr_face->nPoints() );

                DVLOG(2) << "[PointSetQuadrature::constructQROnFace] npoints on face "
                         << __f << " : "
                         << __qr_face->nPoints() << "\n";
                // transform the quad nodes on the boundary _reference_
                // element to the face of the reference element face

                for ( uint16_type __ip = 0; __ip < __qr_face->nPoints(); ++__ip )
                {
                    ublas::column( M_n_face[__f][__p.value()], __ip ) = __c->xReal( __ip );

                    // w = w_face * ||B*n|| * J
                    M_w_face[ __f][__p.value()]( __ip ) = __qr_face->weight( __ip )*__c->J( __ip );

                    __len += M_w_face[ __f][__p.value()]( __ip );


                    DVLOG(2) << "face " << __f << " ip = " << __ip << "       J =" << __c->J( __ip ) << "\n";
                    DVLOG(2) << "face " << __f << " ip = " << __ip << "  weight =" << __qr_face->weight( __ip ) << "\n";
                    DVLOG(2) << "face " << __f << " ip = " << __ip << "  weight =" << M_w_face[ __f][__p.value()]( __ip ) << "\n";
                    //            DVLOG(2) << "face " << __f << " ip = " << __ip << "  x  ref =" << __c->xRef( __ip ) << "\n";
                    //            DVLOG(2) << "face " << __f << " ip = " << __ip << "  x real =" << __c->xReal( __ip ) << "\n";
                }

                DVLOG(2) << "length = " << __len << "\n";
            }
    }
}


template<class Convex, typename T, typename IndexT>
template<typename Elem, typename GM, typename IM>
void
PointSetQuadrature<Convex,T, IndexT>::constructQROnEdge( Elem const& ref_convex,
                                                         std::shared_ptr<GM> const& __gm,
                                                         std::shared_ptr<IM> const& __qr_edge,
                                                         mpl::bool_<true>  )
{
    M_n_edge.resize( Elem::numEdges );
    M_w_edge.resize( Elem::numEdges );

    using permutation_type = typename Elem::edge_permutation_type;
    using edge_pc_type =  typename GM::edge_gm_type::precompute_type;
    using edge_pc_ptrtype =  typename GM::edge_gm_type::precompute_ptrtype;

    edge_pc_ptrtype __geopc = std::make_shared<edge_pc_type>( __gm->edgeMap(),__qr_edge->points() );

    for ( uint16_type __f = 0; __f < Elem::numEdges; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                using element_type =  typename Elem::edge_type;
                
                // get ref convex of edge with a permutation given
                element_type ref_convex_edge = ref_convex.edge( __f, __p.value() );
                DVLOG(2) << "[quadpt] ref_convex_edge "  << __f << "=" << ref_convex_edge.points() << "\n";
                //toPython( ref_convex_edge );

                auto ctx = __gm->context<(vm::JACOBIAN|vm::POINT|vm::KB)>( __gm->edgeMap(), ref_convex_edge, __geopc );
                ctx.template update<vm::JACOBIAN|vm::POINT|vm::KB>( ref_convex_edge );
                DVLOG(2) << "[quadpt] ref_convex_edge "  << __f << " xref" << ctx.xRefs() << "\n";
                DVLOG(2) << "[quadpt] ref_convex_edge "  << __f << " xreal" << ctx.xReal() << "\n";

                value_type __len = 0.0;
                M_n_edge[__f][__p.value()].resize( Elem::nDim, __qr_edge->nPoints() );
                M_w_edge[__f][__p.value()].resize( __qr_edge->nPoints() );

                DVLOG(2) << "[PointSetQuadrature::constructQROnEdge] npoints on edge "
                         << __f << " : "
                         << __qr_edge->nPoints() << "\n";
                // transform the quad nodes on the boundary _reference_
                // element to the edge of the reference element edge

                for ( uint16_type __ip = 0; __ip < __qr_edge->nPoints(); ++__ip )
                {
                    ublas::column( M_n_edge[__f][__p.value()], __ip ) = ctx.xReal( __ip );

                    // w = w_edge * ||B*n|| * J
                    M_w_edge[ __f][__p.value()]( __ip ) = __qr_edge->weight( __ip )*ctx.J( __ip );

                    __len += M_w_edge[ __f][__p.value()]( __ip );


                    DVLOG(2) << "edge " << __f << " ip = " << __ip << "       J =" << ctx.J( __ip ) << "\n";
                    DVLOG(2) << "edge " << __f << " ip = " << __ip << "  weight =" << __qr_edge->weight( __ip ) << "\n";
                    DVLOG(2) << "edge " << __f << " ip = " << __ip << "  weight =" << M_w_edge[ __f][__p.value()]( __ip ) << "\n";
                    //            DVLOG(2) << "edge " << __f << " ip = " << __ip << "  x  ref =" << ctx.xRef( __ip ) << "\n";
                    //            DVLOG(2) << "edge " << __f << " ip = " << __ip << "  x real =" << ctx.xReal( __ip ) << "\n";
                }

                DVLOG(2) << "length = " << __len << "\n";
            }
    }
}


} // Feel


namespace std
{
template<class Convex, typename T, typename IndexT>
std::ostream& operator<<( std::ostream& os, Feel::PointSetQuadrature<Convex,T, IndexT> const& qps )
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

#endif /* FEELPP_POINTSETQUADRATURE_HPP */
