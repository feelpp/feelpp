/*
  This file is part of the Feel library

  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Université Joseph Fourier Grenoble 1

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
#ifndef _GEOMAP_H
#define _GEOMAP_H

#include <boost/version.hpp>
#if ( BOOST_VERSION >= 103400 )
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */


#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Geometry>
#include <Eigen/src/Core/util/ReenableStupidWarnings.h>

#include <boost/functional/hash.hpp>
#include <memory>
#include <boost/mpl/vector.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>

//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/atlas/cblas3.hpp>

#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/iteration.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/svd.hpp>

#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/entitymarkers.hpp>
#include <feel/feelpoly/context.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feelpoly/fekete.hpp>

#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/lagrange.hpp>


namespace Feel
{
/**
 * \enum type of geomap strategy
 */
enum class GeomapStrategyType
{
    GEOMAP_OPT = 0,
    GEOMAP_O1 = 1,
    GEOMAP_HO = 2
};
namespace detail
{
template <typename RangeType>
GeomapStrategyType
geomapStrategy( RangeType const& /**/, GeomapStrategyType gs )
{
    constexpr uint16_type geoOrder = entity_range_t<RangeType>::nOrder;
    if ( geoOrder == 1 )
        return GeomapStrategyType::GEOMAP_HO;
    else
        return gs;
}
}

//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;

struct GeomapInverse
{
    //GeomapInverse
};

template <uint16_type O,
          template <uint16_type Dim> class PolySetType,
          typename ContinuityType,
          template <class, uint16_type, class> class Pts,
          uint16_type TheTag>
class Lagrange;

/**
 * \class GeoMap
 * \brief Structure for the geometrical mapping
 * \author C. Prud'homme
 *
 * This class contains the geometrical transformation that maps the reference
 * element on the current element, and its values on integration points
 *
 */
template <uint16_type Dim,
          uint16_type Order,
          uint16_type RealDim,
          typename T = double,
          template <uint16_type, uint16_type, uint16_type> class Entity = Simplex,
          template <uint16_type, template <uint16_type RDim> class PolySetType, typename ContinuityType,
                    template <class, uint16_type, class> class Pts, uint16_type> class PP = Lagrange>
class GeoMap
    : public PP<Order, Scalar, Continuous, PointSetEquiSpaced, 0>::template apply<Dim, RealDim /*Dim*/, T, Entity<Dim, Order, /*RealDim*/ Dim>>::result_type //,
                                                                                                                                                             //public boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > >
//public PP<Order,Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type
{
    //typedef typename PP<Order, Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type super;
    typedef typename PP<Order, Scalar, Continuous, PointSetEquiSpaced, 0>::template apply<Dim, RealDim /*Dim*/, T, Entity<Dim, Order, /*RealDim*/ Dim>>::result_type super;

    //typedef boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > > super_enable_this;

    static inline const uint16_type nRealDimCheck2d = mpl::if_<mpl::less_equal<mpl::int_<2>, mpl::int_<RealDim>>,
                                                        mpl::int_<RealDim>,
                                                        mpl::int_<Dim>>::type::value;

    typedef mpl::vector<boost::none_t, boost::none_t, boost::none_t,
                        GeoMap<1, Order, RealDim, T, Entity, PP>>
        geomap_edges_t;

    typedef mpl::vector<boost::none_t, boost::none_t,
                        GeoMap<1, Order, RealDim, T, Entity, PP>,
                        GeoMap<2, Order, nRealDimCheck2d /*RealDim*/, T, Entity, PP>>
        geomap_faces_t;

    typedef mpl::vector<GeoMap<1, Order, RealDim, T, Entity, PP>,
                        GeoMap<2, Order, nRealDimCheck2d /*RealDim*/, T, Entity, PP>,
                        GeoMap<3, Order, 3, T, Entity, PP>,
                        boost::none_t>
        geomap_elements_t;

  public:
    static constexpr uint16_type nDim = super::nDim;
    static constexpr uint16_type nRealDim = super::nRealDim;
    static constexpr uint16_type nDof = super::nDof;
    static constexpr uint16_type nOrder = Order;
    static constexpr uint16_type nNodes = super::nNodes;
    static constexpr fem::transformation_type trans = super::trans;
    static constexpr bool is_linear = ( trans == fem::LINEAR );
    
    typedef typename super::value_type value_type;

    typedef typename super::PreCompute precompute_type;
    typedef std::shared_ptr<precompute_type> precompute_ptrtype;

    typedef typename super::convex_type convex_type;
    typedef Reference<convex_type, nDim, Order, nDim /*nRealDim*/> reference_convex_type;

    typedef GeoMap<Dim, Order, RealDim, T, Entity, PP> self_type;
    typedef self_type geometric_mapping_type;
    typedef std::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

    typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim>>::type element_gm_type;
    typedef std::shared_ptr<element_gm_type> element_gm_ptrtype;
    typedef element_gm_ptrtype gm_ptrtype;
    using gm_type = element_gm_type;

    typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim>>::type face_gm_type;
    typedef std::shared_ptr<face_gm_type> face_gm_ptrtype;

    template <int N>
    struct face_gm
    {
        typedef typename mpl::at<geomap_faces_t, mpl::int_<N>>::type type;
    };

    typedef typename mpl::at<geomap_edges_t, mpl::int_<nDim>>::type edge_gm_type;
    typedef std::shared_ptr<edge_gm_type> edge_gm_ptrtype;

    template <int N>
    struct edge_gm
    {
        typedef typename mpl::at<geomap_edges_t, mpl::int_<N>>::type type;
    };

    typedef typename node<value_type>::type node_t_type;
    typedef typename matrix_node<value_type>::type matrix_node_t_type;

    typedef node_t_type normal_type;
    typedef ublas::vector<normal_type> normals_type;
    typedef typename normals_type::const_iterator normal_const_iterator;
    typedef node_t_type tangent_type;
    typedef ublas::vector<tangent_type> tangents_type;
    typedef typename tangents_type::const_iterator tangent_const_iterator;

    typedef typename ublas::vector<value_type> vector_type;
    typedef typename ublas::matrix<value_type, ublas::column_major> matrix_type;
    
    using eigen_vector_p_type = eigen_vector_type<Dim,value_type>;
    using vector_eigen_vector_p_type = vector_eigen_vector_type<Dim,value_type>;

    using eigen_vector_n_type = eigen_vector_type<RealDim,value_type>;
    using vector_eigen_vector_n_type = vector_eigen_vector_type<RealDim,value_type>;

    using eigen_map_matrix_type = Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>;
    using eigen_map_vector_type = Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,1,Eigen::ColMajor>>;

    using hessian_basis_type = Eigen::Tensor<value_type,3>;
    
    /** default constructor */
    GeoMap()
        : super(),
          M_is_cached( false ),
          _elementMap(),
          _boundaryMap(),
          M_g_linear( nNodes, nDim ),
          M_refconvex()
    {
        if ( trans == fem::LINEAR )
        {
            //M_g_linear.resize( nNodes, nDim );
            matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );

            ublas::vector<ublas::matrix<value_type>> m = super::derivate( __dummy_pts );

            FEELPP_ASSERT( M_g_linear.size2() == m.size() )
            ( M_g_linear.size2() )( m.size() ).error( "invalid dimension" );
            FEELPP_ASSERT( m( 0 ).size2() == 1 )
            ( m( 0 ).size2() ).error( "Invalid number of points" );

            FEELPP_ASSERT( M_g_linear.size1() == m( 0 ).size1() )
            ( M_g_linear.size1() )( m( 0 ).size1() ).error( "invalid number of DOF" );

            //std::cout << "nNodes= " << nNodes << "\n"
            //<< "nDim= " << nDim << "\n";
            //std::cout << "M_g_linear = " << M_g_linear << "\n"
            //<< "m(0) = " << m( 0 ) << "\n";

            for ( uint16_type i = 0; i < nNodes; ++i )
            {
                for ( uint16_type n = 0; n < nDim; ++n )
                {
                    //std::cout << "m(n)= " << m( n ) << "\n";
                    M_g_linear( i, n ) = m( n )( i, 0 );
                }
            }

#if 0

                       for ( uint16_type i = 0; i < nNodes; ++i )
                       {
                           for ( uint16_type n = 0; n < nDim; ++n )
                           {
                               M_g_linear( i, n ) = this->dPhi( i, n, __dummy_pt );
                           }
                       }

#endif // 0
        }
    }
    /** default constructor */
    GeoMap( element_gm_ptrtype const& e, face_gm_ptrtype const& f )
        : super(),
          M_is_cached( false ),
          _elementMap( e ),
          _boundaryMap( f ),
          M_g_linear( nNodes, nDim ),
          M_refconvex()
    {
        if ( trans == fem::LINEAR )
        {
            //M_g_linear.resize( nNodes, nDim );
            node_t_type __dummy_pt( nDim );

            matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );

            //std::cout << "geomap::derivate<> pts=" << __dummy_pts << "\n";
            //std::cout << "geomap::derivate<> m=" << super::derivate( __dummy_pts ) << "\n";
            ublas::vector<ublas::matrix<value_type>> m = super::derivate( __dummy_pts );
            //std::cout << "nNodes= " << nNodes << "\n"
            //<< "nDim= " << nDim << "\n";
            //std::cout << "M_g_linear = " << M_g_linear << "\n"
            //<< "m(0) = " << m( 0 ) << "\n";

            for ( uint16_type i = 0; i < nNodes; ++i )
            {
                for ( uint16_type n = 0; n < nDim; ++n )
                {
                    //std::cout << "m(n)= " << m( n ) << "\n";
                    M_g_linear( i, n ) = m( n )( i, 0 );
                }
            }
        }
    }
    /**
            destructor
            */
    ~GeoMap() override
    {
    }

    /**
            \return the dimension of the underlying element
            */
    constexpr uint16_type dim() const noexcept
    {
        return nDim;
    }

    constexpr uint16_type realDim() const noexcept
    {
        return nRealDim;
    }

    /**
   \return true if the geometric mapping is linear, false otherwise
*/
    bool isLinear() const
    {
        return trans == fem::LINEAR;
    }

    /**
 *  \return the natural mapping for the element
 */
    const element_gm_ptrtype elementMap() const
    {
        return _elementMap;
    }

    /**
 *  \return the natural mapping for the boundary of the element
 */
    const face_gm_ptrtype& boundaryMap() const
    {
        return _boundaryMap;
    }

    /**
 * \return the i-th reference node associated with the geometric
 * mapping construction
 */
    ublas::vector<value_type> refNode( uint16_type i ) const
    {
        return ublas::column( this->points(), i );
    }

    /**
 * @return the reference convex
 */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    /**
 * @return a positive number if the point is within the reference
 * convex
 */
    boost::tuple<bool, value_type> isIn( typename node<value_type>::type const& pt ) const
    {
        return M_refconvex.isIn( pt );
    }

    /**
     *  apply the geometric mapping to the point \c pt given the real
     *  geometric nodes stored in a NxNg matrix \c G
     */
    node_t_type transform(const node_t_type& ref_p,
                          const matrix_node_t_type& G) const
    {
        // Create a real_p node of the appropriate size and clear it
        typename node<value_type>::type real_p(G.size1());
        real_p.clear();

        // Loop over all nodes
        for (uint16_type i = 0; i < nNodes; ++i)
        {
            // Evaluate phi at the current point and assign it to phi_at_pt
            value_type phi_at_pt = super::evaluate(i, ref_p)(0);

            // Perform the transformation by accumulating the product of phi_at_pt and the i-th column of G
            real_p.plus_assign(phi_at_pt * ublas::column(G, i));
        }

        return real_p;
    }

    /**
     * apply the geometric mapping to the point index \c id_pt given the real
     * geometric nodes stored in a NxNg matrix \c G
     */
    node_t_type transform( uint16_type id_ref,
                           const matrix_node_t_type& G,
                           const precompute_type* pc) const
    {
        // Initialize the real_p node with the appropriate size and clear it
        node_t_type real_p(G.size1());
        real_p.clear();

        // Loop over all nodes
        for (uint16_type i = 0; i < nNodes; ++i)
        {
            // Evaluate the transformation at the given point
            value_type phi_at_pt = pc->phi(id_ref, i);

            // Accumulate the result of phi_at_pt multiplied by the i-th column of G
            real_p.plus_assign(phi_at_pt * ublas::column(G, i));
        }

        return real_p;
}

    /**
     * compute real coordinates from a matrix of ref coordinates
     */
    void transform( matrix_node_t_type const& G,
                    precompute_type const* pc,
                    matrix_type& x ) const
    {
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> P( G.data(), G.rows(), G.cols() );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> Phi( pc->phi().data(), pc->phi().rows(), pc->phi().cols() );
        Eigen::Map<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> X( x.data(), x.size() );
        X = P * Phi;
        //ublas::axpy_prod( G, pc->phi(), x, true );
    }

    /**
     *  compute the gradient of the transformation in the reference
     *  element
     *
     *  Compute the gradient at node \c x, pc is resized to
     *  [nbNodes() x dim()] if the transformation is linear, \c x is
     *  not used at all
     */
    void gradient( const node_t_type& __pt,
                   matrix_type& __g ) const
        {
            if ( trans == fem::LINEAR )
            {
                __g = M_g_linear;
            }

            else
            {
                FEELPP_ASSERT( __pt.size() == dim() )
                    ( __pt.size() )( dim() ).error( "invalid dimension" );

                matrix_node_t_type __pts( nDim, 1 );
                ublas::column( __pts, 0 ) = __pt;

                ublas::vector<ublas::matrix<value_type>> m = super::derivate( __pts );

                for ( uint16_type n = 0; n < nDim; ++n )
                {
                    ublas::column( __g, n ) = ublas::column( m( n ), 0 );
                }
            }
        }

    /**
   compute the gradient of the transformation in the reference
   element

   Compute the gradient at node \c x, pc is resized to
   [nbNodes() x dim()] if the transformation is linear, \c x is
   not used at all
*/
    void gradient( uint16_type __idref,
                   matrix_type& __g,
                   precompute_type const* __pc ) const
    {
        FEELPP_ASSERT( __pc )
        ( __idref ).error( "a PreCompute must be set first before using this function" );

        if ( trans == fem::LINEAR )
        {
            __g = M_g_linear;
        }

        else
        {
            FEELPP_ASSERT( __pc->dim() == dim() )
            ( __pc->dim() )( dim() ).error( "invalid dimension" );

            for ( size_type i = 0; i < nNodes; ++i )
            {
                for ( uint16_type n = 0; n < nDim; ++n )
                {
                    __g( i, n ) = __pc->grad( i, 0, n, __idref );
                }
            }
        }
    }

    /**
     compute the gradient of the transformation in the reference
     element
     
     Compute the gradient at node \c x, pc is resized to
     [nbNodes() x dim()] if the transformation is linear, \c x is
     not used at all
     */
    void hessianBasisAtPoint( uint16_type __idref,
                              hessian_basis_type& _hessian,
                              precompute_type const* __pc ) const
        {
            DCHECK( __pc ) << "a PreCompute must be set first before using this function:"  << __idref;
            
            for ( size_type i = 0; i < nNodes; ++i )
            {
                _hessian.chip( i, 0 ) = __pc->hessian( i, __idref ).chip( 0, 2 );
            }
        }

    /**
 * get an estimate of the radius of the element defined by G
 *
 * @param G matrix of nodes defining the element
 * @return the estimate of the radius of the element defined by G
 */
    value_type radiusEstimate( matrix_node_t_type const& G ) const
    {
        size_type N = G.size1();
        size_type Nref = this->nbNodes();
        //size_type n = ( this->isLinear()) ? 1 : Nref;

        // compute K
        matrix_node_t_type __g( Nref, dim() );
        // xref is not defined at this point: it serves only as a dummy point
        // in the case of linear inversion (__g is a constant matrix in this case)
        this->gradient( this->refNode( 0 ), __g );

        FEELPP_ASSERT( __g.size1() == G.size2() )
        ( __g.size1() )( G.size2() ).error( "invalid sizes" );

        matrix_node_t_type K( G.size1(), __g.size2() );
        ublas::axpy_prod( G, __g, K );

        SVD<matrix_node_t_type> __svd( K );
        value_type __max;
        value_type __min;
        __svd.conditionNumber( __max, __min );
        return __max * sqrt( value_type( N ) ) / value_type( N );
    }

    bool M_is_cached;
    bool isCached() const
    {
        return M_is_cached;
    }

    std::unordered_map<uint32_type,value_type> M_J;
    using element_face_pair_t = std::pair<uint32_type,uint16_type>;
    std::unordered_map<element_face_pair_t,value_type,boost::hash<element_face_pair_t>> M_JFace, M_n_norm;
    using allocator_vector_n_t = Eigen::aligned_allocator<std::pair<const element_face_pair_t, eigen_vector_n_type> >;
    std::unordered_map<element_face_pair_t,eigen_vector_n_type,boost::hash<element_face_pair_t>, std::equal_to<element_face_pair_t>, allocator_vector_n_t> M_un_real;
    std::unordered_map<element_face_pair_t,int,boost::hash<element_face_pair_t>> M_permutation_element_face_neighbor;
    unordered_map_eigen_matrix_type<uint32_type,nRealDim,nDim,value_type> M_K;
    unordered_map_eigen_matrix_type<uint32_type,nRealDim,nDim,value_type> M_B;


    template <typename MeshType>
    void initCache( MeshType const* mesh )
    {
#if 0
        size_type nelts = mesh->numElements();
        LOG( INFO ) << "[Geomap] start caching J,K,B for " << nelts << " elements\n";

        M_J.reserve( nelts );
        M_K.clear();
        M_K.reserve( nelts );
        M_B.reserve( nelts );
        M_JFace.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_n_norm.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_un_real.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_permutation_element_face_neighbor.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_is_cached = true;
#endif
    }
    bool isCacheValid() const
    {
        if ( !( M_J.size() > 0 &&
                M_K.size() > 0 &&
                M_B.size() > 0 ) )
        {
            return false;
        }

        return true;
    }
    //! @return true if geomap data are cache for element @p e, false otherwise
    bool cachedK( int e ) const
    {
        return M_K.count( e ) > 0;
    }
    bool cachedB( int e ) const
    {
        return M_B.count( e ) > 0;
    }
    bool cachedJ( int e ) const
    {
        return M_J.count( e ) > 0;
    }
    //! return jacobian at element @p e
    value_type J( int e ) const
    {
        return M_J.at(e);
    }
    void setJacobian( int e, value_type v ) 
        {
            M_J[e]=v;
        }
    bool hasJacobian( int e ) const
        {
            return M_J.count(e);
        }
    value_type jacobianAtFace( int e, int f ) const
        {
            return M_JFace.at({e,f});
        }
    eigen_matrix_type<nDim,nRealDim,value_type> const& B( int e ) const
    {
        return M_B.at(e);
    }
    eigen_matrix_type<nDim,nRealDim,value_type> const& K( int e ) const
    {
        return M_K.at(e);
    }
    void setK( int e, eigen_matrix_type<nDim,nRealDim,value_type> const& K ) 
        {
            M_K[e]=K;
        }
    template<typename G = self_type>
    bool cacheK( int e,
                 eigen_matrix_type<nRealDim,nDim,value_type>& K )
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                if ( !cachedK(e) ) return false;
                K.noalias() = M_K.at(e);
                return true;
            }
            else
                return false;
        }
    template<typename G = self_type>
    bool cacheB( int e,
                 eigen_matrix_type<nRealDim,nDim,value_type>& B )
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                if ( !cachedB(e) ) return false;
                B.noalias() = M_B.at(e);
                return true;
            }
            else
                return false;
        }
    template<typename G = self_type>
    bool cacheJ( int e, value_type& J )
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                if ( !cachedJ(e) ) return false;
                J = M_J.at(e);
                return true;
            }
            else
                return false;
        }
    template<typename G = self_type>
    void updateCacheK( int e,
                       eigen_matrix_type<nRealDim,nDim,value_type> const& K )
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                M_K[e] = K;
            }
        }
    template<typename G = self_type>
    void updateCacheB( int e,
                       eigen_matrix_type<nRealDim,nDim,value_type> const & B)
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                M_B[e] = B;
            }
        }
    template<typename G = self_type>
    void updateCacheJ( int e,
                       value_type const& J )
        {
            if constexpr ( G::is_linear && (nDim==nRealDim) )
            {
                M_J[e] = J;
            }
        }
    void addJ( int e, value_type v )
    {
        DCHECK( e >= 0 && e < M_J.size() ) << "invalid element id " << e << "( " << M_J.size() << ") for geomap cache, are you using the proper mesh\n";
        M_J[e] = v;
    }
    bool hasJacobianAtFace( int e, int f ) const
        {
            return M_JFace.count({e,f});
        }
    void setJacobianAtFace( int e, int f, value_type v )
        {
            M_JFace[{e,f}] = v;
        }
    value_type normalNormAtFace( int e, int f ) const
        {
            return M_n_norm.at({e,f});
        }
    bool hasNormalNormAtFace( int e, int f ) const
        {
            return M_n_norm.count({e,f});
        }
    void setNormalNormAtFace( int e, int f, value_type v )
        {
            M_n_norm[{e,f}] = v;
        }
    eigen_vector_n_type const& unitNormalAtFace( int e, int f ) const
        {
            return M_un_real.at({e,f});
        }
    bool hasUnitNormalAtFace( int e, int f ) const
        {
            return M_un_real.count({e,f});
        }
    void setUnitNormalAtFace( int e, int f, eigen_vector_n_type const& v )
        {
            M_un_real[{e,f}] = v;
        }
    int permutationWithNeighborFace( int e, int f ) const
        {
            return M_permutation_element_face_neighbor.at({e,f});
        }
    bool hasPermutationWithNeighborFace( int e, int f ) const
        {
            return M_permutation_element_face_neighbor.count({e,f});
        }
    void setPermutationWithNeighborFace( int e, int f, int v )
        {
            M_permutation_element_face_neighbor[{e,f}] = v;
        }
    template <typename Derived>
    void setB( int e, const Eigen::MatrixBase<Derived>& B )
    {
        M_B[e].noalias() = B;
    }

    /**
     * \class Context
     *
     * Context for the geometric mapping depend on a node in the
     * reference element
     */
    template <typename ElementType, int SubEntityCoDim = 0>
    class Context
    {
      public:
        static const int subEntityCoDim = SubEntityCoDim;
        static const int subEntityCoDimFix = SubEntityCoDim > 0 ? SubEntityCoDim : 1;
        // reference space dimension
        static inline const uint16_type PDim = ElementType::nDim;
        // real space dimension
        static constexpr uint16_type NDim = ElementType::nRealDim;
        static inline const uint16_type nDim = NDim;
        // type of transformation (linear or not)
        static const fem::transformation_type trans = geometric_mapping_type::trans;
        static const bool is_linear = ( trans == fem::LINEAR );

        static const bool condition = ( ( PDim == NDim ) || ( ( NDim >= 1 ) && ( PDim == NDim - 1 ) ) );
        //BOOST_MPL_ASSERT_MSG( condition, INVALID_DIM, (mpl::int_<NDim>, mpl::int_<PDim>, ElementType ) );
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim>>,
                                  mpl::identity<GeoMap<Dim, Order, NDim, T, Entity, PP>>,
                                  typename mpl::if_<mpl::and_<mpl::greater_equal<mpl::int_<NDim>, mpl::int_<1>>,
                                                              mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim - 1>>>,
                                                    //typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
                                                    // mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
                                                    mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP>::face_gm_type>,
                                                    typename mpl::if_<mpl::and_<mpl::equal_to<mpl::int_<NDim>, mpl::int_<3>>,
                                                                                mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim - 2>>>,
                                                                      //typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
                                                                      // mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
                                                                      mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP>::edge_gm_type>,
                                                                      mpl::identity<boost::none_t>>::type>::type>::type::type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;

        typedef typename gm_type::value_type value_type;

        typedef typename gm_type::precompute_ptrtype precompute_ptrtype;

        typedef Context<ElementType, SubEntityCoDim> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;

        typedef node_t_type normal_type;
        typedef ublas::vector<normal_type> normals_type;
        typedef typename normals_type::const_iterator normal_const_iterator;
        typedef node_t_type tangent_type;
        typedef ublas::vector<tangent_type> tangents_type;
        typedef typename tangents_type::const_iterator tangent_const_iterator;

        typedef ElementType element_type;
        //typedef typename element_type::permutation_type permutation_type;
        typedef typename element_type::template PermutationSubEntity<subEntityCoDimFix> permutation_type;
        using mesh_marker_type = typename element_type::marker_type;

        using eigen_matrix_nx_type = eigen_matrix_type<NDim,Eigen::Dynamic,value_type>;
        using eigen_matrix_xn_type = eigen_matrix_type<Eigen::Dynamic,NDim,value_type>;
        using eigen_matrix_px_type = eigen_matrix_type<PDim,Eigen::Dynamic,value_type>;
        using eigen_matrix_xp_type = eigen_matrix_type<Eigen::Dynamic,PDim,value_type>;

        using eigen_matrix_nn_type = eigen_matrix_type<NDim,NDim,value_type>;
        using vector_eigen_matrix_pp_type = vector_eigen_matrix_type<PDim,PDim,value_type>;
        using vector_eigen_matrix_nn_type = vector_eigen_matrix_type<NDim,NDim,value_type>;
        using eigen_matrix_np_type = eigen_matrix_type<NDim,PDim,value_type>;
        using eigen_matrix_pp_type = eigen_matrix_type<PDim,PDim,value_type>;
        using vector_eigen_matrix_np_type = vector_eigen_matrix_type<NDim,PDim,value_type>;

        using hessian_type = tensor3_fixed_size_t<NDim,PDim,PDim,value_type>;
        using vector_hessian_type = vector_tensor3_fixed_size_t<NDim,PDim,PDim,value_type>;

        static constexpr bool is_on_face = subEntityCoDim == 1;

        Context() = default;

        template <int CTX>
        Context( gm_ptrtype __gm,
                 element_type const& __e,
                 std::integral_constant<int, CTX> const& statctx,
                 precompute_ptrtype const& __pc = precompute_ptrtype(),
                 uint16_type __f = invalid_uint16_type_value,
                 size_type dynctx = 0 )
            : M_gm( __gm ),
              M_element( nullptr ),
              M_pc( __pc ),
              M_pc_faces(),
              M_npoints( ( M_pc ) ? M_pc->nPoints() : 0 ),
              M_G( NDim, gm_type::nNodes ),
              M_g( M_G.cols(), PDim ),
              M_hessian_basis_at_pt( M_G.cols(), PDim, PDim ),
              M_CS(),
              M_CSi(),
              M_B3( boost::extents[NDim][NDim][PDim][PDim] ),
              M_id( __e.id() ),
              M_e_markers( __e.markers() ),
              M_elem_id_1( invalid_v<size_type> ),          // __e.ad_first() ),
              M_pos_in_elem_id_1( invalid_uint16_type_value ), //__e.pos_first() ),
              M_elem_id_2( invalid_v<size_type> ),          //__e.ad_second() ),
              M_pos_in_elem_id_2( invalid_uint16_type_value ), //__e.pos_second() ),
              M_face_id( __f ),
              M_h( 0 ),
              M_h_min( 0 ),
              M_h_face( 0 ),
              M_meas( 0 ),
              M_measface( 0 ),
              M_perm(),
              M_dynamic_context( dynctx ),
              M_context_updated( 0 )
        {
            this->setElement( __e );

            if constexpr ( subEntityCoDim == 0 )
            {
                if ( M_pc )
                {
                    resize<CTX,true>();
                    update<CTX>( __e );
                }
            }
        }

        template <int CTX>
        Context( gm_ptrtype __gm,
                 element_type const& __e,
                 std::integral_constant<int, CTX> const& statctx,
                 std::vector<std::map<permutation_type, precompute_ptrtype> > & __pc,
                 uint16_type __f,
                 size_type dynctx = 0 )
            :
            Context( __gm, __e, statctx, __pc.empty()?precompute_ptrtype{}:__pc[__f][__e.permutation(__f, mpl::int_<subEntityCoDimFix>() )], __f, dynctx )
            {
                if ( !__pc.empty() )
                    M_pc_faces =  __pc;
                if ( M_pc )
                {
                    resize<CTX,true>();
                    update<CTX>( __e, M_face_id );
                }
            }

        Context( gmc_ptrtype& p )
            :
            Context( *p )
            {}

        Context( Context const& ) = default;
        Context( Context && ) = default;

        /**
         * clone the context
         */
        gmc_ptrtype clone()
        {
            return std::make_shared<gmc_type>( *this );
        }

        /**
         * update information on this context
         *
         *  - update the coordinate of the real points
         *  - update the pseudo-inverse of the gradient of the transformation
         *
         *  compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         *
         *  compute \f$ B(x_{\mathrm{ref}}) = K ( K^T K )^{-1} \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         */

        template <size_type CTX, int TheSubEntityCoDim = SubEntityCoDim>
        void update( element_type const& __e,
                     std::enable_if_t< TheSubEntityCoDim == 0 >* = nullptr )
        {
            this->setElement( __e );
            this->updateImpl<CTX>();
        }

        template<size_type CTX,int TheSubEntityCoDim = SubEntityCoDim>
        void update( element_type const& __e, uint16_type __f,
                     std::enable_if_t< TheSubEntityCoDim != 0 >* = nullptr)
            {
                this->setElement( __e );
                this->setSubEntity( __f );
                this->updateImpl<CTX>();
            }

        template<size_type CTX,int TheSubEntityCoDim = SubEntityCoDim>
        void update( element_type const& __e, uint16_type __f, permutation_type __perm,
                     std::enable_if_t< TheSubEntityCoDim != 0 >* = nullptr)
            {
                this->setElement( __e );
                this->setSubEntity( __f,__perm );
                this->updateImpl<CTX>();
            }

        /**
        * update information on this context
        *
        *  - update the coordinate of the real points
        *  - update the pseudo-inverse of the gradient of the transformation
        *
        *  compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
        *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
        *
        *  compute \f$ B(x_{\mathrm{ref}}) = K ( K^T K )^{-1} \f$
        *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
        */

        template<size_type CTX>
        void update( element_type const& __e,
                     precompute_ptrtype const& __pc,
                     uint16_type __f = invalid_uint16_type_value,
                     bool updatePc = true )
        {
            M_pc = __pc;

            if ( M_npoints != M_pc->nPoints() )
            {
                M_npoints = M_pc.get()->nPoints();
                resize<CTX>();
            }

            if constexpr ( subEntityCoDim > 0 )
                update<CTX>( __e, __f/*, updatePc*/ );
            else
                update<CTX>( __e );
        }


        template<size_type CTX>
        void updateContext( size_type dynctx = 0 )
            {
                M_dynamic_context = M_dynamic_context | dynctx;

                if ( M_dynamic_context != 0 )
                {
                    this->resize<CTX|vm::DYNAMIC,true>();
                    this->updateImpl<CTX|vm::DYNAMIC>();
                }
                else
                {
                    this->resize<CTX,true>();
                    this->updateImpl<CTX>();
                }
            }


        FEELPP_STRONG_INLINE
        void setElement( element_type const& __e )
            {
                element_type const* newElement = std::addressof( __e );
                bool isSameElt = newElement == M_element; // ( dynamic_cast<void const*>( newElement ) == dynamic_cast<void const*>( M_element ) )
                if ( !isSameElt )
                {
                    M_context_updated = 0;
                    M_element = std::addressof( __e );
                    if constexpr ( !is_reference_convex_v<element_type> )
                    {
                        if constexpr ( gm_type::nNodes == element_type::numVertices )
                            __e.updateVertices( M_G );
                        else
                            __e.updateG( M_G );
                    }
                    else
                    {
                        if constexpr( std::is_same_v<value_type,typename element_type::value_type> )
                            M_G = em_cmatrix_col_type<value_type>( ( ( gm_type::nNodes == element_type::numVertices ) ? __e.vertices() : __e.G() ).data().begin(), M_G.rows(), M_G.cols() );
                        else
                            M_G = em_cmatrix_col_type<typename element_type::value_type>( ( ( gm_type::nNodes == element_type::numVertices ) ? __e.vertices() : __e.G() ).data().begin(), M_G.rows(), M_G.cols() ).template cast<value_type>();
                    }

                    M_id = __e.id();
                    M_e_markers = __e.markers();
                }

            }

        template<int TheSubEntityCoDim = SubEntityCoDim >
        void setSubEntity( uint16_type __f,
                           std::enable_if_t< TheSubEntityCoDim != 0 >* = nullptr )
            {
                M_face_id = __f;
                CHECK( M_face_id != invalid_v<uint16_type> ) << "sub entity id not defined";
                M_perm = this->element().permutation( M_face_id, mpl::int_<subEntityCoDim>() );
                M_pc = M_pc_faces[__f][M_perm];
                M_f_markers = entityMarkers<subEntityCoDim>( this->element(), __f );
            }
        template<int TheSubEntityCoDim = SubEntityCoDim >
        void setSubEntity( uint16_type __f, permutation_type __perm,
                           std::enable_if_t< TheSubEntityCoDim != 0 >* = nullptr )
            {
                M_face_id = __f;
                CHECK( M_face_id != invalid_v<uint16_type> ) << "sub entity id not defined";
                M_perm = __perm;
                M_pc = M_pc_faces[__f][M_perm];
                M_f_markers = entityMarkers<subEntityCoDim>( this->element(), __f );
            }

        ~Context()
        {
        }

        /** @name Accessors
         */
        //@{

        /**
         \return the geometric mapping associated with the context
         */
        gm_ptrtype const& geometricMapping() const
            {
                return M_gm;
            }

        /**
         \return the dimension of the space of the real element
         */
        uint16_type N() const
            {
                return NDim;
            }

        /**
         \return the dimension of the space of the reference element
         */
        uint16_type P() const
            {
                return PDim;
        }

        /**
         * \return the real element to which this context is associated
         */
        element_type const& element() const
            {
            return *M_element;
            }

        element_type const& element_c() const
            {
                return *M_element;
            }

        //!
        //! @return true if geomap associated to a face, false otherwise
        //!
        bool isOnFace() const
            {
                return (subEntityCoDim == 1) && (M_face_id != invalid_uint16_type_value);
            }
        //!
        //! @return true if geomap associated to a face, edge or point, false otherwise
        //!
        bool isOnSubEntity() const
            {
                return (subEntityCoDim > 0) && (M_face_id != invalid_uint16_type_value);
                //return subEntityCoDim > 1 || ( (subEntityCoDim == 1) && (M_face_id != invalid_uint16_type_value) );
            }

        //!
        //! @return the number of points to transfer from reference to real element
        //!
        uint16_type nPoints() const
            {
                return M_npoints;
            }

        //!
        //! in the case of linear transformation, the jacobian is constant and
        //! we don't have to compute the information at all points, only at one
        //!
        //! @return the number of points at which the geometric mapping is really computed. 
        //!
        constexpr uint16_type nComputedPoints() const
        {
            //return is_linear_polynomial_v<gm_type>?1:M_npoints;
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return 1;
            else
                return M_npoints;
        }

        /**
         * \return the set of points in the reference convex
         */
        matrix_node_t_type const& xRefs() const
        {
            return M_pc->nodes();
        }

        /**
         * \return the q-th point in the reference convex
         */
        ublas::matrix_column<matrix_node_t_type const> xRef( int q ) const
        {
            return M_pc->node( q );
        }

        /**
         * \return the node in the real element
         */
        //matrix_node_t_type const& xReal() const
        matrix_type const& xReal() const
        {
            //BOOST_STATIC_ASSERT( vm::has_point<context>::value );
            return M_xrealq;
        }

        /**
         * \return the node in the real element
         */
        ublas::matrix_column<matrix_type const> xReal( int q ) const
        {
            // BOOST_STATIC_ASSERT( vm::has_point<context>::value );
            return ublas::column( M_xrealq, q );
        }

        /**
         * Get the jacobian of the transformation at the \p q -th point
         *
         * \param q index of the point where the jacobian is requested
         * \return the jacobian of the transformation
         */
        value_type J( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_J[0];
            else
                return M_J[q];
        }
        eigen_vector_type<Eigen::Dynamic,value_type> const& J() const
        {
            return M_J;
        }
        /**
         * \return the matrix associated with the geometric nodes
         */
        //matrix_node_t_type const& G() const { return M_G; }
        //matrix_type const& G() const
        eigen_matrix_nx_type const& G() const
        {
            return M_G;
        }

        /**
         * Get the inverse of the transformation at the \p i -th point
         *
         * \param q the index of point where the pseudo-inverse is requested
         * \return the pseudo inverse of the transformation
         */
        eigen_matrix_np_type const& B( int q  ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_B[0];
            else
                return M_B[q];
        }

        eigen_matrix_np_type const& B( int c1, int c2, int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_B[0](c1,c2);
            else
                return M_B[q]( c1, c2 );
        }
        vector_eigen_matrix_np_type const& B() const
        {
            return M_B;
        }
        eigen_matrix_np_type const& localBasis( int q  ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_local_basis_real[0];
                else
                   return M_local_basis_real[q]; 
            }

        typename eigen_matrix_np_type::ColXpr const& basisN( int q ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_local_basis_real[0].col(0);
                else
                    return M_local_basis_real[q].col(0);
            }


        value_type const& localBasis( int c1, int c2, int q ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_local_basis_real[0]( c1, c2 );
                else
                    return M_local_basis_real[q](c1,c2);
            }

        value_type const& basisN( int c1, int c2, int q ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_local_basis_real[0](c1,0);
                else
                    return M_local_basis_real[q]( c1, 0 );
            }

        eigen_matrix_nn_type const& projectorTangent( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_Ptangent[0];
            else
                return M_Ptangent[q];
        }

        value_type const& projectorTangent( int c1, int c2, int q ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_Ptangent[0]( c1, c2 );
                else
                    return M_Ptangent[q](c1,c2);
            }

        eigen_matrix_np_type const& K( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_K[0];
            else
                return M_K[q];
        }

        value_type const& K( int c1, int c2, int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_K[0]( c1, c2 );
            else
                return M_K[q]( c1, c2 );
        }

        vector_eigen_matrix_np_type const& K() const
            {
                return M_K;
            }

        hessian_type const& hessian( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_hessian[0];
            else
                return M_hessian[q];
        }

        /**
         * the tensor of rank 4 for the transformation of 2nd
         * order derivatives from the reference element to the real
         * element. The tensor has the shape \f$[N,N,P,P]\f$.
         *
         * \return the tensor of rank 4
         */
        boost::multi_array<value_type, 4> const& B3() const
        {
            return M_B3;
        }

        /**
         * \return the barycenter of the reference nodes
         */
        node_t_type barycenterRef() const
        {
            node_t_type __barycenter( M_gm->dim() );
            __barycenter.clear();

            for ( uint16_type __c = 0; __c < M_gm->refNodes().size(); ++__c )
            {
                __barycenter += M_gm->refNode( __c );
            }

            __barycenter /= M_gm->refNodes().size();
            return __barycenter;
        }

        /**
         \return the barycenter of the geometric nodes
         */
        eigen_vector_type<NDim,value_type> barycenterReal() const
            {
                int nPts = M_G.cols();
                return (1./nPts)*M_G*eigen_vector_type<Eigen::Dynamic,value_type>::Ones( nPts );
            }
        /**
         * tell whether the point is on the surface of the convex
         * @return true if the point is on the surface, false otherwise
         */
        bool isOnConvexSurface() const
        {
            if constexpr ( trans == fem::LINEAR )
            {
#if 0
                // x -x0 - K(0)\bar{x}
                return std::abs( ublas::norm_2( xReal() - ublas::column( M_G, 0 ) -
                                                ublas::prod( M_K, xRef() ) ) ) < 1e-10;
#else
                CHECK( false ) << "TODO";
                return false;
#endif
            }
            else
                return false;
        }

        node_t_type const& refNormal( int /*q*/ ) const
        {
            return M_gm->referenceConvex().normal( M_face_id );
        }

        /**
         * get the norm_2 of normal of the real element
         *
         * @return the norm_2 of the normal of the real element
         */
        value_type normalNorm( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_normal_norms[0];
            else
                return M_normal_norms[q];
        }
        /**
         * get the norm_2 of normal in the real element
         *
         * @return the norm_2 of the normal in the real element
         */
        eigen_vector_type<Eigen::Dynamic,value_type> const& normalNorms() const
        {
            return M_normal_norms;
        }
        /**
         * get the normal of the real element
         *
         * @return the normal of the real element
         */
        eigen_vector_n_type const& normal( int q ) const
            {
                if constexpr ( is_linear_polynomial_v<gm_type> )
                    return M_normals[0];
                else
                    return M_normals[q];
            }

        //!
        //! @return the unit normal of the real element
        //!
        vector_eigen_vector_n_type const& unitNormal() const
        {
            return M_unit_normals;
        }

        //!
        //! @return the unit_normal at point @c q
        //!
        eigen_vector_n_type const& unitNormal( int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_unit_normals[0];
            else
                return M_unit_normals[q];
        }

        //!
        //! @return the unit_normal component @c n at point @c q
        //!
        value_type const& unitNormal( int n, int q ) const
        {
            if constexpr ( is_linear_polynomial_v<gm_type> )
                return M_unit_normals[0]( n );
            else
                return M_unit_normals[q]( n );
        }

        //!
        //! @return the scaled tangent to the current face
        //!
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& tangent(std::enable_if_t<is_linear_polynomial_v<GeoMapT>>* = nullptr) const
        {
            DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
            DCHECK( is_linear ) << "Invalid call to unitTangent, the geometric mapping is not linear";
            return M_tangents[0];
        }

        /**
         * @brief scaled tangent at point
         
         * @param q index of the point at which the tangent is evaluated
         * @return the scaled tamgent at the point
         */
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& tangent( int q, std::enable_if_t<is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
        {
            DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
            return M_tangents[0];
        }
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& tangent( int q, std::enable_if_t<!is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
            {
                DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
                return M_tangents[q];
            }
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& unitTangent( int q, std::enable_if_t<is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
            {
                DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
                return M_unit_tangents[0];
            }
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& unitTangent( int q, std::enable_if_t<!is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
            {
                DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
                return M_unit_tangents[q];
            }
        template<typename GeoMapT=gm_type>
        eigen_vector_n_type const& unitTangent( int c1, int c2, int q, std::enable_if_t<is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
            {
                DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
                return M_unit_tangents[0](c1, c2);
            }
        template<typename GeoMapT=gm_type>
        value_type const& unitTangent( int c1, int c2, int q, std::enable_if_t<!is_linear_polynomial_v<GeoMapT>>* = nullptr ) const
            {
                DCHECK( nDim == 2 ) << "Tangent is available only in 2D";
                return M_unit_tangents[q]( c1, c2 );
            }


        /**
         * get the id of the element
         *
         * @return the id of the element
         */
        size_type id() const
        {
                return M_id;
        }

        /*
         * \return the face id
         */
        uint16_type faceId() const
        {
            return M_face_id;
        }

        /**
         * \return true if the element is a face, false otherwise
         */
        bool elementIsAFace() const
        {
            return M_face_id != invalid_uint16_type_value;
        }

        /**
         * get the marker of the element
         *
         * @return the marker of the element
         */
        mesh_marker_type marker( uint16_type k ) const
            {
                auto itFindMarker = M_e_markers.find( k );
                if ( itFindMarker!= M_e_markers.end() )
                    return itFindMarker->second;
                else
                    return mesh_marker_type{};
            }
        /**
         * get the marker of the element
         *
         * @return the marker of the element
         */
        mesh_marker_type marker() const
        {
            auto itFindMarker = M_e_markers.find( 1 );
            if ( itFindMarker!= M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker2 of the element
         *
         * @return the marker2 of the element
         */
        mesh_marker_type marker2() const
        {
            auto itFindMarker = M_e_markers.find( 2 );
            if ( itFindMarker!= M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker3 of the element
         *
         * @return the marker3 of the element
         */
        mesh_marker_type marker3() const
        {
            auto itFindMarker = M_e_markers.find( 3 );
            if ( itFindMarker!= M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker of the face of the element
         *
         * @return the marker of the face of the  element
         */
        mesh_marker_type entityMarker( uint16_type k = 1 ) const
            {
                if ( !isOnSubEntity() || !M_f_markers )
                    return mesh_marker_type{};
                auto itFindMarker = M_f_markers->find( k );
                if ( itFindMarker!= M_f_markers->end() )
                    return itFindMarker->second;
                else
                    return mesh_marker_type{};
            }

        /**
         * get the id of the first element containing the element
         *
         * @return the id of the first element containing the element
         */
        size_type id1() const
        {
            return M_elem_id_1;
        }

        /**
         * get the local id of the element in the first element containing the element
         *
         * @return the local id of element in the first element containing the element
         */
        uint16_type idIn1() const
        {
            return M_pos_in_elem_id_1;
        }

        /**
         * get the id of the second element containing the element
         *
         * @return the id of the second element containing the element
         */
        size_type id2() const
        {
            return M_elem_id_2;
        }

        /**
         * get the local id of the element in the second element containing the element
         *
         * @return the local id of element in the second element containing the element
         */
        uint16_type idIn2() const
        {
            return M_pos_in_elem_id_2;
        }

        /**
         * get an estimate of the radius of the current element
         *
         * @return the radius estimate of the current element
         */
        value_type radiusEstimate() const
            {
#if 0
            return M_gm->radiusEstimate( M_G );
#else
            CHECK( false ) << "TODO";
            return value_type(0);
#endif
        }

        /**
         * @return the max length of the edges of the element
         */
        value_type h() const
        {
            return M_h;
        }
        /**
     * @return the min length of the edges of the element
     */
        value_type hMin() const
        {
            return M_h_min;
        }
        /**
     * Get max length of the edge of the face of the element
     *
     * @return the max length of the edge of the face of the element
     */
        value_type hFace() const
        {
            return M_h_face;
        }

        /*
     * @return the measure of the element
     */
        value_type meas() const
        {
            return M_meas;
        }

        /*
     * @return the measure of the set of elements which share a vertex with \p M_elements including himself
     */
        value_type measurePointElementNeighbors() const
        {
            return M_element->measurePointElementNeighbors();
        }

        /*
     * @return the measure of the (current) face of the element
     */
        value_type measFace() const
        {
            return M_measface;
        }

        /**
     * \return the permutation associated with the face
     */
        permutation_type permutation() const
        {
            return permutation( mpl::bool_<( nDim >= 2 )>() );
        }

        /**
     * \return the precompute type
     */
        precompute_ptrtype const& pc() const
        {
            return M_pc;
        }

        /**
     * \return the precompute type for the faces
     */
        std::vector<std::map<permutation_type, precompute_ptrtype>> const& pcFaces() const
        {
            return M_pc_faces;
        }

        //! @return the dynamic context associated to the geomap
        size_type dynamicContext() const { return M_dynamic_context; }
        
        //@}

        /** @name  Mutators
     */
        //@{

        /**
     set some precomputed data on the reference element
     */
        void setPc( precompute_ptrtype const& __pc )
        {
            M_pc = __pc;
        }

        void setPcFaces( std::vector<std::map<permutation_type, precompute_ptrtype>> const& __pcfaces )
        {
            M_pc_faces = __pcfaces;
        }

        void
        edgeTangent( int edgeId, ublas::vector<value_type>& t, bool scaled = false ) const
        {
#if 0
            auto const& K = this->K( 0 );

            ublas::axpy_prod( K,
                              this->geometricMapping()->referenceConvex().tangent( edgeId ),
                              t,
                              true );
            if ( scaled )
                t *= this->element().hEdge( edgeId ) / ublas::norm_2( t );
            else
                t /= ublas::norm_2( t );
#else
            Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,1>> e_t( t.data().begin(), t.size() );
            e_t=M_K[0] * M_gm->referenceConvex().tangent( edgeId );
            e_t.normalize();
            if (scaled)
                e_t*=this->element().hEdge( edgeId );
#endif
        }

        void
        faceNormal( int faceId, ublas::vector<value_type>& n, bool scaled = false ) const
        {
#if 0
            auto const& K = this->K( 0 );
            auto const& B = this->B( 0 );

            ublas::axpy_prod( B,
                              this->geometricMapping()->referenceConvex().normal( faceId ),
                              n,
                              true );
            if ( scaled )
                n *= this->element().faceMeasure( faceId ) / ublas::norm_2( n );
            else
                n /= ublas::norm_2( n );
#else
            Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,1>> e_n( n.data().begin(), n.size() );
            e_n=M_B[0] * M_gm->referenceConvex().normal( faceId );
            e_n.normalize();
            if (scaled)
                e_n*=this->element().faceMeasure( faceId );

#endif
        }

        /**
         * update the geomap through a neighbor matching face.
         *
         * Matching means that the points coordinates from both geomap are the same
         * up to a permutation. We compute here the permutation that ensures this
         * and update the geomap accordingly
         *
         * @return true if the permutation has been found and the geomap updated,
         * false otherwise.
         */
        template <int CTX,typename EltType, typename NeighborGeoType,int TheSubEntityCoDim = SubEntityCoDim>
        bool updateFromNeighborMatchingFace( EltType const& elt, uint16_type face_in_elt, std::shared_ptr<NeighborGeoType> const& gmc,
                                             std::enable_if_t< TheSubEntityCoDim == 1 >* = nullptr )
            {
                if ( M_gm->hasPermutationWithNeighborFace( elt.id(), face_in_elt ) == false )
                {
                    auto [found_permutation,perm] = this->updateFromMatchingNodes<CTX>(elt, face_in_elt, gmc );
#if 0 // NOT WORK IN ALL CASES
                    if ( found_permutation )
                        M_gm->setPermutationWithNeighborFace( elt.id(), face_in_elt, perm.value() );
#endif
                    return found_permutation;
                }
                this->update<CTX>( elt, face_in_elt, permutation_type(M_gm->permutationWithNeighborFace(elt.id(), face_in_elt)) );
                return true;
            }

        template <int CTX,typename EltType, typename MatchingGmcType,int TheSubEntityCoDim = SubEntityCoDim>
        std::pair<bool,permutation_type> updateFromMatchingNodes( EltType const& elt, uint16_type subentity_in_elt, std::shared_ptr<MatchingGmcType> const& gmc,
                                                                  std::enable_if_t< TheSubEntityCoDim != 0 >* = nullptr )
            {
                this->setElement( elt );
                auto gmcptr = gmc.get();
                for ( permutation_type __p( permutation_type::IDENTITY );
                      __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    // update only xReal in current geomap
                    this->setSubEntity( subentity_in_elt,__p );
                    this->updatePoints();

                    bool check = true;
                    for ( uint16_type i = 0; i < gmc->nPoints() && check; ++i )
                    {
                        for ( uint16_type d = 0; (d < NDim); ++d )
                        {
                            check = check && ( std::abs( gmcptr->xReal( i )[d] - this->xReal( i )[d] ) < 1e-8 );
                        }
                    }
                    // if check update full gmc context with the good permutation
                    if ( check )
                    {
                        this->update<CTX>( elt, subentity_in_elt, __p );
                        return std::make_pair(true, __p);
                    }
                }
                return std::make_pair(false,permutation_type( permutation_type::IDENTITY ));
            }

        //@}
      private:
        /**
         * \return the permutation associated with the face
         */
        permutation_type permutation( mpl::bool_<false> ) const
        {
            return M_perm;
        }

        /**
         * \return the permutation associated with the face
         */
        permutation_type permutation( mpl::bool_<true> ) const
        {
            FEELPP_ASSERT( M_face_id == invalid_uint16_type_value ||
                           ( M_face_id != invalid_uint16_type_value &&
                             M_perm != permutation_type( permutation_type::NO_PERMUTATION ) ) )
            ( M_face_id ).error( "invalid permutation" );
            return M_perm;
        }

        template<size_type CTX>
        void updateImpl()
            {
                DCHECK( M_G.cols() == M_gm->nbPoints()  ) << "Invalid number of points got " << M_G.cols() << " expected: " << M_gm->nbPoints();
                DCHECK( M_pc ) << "invalid precompute data structure";

                if constexpr ( vm::has_point_v<CTX> )
                {
                    this->updatePoints();
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasPOINT( M_dynamic_context ) )
                        this->updatePoints();
                }
                if constexpr ( vm::has_measure_v<CTX> )
                {
                    this->updateMeasures();
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasMEASURE( M_dynamic_context ) )
                        this->updateMeasures();
                }

                if constexpr ( subEntityCoDim == 1 )//this->isOnFace()
                {
                    if constexpr ( vm::has_measure_v<CTX> || vm::has_tangent_v<CTX> )
                        this->updateFaceMeasures();
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasMEASURE( M_dynamic_context ) || hasTANGENT( M_dynamic_context ) )
                            this->updateFaceMeasures();
                    }
                }

                updateJKBN<CTX>();
                M_context_updated = CTX | M_dynamic_context;
            }

        FEELPP_STRONG_INLINE
        void updateGradient( int q ) noexcept
            {
                M_gm->gradient( q, M_g, M_pc.get() );
            }

        FEELPP_STRONG_INLINE
        void updateMeasures() noexcept
            {
                M_h = M_element->h();
                M_h_min = M_element->hMin();
                M_meas = M_element->measure();
            }

        FEELPP_STRONG_INLINE
        void updateFaceMeasures() noexcept
            {
                M_measface = M_element->faceMeasure( M_face_id );
                M_h_face = M_element->hFace( M_face_id );
            }


        FEELPP_STRONG_INLINE
        void updatePoints() noexcept
            {
                em_matrix_col_type<value_type> Xreal( M_xrealq.data().begin(), M_xrealq.size1(), M_xrealq.size2() );
                //em_matrix_col_type<value_type> Pts( M_G.data().begin(), M_G.size1(), M_G.size2() );
                Xreal.noalias() = M_G * M_pc->phiEigen();
            }

        template<typename ConvexType = ElementType>
        void updateJacobian( eigen_matrix_np_type const & K, value_type& J,
                             std::enable_if_t<dimension_v<ConvexType> == real_dimension_v<ConvexType>>* = nullptr )
        {
            J = math::abs( K.determinant() );
        }
        template<typename ConvexType = ElementType>
        void updateJacobian( eigen_matrix_np_type const& K, value_type& J,
                             std::enable_if_t<dimension_v<ConvexType> != real_dimension_v<ConvexType>>* = nullptr )
        {
            // CS = K^T K
            //M_CSi.noalias() = K.transpose()*K;
            J = math::sqrt( math::abs( M_CSi.determinant() ) );
        }
        template<typename ConvexType = ElementType>
        void updateB( eigen_matrix_np_type const & K, eigen_matrix_np_type& B,
                      std::enable_if_t<dimension_v<ConvexType> == real_dimension_v<ConvexType>>* = nullptr )
        {
            M_CS.noalias() = K.inverse();
            B.noalias() = M_CS.transpose();
        }
        template<typename ConvexType = ElementType>
        void updateB( eigen_matrix_np_type const& K, eigen_matrix_np_type& B,
                             std::enable_if_t<dimension_v<ConvexType> != real_dimension_v<ConvexType>>* = nullptr )
        {
            // CS = K^T K
            // M_CSi.noalias() = K.transpose()*K;
            M_CS=M_CSi.inverse();
            // B = K CS
            B.noalias() = K*M_CS;
            //std::cout << "2.B=" << B << std::endl;
        }
        //!
        //! update hessian information
        //!
        template<int CTX>
        void updateHessian( eigen_matrix_np_type const& B )
        {
            if constexpr ( vm::has_hessian_v<CTX> || vm::has_laplacian_v<CTX> )
            {
#if 0
                for ( uint16_type k = 0; k < NDim; ++k )
                    for ( uint16_type l = 0; l < NDim; ++l )
                        for ( uint16_type i = 0; i < PDim; ++i )
                            for ( uint16_type j = 0; j < PDim; ++j )
                                M_B3[k][l][i][j] = B( k, i ) * B( l, j );
#endif
            }
        }
        //!
        //! update normal in the real element
        //!

        FEELPP_STRONG_INLINE
        void updateNormals( eigen_matrix_np_type const& K, eigen_matrix_np_type const& B,
                            eigen_vector_n_type& N, eigen_vector_n_type& unitN, value_type& Nnorm ) noexcept
        {
            N.noalias() = B * M_gm->referenceConvex().normal( M_face_id );
            #if 0
            if constexpr ( PDim == NDim-1 && NDim == 3)
            {
                N.noalias() = K.col(0).cross(K.col(1));
            }
            else
            {
                N.noalias() = B * M_gm->referenceConvex().normal( M_face_id );
            }
            #endif
            Nnorm = N.norm();
            unitN = N/Nnorm;

            //M_gm->setNormalNormAtFace( M_id, M_face_id, M_n_norm );
            //M_gm->setUnitNormalAtFace( M_id, M_face_id, M_un_real );

        }

        FEELPP_STRONG_INLINE
        void updateTangents( eigen_matrix_np_type const& K, eigen_vector_n_type const& unitN,
                             eigen_vector_n_type& Ta, eigen_vector_n_type& unitT, value_type& Tnorm,
                             eigen_matrix_nn_type& P ) noexcept
        {
            if constexpr ( NDim == 2 )
            {
                // t = |\hat{e}|*o_K*(K*t_ref)/|e| where o_K is the sign(e*x_K(\hat{e}))
                Ta.noalias() = K * M_gm->referenceConvex().tangent( M_face_id );
                Tnorm = Ta.norm();
                unitT = Ta/Tnorm;
            }

            // compute projector on tangent plane
            P.setIdentity();
            P.noalias() -= unitN * unitN.transpose();
        }
        //!
        //! update normal in the real element
        //!
        FEELPP_STRONG_INLINE
        void updateLocalBasis( eigen_matrix_np_type const& B,
                               eigen_matrix_pp_type& local_basis_ref, eigen_matrix_np_type& local_basis_real ) noexcept
        {
#if 0
            if  ( !this->isOnFace() )
                return;//throw std::logic_error("Local basis defined only on faces ");
#endif
            local_basis_ref = eigen_matrix_pp_type::Identity();
            eigen_vector_p_type Np = M_gm->referenceConvex().normal( M_face_id );
            int max_col;
            Np.array().abs().maxCoeff( &max_col );
            if ( max_col != 0 )
                local_basis_ref.col( max_col ) = local_basis_ref.col( 0 );
            local_basis_ref.col(0)=Np;
            local_basis_real = B*local_basis_ref;
            // orthogonalize columns using the Gram-Schmidt algorithm
            for (int col = 0; col < PDim; ++col)
            {
                typename eigen_matrix_np_type::ColXpr colVec = local_basis_real.col(col);
                for (int prevCol = 0; prevCol < col; ++prevCol)
                {
                    typename eigen_matrix_np_type::ColXpr prevColVec = local_basis_real.col(prevCol);
                    colVec -= colVec.dot(prevColVec)*prevColVec;
                }
                local_basis_real.col(col) = colVec.normalized();
            }

            // Ensure basis_real is direct
            if ( (NDim == PDim) && (NDim>1) && ( local_basis_real.determinant() < 0 ) )
            {
                local_basis_real.col(1) *= -1;
            }
        }


        template<size_type CTX,bool ForceResizeJKBN=false>
        void resize()
        {
            if constexpr ( vm::has_point_v<CTX> )
            {
                M_xrealq.resize( NDim, nPoints() );
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasPOINT( M_dynamic_context ) )
                    M_xrealq.resize( NDim, nPoints() );
            }
            if constexpr ( !is_linear_polynomial_v<gm_type> || ForceResizeJKBN )
                 resizeJKBN<CTX>();
        }
        template<size_type CTX>
        void resizeJKBN()
            {
                static const bool IsOnFace = subEntityCoDim == 1;
                if constexpr ( !vm::has_dynamic_v<CTX> )
                {
                    if constexpr ( !IsOnFace )
                    {
                        if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> ) )
                            return;
                    }
                    else
                    {
                        if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> ) )
                            return;
                    }
                }

                static const bool hasStaticB = IsOnFace?
                    vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> :
                    vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX>;
                bool hasJ = vm::has_jacobian_v<CTX>;
                bool hasLocalBasis = IsOnFace && vm::has_local_basis_v<CTX>;
                bool hasTangent = IsOnFace && vm::has_tangent_v<CTX>;
                bool hasNormal = IsOnFace && (vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>);
                bool hasB = hasStaticB;//vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>;
                bool hasSecondDerivative = vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX>;
                if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if constexpr ( !vm::has_jacobian_v<CTX> )
                        hasJ = hasJACOBIAN( M_dynamic_context );
                    if constexpr ( !vm::has_second_derivative_v<CTX> && !vm::has_hessian_v<CTX> )
                        hasSecondDerivative = hasSECOND_DERIVATIVE( M_dynamic_context ) || hasHESSIAN( M_dynamic_context );
                    if constexpr ( IsOnFace )
                    {
                        hasLocalBasis = hasLOCAL_BASIS( M_dynamic_context );
                        hasTangent = hasTANGENT( M_dynamic_context );
                        hasNormal = hasTangent || hasNORMAL( M_dynamic_context );
                        if constexpr ( !hasStaticB )
                             hasB = hasNormal || hasKB( M_dynamic_context );
                    }
                    else
                    {
                        if constexpr ( !hasStaticB )
                             hasB =  hasKB( M_dynamic_context );
                    }

                    if constexpr ( !(vm::has_jacobian_v<CTX> || hasStaticB) )
                           if ( !(hasJ || hasB) )
                               return;
                }

                const uint16_type _nComputedPoints = nComputedPoints();
                // K
                M_K.resize( _nComputedPoints );

                // jacobian
                if constexpr ( vm::has_jacobian_v<CTX> )
                {
                    M_J.resize( _nComputedPoints );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasJ )
                        M_J.resize( _nComputedPoints );
                }

                // B
                if constexpr ( hasStaticB )
                {
                    M_B.resize( _nComputedPoints );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasB )
                        M_B.resize( _nComputedPoints );
                }

                // second derivative
                if constexpr ( !gmc_type::is_linear )
                {
                    if constexpr ( vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> /*|| vm::has_laplacian_v<CTX>*/ )
                    {
                        M_hessian.resize( _nComputedPoints );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasSecondDerivative )
                            M_hessian.resize( _nComputedPoints );
                    }
                }

                if constexpr ( IsOnFace )
                {
                    // normal
                    if constexpr ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> )
                    {
                        M_normals.resize( _nComputedPoints );
                        M_unit_normals.resize( _nComputedPoints );
                        M_normal_norms.resize( _nComputedPoints );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasNormal )
                        {
                            M_normals.resize( _nComputedPoints );
                            M_unit_normals.resize( _nComputedPoints );
                            M_normal_norms.resize( _nComputedPoints );
                        }
                    }
                    // tangent
                    if constexpr ( vm::has_tangent_v<CTX> )
                    {
                        M_tangents.resize( _nComputedPoints );
                        M_unit_tangents.resize( _nComputedPoints );
                        M_tangent_norms.resize( _nComputedPoints );
                        M_Ptangent.resize( _nComputedPoints );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasTangent )
                        {
                            M_tangents.resize( _nComputedPoints );
                            M_unit_tangents.resize( _nComputedPoints );
                            M_tangent_norms.resize( _nComputedPoints );
                            M_Ptangent.resize( _nComputedPoints );
                        }
                    }
                    // local basis
                    if constexpr ( vm::has_local_basis_v<CTX> )
                    {
                        M_local_basis_ref.resize( _nComputedPoints );
                        M_local_basis_real.resize( _nComputedPoints );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasLocalBasis )
                        {
                            M_local_basis_ref.resize( _nComputedPoints );
                            M_local_basis_real.resize( _nComputedPoints );
                        }
                    }
                }

            }



        template<int CTX>
        void updateJKBN() noexcept
            {
                static const bool IsOnFace = subEntityCoDim == 1;
                if constexpr ( !vm::has_dynamic_v<CTX> )
                {
                    if constexpr ( !IsOnFace )
                    {
                        if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX>) )
                            return;
                    }
                    else
                    {
                        if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> ) )
                            return;
                    }
                }

                static const bool hasStaticB = IsOnFace?
                    vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> :
                    vm::has_kb_v<CTX> || vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX>;
                bool hasJ = vm::has_jacobian_v<CTX>;
                bool hasLocalBasis = IsOnFace && vm::has_local_basis_v<CTX>;
                bool hasTangent = IsOnFace && vm::has_tangent_v<CTX>;
                bool hasNormal = IsOnFace && (vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>);
                bool hasB = hasStaticB;//vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>;
                bool hasSecondDerivative = vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX>;
                if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if constexpr ( !vm::has_jacobian_v<CTX> )
                        hasJ = hasJACOBIAN( M_dynamic_context );
                    if constexpr ( !vm::has_second_derivative_v<CTX> && !vm::has_hessian_v<CTX> )
                        hasSecondDerivative = hasSECOND_DERIVATIVE( M_dynamic_context ) || hasHESSIAN( M_dynamic_context );
                    if constexpr ( IsOnFace )
                    {
                        hasLocalBasis = hasLOCAL_BASIS( M_dynamic_context );
                        hasTangent = hasTANGENT( M_dynamic_context );
                        hasNormal = hasTangent || hasNORMAL( M_dynamic_context );
                        if constexpr ( !hasStaticB )
                             hasB = hasNormal || hasKB( M_dynamic_context );
                    }
                    else
                    {
                        if constexpr ( !hasStaticB )
                             hasB =  hasKB( M_dynamic_context );
                    }

                    if constexpr ( !(vm::has_jacobian_v<CTX> || hasStaticB) )
                           if ( !(hasJ || hasB) )
                               return;
                }

                //em_matrix_col_type<value_type> Pts( M_G.data().begin(), M_G.size1(), M_G.size2() );

                //bool _gradient_needs_update = resizeGradient<CTX>();
                em_matrix_col_type<value_type> GradPhi( M_g.data().begin(),
                                                        M_G.cols(), PDim );

                tensor_map_t<2,value_type> TPts( M_G.data()/*.begin()*/, M_G.rows(), M_G.cols() );
                Eigen::array<dimpair_t, 1> dims = {{dimpair_t(1, 0)}};

                //if constexpr ( is_linear_polynomial_v<gm_type> )
                const uint16_type nPts = nComputedPoints();
                for ( uint16_type q = 0; q < nPts; ++q )
                {
                    //uint16_type q = 0;
                    updateGradient( q );

                    auto & K = M_K[q];

                    K.noalias() = /*Pts*/ M_G* GradPhi;

                    if constexpr ( dimension_v<ElementType> != real_dimension_v<ElementType> )
                    {
                        // CS = K^T K
                        M_CSi.noalias() = K.transpose()*K;
                    }

                    // jacobian
                    if constexpr ( vm::has_jacobian_v<CTX> )
                    {
                        updateJacobian( K, M_J[q] );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasJ )
                            updateJacobian( K, M_J[q] );
                    }

                    // B
                    if constexpr ( hasStaticB )
                    {
                        updateB( K, M_B[q] );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasB )
                            updateB( K, M_B[q] );
                    }

                    // second derivative
                    if constexpr ( !gmc_type::is_linear )
                    {
                        if constexpr ( vm::has_second_derivative_v<CTX> || vm::has_hessian_v<CTX> /*|| vm::has_laplacian_v<CTX>*/ )
                        {
                            //M_h.resize( { NDim, NDim } );
                            M_gm->hessianBasisAtPoint( q, M_hessian_basis_at_pt, M_pc.get() );
                            M_hessian[q] = TPts.contract(M_hessian_basis_at_pt,dims);
                            //std::cout << "M_hessian[" << q << "]=" << M_hessian[q] << std::endl;
                        }
                        else if constexpr ( vm::has_dynamic_v<CTX> )
                        {
                            if ( hasSecondDerivative )
                            {
                                M_gm->hessianBasisAtPoint( q, M_hessian_basis_at_pt, M_pc.get() );
                                M_hessian[q] = TPts.contract(M_hessian_basis_at_pt,dims);
                            }
                        }
                    }


                    if constexpr ( IsOnFace )
                    {
                        // normal
                        if constexpr ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> )
                        {
                            updateNormals( M_K[q], M_B[q], M_normals[q], M_unit_normals[q], M_normal_norms[q] );
                        }
                        else if constexpr ( vm::has_dynamic_v<CTX> )
                        {
                            if ( hasNormal )
                                updateNormals( M_K[q], M_B[q], M_normals[q], M_unit_normals[q], M_normal_norms[q] );
                        }
                        // tangent
                        if constexpr ( vm::has_tangent_v<CTX> )
                        {
                            updateTangents( K, M_unit_normals[q], M_tangents[q], M_unit_tangents[q], M_tangent_norms[q], M_Ptangent[q] );
                        }
                        else if constexpr ( vm::has_dynamic_v<CTX> )
                        {
                            if ( hasTangent )
                                updateTangents( K, M_unit_normals[q], M_tangents[q], M_unit_tangents[q], M_tangent_norms[q], M_Ptangent[q] );
                        }
                        // local basis
                        if constexpr ( vm::has_local_basis_v<CTX> )
                        {
                            updateLocalBasis( M_B[q], M_local_basis_ref[q], M_local_basis_real[q] );
                        }
                        else if constexpr ( vm::has_dynamic_v<CTX> )
                        {
                            if ( hasLocalBasis )
                                updateLocalBasis( M_B[q], M_local_basis_ref[q], M_local_basis_real[q] );
                        }
                    }
                }
            }

        friend class boost::serialization::access;

        template <class Archive>
        void serialize( Archive& ar, const unsigned int version )
        {
            if constexpr ( subEntityCoDim == 0 )
            {
                ar& BOOST_SERIALIZATION_NVP( M_dynamic_context );
                ar& BOOST_SERIALIZATION_NVP( M_context_updated );
                if constexpr ( Archive::is_saving::value )
                {
                    ar & boost::serialization::make_nvp( "xRefs", this->xRefs() );
                }
                else if constexpr ( Archive::is_loading::value )
                {
                    size_type saveDynamicContext = M_dynamic_context;
                    matrix_node_t_type _xRefs;
                    ar & boost::serialization::make_nvp( "xRefs", _xRefs );
                    M_pc = M_gm->preCompute( _xRefs );
                    M_dynamic_context = M_context_updated;
                    M_context_updated = 0;
                    M_npoints = M_pc->nPoints();
                    resize<vm::DYNAMIC,true>();
                    update<vm::DYNAMIC>( this->element()/*, __f, updatePc*/ );
                    M_dynamic_context = saveDynamicContext;
                }
            }
            else // TODO apply same procedure as above
            {

            ar& BOOST_SERIALIZATION_NVP( M_npoints );
            ar& BOOST_SERIALIZATION_NVP( M_J );
            ar& BOOST_SERIALIZATION_NVP( M_G );
            //ar& BOOST_SERIALIZATION_NVP( M_ref_normals );
            ar& BOOST_SERIALIZATION_NVP( M_normals );
            ar& BOOST_SERIALIZATION_NVP( M_unit_normals );
            ar& BOOST_SERIALIZATION_NVP( M_normal_norms );
            ar& BOOST_SERIALIZATION_NVP( M_tangents );
            ar& BOOST_SERIALIZATION_NVP( M_unit_tangents );
            ar& BOOST_SERIALIZATION_NVP( M_tangent_norms );
            ar& BOOST_SERIALIZATION_NVP( M_xrealq );
            //ar& BOOST_SERIALIZATION_NVP( M_g_linear );
            ar& BOOST_SERIALIZATION_NVP( M_g );
            ar& BOOST_SERIALIZATION_NVP( M_K );
            ar& BOOST_SERIALIZATION_NVP( M_CS );
            ar& BOOST_SERIALIZATION_NVP( M_CSi );
            ar& BOOST_SERIALIZATION_NVP( M_B );
            ar& BOOST_SERIALIZATION_NVP( M_Ptangent );
            ar& BOOST_SERIALIZATION_NVP( M_B3 ); //
            ar& BOOST_SERIALIZATION_NVP( M_id );
            ar& BOOST_SERIALIZATION_NVP( M_e_markers );
            ar& BOOST_SERIALIZATION_NVP( M_elem_id_1 );
            ar& BOOST_SERIALIZATION_NVP( M_pos_in_elem_id_1 );
            ar& BOOST_SERIALIZATION_NVP( M_elem_id_2 );
            ar& BOOST_SERIALIZATION_NVP( M_pos_in_elem_id_2 );
            ar& BOOST_SERIALIZATION_NVP( M_face_id );
            ar& BOOST_SERIALIZATION_NVP( M_h );
            ar& BOOST_SERIALIZATION_NVP( M_h_min );
            ar& BOOST_SERIALIZATION_NVP( M_h_face );
            ar& BOOST_SERIALIZATION_NVP( M_meas );
            ar& BOOST_SERIALIZATION_NVP( M_measface );
            ar& BOOST_SERIALIZATION_NVP( M_perm );
            }
        }

      private:
        gm_ptrtype M_gm;

        element_type const* M_element;
        //std::shared_ptr<element_type const> M_element_c;
        //element_type M_element_c;

        precompute_ptrtype M_pc;
        std::vector<std::map<permutation_type, precompute_ptrtype>> M_pc_faces;
        uint16_type M_npoints;

        eigen_vector_type<Eigen::Dynamic,value_type> M_J;

        //matrix_type M_G;
        eigen_matrix_nx_type M_G;
        //vector_eigen_vector_p_type M_ref_normals;
        vector_eigen_vector_n_type M_normals, M_unit_normals;
        eigen_vector_type<Eigen::Dynamic,value_type> M_normal_norms;
        vector_eigen_vector_n_type M_tangents, M_unit_tangents;
        std::vector<value_type> M_tangent_norms;
        vector_eigen_matrix_np_type M_local_basis_real;
        vector_eigen_matrix_pp_type M_local_basis_ref;

        matrix_type M_xrealq;

        //matrix_type M_g_linear;
        matrix_type M_g;
        hessian_basis_type M_hessian_basis_at_pt;
        
        vector_eigen_matrix_np_type M_K;
        eigen_matrix_pp_type M_CS;
        eigen_matrix_pp_type M_CSi;
        vector_eigen_matrix_np_type M_B;
        vector_hessian_type M_hessian;
        vector_eigen_matrix_nn_type M_Ptangent;
        
        boost::multi_array<value_type, 4> M_B3;

        size_type M_id;
        std::map<uint16_type,mesh_marker_type> M_e_markers;
        std::optional<std::map<uint16_type,mesh_marker_type>> M_f_markers;
        size_type M_elem_id_1;
        uint16_type M_pos_in_elem_id_1;
        size_type M_elem_id_2;
        uint16_type M_pos_in_elem_id_2;

        uint16_type M_face_id;

        value_type M_h;
        value_type M_h_min;
        value_type M_h_face;
        value_type M_meas;
        value_type M_measface;

        permutation_type M_perm;

        size_type M_dynamic_context;
        size_type M_context_updated;

    }; // Context

    template <size_type Context_v, int SubEntityCoDim_v = 0, typename ElementType_t>
        std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( geometric_mapping_ptrtype gm, ElementType_t const& e, precompute_ptrtype const& pc, size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>( gm,
                                                                           e, std::integral_constant<int, Context_v>{},
                                                                           pc, invalid_uint16_type_value, dynctx );
    }

    template <size_type Context_v, int SubEntityCoDim_v = 0, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t,SubEntityCoDim_v>>
    context( ElementType_t const& e, precompute_ptrtype const& pc, size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t,SubEntityCoDim_v>>(
            //super_enable_this::shared_from_this(),
            std::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP>>( this->shared_from_this() ),
            e,  std::integral_constant<int, Context_v>{},
            pc, invalid_uint16_type_value, dynctx );
    }

    template <size_type Context_v, int SubEntityCoDim_v = 1, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( geometric_mapping_ptrtype gm,
             ElementType_t const& e,
             std::vector<std::map<typename ElementType_t::permutation_type, precompute_ptrtype>>& pc,
             uint16_type f,
             size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>( gm,
                                                                           e, std::integral_constant<int, Context_v>{},
                                                                           pc,
                                                                           f, dynctx );
    }
    template <size_type Context_v, int SubEntityCoDim_v = 1, typename ElementType_t>
        std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( ElementType_t const& e,
             std::vector<std::map<typename ElementType_t::template PermutationSubEntity<SubEntityCoDim_v>, precompute_ptrtype>>& pc,
             uint16_type f, size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>(
                //super_enable_this::shared_from_this(),
                std::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP>>( this->shared_from_this() ),
                e, std::integral_constant<int, Context_v>{},
                pc,
                f, dynctx );
    }

    /**
 * Inverse of the geometric mapping for a given context
 */
    class Inverse
    {
      public:
        typedef GeoMap geometric_mapping_type;
        typedef std::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

        static const fem::transformation_type trans = geometric_mapping_type::trans;

        typedef typename geometric_mapping_type::node_t_type node_t_type;
        typedef typename geometric_mapping_type::node_t_type node_type;
        typedef typename geometric_mapping_type::matrix_node_t_type matrix_node_t_type;
        typedef typename geometric_mapping_type::matrix_node_t_type matrix_node_type;

        typedef ublas::vector<double> dense_vector_type;
        typedef ublas::matrix<double> dense_matrix_type;

        template <typename GeoElem>
        Inverse( geometric_mapping_ptrtype __gm, GeoElem const& __ge,
                 worldcomm_ptr_t const& worldComm = Environment::worldCommSeqPtr() )
            : M_gm( __gm ),
              M_xref( __gm->dim() ),
              M_xreal( __ge.G().size1() ),
              M_is_in( false ),
              M_G( __ge.G() ),
              M_K( N(), __gm->dim() ),
              M_B( N(), __gm->dim() ),
              M_CS( __gm->dim(), __gm->dim() ),
              M_g( M_gm->nbPoints(), __gm->dim() ),
#if defined( FEELPP_HAS_PETSC )
              //M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC, worldComm ) )
              M_nlsolver( SolverNonLinear<double>::build( "petsc", "", worldComm ) )
#else
              M_nlsolver( SolverNonLinear<double>::build( "eigen_dense", "", worldComm ) )
#endif
        {
            this->init();
        }

        template <typename GeoElem>
        Inverse( geometric_mapping_ptrtype __gm, GeoElem const& __ge, mpl::int_<1> /**/,
                 worldcomm_ptr_t const& worldComm = Environment::worldCommSeqPtr() )
            : M_gm( __gm ),
              M_xref( __gm->dim() ),
              M_xreal( __ge.vertices().size1() ),
              M_is_in( false ),
              M_G( __ge.vertices() ),
              M_K( N(), __gm->dim() ),
              M_B( N(), __gm->dim() ),
              M_CS( __gm->dim(), __gm->dim() ),
              M_g( M_gm->nbPoints(), __gm->dim() ),
#if defined( FEELPP_HAS_PETSC )
              //M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC,worldComm) )
              M_nlsolver( SolverNonLinear<double>::build( "petsc", "", worldComm ) )
#else
              M_nlsolver( SolverNonLinear<double>::build( "eigen_dense", "", worldComm ) )
#endif
        {
            this->init();
        }

      private:
        void init()
        {
            FEELPP_ASSERT( M_G.size2() == M_gm->nbPoints() )
            ( M_G.size2() )( M_gm->nbPoints() ).error( "invalid dimensions" );

            if ( M_gm->isLinear() )
            {
                update();
            }
            else
            {
#if defined( FEELPP_HAS_PETSC )
                M_nlsolver->dense_residual = std::bind( &Inverse::updateResidual, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
                M_nlsolver->dense_jacobian = std::bind( &Inverse::updateJacobian, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
                // find xref by solving the non linear equation
                //M_nlsolver->setType( TRUST_REGION );
                M_nlsolver->setKspSolverType( SolverType::PREONLY );
                M_nlsolver->setPreconditionerType( PreconditionerType::LU_PRECOND );
                M_nlsolver->setRelativeResidualTol( 1e-12 /*1e-16*/ );

                M_J.reset( new dense_matrix_type( M_xref.size(), M_xref.size() ) );
                M_R.reset( new dense_vector_type( M_xref.size() ) );
#endif
            }
        }

      public:
        template <typename GeoElem>
        void
        update( GeoElem const& __ge )
        {
            //M_G = ( geometric_mapping_type::nNodes == GeoElem::numVertices ) ?__ge.vertices() : __ge.G();
            M_G = __ge.G();
            if ( M_gm->isLinear() )
            {
                update();
            }
        }
        template <typename GeoElem>
        void
        update( GeoElem const& __ge, mpl::int_<1> /**/ )
        {
            M_G = __ge.vertices();
            if ( M_gm->isLinear() )
            {
                update();
            }
        }

        /** @name Accessors
 */
        //@{

        /**
   \return the geometric mapping associated with the context
*/
        geometric_mapping_ptrtype const& geometricMapping() const
        {
            return M_gm;
        }

        /**
   \return the dimension of the space of the real element
*/
        uint16_type N() const
        {
            return G().size1();
        }
        uint16_type P() const
        {
            return M_gm->dim();
        }

        /**
   \return the node in the reference element
*/
        node_t_type const& xRef() const
        {
            return M_xref;
        }

        /**
   \return the node in the real element
*/
        node_t_type const& xReal() const
        {
            return M_xreal;
        }

        /**
   \return the matrix associated with the geometric nodes
*/
        matrix_type const& G() const
        {
            return M_G;
        }

        /**
   \return the gradient of the transformation at reference node
*/
        matrix_type const& K() const
        {
            return M_K;
        }

        value_type J() const
        {
            return math::abs( det<Dim>( M_K ) );
        }

        /**
   \return the pseudo-inverse of the gradient of the
   transformation at reference node
*/
        matrix_type const& B() const
        {
            return M_B;
        }

        /**
   \return the barycenter of the reference nodes
*/
        node_type barycenterRef() const
        {
            node_type __barycenter( M_gm->dim() );
            __barycenter.clear();

            for ( uint16_type __c = 0; __c < M_gm->referenceConvex().nPoints(); ++__c )
            {
                __barycenter += M_gm->refNode( __c );
            }

            __barycenter /= M_gm->referenceConvex().nPoints();
            return __barycenter;
        }

        /**
   \return the barycenter of the geometric nodes
*/
        node_type barycenterReal() const;

        /**
   \return \c true if the node is in the convex defined by G,
   \c false otherwise
*/
        bool isIn() const
        {
            return M_is_in;
        }

        /**
 * tell whether the point is on the surface of the convex
 * @return true if the point is on the surface, false otherwise
 */
        bool isOnConvexSurface() const
        {
            if ( M_gm->isLinear() )
            {
                // x -x0 - K(0)\bar{x}
                return std::abs( ublas::norm_2( xReal() - ublas::column( M_G, 0 ) -
                                                ublas::prod( M_K, xRef() ) ) ) < 1e-10;
            }

            return false;
        }

        //@}

        /** @name  Mutators
 */
        //@{

        /**
   set the real node
*/
        void setXReal( node_type const& __xreal )
        {
            M_xreal = __xreal;

            if ( M_gm->isLinear() )
            {
                //update();
                M_is_in = linearInverse();
            }

            else
            {
                M_is_in = nonLinearInversePetsc();

                //bool isin = nonLinearInverse();
            }
        }

        matrix_node_t_type operator()( matrix_node_t_type const& real_pts, bool allow_extrapolation = false ) const
        {
            if ( trans == fem::LINEAR )
                return linearInversePoints( real_pts, allow_extrapolation );

            else
                return matrix_node_t_type();
        }

        //@}
      private:
        class GeomapInverseConvex
        {
          public:
            GeomapInverseConvex( Inverse& __gmi,
                                 const node_t_type& __xr )
                : gmi( __gmi ),
                  xreal( __xr )
            {
            }
            // evaluate
            scalar_type operator()( const node_t_type& x ) const
            {
                node_type r = gmi.geometricMapping()->transform( x, gmi.G() ) - xreal;
                return ublas::inner_prod( r, r ) / 2.;
            }
            void operator()( const node_t_type& x, node_t_type& gr ) const
            {
                //             gmi.setXreal( xreal );
                gmi.M_xref.assign( x );
                gmi.update();
                node_type r = gmi.geometricMapping()->transform( x, gmi.G() ) - xreal;
                gr.resize( x.size() );
                ublas::prod( ublas::trans( gmi.K() ), r, gr );
            }

          private:
            Inverse& gmi;
            node_type xreal;
        };

        void updateResidual( dense_vector_type const& x, dense_vector_type& r )
        {
            dense_vector_type y = M_gm->transform( x, M_G );

            if ( N() == P() )
                r = y - M_xreal;

            else
            {
                M_gm->gradient( x, M_g );
                ublas::axpy_prod( M_G, M_g, M_K );
                ublas::prod( ublas::trans( M_K ), y - M_xreal, r );
            }

#if 0
    LOG(INFO) << "[geomap::residual] begin ------------------------------\n";
    LOG(INFO) << "[geomap::residual] x =" << x << "\n";
    LOG(INFO) << "[geomap::residual] M_G =" << M_G << "\n";

    LOG(INFO) << "[geomap::residual] y =" << y << "\n";
    LOG(INFO) << "[geomap::residual] xreal =" << M_xreal << "\n";
    LOG(INFO) << "[geomap::residual] r(xreal-y) =" << r << "\n";
    LOG(INFO) << "[geomap::residual] end   ------------------------------\n";
#endif // 0
        }

        void updateJacobian( dense_vector_type const& x, dense_matrix_type& j )
        {
            M_gm->gradient( x, M_g );

            if ( N() == P() )
                ublas::axpy_prod( M_G, M_g, j );

            else
            {
                ublas::axpy_prod( M_G, M_g, M_K );
                ublas::prod( ublas::trans( M_K ), M_K, j );
            }

#if 0
    LOG(INFO) << "[geomap::jacobian] begin ------------------------------\n";
    LOG(INFO) << "[geomap::jacobian] x =" << x << "\n";
    LOG(INFO) << "[geomap::jacobian] j =" << j << "\n";
    LOG(INFO) << "[geomap::jacobian] end   ------------------------------\n";
#endif
        }

        /**
   update information on this context

   -# update the coordinate of the real points
   -# update the pseudo-inverse of the gradient of the transformation

   compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
   where \f$G \f$ is the matrix representing the geometric nodes nDof x dim

   compute \f$ B(x_{\mathrm{ref}}) = K ( K^T K )^{-1} \f$
   where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
*/
        void update()
        {
            // xref is not defined at this point: it serves only as a dummy point
            // in the case of linear inversion (__g is a constant matrix in this case)
            // in the non linear case , xRef() is update in the newton iterations
            M_gm->gradient( xRef(), M_g );
            DVLOG( 2 ) << "[update] g = " << M_g << "\n";

            checkInvariant();

            ublas::axpy_prod( M_G, M_g, M_K );
            DVLOG( 2 ) << "[update] K(0) = " << M_K << "\n";

            // compute B
            if ( M_gm->dim() != N() )
            {
                ublas::prod( ublas::trans( M_K ), M_K, M_CS );
                LU<matrix_type> __lu( M_CS );
                __lu.inverse( M_CS );
                ublas::axpy_prod( M_K, M_CS, M_B );
            }

            else
            {
                LU<matrix_type> __lu( M_K );
                __lu.inverse( M_CS );
                M_B = ublas::trans( M_CS );
            }

            DVLOG( 2 ) << "[update] B(0) = " << M_B << "\n";
        }

        void
        checkInvariant() const
        {
            FEELPP_ASSERT( M_G.size2() == M_g.size1() )
            ( M_G.size2() )( M_g.size1() ).error( "G,g invalid dimensions" );
            FEELPP_ASSERT( M_G.size1() == M_K.size1() )
            ( M_G.size1() )( M_K.size1() ).error( "G,K invalid dimensions" );
            FEELPP_ASSERT( M_g.size2() == M_K.size2() )
            ( M_g.size2() )( M_K.size2() ).error( "g,K invalid dimensions" );
            FEELPP_ASSERT( M_B.size2() == M_gm->dim() )
            ( M_B.size1() )( N() ).error( "B,gm invalid dimensions" );
        }

        bool linearInverse()
        {
            checkInvariant();

            size_type N = M_xreal.size();
            size_type P = M_xref.size();

            node_type y( M_xreal );

            DVLOG( 2 ) << "y = xreal = " << y << "\n";
            //DVLOG(2) << "G(0)  = " << node_type( M_x0 << "\n";
            y.minus_assign( ublas::column( M_G, 0 ) );
            DVLOG( 2 ) << "y - G(0) = " << y << "\n";

            DVLOG( 2 ) << "B(0) = " << M_B << "\n";
            DVLOG( 2 ) << "xref = " << ublas::prod( ublas::trans( M_B ), y ) << "\n";

            // xref = B^T * y = B^T * ( x_real - x_0)
            M_xref.assign( ublas::prod( ublas::trans( M_B ), y ) - ublas::scalar_vector<value_type>( P, 1.0 ) );

            DVLOG( 2 ) << "[GeoMap::Inverse::linearInverse] xref : " << M_xref << "\n";

            bool __isin;
            double vmin;
            boost::tie( __isin, vmin ) = M_gm->isIn( M_xref );
            DVLOG( 2 ) << "[GeoMap::Inverse::linearInverse] isIn : " << __isin << "\n";

            ///if ( __isin < 1e-10 )
            if ( __isin )
            {
                if ( N == P )
                    return true;

                else
                {
                    // y = y - K * x_ref
                    ublas::axpy_prod( M_K, -M_xref, y );

                    if ( ublas::norm_2( y ) < 1e-10 )
                        return true;
                }
            }

            return false;
        }

        matrix_node_t_type linearInversePoints( matrix_node_t_type const& real_pts, bool /*allow_extrapolation*/ = false ) const
        {
            return ublas::prod( ublas::trans( M_B ),
                                real_pts - ublas::outer_prod( ublas::column( M_G, 0 ),
                                                              ublas::scalar_vector<value_type>( real_pts.size2(), value_type( 1 ) ) ) ) -
                   ublas::scalar_matrix<value_type>( M_B.size2(), real_pts.size2(), value_type( 1 ) );
        }

        /*
 *inversion for non-linear geometric transformations
 *  (Newton on Grad(pgt)(y - pgt(x)) = 0 )
 */
        bool nonLinearInversePetsc()
            {
                //LOG(INFO) << "starting new nonlinear inverse\n";
                //const double EPS = 1e-10;
                const double IN_EPS = 1e-10;
                size_type N = M_xreal.size();
                size_type P = M_xref.size();

#if 0
                /*
                 find an initial guess: closest geometric node to M_xreal
                 */
                node_type x0 = M_gm->refNode( 0 );
                node_type y = ublas::column( M_G, 0 );
                scalar_type d = ublas::inner_prod( y, M_xreal );

                for ( size_type j = 1; j < M_gm->nbPoints(); ++j )
                {
                    scalar_type d2 = ublas::inner_prod( ublas::column( M_G, j ), M_xreal );

                    if ( d2 < d )
                    {
                        d = d2;
                        x0 = M_gm->refNode( j );
                        y = ublas::column( M_G, j );
                    }
                }

                M_xref = x0;
                node_type x0 = barycenterRef();

#else
                M_xref = barycenterRef();
#endif

#if 0
                dense_matrix_type J( P, P );
                dense_vector_type R( P );
                updateResidual( M_xref, R );
                updateJacobian( M_xref, J );

                // find xref by solving the non linear equation
                //M_nlsolver->setType( TRUST_REGION );
                M_nlsolver->setKspSolverType( SolverType::PREONLY );
                M_nlsolver->setPreconditionerType( PreconditionerType::LU_PRECOND );
                M_nlsolver->setRelativeResidualTol( 1e-16 );
#endif
                M_nlsolver->solve( *M_J, M_xref, *M_R, 1e-10, 10 );

                // compute the location of xref: inside or outside the element
                bool __isin;
                double vmin;
                this->updateResidual( M_xref, *M_R );
                boost::tie( __isin, vmin ) = M_gm->isIn( M_xref );

                if ( __isin &&
                     ( P == N || ublas::norm_2( *M_R ) < IN_EPS ) )
                {
                    //LOG(INFO) << "point " << M_xref << "in IN (" << vmin << ") residual = " << ublas::norm_2(R) << "\n";
                    return true;
                }

                else
                {
                    //LOG(INFO) << "point " << M_xref << "in OUT (" << vmin << ") residual = " << ublas::norm_2(R) << "\n";
                }

                //LOG(INFO) << "done in new nonlinear inverse\n";
                return false;
            }
        /*
 *inversion for non-linear geometric transformations
 *  (Newton on Grad(pgt)(y - pgt(x)) = 0 )
 */
        bool nonLinearInverse()
            {
                const double EPS = 1e-10;
                const double IN_EPS = 1e-10;
                size_type N = M_xreal.size();
                size_type P = M_xref.size();

                /*
                 find an initial guess: closest geometric node to M_xreal
                 */
                node_type x0 = M_gm->refNode( 0 );
                node_type y = ublas::column( M_G, 0 );
                scalar_type d = ublas::inner_prod( y, M_xreal );

                for ( size_type j = 1; j < M_gm->nbPoints(); ++j )
                {
                    scalar_type d2 = ublas::inner_prod( ublas::column( M_G, j ), M_xreal );

                    if ( d2 < d )
                    {
                        d = d2;
                        x0 = M_gm->refNode( j );
                        y = ublas::column( M_G, j );
                    }
                }

                M_xref = x0;

                node_type vres( N );
                node_type rn( M_xreal );
                rn.minus_assign( y );

                this->update();

                // vres = K^T rn
                ublas::axpy_prod( rn, M_K, vres );
                scalar_type res = ublas::norm_2( vres );

                //     std::cerr << "DEBUT: res0=" << res << ", X=" << M_xreal << "\nB=" << M_B << ", K=" << M_K << "\n";
                unsigned cnt = 50;

                while ( res > EPS / 10 && cnt > 0 )
                {
                    node_type xn( P );

                    // xn = B^T rn
                    ublas::axpy_prod( rn, M_B, xn );

                    scalar_type newres;

                    for ( int16_type i = 1; i <= 256; i *= 2 )
                    {
                        M_xref.plus_assign( xn / scalar_type( i ) );
                        y = M_gm->transform( M_xref, G() );

                        rn = M_xreal - y;

                        this->update();

                        if ( P != N )
                        {
                            // vres = K^T rn
                            ublas::axpy_prod( rn, M_K, vres );
                            newres = ublas::norm_2( vres );
                        }

                        else
                        {
                            newres = ublas::norm_2( rn );
                        }

                        if ( newres < 1.5 * res )
                            break;
                    }

                    res = newres;
                    //         std::cout << "cnt=" << cnt << ", x=" << M_xref << ", res=" << res << "\n";
                    --cnt;
                }

                //     std::cerr << " invert_nonlin done\n";
                //     std::cerr << "cnt=" << cnt << ", P=" << P << ", N=" << N
                //               << "\nX=" << M_xreal << " Xref=" << M_xref << "\nresidu=" << res << "\n";
                //<< ", G=" << G << "\nX=" << M_xreal << " Xref=" << x << "\nresidu=" << res << "\nB=" << B << ", K=" << K << "\n"  << "\n-------------------^^^^^^^^\n";
                if ( cnt == 0 )
                {
#if 0
                    std::cerr << "BFGS in geotrans_inv_convex!\n";
                    GeomapInverseConvex b( *this, M_xreal );

                    iteration_ptrtype iter( Iteration<double>::New() );
                    iter->setMaximumNumberOfIterations( 50 );
                    iter->setRelativePrecision( 1e-8 );

                    node_type x( x0 );
                    bfgs( b,b,x,10,*iter );
                    rn = M_gm->transform( x,G() ) - M_xreal;

                    if ( M_gm->isIn( x ) < IN_EPS &&
                         N==P && ublas::norm_2( rn ) > IN_EPS )
                        throw "inversion of non-linear geometric transformation "
                            "failed ! (too much iterations)";

#endif
                }

                bool __isin;
                double vmin;
                boost::tie( __isin, vmin ) = M_gm->isIn( M_xref );

                if ( __isin &&
                     ( P == N || ublas::norm_2( rn ) < IN_EPS ) )
                {

                    return true;
                }

                else
                {
                    //LOG(INFO) << "point " << M_xref << "in OUT (" << vmin << ")\n";
                }

                return false;
            }

      private:
        geometric_mapping_ptrtype M_gm;
        node_type M_xref;
        node_type M_xreal;

        bool M_is_in;

        //matrix_type const& M_G;
        matrix_type M_G;
        matrix_type M_K;
        matrix_type M_B;
        matrix_type M_CS;
        matrix_type M_g;

        std::shared_ptr<SolverNonLinear<double>> M_nlsolver;
        std::shared_ptr<dense_matrix_type> M_J;
        std::shared_ptr<dense_vector_type> M_R;

    }; // Inverse

  private:
    element_gm_ptrtype _elementMap;
    face_gm_ptrtype _boundaryMap;
    matrix_type M_g_linear;

    friend class Inverse;

    reference_convex_type M_refconvex;
};

#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>

template <int Dim, int Order, int RealDim, template <uint16_type, uint16_type, uint16_type> class Entity = Simplex, typename T = double>
struct GT_Lagrange
{
};

template <int Dim, int Order, int RealDim, typename T = double>
struct GT_QK
{};

# /* List of dims. */
# define FEELPP_GEOMAP                                    \
    BOOST_PP_TUPLE_TO_LIST(                               \
        1,                                                \
        (                                                 \
            Lagrange                                      \
                                                          ) \
                                                          ) \
    /**/
# /* List of dims. */
# define FEELPP_DIMS                                      \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           4,                           \
                           (                            \
                               0,1,2,3                  \
                                                    )   \
                                                    )   \
    /**/
# /* List of real dims. */
# define FEELPP_REALDIMS                                  \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           3,                           \
                           (                            \
                            1,2,3                       \
                                                    )   \
                                                    )   \

# /* List of real dims. */

# define FEELPP_NEWDIMS                                           \
    BOOST_PP_TUPLE_TO_LIST(                                     \
                           9,                                   \
                           (                                    \
                            (0,1),(0,2),(0,3),                  \
                            (1,1),(1,2),(1,3),                  \
                            (2,2),(2,3),                        \
                            (3,3)                                   \
                                                                )   \
                                                            )   \
    /**/
# /* List of orders. */
# define FEELPP_ORDERS                                    \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           5,                           \
                           (                            \
                            1,2,3,4,5                   \
                                                    )   \
                                                    )   \
    /**/
#
#define FEELPP_ENTITY BOOST_PP_TUPLE_TO_LIST( 2, ( Simplex, Hypercube ) )
/**/
#
#define FEELPP_GEN_GT( GEOM, LDIM, LORDER ) \
    "GT_" #GEOM "(" #LDIM "," #LORDER ")" /**/
#
#/* Generates code for all dim and order. */
#define FEELPP_GT_FACTORY_OP( _, GDO ) \
    FEELPP_GT_FACTORY GDO /**/
#
#
#define FEELPP_GT_DIM( T )         \
    BOOST_PP_TUPLE_ELEM( 2, 0, T ) \
/**/
#define FEELPP_GT_REALDIM( T )     \
    BOOST_PP_TUPLE_ELEM( 2, 1, T ) \
/**/
#
#
#define FEELPP_GT_FACTORY( GEOM, LDIMS, LORDER, ENTITY )                                                        \
    template <typename T>                                                                                       \
    struct BOOST_PP_CAT( GT_, GEOM )<FEELPP_GT_DIM( LDIMS ), LORDER, FEELPP_GT_REALDIM( LDIMS ), ENTITY, T>     \
        : public GeoMap<FEELPP_GT_DIM( LDIMS ), LORDER, FEELPP_GT_REALDIM( LDIMS ), T, ENTITY, GEOM>            \
    {                                                                                                           \
        static inline const uint16_type nDim = FEELPP_GT_DIM( LDIMS );                                                 \
        static inline const uint16_type order = LORDER;                                                                \
        static inline const uint16_type nRealDim = FEELPP_GT_REALDIM( LDIMS );                                         \
        static inline const uint16_type nRealDimCheck2d = mpl::if_<mpl::less_equal<mpl::int_<2>, mpl::int_<nRealDim>>, \
                                                            mpl::int_<nRealDim>,                                \
                                                            mpl::int_<nDim>>::type::value;                      \
                                                                                                                \
        typedef mpl::vector<boost::none_t,                                                                      \
                            boost::none_t,                                                                      \
                            GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>,                                       \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM>>                                \
            geomap_faces_t;                                                                                     \
        typedef mpl::vector<GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>,                                       \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM>,                                \
                            GeoMap<3, LORDER, 3, T, ENTITY, GEOM>,                                              \
                            boost::none_t>                                                                      \
            geomap_elements_t;                                                                                  \
        typedef typename type_traits<T>::value_type value_type;                                                 \
        typedef GeoMap<nDim, LORDER, nRealDim, T, ENTITY, GEOM> super;                                          \
        typedef BOOST_PP_CAT( GT_, GEOM )<nDim - 1, LORDER, nRealDim, ENTITY, T> face_geo_type;                 \
                                                                                                                \
        static inline const uint16_type nDof = super::nDof;                                                            \
        static inline const uint16_type nNodes = super::nNodes;                                                        \
        typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim>>::type element_gm_type;                     \
        typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim>>::type face_gm_type;                           \
        template <int N>                                                                                        \
        struct face_gm                                                                                          \
        {                                                                                                       \
            typedef typename mpl::at<geomap_faces_t, mpl::int_<N>>::type type;                                  \
        };                                                                                                      \
                                                                                                                \
        BOOST_PP_CAT( GT_, GEOM )                                                                               \
        ()                                                                                                      \
            : super( std::make_shared<element_gm_type>(), std::make_shared<face_gm_type>() )                \
        {                                                                                                       \
        }                                                                                                       \
    };                                                                                                          \
/**/

template<typename Convex,typename T>
typename mpl::if_<Feel::is_simplex<Convex>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim,1,Convex::nDim,Simplex,T>>>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim,1,Convex::nDim,Hypercube,T>>> >::type::type
makeGeometricTransformation()
{
    using gmt = typename mpl::if_<is_simplex<Convex>,
                                  mpl::identity<GT_Lagrange<Convex::nDim,1,Convex::nDim,Simplex,T>>,
                                  mpl::identity<GT_Lagrange<Convex::nDim,1,Convex::nDim,Hypercube,T>> >::type::type;
    return std::make_shared<gmt>();
}

#if 0
#define FEELPP_GT_FACTORY( GEOM, LDIM, LORDER, LREALDIM, ENTITY ) \
    FEELPP_GT_FACTORY_BIS(GEOM,LDIM,LORDER,LREALDIM,ENTITY)),             \
                 BOOST_PP_EMPTY )                                       \
    /**/
#endif

//BOOST_PP_LIST_FOR_EACH_PRODUCT(FEELPP_GT_FACTORY_OP, 5, (FEELPP_GEOMAP, FEELPP_DIMS, FEELPP_ORDERS, FEELPP_REALDIMS, FEELPP_ENTITY))
BOOST_PP_LIST_FOR_EACH_PRODUCT( FEELPP_GT_FACTORY_OP, 4, ( FEELPP_GEOMAP, FEELPP_NEWDIMS, FEELPP_ORDERS, FEELPP_ENTITY ) )

#undef FEELPP_DIMS
#undef FEELPP_ORDERS
#undef FEELPP_REALDIMS
#undef FEELPP_NEWDIMS

template <typename Elem, template <uint16_type, uint16_type, uint16_type> class Entity = Simplex, typename T = double>
class RealToReference
{
  public:
    static inline const uint16_type nDim = Elem::nDim;
    static inline const uint16_type nRealDim = Elem::nRealDim;

    typedef T value_type;
    typedef typename matrix_node<T>::type points_type;
    typedef GT_Lagrange<nDim, 1, nRealDim, Entity, value_type> gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::Inverse inverse_gm_type;

    RealToReference( Elem const& elem )
        : M_gm( new gm_type ),
          M_igm( M_gm, elem )
    {
    }

    points_type operator()( points_type const& pts ) const
    {
        return M_igm( pts );
    }

    value_type J() const
    {
        return M_igm.J();
    }

  private:
    gm_ptrtype M_gm;
    inverse_gm_type M_igm;
};

} // Feel
#endif
