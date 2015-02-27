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
#if (BOOST_VERSION >= 103400)
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <Eigen/Core>

#include <boost/mpl/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/enable_shared_from_this.hpp>

//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/atlas/cblas3.hpp>

#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/svd.hpp>
#include <feel/feelalg/solvernonlinear.hpp>
#include <feel/feelalg/lu.hpp>
#include <feel/feelalg/iteration.hpp>

#include <feel/feelpoly/context.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feelmesh/marker.hpp>



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

//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;

struct GeomapInverse
{
    //GeomapInverse
};

template<uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTag > class Lagrange;


/**
 * \class GeoMap
 * \brief Structure for the geometrical mapping
 * \author C. Prud'homme
 *
 * This class contains the geometrical transformation that maps the reference
 * element on the current element, and its values on integration points
 *
 */
template<uint16_type Dim,
         uint16_type Order,
         uint16_type RealDim,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex,
         template<uint16_type, template<uint16_type RDim> class PolySetType, typename ContinuityType,
         template<class, uint16_type, class> class Pts, uint16_type > class PP = Lagrange>
class GeoMap
    :
    public PP<Order,Scalar, Continuous,PointSetEquiSpaced, 0>::template apply<Dim,RealDim/*Dim*/, T, Entity<Dim,Order,/*RealDim*/Dim> >::result_type//,
//public boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > >
       //public PP<Order,Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type
       {
           //typedef typename PP<Order, Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type super;
           typedef typename PP<Order, Scalar, Continuous, PointSetEquiSpaced, 0>::template apply<Dim, RealDim/*Dim*/, T, Entity<Dim,Order,/*RealDim*/Dim> >::result_type super;

           //typedef boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > > super_enable_this;

           static const uint16_type nRealDimCheck2d = mpl::if_< mpl::less_equal<mpl::int_<2>,mpl::int_<RealDim> >,
           mpl::int_<RealDim>,
           mpl::int_<Dim> >::type::value;


           typedef mpl::vector<boost::none_t, boost::none_t, boost::none_t,
           GeoMap<1, Order, RealDim, T, Entity, PP> > geomap_edges_t;

           typedef mpl::vector<boost::none_t, boost::none_t,
           GeoMap<1, Order, RealDim, T, Entity, PP>,
           GeoMap<2, Order, nRealDimCheck2d/*RealDim*/, T, Entity, PP> > geomap_faces_t;

           typedef mpl::vector<GeoMap<1, Order, RealDim,T, Entity, PP>,
           GeoMap<2, Order, nRealDimCheck2d/*RealDim*/,T, Entity, PP>,
           GeoMap<3, Order, 3,T, Entity, PP>,
           boost::none_t> geomap_elements_t;
           public:

           static const uint16_type nDim = super::nDim;
           static const uint16_type nRealDim = super::nRealDim;
           static const uint16_type nDof = super::nDof;
           static const uint16_type nNodes = super::nNodes;
           static const fem::transformation_type trans = super::trans;

           typedef typename super::value_type value_type;

           typedef typename super::PreCompute precompute_type;
           typedef boost::shared_ptr<precompute_type> precompute_ptrtype;

           typedef typename super::convex_type convex_type;
           typedef Reference<convex_type,nDim,Order,nDim/*nRealDim*/> reference_convex_type;

           typedef GeoMap<Dim, Order, RealDim,T, Entity, PP > self_type;
           typedef self_type geometric_mapping_type;
           typedef boost::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

           typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim> >::type element_gm_type;
           typedef boost::shared_ptr<element_gm_type> element_gm_ptrtype;
           typedef element_gm_ptrtype gm_ptrtype;

           typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim> >::type face_gm_type;
           typedef boost::shared_ptr<face_gm_type> face_gm_ptrtype;

           template<int N>
           struct face_gm
           {
               typedef typename mpl::at<geomap_faces_t, mpl::int_<N> >::type type;
           };


           typedef typename mpl::at<geomap_edges_t, mpl::int_<nDim> >::type edge_gm_type;
           typedef boost::shared_ptr<edge_gm_type> edge_gm_ptrtype;

           template<int N>
           struct edge_gm
           {
               typedef typename mpl::at<geomap_edges_t, mpl::int_<N> >::type type;
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
           typedef typename ublas::matrix<value_type> matrix_type;



           /** default constructor */
           GeoMap()
               :
               super(),
               M_is_cached( false ),
               _elementMap( ),
               _boundaryMap(),
               M_g_linear( nNodes, nDim ),
               M_refconvex()
{
    if ( trans == fem::LINEAR )
    {
        //M_g_linear.resize( nNodes, nDim );
        matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );

        ublas::vector<ublas::matrix<value_type> > m = super::derivate( __dummy_pts );

        FEELPP_ASSERT( M_g_linear.size2() == m.size() )( M_g_linear.size2() )(  m.size() ).error( "invalid dimension" );
        FEELPP_ASSERT( m( 0 ).size2() == 1 )( m( 0 ).size2() ).error( "Invalid number of points" );

        FEELPP_ASSERT( M_g_linear.size1() == m( 0 ).size1() )( M_g_linear.size1() )( m( 0 ).size1() ).error( "invalid number of DOF" );

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
GeoMap( element_gm_ptrtype const& e,  face_gm_ptrtype const& f )
    :
    super(),
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
        ublas::vector<ublas::matrix<value_type> > m = super::derivate( __dummy_pts );
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
/**
   destructor
*/
~GeoMap()
{}

/**
   \return the dimension of the underlying element
*/
uint16_type dim() const
{
    return nDim;
}

uint16_type realDim() const
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
ublas::vector<value_type>  refNode( uint16_type i ) const
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
boost::tuple<bool,value_type> isIn( typename node<value_type>::type const& pt ) const
{
    return M_refconvex.isIn( pt );
}


/**
 *  apply the geometric mapping to the point \c pt given the real
 *  geometric nodes stored in a NxNg matrix \c G
 */
node_t_type transform( const node_t_type &__ref_p,
                       matrix_node_t_type const& __G ) const
{
    namespace lambda = boost::lambda;

    typename node<value_type>::type __real_p(  __G.size1() );
    __real_p.clear();

    for ( uint16_type __i = 0; __i < nNodes; ++__i )
    {
        //value_type __phi_at_pt = super::phi( __i, __ref_p );
        value_type __phi_at_pt = super::evaluate( __i, __ref_p )( 0 );
        __real_p.plus_assign( __phi_at_pt*ublas::column( __G, __i ) );
    }

    return __real_p;
}

/**
 * apply the geometric mapping to the point index \c id_pt given the real
 * geometric nodes stored in a NxNg matrix \c G
 */
node_t_type transform( uint16_type __idref,
                       matrix_node_t_type const& __G,
                       precompute_type const* __pc ) const
{
    node_t_type __real_p( __G.size1() );
    __real_p.clear();

    for ( uint16_type __i = 0; __i < nNodes; ++__i )
    {
        // evaluate transformation at point pt
        value_type __phi_at_pt = __pc->phi( __idref, __i );
        __real_p.plus_assign(  __phi_at_pt*ublas::column( __G, __i ) );
    }

    return __real_p;
}

/**
 * compute real coordinates from a matrix of ref coordinates
 */
void transform( matrix_node_t_type const& G,
                precompute_type const* pc,
                matrix_type & x ) const
{
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> P( G.data(), G.rows(), G.cols() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> Phi( pc->phi().data(), pc->phi().rows(), pc->phi().cols() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> X( x.data(), x.size() );
        X  = P*Phi;
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
    namespace lambda = boost::lambda;

    if ( trans == fem::LINEAR )
    {
        __g = M_g_linear;
    }

    else
    {
        FEELPP_ASSERT( __pt.size() == dim() )( __pt.size() )( dim() ).error( "invalid dimension" );

        matrix_node_t_type __pts( nDim, 1 );
        ublas::column( __pts, 0 ) = __pt;

        ublas::vector<ublas::matrix<value_type> > m = super::derivate( __pts );

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
    FEELPP_ASSERT( __pc )( __idref ).error( "a PreCompute must be set first before using this function" );

    if ( trans == fem::LINEAR )
    {
        __g = M_g_linear;
    }

    else
    {
        FEELPP_ASSERT( __pc->dim() == dim() )( __pc->dim() )( dim() ).error( "invalid dimension" );

        for ( size_type i = 0; i < nNodes; ++i )
        {
            for ( uint16_type n = 0; n < nDim; ++n )
            {
                __g( i, n ) = __pc->grad( i,  0, n, __idref );
            }
        }
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

    FEELPP_ASSERT(  __g.size1() == G.size2() )( __g.size1() )( G.size2() ).error( "invalid sizes" );

    matrix_node_t_type K( G.size1(), __g.size2() );
    ublas::axpy_prod( G, __g, K );

    SVD<matrix_node_t_type> __svd( K );
    value_type __max;
    value_type __min;
    __svd.conditionNumber( __max, __min );
    return __max*sqrt( value_type( N ) ) / value_type( N );
}


bool M_is_cached;
bool isCached() const
{
    return M_is_cached;
}

std::vector<bool> M_cached;
std::vector<double> M_J;
//boost::multi_array<double,3> M_K;
//boost::multi_array<double,3> M_B;
std::vector<matrix_type> M_K;
std::vector<matrix_type> M_B;

template<typename MeshType>
void initCache( MeshType const* mesh  )
{
    size_type nelts = mesh->numElements();
    LOG(INFO) << "[Geomap] start caching J,K,B for " << nelts << " elements\n";
    M_cached.resize( nelts );
    std::fill( M_cached.begin(), M_cached.end(), false );

    M_J.resize( nelts );
    M_K.resize( nelts );
    M_B.resize( nelts );
    M_is_cached = true;
}
bool isCacheValid() const
{
    if ( !( M_cached.size() > 0 &&
            M_J.size() > 0 &&
            M_K.size() > 0 &&
            M_B.size() > 0 ) )
    {
        LOG(INFO) << "invalid cache\n";
        return false;
    }

    return true;
}
bool cached( int e ) const
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    return M_cached[e];
}
void setCached( int e, bool v )
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    M_cached[e] = v;
}
double J( int e ) const
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_J.size() ) << "invalid element id " << e << "( " << M_J.size() << ") for geomap cache, are you using the proper mesh\n";
    return M_J[e];
}
matrix_type const& B( int e ) const
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_B.size() ) << "invalid element id " << e << "( " << M_B.size() << ") for geomap cache, are you using the proper mesh\n";
    return M_B[e];
}
matrix_type const& K( int e ) const
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_K.size() ) << "invalid element id " << e << "( " << M_K.size() << ") for geomap cache, are you using the proper mesh\n";
    return M_K[e];
}
void addJ( int e, double v )
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_J.size() ) << "invalid element id " << e << "( " << M_J.size() << ") for geomap cache, are you using the proper mesh\n";
    M_J[e] = v;
}
void addK( int e, matrix_type const& K )
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_K.size() ) << "invalid element id " << e << "( " << M_K.size() << ") for geomap cache, are you using the proper mesh\n";
    M_K[e].resize( K.size1(), K.size2() );

    for ( size_type i = 0; i < K.size1(); ++i )
        for ( size_type j = 0; j < K.size2(); ++j )
            M_K[e]( i, j ) = K( i, j );
}
void addB( int e, matrix_type const& B )
{
    FEELPP_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
    DCHECK( e >= 0 && e < M_B.size() ) << "invalid element id " << e << "( " << M_B.size() << ") for geomap cache, are you using the proper mesh\n";
    M_B[e].resize( B.size1(), B.size2() );

    for ( size_type i = 0; i < B.size1(); ++i )
        for ( size_type j = 0; j < B.size2(); ++j )
            M_B[e]( i, j ) = B( i, j );
}

/**
 * \class Context
 *
 * Context for the geometric mapping depend on a node in the
 * reference element
 */
template<size_type context_v, typename ElementType>
class Context
{
    public:
    static const size_type contextv = context_v;
    static const size_type context = context_v;
    // reference space dimension
    static const uint16_type PDim = ElementType::nDim;
    // real space dimension
    static const uint16_type NDim = ElementType::nRealDim;
    static const uint16_type nDim = NDim;
    // type of transformation (linear or not)
    static const fem::transformation_type trans = geometric_mapping_type::trans;
    static const bool is_linear = ( trans == fem::LINEAR );

    static const bool condition = ( ( PDim==NDim )||( ( NDim>=1 )&&( PDim==NDim-1 ) ) );
    //BOOST_MPL_ASSERT_MSG( condition, INVALID_DIM, (mpl::int_<NDim>, mpl::int_<PDim>, ElementType ) );
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim> >,
    mpl::identity<GeoMap<Dim, Order, NDim, T, Entity, PP > >,
    typename mpl::if_<mpl::and_<mpl::greater_equal<mpl::int_<NDim>, mpl::int_<1> >,
    mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> > >,
    //typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
    // mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
    mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP >::face_gm_type>,
    typename mpl::if_<mpl::and_<mpl::equal_to<mpl::int_<NDim>, mpl::int_<3> >,
    mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-2> > >,
    //typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
    // mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
    mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP >::edge_gm_type>,
    mpl::identity<boost::none_t> >::type>::type>::type::type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;

    typedef typename gm_type::value_type value_type;

    typedef typename gm_type::precompute_ptrtype precompute_ptrtype;

    typedef Context<contextv,ElementType> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    typedef node_t_type normal_type;
    typedef ublas::vector<normal_type> normals_type;
    typedef typename normals_type::const_iterator normal_const_iterator;
    typedef node_t_type tangent_type;
    typedef ublas::vector<tangent_type> tangents_type;
    typedef typename tangents_type::const_iterator tangent_const_iterator;

    typedef ElementType element_type;
    typedef typename element_type::permutation_type permutation_type;

    Context( gm_ptrtype __gm,
             element_type const& __e,
             precompute_ptrtype const& __pc )
        :
        M_gm( __gm ),
        M_element( boost::addressof( __e ) ),
        M_pc( __pc ),
        M_pc_faces(),
        M_npoints( M_pc->nPoints() ),

        //M_xref( PDim ),
        //M_xreal( NDim ),
        //M_x0( NDim ),
        M_J( 0 ),
        M_G( ( gm_type::nNodes == element_type::numVertices ) ?__e.vertices() : __e.G() ),
        M_n( M_gm->referenceConvex().normals() ),
        M_n_real( NDim ),
        M_un_real( NDim ),
        M_n_norm( 0 ),
        M_t_real( NDim ),
        M_ut_real( NDim ),
        M_t_norm( 0 ),
        M_xrefq( PDim, nPoints() ),
        M_xrealq( NDim, nPoints() ),
        M_n_realq( NDim, nPoints() ),
        M_un_realq( NDim, nPoints() ),
        M_n_normq( nPoints() ),

        M_g( M_G.size2(), PDim ),
        M_K( NDim, PDim ),
        M_CS( PDim, PDim ),
        M_CSi( PDim, PDim ),
        M_B( NDim, PDim ),
        M_Ptangent( NDim, NDim ),
        M_B3( boost::extents[NDim][NDim][PDim][PDim] ),
        M_id( __e.id() ),
        M_e_marker( __e.marker() ),
        M_e_marker2( __e.marker2() ),
        M_e_marker3( __e.marker3() ),
        M_elem_id_1( invalid_size_type_value ),// __e.ad_first() ),
        M_pos_in_elem_id_1( invalid_uint16_type_value ),  //__e.pos_first() ),
        M_elem_id_2( invalid_size_type_value ),  //__e.ad_second() ),
        M_pos_in_elem_id_2( invalid_uint16_type_value ),  //__e.pos_second() ),
        M_face_id( invalid_uint16_type_value ),
        M_h( __e.h() ),
        M_h_min( __e.hMin() ),
        M_h_face( 0 ),
        M_meas( __e.measure() ),
        M_measface( 0 ),
        M_Jt(),
        M_Bt(),
        M_perm( )
{

    if ( is_linear )
    {
        M_gm->gradient( node_t_type(), M_g_linear );
    }

    else
    {
        M_Jt.resize( nPoints() );
        M_Bt.resize( nPoints() );
    }

    update( __e );
}

Context( gm_ptrtype __gm,
         element_type const& __e,
         std::vector<std::map<permutation_type, precompute_ptrtype> > & __pc,
         uint16_type __f )
    :
    M_gm( __gm ),
    M_element( boost::addressof( __e ) ),
    M_pc(),
    M_pc_faces( __pc ),
    M_npoints( __pc[__f][__e.permutation( __f )]->nPoints() ),

    //M_xref( PDim ),
    //M_xreal( NDim ),
    //M_x0( NDim ),
    M_J( 0 ),
    M_G( ( gm_type::nNodes == element_type::numVertices ) ?__e.vertices() : __e.G() ),
    M_n( M_gm->referenceConvex().normals() ),
    M_n_real( NDim ),
    M_un_real( NDim ),
    M_n_norm( 0 ),
    M_t_real( NDim ),
    M_ut_real( NDim ),
    M_t_norm( 0 ),
    M_xrefq( PDim, nPoints() ),
    M_xrealq( NDim, nPoints() ),
    M_n_realq( NDim, nPoints() ),
    M_un_realq( NDim, nPoints() ),
    M_n_normq( nPoints() ),

    M_g( M_G.size2(), PDim ),
    M_K( NDim, PDim ),
    M_CS( PDim, PDim ),
    M_CSi( PDim, PDim ),
    M_B( NDim, PDim ),
    M_Ptangent( NDim, NDim ),
    M_B3( boost::extents[NDim][NDim][PDim][PDim] ),
    M_id( __e.id() ),
    M_e_marker( __e.marker() ),
    M_e_marker2( __e.marker2() ),
    M_e_marker3( __e.marker3() ),
    M_elem_id_1( invalid_size_type_value ),// __e.ad_first() ),
    M_pos_in_elem_id_1( invalid_uint16_type_value ),  //__e.pos_first() ),
    M_elem_id_2( invalid_size_type_value ),  //__e.ad_second() ),
    M_pos_in_elem_id_2( invalid_uint16_type_value ),  //__e.pos_second() ),
    M_face_id( __f ),
    M_h( __e.h() ),
    M_h_min( __e.hMin() ),
    M_h_face( 0 ),
    M_meas( __e.measure() ),
    M_measface( __e.faceMeasure( __f ) ),
    M_Jt(),
    M_Bt(),
    M_perm( )
{

    if ( is_linear )
    {
        M_gm->gradient( node_t_type(), M_g_linear );
    }

    else
    {
        M_Jt.resize( nPoints() );
        M_Bt.resize( nPoints() );
    }

    update( __e, __f );
}
Context( gmc_ptrtype& p )
    :
    M_gm( p->M_gm ),
    M_element( p->M_element ),
    M_pc( p->M_pc ),
    M_pc_faces( p->M_pc_faces ),
    M_npoints( M_pc->nPoints() ),
    //M_xref( PDim ),
    //M_xreal( NDim ),
    //M_x0( NDim ),
    M_J( p->M_J ),
    M_G( ( gm_type::nNodes == element_type::numVertices ) ?M_element->vertices() : M_element->G() ),
    M_n( p->M_n ),
    M_n_real( p->M_n_real ),
    M_un_real( p->M_un_real ),
    M_n_norm( p->M_n_norm ),
    M_t_real( p->M_t_real ),
    M_ut_real( p->M_ut_real ),
    M_t_norm( p->M_t_norm ),
    M_xrefq( p->M_xrefq ),
    M_xrealq( p->M_xrealq ),
    M_n_realq( p->M_n_realq ),
    M_un_realq( p->M_un_realq ),
    M_n_normq( p->M_n_normq ),
    M_g( p->M_g ),
    M_K( p->M_K ),
    M_CS( p->M_CS ),
    M_CSi( p->M_CSi ),
    M_B( p->M_B ),
    M_Ptangent( p->M_Ptangent ),
    M_B3( p->M_B3 ),
    M_id( p->M_id ),
    M_e_marker( p->M_e_marker ),
    M_e_marker2( p->M_e_marker2 ),
    M_e_marker3( p->M_e_marker3 ),
    M_elem_id_1( invalid_size_type_value ),// M_element.ad_first() ),
    M_pos_in_elem_id_1( invalid_uint16_type_value ),  //M_element.pos_first() ),
    M_elem_id_2( invalid_size_type_value ),  //M_element.ad_second() ),
    M_pos_in_elem_id_2( invalid_uint16_type_value ),  //M_element.pos_second() ),
    M_face_id( invalid_uint16_type_value ),
    M_h( p->M_h ),
    M_h_min( p->M_h_min ),
    M_h_face( p->M_h_face ),
    M_meas( p->M_meas ),
    M_measface( p->M_measface ),
    M_Jt(),
    M_Bt(),
    M_perm( p->M_perm )
{
    if ( is_linear )
    {
        M_gm->gradient( node_t_type(), M_g_linear );
    }

    else
    {
        M_Jt.resize( nPoints() );
        M_Bt.resize( nPoints() );
    }

    update( *M_element );
}

/**
 * clone this context
 */
gmc_ptrtype clone()
{
    return gmc_ptrtype( new gmc_type( *this ) );
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
void update( element_type const& __e, uint16_type __f )
{
    //M_element_c = boost::shared_ptr<element_type const>(&__e);
    M_element = boost::addressof( __e );
    M_face_id = __f;

    M_perm = __e.permutation( M_face_id );

    M_h_face = __e.hFace( M_face_id );
    //M_h_edge = __e.hEdge( M_face_id );

    M_pc = M_pc_faces[__f][M_perm];
    //M_G = __e.G();
    M_G = ( gm_type::nNodes == element_type::numVertices ) ?__e.vertices() : __e.G();
    M_id = __e.id();
    M_e_marker = __e.marker();
    M_e_marker2 = __e.marker2();
    M_e_marker3 = __e.marker3();
    M_h = __e.h();
    M_h_min = __e.hMin();
    M_meas = __e.measure();
    M_measface = __e.faceMeasure( __f );
    M_xrefq = M_pc->nodes();

    FEELPP_ASSERT( M_G.size2() == M_gm->nbPoints() )( M_G.size2() )( M_gm->nbPoints() ).error( "invalid dimensions" );
    FEELPP_ASSERT( M_pc ).error( "invalid precompute data structure" );

    if ( vm::has_point<context>::value )
    {

        //ublas::axpy_prod( M_G, pc->phi(), M_xrealq, true );
        std::fill( M_xrealq.data().begin(), M_xrealq.data().end(), value_type( 0 ) );
        const uint16_type size1 = M_G.size1();
        const uint16_type size3 = M_G.size2();
        const uint16_type size2 = M_pc->nPoints();

        for ( uint16_type i = 0; i < size1; ++i )
            for ( uint16_type j = 0; j < size2; ++j )
            {
                for ( uint16_type k = 0; k < size3; ++k )
                    M_xrealq( i, j ) += M_G( i, k ) * M_pc->phi()[k][j]( 0,0 );
            }
    }

    if ( vm::has_jacobian<context>::value )
    {
        updateJKBN( mpl::bool_<is_linear>() );
    }

}

void update( element_type const& __e, uint16_type __f, permutation_type __perm, bool __updateJacobianCtx=true )
{
    //M_element_c = boost::shared_ptr<element_type const>(&__e);
    M_element = boost::addressof( __e );
    M_face_id = __f;

    M_perm = __perm;

    M_h_face = __e.hFace( M_face_id );
    //M_h_edge = __e.hEdge( M_face_id );

    M_pc = M_pc_faces[__f][M_perm];
    //M_G = __e.G();
    M_G = ( gm_type::nNodes == element_type::numVertices ) ?__e.vertices() : __e.G();
    M_id = __e.id();
    M_e_marker = __e.marker();
    M_e_marker2 = __e.marker2();
    M_e_marker3 = __e.marker3();
    M_h = __e.h();
    M_h_min = __e.hMin();
    M_meas = __e.measure();
    M_measface = __e.faceMeasure( __f );
    M_xrefq = M_pc->nodes();

    FEELPP_ASSERT( M_G.size2() == M_gm->nbPoints() )( M_G.size2() )( M_gm->nbPoints() ).error( "invalid dimensions" );
    FEELPP_ASSERT( M_pc ).error( "invalid precompute data structure" );

    if ( vm::has_point<context>::value )
    {

        //ublas::axpy_prod( M_G, pc->phi(), M_xrealq, true );
        std::fill( M_xrealq.data().begin(), M_xrealq.data().end(), value_type( 0 ) );
        const uint16_type size1 = M_G.size1();
        const uint16_type size3 = M_G.size2();
        const uint16_type size2 = M_pc->nPoints();

        for ( uint16_type i = 0; i < size1; ++i )
            for ( uint16_type j = 0; j < size2; ++j )
            {
                for ( uint16_type k = 0; k < size3; ++k )
                    M_xrealq( i, j ) += M_G( i, k ) * M_pc->phi()[k][j]( 0,0 );
            }
    }

    if ( vm::has_jacobian<context>::value && __updateJacobianCtx )
    {
        updateJKBN( mpl::bool_<is_linear>() );
    }

}

void update( element_type const& __e,
             precompute_ptrtype const& __pc )
{
    M_pc = __pc;


    if ( M_npoints != M_pc->nPoints() )
    {
        M_npoints = M_pc.get()->nPoints();

        M_xrefq.resize( PDim, nPoints() );
        M_xrealq.resize( NDim, nPoints() );
        M_n_realq.resize( NDim, nPoints() );
        M_un_realq.resize( NDim, nPoints() );
        M_n_normq.resize( nPoints() );

        if ( is_linear )
        {
            M_gm->gradient( node_t_type(), M_g_linear );
        }

        else
        {
            M_Jt.resize( nPoints() );
            M_Bt.resize( nPoints() );
        }
    }

    update( __e );
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
void update( element_type const& __e )
{
    M_G = ( gm_type::nNodes == element_type::numVertices ) ?__e.vertices() : __e.G();
    //M_G = __e.G();
    M_g.resize( M_G.size2(), PDim );
    //M_element_c = boost::shared_ptr<element_type const>(&__e);
    M_element = boost::addressof( __e );
    M_id = __e.id();
    M_e_marker = __e.marker();
    M_e_marker2 = __e.marker2();
    M_e_marker3 = __e.marker3();
    M_face_id = invalid_uint16_type_value;
    M_h = __e.h();
    M_h_min = __e.hMin();
    M_meas = __e.measure();
    M_xrefq = M_pc->nodes();

    FEELPP_ASSERT( M_G.size2() == M_gm->nbPoints() )( M_G.size2() )( M_gm->nbPoints() ).error( "invalid dimensions" );
    FEELPP_ASSERT( M_pc ).error( "invalid precompute data structure" );

    if ( vm::has_point<context>::value )
    {
        std::fill( M_xrealq.data().begin(), M_xrealq.data().end(), value_type( 0 ) );
        const uint16_type size1 = M_G.size1();
        const uint16_type size3 = M_G.size2();
        const uint16_type size2 = M_pc->nPoints();

        for ( uint16_type i = 0; i < size1; ++i )
            for ( uint16_type j = 0; j < size2; ++j )
            {
                for ( uint16_type k = 0; k < size3; ++k )
                    M_xrealq( i, j ) += M_G( i, k ) * M_pc->phi()[k][j]( 0,0 );
            }
    }

    if ( vm::has_jacobian<context>::value )
    {
        updateJKBN( mpl::bool_<is_linear>() );
    }

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

/**
 *
 */
uint16_type nPoints() const
{
    return M_npoints;
}

/**
 * \return the set of points in the reference convex
 */
matrix_node_t_type const& xRefs() const
{
    return M_xrefq;
}

/**
 * \return the q-th point in the reference convex
 */
ublas::matrix_column<matrix_node_t_type const> xRef( int q ) const
{
    return ublas::column( M_xrefq, q );
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
#if 0
value_type J( int q ) const
{
    //BOOST_STATIC_ASSERT( vm::has_jacobian<context>::value );
    //return J( q, mpl::bool_<is_linear>() );
    if ( is_linear )
        return M_J;

    else
        return M_Jt[q];
}
#else
value_type J( int q ) const
{
    if ( is_linear )
        return M_J;

    else
        return M_Jt[q];
}
#endif
#if 0
/**
 * \internal
 */
value_type J( int /*q*/, mpl::bool_<true> ) const
{
    return M_J;
}

/**
 * \internal
 */
value_type J( int q, mpl::bool_<false> ) const
{
    return M_Jt[q];
}
#endif

/**
 * \return the matrix associated with the geometric nodes
 */
//matrix_node_t_type const& G() const { return M_G; }
matrix_type const& G() const
{
    return M_G;
}

/**
 * Get the inverse of the transformation at the \p i -th point
 *
 * \param i the index of point where the pseudo-inverse is requested
 * \return the pseudo inverse of the transformation
 */
matrix_type const& B( int i ) const
{
    //BOOST_STATIC_ASSERT( vm::has_kb<context>::value );
    return B( i, mpl::bool_<is_linear>() );
}
value_type B( int c1, int c2, int q ) const
{
    //BOOST_STATIC_ASSERT( vm::has_kb<context>::value );
    return B( q, mpl::bool_<is_linear>() )( c1, c2 );
}

matrix_type const& projectorTangent( int i ) const
{
    return M_Ptangent;
}
value_type projectorTangent( int c1, int c2, int q ) const
{
    return M_Ptangent(c1,c2);
}

matrix_type const& K( int i ) const
{
    return M_K;
}
value_type K( int c1, int c2, int q ) const
{
    return M_K( c1, c2 );
}
/**
 * \internal
 */
matrix_type const& B( int /*i*/, mpl::bool_<true> ) const
{
    return M_B;
}
/**
 * \internal
 */
matrix_type const& B( int i, mpl::bool_<false> ) const
{
    return M_Bt[i];
}

/**
 * the tensor of rank 4 for the transformation of 2nd
 * order derivatives from the reference element to the real
 * element. The tensor has the shape \f$[N,N,P,P]\f$.
 *
 * \return the tensor of rank 4
 */
boost::multi_array<value_type,4> const& B3() const
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
node_t_type barycenterReal() const
{
    node_t_type __barycenter( M_G.size1() );
    __barycenter.clear();


    for ( uint16_type __c = 0; __c < M_G.size2(); ++__c )
    {
        __barycenter += ublas::column( M_G, __c );
    }

    __barycenter /= M_G.size2();
    return __barycenter;
}


/**
 * tell whether the point is on the surface of the convex
 * @return true if the point is on the surface, false otherwise
 */
bool isOnConvexSurface() const
{
    if ( trans == fem::LINEAR )
    {
        // x -x0 - K(0)\bar{x}
        return std::abs( ublas::norm_2( xReal()-ublas::column( M_G, 0 )-
                                        ublas::prod( M_K, xRef() ) ) ) < 1e-10;
    }

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
    if ( is_linear )
        return M_n_norm;

    else
        return M_n_normq[q];
}

/**
 * get the normal of the real element
 *
 * @return the normal of the real element
 */
node_t_type const& normal() const
{
    //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
    return M_n_real;
}

/**
 * normal getter
 * @param q the index of the normal to be returned
 * @return the \p q -th normal
 */
ublas::matrix_column<matrix_node_t_type const> normal( int q ) const
{
    //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
    if ( is_linear )
        return M_n_real;

    else
    {
        return ublas::column( M_n_realq
                     , q );
    }
}

/**
 * get the unit normal of the real element
 *
 * @return the unit normal of the real element
 */
node_t_type const& unitNormal() const
{
    //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
    return M_un_real;
}

//ublas::matrix_column<matrix_node_t_type const> unitNormal( int q ) const
// node_t_type const& unitNormal( int q ) const
node_t_type unitNormal( int q ) const
{
    //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
    if ( is_linear )
        return M_un_real;

    else
        return ublas::column( M_un_realq, q );
}


value_type const& unitNormal( int n, int q ) const
{
    //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );

    if ( is_linear )
        return M_un_real( n );

    else
        return M_un_realq( n, q );
}
/**
 * @brief the scaled tangent to the current face
 * @return the scaled tangent
 */
node_t_type const& tangent() const
{
    DCHECK( nDim == 2 ) <<  "Tangent is available only in 2D";
    DCHECK( is_linear ) << "Invalid call to unitTangent, the geometric mapping is not linear";
    return M_t_real;
}

/**
 * @brief scaled tangent at point

 * @param q index of the point at which the tangent is evaluated
 * @return the scaled tamgent at the point
 */
node_t_type const& tangent( int q ) const
{
    DCHECK( nDim == 2 ) <<  "Tangent is available only in 2D";
    if ( is_linear )
        return M_t_real;
    return M_t_realq(q);
}
/**
 * @brief the unit Tangent (linear case)
 *
 * @return the unit tangent in the linear case
 */
node_t_type const& unitTangent() const
{
    DCHECK( nDim == 2 ) <<  "Tangent is available only in 2D";
    DCHECK( is_linear ) << "Invalid call to unitTangent, the geometric mapping is not linear";
    return M_ut_real;
}

/**
 * @brief the unit tangent a point \p q
 *
 * @param q point index
 * @return unit tangent
 */
node_t_type const& unitTangent( int q ) const
{
    DCHECK( nDim == 2 ) <<  "Tangent is available only in 2D";
    DCHECK( is_linear ) << "Invalid call to unitTangent, the geometric mapping is not linear";
    return M_ut_real;
}

/**
 * @brief unit tangent coordinate \p n at point \p q
 *
 * @param n coordinate
 * @param q point index
 *
 * @return [description]
 */
value_type const& unitTangent( int n, int q ) const
{
    DCHECK( nDim == 2 ) <<  "Tangent is available only in 2D";
    if ( is_linear )
        return M_ut_real( n );
    return M_ut_realq( n, q );
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
Marker1 marker() const
{
    return M_e_marker;
}

/**
 * get the marker2 of the element
 *
 * @return the marker2 of the element
 */
Marker2 marker2() const
{
    return M_e_marker2;
}

/**
 * get the marker3 of the element
 *
 * @return the marker3 of the element
 */
Marker2 marker3() const
{
    return M_e_marker3;
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
    return M_gm->radiusEstimate( M_G );
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
    return permutation( mpl::bool_<( nDim>=2 )>() );
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
std::vector<std::map<permutation_type, precompute_ptrtype> > const & pcFaces() const
{
    return M_pc_faces;
}
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

void setPcFaces( std::vector<std::map<permutation_type, precompute_ptrtype> > const& __pcfaces )
{
    M_pc_faces = __pcfaces;
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

void updateJKBN( mpl::bool_<true>, mpl::true_  )
{
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> P( M_G.data().begin(), M_G.size1(), M_G.size2() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> GradPhi( M_g_linear.data().begin(), M_g_linear.size1(), M_g_linear.size2() );
        Eigen::Map<Eigen::Matrix<value_type,NDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MK( M_K.data().begin(), M_K.size1(), M_K.size2() );
        Eigen::Map<Eigen::Matrix<value_type,PDim,PDim,((NDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MCS( M_CS.data().begin(), M_CS.size1(), M_CS.size2() );
        Eigen::Map<Eigen::Matrix<value_type,PDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MCSi( M_CSi.data().begin(), M_CSi.size1(), M_CSi.size2() );
        Eigen::Map<Eigen::Matrix<value_type,NDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MB( M_B.data().begin(), M_B.size1(), M_B.size2() );
        MK.noalias() =  P*GradPhi;
        M_J = math::abs( MK.determinant() );
        MCS = MK.inverse();
        MB.noalias() = MCS.transpose();

}
void updateJKBN( mpl::bool_<true>, mpl::false_  )
{
#if 0
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> P( M_G.data().begin(), M_G.size1(), M_G.size2() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> GradPhi( M_g_linear.data().begin(), M_g_linear.size1(), M_g_linear.size2() );
        Eigen::Map<Eigen::Matrix<value_type,NDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MK( M_K.data().begin(), M_K.size1(), M_K.size2() );
        Eigen::Map<Eigen::Matrix<value_type,PDim,PDim,((NDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MCS( M_CS.data().begin(), M_CS.size1(), M_CS.size2() );
        Eigen::Map<Eigen::Matrix<value_type,PDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MCSi( M_CSi.data().begin(), M_CSi.size1(), M_CSi.size2() );
        Eigen::Map<Eigen::Matrix<value_type,NDim,PDim,((PDim==1)?Eigen::ColMajor:Eigen::RowMajor)>> MB( M_B.data().begin(), M_B.size1(), M_B.size2() );
        MK.noalias() =  P*GradPhi;
        //ublas::axpy_prod( M_G, M_g_linear, M_K, true );

        // CS = K^T K
        MCSi.noalias() = MK.transpose()*MK;
        M_J = math::sqrt( math::abs( MCSi.determinant() ) );
        //if ( vm::has_kb<context>::value )
        {
                MCS=MCSi.inverse();
                // B = K CS
                MB.noalias() = MK*MCS;
        }
#else
        ublas::axpy_prod( M_G, M_g_linear, M_K, true );
        ublas::noalias( M_CS ) = ublas::prod( ublas::trans( M_K ), M_K );
        M_J = math::sqrt( math::abs( det<PDim>( M_CS ) ) );
        inverse<PDim>( M_CS, M_CSi );
        ublas::axpy_prod( M_K, M_CSi, M_B, true );

#endif

}
/**
 * update Jacobian data : linear case
 */
void updateJKBN( mpl::bool_<true>  )
{
    if ( !M_gm->isCached() ||
            ( M_gm->isCached() && M_gm->cached( M_id ) == false ) )
    {
        updateJKBN( mpl::true_(), mpl::bool_<NDim==PDim>() );
        if ( M_gm->isCached() )
        {
            // cache J, K and B
            M_gm->addJ( M_id, M_J );
            //if ( vm::has_kb<context>::value )
            {
                M_gm->addK( M_id, M_K );
                M_gm->addB( M_id, M_B );
            }
            M_gm->setCached( M_id, true );

            //LOG(INFO) << "(add to cache) J[" << M_id << "]=" <<  M_J << "\n";
            //LOG(INFO) << "(add to cache) B[" << M_id << "]=" <<  M_B << "\n";
        }
    }

    else
    {
        M_J = M_gm->J( M_id );
        //LOG(INFO) << "(use cache) J[" << M_id << "]=" <<  M_J << "\n";
        //if ( vm::has_kb<context>::value )
        {
            M_K = M_gm->K( M_id );
            M_B = M_gm->B( M_id );
            //LOG(INFO) << "(use cache) B[" << M_id << "]=" <<  M_B << "\n";
        }


    }

    if ( vm::has_hessian<context>::value || vm::has_laplacian<context>::value )
    {

        for ( uint16_type k = 0; k < NDim; ++k )
            for ( uint16_type l = 0; l < NDim; ++l )
                for ( uint16_type i = 0; i < PDim; ++i )
                    for ( uint16_type j = 0; j < PDim; ++j )
                        M_B3[k][l][i][j] = M_B( k,i )*M_B( l,j );

#if 0
        //B3(k + N_*l, i + P*j) = BB(k, i) * BB(l, j);
        std::cout << "element " << this->id() << " B3 = " << "\n";

        for ( int i = 0; i < nDim; ++i )
            for ( int j = 0; j < nDim; ++j )
                for ( typename boost::multi_array<value_type,4>::index k = 0; k < nDim; ++k )
                    for ( typename boost::multi_array<value_type,4>::index l = 0; l < nDim; ++l )
                        std::cout << "B3[" << i << "][" << j << "][" << k << "][" << l << "]="
                                  << this->B3()[i][j][k][l] << "\n";

#endif
    }

    if ( ( ( NDim != PDim ) || ( vm::has_normal<context>::value ) ) && ( M_face_id != invalid_uint16_type_value ) )
    {
        ublas::axpy_prod( M_B,
                          M_gm->referenceConvex().normal( M_face_id ),
                          M_n_real,
                          true );
        M_n_norm = ublas::norm_2( M_n_real );
        M_un_real = M_n_real/M_n_norm;

    }

    if( vm::has_tangent<context>::value  && ( M_face_id != invalid_uint16_type_value ) )
    {
        if ( NDim == 2 )
        {
            // t = |\hat{e}|*o_K*(K*t_ref)/|e| where o_K is the sign(e*x_K(\hat{e}))
            ublas::axpy_prod( M_K,
                              M_gm->referenceConvex().tangent( M_face_id ),
                              M_t_real,
                              true );
            double ratio = M_gm->referenceConvex().h( M_face_id )/M_h_face;

            M_t_norm = ublas::norm_2( M_t_real );
            M_ut_real = M_t_real/M_t_norm;
        }

        // compute projector on tangent plane
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > P(M_Ptangent.data().begin(), NDim, NDim);
        P.setIdentity(NDim, NDim);
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > N(M_un_real.data().begin(), NDim,1);
        P.noalias() -= N*N.transpose();
    }


}

/**
 * update Jacobian data : nonlinear case
 */
void updateJKBN( mpl::bool_<false>  )
{

     const bool doComputeNormal = ( ( NDim != PDim ) || ( vm::has_normal<context>::value ) ) && ( M_face_id != invalid_uint16_type_value );

    //std::cout << "nPoints() =" << nPoints() << "\n";
    //VLOG(1) << "[geomap] G = "<< M_G << "\n";
    //double res = 0;
    for ( int q = 0; q < nPoints(); ++q )
    {
        //std::cout << "q =" << q << "\n";
        M_gm->gradient( q, M_g, M_pc.get() );
        //VLOG(1) << "[geomap] g[" << q << "] = "<< M_g << "\n";

#if 0
        blas::gemm( M_G, M_g, M_K );
#else
        ublas::axpy_prod( M_G, M_g, M_K, true );
#endif

        if ( NDim == PDim )
        {
            M_J = math::abs( det<NDim>( M_K ) );

            if ( vm::has_kb<context>::value )
            {
                inverse<NDim>( M_K, M_CS );
                //ublas::noalias(M_B) = ublas::trans( M_CS );
                M_B = ublas::trans( M_CS );
                //std::cout << "========== B[" << q << "]=" << M_B << "\n";
            }
        }

        else // N != P
        {
            // CS = K^T K
#if 1
            ublas::prod( ublas::trans( M_K ), M_K, M_CS );
#else
            blas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                        1.0, M_K, M_K,
                        0.0, M_CS );
#endif
            M_J = math::sqrt( math::abs( det<PDim>( M_CS ) ) );

            if ( vm::has_kb<context>::value )
            {
                inverse<PDim>( M_CS, M_CSi );
                // B = K CS
#if 1
                ublas::axpy_prod( M_K, M_CSi, M_B );
#else
                blas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                            1.0, M_K, M_CSi,
                            0.0, M_B );
#endif
            }

        }

        //VLOG(1) << "[geomap] J[" << q << "]= "<< M_J << "\n";
        //res += M_J;
        // store q-th jacobian entry

        M_Jt[q] = M_J;

        if ( vm::has_kb<context>::value || doComputeNormal )
        {
            //VLOG(1) << "[geomap] B[" << q << "]= "<< M_B << "\n";
            M_Bt[q].resize( M_B.size1(), M_B.size2() );
            M_Bt[q] = M_B;
        }


    }


    //VLOG(1) << "[geomap] res(sum J) = " << res << "\n";
    if ( doComputeNormal )
    {
        //std::cout << "has normal\n";
        for ( int q = 0; q < nPoints(); ++q )
        {
#if 0
            blas::gemv( traits::NO_TRANSPOSE,
                        1.0, M_Bt[ q ], M_n[M_face_id],
                        0.0, M_n_real );
#else

            if ( 0 ) //trans == fem::LINEAR )
            {
                ublas::axpy_prod( M_B,
                                  M_gm->referenceConvex().normal( M_face_id ),
                                  M_n_real,
                                  true );
            }

            else
                ublas::axpy_prod( M_Bt[ q ],
                                  M_gm->referenceConvex().normal( M_face_id ),
                                  M_n_real,
                                  true );

#endif
            //std::cout << "[geomap] point " << q << " n_real = " << M_n_real << "\n";
            M_n_norm = ublas::norm_2( M_n_real );
            M_un_real = M_n_real/M_n_norm;
            ublas::column( M_n_realq
                 , q ) = M_n_real;
            ublas::column( M_un_realq, q ) = M_un_real;
            M_n_normq[q] = M_n_norm;

            if ( NDim != PDim )
                M_Jt[q] *= M_n_norm;

            if ( vm::has_tangent<context>::value )
            {

            }
        }

#if 0
        std::cout << "[geomap] face id = " << M_face_id << "\n"
                  << "[geomap] ref normal = " << M_gm->referenceConvex().normal( M_face_id ) << "\n"
                  << "[geomap] M_n_real = " << M_n_realq
          << "\n"
                  << "[geomap] M_un_realq = " << M_un_realq << "\n"
                  << "[geomap] M_n_normq = " << M_n_normq << "\n";
#endif
    }
}

Context();

private:

gm_ptrtype M_gm;

element_type const* M_element;
//boost::shared_ptr<element_type const> M_element_c;
//element_type M_element_c;

precompute_ptrtype M_pc;
std::vector<std::map<permutation_type, precompute_ptrtype> > M_pc_faces;
uint16_type M_npoints;

value_type M_J;

matrix_type  M_G;
ublas::vector<node_t_type> M_n;
node_t_type M_n_real;
node_t_type M_un_real;
value_type M_n_norm;
node_t_type M_t_real;
node_t_type M_ut_real;
value_type M_t_norm;

matrix_node_t_type  M_xrefq;
matrix_type M_xrealq;

// normals
matrix_node_t_type M_n_realq;
matrix_node_t_type M_un_realq;
ublas::vector<value_type> M_n_normq;
// tangents
matrix_node_t_type M_t_realq;
matrix_node_t_type M_ut_realq;
ublas::vector<value_type> M_t_normq;

matrix_type M_g_linear;
matrix_type M_g;
matrix_type M_K;
matrix_type M_CS;
matrix_type M_CSi;
matrix_type M_B;
matrix_type M_Ptangent;

boost::multi_array<value_type,4> M_B3;

size_type M_id;
Marker1 M_e_marker;
Marker2 M_e_marker2;
Marker3 M_e_marker3;
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

vector_type M_Jt;
std::vector<matrix_type> M_Bt;

permutation_type M_perm;
}; // Context




template<size_type context_v, typename ElementType>
boost::shared_ptr<Context<context_v,ElementType> >
context( geometric_mapping_ptrtype gm, ElementType const& e, precompute_ptrtype const& pc )
{
    return boost::shared_ptr<Context<context_v,ElementType> >(
               new Context<context_v, ElementType>( gm,
                       e,
                       pc ) );
}

template<size_type context_v, typename ElementType>
boost::shared_ptr<Context<context_v,ElementType> >
context( ElementType const& e, precompute_ptrtype const& pc )
{
    return boost::shared_ptr<Context<context_v,ElementType> >(
                new Context<context_v, ElementType>
                (
                 //super_enable_this::shared_from_this(),
                 boost::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP > >(this->shared_from_this()),
                e,
                pc ) );
}

template<size_type context_v, typename ElementType>
boost::shared_ptr<Context<context_v,ElementType> >
context( geometric_mapping_ptrtype gm,
         ElementType const& e,
         std::vector<std::map<typename ElementType::permutation_type,precompute_ptrtype> > & pc,
         uint16_type f )
{
    return boost::shared_ptr<Context<context_v,ElementType> >(
               new Context<context_v, ElementType>( gm,
                       e,
                       pc,
                       f ) );
}
template<size_type context_v, typename ElementType>
boost::shared_ptr<Context<context_v,ElementType> >
context( ElementType const& e,
         std::vector<std::map<typename ElementType::permutation_type,precompute_ptrtype> > & pc,
         uint16_type f )
{
    return boost::shared_ptr<Context<context_v,ElementType> >(
               new Context<context_v, ElementType>
               (
                //super_enable_this::shared_from_this(),
                boost::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP > >(this->shared_from_this()),
               e,
               pc,
               f ) );
}


/**
 * Inverse of the geometric mapping for a given context
 */
class Inverse
{
    public:

    typedef GeoMap geometric_mapping_type;
    typedef boost::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

    static const fem::transformation_type trans = geometric_mapping_type::trans;

    typedef typename geometric_mapping_type::node_t_type node_t_type;
    typedef typename geometric_mapping_type::node_t_type node_type;
    typedef typename geometric_mapping_type::matrix_node_t_type matrix_node_t_type;
    typedef typename geometric_mapping_type::matrix_node_t_type matrix_node_type;

    typedef ublas::vector<double> dense_vector_type;
    typedef ublas::matrix<double> dense_matrix_type;

    template<typename GeoElem>
    Inverse( geometric_mapping_ptrtype __gm, GeoElem const& __ge,
             WorldComm const& worldComm = Environment::worldCommSeq() )
        :
        M_gm( __gm ),
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
        M_nlsolver( SolverNonLinear<double>::build( "petsc","", worldComm ) )
#else
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_GMM, worldComm ) )
#endif
{
    this->init();
}

    template<typename GeoElem>
    Inverse( geometric_mapping_ptrtype __gm, GeoElem const& __ge, mpl::int_<1>/**/ ,
             WorldComm const& worldComm = Environment::worldCommSeq() )
    :
        M_gm( __gm ),
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
        M_nlsolver( SolverNonLinear<double>::build( "petsc","", worldComm ) )
#else
        M_nlsolver( SolverNonLinear<double>::build( SOLVERS_GMM,worldComm ) )
#endif
{
    this->init();
}


private :
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
            M_nlsolver->setRelativeResidualTol( 1e-12/*1e-16*/ );

            M_J.reset( new dense_matrix_type( M_xref.size(),M_xref.size() ) );
            M_R.reset( new dense_vector_type( M_xref.size() ) );
#endif
        }
}

public :
    template<typename GeoElem>
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
    template<typename GeoElem>
    void
    update( GeoElem const& __ge, mpl::int_<1>/**/ )
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
        return std::abs( ublas::norm_2( xReal()-ublas::column( M_G, 0 )-
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

matrix_node_t_type operator() ( matrix_node_t_type const& real_pts, bool allow_extrapolation = false ) const
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
        :
        gmi( __gmi ),
        xreal( __xr )
{}
// evaluate
scalar_type operator()( const node_t_type& x ) const
{
    node_type r = gmi.geometricMapping()->transform( x, gmi.G() ) - xreal;
    return ublas::inner_prod( r, r )/2.;
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
        r = y - M_xreal ;

    else
    {
        M_gm->gradient( x, M_g );
        ublas::axpy_prod( M_G, M_g, M_K );
        ublas::prod( ublas::trans( M_K ), y-M_xreal, r );
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
    DVLOG(2) << "[update] g = " << M_g << "\n";

    checkInvariant();

    ublas::axpy_prod( M_G, M_g, M_K );
    DVLOG(2) << "[update] K(0) = " << M_K << "\n";

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

    DVLOG(2) << "[update] B(0) = " << M_B << "\n";
}

void
checkInvariant() const
{
    FEELPP_ASSERT( M_G.size2() == M_g.size1() )( M_G.size2() )( M_g.size1() ).error( "G,g invalid dimensions" );
    FEELPP_ASSERT( M_G.size1() == M_K.size1() )( M_G.size1() )( M_K.size1() ).error( "G,K invalid dimensions" );
    FEELPP_ASSERT( M_g.size2() == M_K.size2() )( M_g.size2() )( M_K.size2() ).error( "g,K invalid dimensions" );
    FEELPP_ASSERT( M_B.size2() == M_gm->dim() )( M_B.size1() )( N() ).error( "B,gm invalid dimensions" );
}


bool linearInverse()
{
    checkInvariant();

    size_type N = M_xreal.size();
    size_type P = M_xref.size();

    node_type y( M_xreal );

    DVLOG(2) << "y = xreal = " << y << "\n";
    //DVLOG(2) << "G(0)  = " << node_type( M_x0 << "\n";
    y.minus_assign(  ublas::column( M_G, 0 ) );
    DVLOG(2) << "y - G(0) = " << y << "\n";

    DVLOG(2) << "B(0) = " << M_B << "\n";
    DVLOG(2) << "xref = " << ublas::prod( ublas::trans( M_B ), y ) << "\n";

    // xref = B^T * y = B^T * ( x_real - x_0)
    M_xref.assign( ublas::prod( ublas::trans( M_B ), y )-ublas::scalar_vector<value_type>( P, 1.0 ) );

    DVLOG(2) << "[GeoMap::Inverse::linearInverse] xref : " << M_xref << "\n";

    bool __isin;
    double vmin;
    boost::tie( __isin, vmin ) = M_gm->isIn( M_xref );
    DVLOG(2) << "[GeoMap::Inverse::linearInverse] isIn : " << __isin << "\n";

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
                        real_pts-ublas::outer_prod( ublas::column( M_G, 0 ),
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

    if ( __isin  &&
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

    while ( res > EPS/10 && cnt > 0 )
    {
        node_type xn( P );

        // xn = B^T rn
        ublas::axpy_prod( rn, M_B, xn );

        scalar_type newres;

        for ( int16_type i=1; i<=256; i*=2 )
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

            if ( newres < 1.5*res )
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

    if ( __isin  &&
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

boost::shared_ptr<SolverNonLinear<double> > M_nlsolver;
boost::shared_ptr<dense_matrix_type> M_J;
boost::shared_ptr<dense_vector_type> M_R;

}; // Inverse

private:

element_gm_ptrtype _elementMap;
face_gm_ptrtype _boundaryMap;
matrix_type M_g_linear;

friend class Inverse;

reference_convex_type M_refconvex;

       };

#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/facilities/identity.hpp>

template<int Dim, int Order,int RealDim,  template<uint16_type,uint16_type,uint16_type> class Entity = Simplex, typename T = double>
struct GT_Lagrange
{};

template<int Dim, int Order, int RealDim, typename T = double>
struct GT_QK
{};

# /* List of dims. */
# define FEELPP_GEOMAP                                    \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           1,                           \
                           (                            \
                            Lagrange                    \
                                                    )   \
                                                    )   \
    /**/
# /* List of dims. */
# define FEELPP_DIMS                                      \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           3,                           \
                           (                            \
                            1,2,3                       \
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
                           6,                                   \
                           (                                    \
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
# define FEELPP_ENTITY BOOST_PP_TUPLE_TO_LIST( 2, ( Simplex,Hypercube ) )
/**/
#
# define FEELPP_GEN_GT(GEOM,LDIM,LORDER)         \
    "GT_" #GEOM "(" #LDIM "," #LORDER ")"       \
    /**/
#
# /* Generates code for all dim and order. */
# define FEELPP_GT_FACTORY_OP(_, GDO)            \
    FEELPP_GT_FACTORY GDO                        \
    /**/
#
#
# define FEELPP_GT_DIM(T)  BOOST_PP_TUPLE_ELEM(2, 0, T) \
/**/
# define FEELPP_GT_REALDIM(T)  BOOST_PP_TUPLE_ELEM(2, 1, T)   \
/**/
#
#
#define FEELPP_GT_FACTORY(GEOM,LDIMS,LORDER,ENTITY)                       \
    template<typename T>                                                \
    struct BOOST_PP_CAT(GT_, GEOM)<FEELPP_GT_DIM(LDIMS), LORDER, FEELPP_GT_REALDIM(LDIMS), ENTITY, T> \
        :                                                               \
        public GeoMap<FEELPP_GT_DIM(LDIMS), LORDER, FEELPP_GT_REALDIM(LDIMS), T, ENTITY, GEOM> \
    {                                                                   \
        static const uint16_type nDim = FEELPP_GT_DIM(LDIMS);             \
        static const uint16_type order = LORDER;                        \
        static const uint16_type nRealDim = FEELPP_GT_REALDIM(LDIMS);     \
        static const uint16_type nRealDimCheck2d = mpl::if_< mpl::less_equal<mpl::int_<2>,mpl::int_<nRealDim> >, \
                                                             mpl::int_<nRealDim>, \
                                                             mpl::int_<nDim> >::type::value; \
                                                                        \
        typedef mpl::vector<boost::none_t,                              \
                            boost::none_t,                              \
                            GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>, \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM> > geomap_faces_t; \
        typedef mpl::vector<GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>, \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM>, \
                            GeoMap<3, LORDER, 3, T, ENTITY, GEOM>, \
                            boost::none_t> geomap_elements_t;           \
        typedef typename type_traits<T>::value_type value_type;         \
        typedef GeoMap<nDim, LORDER, nRealDim, T, ENTITY, GEOM> super;  \
        typedef BOOST_PP_CAT(GT_,GEOM)<nDim-1, LORDER, nRealDim,ENTITY, T> face_geo_type; \
                                                                        \
        static const uint16_type nDof = super::nDof;                    \
        static const uint16_type nNodes = super::nNodes;                \
        typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim> >::type element_gm_type; \
        typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim> >::type face_gm_type; \
        template<int N>                                                 \
            struct face_gm                                              \
        {                                                               \
            typedef typename mpl::at<geomap_faces_t, mpl::int_<N> >::type type; \
        };                                                              \
                                                                        \
        BOOST_PP_CAT(GT_,GEOM)()                                        \
            :                                                           \
            super( boost::shared_ptr<element_gm_type>(new element_gm_type()), boost::shared_ptr<face_gm_type>(new face_gm_type() )) \
            {}                                                          \
    };                                                                  \
    /**/
#if 0
#define FEELPP_GT_FACTORY(GEOM,LDIM,LORDER,LREALDIM,ENTITY)               \
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


template<typename Elem, template<uint16_type,uint16_type,uint16_type> class Entity = Simplex, typename T = double>
class RealToReference
{
public:
    static const uint16_type nDim = Elem::nDim;
    static const uint16_type nRealDim = Elem::nRealDim;

    typedef T value_type;
    typedef typename matrix_node<T>::type points_type;
    typedef GT_Lagrange<nDim,1, nRealDim, Entity, value_type> gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::Inverse inverse_gm_type;

    RealToReference( Elem const& elem )
        :
        M_gm( new gm_type ),
        M_igm( M_gm, elem )
    {}

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
