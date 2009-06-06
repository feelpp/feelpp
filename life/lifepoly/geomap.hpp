/*
  This file is part of the Life library

  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007,2008 Université Joseph Fourier Grenoble 1

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

#include <boost/mpl/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/signal.hpp>

//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/atlas/cblas3.hpp>

#include <life/lifecore/life.hpp>

#include <life/lifealg/svd.hpp>
#include <life/lifealg/solvernonlinear.hpp>
#include <life/lifealg/lu.hpp>
#include <life/lifealg/iteration.hpp>

#include <life/lifepoly/context.hpp>
#include <life/lifepoly/lagrange.hpp>
#include <life/lifemesh/marker.hpp>




namespace Life
{
//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;

struct GeomapInverse
{
    //GeomapInverse
};

template<uint16_type O,
         template<uint16_type Dim> class PolySetType,
         template<class, uint16_type, class> class Pts> class Lagrange;


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
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Entity = Simplex,
         template<uint16_type, template<uint16_type Dim> class PolySetType,
                  template<class, uint16_type, class> class Pts> class PP = Lagrange>
class GeoMap
    :
    public PP<Order,Scalar, PointSetEquiSpaced>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type
{
    typedef typename PP<Order, Scalar, PointSetEquiSpaced>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type super;

    typedef mpl::vector<boost::none_t, boost::none_t, GeoMap<1, Order, T, Entity, PP>, GeoMap<2, Order, T, Entity, PP> > geomap_faces_t;
    typedef mpl::vector<GeoMap<1, Order, T, Entity, PP>, GeoMap<2, Order, T, Entity, PP>, GeoMap<3, Order, T, Entity, PP>, boost::none_t> geomap_elements_t;
public:

    static const uint16_type nDim = super::nDim;
    static const uint16_type nDof = super::nDof;
    static const uint16_type nNodes = super::nNodes;
    static const fem::transformation_type trans = super::trans;

    typedef typename super::value_type value_type;

    typedef typename super::PreCompute precompute_type;
    typedef boost::shared_ptr<precompute_type> precompute_ptrtype;

    typedef typename super::convex_type convex_type;
    typedef Reference<convex_type,nDim,Order,nDim> reference_convex_type;

    typedef GeoMap<Dim, Order, T, Entity, PP > self_type;
    typedef self_type geometric_mapping_type;
    typedef boost::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

    typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim> >::type element_gm_type;
    typedef boost::shared_ptr<element_gm_type> element_gm_ptrtype;

    typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim> >::type face_gm_type;
    typedef boost::shared_ptr<face_gm_type> face_gm_ptrtype;

    template<int N>
    struct face_gm
    {
        typedef typename mpl::at<geomap_faces_t, mpl::int_<N> >::type type;
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
        _M_g_linear( nNodes, nDim ),
        M_refconvex()
    {
        LIFE_CONSTRUCTOR_END( 5046, "GeoMap" );

        if ( trans == fem::LINEAR )
            {
                //_M_g_linear.resize( nNodes, nDim );
                matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );

                ublas::vector<ublas::matrix<value_type> > m = super::derivate( __dummy_pts );

                LIFE_ASSERT( _M_g_linear.size2() == m.size() )( _M_g_linear.size2() )(  m.size() ).error( "invalid dimension" );
                LIFE_ASSERT( m( 0 ).size2() == 1 )( m( 0 ).size2() ).error( "Invalid number of points" );

                LIFE_ASSERT( _M_g_linear.size1() == m( 0 ).size1() )( _M_g_linear.size1() )( m( 0 ).size1() ).error( "invalid number of DOF" );

                //std::cout << "nNodes= " << nNodes << "\n"
                //<< "nDim= " << nDim << "\n";
                //std::cout << "_M_g_linear = " << _M_g_linear << "\n"
                //<< "m(0) = " << m( 0 ) << "\n";

                for ( uint16_type i = 0; i < nNodes; ++i)
                    {
                        for ( uint16_type n = 0; n < nDim; ++n)
                            {
                                //std::cout << "m(n)= " << m( n ) << "\n";
                                _M_g_linear( i, n) = m(n)( i, 0 );
                            }
                    }
#if 0
                for ( uint16_type i = 0; i < nNodes; ++i)
                    {
                        for ( uint16_type n = 0; n < nDim; ++n)
                            {
                                _M_g_linear( i, n ) = this->dPhi( i, n, __dummy_pt );
                            }
                    }
#endif // 0
            }
        LIFE_CONSTRUCTOR_END( 5046, "GeoMap" );
    }
    /** default constructor */
    GeoMap( element_gm_ptrtype const& e,  face_gm_ptrtype const& f )
        :
        super(),
        M_is_cached( false ),
        _elementMap( e ),
        _boundaryMap( f ),
        _M_g_linear( nNodes, nDim ),
        M_refconvex()
    {
        LIFE_CONSTRUCTOR_END( 5046, "GeoMap" );

        if ( trans == fem::LINEAR )
            {
                //_M_g_linear.resize( nNodes, nDim );
                node_t_type __dummy_pt( nDim );

                matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );

                //std::cout << "geomap::derivate<> pts=" << __dummy_pts << "\n";
                //std::cout << "geomap::derivate<> m=" << super::derivate( __dummy_pts ) << "\n";
                ublas::vector<ublas::matrix<value_type> > m = super::derivate( __dummy_pts );
                //std::cout << "nNodes= " << nNodes << "\n"
                //<< "nDim= " << nDim << "\n";
                //std::cout << "_M_g_linear = " << _M_g_linear << "\n"
                //<< "m(0) = " << m( 0 ) << "\n";

                for ( uint16_type i = 0; i < nNodes; ++i)
                    {
                        for ( uint16_type n = 0; n < nDim; ++n)
                            {
                                //std::cout << "m(n)= " << m( n ) << "\n";
                                _M_g_linear( i, n) = m(n)( i, 0 );
                            }
                    }
#if 0
                for ( uint16_type i = 0; i < nNodes; ++i)
                    {
                        for ( uint16_type n = 0; n < nDim; ++n)
                            {
                                _M_g_linear( i, n ) = this->dPhi( i, n, __dummy_pt );
                            }
                    }
#endif // 0
            }

        LIFE_CONSTRUCTOR_END( 5046, "GeoMap" );
    }
    /**
       destructor
    */
    ~GeoMap()
    {}

    /**
       \return the dimension of the underlying element
    */
    uint16_type dim() const { return nDim; }

    /**
       \return true if the geometric mapping is linear, false otherwise
    */
    bool isLinear() const { return trans == fem::LINEAR; }


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
    ublas::vector<value_type>  refNode( uint16_type i ) const { return ublas::column( this->points(), i ); }

    /**
     * @return the reference convex
     */
    reference_convex_type const& referenceConvex() const { return M_refconvex;}

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
                value_type __phi_at_pt = super::evaluate( __i, __ref_p )(0);
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
#if 0
        blas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                    1.0, G, pc->phi(),
                    0.0, x );
#else
        ublas::axpy_prod( G, pc->phi(), x, true );
#endif
    }


    /**
     *  compute the gradient of the transformation in the reference
     *  element
     *
     *  Compute the gradient at node \c x, pc is resized to
     *  [nbNodes() x dim()] if the transformation is linear, \c x is
     *  not used at all
     */
    void gradient(const node_t_type& __pt,
                  matrix_type& __g ) const
    {
        namespace lambda = boost::lambda;

        if ( trans == fem::LINEAR )
            {
                __g = _M_g_linear;
            }
        else
            {
                LIFE_ASSERT( __pt.size() == dim() )( __pt.size() )( dim() ).error( "invalid dimension" );

                matrix_node_t_type __pts( nDim, 1 );
                ublas::column( __pts, 0 ) = __pt;

                ublas::vector<ublas::matrix<value_type> > m = super::derivate( __pts );
                for ( uint16_type n = 0; n < nDim; ++n)
                    {
                        ublas::column( __g, n) = ublas::column( m(n), 0 );
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
        LIFE_ASSERT( __pc )( __idref ).error( "a PreCompute must be set first before using this function" );

        if ( trans == fem::LINEAR )
            {
                __g = _M_g_linear;
            }
        else
            {
                LIFE_ASSERT( __pc->dim() == dim() )( __pc->dim() )( dim() ).error( "invalid dimension" );
                for ( size_type i = 0; i < nNodes; ++i)
                    {
                        for ( uint16_type n = 0; n < nDim; ++n)
                            {
                                __g( i, n ) = __pc->grad(i,  0, n, __idref );
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

        LIFE_ASSERT(  __g.size1() == G.size2() )( __g.size1() )( G.size2() ).error( "invalid sizes" );

        matrix_node_t_type K( G.size1(), __g.size2() );
        ublas::axpy_prod( G, __g, K );

        SVD<matrix_node_t_type> __svd( K );
        value_type __max;
        value_type __min;
        __svd.conditionNumber( __max, __min );
        return __max*sqrt(value_type(N)) / value_type(N);
    }


    bool M_is_cached;
    bool isCached() const { return M_is_cached; }

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
        Log() << "[Geomap] start caching J,K,B for " << nelts << " elements\n";
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
                Log() << "invalid cache\n";
                return false;
            }
        return true;
    }
    bool cached( int e ) const {  LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); return M_cached[e]; }
    void setCached( int e, bool v ) { LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); M_cached[e] = v; }
    double J( int e ) const { LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); return M_J[e]; }
    matrix_type const& B( int e ) const { LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); return M_B[e]; }
    matrix_type const& K( int e ) const { LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); return M_K[e]; }
    void addJ( int e, double v )
    {
        LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
        M_J[e] = v;
    }
    void addK( int e, matrix_type const& K )
    {
        LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" ); M_K[e].resize( K.size1(), K.size2() );
        for( size_type i = 0; i < K.size1(); ++i )
            for( size_type j = 0; j < K.size2(); ++j )
                M_K[e]( i, j ) = K( i, j );
    }
    void addB( int e, matrix_type const& B )
    {
        LIFE_ASSERT( this->isCacheValid() )( e ).error( "invalid cache" );
        M_B[e].resize( B.size1(), B.size2() );
        for( size_type i = 0; i < B.size1(); ++i )
            for( size_type j = 0; j < B.size2(); ++j )
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

        static const size_type context = context_v;
        // reference space dimension
        static const uint16_type PDim = ElementType::nDim;
        // real space dimension
        static const uint16_type NDim = ElementType::nRealDim;
        static const uint16_type nDim = NDim;
        // type of transformation (linear or not)
        static const fem::transformation_type trans = geometric_mapping_type::trans;
        static const bool is_linear = (trans == fem::LINEAR);

        static const bool condition = ((PDim==NDim)||((NDim>=1)&&(PDim==NDim-1)));
        BOOST_MPL_ASSERT_MSG( condition, INVALID_DIM, (mpl::int_<NDim>, mpl::int_<PDim>, ElementType ) );
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim> >,
                                  mpl::identity<GeoMap<Dim, Order, T, Entity, PP > >,
                                  typename mpl::if_<mpl::and_<mpl::greater_equal<mpl::int_<NDim>, mpl::int_<1> >,
                                                              mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> > >,
                                                    //typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
                                                    // mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
                                                    mpl::identity<typename GeoMap<NDim, Order, T, Entity, PP >::face_gm_type>,
                                                    mpl::identity<boost::none_t> >::type>::type::type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;

    typedef typename gm_type::value_type value_type;

    typedef typename gm_type::precompute_ptrtype precompute_ptrtype;

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
             precompute_ptrtype const& __pc,
             uint16_type __f = invalid_uint16_type_value )
        :
        _M_gm( __gm ),
        _M_element( __e ),
        _M_pc( __pc ),
        _M_npoints( _M_pc.get()->nPoints() ),

        //_M_xref( PDim ),
        //_M_xreal( NDim ),
        //_M_x0( NDim ),
        _M_J( 0 ),
        _M_G( __e.G() ),
        _M_n( _M_gm->referenceConvex().normals() ),
        _M_n_real( NDim ),
        _M_u_n_real( NDim ),
        _M_n_norm( 0 ),
        _M_xrefq( PDim, nPoints() ),
        _M_xrealq( NDim, nPoints() ),
        _M_nrealq( NDim, nPoints() ),
        _M_unrealq( NDim, nPoints() ),
        _M_nnormq( nPoints() ),

        _M_g( _M_G.size2(), PDim ),
        _M_K( NDim, PDim ),
        _M_CS( PDim, PDim ),
        _M_CSi( PDim, PDim ),
        _M_B( NDim, PDim ),
        _M_B3( boost::extents[NDim][NDim][PDim][PDim] ),
        _M_id( __e.id() ),
        _M_e_marker( __e.marker() ),
        _M_e_marker2( __e.marker2() ),
        _M_e_marker3( __e.marker3() ),
        _M_elem_id_1( invalid_size_type_value ),// __e.ad_first() ),
        _M_pos_in_elem_id_1( invalid_uint16_type_value ),  //__e.pos_first() ),
        _M_elem_id_2( invalid_size_type_value ),  //__e.ad_second() ),
        _M_pos_in_elem_id_2( invalid_uint16_type_value ),  //__e.pos_second() ),
        _M_face_id( __f ),
        _M_h( __e.h() ),
        _M_h_face(0),
        _M_Jt(),
        _M_Bt(),
        _M_perm( )
    {

        if ( is_linear )
            {
                _M_gm->gradient( node_t_type(), _M_g_linear );
            }
        else
            {
                _M_Jt.resize( nPoints() );
                _M_Bt.resize( nPoints() );
            }
        update( __e, __pc, __f );
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
     *  compute \f$ B(x_{\mathrm{ref}}) = K ( K^{-T} K )^{-1} \f$
     *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
     */
    void update( element_type const& __e,
                 precompute_ptrtype const& pc,
                 uint16_type __f = invalid_uint16_type_value )
    {
        _M_pc = pc;
        _M_G = __e.G();
        _M_id = __e.id();
        _M_e_marker = __e.marker();
        _M_e_marker2 = __e.marker2();
        _M_e_marker3 = __e.marker3();
        _M_face_id = __f;
        _M_h = __e.h();
        _M_xrefq = pc->nodes();

        LIFE_ASSERT( __e.G().size2() == _M_gm->nbPoints() )( __e.G().size2() )( _M_gm->nbPoints() ).error( "invalid dimensions" );
        LIFE_ASSERT( _M_pc ).error( "invalid precompute data structure" );

        if (  _M_face_id != invalid_uint16_type_value )
            {
                _M_perm = __e.permutation( _M_face_id );

                _M_h_face = __e.hFace( _M_face_id );
#if 0
                std::cout << "[geomap] _M_face_id = " << _M_face_id << "\n"
                          << "[geomap] _M_h = " << _M_h << "\n"
                          << "[geomap] _M_h_face = " << _M_h_face << "\n";
#endif

            }
        if ( vm::has_point<context>::value )
            {

                //ublas::axpy_prod( _M_G, pc->phi(), _M_xrealq, true );
                std::fill( _M_xrealq.data().begin(), _M_xrealq.data().end(), value_type(0));
                const uint16_type size1 = _M_G.size1();
                const uint16_type size3 = _M_G.size2();
                const uint16_type size2 = pc->nPoints();

                for( uint16_type i = 0; i < size1; ++i )
                    for( uint16_type j = 0; j < size2; ++j )
                        {
                            for( uint16_type k = 0; k < size3; ++k )
                                _M_xrealq( i, j ) += _M_G( i, k ) * pc->phi()[k][0][0][j];
                        }
#if 0
                std::cerr << "calling gemm\n";
                atlas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                             1.0, _M_G, __pc->phi(),
                             0.0, _M_xrealq );
                std::cerr << "done calling gemm\n";

                std::cout <<"xref = " << _M_xrefq << "\n"
                          << "xrealq = " << _M_xrealq << "\n"
                          << "G=" << _M_G << "\n";
#endif
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
    gm_ptrtype const& geometricMapping() const { return _M_gm; }

    /**
       \return the dimension of the space of the real element
    */
    uint16_type N() const { return NDim; }

    /**
       \return the dimension of the space of the reference element
    */
    uint16_type P() const { return PDim; }

    /**
     * \return the real element to which this context is associated
     */
    element_type const& element() const { return _M_element; }

    /**
     *
     */
    uint16_type nPoints() const { return _M_npoints; }

    /**
     * \return the set of points in the reference convex
     */
    matrix_node_t_type const& xRefs() const { return _M_xrefq; }

    /**
     * \return the q-th point in the reference convex
     */
    ublas::matrix_column<matrix_node_t_type const> xRef( int q ) const { return ublas::column( _M_xrefq, q ); }

    /**
     * \return the node in the real element
     */
    //matrix_node_t_type const& xReal() const
    matrix_type const& xReal() const
    {
        //BOOST_STATIC_ASSERT( vm::has_point<context>::value );
        return _M_xrealq;
    }

    /**
     * \return the node in the real element
     */
    ublas::matrix_column<matrix_type const> xReal( int q ) const
    {
        // BOOST_STATIC_ASSERT( vm::has_point<context>::value );
        return ublas::column( _M_xrealq, q );
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
            return _M_J;
        else
            return _M_Jt[q];
    }
#else
    value_type J( int q ) const
    {
        if ( is_linear )
            return _M_J;
        else
            return _M_Jt[q];
    }
#endif
#if 0
    /**
     * \internal
     */
    value_type J( int /*q*/, mpl::bool_<true> ) const  { return _M_J; }

    /**
     * \internal
     */
    value_type J( int q, mpl::bool_<false> ) const  { return _M_Jt[q]; }
#endif

    /**
     * \return the matrix associated with the geometric nodes
     */
    //matrix_node_t_type const& G() const { return _M_G; }
    matrix_type const& G() const { return _M_G; }

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

    matrix_type const& K( int i ) const
    {
        return _M_K;
    }
    /**
     * \internal
     */
    matrix_type const& B( int /*i*/, mpl::bool_<true> ) const { return _M_B; }
    /**
     * \internal
     */
    matrix_type const& B( int i, mpl::bool_<false> ) const { return _M_Bt[i]; }

    /**
     * the tensor of rank 4 for the transformation of 2nd
     * order derivatives from the reference element to the real
     * element. The tensor has the shape \f$[N,N,P,P]\f$.
     *
     * \return the tensor of rank 4
     */
    boost::multi_array<value_type,4> const& B3() const
    {
        return _M_B3;
    }

    /**
     * \return the barycenter of the reference nodes
     */
    node_t_type barycenterRef() const
    {
        node_t_type __barycenter( _M_gm->dim() );
        __barycenter.clear();

        for ( uint16_type __c = 0;__c < _M_gm->refNodes().size(); ++__c )
            {
                __barycenter += _M_gm->refNode( __c );
            }
        __barycenter /= _M_gm->refNodes().size();
        return __barycenter;
    }

    /**
       \return the barycenter of the geometric nodes
    */
    node_t_type barycenterReal() const
    {
        node_t_type __barycenter( _M_G.size1() );
        __barycenter.clear();


        for ( uint16_type __c = 0;__c < _M_G.size2(); ++__c )
            {
                __barycenter += ublas::column( _M_G, __c );
            }
        __barycenter /= _M_G.size2();
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
                return std::abs( ublas::norm_2( xReal()-ublas::column( _M_G, 0 )-
                                                ublas::prod( _M_K, xRef() ) ) ) < 1e-10;
            }
        return false;
    }

    /**
     * get the norm_2 of normal of the real element
     *
     * @return the norm_2 of the normal of the real element
     */
    value_type normalNorm( int q ) const
    {
        if ( is_linear )
            return _M_n_norm;
        else
            return _M_nnormq[q];
    }

    /**
     * get the normal of the real element
     *
     * @return the normal of the real element
     */
    node_t_type const& normal() const
    {
        //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
        return _M_n_real;
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
            return _M_n_real;
        else
            {
                return ublas::column( _M_nrealq, q );
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
        return _M_u_n_real;
    }

    //ublas::matrix_column<matrix_node_t_type const> unitNormal( int q ) const
    // node_t_type const& unitNormal( int q ) const
    node_t_type unitNormal( int q ) const
    {
        //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
        if ( is_linear )
            return _M_u_n_real;
        else
            return ublas::column( _M_unrealq, q );
    }


    value_type const& unitNormal( int n, int q ) const
    {
        //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );

        if ( is_linear )
            return _M_u_n_real( n );
        else
            return _M_unrealq( n, q );
    }

#if 0
    //ublas::matrix_column<matrix_node_t_type const> unitNormal( int q ) const
    // node_t_type const& unitNormal( int q ) const
    node_t_type unitTangent( int q ) const
    {
        //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );
        if ( is_linear )
            return _M_u_t_real;
        else
            return ublas::column( _M_utrealq, q );
    }


    value_type const& unitTangent( int n, int q ) const
    {
        //BOOST_STATIC_ASSERT( vm::has_normal<context>::value );

        if ( is_linear )
            return _M_u_t_real( n );
        else
            return _M_utrealq( n, q );
    }
#endif
    /**
     * get the id of the element
     *
     * @return the id of the element
     */
    size_type id() const { return _M_id; }

    /*
     * \return the face id
     */
    uint16_type faceId() const { return _M_face_id; }

    /**
     * \return true if the element is a face, false otherwise
     */
    bool elementIsAFace() const { return _M_face_id != invalid_uint16_type_value; }

    /**
     * get the marker of the element
     *
     * @return the marker of the element
     */
    Marker1 marker() const { return _M_e_marker; }

    /**
     * get the marker2 of the element
     *
     * @return the marker2 of the element
     */
    Marker2 marker2() const { return _M_e_marker2; }

    /**
     * get the marker3 of the element
     *
     * @return the marker3 of the element
     */
    Marker2 marker3() const { return _M_e_marker3; }

    /**
     * get the id of the first element containing the element
     *
     * @return the id of the first element containing the element
     */
    size_type id1() const { return _M_elem_id_1; }

    /**
     * get the local id of the element in the first element containing the element
     *
     * @return the local id of element in the first element containing the element
     */
    uint16_type idIn1() const { return _M_pos_in_elem_id_1; }

    /**
     * get the id of the second element containing the element
     *
     * @return the id of the second element containing the element
     */
    size_type id2() const { return _M_elem_id_2; }

    /**
     * get the local id of the element in the second element containing the element
     *
     * @return the local id of element in the second element containing the element
     */
    uint16_type idIn2() const { return _M_pos_in_elem_id_2; }

    /**
     * get an estimate of the radius of the current element
     *
     * @return the radius estimate of the current element
     */
    value_type radiusEstimate() const
    {
        return _M_gm->radiusEstimate( _M_G );
    }

    /**
     *
     *
     *
     * @return the max length of the edges of the element
     */
    value_type h() const { return _M_h; }

    /**
     * Get max length of the edge of the face of the element
     *
     * @return the max length of the edge of the face of the element
     */
    value_type hFace() const { return _M_h_face; }

    /**
     * \return the permutation associated with the face
     */
    permutation_type permutation() const
    {
        return permutation( mpl::bool_<(nDim>=2)>() );
    }

    /**
     * \return the precompute type
     */
    precompute_ptrtype const& pc() { return _M_pc; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
       set some precomputed data on the reference element
    */
    void setPc( precompute_ptrtype const& __pc )
    {
        _M_pc = __pc;
    }

    //@}
private:

    /**
     * \return the permutation associated with the face
     */
    permutation_type permutation( mpl::bool_<false> ) const
    {
        return _M_perm;
    }

    /**
     * \return the permutation associated with the face
     */
    permutation_type permutation( mpl::bool_<true> ) const
    {
        LIFE_ASSERT( _M_face_id == invalid_uint16_type_value ||
                     (_M_face_id != invalid_uint16_type_value &&
                      _M_perm != permutation_type( permutation_type::NO_PERMUTATION ) ) )
            ( _M_face_id ).error( "invalid permutation" );
        return _M_perm;
    }

    /**
     * update Jacobian data : linear case
     */
    void updateJKBN( mpl::bool_<true>  )
    {

        if ( !_M_gm->isCached() ||
             (_M_gm->isCached() && _M_gm->cached( _M_id ) == false) )
            {
#if 0
                if ( boost::is_arithmetic<value_type>::value )
                    atlas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                                 1.0, _M_G, _M_g_linear,
                                 0.0, _M_K );
                else
#endif
                    ublas::axpy_prod( _M_G, _M_g_linear, _M_K, true );

                if ( NDim == PDim )
                    {
                        _M_J = math::abs( det<NDim>( _M_K ) );
                        //if ( vm::has_kb<context>::value )
                            {
#if 0
                                inverse<NDim>( _M_K, _M_CS, _M_J );
                                ublas::noalias(_M_B) = ublas::trans( _M_CS );
#else
                                inverse<NDim>( _M_K, _M_CS );
                                ublas::noalias(_M_B) = ublas::trans( _M_CS );
#endif
                            }
                    }
                else // N != P
                    {
                        // CS = K^T K
#if 0
                        if ( boost::is_arithmetic<value_type>::value )
                            atlas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                                         1.0, _M_K, _M_K,
                                         0.0, _M_CS );
                        else
#endif

                            ublas::noalias(_M_CS) = ublas::prod( ublas::trans( _M_K ), _M_K );

                        _M_J = math::sqrt( math::abs( det<PDim>( _M_CS ) ) );
                        //if ( vm::has_kb<context>::value )
                            {
                                inverse<PDim>( _M_CS, _M_CSi );
                                // B = K CS
#if 0
                                if ( boost::is_arithmetic<value_type>::value )
                                    atlas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                                                 1.0, _M_K, _M_CSi,
                                                 0.0, _M_B );
                                else
#endif
                                    ublas::axpy_prod( _M_K, _M_CSi, _M_B, true );
                            }

                    }

                if ( _M_gm->isCached() )
                    {
                        // cache J, K and B
                        _M_gm->addJ( _M_id, _M_J );
                        //if ( vm::has_kb<context>::value )
                            {
                                _M_gm->addK( _M_id, _M_K );
                                _M_gm->addB( _M_id, _M_B );
                            }
                        _M_gm->setCached( _M_id, true );

                        //Log() << "(add to cache) J[" << _M_id << "]=" <<  _M_J << "\n";
                        //Log() << "(add to cache) B[" << _M_id << "]=" <<  _M_B << "\n";
                    }
            }
        else
            {
                _M_J = _M_gm->J( _M_id );
                //Log() << "(use cache) J[" << _M_id << "]=" <<  _M_J << "\n";
                //if ( vm::has_kb<context>::value )
                    {
                        _M_K = _M_gm->K( _M_id );
                        _M_B = _M_gm->B( _M_id );
                        //Log() << "(use cache) B[" << _M_id << "]=" <<  _M_B << "\n";
                    }


            }
        if ( vm::has_hessian<context>::value )
            {

                for (uint16_type k = 0; k < NDim; ++k)
                    for (uint16_type l = 0; l < NDim; ++l)
                        for (uint16_type i = 0; i < PDim; ++i)
                            for (uint16_type j = 0; j < PDim; ++j)
                                _M_B3[k][l][i][j] = _M_B(k,i)*_M_B(l,j);
#if 0
                //B3(k + N_*l, i + P*j) = BB(k, i) * BB(l, j);
                std::cout << "element " << this->id() << " B3 = " << "\n";
                for( int i = 0; i < nDim; ++i )
                    for( int j = 0; j < nDim; ++j )
                        for( typename boost::multi_array<value_type,4>::index k = 0; k < nDim; ++k )
                            for( typename boost::multi_array<value_type,4>::index l = 0; l < nDim; ++l )
                                std::cout << "B3[" << i << "][" << j << "][" << k << "][" << l << "]="
                                          << this->B3()[i][j][k][l] << "\n";
#endif
            }
        if ( ((NDim != PDim) || ( vm::has_normal<context>::value) ) && (_M_face_id != invalid_uint16_type_value ) )
            {
#if 0
                blas::gemv( traits::NO_TRANSPOSE,
                            1.0, _M_Bt[ q ], _M_n[_M_face_id],
                            0.0, _M_n_real );
#else
                ublas::axpy_prod( _M_B,
                                  _M_gm->referenceConvex().normal( _M_face_id ),
                                  _M_n_real,
                                  true );
#endif
                _M_n_norm = ublas::norm_2( _M_n_real );
                _M_u_n_real = _M_n_real/_M_n_norm;

                if ( vm::has_tangent<context>::value )
                    {

                    }
            }
    }

    /**
     * update Jacobian data : nonlinear case
     */
    void updateJKBN( mpl::bool_<false>  )
    {
        //std::cout << "nPoints() =" << nPoints() << "\n";
        //Debug() << "[geomap] G = "<< _M_G << "\n";
        //double res = 0;
        for ( int q = 0; q < nPoints(); ++q )
            {
                //std::cout << "q =" << q << "\n";
                _M_gm->gradient( q, _M_g, _M_pc.get().get() );
                //Debug() << "[geomap] g[" << q << "] = "<< _M_g << "\n";

#if 0
                blas::gemm( _M_G, _M_g, _M_K );
#else
                ublas::axpy_prod( _M_G, _M_g, _M_K, true );
#endif

                if ( NDim == PDim )
                    {
                        _M_J = math::abs( det<NDim>( _M_K ) );
                        if ( vm::has_kb<context>::value )
                            {
                                inverse<NDim>( _M_K, _M_CS );
                                //ublas::noalias(_M_B) = ublas::trans( _M_CS );
                                _M_B = ublas::trans( _M_CS );
                                //std::cout << "========== B[" << q << "]=" << _M_B << "\n";
                            }
                    }
                else // N != P
                    {
                        // CS = K^T K
#if 1
                        ublas::prod( ublas::trans( _M_K ), _M_K, _M_CS );
#else
                        blas::gemm( traits::TRANSPOSE, traits::NO_TRANSPOSE,
                                    1.0, _M_K, _M_K,
                                    0.0, _M_CS );
#endif
                        _M_J = math::sqrt( math::abs( det<PDim>( _M_CS ) ) );
                        if ( vm::has_kb<context>::value )
                            {
                                inverse<PDim>( _M_CS, _M_CSi );
                                // B = K CS
#if 1
                                ublas::axpy_prod( _M_K, _M_CSi, _M_B );
#else
                                blas::gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE,
                                            1.0, _M_K, _M_CSi,
                                            0.0, _M_B );
#endif
                            }

                    }
                //Debug() << "[geomap] J[" << q << "]= "<< _M_J << "\n";
                //res += _M_J;
                // store q-th jacobian entry

                _M_Jt[q] = _M_J;

                if ( vm::has_kb<context>::value )
                    {
                        //Debug() << "[geomap] B[" << q << "]= "<< _M_B << "\n";
                        _M_Bt[q].resize( _M_B.size1(), _M_B.size2() );
                        _M_Bt[q] = _M_B;
                    }


            }
        //Debug() << "[geomap] res(sum J) = " << res << "\n";
        if ( ((NDim != PDim) || (vm::has_normal<context>::value)) && (_M_face_id != invalid_uint16_type_value ) )
            {
                //std::cout << "has normal\n";
                for ( int q = 0;q < nPoints(); ++q )
                    {
#if 0
                        blas::gemv( traits::NO_TRANSPOSE,
                                    1.0, _M_Bt[ q ], _M_n[_M_face_id],
                                    0.0, _M_n_real );
#else
                        if ( 0 ) //trans == fem::LINEAR )
                            {
                                ublas::axpy_prod( _M_B,
                                                  _M_gm->referenceConvex().normal( _M_face_id ),
                                                  _M_n_real,
                                                  true );
                            }
                        else
                            ublas::axpy_prod( _M_Bt[ q ],
                                              _M_gm->referenceConvex().normal( _M_face_id ),
                                              _M_n_real,
                                              true );
#endif
                        //std::cout << "[geomap] point " << q << " n_real = " << _M_n_real << "\n";
                        _M_n_norm = ublas::norm_2( _M_n_real );
                        _M_u_n_real = _M_n_real/_M_n_norm;
                        ublas::column( _M_nrealq, q ) = _M_n_real;
                        ublas::column( _M_unrealq, q ) = _M_u_n_real;
                        _M_nnormq[q] = _M_n_norm;

                        if ( NDim != PDim )
                            _M_Jt[q] *= _M_n_norm;
                        if ( vm::has_tangent<context>::value )
                            {

                            }
                    }

#if 0
                std::cout << "[geomap] face id = " << _M_face_id << "\n"
                          << "[geomap] ref normal = " << _M_gm->referenceConvex().normal( _M_face_id ) << "\n"
                          << "[geomap] _M_n_real = " << _M_nrealq << "\n"
                          << "[geomap] _M_unrealq = " << _M_unrealq << "\n"
                          << "[geomap] _M_nnormq = " << _M_nnormq << "\n";
#endif
            }
    }

    Context();
    Context( Context const& );
private:

    gm_ptrtype _M_gm;

    element_type const& _M_element;

    boost::optional<precompute_ptrtype> _M_pc;
    uint16_type _M_npoints;

    value_type _M_J;

    //matrix_node_t_type  _M_G;
    matrix_type  _M_G;
    ublas::vector<node_t_type> _M_n;
    node_t_type _M_n_real;
    node_t_type _M_u_n_real;
    value_type _M_n_norm;

    //matrix_node_t_type const& _M_xrefq;
    matrix_node_t_type  _M_xrefq;
    //matrix_node_t_type _M_xrealq;
    matrix_type _M_xrealq;
    matrix_node_t_type _M_nrealq;
    matrix_node_t_type _M_unrealq;
    ublas::vector<value_type> _M_nnormq;

    matrix_type _M_g_linear;
    matrix_type _M_g;
    matrix_type _M_K;
    matrix_type _M_CS;
    matrix_type _M_CSi;
    matrix_type _M_B;
    boost::multi_array<value_type,4> _M_B3;

    size_type _M_id;
    Marker1 _M_e_marker;
    Marker2 _M_e_marker2;
    Marker3 _M_e_marker3;
    size_type _M_elem_id_1;
    uint16_type _M_pos_in_elem_id_1;
    size_type _M_elem_id_2;
    uint16_type _M_pos_in_elem_id_2;

    uint16_type _M_face_id;

    value_type _M_h;
    value_type _M_h_face;

    vector_type _M_Jt;
    std::vector<matrix_type> _M_Bt;

    permutation_type _M_perm;
}; // Context


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


    template<typename GeoElem>
    Inverse( geometric_mapping_ptrtype __gm, GeoElem const& __ge )
        :
        _M_gm( __gm ),
        _M_xref( __gm->dim() ),
        _M_xreal( __ge.G().size1() ),
        _M_is_in( false ),
        _M_G( __ge.G() ),
        _M_K( N(), __gm->dim() ),
        _M_B( N(), __gm->dim() ),
        _M_CS( __gm->dim(), __gm->dim() ),
        _M_g( _M_gm->nbPoints(), __gm->dim() ),
#if defined( HAVE_PETSC )
        _M_nlsolver( SolverNonLinear<double>::build( SOLVERS_PETSC ) )
#else
        _M_nlsolver( SolverNonLinear<double>::build( SOLVERS_GMM ) )
#endif
    {
        LIFE_ASSERT( _M_G.size2() == __gm->nbPoints() )
            ( _M_G.size2() )( __gm->nbPoints() ).error( "invalid dimensions" );

        if ( _M_gm->isLinear() )
            update();
        else
            {
#if defined( HAVE_PETSC )
                _M_nlsolver->dense_residual = boost::bind( &Inverse::updateResidual, boost::ref( *this ), _1, _2 );
                _M_nlsolver->dense_jacobian = boost::bind( &Inverse::updateJacobian, boost::ref( *this ), _1, _2 );
#else

#endif
            }
    }



    /** @name Accessors
     */
    //@{

    /**
       \return the geometric mapping associated with the context
    */
    geometric_mapping_ptrtype const& geometricMapping() const { return _M_gm; }

    /**
       \return the dimension of the space of the real element
    */
    uint16_type N() const { return G().size1(); }

    /**
       \return the node in the reference element
    */
    node_t_type const& xRef() const { return _M_xref; }

    /**
       \return the node in the real element
    */
    node_t_type const& xReal() const { return _M_xreal; }

    /**
       \return the matrix associated with the geometric nodes
    */
    matrix_type const& G() const { return _M_G; }

    /**
       \return the gradient of the transformation at reference node
    */
    matrix_type const& K() const { return _M_K; }

    value_type J() const { return math::abs( det<Dim>( _M_K ) ); }

    /**
       \return the pseudo-inverse of the gradient of the
       transformation at reference node
    */
    matrix_type const& B() const { return _M_B; }

    /**
       \return the barycenter of the reference nodes
    */
    node_type barycenterRef() const;

    /**
       \return the barycenter of the geometric nodes
    */
    node_type barycenterReal() const;

    /**
       \return \c true if the node is in the convex defined by G,
       \c false otherwise
    */
    bool isIn() const { return _M_is_in; }

    /**
     * tell whether the point is on the surface of the convex
     * @return true if the point is on the surface, false otherwise
     */
    bool isOnConvexSurface() const
    {
        if ( _M_gm->isLinear() )
            {
                // x -x0 - K(0)\bar{x}
                return std::abs( ublas::norm_2( xReal()-ublas::column( _M_G, 0 )-
                                                ublas::prod( _M_K, xRef() ) ) ) < 1e-10;
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
        _M_xreal = __xreal;

        if ( _M_gm->isLinear() )
            {
                update();
                _M_is_in = linearInverse();
            }
        else
            {
                _M_is_in = nonLinearInversePetsc();

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
        scalar_type operator()(const node_t_type& x) const
        {
            node_type r = gmi.geometricMapping()->transform( x, gmi.G() ) - xreal;
            return ublas::inner_prod( r, r )/2.;
        }
        void operator()(const node_t_type& x, node_t_type& gr) const
        {
            //             gmi.setXreal( xreal );
            gmi._M_xref.assign( x );
            gmi.update();
            node_type r = gmi.geometricMapping()->transform( x, gmi.G() ) - xreal;
            gr.resize(x.size());
            ublas::prod( ublas::trans( gmi.K() ), r, gr );
        }
    private:

        Inverse& gmi;
        node_type xreal;

    };
    typedef ublas::vector<double> dense_vector_type;
    typedef ublas::matrix<double> dense_matrix_type;

    void updateResidual( dense_vector_type const& x, dense_vector_type& r )
    {
        dense_vector_type y = _M_gm->transform( x, _M_G );
        r = y - _M_xreal ;
#if 0
        Log() << "[geomap::residual] begin ------------------------------\n";
        Log() << "[geomap::residual] x =" << x << "\n";
        Log() << "[geomap::residual] _M_G =" << _M_G << "\n";

        Log() << "[geomap::residual] y =" << y << "\n";
        Log() << "[geomap::residual] xreal =" << _M_xreal << "\n";
        Log() << "[geomap::residual] r(xreal-y) =" << r << "\n";
        Log() << "[geomap::residual] end   ------------------------------\n";
#endif // 0

    }

    void updateJacobian( dense_vector_type const& x, dense_matrix_type& j )
    {
        _M_gm->gradient( x, _M_g );
        ublas::axpy_prod( _M_G, _M_g, j );
#if 0
        Log() << "[geomap::jacobian] begin ------------------------------\n";
        Log() << "[geomap::jacobian] x =" << x << "\n";
        Log() << "[geomap::jacobian] j =" << j << "\n";
        Log() << "[geomap::jacobian] end   ------------------------------\n";
#endif
    }

    /**
       update information on this context

       -# update the coordinate of the real points
       -# update the pseudo-inverse of the gradient of the transformation

       compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
       where \f$G \f$ is the matrix representing the geometric nodes nDof x dim

       compute \f$ B(x_{\mathrm{ref}}) = K ( K^{-T} K )^{-1} \f$
       where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
    */
    void update()
    {
        // xref is not defined at this point: it serves only as a dummy point
        // in the case of linear inversion (__g is a constant matrix in this case)
        // in the non linear case , xRef() is update in the newton iterations
        _M_gm->gradient( xRef(), _M_g );
        Debug( 5046 ) << "[update] g = " << _M_g << "\n";

        checkInvariant();

        ublas::axpy_prod( _M_G, _M_g, _M_K );
        Debug( 5046 ) << "[update] K(0) = " << _M_K << "\n";

        // compute B
        if ( _M_gm->dim() != N() )
            {
                ublas::prod( ublas::trans( _M_K ), _M_K, _M_CS );
                LU<matrix_type> __lu( _M_CS );
                __lu.inverse( _M_CS );
                ublas::axpy_prod( _M_K, _M_CS, _M_B );
            }
        else
            {
                LU<matrix_type> __lu( _M_K );
                __lu.inverse( _M_CS );
                _M_B = ublas::trans( _M_CS );
            }

        Debug( 5046 ) << "[update] B(0) = " << _M_B << "\n";
    }

    void
    checkInvariant() const
    {
        LIFE_ASSERT( _M_G.size2() == _M_g.size1() )( _M_G.size2() )( _M_g.size1() ).error( "G,g invalid dimensions" );
        LIFE_ASSERT( _M_G.size1() == _M_K.size1() )( _M_G.size1() )( _M_K.size1() ).error( "G,K invalid dimensions" );
        LIFE_ASSERT( _M_g.size2() == _M_K.size2() )( _M_g.size2() )( _M_K.size2() ).error( "g,K invalid dimensions" );
        LIFE_ASSERT( _M_B.size2() == _M_gm->dim() )( _M_B.size1() )( N() ).error( "B,gm invalid dimensions" );
    }


    bool linearInverse()
    {
        checkInvariant();

        size_type N = _M_xreal.size();
        size_type P = _M_xref.size();

        node_type y(_M_xreal);

        Debug( 5046 ) << "y = xreal = " << y << "\n";
        //Debug( 5046 ) << "G(0)  = " << node_type( _M_x0 << "\n";
        y.minus_assign(  ublas::column( _M_G, 0 ) );
        Debug( 5046 ) << "y - G(0) = " << y << "\n";

        Debug( 5046 ) << "B(0) = " << _M_B << "\n";
        Debug( 5046 ) << "xref = " << ublas::prod( ublas::trans( _M_B ), y ) << "\n";

        // xref = B^T * y = B^T * ( x_real - x_0)
        _M_xref.assign( ublas::prod( ublas::trans( _M_B ), y )-ublas::scalar_vector<value_type>( P, 1.0 ) );

        Debug( 5046 ) << "[GeoMap::Inverse::linearInverse] xref : " << _M_xref << "\n";

        bool __isin;
        double vmin;
        boost::tie( __isin, vmin ) = _M_gm->isIn( _M_xref );
        Debug( 5046 ) << "[GeoMap::Inverse::linearInverse] isIn : " << __isin << "\n";
        ///if ( __isin < 1e-10 )
        if ( __isin )
            {
                if ( N == P )
                    return true;
                else
                    {
                        // y = y - K * x_ref
                        ublas::axpy_prod( _M_K, -_M_xref, y );
                        if ( ublas::norm_2( y ) < 1e-10 )
                            return true;
                    }
            }
        return false;
    }

    matrix_node_t_type linearInversePoints( matrix_node_t_type const& real_pts, bool /*allow_extrapolation*/ = false ) const
    {
        return ublas::prod( ublas::trans( _M_B ),
                            real_pts-ublas::outer_prod( ublas::column( _M_G, 0 ),
                                                        ublas::scalar_vector<value_type>( real_pts.size2(), value_type( 1 ) ) ) ) -
            ublas::scalar_matrix<value_type>( _M_B.size2(), real_pts.size2(), value_type( 1 ) );

    }

    /*
     *inversion for non-linear geometric transformations
     *  (Newton on Grad(pgt)(y - pgt(x)) = 0 )
     */
    bool nonLinearInversePetsc()
    {
        //Log() << "starting new nonlinear inverse\n";
        //const double EPS = 1e-10;
        const double IN_EPS = 1e-10;
        size_type N = _M_xreal.size();
        size_type P = _M_xref.size();

        /*
          find an initial guess: closest geometric node to _M_xreal
        */
        node_type x0 = _M_gm->refNode( 0 );
        node_type y = ublas::column( _M_G, 0 );
        scalar_type d = ublas::inner_prod( y, _M_xreal );
        for (size_type j = 1; j < _M_gm->nbPoints(); ++j)
            {
                scalar_type d2 = ublas::inner_prod( ublas::column( _M_G, j ), _M_xreal );
                if (d2 < d)
                    {
                        d = d2;
                        x0 = _M_gm->refNode(j);
                        y = ublas::column( _M_G, j );
                    }
            }
        _M_xref = x0;

        dense_matrix_type J( N, N );
        dense_vector_type R( N );

        updateResidual( _M_xref, R );
        updateJacobian( _M_xref, J );

        // find xref by solving the non linear equation
        _M_nlsolver->solve( J, _M_xref, R, 1e-10, 10 );


        // compute the location of xref: inside or outside the element
        bool __isin;
        double vmin;
        this->updateResidual( _M_xref, R );
        boost::tie( __isin, vmin ) = _M_gm->isIn( _M_xref );
        if ( __isin  &&
             ( P == N || ublas::norm_2(R) < IN_EPS) )
            {
                //Log() << "point " << _M_xref << "in IN (" << vmin << ") residual = " << ublas::norm_2(R) << "\n";
                return true;
            }
        else
            {
                //Log() << "point " << _M_xref << "in OUT (" << vmin << ") residual = " << ublas::norm_2(R) << "\n";
            }
        //Log() << "done in new nonlinear inverse\n";
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
        size_type N = _M_xreal.size();
        size_type P = _M_xref.size();

        /*
          find an initial guess: closest geometric node to _M_xreal
        */
        node_type x0 = _M_gm->refNode( 0 );
        node_type y = ublas::column( _M_G, 0 );
        scalar_type d = ublas::inner_prod( y, _M_xreal );
        for (size_type j = 1; j < _M_gm->nbPoints(); ++j)
            {
                scalar_type d2 = ublas::inner_prod( ublas::column( _M_G, j ), _M_xreal );
                if (d2 < d)
                    {
                        d = d2;
                        x0 = _M_gm->refNode(j);
                        y = ublas::column( _M_G, j );
                    }
            }
        _M_xref = x0;

        node_type vres(N);
        node_type rn( _M_xreal );
        rn.minus_assign( y );

        this->update();

        // vres = K^T rn
        ublas::axpy_prod( rn, _M_K, vres );
        scalar_type res = ublas::norm_2(vres);

        //     std::cerr << "DEBUT: res0=" << res << ", X=" << _M_xreal << "\nB=" << _M_B << ", K=" << _M_K << "\n";
        unsigned cnt = 50;
        while (res > EPS/10 && cnt > 0 )
            {
                node_type xn( P );

                // xn = B^T rn
                ublas::axpy_prod( rn, _M_B, xn );

                scalar_type newres;
                for (int16_type i=1; i<=256; i*=2 )
                    {
                        _M_xref.plus_assign( xn / scalar_type(i) );
                        y = _M_gm->transform( _M_xref, G() );

                        rn = _M_xreal - y;

                        this->update();

                        if (P != N)
                            {
                                // vres = K^T rn
                                ublas::axpy_prod( rn, _M_K, vres );
                                newres = ublas::norm_2( vres );
                            }
                        else
                            {
                                newres = ublas::norm_2( rn );
                            }

                        if (newres < 1.5*res)
                            break;
                    }
                res = newres;
                //         std::cout << "cnt=" << cnt << ", x=" << _M_xref << ", res=" << res << "\n";
                --cnt;
            }
        //     std::cerr << " invert_nonlin done\n";
        //     std::cerr << "cnt=" << cnt << ", P=" << P << ", N=" << N
        //               << "\nX=" << _M_xreal << " Xref=" << _M_xref << "\nresidu=" << res << "\n";
        //<< ", G=" << G << "\nX=" << _M_xreal << " Xref=" << x << "\nresidu=" << res << "\nB=" << B << ", K=" << K << "\n"  << "\n-------------------^^^^^^^^\n";
        if (cnt == 0)
            {
#if 0
                std::cerr << "BFGS in geotrans_inv_convex!\n";
                GeomapInverseConvex b( *this, _M_xreal );

                iteration_ptrtype iter( Iteration<double>::New() );
                iter->setMaximumNumberOfIterations( 50 );
                iter->setRelativePrecision( 1e-8 );

                node_type x( x0 );
                bfgs(b,b,x,10,*iter);
                rn = _M_gm->transform(x,G()) - _M_xreal;

                if ( _M_gm->isIn(x) < IN_EPS &&
                     N==P && ublas::norm_2(rn) > IN_EPS)
                    throw "inversion of non-linear geometric transformation "
                        "failed ! (too much iterations)";
#endif
            }

        bool __isin;
        double vmin;
        boost::tie( __isin, vmin ) = _M_gm->isIn( _M_xref );
        if ( __isin  &&
             ( P == N || ublas::norm_2(rn) < IN_EPS) )
            {

                return true;
            }
        else
            {
                //Log() << "point " << _M_xref << "in OUT (" << vmin << ")\n";
            }
        return false;
    }



private:

    geometric_mapping_ptrtype _M_gm;
    node_type _M_xref;
    node_type _M_xreal;

    bool _M_is_in;

    //matrix_type const& _M_G;
    matrix_type _M_G;
    matrix_type _M_K;
    matrix_type _M_B;
    matrix_type _M_CS;
    matrix_type _M_g;

    boost::shared_ptr<SolverNonLinear<double> > _M_nlsolver;
}; // Inverse

private:

element_gm_ptrtype _elementMap;
face_gm_ptrtype _boundaryMap;
matrix_type _M_g_linear;

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

template<int Dim, int Order,  template<uint16_type,uint16_type,uint16_type> class Entity = Simplex, typename T = double>
struct GT_Lagrange
{};

template<int Dim, int Order,  typename T = double>
struct GT_QK
{};

# /* List of dims. */
# define LIFE_GEOMAP                                   \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           1,                           \
                           (                            \
                            Lagrange                    \
                                                    )   \
                                                    )   \
    /**/
# /* List of dims. */
# define LIFE_DIMS                                     \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           3,                           \
                           (                            \
                            1,2,3                       \
                                                    )   \
                                                    )   \
    /**/
# /* List of orders. */
# define LIFE_ORDERS                                   \
    BOOST_PP_TUPLE_TO_LIST(                             \
                           5,                           \
                           (                            \
                            1,2,3,4,5                   \
                                                    )   \
                                                    )   \
    /**/
#
# define LIFE_ENTITY BOOST_PP_TUPLE_TO_LIST( 2, ( Simplex,SimplexProduct ) )
    /**/
#
# define LIFE_GEN_GT(GEOM,LDIM,LORDER)         \
    "GT_" #GEOM "(" #LDIM "," #LORDER ")"       \
    /**/
#
# /* Generates code for all dim and order. */
# define LIFE_GT_FACTORY_OP(_, GDO)            \
    LIFE_GT_FACTORY GDO                        \
    /**/
#
#
#define LIFE_GT_FACTORY(GEOM,LDIM,LORDER,ENTITY)                       \
    template<typename T>                                                \
    struct BOOST_PP_CAT(GT_, GEOM)<LDIM, LORDER, ENTITY, T>             \
        :                                                               \
        public GeoMap<LDIM, LORDER, T, ENTITY, GEOM>                    \
    {                                                                   \
        typedef mpl::vector<boost::none_t,                              \
            boost::none_t,                                              \
            GeoMap<1, LORDER, T, ENTITY, GEOM>,                    \
            GeoMap<2, LORDER, T, ENTITY, GEOM> > geomap_faces_t;   \
    typedef mpl::vector<GeoMap<1, LORDER, T, ENTITY, GEOM>,        \
            GeoMap<2, LORDER, T, ENTITY, GEOM>,                    \
            GeoMap<3, LORDER, T, ENTITY, GEOM>,                    \
            boost::none_t> geomap_elements_t;                           \
    typedef typename type_traits<T>::value_type value_type;             \
    typedef GeoMap<LDIM, LORDER, T, ENTITY, GEOM> super;           \
    typedef BOOST_PP_CAT(GT_,GEOM)<LDIM-1, LORDER, ENTITY, T> face_geo_type; \
                                                                        \
    static const uint16_type nDim = LDIM;                                   \
    static const uint16_type order = LORDER;                                \
    static const uint16_type nDof = super::nDof;                            \
    static const uint16_type nNodes = super::nNodes;                        \
    typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim> >::type element_gm_type; \
    typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim> >::type face_gm_type; \
    template<int N>                                                     \
        struct face_gm                                                  \
    {                                                                   \
        typedef typename mpl::at<geomap_faces_t, mpl::int_<N> >::type type; \
    };                                                                  \
                                                                        \
    BOOST_PP_CAT(GT_,GEOM)()                                            \
        :                                                               \
        super( boost::shared_ptr<element_gm_type>(new element_gm_type()), boost::shared_ptr<face_gm_type>(new face_gm_type() )) \
        {}                                                              \
    };                                                                  \
    /**/
BOOST_PP_LIST_FOR_EACH_PRODUCT(LIFE_GT_FACTORY_OP, 4, (LIFE_GEOMAP, LIFE_DIMS, LIFE_ORDERS, LIFE_ENTITY))

#undef LIFE_DIMS
#undef LIFE_ORDERS


template<typename Elem, template<uint16_type,uint16_type,uint16_type> class Entity = Simplex, typename T = double>
class RealToReference
{
public:
    static const uint16_type nDim = Elem::nDim;
    typedef T value_type;
    typedef typename matrix_node<T>::type points_type;
    typedef GT_Lagrange<nDim,1, Entity, value_type> gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::Inverse inverge_gm_type;

    RealToReference( Elem const& elem )
        :
        _M_gm( new gm_type ),
        _M_igm( _M_gm, elem )
    {}

    points_type operator()( points_type const& pts ) const
    {
        return _M_igm( pts );
    }

    value_type J() const { return _M_igm.J(); }

private:

    gm_ptrtype _M_gm;
    inverge_gm_type _M_igm;

};

} // Life
#endif
