/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-20

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   \file integrators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-20
 */
#ifndef FEELPP_INTEGRATORS_HPP
#define FEELPP_INTEGRATORS_HPP 1

#include <cxxabi.h>
#include <typeinfo>
#include <boost/timer.hpp>

#include <Eigen/Eigen>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/quadmapped.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/detail/clean.hpp>
#include <feel/feelvf/block.hpp>

#include <feel/feelvf/formcontextbase.hpp>
#include <feel/feelvf/bilinearform.hpp>
#include <feel/feelvf/linearform.hpp>
#include <feel/feeldiscr/quadptlocalization.hpp>
#if defined( FEELPP_HAS_GOOGLE_PROFILER_H )
#include <google/profiler.h>
#endif

#include <feel/feelvf/detail/integrator.hpp>

#if defined(FEELPP_HAS_HARTS)
#include "RunTimeSystem/Model/RunTimeSysEnv.h"
#include "RunTimeSystem/DataMng/DataHandler.h"
#include "RunTimeSystem/TaskMng/TaskMng.h"
#include "RunTimeSystem/TaskMng/AsynchTask.h"
#include "RunTimeSystem/TaskMng/StdScheduler.h"
#include "RunTimeSystem/TaskMng/StdDriver.h"
#include "RunTimeSystem/TaskMng/TBBScheduler.h"
#include "RunTimeSystem/TaskMng/TBBDriver.h"
#include "RunTimeSystem/TaskMng/PTHScheduler.h"
#include "RunTimeSystem/TaskMng/PTHDriver.h"

#include "RunTimeSystem/DataMng/DataArgs.h"

#include "Utils/PerfTools/PerfCounterMng.h"
#endif //defined(FEELPP_HAS_HARTS)

namespace Feel
{
namespace vf
{

/// \cond detail
enum IntegratorType
{
    INT_0D = -1,                 /**< 0D integrator */
    INT_1D = 1,                 /**< 1D integrator */
    INT_2D = 2,                 /**< 2D integrator */
    INT_3D = 3,                 /**< 3D integrator */
    INT_ELEMENTS = 4,           /**< integration over elements */
    INT_ELEMENT_FACES = 5,      /**< integration over face elements */
    INT_INTERNAL_FACES = 6      /**< integration over internal faces */

};

/**
 * \class Integrator
 * \brief base class for integrators
 *
 * @author Christophe Prud'homme
 * @see IntegratorOn
 */
template<typename Elements, typename Im, typename Expr, typename Im2=Im>
class Integrator: public IntegratorBase
{
public:

    /** @name Constants
     */
    //@{

    static const size_type context = Expr::context;
    static const bool is_terminal = false;

    static const uint16_type imorder = 0;
    static const bool imIsPoly = true;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = Expr::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = Expr::template HasTrialFunction<Func>::result;
    };

    static const size_type iDim = boost::tuples::template element<0, Elements>::type::value;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef Integrator<Elements, Im, Expr, Im2> self_type;

    typedef typename boost::tuples::template element<1, Elements>::type element_iterator;

    typedef Expr expression_type;
    typedef typename expression_type::value_type expression_value_type;
    typedef ublas::vector<int> vector_dof_type;

    struct eval
    {
        //
        // some typedefs
        //
        typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
        typedef typename boost::remove_const<const_t>::type the_face_element_type;
        typedef typename the_face_element_type::super2::template Element<the_face_element_type>::type the_element_type;

        typedef typename mpl::if_<mpl::bool_<the_element_type::is_simplex>,
                mpl::identity<typename Im::template apply<the_element_type::nDim, expression_value_type, Simplex>::type >,
                    mpl::identity<typename Im::template apply<the_element_type::nDim, expression_value_type, Hypercube>::type >
        >::type::type im_type;

        typedef typename mpl::if_<mpl::bool_<the_element_type::is_simplex>,
                mpl::identity<typename Im2::template apply<the_element_type::nDim, expression_value_type, Simplex>::type >,
                    mpl::identity<typename Im2::template apply<the_element_type::nDim, expression_value_type, Hypercube>::type >
        >::type::type im2_type;

        typedef the_element_type element_type;
        typedef typename the_element_type::gm_type gm_type;
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename the_element_type::gm1_type gm1_type;
        typedef boost::shared_ptr<gm1_type> gm1_ptrtype;
        //typedef typename gm_type::template Context<expression_type::context, the_element_type, im_type::numPoints> gmc_type;
        typedef typename gm_type::template Context<expression_type::context|vm::JACOBIAN, the_element_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm1_type::template Context<expression_type::context|vm::JACOBIAN, the_element_type> gmc1_type;
        typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;
#if 0
        typedef typename gm_type::template precompute<im_type::numPoints>::type gmpc_type;
        typedef typename gm_type::template precompute<im_type::numPoints>::ptrtype gmpc_ptrtype;
#else
        typedef typename gm_type::PreCompute gmpc_type;
        typedef boost::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef typename gm1_type::PreCompute gmpc1_type;
        typedef boost::shared_ptr<gmpc1_type> gmpc1_ptrtype;
#endif
        //typedef typename eval_expr_type::value_type value_type;
        //typedef typename strongest_numeric_type<typename Im::value_type, typename expression_type::value_type>::type value_type;
        //typedef typename expression_type::value_type value_type;

        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef typename eval_expr_type::shape shape;

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
        typedef typename expression_type::template tensor<map_gmc1_type> eval_expr1_type;
        typedef typename eval_expr1_type::shape shape1; // should be the same as shape
        //typedef typename shape_type::storage<value_type> storage_type;
        /*
          typedef mpl::if_<mpl::bool_<shape_type::is_scalar>,
          mpl::identity<value_type>,
          mpl::identity<ublas::vector<value_type,storage_type> > >::type::type ret_type;*/
        //typedef ublas::matrix<value_type> ret_type;
#if 0
        typedef typename mpl::if_<mpl::and_<mpl::equal_to<mpl::int_<shape::M>, mpl::int_<1> >,
                mpl::equal_to<mpl::int_<shape::N>, mpl::int_<1> > >,
                mpl::identity<expression_value_type>,
                mpl::identity<Eigen::Matrix<expression_value_type, shape::M, shape::N> > >::type::type value_type;
#else
        typedef Eigen::Matrix<expression_value_type, shape::M, shape::N> value_type;
        typedef Eigen::Matrix<expression_value_type, shape::M, shape::N> evaluate_type;
#endif
        typedef Eigen::Matrix<expression_value_type, shape::M, shape::N> matrix_type;
        static value_type zero( mpl::bool_<false> )
        {
            return value_type::Zero();
        }
        static value_type zero( mpl::bool_<true> )
        {
            return 0;
        }
        static value_type zero()
        {
            return zero( boost::is_scalar<value_type>() );
        }

    };

    typedef typename eval::im_type im_type;
    typedef typename eval::im2_type im2_type;
    typedef typename im_type::face_quadrature_type im_face_type;
    typedef typename im2_type::face_quadrature_type im2_face_type;
    //typedef typename eval::value_type value_type;
    typedef typename eval::matrix_type matrix_type;
    typedef typename eval::matrix_type value_type;
    typedef typename eval::matrix_type evaluate_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Integrator( Elements const& elts, Im const& /*__im*/, expression_type const& __expr, GeomapStrategyType gt, Im2 const& /*__im2*/, bool use_tbb, bool use_harts, int grainsize, std::string const& partitioner,
                boost::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > qpl )
        :
        M_elts(),
        M_eltbegin( elts.template get<1>() ),
        M_eltend( elts.template get<2>() ),
        M_im( ),
        M_im2( ),
        M_expr( __expr ),
        M_gt( gt ),
        M_use_tbb( use_tbb ),
        M_use_harts( use_harts ),
        M_grainsize( grainsize ),
        M_partitioner( partitioner ),
        M_QPL( qpl )
    {
        M_elts.push_back( elts );
        DLOG(INFO) << "Integrator constructor from expression\n";
    }

    Integrator( std::list<Elements> const& elts, Im const& /*__im*/, expression_type const& __expr,
                GeomapStrategyType gt, Im2 const& /*__im2*/, bool use_tbb, bool use_harts, int grainsize, std::string const& partitioner,
                boost::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > qpl )
        :
        M_elts( elts ),
        M_im( ),
        M_im2( ),
        M_expr( __expr ),
        M_gt( gt ),
        M_use_tbb( use_tbb ),
        M_use_harts( use_harts ),
        M_grainsize( grainsize ),
        M_partitioner( partitioner ),
        M_QPL( qpl )
    {
        DLOG(INFO) << "Integrator constructor from expression\n";
        if ( elts.size() )
        {
            M_eltbegin = elts.begin()->template get<1>();
            M_eltend = elts.begin()->template get<2>();
        }
    }

    Integrator( Integrator const& __vfi )
        :
        M_elts( __vfi.M_elts) ,
        M_eltbegin( __vfi.M_eltbegin ),
        M_eltend( __vfi.M_eltend ),
        M_im( __vfi.M_im ),
        M_im2( __vfi.M_im2 ),
        M_expr( __vfi.M_expr ),
        M_gt( __vfi.M_gt ),
        M_use_tbb( __vfi.M_use_tbb ),
        M_use_harts( __vfi.M_use_harts ),
        M_grainsize( __vfi.M_grainsize ),
        M_partitioner( __vfi.M_partitioner ),
        M_QPL( __vfi.M_QPL )
    {
        DLOG(INFO) << "Integrator copy constructor\n";
    }

    virtual ~Integrator() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    template<typename... TheExpr>
    struct Lambda
    {
        typedef typename expression_type::template Lambda<TheExpr...>::type expr_type;
        typedef _Q< ExpressionOrder<Elements,expr_type>::value > quad_type;
        typedef _Q< ExpressionOrder<Elements,expr_type>::value_1 > quad1_type;
        typedef Integrator<Elements, quad_type, expr_type, quad1_type> type;
    };

    template<typename... ExprT>
    typename Lambda<ExprT...>::type
    operator()( ExprT...  e )
        {
#if 0
            typedef decltype(expr(M_expr(e...))) t1;
            typedef typename Lambda<ExprT>::expr_type t2;
            BOOST_MPL_ASSERT_MSG( (boost::is_same<t1,t2>), INVALID_TYPE_IN_META_EXPRESSION,
                                  (decltype(expr(M_expr(e...))), typename Lambda<ExprT>::expr_type ) );
#endif
            auto new_expr = M_expr(e...);
            typedef decltype(new_expr) expr_type;
            typedef typename Lambda<ExprT...>::expr_type e_type;
            typedef _Q< ExpressionOrder<Elements,e_type>::value > quad_type;
            typedef _Q< ExpressionOrder<Elements,e_type>::value_1 > quad1_type;
            typedef boost::shared_ptr<QuadPtLocalization<Elements,quad_type,expr_type > > quadptloc_ptrtype;
            quad_type quad;
            quad1_type quad1;
            //BOOST_STATIC_ASSERT( ( boost::is_same<expr_type,e_type> ) );
            auto i = Integrator<Elements, quad_type, expr_type, quad1_type>( M_elts, quad, new_expr, M_gt, quad1, M_use_tbb, M_use_harts, M_grainsize, M_partitioner, quadptloc_ptrtype() );
            DLOG(INFO) << " -- M_elts size=" << M_elts.size() << "\n";
            DLOG(INFO) << " -- nelts=" << std::distance( M_eltbegin, M_eltend ) << "\n";
            DLOG(INFO) << " -- integrate: quad = " << i.im().nPoints() << "\n";
            DLOG(INFO) << " -- integrate: quad1 = " << i.im2().nPoints() << "\n";
            i.setBeginElement( M_eltbegin );
            i.setEndElement( M_eltend );
            return i;
        }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get the integration method
     *
     *
     * @return the integration method
     */
    im_type const& im() const
    {
        return M_im;
    }

    /**
     * get the integration method
     *
     *
     * @return the integration method
     */
    im2_type const& im2() const
    {
        return M_im2;
    }

    /**
     * get the integration method on face f
     *
     *
     * @return the integration method on face f
     */
    im_face_type  im( uint16_type f ) const
    {
        return M_im.face( f );
    }

    /**
     * get the integration method on face f
     *
     *
     * @return the integration method on face f
     */
    im2_face_type  im2( uint16_type f ) const
    {
        return M_im2.face( f );
    }

    /**
     * get the variational expression
     *
     *
     * @return the variational expression
     */
    expression_type const& expression() const
    {
        return M_expr;
    }


    /**
     * get the geometric mapping integrator type
     *
     * @return the geometric mapping integrator type
     */
    GeomapStrategyType geomapIntegratorType() const
    {
        return M_gt;
    }

    /**
     * iterator that points at the beginning of the container that
     * holds the data that will apply the integration upon
     */
    element_iterator beginElement() const
    {
        return M_eltbegin;
    }

    /**
     * iterator that points at the end of the container that
     * holds the data that will apply the integration upon
     */
    element_iterator endElement() const
    {
        return M_eltend;
    }

    /**
     * tell whether the expression is symmetric or not
     *
     *
     * @return true if symmetric, false otherwise
     */
    bool isSymmetric() const
    {
        return M_expr.isSymmetric();
    }

    //@}

    /** @name  Mutators
     */
    //@{
    void setBeginElement( element_iterator it ) { M_eltbegin = it; }
    void setEndElement( element_iterator en ) { M_eltend = en; }

    //@}

    /** @name  Methods
     */
    //@{


    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __u,
                   boost::shared_ptr<Elem2> const& __v,
                   FormType& __form ) const;

    template<typename Elem1, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __v,
                   FormType& __form ) const;

    //typename expression_type::template tensor<Geo_t>::value_type
    template<typename P0hType>
    typename P0hType::element_type  broken( boost::shared_ptr<P0hType>& P0h ) const
    {
        auto p0 = broken( P0h, mpl::int_<iDim>() );
        return p0;
    }
#if 1
    matrix_type
    evaluate( bool parallel=true,
              WorldComm const& worldcomm = Environment::worldComm() ) const
#else
    //typename expression_type::template tensor<Geo_t>::value_type
    BOOST_PARAMETER_MEMBER_FUNCTION( ( matrix_type ),
                                     evaluate,
                                     tag,
                                     /*(required
                                     //(h,*(double)))*/
                                     ( optional
                                       ( parallel,*( bool ), true ) ) )
#endif
    {
        typename eval::matrix_type loc =  evaluate( mpl::int_<iDim>() );

        if ( !parallel )
            return loc;

        else // parallel
        {
            typename eval::matrix_type glo( loc );
            // maybe better to create anoter worldcomm which split the mesh worldcomm
            // with only partition that contains at least one element (Vincent C.)
            // and thus argument worldComm can be remove
            // auto const& worldcomm = const_cast<MeshBase*>( this->beginElement()->mesh() )->worldComm();

            if ( worldcomm.localSize() > 1 )
            {
                mpi::all_reduce( worldcomm.localComm(),
                                 loc,
                                 glo,
                                 [] ( matrix_type const& x, matrix_type const& y )
                {
                    return x + y;
                } );
            }

            return glo;
        }
    }

#if defined( FEELPP_HAS_TBB )
    template<typename FormType, typename ExprType, typename IMType, typename EltType>
    class Context
    {
    public:
        //
        // some typedefs
        //
        typedef typename eval::gm_type gm_type;
        typedef typename eval::gmc_type gmc_type;
        typedef typename eval::gmpc_type gmpc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef boost::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
        //typedef vf::detail::FormContextBase<map_gmc_type,im_type> fcb_type;
        typedef form_context_type fcb_type;
        typedef boost::shared_ptr<fcb_type> fcb_ptrtype;


        Context( FormType& _form,
                 ExprType const& _expr,
                 IMType const& _im,
                 EltType const& _elt )
            :
            M_geopc( new gmpc_type( _form.gm(), _im.points() ) ),
            M_c( new gmc_type( _form.gm(), _elt, M_geopc ) ),
            M_formc( new form_context_type( _form,
                                            fusion::make_pair<vf::detail::gmc<0> >( M_c ),
                                            fusion::make_pair<vf::detail::gmc<0> >( M_c ),
                                            _expr,
                                            _im ) )

        {
        }
        Context( Context const& c )
            :
            M_geopc( new gmpc_type( *c.M_geopc ) ),
            M_c( new gmc_type( *c.M_c ) ),
            M_formc( new form_context_type( *c.M_formc ) )
        {
        }
        //std::vector<boost::reference_wrapper<const typename mesh_type::element_type> > _v;
        typedef typename std::vector<boost::reference_wrapper<const typename eval::element_type> >::iterator elt_iterator;
        void operator() ( const tbb::blocked_range<elt_iterator>& r ) const
        {
#if 1
            tbb::mutex m;
            tbb::mutex::scoped_lock lock( m  );
            lock.release();
            auto  mapgmc = fusion::make_pair<vf::detail::gmc<0> >( M_c );

            for ( auto _elt = r.begin(); _elt != r.end(); ++_elt )
            {
                M_c->update( *_elt );
                M_formc->update( mapgmc, mapgmc );
                M_formc->integrate();
                lock.acquire( m );
                M_formc->assemble();
                lock.release();
            }

#else

            for ( auto _elt = r.begin(); _elt != r.end(); ++_elt )
            {
                M_c->update( *_elt );
                M_formc->update( fusion::make_pair<vf::detail::gmc<0> >( M_c ),
                                 fusion::make_pair<vf::detail::gmc<0> >( M_c ) );
                M_formc->integrate();
            }

#endif
        }
        mutable gmpc_ptrtype M_geopc;
        mutable gmc_ptrtype M_c;
        mutable fcb_ptrtype M_formc;

    };

    template<typename ExprType, typename IMType, typename EltType>
    class ContextEvaluate
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //
        // some typedefs
        //
        typedef ExprType expression_type;
        typedef typename eval::gm_type gm_type;
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename eval::gmc_type gmc_type;
        typedef typename eval::gmpc_type gmpc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef boost::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef typename eval_expr_type::shape shape;
        typedef IMType im_type;
        typedef typename eval::matrix_type value_type;
        ContextEvaluate( ExprType const& _expr,
                         IMType const& _im,
                         EltType const& _elt )
            :
            M_gm( new gm_type( *_elt.gm() ) ),
            M_geopc( new gmpc_type( M_gm, _im.points() ) ),
            M_c( new gmc_type( M_gm, _elt, M_geopc ) ),
            M_expr( _expr, map_gmc_type( fusion::make_pair<vf::detail::gmc<0> >( M_c ) ) ),
            M_im( _im ),
            M_ret( eval::matrix_type::Zero() )
        {
        }
        ContextEvaluate( ContextEvaluate& c, tbb::split )
            :
            M_gm( new gm_type( *c.M_gm ) ),
            //M_geopc( new gmpc_type( M_gm, c.M_im.points() ) ),
            M_geopc( c.M_geopc ),
            M_c( new gmc_type( M_gm, c.M_c->element(), M_geopc ) ),
            M_expr( c.M_expr ),
            M_im( c.M_im ),
            M_ret( eval::matrix_type::Zero() )
        {}

        ContextEvaluate( ContextEvaluate const& c )
            :
            M_gm( new gm_type( *c.M_gm ) ),
            //M_geopc( new gmpc_type( M_gm, c.M_im.points() ) ),
            M_geopc( c.M_geopc ),
            M_c( new gmc_type( M_gm, c.M_c->element(), M_geopc ) ),
            M_expr( c.M_expr ),
            M_im( c.M_im ),
            M_ret( c.M_ret )
        {
        }
        //std::vector<boost::reference_wrapper<const typename mesh_type::element_type> > _v;
        typedef typename std::vector<boost::reference_wrapper<const typename eval::element_type> >::iterator elt_iterator;
        void operator() ( const tbb::blocked_range<elt_iterator>& r )
        {
#if 1

            for ( auto _elt = r.begin(); _elt != r.end(); ++_elt )
            {
                M_c->update( *_elt );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( M_c ) );

                M_expr.update( mapgmc );
                M_im.update( *M_c );

#if 1

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    {
                        M_ret( c1,c2 ) += M_im( M_expr, c1, c2 );
                    }

#endif
            }

#else
#if 0

            for ( int i = 0; i < 10000; ++i )
                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    {
                        M_ret( c1,c2 ) += i*i; //M_im( M_expr, c1, c2 );
                    }

#endif
#endif
        }
        void join( ContextEvaluate const& other )
        {
            M_ret += other.M_ret;
        }

        value_type result() const
        {
            return M_ret;
        }

        gm_ptrtype M_gm;
        gmpc_ptrtype M_geopc;
        gmc_ptrtype M_c;
        eval_expr_type M_expr;
        im_type M_im;
        value_type M_ret;
    };
#endif // FEELPP_HAS_TBB

    //@}


private:

    template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
    bool useSameMesh( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                      FaceRangeType const& faceRange ) const;
    template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
    bool useSameMesh( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                      FaceRangeType const& faceRange ) const;

    template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
    boost::tuple<size_type,rank_type,uint16_type>
    testElt0IdFromFaceRange( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                             FaceRangeType const& faceRange ) const;
    template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
    boost::tuple<size_type,rank_type,uint16_type>
    testElt0IdFromFaceRange( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                             FaceRangeType const& faceRange ) const;

    template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
    bool faceIntegratorUseTwoConnections( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                          FaceRangeType const& faceRange,
                                          typename vf::detail::BilinearForm<FE1,FE2,ElemContType>::mesh_1_type::element_type const& eltTest,
                                          uint16_type faceIdInElt ) const;
    template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
    bool faceIntegratorUseTwoConnections( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                                          FaceRangeType const& faceRange,
                                          typename vf::detail::LinearForm<FE,VectorType,ElemContType>::mesh_type::element_type const& eltTest,
                                          uint16_type faceIdInElt ) const;



    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool hasRelation  ) const
        {
            assemble( __form, mpl::int_<MESH_ELEMENTS>() /**/, mpl::bool_<true>() /**/, hasRelation,
                      mpl::bool_<FormType::test_space_type::is_mortar||FormType::trial_space_type::is_mortar>()  );
        }
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool hasRelation, mpl::false_  ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool hasRelation, mpl::true_  ) const;

    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<false> /**/, bool hasRelation ) const;

    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<true> /**/, bool hasRelation ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<false> /**/, bool hasRelation ) const;

    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
        {
            LOG(INFO) << "Intergrator::assembleWithRelationDifferentMeshType mortar case:" << FE1::is_mortar||FE2::is_mortar;
            assembleWithRelationDifferentMeshType( __form, mpl::int_<MESH_ELEMENTS>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::true_ ) const;

    template<typename FE,typename VectorType,typename ElemContType>
    void assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
        {
            LOG(INFO) << "Intergrator::assembleWithRelationDifferentMeshType mortar case:" << FE1::is_mortar||FE2::is_mortar;
            assembleWithRelationDifferentMeshType( __form, mpl::int_<MESH_FACES>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/, mpl::false_ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/, mpl::true_ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;

    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
        {
            LOG(INFO) << "Intergrator::assembleInCaseOfInterpolate :" << FE1::is_mortar||FE2::is_mortar;
            assembleInCaseOfInterpolate( __form, mpl::int_<MESH_ELEMENTS>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }

    // bilinear interpolation case (no relation)
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
        {
            LOG(INFO) << "Intergrator::assembleInCaseOfInterpolate mortar case:" << FE1::is_mortar||FE2::is_mortar;
            assembleInCaseOfInterpolate( __form, mpl::int_<MESH_FACES>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::true_ ) const;


    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/, mpl::false_ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/, mpl::true_ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;


    template<typename P0hType>
    typename P0hType::element_type  broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const;
    template<typename P0hType>
    typename P0hType::element_type  broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const;

    typename eval::matrix_type evaluate( mpl::int_<MESH_ELEMENTS> ) const;
    typename eval::matrix_type evaluate( mpl::int_<MESH_FACES> ) const;
    typename eval::matrix_type evaluate( mpl::int_<MESH_POINTS> ) const;

private:

    std::list<Elements> M_elts;
    element_iterator M_eltbegin;
    element_iterator M_eltend;
    mutable im_type M_im;
    mutable im2_type M_im2;
    expression_type   M_expr;
    GeomapStrategyType M_gt;
    bool M_use_tbb;
    bool M_use_harts;
    int M_grainsize;
    std::string M_partitioner;
    mutable boost::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > M_QPL;
    //     mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, double, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, double> > > M_profile_local_assembly;

    //     mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, double, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, double> > > M_profile_global_assembly;
};

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename Elem1, typename Elem2, typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( boost::shared_ptr<Elem1> const& __u,
        boost::shared_ptr<Elem2> const& __v,
        FormType& __form ) const
{
#if 0
    details::GlobalMatrixAssembler<iDim, self_type, Elem1, Elem2, FormType>( *this,
            __u,
            __v,
            __form );
#endif

    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same1_mesh_type;
    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem2::mesh_type::element_type>::type same2_mesh_type;
    typedef typename boost::mpl::and_< same1_mesh_type,same2_mesh_type>::type same_mesh_type;

    element_iterator it, en;
    // get one elt for init
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        it = lit->template get<1>();
        en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( it == en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;

    // run assemble process according to isRelated meshes
    const bool test_related_to_trial = __v->mesh()->isRelatedTo( __u->mesh() ) && ( __u->mesh()->isRelatedTo( it->mesh() ) || __v->mesh()->isRelatedTo( it->mesh() ) );
    LOG(INFO) << "[integrator::assemble bilinear form] with_relation_mesh (same_mesh: " << same_mesh_type::value
              << " test_related_to_trial: " << test_related_to_trial << ")\n";
    if ( test_related_to_trial )
    {
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<same_mesh_type::value>(), test_related_to_trial );
    }
    else
    {
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<false>(), test_related_to_trial );
    }

}


template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename Elem1, typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( boost::shared_ptr<Elem1> const& __v,
        FormType& __form ) const
{
#if 0
    details::GlobalVectorAssembler<iDim, self_type, Elem1, FormType>( *this,
            __v,
            __form );
#endif

    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same_mesh_type;

    element_iterator it, en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        it = lit->template get<1>();
        en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( it == en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;


    const bool test_related_to_range = __v->mesh()->isRelatedTo( it->mesh() );

    //if ( dynamic_cast<void*>( const_cast<MeshBase*>( it->mesh() ) ) == dynamic_cast<void*>( __v->mesh().get() ) )
    if ( test_related_to_range )
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<same_mesh_type::value>(), test_related_to_range/*true*/ );

    else
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<false>(), test_related_to_range/*false*/ );

    //assemble( __form, mpl::int_<iDim>(), mpl::bool_<true>() );
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool /*hasRelation*/, mpl::false_ ) const
{
    boost::timer __timer;
    LOG(INFO) << "[integrator::assemble FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true>\n";

#if defined(FEELPP_HAS_TBB)

    //std::cout << "Integrator Uses TBB: " << M_use_tbb << "\n";
    if ( !M_use_tbb )
#else
        if ( 1 )
#endif
        {
            //
            // some typedefs
            //
            typedef typename eval::gm_type gm_type;
            typedef typename eval::gm1_type gm1_type;
            typedef typename eval::gmc_type gmc_type;
            typedef typename eval::gmc1_type gmc1_type;
            typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
            typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;

            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
            typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
            typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type> form1_context_type;
            //typedef vf::detail::FormContextBase<map_gmc_type,im_type> fcb_type;
            typedef form_context_type fcb_type;
            typedef form1_context_type fcb1_type;
            typedef fcb_type* focb_ptrtype;
            typedef fcb1_type* focb1_ptrtype;

            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //
            typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( __form.gm(), this->im().points() ) );
            typename eval::gmpc1_ptrtype __geopc1( new typename eval::gmpc1_type( __form.gm1(), this->im2().points() ) );


            for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
            {
                element_iterator it = lit->template get<1>();
                element_iterator en = lit->template get<2>();
                DLOG(INFO) << "integrating over " << std::distance( it, en )  << " elements\n";

                // check that we have elements to iterate over
                if ( it == en )
                    continue;

                size_type idEltTestInit = it->id();
                if ( it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                    idEltTestInit = it->mesh()->subMeshToMesh( idEltTestInit );
                else if ( __form.testSpace()->mesh()->isSubMeshFrom( it->mesh() ) )
                    idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltTestInit );

                auto const& eltTestInit = __form.testSpace()->mesh()->element( idEltTestInit );


                gmc_ptrtype __c( new gmc_type( __form.gm(), eltTestInit, __geopc ) );
                gmc1_ptrtype __c1( new gmc1_type( __form.gm1(), eltTestInit, __geopc1 ) );


                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

                focb_ptrtype formc( new form_context_type( __form,
                                                           mapgmc,
                                                           mapgmc,
                                                           mapgmc,
                                                           this->expression(),
                                                           this->im() ) );
                focb1_ptrtype formc1( new form1_context_type( __form,
                                                              mapgmc1,
                                                              mapgmc1,
                                                              mapgmc1,
                                                              this->expression(),
                                                              this->im2() ) );

                //int nelt = std::distance( this->beginElement(), this->endElement() );
                boost::timer ti0,ti1, ti2, ti3;

                //double t0 = 0, t1 = 0,t2 = 0,t3 = 0;
                //
                // start the real intensive job:
                // -# iterate over all elements to integrate over
                // -# construct the associated geometric mapping with the reference element
                // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
                // -# assemble the local contribution in the global representation of the bilinear form
                //
                for ( ; it != en; ++it )
                {
                    size_type idElt = it->id();
                    if ( it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                        idElt = it->mesh()->subMeshToMesh( idElt );
                    else if ( __form.testSpace()->mesh()->isSubMeshFrom( it->mesh() ) )
                        idElt = __form.testSpace()->mesh()->meshToSubMesh( idElt );

                    auto const& eltTest = __form.testSpace()->mesh()->element( idElt );

                    if ( formc->isZero( idElt ) )
                        continue;

                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        //ti0.restart();
                        __c->update( eltTest );
                        //t0+=ti0.elapsed();
#if 0
                        std::cout << "Element: " << it->id() << "\n"
                                  << " o - points : " << it->G() << "\n"
                                  << " o - quadrature :\n"
                                  << "     ref : " << this->im().points() << "\n"
                                  << "     real : " << __c->xReal() << "\n";
#endif
                        //ti1.restart();
                        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                        formc->update( mapgmc,mapgmc,mapgmc );
                        //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                        //t1+=ti1.elapsed();

                        //ti2.restart();
                        formc->integrate();
                        //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                        //t2+=ti2.elapsed();

                        //ti3.restart();
                        formc->assemble();
                        //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                        //t3+=ti3.elapsed();
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_O1:
                    {
                        //ti0.restart();
                        __c1->update( eltTest );
                        //t0+=ti0.elapsed();
#if 0
                        DLOG(INFO) << "Element: " << it->id() << "\n"
                                      << " o - points : " << it->G() << "\n"
                                      << " o - quadrature :\n"
                                      << "     ref : " << this->im().points() << "\n"
                                      << "     real : " << __c->xReal() << "\n";
#endif
                        //ti1.restart();
                        map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                        formc1->update( mapgmc1,mapgmc1,mapgmc1 );
                        //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                        //t1+=ti1.elapsed();

                        //ti2.restart();
                        formc1->integrate();
                        //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                        //t2+=ti2.elapsed();

                        //ti3.restart();
                        formc1->assemble();
                        //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                        //t3+=ti3.elapsed();
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_OPT:
                    {
                        if ( it->isOnBoundary() )
                        {
                            //ti0.restart();
                            __c->update( eltTest );
                            //t0+=ti0.elapsed();
#if 0
                            DLOG(INFO) << "Element: " << it->id() << "\n"
                                          << " o - points : " << it->G() << "\n"
                                          << " o - quadrature :\n"
                                          << "     ref : " << this->im().points() << "\n"
                                          << "     real : " << __c->xReal() << "\n";
#endif
                            //ti1.restart();
                            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                            formc->update( mapgmc,mapgmc,mapgmc );

                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formc->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formc->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }

                        else
                        {
                            //ti0.restart();
                            __c1->update( eltTest );
                            //t0+=ti0.elapsed();
#if 0
                            DLOG(INFO) << "Element: " << it->id() << "\n"
                                          << " o - points : " << it->G() << "\n"
                                          << " o - quadrature :\n"
                                          << "     ref : " << this->im().points() << "\n"
                                          << "     real : " << __c->xReal() << "\n";
#endif
                            //ti1.restart();
                            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                            formc1->update( mapgmc1,mapgmc1,mapgmc1 );

                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formc1->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formc1->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }
                    }
                    break;
                    }
                } // end loop on elements
                delete formc;
                delete formc1;
            }// end loop on list of elements
#if 0
            DLOG(INFO) << "[elements] Overall geometric mapping update time : " << ( t0+t1+t2 ) << " per element:" << ( t0+t1+t2 )/std::distance( this->beginElement(), this->endElement() ) << "\n";
            DLOG(INFO) << "[elements] Overall geometric mapping update time : " << t0 << "\n";
            DLOG(INFO) << "[elements] Overall form update time : " << t1 << "\n";
            DLOG(INFO) << "[elements] Overall local assembly time : " << t2 << "\n";
            DLOG(INFO) << "[elements] Overall global assembly time : " << t3 << "\n";
#endif
            DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";


        }

#if defined( FEELPP_HAS_TBB )

        else
        {
            element_iterator it = this->beginElement();
            element_iterator en = this->endElement();

            if ( it == en )
                return;

            std::vector<boost::reference_wrapper<const typename eval::element_type> > _v;

            for ( auto _it = it; _it != en; ++_it )
                _v.push_back( boost::cref( *_it ) );

            //tbb::blocked_range<decltype(_v.begin())> r( _v.begin(), _v.end(), M_grainsize );
            tbb::blocked_range<decltype( _v.begin() )> r( _v.begin(), _v.end(), std::distance( it, en ) );
            Context<FormType,expression_type, im_type, typename eval::the_element_type> thecontext ( __form,
                                                                                                     this->expression(),
                                                                                                     this->im(),
                                                                                                     *it );
            tbb::parallel_for( r,  thecontext );
        }

#endif // FEELPP_HAS_TBB
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/,
                                               mpl::bool_<true> /**/, bool /*hasRelation*/, mpl::true_ ) const
{
    boost::timer __timer;
    LOG(INFO) << "[integrator::assemble FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> mortar case\n";

#if defined(FEELPP_HAS_TBB)

    //std::cout << "Integrator Uses TBB: " << M_use_tbb << "\n";
    if ( !M_use_tbb )
#else
        if ( 1 )
#endif
        {
            //
            // some typedefs
            //
            typedef typename eval::gm_type gm_type;
            typedef typename eval::gm1_type gm1_type;
            typedef typename eval::gmc_type gmc_type;
            typedef typename eval::gmc1_type gmc1_type;
            typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
            typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;

            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
            typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
            typedef typename FormType::template Context<map_gmc_type, expression_type, im_type,map_gmc_type,map_gmc_type,true> form_mortar_context_type;
            typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type> form1_context_type;
            typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type,map_gmc1_type,map_gmc1_type,true> form1_mortar_context_type;
            //typedef vf::detail::FormContextBase<map_gmc_type,im_type> fcb_type;
            typedef form_context_type fcb_type;
            typedef form_mortar_context_type mortar_fcb_type;
            typedef form1_context_type fcb1_type;
            typedef form1_mortar_context_type mortar_fcb1_type;
            typedef fcb_type* focb_ptrtype;
            typedef mortar_fcb_type* mortar_focb_ptrtype;
            typedef fcb1_type* focb1_ptrtype;
            typedef mortar_fcb1_type* mortar_focb1_ptrtype;

            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //
            typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( __form.gm(), this->im().points() ) );
            typename eval::gmpc1_ptrtype __geopc1( new typename eval::gmpc1_type( __form.gm1(), this->im2().points() ) );


            for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
            {
                element_iterator it = lit->template get<1>();
                element_iterator en = lit->template get<2>();
                DLOG(INFO) << "integrating over " << std::distance( it, en )  << " elements\n";

                // check that we have elements to iterate over
                if ( it == en )
                    continue;

                size_type idEltTestInit = it->id();
                if ( it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                    idEltTestInit = it->mesh()->subMeshToMesh( idEltTestInit );
                else if ( __form.testSpace()->mesh()->isSubMeshFrom( it->mesh() ) )
                    idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltTestInit );

                auto const& eltTestInit = __form.testSpace()->mesh()->element( idEltTestInit );


                gmc_ptrtype __c( new gmc_type( __form.gm(), eltTestInit, __geopc ) );
                gmc1_ptrtype __c1( new gmc1_type( __form.gm1(), eltTestInit, __geopc1 ) );


                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

                focb_ptrtype formc( new form_context_type( __form,
                                                           mapgmc,
                                                           mapgmc,
                                                           mapgmc,
                                                           this->expression(),
                                                           this->im() ) );
                mortar_focb_ptrtype formcm( new form_mortar_context_type( __form,
                                                                           mapgmc,
                                                                           mapgmc,
                                                                           mapgmc,
                                                                           this->expression(),
                                                                           this->im() ) );
                focb1_ptrtype formc1( new form1_context_type( __form,
                                                              mapgmc1,
                                                              mapgmc1,
                                                              mapgmc1,
                                                              this->expression(),
                                                              this->im2() ) );
                mortar_focb1_ptrtype formc1m( new form1_mortar_context_type( __form,
                                                                             mapgmc1,
                                                                             mapgmc1,
                                                                             mapgmc1,
                                                                             this->expression(),
                                                                             this->im2() ) );

                //int nelt = std::distance( this->beginElement(), this->endElement() );
                boost::timer ti0,ti1, ti2, ti3;

                //double t0 = 0, t1 = 0,t2 = 0,t3 = 0;
                //
                // start the real intensive job:
                // -# iterate over all elements to integrate over
                // -# construct the associated geometric mapping with the reference element
                // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
                // -# assemble the local contribution in the global representation of the bilinear form
                //
                for ( ; it != en; ++it )
                {
                    size_type idElt = it->id();
                    if ( it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                        idElt = it->mesh()->subMeshToMesh( idElt );
                    else if ( __form.testSpace()->mesh()->isSubMeshFrom( it->mesh() ) )
                        idElt = __form.testSpace()->mesh()->meshToSubMesh( idElt );

                    auto const& eltTest = __form.testSpace()->mesh()->element( idElt );

                    if ( formc->isZero( idElt ) )
                        continue;

                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        //ti0.restart();
                        __c->update( eltTest );
                        //t0+=ti0.elapsed();
#if 0
                        std::cout << "Element: " << it->id() << "\n"
                                  << " o - points : " << it->G() << "\n"
                                  << " o - quadrature :\n"
                                  << "     ref : " << this->im().points() << "\n"
                                  << "     real : " << __c->xReal() << "\n";
#endif
                        //ti1.restart();
                        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

                        if ( eltTest.isOnBoundary() )
                        {
                            formcm->update( mapgmc,mapgmc,mapgmc );
                            formcm->integrate();
                            formcm->assemble();
                        }
                        else
                        {
                            formc->update( mapgmc,mapgmc,mapgmc );
                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formc->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formc->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_O1:
                    {
                        //ti0.restart();
                        __c1->update( eltTest );
                        //t0+=ti0.elapsed();
#if 0
                        DLOG(INFO) << "Element: " << it->id() << "\n"
                                      << " o - points : " << it->G() << "\n"
                                      << " o - quadrature :\n"
                                      << "     ref : " << this->im().points() << "\n"
                                      << "     real : " << __c->xReal() << "\n";
#endif
                        //ti1.restart();
                        map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

                        if ( eltTest.isOnBoundary() )
                        {
                            formc1m->update( mapgmc,mapgmc,mapgmc );
                            formc1m->integrate();
                            formc1m->assemble();
                        }
                        else
                        {
                            formc1->update( mapgmc1,mapgmc1,mapgmc1 );
                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formc1->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formc1->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_OPT:
                    {
                        if ( it->isOnBoundary() )
                        {
                            //ti0.restart();
                            __c->update( eltTest );
                            //t0+=ti0.elapsed();
#if 0
                            DLOG(INFO) << "Element: " << it->id() << "\n"
                                          << " o - points : " << it->G() << "\n"
                                          << " o - quadrature :\n"
                                          << "     ref : " << this->im().points() << "\n"
                                          << "     real : " << __c->xReal() << "\n";
#endif
                            //ti1.restart();
                            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                            formcm->update( mapgmc,mapgmc,mapgmc );

                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formcm->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formcm->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }

                        else
                        {
                            //ti0.restart();
                            __c1->update( eltTest );
                            //t0+=ti0.elapsed();
#if 0
                            DLOG(INFO) << "Element: " << it->id() << "\n"
                                          << " o - points : " << it->G() << "\n"
                                          << " o - quadrature :\n"
                                          << "     ref : " << this->im().points() << "\n"
                                          << "     real : " << __c->xReal() << "\n";
#endif
                            //ti1.restart();
                            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                            formc1->update( mapgmc1,mapgmc1,mapgmc1 );

                            //DLOG(INFO)  << "update gmc : " << ti1.elapsed() << "\n";
                            //t1+=ti1.elapsed();

                            //ti2.restart();
                            formc1->integrate();
                            //DLOG(INFO)  << "integrate : " << ti2.elapsed() << "\n";
                            //t2+=ti2.elapsed();

                            //ti3.restart();
                            formc1->assemble();
                            //DLOG(INFO)  << "assemble : " << ti3.elapsed() << "\n";
                            //t3+=ti3.elapsed();
                        }
                    }
                    break;
                    }
                } // end loop on elements
                delete formc;
                delete formc1;
                delete formcm;
                delete formc1m;
            }// end loop on list of elements
#if 0
            DLOG(INFO) << "[elements] Overall geometric mapping update time : " << ( t0+t1+t2 ) << " per element:" << ( t0+t1+t2 )/std::distance( this->beginElement(), this->endElement() ) << "\n";
            DLOG(INFO) << "[elements] Overall geometric mapping update time : " << t0 << "\n";
            DLOG(INFO) << "[elements] Overall form update time : " << t1 << "\n";
            DLOG(INFO) << "[elements] Overall local assembly time : " << t2 << "\n";
            DLOG(INFO) << "[elements] Overall global assembly time : " << t3 << "\n";
#endif
            DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";


        }

#if defined( FEELPP_HAS_TBB )

        else
        {
            element_iterator it = this->beginElement();
            element_iterator en = this->endElement();

            if ( it == en )
                return;

            std::vector<boost::reference_wrapper<const typename eval::element_type> > _v;

            for ( auto _it = it; _it != en; ++_it )
                _v.push_back( boost::cref( *_it ) );

            //tbb::blocked_range<decltype(_v.begin())> r( _v.begin(), _v.end(), M_grainsize );
            tbb::blocked_range<decltype( _v.begin() )> r( _v.begin(), _v.end(), std::distance( it, en ) );
            Context<FormType,expression_type, im_type, typename eval::the_element_type> thecontext ( __form,
                                                                                                     this->expression(),
                                                                                                     this->im(),
                                                                                                     *it );
            tbb::parallel_for( r,  thecontext );
        }

#endif // FEELPP_HAS_TBB
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<false> /**/, bool hasRelation ) const
{
    if ( hasRelation )
    {
        assembleWithRelationDifferentMeshType( __form,mpl::int_<MESH_ELEMENTS>() );
    }
    else
    {
        assembleInCaseOfInterpolate( __form,mpl::int_<MESH_ELEMENTS>() );
    }
}

namespace detail
{
template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::true_ )
{
    return gmcExpr;
}

template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::false_ )
{
    CHECK( false ) << "not allowed\n";
    return boost::shared_ptr<GmcType>();
}

template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/ )
{
    return buildGmcWithRelationDifferentMeshType<SpaceType,ImType,GmcType,GmcExprType>( space,gm,im, idElt, gmcExpr,
                                                                                        mpl::int_<0>(),
                                                                                        mpl::bool_<boost::is_same<GmcType,GmcExprType>::type::value>() );
}
template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                       size_type idElt,boost::shared_ptr<GmcExprType> /**/ ,mpl::int_<1> /**/ )
{
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;
    typedef GmcType gmc_type;
#if 0
    typedef typename SpaceType::gm_type::precompute_type geopc_type;
    typedef typename SpaceType::gm_type::precompute_ptrtype geopc_ptrtype;
#else
    typedef typename gmc_type::gm_type::precompute_type geopc_type;
    typedef typename gmc_type::gm_type::precompute_ptrtype geopc_ptrtype;
#endif
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( im.nFaces() );

    for ( uint16_type __f = 0; __f < im.nFaces(); ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = geopc_ptrtype( new geopc_type( gm, im.fpoints(__f, __p.value() ) ) );
        }
    }

    auto const& eltInit = space->mesh()->face( idElt );
    gmc_ptrtype gmc( new gmc_type( gm,  eltInit.element0(), __geopc, eltInit.pos_first() /*face_id_in_elt_0*/ ) );

    return gmc;
}
template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& /*space*/,
                                        boost::shared_ptr<GmcType> /*gmc*/, boost::shared_ptr<GmcExprType> /*gmcExpr*/,
                                        size_type /*idElt*/, mpl::int_<0> /**/ )
{
    // nothing to do!
}
template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType( boost::shared_ptr<SpaceType> const& space,
                                        boost::shared_ptr<GmcType> gmc, boost::shared_ptr<GmcExprType> gmcExpr,
                                        size_type idElt, mpl::int_<1> /**/ )
{
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;

    auto const& theface = space->mesh()->face( idElt );
    bool findPermutation=false;
    for ( permutation_type __p( permutation_type::IDENTITY );
          __p < permutation_type( permutation_type::N_PERMUTATIONS ) && !findPermutation; ++__p )
    {
        // update only xReal in gmc
        gmc->update( theface.element0(), theface.pos_first(), __p, false );

        bool check=true;
        for ( uint16_type i=0;i<gmc->nPoints() && check;++i )
            for (uint16_type d=0;d<GmcType::NDim;++d)
                check = check && ( std::abs(gmc->xReal(i)[d]-gmcExpr->xReal(i)[d])<1e-8 );

        // if check compute full gmc context with the good permutation
        if (check) { gmc->update( theface.element0(), theface.pos_first(), __p ); findPermutation=true; }
    }
    CHECK(findPermutation) << "the permutation of quad point is not find\n";

}


template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( boost::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::true_ )
{
    return gmcExpr;
}

template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( boost::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::false_ )
{
    CHECK( false ) << "not allowed\n";
    return boost::shared_ptr<GmcType>();
}

template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( boost::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt ,boost::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/ )
{
    return buildGmcWithRelationDifferentMeshType2<SpaceType,ImType,GmcType,GmcExprType>( space,gm,im, idElt, gmcExpr,
                                                                                         mpl::int_<0>(),
                                                                                         mpl::bool_<boost::is_same<GmcType,GmcExprType>::type::value>() );
}
template<typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
boost::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( boost::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt,boost::shared_ptr<GmcExprType> /**/ ,mpl::int_<1> /**/ )
{
    typedef GmcType gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gmc_type::gm_type::precompute_type geopc_type;
    typedef typename gmc_type::gm_type::precompute_ptrtype geopc_ptrtype;

    auto const& eltInit = space->mesh()->element( idElt );

    geopc_ptrtype geopc( new geopc_type( gm, im.points() ) );
    gmc_ptrtype gmc( new gmc_type( gm, eltInit, geopc ) );

    return gmc;
}

template<typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType2( FaceType const& theface, boost::shared_ptr<SpaceType> const& /*space*/,
                                         boost::shared_ptr<GmcType> /*gmc*/, boost::shared_ptr<GmcExprType> /*gmcExpr*/,
                                         size_type /*idElt*/, mpl::int_<0> /**/ )
{
    // nothing to do!
}
template<typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType2( FaceType const& theface, boost::shared_ptr<SpaceType> const& space,
                                         boost::shared_ptr<GmcType> gmc, boost::shared_ptr<GmcExprType> gmcExpr,
                                         size_type idElt, mpl::int_<1> /**/ )
{
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;

    auto const& theelt = space->mesh()->element( idElt );
    gmc->update(theelt);

    bool findPermutation=false;
    for ( permutation_type __p( permutation_type::IDENTITY );
          __p < permutation_type( permutation_type::N_PERMUTATIONS ) && !findPermutation; ++__p )
    {
        // update only xReal in gmc
        gmcExpr->update( theface.element0(), theface.pos_first(), __p, false );

        bool check=true;
        for ( uint16_type i=0;i<gmc->nPoints() && check;++i )
            for (uint16_type d=0;d<GmcType::NDim;++d)
                check = check && ( std::abs(gmc->xReal(i)[d]-gmcExpr->xReal(i)[d])<1e-8 );

        // if check compute full gmc context with the good permutation
        if (check) { gmcExpr->update( theface.element0(), theface.pos_first(), __p ); findPermutation=true; }
    }
    CHECK(findPermutation) << "the permutation of quad point is not find\n";
}

} // namespace detail

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                           mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const
{
    LOG(INFO) << "[integrator::assembleWithRelationDifferentMeshType] vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS>\n";

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc1_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef boost::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::gm1_1_type gm1_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc1_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef boost::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc1_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef boost::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //static const bool gmTestIsGmExpr = boost::is_same<gm_expr_type,gm_formTest_type>::type::value;
    //static const bool gmTrialIsGmExpr = boost::is_same<gm_expr_type,gm_formTrial_type>::type::value;
    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    static const uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_range_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtest_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtest_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtrial_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtrial_type;

    im_range_type imRange;
    im1_range_type im1Range;
    im_formtest_type imTest;
    im_formtrial_type imTrial;


    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, im1_range_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef form_context_type fcb_type;
    typedef form1_context_type fcb1_type;
    typedef fcb_type* focb_ptrtype;
    typedef fcb1_type* focb1_ptrtype;

    //-----------------------------------------------//

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
        {
            auto elt_it = lit->template get<1>();
            auto const elt_en = lit->template get<2>();
            DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

            // check that we have elements to iterate over
            if ( elt_it == elt_en )
                continue;

            const bool rangeIsSubMeshFromTest = elt_it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
            const bool rangeIsSubMeshFromTrial = elt_it->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
            const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );
            const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );

            const size_type idEltRangeInit = elt_it->id();
            size_type idEltTestInit  = idEltRangeInit;
            if ( rangeIsSubMeshFromTest )
                idEltTestInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
            else if ( testIsSubMeshFromRange )
                idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            size_type idEltTrialInit = idEltRangeInit;
            if ( rangeIsSubMeshFromTrial )
                idEltTrialInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
            else if ( trialIsSubMeshFromRange )
                idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            //-----------------------------------------------//
            pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), imRange.points() ) );
            pc1_expr_ptrtype geopc1Expr( new pc1_expr_type( elt_it->gm1(), im1Range.points() ) );
            gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->gm(),*elt_it, geopcExpr ) );
            gmc1_expr_ptrtype gmc1Expr( new gmc1_expr_type( elt_it->gm1(),*elt_it, geopc1Expr ) );
            //-----------------------------------------------//
            auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType< typename FormType::space_1_type,im_formtest_type,
                                                                              gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                                idEltTestInit, gmcExpr, mpl::int_<gmTestRangeRelation>() );
            auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType< typename FormType::space_2_type,im_formtrial_type,
                                                                               gmc_formTrial_type,gmc_expr_type>( __form.trialSpace(), __form.trialSpace()->gm(), imTrial,
                                                                                                                  idEltTrialInit, gmcExpr, mpl::int_<gmTrialRangeRelation>() );
            //-----------------------------------------------//
            map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
            map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
            map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
            //-----------------------------------------------//
            focb_ptrtype formc( new form_context_type( __form,
                                                       mapgmcFormTest,
                                                       mapgmcFormTrial,
                                                       mapgmcExpr,
                                                       this->expression(),
                                                       imRange,imTest,imTrial ) );

            //-----------------------------------------------//

            for ( ; elt_it != elt_en; ++elt_it )
                {
                    const size_type idEltRange = elt_it->id();
                    size_type idEltTest = idEltRange;
                    if ( rangeIsSubMeshFromTest )
                        idEltTest = elt_it->mesh()->subMeshToMesh( idEltRange );
                    else if ( testIsSubMeshFromRange )
                        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
                    size_type idEltTrial = idEltRange;
                    if ( rangeIsSubMeshFromTrial )
                        idEltTrial = elt_it->mesh()->subMeshToMesh( idEltRange );
                    else if ( trialIsSubMeshFromRange )
                        idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltRange );

#if 0
                    if ( formc->isZero( eltTest.element0() /*idElt*/ ) )
                        continue;
#endif

                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        gmcExpr->update(*elt_it);
                        detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_1_type,im_formtest_type,
                                                                       gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                        idEltTest, mpl::int_<gmTestRangeRelation>() );
                        detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_2_type,im_formtrial_type,
                                                                       gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                         idEltTrial, mpl::int_<gmTrialRangeRelation>() );

                        formc->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                        formc->integrate();
                        formc->assemble();
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_O1:
                    {
                        CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_ELEMENTS>] : GEOMAP_O1 not implement\n";
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_OPT:
                    {
                        if ( ( rangeIsSubMeshFromTest && __form.testSpace()->mesh()->face( idEltTest ).isOnBoundary() ) ||
                             ( rangeIsSubMeshFromTrial && __form.trialSpace()->mesh()->face( idEltTrial ).isOnBoundary() ) )
                        {
                            gmcExpr->update(*elt_it);
                            detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                            idEltTest, mpl::int_<gmTestRangeRelation>() );
                            detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_2_type,im_formtrial_type,
                                                                           gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                             idEltTrial, mpl::int_<gmTrialRangeRelation>() );

                            formc->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                            formc->integrate();
                            formc->assemble();
                        }
                        else
                        {
                            CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_ELEMENTS>] : GEOMAP_O1 not implement\n";
                        }
                    }
                    break;

                    }
                } // end loop on elements
            delete formc;
            //delete formc1;
        } // end loop on list of elements
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                           mpl::int_<MESH_ELEMENTS> /**/, mpl::true_ ) const
{
    LOG(INFO) << "[integrator::assembleWithRelationDifferentMeshType(ELEMENTS)] mortar case\n";

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc1_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef boost::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::gm1_1_type gm1_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc1_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef boost::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc1_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef boost::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //static const bool gmTestIsGmExpr = boost::is_same<gm_expr_type,gm_formTest_type>::type::value;
    //static const bool gmTrialIsGmExpr = boost::is_same<gm_expr_type,gm_formTrial_type>::type::value;
    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    static const uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_range_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtest_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtest_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtrial_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtrial_type;

    im_range_type imRange;
    im1_range_type im1Range;
    im_formtest_type imTest;
    im_formtrial_type imTrial;


    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type,true> form_mortar_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, im1_range_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, im1_range_type, map_gmc1_expr_type, map_gmc1_formTrial_type,true> form1_mortar_context_type;
    typedef form_context_type fcb_type;
    typedef form_mortar_context_type mortar_fcb_type;
    typedef form1_context_type fcb1_type;
    typedef form1_mortar_context_type mortar_fcb1_type;
    typedef fcb_type* focb_ptrtype;
    typedef mortar_fcb_type* mortar_focb_ptrtype;
    typedef fcb1_type* focb1_ptrtype;
    typedef mortar_fcb1_type* mortar_focb1_ptrtype;

    //-----------------------------------------------//

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
        {
            auto elt_it = lit->template get<1>();
            auto const elt_en = lit->template get<2>();
            DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

            // check that we have elements to iterate over
            if ( elt_it == elt_en )
                continue;

            const bool rangeIsSubMeshFromTest = elt_it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
            const bool rangeIsSubMeshFromTrial = elt_it->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
            const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );
            const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );

            const size_type idEltRangeInit = elt_it->id();
            size_type idEltTestInit  = idEltRangeInit;
            if ( rangeIsSubMeshFromTest )
                idEltTestInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
            else if ( testIsSubMeshFromRange )
                idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            size_type idEltTrialInit = idEltRangeInit;
            if ( rangeIsSubMeshFromTrial )
                idEltTrialInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
            else if ( trialIsSubMeshFromRange )
                idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            //-----------------------------------------------//
            pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), imRange.points() ) );
            pc1_expr_ptrtype geopc1Expr( new pc1_expr_type( elt_it->gm1(), im1Range.points() ) );
            gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->gm(),*elt_it, geopcExpr ) );
            gmc1_expr_ptrtype gmc1Expr( new gmc1_expr_type( elt_it->gm1(),*elt_it, geopc1Expr ) );
            //-----------------------------------------------//
            auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType< typename FormType::space_1_type,im_formtest_type,
                                                                              gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                                idEltTestInit, gmcExpr, mpl::int_<gmTestRangeRelation>() );
            auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType< typename FormType::space_2_type,im_formtrial_type,
                                                                               gmc_formTrial_type,gmc_expr_type>( __form.trialSpace(), __form.trialSpace()->gm(), imTrial,
                                                                                                                  idEltTrialInit, gmcExpr, mpl::int_<gmTrialRangeRelation>() );
            //-----------------------------------------------//
            map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
            map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
            map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
            //-----------------------------------------------//
            focb_ptrtype formc( new form_context_type( __form,
                                                       mapgmcFormTest,
                                                       mapgmcFormTrial,
                                                       mapgmcExpr,
                                                       this->expression(),
                                                       imRange,imTest,imTrial ) );
            mortar_focb_ptrtype formcm( new form_mortar_context_type( __form,
                                                                      mapgmcFormTest,
                                                                      mapgmcFormTrial,
                                                                      mapgmcExpr,
                                                                      this->expression(),
                                                                      imRange,imTest,imTrial ) );

            //-----------------------------------------------//

            for ( ; elt_it != elt_en; ++elt_it )
                {
                    const size_type idEltRange = elt_it->id();
                    size_type idEltTest = idEltRange;
                    if ( rangeIsSubMeshFromTest )
                        idEltTest = elt_it->mesh()->subMeshToMesh( idEltRange );
                    else if ( testIsSubMeshFromRange )
                        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
                    size_type idEltTrial = idEltRange;
                    if ( rangeIsSubMeshFromTrial )
                        idEltTrial = elt_it->mesh()->subMeshToMesh( idEltRange );
                    else if ( trialIsSubMeshFromRange )
                        idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltRange );

#if 0
                    if ( formc->isZero( eltTest.element0() /*idElt*/ ) )
                        continue;
#endif
                    LOG(INFO) << "elt "  << elt_it->id() << " bdy = " << elt_it->isOnBoundary();
                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        gmcExpr->update(*elt_it);
                        detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_1_type,im_formtest_type,
                                                                       gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                        idEltTest, mpl::int_<gmTestRangeRelation>() );
                        detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_2_type,im_formtrial_type,
                                                                       gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                         idEltTrial, mpl::int_<gmTrialRangeRelation>() );
                        if ( elt_it->isOnBoundary()  )
                        {
                            formcm->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                            formcm->integrate();
                            formcm->assemble();
                        }
                        else
                        {
                            formc->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                            formc->integrate();
                            formc->assemble();
                        }
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_O1:
                    {
                        CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_ELEMENTS>] : GEOMAP_O1 not implement\n";
                    }
                    break;

                    case GeomapStrategyType::GEOMAP_OPT:
                    {
                        if ( ( rangeIsSubMeshFromTest && __form.testSpace()->mesh()->face( idEltTest ).isOnBoundary() ) ||
                             ( rangeIsSubMeshFromTrial && __form.trialSpace()->mesh()->face( idEltTrial ).isOnBoundary() ) )
                        {
                            gmcExpr->update(*elt_it);
                            detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                            idEltTest, mpl::int_<gmTestRangeRelation>() );
                            detail::updateGmcWithRelationDifferentMeshType<typename FormType::space_2_type,im_formtrial_type,
                                                                           gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                             idEltTrial, mpl::int_<gmTrialRangeRelation>() );

                            if ( elt_it->isOnBoundary()  )
                            {
                                formcm->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                                formcm->integrate();
                                formcm->assemble();
                            }
                            else
                            {
                                formc->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr );
                                formc->integrate();
                                formc->assemble();
                            }
                        }
                        else
                        {
                            CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_ELEMENTS>] : GEOMAP_O1 not implement\n";
                        }
                    }
                    break;

                    }
                } // end loop on elements
            delete formc;
            //delete formc1;
        } // end loop on list of elements
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
{
    CHECK ( false ) << "[assembleWithRelationDifferentMeshType<LinearForm,MESH_ELEMENTS>] : not implement\n";
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                 mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const
{
    LOG(INFO) << "[integrator::assembleInCaseOfInterpolate] vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS>\n";

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;


    //-----------------------------------------------//
    im_range_type imRange;
    pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->gm(),*elt_it, geopcExpr ) );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
    //-----------------------------------------------//
    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
    //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               imRange ) );

    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest,meshTrial );
    }

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->result().begin();
        auto const res_en = M_QPL->result().end();
        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& map = res_it->template get<1>();
            auto map_it = map.begin();
            auto const map_en = map.end();
            for ( ; map_it != map_en ; ++map_it )
            {
                auto const idEltTrial = map_it->first;
                auto const& eltTrial = meshTrial->element( idEltTrial );
                auto const& eltTest = meshTest->element( idEltTest );

                auto const& ptRefTest = map_it->second.template get<1>();
                auto const& ptRefTrial = map_it->second.template get<2>();
                auto const& themapQuad = map_it->second.template get<0>();
                auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
                auto gmcExpr_it = vec_gmcExpr.begin();
                auto const gmcExpr_en = vec_gmcExpr.end();
                bool isFirstExperience = true;
                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcFormTest->update( gmcExpr_it->template get<2>() );
                    geopcFormTrial->update( gmcExpr_it->template get<3>() );

                    gmcFormTest->update( eltTest,geopcFormTest );
                    gmcFormTrial->update( eltTrial,geopcFormTrial );

                    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );

                    formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->template get<0>() );

                    formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                    isFirstExperience = false;
                }

                formc->assembleInCaseOfInterpolate();

            }
        }
    } // if (!M_QPL->hasPrecomputeBF())
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcFormTest*/ ) );
                map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<3>()/*gmcFormTrial*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );

                formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assembleInCaseOfInterpolate();
        }

    } //else !M_QPL->hasPrecomputeBF()

    delete formc;

} // assembleInCaseOfInterpolate

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                 mpl::int_<MESH_ELEMENTS> /**/, mpl::true_ ) const
{
    LOG(INFO) << "[integrator::assembleInCaseOfInterpolate] (elements) mortar case";

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type,true> form_mortar_context_type;
    typedef form_context_type fcb_type;
    typedef form_mortar_context_type mortar_fcb_type;
    typedef fcb_type* focb_ptrtype;
    typedef mortar_fcb_type* mortar_focb_ptrtype;

    //-----------------------------------------------//

    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;


    //-----------------------------------------------//
    im_range_type imRange;
    pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->gm(),*elt_it, geopcExpr ) );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
    //-----------------------------------------------//
    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), __form.template testFiniteElement<false>()->points() ) );
    VLOG(2) << "pts non mortar : " << __form.template testFiniteElement<false>()->points();
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );

    // mortar
    pc_formTest_ptrtype geopcFormTestMortar( new pc_formTest_type( __form.gm(), __form.template testFiniteElement<true>()->points() ) );
    VLOG(2) << "pts mortar : " << __form.template testFiniteElement<true>()->points();
    gmc_formTest_ptrtype gmcFormTestMortar( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTestMortar ) );
    map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );
    //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               imRange ) );
    mortar_focb_ptrtype formcm( new form_mortar_context_type( __form,
                                                              mapgmcFormTestMortar,
                                                              mapgmcFormTrial,
                                                              mapgmcExpr,
                                                              this->expression(),
                                                              imRange ) );

    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest,meshTrial );
    }

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->result().begin();
        auto const res_en = M_QPL->result().end();
        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& map = res_it->template get<1>();
            auto map_it = map.begin();
            auto const map_en = map.end();

            for ( ; map_it != map_en ; ++map_it )
            {
                auto const idEltTrial = map_it->first;
                auto const& eltTrial = meshTrial->element( idEltTrial );
                auto const& eltTest = meshTest->element( idEltTest );

                auto const& ptRefTest = map_it->second.template get<1>();
                auto const& ptRefTrial = map_it->second.template get<2>();
                auto const& themapQuad = map_it->second.template get<0>();
                auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
                auto gmcExpr_it = vec_gmcExpr.begin();
                auto const gmcExpr_en = vec_gmcExpr.end();
                bool isFirstExperience = true;
                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcFormTrial->update( gmcExpr_it->template get<3>() );
                    gmcFormTrial->update( eltTrial,geopcFormTrial );
                    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );

                    if ( eltTest.isOnBoundary() )
                    {
                        geopcFormTestMortar->update( gmcExpr_it->template get<2>() );
                        gmcFormTestMortar->update( eltTest,geopcFormTestMortar );

                        map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );

                        formcm->updateInCaseOfInterpolate( mapgmcFormTestMortar, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->template get<0>() );
                        formcm->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                    }
                    else
                    {
                        geopcFormTest->update( gmcExpr_it->template get<2>() );
                        gmcFormTest->update( eltTest,geopcFormTest );

                        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );

                        formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->template get<0>() );
                        formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                    }
                    isFirstExperience = false;
                }
                if ( eltTest.isOnBoundary() )
                {
                    formcm->assembleInCaseOfInterpolate();
                }
                else
                {
                    formc->assembleInCaseOfInterpolate();
                }

            }
        }
    } // if (!M_QPL->hasPrecomputeBF())
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcFormTest*/ ) );
                map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<3>()/*gmcFormTrial*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );

                formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assembleInCaseOfInterpolate();
        }

    } //else !M_QPL->hasPrecomputeBF()

    delete formc;

} // assembleInCaseOfInterpolate


template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
{
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT|vm::JACOBIAN,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_form_type, expression_type, im_range_type, map_gmc_expr_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;

    //-----------------------------------------------//

    im_range_type imRange;
    pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->gm(),*elt_it, geopcExpr ) );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcForm ) );
    map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( gmcForm ) );

    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcForm,
                                               mapgmcForm,
                                               mapgmcExpr,
                                               this->expression(),
                                               imRange ) );

    //-----------------------------------------------//

    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest );
    }

    //-----------------------------------------------//

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->resultLinear().begin();
        auto const res_en = M_QPL->resultLinear().end();

        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& eltTest = meshTest->element( idEltTest );
            auto const& ptRefTest = res_it->template get<2>();
            auto const& themapQuad = res_it->template get<1>();
            auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest );
            auto gmcExpr_it = vec_gmcExpr.begin();
            auto const gmcExpr_en = vec_gmcExpr.end();
            bool isFirstExperience = true;

            for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
            {
                geopcForm->update( gmcExpr_it->template get<2>() );
                gmcForm->update( eltTest,geopcForm );
                map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( gmcForm ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,gmcExpr_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assemble();
        }
    } //if (!M_QPL->hasPrecompute())
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcForm*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );
                formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }
            formc->assemble();
        }
    } //else

    delete formc;

}


template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
bool
Integrator<Elements, Im, Expr, Im2>::useSameMesh( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                  FaceRangeType const& faceRange ) const
{
    bool res = faceRange.mesh()->isSameMesh( __form.testSpace()->mesh() ) && faceRange.mesh()->isSameMesh( __form.trialSpace()->mesh() );
    return res;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
bool
Integrator<Elements, Im, Expr, Im2>::useSameMesh( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                                                  FaceRangeType const& faceRange ) const
{
    bool res = faceRange.mesh()->isSameMesh( __form.testSpace()->mesh() );
    return res;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
boost::tuple<size_type,rank_type,uint16_type>
Integrator<Elements, Im, Expr, Im2>::testElt0IdFromFaceRange( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                              FaceRangeType const& faceRange ) const
{
    uint16_type __face_id_in_elt_0 = faceRange.pos_first();
    rank_type procIdElt0 = faceRange.proc_first();
    //uint16_type idEltConnectedFaceRange = 0;
    size_type idEltTest = faceRange.element( 0 ).id();
    bool trialEltIsOk = false;

    if ( !faceRange.element0().isGhostCell() )
    {


    if ( faceRange.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
    {
        idEltTest = faceRange.mesh()->subMeshToMesh( idEltTest );
    }
    else if ( __form.testSpace()->mesh()->isSubMeshFrom( faceRange.mesh() ) )
    {
        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltTest );
        if ( idEltTest == invalid_size_type_value )
        {
            if ( faceRange.isConnectedTo1() && !faceRange.element1().isGhostCell() )
            {
                __face_id_in_elt_0 = faceRange.pos_second();
                procIdElt0 = faceRange.proc_second();
                //idEltConnectedFaceRange = 1;
                idEltTest = __form.testSpace()->mesh()->meshToSubMesh( faceRange.element( 1 ).id() );
            }
        }
    }

    if ( idEltTest != invalid_size_type_value )
    {
        size_type idEltTrial = idEltTest;
        if ( __form.trialSpace()->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
            idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltTest );
        else if ( __form.testSpace()->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() ) )
            idEltTrial = __form.testSpace()->mesh()->subMeshToMesh( idEltTest );

        if (idEltTrial != invalid_size_type_value)
            trialEltIsOk=true;
     }

    }
    else // ghost cell
    {
        if ( !faceRange.isConnectedTo1() || faceRange.element1().isGhostCell() )
            return boost::make_tuple( invalid_size_type_value, procIdElt0, __face_id_in_elt_0 );
        else
            idEltTest = invalid_size_type_value; // continue algo
    }

    // if first pass not good, restart with faceRange.element1
    if ( ( idEltTest == invalid_size_type_value || !trialEltIsOk )  && faceRange.isConnectedTo1() && !faceRange.element1().isGhostCell() )
    {
        __face_id_in_elt_0 = faceRange.pos_second();
        procIdElt0 = faceRange.proc_second();
        idEltTest = faceRange.element( 1 ).id();
        if ( faceRange.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
            idEltTest = faceRange.mesh()->subMeshToMesh( idEltTest );
        else if ( __form.testSpace()->mesh()->isSubMeshFrom( faceRange.mesh() ) )
            idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltTest );

        trialEltIsOk = false;
        if ( idEltTest != invalid_size_type_value )
        {
            size_type idEltTrial = idEltTest;
            if ( __form.trialSpace()->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltTest );
            else if ( __form.testSpace()->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() ) )
                idEltTrial = __form.testSpace()->mesh()->subMeshToMesh( idEltTest );

            if (idEltTrial != invalid_size_type_value)
                {

                trialEltIsOk=true;
                }
        }
    }

    // if test or trial id not find, return invalid value
    if ( idEltTest == invalid_size_type_value || !trialEltIsOk )
        return boost::make_tuple( invalid_size_type_value, procIdElt0, __face_id_in_elt_0 );
    else
        return boost::make_tuple( idEltTest,procIdElt0, __face_id_in_elt_0 );
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
boost::tuple<size_type,rank_type,uint16_type>
Integrator<Elements, Im, Expr, Im2>::testElt0IdFromFaceRange( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                                                              FaceRangeType const& faceRange ) const
{
    uint16_type __face_id_in_elt_0 = faceRange.pos_first();
    rank_type procIdElt0 = faceRange.proc_first();
    size_type idEltTest = faceRange.element( 0 ).id();
    if ( faceRange.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
    {
        idEltTest = faceRange.mesh()->subMeshToMesh( idEltTest );
    }
    else if ( __form.testSpace()->mesh()->isSubMeshFrom( faceRange.mesh() ) )
    {
        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltTest );
        if ( idEltTest == invalid_size_type_value )
        {
            if ( faceRange.isConnectedTo1() )
            {
                __face_id_in_elt_0 = faceRange.pos_second();
                procIdElt0 = faceRange.proc_second();
                idEltTest = __form.testSpace()->mesh()->meshToSubMesh( faceRange.element( 1 ).id() );
            }
        }
    }
    return boost::make_tuple( idEltTest, procIdElt0, __face_id_in_elt_0 );
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType,typename FaceRangeType>
bool
Integrator<Elements, Im, Expr, Im2>::faceIntegratorUseTwoConnections( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                      FaceRangeType const& faceRange,
                                                                      typename vf::detail::BilinearForm<FE1,FE2,ElemContType>::mesh_1_type::element_type const& eltTest,
                                                                      uint16_type faceIdInElt ) const
{
    bool res = faceRange.isConnectedTo0() && faceRange.isConnectedTo1() &&
        eltTest.face(faceIdInElt).isConnectedTo0() && eltTest.face(faceIdInElt).isConnectedTo1();

    // if possible, search corresponding id in trial space
    if ( res )
    {
        size_type idEltTrial = eltTest.id();
        if ( __form.trialSpace()->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
            idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltTrial );
        else if ( __form.testSpace()->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() ) )
            idEltTrial = __form.testSpace()->mesh()->subMeshToMesh( idEltTrial );

        if ( idEltTrial == invalid_size_type_value )
            res=false;
        else
        {
            auto const& faceTrial = __form.trialSpace()->mesh()->element(idEltTrial, eltTest.processId() ).face( faceIdInElt );
            res = faceTrial.isConnectedTo0() && faceTrial.isConnectedTo1();
        }
    }

    return res;
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
bool
Integrator<Elements, Im, Expr, Im2>::faceIntegratorUseTwoConnections( vf::detail::LinearForm<FE,VectorType,ElemContType>& __form,
                                                                      FaceRangeType const& faceRange,
                                                                      typename vf::detail::LinearForm<FE,VectorType,ElemContType>::mesh_type::element_type const& eltTest,
                                                                      uint16_type faceIdInElt ) const
{
    bool res = faceRange.isConnectedTo0() && faceRange.isConnectedTo1() &&
        eltTest.face(faceIdInElt).isConnectedTo0() && eltTest.face(faceIdInElt).isConnectedTo1();
    return res;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<true> /**/, bool /*hasRelation*/ ) const
{
    //DLOG(INFO) << "integrating over "
    //              << std::distance( this->beginElement(), this->endElement() )  << " faces\n";
    boost::timer __timer;

    //
    // some typedefs
    //
    typedef typename eval::gm_type gm_type;
    typedef typename eval::gm1_type gm1_type;
    //typedef typename FormType::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef boost::shared_ptr<gm1_type> gm1_ptrtype;
    typedef typename gm_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_type;
    typedef typename gm1_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc1_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;
    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename eval::gmpc_type pc_type;
    typedef typename eval::gmpc1_type pc1_type;
    typedef typename eval::gmpc_ptrtype pc_ptrtype;
    typedef typename eval::gmpc1_ptrtype pc1_ptrtype;
    //typedef typename mpl::if_<mpl::equal_to<mpl::int_<FormType::nDim>, mpl::int_<2> >, mpl::identity<typename eval::element_type::edge_permutation_type>, mpl::identity<typename eval::element_type::face_permutation_type> >::type::type permutation_type;

    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;
    //typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );


    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( this->im().nFaces() );
    std::vector<std::map<permutation_type, pc1_ptrtype> > __geopc1( this->im2().nFaces() );
    typedef typename im_type::face_quadrature_type face_im_type;
    typedef typename im2_type::face_quadrature_type face_im2_type;

    //__typeof__(im(__face_id_in_elt_0 ) ) im_face ( im(__face_id_in_elt_0 ) );
    std::vector<face_im_type> face_ims( im().nFaces() );
    std::vector<face_im2_type> face_ims2( im2().nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = this->im( __f );
        face_ims2[__f] = this->im2( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopc[__f][__p] = pc_ptrtype(  new pc_type( __form.gm(), this->im().fpoints(__f, __p.value() ) ) );
            __geopc1[__f][__p] = pc1_ptrtype(  new pc1_type( __form.gm1(), this->im2().fpoints(__f, __p.value() ) ) );
        }
    }

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        DLOG(INFO) << "integrating over "
                   << std::distance( it, en )  << " faces\n";

        // check that we have elements to iterate over
        if ( (it == en) || (it->isConnectedTo0() == false) )
            continue;

        // true if range/test/trial are same mesh
        bool useSameMesh = this->useSameMesh( __form,*it );

        uint16_type __face_id_in_elt_0 = it->pos_first();
        rank_type procIdElt0 = it->proc_first();
        size_type idEltTestInit = it->element( 0 ).id();
        if ( !useSameMesh )
        {
            bool hasFindEltToInit = false;
            while( !hasFindEltToInit)
            {
                if ( it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                {
                    idEltTestInit = it->mesh()->subMeshToMesh( idEltTestInit );
                }
                else if ( __form.testSpace()->mesh()->isSubMeshFrom( it->mesh() ) )
                {
                    idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltTestInit );
                    if ( idEltTestInit == invalid_size_type_value )
                    {
                        if ( it->isConnectedTo1() )
                        {
                            __face_id_in_elt_0 = it->pos_second();
                            procIdElt0 = it->proc_second();
                            idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( it->element( 1 ).id() );
                        }
                    }
                }
                if ( idEltTestInit == invalid_size_type_value )
                {
                    ++it;
                    if (it==en) break;
                }
                else
                    hasFindEltToInit=true;
            }
            if ( !hasFindEltToInit )
                continue;
            else
                CHECK( idEltTestInit != invalid_size_type_value ) << "mesh relation fail : no find a corresponding element\n";

        }

        auto const& elt0TestInit = __form.testSpace()->mesh()->element( idEltTestInit,procIdElt0 );
        //auto const& faceTestInit = elt0TestInit.face( __face_id_in_elt_0 );

        // get the geometric mapping associated with element 0
        //DLOG(INFO) << "element " << it->element(0)  << "face " << __face_id_in_elt_0 << " permutation " << it->element(0).permutation( __face_id_in_elt_0 ) << "\n";
        gm_ptrtype __gm = elt0TestInit.gm();
        gm1_ptrtype __gm1 = elt0TestInit.gm1();
        //DLOG(INFO) << "[integrator] evaluate(faces), gm is cached: " << __gm->isCached() << "\n";
        gmc_ptrtype __c0( new gmc_type( __gm, elt0TestInit, __geopc, __face_id_in_elt_0 ) );
        gmc1_ptrtype __c01( new gmc1_type( __gm1, elt0TestInit, __geopc1, __face_id_in_elt_0 ) );


        //
        // the case where the face is connected only to one element
        //
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
        typedef typename FormType::template Context<map_gmc_type, expression_type, face_im_type> form_context_type;
        typedef typename FormType::template Context<map_gmc1_type, expression_type, face_im2_type> form1_context_type;
        typedef boost::shared_ptr<form_context_type> form_context_ptrtype;
        typedef boost::shared_ptr<form1_context_type> form1_context_ptrtype;
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
        map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c01 ) );
        form_context_ptrtype form;
        form1_context_ptrtype form1;

        //
        // the case where the face is connected only to two elements
        //
        // get the geometric mapping associated with element 1
        gmc_ptrtype __c1;
        gmc1_ptrtype __c11;

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc1_ptrtype> > map21_gmc_type;
        typedef typename FormType::template Context<map2_gmc_type, expression_type, face_im_type> form2_context_type;
        typedef typename FormType::template Context<map21_gmc_type, expression_type, face_im2_type> form21_context_type;
        typedef boost::shared_ptr<form2_context_type> form2_context_ptrtype;
        typedef boost::shared_ptr<form21_context_type> form21_context_ptrtype;
        form2_context_ptrtype form2;
        form21_context_ptrtype form21;

        bool isInitConnectionTo0=false;
        bool isInitConnectionTo1=false;

        // true if connected to another element, false otherwise
        //if ( it->isConnectedTo1() )
        if ( this->faceIntegratorUseTwoConnections(__form,*it,elt0TestInit,__face_id_in_elt_0) )
        {
            uint16_type __face_id_in_elt_1 = it->pos_second();
            rank_type procIdElt1 = it->proc_second();
            size_type idElt1TestInit = it->element( 1 ).id();
            if ( !useSameMesh )
            {
                // search other connection
                if ( elt0TestInit.face( __face_id_in_elt_0 ).element( 0 ).id() == elt0TestInit.id() )
                {
                    __face_id_in_elt_1 = elt0TestInit.face( __face_id_in_elt_0 ).pos_second();
                    procIdElt1 = elt0TestInit.face( __face_id_in_elt_0 ).proc_second();
                    idElt1TestInit = elt0TestInit.face( __face_id_in_elt_0 ).element( 1 ).id();
                }
                else
                {
                    __face_id_in_elt_1 = elt0TestInit.face( __face_id_in_elt_0 ).pos_first();
                    procIdElt1 = elt0TestInit.face( __face_id_in_elt_0 ).proc_first();
                    idElt1TestInit = elt0TestInit.face( __face_id_in_elt_0 ).element( 0 ).id();
                }
            }
            //CHECK( idElt1TestInit != invalid_size_type_value ) << "mesh relation fail : no find a corresponding element\n";
            // get element1
            auto const& elt1TestInit = __form.testSpace()->mesh()->element( idElt1TestInit,procIdElt1 );

            // init linear/bilinear form for two connections
            __c1 = gmc_ptrtype( new gmc_type( __gm, elt1TestInit, __geopc, __face_id_in_elt_1 ) );
            __c11 = gmc1_ptrtype( new gmc1_type( __gm1, elt1TestInit, __geopc1, __face_id_in_elt_1 ) );
            map2_gmc_type mapgmc2( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                   fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );
            map21_gmc_type mapgmc21( fusion::make_pair<vf::detail::gmc<0> >( __c01 ),
                                     fusion::make_pair<vf::detail::gmc<1> >( __c11 ) );

            form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
            form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
            isInitConnectionTo1=true;
        }

        else
        {
            // init linear/bilinear form for one connection
            form = form_context_ptrtype( new form_context_type( __form, mapgmc, mapgmc, mapgmc, expression(), face_ims[__face_id_in_elt_0], this->im() ) );
            form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
            isInitConnectionTo0=true;
        }

        boost::timer ti0,ti1, ti2, ti3;
        //double t0 = 0, t1 = 0,t2 = 0,t3 = 0;
        DLOG(INFO) << "[Integrator::faces/forms] starting...\n";


        //
        // start the real intensive job:
        // -# iterate over all elements to integrate over
        // -# construct the associated geometric mapping with the reference element
        // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
        // -# assemble the local contribution in the global representation of the bilinear form
        //
        for ( ; it != en; ++it )
        {
            // get some info about test mesh element0 connected to the current face
            uint16_type __face_id_in_elt_0 = it->pos_first();
            rank_type procIdElt0 = it->proc_first();
            size_type idEltTest = it->element( 0 ).id();
            bool swapElt0Elt1WithSameMesh = false;
            if ( !useSameMesh)
            {
                boost::tie( idEltTest, procIdElt0, __face_id_in_elt_0) = this->testElt0IdFromFaceRange(__form,*it);
                if ( idEltTest == invalid_size_type_value )
                    continue;
            }
            else if ( it->element0().isGhostCell() && it->isConnectedTo1() )
            {
                __face_id_in_elt_0 = it->pos_second();
                procIdElt0 = it->proc_second();
                idEltTest = it->element( 1 ).id();
                swapElt0Elt1WithSameMesh = true;
            }

            // element0 (from test mesh) used in integration
            auto const& elt0Test = __form.testSpace()->mesh()->element( idEltTest,procIdElt0 );
            CHECK( !elt0Test.isGhostCell() ) << "elt0 can't be a ghost element";
            //auto const& faceTest = elt0Test.face( __face_id_in_elt_0 );


            //if ( it->isConnectedTo1() )
            if ( this->faceIntegratorUseTwoConnections(__form,*it,elt0Test,__face_id_in_elt_0) )
            {

                if ( it->isGhostFace() )
                {
                    LOG(WARNING) << "face id : " << it->id() << " is a ghost face" << it->G();
                    continue;
                }
                // if is a interprocess faces, only integrate in one process
                if ( it->isInterProcessDomain() && it->partition1() > it->partition2() )
                    continue;


                // get some info about test mesh element1 connected to the current face
                uint16_type __face_id_in_elt_1 = it->pos_second();
                rank_type procIdElt1 = it->proc_second();
                size_type idElt1Test = it->element( 1 ).id();
                if ( !useSameMesh )
                {
                    // search other connection
                    if ( elt0Test.face( __face_id_in_elt_0 ).element( 0 ).id() == elt0Test.id() )
                    {
                        __face_id_in_elt_1 = elt0Test.face( __face_id_in_elt_0 ).pos_second();
                        procIdElt1 = elt0Test.face( __face_id_in_elt_0 ).proc_second();
                        idElt1Test = elt0Test.face( __face_id_in_elt_0 ).element( 1 ).id();
                    }
                    else
                    {
                        __face_id_in_elt_1 = elt0Test.face( __face_id_in_elt_0 ).pos_first();
                        procIdElt1 = elt0Test.face( __face_id_in_elt_0 ).proc_first();
                        idElt1Test = elt0Test.face( __face_id_in_elt_0 ).element( 0 ).id();
                    }
                }
                else if ( swapElt0Elt1WithSameMesh )
                {
                    uint16_type __face_id_in_elt_1 = it->pos_first();
                    rank_type procIdElt1 = it->proc_first();
                    size_type idElt1Test = it->element( 0 ).id();
                }
                CHECK( idElt1Test != invalid_size_type_value ) << "mesh relation fail : no find a corresponding element\n";

                // element1 (from test mesh) used in integration
                auto const& elt1Test = __form.testSpace()->mesh()->element( idElt1Test,procIdElt1 );

                if ( !isInitConnectionTo1 )
                {
                    // init linear/bilinear form for element1
                    __c1 = gmc_ptrtype( new gmc_type( __gm, elt1Test, __geopc, __face_id_in_elt_1 ) );
                    __c11 = gmc1_ptrtype( new gmc1_type( __gm1, elt1Test, __geopc1, __face_id_in_elt_1 ) );
                    map2_gmc_type mapgmc2( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                           fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );
                    map21_gmc_type mapgmc21( fusion::make_pair<vf::detail::gmc<0> >( __c01 ),
                                             fusion::make_pair<vf::detail::gmc<1> >( __c11 ) );

                    form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
                    form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
                    isInitConnectionTo1=true;
                }

                switch ( M_gt )
                {
                default:
                case GeomapStrategyType::GEOMAP_HO:
                {
                    FEELPP_ASSERT( it->isOnBoundary() == false  )
                        ( it->id() ).error( "face on boundary but connected on both sides" );
                    //ti0.restart();
                    __c0->update( elt0Test, __face_id_in_elt_0 );
#if 0
                    __c1->update( elt1Test, __face_id_in_elt_1 );
                    bool check=true;
                    for ( uint16_type i=0;i<__c0->nPoints() && check;++i )
                        for (uint16_type d=0;d<gmc1_type::NDim;++d)
                            check = check && ( std::abs(__c0->xReal(i)[d] - __c1->xReal(i)[d])<1e-8 );
                    CHECK( check ) << "quad points in each sides are not connected";
#else
                    bool findPermutation=false;
                    for ( permutation_type __p( permutation_type::IDENTITY );
                          __p < permutation_type( permutation_type::N_PERMUTATIONS ) && !findPermutation; ++__p )
                        {
                            // update only xReal in gmc
                            __c1->update( elt1Test, __face_id_in_elt_1, __p, false );

                            bool check=true;
                            for ( uint16_type i=0;i<__c0->nPoints() && check;++i )
                            {
                                //LOG(INFO) << "c0.xreal = " << __c0->xReal(i);
                                //LOG(INFO) << "c1.xreal(" << __p << ") = " << __c1->xReal(i);
                                for (uint16_type d=0;d<gmc1_type::NDim;++d)
                                    check = check && ( std::abs(__c0->xReal(i)[d] - __c1->xReal(i)[d])<1e-8 );
                            }
                            // if check compute full gmc context with the good permutation
                            if (check) { __c1->update( elt1Test, __face_id_in_elt_1, __p ); findPermutation=true; }
                         }
                    CHECK(findPermutation) << "the permutation of quad point is not find\n";
#endif
                    //t0 += ti0.elapsed();

                    //ti1.restart();
                    map2_gmc_type mapgmc2 = map2_gmc_type( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                                           fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );
                    form2->update( mapgmc2, mapgmc2, mapgmc2, face_ims[__face_id_in_elt_0], mpl::int_<2>() );
                    //t1 += ti1.elapsed();

                    //ti2.restart();
                    form2->integrate( );
                    //t2 += ti2.elapsed();

                    //ti3.restart();
                    form2->assemble( elt0Test.id(), elt1Test.id() );
                    //t3 += ti3.elapsed();
                }
                break;

                case GeomapStrategyType::GEOMAP_O1:
                case GeomapStrategyType::GEOMAP_OPT:
                {
                    FEELPP_ASSERT( it->isOnBoundary() == false  )
                        ( it->id() ).error( "face on boundary but connected on both sides" );
                    //ti0.restart();
                    __c01->update( elt0Test, __face_id_in_elt_0 );
#if 0
                    __c11->update( elt1Test, __face_id_in_elt_1 );
#else

                    bool findPermutation=false;
                    for ( permutation_type __p( permutation_type::IDENTITY );
                          __p < permutation_type( permutation_type::N_PERMUTATIONS ) && !findPermutation; ++__p )
                        {
                            // update only xReal in gmc
                            __c11->update( elt1Test, __face_id_in_elt_1, __p, false );

                            bool check=true;
                            for ( uint16_type i=0;i<__c01->nPoints() && check;++i )
                            {
                                //LOG(INFO) << "c0.xreal = " << __c01->xReal(i);
                                //LOG(INFO) << "c1.xreal(" << __p << ") = " << __c11->xReal(i);
                                for (uint16_type d=0;d<gmc1_type::NDim;++d)
                                    check = check && ( std::abs(__c01->xReal(i)[d] - __c11->xReal(i)[d])<1e-8 );
                            }

                            // if check compute full gmc context with the good permutation
                            if (check) { __c11->update( elt1Test, __face_id_in_elt_1, __p ); findPermutation=true; }
                         }
                    CHECK(findPermutation) << "the permutation of quad point is not find\n";
#endif
                    //t0 += ti0.elapsed();

                    //ti1.restart();
                    map21_gmc_type mapgmc21 = map21_gmc_type( fusion::make_pair<vf::detail::gmc<0> >( __c01 ),
                                                              fusion::make_pair<vf::detail::gmc<1> >( __c11 ) );
                    form21->update( mapgmc21, mapgmc21, mapgmc21, face_ims2[__face_id_in_elt_0], mpl::int_<2>() );
                    //t1 += ti1.elapsed();

                    //ti2.restart();
                    form21->integrate( );
                    //t2 += ti2.elapsed();

                    //ti3.restart();
                    form21->assemble( elt0Test.id(), elt1Test.id() );
                    //t3 += ti3.elapsed();
                }
                break;
                }
            }

            else
            {
#if 1
                if ( !isInitConnectionTo0 )
                {
                    form = form_context_ptrtype( new form_context_type( __form, mapgmc, mapgmc, mapgmc, expression(), face_ims[__face_id_in_elt_0], this->im() ) );
                    form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                    isInitConnectionTo0=true;
                }

                //ti0.restart();
                __c0->update( elt0Test,__face_id_in_elt_0 );
                //t0 += ti0.elapsed();

                FEELPP_ASSERT( __face_id_in_elt_0 == __c0->faceId() )
                    ( __face_id_in_elt_0 )
                    ( __c0->faceId() ).warn ( "invalid face id" );

                //ti1.restart();
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                form->update( mapgmc, mapgmc, mapgmc, face_ims[__face_id_in_elt_0] );
                //t1 += ti1.elapsed();

                //ti2.restart();
                form->integrate( );
                //t2 += ti2.elapsed();

                //ti3.restart();
                form->assemble( elt0Test.id() );
                //t3 += ti3.elapsed();
#else
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                form->update( mapgmc, mapgmc, mapgmc, face_ims[__face_id_in_elt_0] );
                form->integrate( );
#endif
            } // end loop on elements
        }
    }// end loop on list of element

#if 0
    DLOG(INFO) << "[faces] Overall integration time : " << ( t0+t1+t2+t3 ) << " per element:" << ( t0+t1+t2+t3 )/std::distance( this->beginElement(), this->endElement() ) << "for " << std::distance( this->beginElement(), this->endElement() ) << "elements\n";
    DLOG(INFO) << "[faces] Overall geometric mapping update time : " << t0 << "\n";
    DLOG(INFO) << "[faces] Overall form update time : " << t1 << "\n";
    DLOG(INFO) << "[faces] Overall local assembly time : " << t2 << "\n";
    DLOG(INFO) << "[faces] Overall global assembly time : " << t3 << "\n";
#endif
    DLOG(INFO) << "integrating over faces done in " << __timer.elapsed() << "s\n";
    //std::cout << "integrating over faces done in " << __timer.elapsed() << "s\n";
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<false> /**/, bool hasRelation ) const
{
    if ( hasRelation )
    {
        assembleWithRelationDifferentMeshType( __form,mpl::int_<MESH_FACES>() );
    }
    else
    {
        assembleInCaseOfInterpolate( __form,mpl::int_<MESH_FACES>() );
    }
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                           mpl::int_<MESH_FACES> /**/, mpl::false_ ) const
{

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc1_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef boost::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::gm1_1_type gm1_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type contextTest = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTest_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT> >::type::value;
    typedef typename gm_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc1_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef boost::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type contextTrial = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTrial_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT> >::type::value;
    typedef typename gm_formTrial_type::template Context<contextTrial,geoelement_formTrial_type> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<contextTrial,geoelement_formTrial_type> gmc1_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef boost::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //-----------------------------------------------------//

    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    static const uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_range_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtest_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtest_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtrial_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtrial_type;

    im_range_type imRange;
    im1_range_type im1Range;
    im_formtest_type imTest;
    im_formtrial_type imTrial;

    //-----------------------------------------------------//

    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    typedef typename im_range_type::face_quadrature_type face_im_type;
    typedef typename im1_range_type::face_quadrature_type face_im1_type;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopc( imRange.nFaces() );
    std::vector<std::map<permutation_type, pc1_expr_ptrtype> > __geopc1( im1Range.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );
    std::vector<face_im1_type> face_ims1( im1Range.nFaces() );
    bool hasInitGeoPc=false;

    //-----------------------------------------------------//
    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef boost::shared_ptr<form_context_type> form_context_ptrtype;
    typedef boost::shared_ptr<form1_context_type> form1_context_ptrtype;


    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto elt_it = lit->template get<1>();
        auto const elt_en = lit->template get<2>();
        DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

        // check that we have elements to iterate over
        if ( (elt_it == elt_en) || (elt_it->isConnectedTo0() == false) )
            continue;

        //-----------------------------------------------//

        if (!hasInitGeoPc)
        {
            for ( uint16_type __f = 0; __f < imRange.nFaces(); ++__f )
            {
                face_ims[__f] = imRange.face( __f );
                face_ims1[__f] = im1Range.face( __f );

                for ( permutation_type __p( permutation_type::IDENTITY );
                      __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                    __geopc[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
                    __geopc1[__f][__p] = pc1_expr_ptrtype(  new pc1_expr_type( elt_it->element( 0 ).gm1(), im1Range.fpoints(__f, __p.value() ) ) );
                }
            }
            hasInitGeoPc=true;
        }

        //-----------------------------------------------//

        const bool rangeIsSubMeshFromTest = elt_it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
        const bool rangeIsSubMeshFromTrial = elt_it->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
        const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );
        const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );

        const size_type idEltRangeInit = elt_it->id();
        size_type idEltTestInit  = idEltRangeInit;
        if ( rangeIsSubMeshFromTest )
            idEltTestInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
        else if ( testIsSubMeshFromRange )
            idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
        size_type idEltTrialInit = idEltRangeInit;
        if ( rangeIsSubMeshFromTrial )
            idEltTrialInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
        else if ( trialIsSubMeshFromRange )
            idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );

        //-----------------------------------------------//
        // get the geometric mapping associated with element 0
        uint16_type __face_id_in_elt_0 = elt_it->pos_first();
        gmc_expr_ptrtype gmcExpr0( new gmc_expr_type( elt_it->element( 0 ).gm(), elt_it->element( 0 ), __geopc, __face_id_in_elt_0 ) );
        //gmc1_expr_ptrtype gmc1Expr0( new gmc1_expr_type( elt_it->element( 0 ).gm1(), elt_it->element( 0 ), __geopc1, __face_id_in_elt_0 ) );

        auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType2< typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                             idEltTestInit, gmcExpr0, mpl::int_<gmTestRangeRelation>() );

        auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType2< typename FormType::space_2_type,im_formtrial_type,
                                                                           gmc_formTrial_type,gmc_expr_type>( __form.trialSpace(), __form.trialSpace()->gm(), imTrial,
                                                                                                              idEltTrialInit, gmcExpr0, mpl::int_<gmTrialRangeRelation>() );

        //-----------------------------------------------//

        map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr0 ) );
        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
        form_context_ptrtype form;
        form1_context_ptrtype form1;


        bool isInitConnectionTo0=false;
        bool isInitConnectionTo1=false;

        // true if connected to another element, false otherwise
        if ( elt_it->isConnectedTo1() )
        {
            CHECK( false ) << "TODO!!!!\n";
#if 0
            uint16_type __face_id_in_elt_1 = it->pos_second();

            __c1 = gmc_ptrtype( new gmc_type( __gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );
            __c11 = gmc1_ptrtype( new gmc1_type( __gm1, it->element( 1 ), __geopc1, __face_id_in_elt_1 ) );
            map2_gmc_type mapgmc2( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                   fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );
            map21_gmc_type mapgmc21( fusion::make_pair<vf::detail::gmc<0> >( __c01 ),
                                     fusion::make_pair<vf::detail::gmc<1> >( __c11 ) );

            form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
            form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
            isInitConnectionTo1=true;
#endif
        }
        else
        {
            form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                this->expression(), face_ims[__face_id_in_elt_0], imRange,imTest,imTrial ) );
            //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
            isInitConnectionTo0=true;
        }



        for ( ; elt_it != elt_en; ++elt_it )
        {
            if ( elt_it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << elt_it->id() << " is a ghost face";
                continue;
            }
            // if is a interprocess faces, only integrate in one process
            if ( elt_it->isInterProcessDomain() && elt_it->partition1() > elt_it->partition2() )
                continue;

            const size_type idEltRange = elt_it->id();
            size_type idEltTest = idEltRange;
            if ( rangeIsSubMeshFromTest )
                idEltTest = elt_it->mesh()->subMeshToMesh( idEltRange );
            else if ( testIsSubMeshFromRange )
                idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
            size_type idEltTrial = idEltRange;
            if ( rangeIsSubMeshFromTrial )
                idEltTrial = elt_it->mesh()->subMeshToMesh( idEltRange );
            else if ( trialIsSubMeshFromRange )
                idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltRange );


            if ( elt_it->isConnectedTo1() )
            {
                CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_FACES>] : isConnectedTo1 not yet implement\n";
            }
            else
            {
                __face_id_in_elt_0 = elt_it->pos_first();

                if ( !isInitConnectionTo0 )
                {
                    form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                        this->expression(), face_ims[__face_id_in_elt_0], imRange,imTest,imTrial ) );
                    //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                    isInitConnectionTo0=true;
                }

                // update gmc
                detail::updateGmcWithRelationDifferentMeshType2<typename element_iterator::reference, typename FormType::space_1_type,im_range_type,
                                                                gmc_formTest_type,gmc_expr_type>(*elt_it,__form.testSpace(), gmcFormTest, gmcExpr0,
                                                                                                 idEltTest, mpl::int_<gmTestRangeRelation>() );

                detail::updateGmcWithRelationDifferentMeshType2<typename element_iterator::reference,typename FormType::space_2_type,im_range_type,
                                                                gmc_formTrial_type,gmc_expr_type>(*elt_it,__form.trialSpace(), gmcFormTrial, gmcExpr0,
                                                                                                  idEltTrial, mpl::int_<gmTrialRangeRelation>() );

                form->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr, face_ims[__face_id_in_elt_0] );
                form->integrate();
                form->assemble( /*elt_it->element( 0 ).id()*/ );
            }

        }

    } // for( auto lit = ... )

}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                           mpl::int_<MESH_FACES> /**/, mpl::true_ ) const
{
    LOG(INFO) << "assembleWithRelationDifferentMeshType Mortar case";
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc1_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef boost::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::gm1_1_type gm1_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type contextTest = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTest_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT> >::type::value;
    typedef typename gm_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc1_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef boost::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type contextTrial = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTrial_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT> >::type::value;
    typedef typename gm_formTrial_type::template Context<contextTrial,geoelement_formTrial_type> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<contextTrial,geoelement_formTrial_type> gmc1_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef boost::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //-----------------------------------------------------//

    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    static const uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_range_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtest_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtest_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtrial_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTrial_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtrial_type;

    im_range_type imRange;
    im1_range_type im1Range;
    im_formtest_type imTest;
    im_formtrial_type imTrial;

    //-----------------------------------------------------//

    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    typedef typename im_range_type::face_quadrature_type face_im_type;
    typedef typename im1_range_type::face_quadrature_type face_im1_type;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopc( imRange.nFaces() );
    std::vector<std::map<permutation_type, pc1_expr_ptrtype> > __geopc1( im1Range.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );
    std::vector<face_im1_type> face_ims1( im1Range.nFaces() );
    bool hasInitGeoPc=false;

    //-----------------------------------------------------//
    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type,true> form_mortar_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type, map_gmc1_formTrial_type,true> form1_mortar_context_type;

    typedef boost::shared_ptr<form_context_type> form_context_ptrtype;
    typedef boost::shared_ptr<form_mortar_context_type> form_mortar_context_ptrtype;
    typedef boost::shared_ptr<form1_context_type> form1_context_ptrtype;
    typedef boost::shared_ptr<form1_mortar_context_type> form1_mortar_context_ptrtype;


    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto elt_it = lit->template get<1>();
        auto const elt_en = lit->template get<2>();
        DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

        // check that we have elements to iterate over
        if ( (elt_it == elt_en) || (elt_it->isConnectedTo0() == false) )
            continue;

        //-----------------------------------------------//

        if (!hasInitGeoPc)
        {
            for ( uint16_type __f = 0; __f < imRange.nFaces(); ++__f )
            {
                face_ims[__f] = imRange.face( __f );
                face_ims1[__f] = im1Range.face( __f );

                for ( permutation_type __p( permutation_type::IDENTITY );
                      __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                    __geopc[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
                    __geopc1[__f][__p] = pc1_expr_ptrtype(  new pc1_expr_type( elt_it->element( 0 ).gm1(), im1Range.fpoints(__f, __p.value() ) ) );
                }
            }
            hasInitGeoPc=true;
        }

        //-----------------------------------------------//

        const bool rangeIsSubMeshFromTest = elt_it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
        const bool rangeIsSubMeshFromTrial = elt_it->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
        const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );
        const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );

        const size_type idEltRangeInit = elt_it->id();
        size_type idEltTestInit  = idEltRangeInit;
        if ( rangeIsSubMeshFromTest )
            idEltTestInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
        else if ( testIsSubMeshFromRange )
            idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
        size_type idEltTrialInit = idEltRangeInit;
        if ( rangeIsSubMeshFromTrial )
            idEltTrialInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
        else if ( trialIsSubMeshFromRange )
            idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );

        //-----------------------------------------------//
        // get the geometric mapping associated with element 0
        uint16_type __face_id_in_elt_0 = elt_it->pos_first();
        gmc_expr_ptrtype gmcExpr0( new gmc_expr_type( elt_it->element( 0 ).gm(), elt_it->element( 0 ), __geopc, __face_id_in_elt_0 ) );
        //gmc1_expr_ptrtype gmc1Expr0( new gmc1_expr_type( elt_it->element( 0 ).gm1(), elt_it->element( 0 ), __geopc1, __face_id_in_elt_0 ) );

        auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType2< typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                             idEltTestInit, gmcExpr0, mpl::int_<gmTestRangeRelation>() );

        auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType2< typename FormType::space_2_type,im_formtrial_type,
                                                                           gmc_formTrial_type,gmc_expr_type>( __form.trialSpace(), __form.trialSpace()->gm(), imTrial,
                                                                                                              idEltTrialInit, gmcExpr0, mpl::int_<gmTrialRangeRelation>() );

        //-----------------------------------------------//

        map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr0 ) );
        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
        form_context_ptrtype form;
        form1_context_ptrtype form1;
        form_mortar_context_ptrtype mortar_form;
        form1_mortar_context_ptrtype mortar_form1;


        bool isInitConnectionTo0=false;
        bool isInitConnectionTo1=false;

        // true if connected to another element, false otherwise
        if ( elt_it->isConnectedTo1() )
        {
            CHECK( false ) << "TODO!!!!\n";
#if 0
            uint16_type __face_id_in_elt_1 = it->pos_second();

            __c1 = gmc_ptrtype( new gmc_type( __gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );
            __c11 = gmc1_ptrtype( new gmc1_type( __gm1, it->element( 1 ), __geopc1, __face_id_in_elt_1 ) );
            map2_gmc_type mapgmc2( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                   fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );
            map21_gmc_type mapgmc21( fusion::make_pair<vf::detail::gmc<0> >( __c01 ),
                                     fusion::make_pair<vf::detail::gmc<1> >( __c11 ) );

            form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
            form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
            isInitConnectionTo1=true;
#endif
        }
        else
        {
            form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                this->expression(), face_ims[__face_id_in_elt_0], imRange,imTest,imTrial ) );
            mortar_form = form_mortar_context_ptrtype( new form_mortar_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                                     this->expression(), face_ims[__face_id_in_elt_0],
                                                                                     imRange,imTest,imTrial ) );
            //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
            isInitConnectionTo0=true;
        }



        for ( ; elt_it != elt_en; ++elt_it )
        {
            if ( elt_it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << elt_it->id() << " is a ghost face";
                continue;
            }
            // if is a interprocess faces, only integrate in one process
            if ( elt_it->isInterProcessDomain() && elt_it->partition1() > elt_it->partition2() )
                continue;

            const size_type idEltRange = elt_it->id();
            size_type idEltTest = idEltRange;
            if ( rangeIsSubMeshFromTest )
                idEltTest = elt_it->mesh()->subMeshToMesh( idEltRange );
            else if ( testIsSubMeshFromRange )
                idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
            size_type idEltTrial = idEltRange;
            if ( rangeIsSubMeshFromTrial )
                idEltTrial = elt_it->mesh()->subMeshToMesh( idEltRange );
            else if ( trialIsSubMeshFromRange )
                idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltRange );


            if ( elt_it->isConnectedTo1() )
            {
                CHECK ( false ) << "[assembleWithRelationDifferentMeshType<BiLinearForm,MESH_FACES>] : isConnectedTo1 not yet implement\n";
            }
            else
            {
                __face_id_in_elt_0 = elt_it->pos_first();

                if ( !isInitConnectionTo0 )
                {
                    if ( !elt_it->isOnBoundary() )
                        form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                            this->expression(), face_ims[__face_id_in_elt_0],
                                                                            imRange,imTest,imTrial ) );
                    else
                        mortar_form = form_mortar_context_ptrtype( new form_mortar_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                                                 this->expression(), face_ims[__face_id_in_elt_0],
                                                                                                 imRange,imTest,imTrial ) );
                    //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                    isInitConnectionTo0=true;
                }

                // update gmc
                detail::updateGmcWithRelationDifferentMeshType2<typename element_iterator::reference, typename FormType::space_1_type,im_range_type,
                                                                gmc_formTest_type,gmc_expr_type>(*elt_it,__form.testSpace(), gmcFormTest, gmcExpr0,
                                                                                                 idEltTest, mpl::int_<gmTestRangeRelation>() );

                detail::updateGmcWithRelationDifferentMeshType2<typename element_iterator::reference,typename FormType::space_2_type,im_range_type,
                                                                gmc_formTrial_type,gmc_expr_type>(*elt_it,__form.trialSpace(), gmcFormTrial, gmcExpr0,
                                                                                                  idEltTrial, mpl::int_<gmTrialRangeRelation>() );
                LOG(INFO) << "elt "  << elt_it->id() << " bdy = " << elt_it->isOnBoundary();
                if ( elt_it->isOnBoundary() )
                {
                    LOG(INFO) << "mortar element";
                    mortar_form->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr, face_ims[__face_id_in_elt_0] );
                    mortar_form->integrate();
                    mortar_form->assemble( /*elt_it->element( 0 ).id()*/ );
                }
                else
                {
                    LOG(INFO) << "internal element";
                    form->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr, face_ims[__face_id_in_elt_0] );
                    form->integrate();
                    form->assemble( /*elt_it->element( 0 ).id()*/ );
                }
            }

        }

    } // for( auto lit = ... )

}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
{
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc1_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef boost::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;

    // typedef on linearform :
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    // test space
    typedef typename FormType::gm_type gm_formTest_type;
    typedef typename FormType::gm1_type gm1_formTest_type;
    typedef typename FormType::mesh_test_element_type geoelement_formTest_type;
    static const size_type contextTest = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTest_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT> >::type::value;
    typedef typename gm_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<contextTest,geoelement_formTest_type> gmc1_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef boost::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    //-----------------------------------------------------//

    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_range_type;

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_formtest_type;
    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im2::template applyIMGeneral<geoelement_formTest_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im1_formtest_type;
    im_range_type imRange;
    im1_range_type im1Range;
    im_formtest_type imTest;

    //-----------------------------------------------------//

    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    typedef typename im_range_type::face_quadrature_type face_im_type;
    typedef typename im1_range_type::face_quadrature_type face_im1_type;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopc( imRange.nFaces() );
    std::vector<std::map<permutation_type, pc1_expr_ptrtype> > __geopc1( im1Range.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );
    std::vector<face_im1_type> face_ims1( im1Range.nFaces() );
    bool hasInitGeoPc=false;

    //-----------------------------------------------------//
    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type> form1_context_type;
    typedef boost::shared_ptr<form_context_type> form_context_ptrtype;
    typedef boost::shared_ptr<form1_context_type> form1_context_ptrtype;


    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto elt_it = lit->template get<1>();
        auto const elt_en = lit->template get<2>();
        DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

        // check that we have elements to iterate over
        if ( (elt_it == elt_en) || (elt_it->isConnectedTo0() == false) )
            continue;

        //-----------------------------------------------//

        if (!hasInitGeoPc)
        {
            for ( uint16_type __f = 0; __f < imRange.nFaces(); ++__f )
            {
                face_ims[__f] = imRange.face( __f );
                face_ims1[__f] = im1Range.face( __f );

                for ( permutation_type __p( permutation_type::IDENTITY );
                      __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                    __geopc[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
                    __geopc1[__f][__p] = pc1_expr_ptrtype(  new pc1_expr_type( elt_it->element( 0 ).gm1(), im1Range.fpoints(__f, __p.value() ) ) );
                }
            }
            hasInitGeoPc=true;
        }

        //-----------------------------------------------//

        const bool rangeIsSubMeshFromTest = elt_it->mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
        const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( elt_it->mesh() );

        const size_type idEltRangeInit = elt_it->id();
        size_type idEltTestInit  = idEltRangeInit;
        if ( rangeIsSubMeshFromTest )
            idEltTestInit = elt_it->mesh()->subMeshToMesh( idEltRangeInit );
        else if ( testIsSubMeshFromRange )
            idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );

        //-----------------------------------------------//
        // get the geometric mapping associated with element 0
        uint16_type __face_id_in_elt_0 = elt_it->pos_first();
        gmc_expr_ptrtype gmcExpr0( new gmc_expr_type( elt_it->element( 0 ).gm(), elt_it->element( 0 ), __geopc, __face_id_in_elt_0 ) );
        //gmc1_expr_ptrtype gmc1Expr0( new gmc1_expr_type( elt_it->element( 0 ).gm1(), elt_it->element( 0 ), __geopc1, __face_id_in_elt_0 ) );

        auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType2< typename FormType::space_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                             idEltTestInit, gmcExpr0, mpl::int_<gmTestRangeRelation>() );

        //-----------------------------------------------//

        map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr0 ) );
        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
        form_context_ptrtype form;
        form1_context_ptrtype form1;


        bool isInitConnectionTo0=false;
        bool isInitConnectionTo1=false;

        // true if connected to another element, false otherwise
        if ( elt_it->isConnectedTo1() )
        {
            CHECK( false ) << "TODO!!!!\n";
        }
        else
        {
            form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcExpr,
                                                                this->expression(), face_ims[__face_id_in_elt_0], imRange,imTest ) );
            //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
            isInitConnectionTo0=true;
        }


        for ( ; elt_it != elt_en; ++elt_it )
        {
            if ( elt_it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << elt_it->id() << " is a ghost face";
                continue;
            }
            // if is a interprocess faces, only integrate in one process
            if ( elt_it->isInterProcessDomain() && elt_it->partition1() > elt_it->partition2() )
                continue;

            const size_type idEltRange = elt_it->id();
            size_type idEltTest = idEltRange;
            if ( rangeIsSubMeshFromTest )
                idEltTest = elt_it->mesh()->subMeshToMesh( idEltRange );
            else if ( testIsSubMeshFromRange )
                idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );

            if ( elt_it->isConnectedTo1() )
            {
                CHECK ( false ) << "[assembleWithRelationDifferentMeshType<LinearForm,MESH_FACES>] : isConnectedTo1 not yet implement\n";
            }
            else
            {
                __face_id_in_elt_0 = elt_it->pos_first();

                if ( !isInitConnectionTo0 )
                {
                    form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcExpr,
                                                                        this->expression(), face_ims[__face_id_in_elt_0], imRange,imTest ) );
                    //form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                    isInitConnectionTo0=true;
                }

                // update gmc
                detail::updateGmcWithRelationDifferentMeshType2<typename element_iterator::reference, typename FormType::space_type,im_range_type,
                                                                gmc_formTest_type,gmc_expr_type>(*elt_it,__form.testSpace(), gmcFormTest, gmcExpr0,
                                                                                                 idEltTest, mpl::int_<gmTestRangeRelation>() );

                form->update( mapgmcFormTest,mapgmcFormTest,mapgmcExpr, face_ims[__face_id_in_elt_0] );
                form->integrate();
                form->assemble();
            }

        } // for ( ; elt_it ... )

    } // for( auto lit = ... )

}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                 mpl::int_<MESH_FACES> /**/,
                                                                 mpl::false_ ) const
{

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    /*typedef typename FormType::gm_type gm_form_type;
      typedef typename FormType::mesh_element_type geoelement_form_type;
      typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
      typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
      typedef typename gm_form_type::precompute_type pc_form_type;
      typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
      typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;*/

    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;



    // typedef on formcontext
    typedef typename im_range_type::face_quadrature_type face_im_type;

    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type,map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;


    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;


    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_range_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );

        }
    }


    uint16_type __face_id_in_elt_0 = elt_it->pos_first();

    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->element( 0 ).gm(),
                                                 elt_it->element( 0 ),
                                                 __geopcExpr,
                                                 __face_id_in_elt_0 ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(),  __form.testSpace()->fe()->points() ) );
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );

    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );

    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               imRange ) );

    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest,meshTrial );
    }

    //-----------------------------------------------//

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->result().begin();
        auto res_en = M_QPL->result().end();

        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& map = res_it->template get<1>();
            auto map_it = map.begin();
            auto const map_en = map.end();

            for ( ; map_it != map_en ; ++map_it )
            {
                auto const idEltTrial = map_it->first;
                auto const& eltTrial = meshTrial->element( idEltTrial );
                auto const& eltTest = meshTest->element( idEltTest );
                auto const& ptRefTest = map_it->second.template get<1>();
                auto const& ptRefTrial = map_it->second.template get<2>();
                auto const& themapQuad = map_it->second.template get<0>();

                auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
                auto gmcExpr_it = vec_gmcExpr.begin();
                auto const gmcExpr_en = vec_gmcExpr.end();
                bool isFirstExperience = true;

                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcFormTest->update( gmcExpr_it->template get<2>() );
                    geopcFormTrial->update( gmcExpr_it->template get<3>() );

                    gmcFormTest->update( eltTest,geopcFormTest );
                    gmcFormTrial->update( eltTrial,geopcFormTrial );
                    //std::cout << "eltTest.id() " << eltTest.id() << " eltTest.G() " << eltTest.G()
                    //          << " eltTrial.id() "<< eltTrial.id() << " eltTrial.G() " << eltTrial.G() << std::endl;
                    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();
                    formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );

                    formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                    isFirstExperience = false;
                }

                formc->assembleInCaseOfInterpolate();
            }

        }
    } // if (!M_QPL->hasPrecompute())
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcFormTest*/ ) );
                map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<3>()/*gmcFormTrial*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );
                __face_id_in_elt_0 = resQPLloc_it->template get<1>()->faceId();
                formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assembleInCaseOfInterpolate();
        }
    } //else

    delete formc;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                 mpl::int_<MESH_FACES> /**/,
                                                                 mpl::true_ ) const
{

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    /*typedef typename FormType::gm_type gm_form_type;
      typedef typename FormType::mesh_element_type geoelement_form_type;
      typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
      typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
      typedef typename gm_form_type::precompute_type pc_form_type;
      typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
      typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;*/

    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;



    // typedef on formcontext
    typedef typename im_range_type::face_quadrature_type face_im_type;

    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type,map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;


    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;


    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_range_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );

        }
    }


    uint16_type __face_id_in_elt_0 = elt_it->pos_first();

    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->element( 0 ).gm(),
                                                 elt_it->element( 0 ),
                                                 __geopcExpr,
                                                 __face_id_in_elt_0 ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(),  __form.testSpace()->fe()->points() ) );
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );

    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );

    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               imRange ) );

    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest,meshTrial );
    }

    //-----------------------------------------------//

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->result().begin();
        auto res_en = M_QPL->result().end();

        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& map = res_it->template get<1>();
            auto map_it = map.begin();
            auto const map_en = map.end();

            for ( ; map_it != map_en ; ++map_it )
            {
                auto const idEltTrial = map_it->first;
                auto const& eltTrial = meshTrial->element( idEltTrial );
                auto const& eltTest = meshTest->element( idEltTest );
                auto const& ptRefTest = map_it->second.template get<1>();
                auto const& ptRefTrial = map_it->second.template get<2>();
                auto const& themapQuad = map_it->second.template get<0>();

                auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
                auto gmcExpr_it = vec_gmcExpr.begin();
                auto const gmcExpr_en = vec_gmcExpr.end();
                bool isFirstExperience = true;

                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcFormTest->update( gmcExpr_it->template get<2>() );
                    geopcFormTrial->update( gmcExpr_it->template get<3>() );

                    gmcFormTest->update( eltTest,geopcFormTest );
                    gmcFormTrial->update( eltTrial,geopcFormTrial );
                    //std::cout << "eltTest.id() " << eltTest.id() << " eltTest.G() " << eltTest.G()
                    //          << " eltTrial.id() "<< eltTrial.id() << " eltTrial.G() " << eltTrial.G() << std::endl;
                    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();
                    formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );

                    formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                    isFirstExperience = false;
                }

                formc->assembleInCaseOfInterpolate();
            }

        }
    } // if (!M_QPL->hasPrecompute())
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcFormTest*/ ) );
                map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<3>()/*gmcFormTrial*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );
                __face_id_in_elt_0 = resQPLloc_it->template get<1>()->faceId();
                formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assembleInCaseOfInterpolate();
        }
    } //else

    delete formc;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
{

    typedef typename mpl::if_<mpl::bool_<eval::the_element_type::is_simplex>,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Simplex>::type >,
                              mpl::identity<typename Im::template applyIMGeneral<eval::the_element_type::nDim, expression_value_type, Hypercube>::type >
                              >::type::type im_range_type;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (test):
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;

    // typedef on formcontext
    typedef typename im_range_type::face_quadrature_type face_im_type;

    typedef typename FormType::template Context<map_gmc_form_type, expression_type, face_im_type,map_gmc_expr_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;


    element_iterator elt_it, elt_en;
    bool findEltForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findEltForInit; ++lit )
    {
        elt_it = lit->template get<1>();
        elt_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( elt_it == elt_en )
            continue;
        else
            findEltForInit = true;
    }
    if (!findEltForInit) return;

    //-----------------------------------------------//

    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange;
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
        }
    }


    uint16_type __face_id_in_elt_0 = elt_it->pos_first();

    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->element( 0 ).gm(),
                                                 elt_it->element( 0 ),
                                                 __geopcExpr,
                                                 __face_id_in_elt_0 ) );

    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcForm ) );
    map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( gmcForm ) );

    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcForm,
                                               mapgmcForm,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               imRange ) );

    //-----------------------------------------------//

    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts ) );
        M_QPL->update( meshTest );
    }

    //-----------------------------------------------//

    if (!M_QPL->hasPrecompute())
    {
        auto res_it = M_QPL->resultLinear().begin();
        auto const res_en = M_QPL->resultLinear().end();
        for ( ; res_it != res_en ; ++res_it )
        {
            auto const idEltTest = res_it->template get<0>();
            auto const& eltTest = meshTest->element( idEltTest );
            auto const& ptRefTest = res_it->template get<2>();
            auto const& themapQuad = res_it->template get<1>();

            auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest );
            auto gmcExpr_it = vec_gmcExpr.begin();
            auto const gmcExpr_en = vec_gmcExpr.end();
            bool isFirstExperience = true;
            for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
            {
                geopcForm->update( gmcExpr_it->template get<2>() );
                gmcForm->update( eltTest,geopcForm );
                map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( gmcForm ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();
                formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }

            formc->assemble();
        }
    }
    else
    {
        auto const& resQPL = M_QPL->getPrecompute(__form);
        auto resQPL_it = resQPL.begin();
        auto const resQPL_en = resQPL.end();
        for ( ; resQPL_it != resQPL_en ; ++resQPL_it)
        {
            auto resQPLloc_it = resQPL_it->begin();
            auto const resQPLloc_en = resQPL_it->end();
            bool isFirstExperience = true;
            for ( ; resQPLloc_it != resQPLloc_en ; ++resQPLloc_it)
            {
                map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<2>()/*gmcForm*/ ) );
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( resQPLloc_it->template get<1>() ) );
                __face_id_in_elt_0 = resQPLloc_it->template get<1>()->faceId();
                formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,face_ims[__face_id_in_elt_0],resQPLloc_it->template get<0>() );

                formc->integrateInCaseOfInterpolate( resQPLloc_it->template get<0>(),isFirstExperience );
                isFirstExperience = false;
            }
            formc->assemble();
        }
    } // else

    delete formc;
}




template<typename Elements, typename Im, typename Expr, typename Im2>
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluate( mpl::int_<MESH_ELEMENTS> ) const
{
    DLOG(INFO)  << "integrating over "
             << std::distance( this->beginElement(), this->endElement() )  << " elements\n";
    boost::timer __timer;

#if defined(FEELPP_HAS_TBB)
    if ( boption(_name="parallel.cpu.enable") && M_use_tbb )
    {
        //std::cout << "Integrator Uses TBB: " << M_use_tbb << "\n";
        element_iterator it = this->beginElement();
        element_iterator en = this->endElement();
        typedef ContextEvaluate<expression_type, im_type, typename eval::the_element_type> context_type;

        if ( it == en )
            return eval::zero();

        std::vector<boost::reference_wrapper<const typename eval::element_type> > _v;

        for ( auto _it = it; _it != en; ++_it )
            _v.push_back( boost::cref( *_it ) );

        tbb::blocked_range<decltype( _v.begin() )> r( _v.begin(), _v.end(), M_grainsize );
        context_type thecontext( this->expression(), this->im(), *it );

        if ( M_partitioner == "auto" )
            tbb::parallel_reduce( r,  thecontext );

        else if ( M_partitioner == "simple" )
            tbb::parallel_reduce( r,  thecontext, tbb::simple_partitioner() );

        //else if ( M_partitioner == "affinity" )
        //tbb::parallel_reduce( r,  thecontext, tbb::affinity_partitioner() );
        return thecontext.result();
    }
    else
#endif // defined(FEELPP_HAS_TBB)
#if defined(FEELPP_HAS_HARTS)
    if ( boption(_name="parallel.cpu.enable") && soption(_name="parallel.cpu.impl").find("harts.") == 0 )
    {
        RunTimeSystem::PerfCounterMng<std::string> perf_mng ;
        perf_mng.init("total") ;
        perf_mng.start("total") ;

        perf_mng.init("init") ;
        perf_mng.start("init") ;

        perf_mng.init("init0") ;
        perf_mng.start("init0") ;

        typename eval::matrix_type res( eval::matrix_type::Zero() );
        typedef Feel::vf::integrator::parallel::HartsContextEvaluate<expression_type, im_type, element_iterator, eval> harts_context_type;

        //std::cout << "Integrator Uses HARTS: " << M_use_harts << "\n";

        // typedef basic types
        typedef RunTimeSystem::PThreadDriver                        PTHDriverType ;
        typedef RunTimeSystem::TaskMng::Task<harts_context_type>    TaskType;
        typedef RunTimeSystem::TaskMng                              TaskMngType;
        typedef RunTimeSystem::TaskMng::ForkJoin<PTHDriverType>     PTHForkJoinTaskType;
        typedef RunTimeSystem::ThreadEnv                            ThreadEnv;

        typedef RunTimeSystem::DataHandler                          DataHandlerType;
        typedef RunTimeSystem::DataMng                              DataMngType;
        typedef RunTimeSystem::NumaAffinityMng                      NumaAffinityMngType;

        // Compute Number of MPI processes
        char * str = NULL;
        int nMPIProc = Environment::worldComm().size();

        // Compute a number of available cores for computation
        // Taking into account the number of MPI processes launched for this app on this node
        //NumaAffinityMngType numaAffMng(RunTimeSystem::NumaAffinityMng::eMode::Interleave);
        NumaAffinityMngType numaAffMng(RunTimeSystem::NumaAffinityMng::eMode::Block);

        // Compute Number of available CPU cores
        int paramNbCores = ioption(_name="parallel.cpu.restrict");
        int nTotalCoresNode = numaAffMng.get_num_cores();
        int nAvailCores = nTotalCoresNode - nMPIProc;
        int coresPerProcess = 0;
        int remainder = 0;

        // TOFIX: Have a repartition taking into account the NUMA architecture
        // Guess the number of threads that can be spawned
        if(paramNbCores <= 0)
        {
            coresPerProcess = nAvailCores / nMPIProc;
            remainder = nAvailCores % nMPIProc;
        }
        // Try to match the nb of cores in parameter
        else
        {
            int nWantedCores = nMPIProc * paramNbCores;
            if(nAvailCores - nWantedCores >= 0)
            {
                coresPerProcess = paramNbCores;
                remainder = nAvailCores - nMPIProc * paramNbCores;
            }
            /* if we can't match the asked number of cores */
            /* Then guess the maximal occupation repartition */
            else
            {
                coresPerProcess = nAvailCores / nMPIProc;
                remainder = nAvailCores % nMPIProc;
            }
        }

        /* create a task and data managers */
        TaskMngType taskMng;
        DataMngType dataMng;

        /* create a task list */
        std::vector<harts_context_type * > hce;
        std::vector<int> taskList;

        /* Using pthread */
        int nbThreads = coresPerProcess;

        //std::cout << "2. HARTS: nMPIProc=" << nMPIProc << ", nTotalCoresNode=" << nTotalCoresNode << ", coresPerProcess=" << coresPerProcess << ", remainder=" << remainder << ", nbElements="<< nbElts <<  std::endl;

        perf_mng.stop("init0") ;

        perf_mng.init("init1") ;
        perf_mng.start("init1") ;

        std::vector<std::vector<std::pair<element_iterator, element_iterator> > > _v;
        _v.resize(nbThreads);

        for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit)
        {
            /* set the iterator at the beginning of the elements */
            auto sit = lit->template get<1>();
            auto eit = lit->template get<2>();
            auto cit = sit;

            /* get number of elements */
            int nbElts = std::distance(sit, eit);
            int nbEltPerRange = nbElts / nbThreads;
            int remainderElt = nbElts % nbThreads;

            int j = 0;
            for(int i = 0; i < nbThreads; i++)
            {
                int nb = nbEltPerRange + (i < remainderElt ? 1 : 0);

                std::cout << Environment::worldComm().rank() << "|" << i << " nbElts=" << nb << std::endl;

                /* save the current iterator position */
                cit = sit;
                /* advance the iterator and save the new position */
                std::advance(sit, nb);
                //std::cout << "T" << i << " adv:" << nb << " " << nbEltPerRange <<  std::endl;
                //std::cout << "T" << i << " adv:" << nb << " " << nbEltPerRange << "(" << j << ", " << (j + nb) << ")" << std::endl;
                j = j + nb;
                _v.at(i).push_back(std::make_pair(cit, sit));
            }
        }

        perf_mng.stop("init1") ;

        perf_mng.init("init2") ;
        perf_mng.start("init2") ;

        ThreadEnv * threadEnv = new ThreadEnv(nbThreads);

        //threadEnv->SetAffinity(numaAffMng);

        PTHDriverType forkjoin(threadEnv);
        PTHForkJoinTaskType* compFkTask = new PTHForkJoinTaskType(forkjoin, taskMng.getTasks());
        taskList.push_back(taskMng.addNew(compFkTask));

        std::ostringstream oss1;
        perf_mng.init("init2.2") ;
        perf_mng.start("init2.2") ;

        perf_mng.init("init2.2.0") ;
        perf_mng.init("init2.2.1") ;
        perf_mng.init("init2.2.2") ;
        perf_mng.init("init2.2.3") ;

        /* getting parameters to avoid problems with const */
        expression_type expr = this->expression();
        im_type im = this->im();
        element_iterator elt_it = M_elts.begin()->template get<1>();

        /* create tasks and associate them data */
        for(int i = 0 ; i< nbThreads; i++)
        {
            DataHandlerType * d_tid = dataMng.getNewData();
            d_tid->set<int>(&i);
            DataHandlerType * d_elts = dataMng.getNewData();
            d_elts->set<std::vector<std::pair<element_iterator, element_iterator> > >(&_v[i]);
            DataHandlerType * d_expr = dataMng.getNewData();
            d_expr->set<expression_type>(&expr);
            DataHandlerType * d_im = dataMng.getNewData();
            d_im->set<im_type>(&im);
            DataHandlerType * d_elt = dataMng.getNewData();
            d_elt->set<element_iterator>(&elt_it);

            if(i == 0) { perf_mng.start("init2.2.0") ; }
            //harts_context_type * t = new harts_context_type(i, this->expression(), this->im(), *(M_elts.begin()->template get<1>()) /* *it */);
            harts_context_type * t = new harts_context_type();
            if(i == 0) { perf_mng.stop("init2.2.0") ; }
            if(i == 0) { perf_mng.start("init2.2.1") ; }
            hce.push_back(t);
            if(i == 0) { perf_mng.stop("init2.2.1") ; }
            if(i == 0) { perf_mng.start("init2.2.2") ; }
            TaskType* task = new TaskType(t);
            if(i == 0) { perf_mng.stop("init2.2.2") ; }
            if(i == 0) { perf_mng.start("init2.2.3") ; }
            task->args().add("threadId", DataHandlerType::R, d_tid);
            task->args().add("elements", DataHandlerType::R, d_elts);
            task->args().add("expr", DataHandlerType::R, d_expr);
            task->args().add("im", DataHandlerType::R, d_im);
            task->args().add("elt", DataHandlerType::R, d_elt);
            if(i == 0) { perf_mng.stop("init2.2.3") ; }

            typename TaskType::FuncType f = &harts_context_type::computeCPU;
            task->set("cpu",f);
            int uid = taskMng.addNew(task);
            compFkTask->add(uid);
        }
        perf_mng.stop("init2.2") ;

        perf_mng.stop("init2") ;

        perf_mng.stop("init") ;

        perf_mng.init("comp") ;
        perf_mng.start("comp") ;

        //Environment::writeCPUData(""); 

        RunTimeSystem::StdScheduler scheduler;
        taskMng.run(scheduler, taskList);

        //Environment::writeCPUData(""); 

        perf_mng.stop("comp") ;

        for(int i = 0 ; i< nbThreads; i++)
        {
            res += hce[i]->result();
        }

        for(int i = 0; i < nbThreads; i++)
        {
            std::cout << Environment::worldComm().rank() << "|" << i << " elapsed=" << hce[i]->elapsed() << std::endl;
            hce[i]->printPerfInfo();
        }

        // Free memory
        #if 1
        delete threadEnv;
        //delete compFkTask;

        for(int i = 0 ; i < nbThreads; i++)
        {
            delete hce[i];
        }
        #endif

        perf_mng.stop("total") ;
        
#if 0
        std::ostringstream oss;
        oss << fs::current_path().string() << "/" << Environment::worldComm().size() << "_" << Environment::worldComm().rank() << "-" << nbThreads << ".dat";
        std::ofstream f;

        f.open(oss.str().c_str(), std::ofstream::out | std::ofstream::trunc);

        for(int i = 0; i < nbThreads; i++)
        {
            f << hce[i]->elapsed() << " ";
            std::cout << Environment::worldComm().rank() << "|" << i << " " << hce[i]->elapsed() << std::endl;
        }

        f << perf_mng.getValueInSeconds("total") << std::endl;

        f.close();
#endif

        std::cout << Environment::worldComm().rank() <<  " Evaluation output: " << res << " in " 
                  << perf_mng.getValueInSeconds("total") << " (" << perf_mng.getValueInSeconds("init") << " ("
                  << perf_mng.getValueInSeconds("init0") << ", " << perf_mng.getValueInSeconds("init1") << ", " << perf_mng.getValueInSeconds("init2") << ") "
                  << ", " << perf_mng.getValueInSeconds("comp") << ")" << std::endl;

        std::cout << Environment::worldComm().rank() <<  " Evaluation " 
                  << perf_mng.getValueInSeconds("init2.2.0") << " "
                  << perf_mng.getValueInSeconds("init2.2.1") << " "
                  << perf_mng.getValueInSeconds("init2.2.2") << " "
                  << perf_mng.getValueInSeconds("init2.2.3") << std::endl;

        DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";
        return res;
    }
    else
#endif // defined(FEELPP_HAS_HARTS)
#if defined(FEELPP_HAS_OPENMP)
    if ( boption(_name="parallel.cpu.enable") && soption(_name="parallel.cpu.impl") == "openmp" )
    {
        // alternative: measure time with omp_get_wtime()
#if defined(FEELPP_HAS_HARTS)
        RunTimeSystem::PerfCounterMng<std::string> perf_mng ;
        perf_mng.init("total") ;
        perf_mng.start("total") ;

        perf_mng.init("init") ;
        perf_mng.start("init") ;

        perf_mng.init("init0") ;
        perf_mng.start("init0") ;
#endif

        typename eval::matrix_type res( eval::matrix_type::Zero() );
        typedef Feel::vf::integrator::parallel::HartsContextEvaluate<expression_type, im_type, element_iterator, eval> harts_context_type;

        //std::cout << "Integrator Uses HARTS: " << M_use_harts << "\n";

        // Compute Number of MPI processes
        char * str = NULL;
        int nMPIProc = Environment::worldComm().size();

        // Compute a number of available cores for computation
        // Taking into account the number of MPI processes launched for this app on this node
        //std::cout << "omp_get_num_procs()=" << omp_get_num_procs() << std::endl;

        // Compute Number of available CPU cores
        int paramNbCores = ioption(_name="parallel.cpu.restrict");
        int nTotalCoresNode = omp_get_num_procs();
        int nAvailCores = nTotalCoresNode - nMPIProc;
        int coresPerProcess = 0;
        int remainder = 0;

        // TOFIX: Have a repartition taking into account the NUMA architecture
        // Guess the number of threads that can be spawned
        if(paramNbCores <= 0)
        {
            coresPerProcess = nAvailCores / nMPIProc;
            remainder = nAvailCores % nMPIProc;
        }
        // Try to match the nb of cores in parameter
        else
        {
            int nWantedCores = nMPIProc * paramNbCores;
            if(nAvailCores - nWantedCores >= 0)
            {
                coresPerProcess = paramNbCores;
                remainder = nAvailCores - nMPIProc * paramNbCores;
            }
            /* if we can't match the asked number of cores */
            /* Then guess the maximal occupation repartition */
            else
            {
                coresPerProcess = nAvailCores / nMPIProc;
                remainder = nAvailCores % nMPIProc;
            }
        }

        /* create a task list */
        std::vector<harts_context_type * > hce;

        /* Using pthread */
        int nbThreads = coresPerProcess;

        std::cout << "2. HARTS: nMPIProc=" << nMPIProc << ", nTotalCoresNode=" << nTotalCoresNode << ", coresPerProcess=" << coresPerProcess << ", remainder=" << remainder << std::endl;

#if defined(FEELPP_HAS_HARTS)
        perf_mng.stop("init0") ;
        perf_mng.init("init1") ;
        perf_mng.start("init1") ;
#endif

        std::vector<std::vector<std::pair<element_iterator, element_iterator> > > _v;
        _v.resize(nbThreads);

        for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit)
        {
            /* set the iterator at the beginning of the elements */
            auto sit = lit->template get<1>();
            auto eit = lit->template get<2>();
            auto cit = sit;

            /* get number of elements */
            int nbElts = std::distance(sit, eit);
            int nbEltPerRange = nbElts / nbThreads;
            int remainderElt = nbElts % nbThreads;

            int j = 0;
            for(int i = 0; i < nbThreads; i++)
            {
                int nb = nbEltPerRange + (i < remainderElt ? 1 : 0);

                std::cout << Environment::worldComm().rank() << "|" << i << " nbElts=" << nb << std::endl;

                /* save the current iterator position */
                cit = sit;
                /* advance the iterator and save the new position */
                std::advance(sit, nb);
                //std::cout << "T" << i << " adv:" << nb << " " << nbEltPerRange <<  std::endl;
                std::cout << "T" << i << " adv:" << nb << " " << nbEltPerRange << "(" << j << ", " << (j + nb) << ")" << std::endl;
                j = j + nb;
                _v.at(i).push_back(std::make_pair(cit, sit));
            }
        }

#if defined(FEELPP_HAS_HARTS)
        perf_mng.stop("init1") ;

        perf_mng.stop("init") ;

        perf_mng.init("comp") ;
        perf_mng.start("comp") ;
#endif

        /* Save the previous count of openmp threads */
        /* in case the use set OMP_NUM_THREADS */
        /* and we don't want to influence the next omp loop */
        int prevThreadCount = omp_get_num_threads();
        omp_set_num_threads(nbThreads);

        value_type * out = new value_type[nbThreads];

        expression_type expr = this->expression();
        im_type im = this->im();
        element_iterator elt_it = M_elts.begin()->template get<1>();


#pragma omp parallel
        {
            int id = omp_get_thread_num();

#if defined(FEELPP_HAS_HARTS)
            RunTimeSystem::PerfCounterMng<std::string> local_perf_mng ;
#endif

            //local_perf_mng.init("Creation") ;
            //local_perf_mng.start("Creation") ;
            //harts_context_type * t = new harts_context_type(id, this->expression(), this->im(), *(M_elts.begin()->template get<1>()) /* *it */);
            harts_context_type * t = new harts_context_type();
            //local_perf_mng.stop("Creation") ;

#if defined(FEELPP_HAS_HARTS)
            std::cout << Environment::worldComm().rank() << "|" << id << " Creation:" << local_perf_mng.getValueInSeconds("Creation") << std::endl;
#endif

#if defined(FEELPP_HAS_HARTS)
        struct timespec ts1;
        struct timespec ts2;
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
#endif

            t->computeCPUOMP(id, &expr, &im, &elt_it, &(_v[id]));
            out[id] = t->result();

#if defined(FEELPP_HAS_HARTS)
        clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
        double t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
        double t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
        std::cout << Environment::worldComm().globalRank() << "|" << id << " elapsed (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;
#endif

            std::cout << Environment::worldComm().rank() << "|" << id << " elapsed=" << t->elapsed() << std::endl;
            t->printPerfInfo();

            delete t;
        }

        /* restore previous openmp thread count */
        omp_set_num_threads(prevThreadCount);

#if defined(FEELPP_HAS_HARTS)
        perf_mng.stop("comp") ;
#endif

        for(int i = 0 ; i< nbThreads; i++)
        {
            res += out[i];
        }
        delete[] out;

#if defined(FEELPP_HAS_HARTS)
        perf_mng.stop("total") ;
        
        std::cout << Environment::worldComm().rank() <<  " Evaluation output: " << res << " in " 
                  << perf_mng.getValueInSeconds("total") << " (" << perf_mng.getValueInSeconds("init") << " ("
                  << perf_mng.getValueInSeconds("init0") << ", " << perf_mng.getValueInSeconds("init1") << ", " << perf_mng.getValueInSeconds("init2") << ") "
                  << ", " << perf_mng.getValueInSeconds("comp") << ")" << std::endl;
#endif

        DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";
        return res;
    }
    else
#endif
    if ( 1 )
    {
#if defined(FEELPP_HAS_HARTS)
        RunTimeSystem::PerfCounterMng<std::string> perf_mng ;
        perf_mng.init("total") ;
        perf_mng.start("total") ;
#endif

        //
        // some typedefs
        //
        typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
        typedef typename boost::remove_const<const_t>::type the_element_type;
        typedef the_element_type element_type;
        typedef typename the_element_type::gm_type gm_type;
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename eval::gmc_type gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

        typedef typename the_element_type::gm1_type gm1_type;
        typedef boost::shared_ptr<gm1_type> gm1_ptrtype;
        typedef typename eval::gmc1_type gmc1_type;
        typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;

        //typedef typename eval_expr_type::value_type value_type;
        //typedef typename Im::value_type value_type;
        typename eval::matrix_type res( eval::matrix_type::Zero() );

        for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
        {
            auto it = lit->template get<1>();
            auto en = lit->template get<2>();

            // make sure that we have elements to iterate over (return 0
            // otherwise)
            if ( it == en )
                continue;

            //std::cout << "nbElts: " << std::distance(it, en) << std::endl;

            //std::cout << "0" << std::endl;

            //
            // Precompute some data in the reference element for
            // geometric mapping and reference finite element
            //
            // warning this is not efficient here, we want to use the geometric mapping
            // from the elements in order to take advantage of the cache if possible
            // this change hsa been made in order to circumvent a bug which is not yet found
            //#warning INEFFICIENT CODE HERE : TO DEBUG
            //gm_ptrtype gm( new gm_type) ;//it->gm();
            gm_ptrtype gm( it->gm() );
            //std::cout << "0.5" << std::endl;
            gm1_ptrtype gm1( new gm1_type ); //it->gm1();
            //std::cout << "0.6:  " << gm1.use_count() << " " << gm.use_count() << std::endl;
            //DDLOG(INFO) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";
            typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( gm,
                                                                               this->im().points() ) );
            //std::cout << "1" << std::endl;
            typename eval::gmpc1_ptrtype __geopc1( new typename eval::gmpc1_type( gm1,
                                                                                  this->im().points() ) );

            //std::cout << "2" << std::endl;
            //it = this->beginElement();

            // wait for all the guys
#ifdef FEELPP_HAS_MPI
            auto const& worldComm = const_cast<MeshBase*>( it->mesh() )->worldComm();
#if 0
            if ( worldComm.localSize() > 1 )
            {
                worldComm.localComm().barrier();
            }
#endif
#endif

            // possibly high order
            gmc_ptrtype __c( new gmc_type( gm, *it, __geopc ) );
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            //std::cout << "3" << std::endl;
            typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
            eval_expr_type expr( expression(), mapgmc );
            typedef typename eval_expr_type::shape shape;
            //std::cout << "4" << std::endl;

            // order 1
            gmc1_ptrtype __c1( new gmc1_type( gm1, *it, __geopc1 ) );
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            //std::cout << "5" << std::endl;
            typedef typename expression_type::template tensor<map_gmc1_type> eval_expr1_type;
            eval_expr1_type expr1( expression(), mapgmc1 );

            //std::cout << "6" << std::endl;


            //value_type res1 = 0;
            for ( ; it != en; ++it )
            {
                switch ( M_gt )
                {
                default:
                case  GeomapStrategyType::GEOMAP_HO :
                {
#if 1
                    __c->update( *it );
                    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                    expr.update( mapgmc );
                    const gmc_type& gmc = *__c;

                    M_im.update( gmc );
#endif


                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        {
                            res( c1,c2 ) += M_im( expr, c1, c2 );
                        }
                }
                break;

                case GeomapStrategyType::GEOMAP_O1:
                {
#if 1
                    //DDLOG(INFO) << "geomap o1" << "\n";
                    __c1->update( *it );
                    map_gmc1_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                    expr1.update( mapgmc );
                    const gmc1_type& gmc = *__c1;

                    M_im.update( gmc );
#endif


                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        {
                            res( c1,c2 ) += M_im( expr1, c1, c2 );
                        }
                }
                break;

                case GeomapStrategyType::GEOMAP_OPT:
                {
                    //DDLOG(INFO) << "geomap opt" << "\n";
                    if ( it->isOnBoundary() )
                    {
#if 1
                        //DDLOG(INFO) << "boundary element using ho" << "\n";
                        __c->update( *it );
                        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                        expr.update( mapgmc );
                        const gmc_type& gmc = *__c;

                        M_im.update( gmc );
#endif


                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                            {
                                res( c1,c2 ) += M_im( expr, c1, c2 );
                            }
                    }

                    else
                    {
#if 1
                        //DDLOG(INFO) << "interior element using order 1" << "\n";
                        __c1->update( *it );
                        map_gmc1_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                        expr1.update( mapgmc );
                        const gmc1_type& gmc = *__c1;

                        M_im.update( gmc );
#endif

                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                            {
                                res( c1,c2 ) += M_im( expr1, c1, c2 );
                            }
                    }
                }

                //break;
                }
            }
        }
#if defined(FEELPP_HAS_HARTS)
        perf_mng.stop("total") ;
        std::cout << Environment::worldComm().rank() <<  " Total: " << perf_mng.getValueInSeconds("total") << std::endl;
#endif

        DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";
        return res;
    }
}
template<typename Elements, typename Im, typename Expr, typename Im2>
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluate( mpl::int_<MESH_FACES> ) const
{
    DLOG(INFO) << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << "faces\n";
    boost::timer __timer;

    //
    // some typedefs
    //
    typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::entity_type the_element_type;
    typedef typename the_element_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, the_element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    //typedef typename eval_expr_type::value_type value_type;
    //typedef typename Im::value_type value_type;

    //BOOST_MPL_ASSERT_MSG( the_element_type::nDim > 1, INVALID_DIM, (mpl::int_<the_element_type::nDim>, mpl::int_<the_face_element_type::nDim>, mpl::identity<the_face_element_type>, mpl::identity<the_element_type> ) );;

    QuadMapped<im_type> qm;
    typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename eval::gmpc_type pc_type;
    typedef typename eval::gmpc_ptrtype pc_ptrtype;
    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( im().nFaces() );
    typedef typename im_type::face_quadrature_type face_im_type;

    std::vector<im_face_type> __integrators;

    typename eval::matrix_type res( eval::matrix_type::Zero() );
    typename eval::matrix_type res0( eval::matrix_type::Zero() );
    typename eval::matrix_type res1( eval::matrix_type::Zero() );


    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        // make sure that we have elements to iterate over (return 0
        // otherwise)
        if ( (it == en) || (it->isConnectedTo0()==false) )
            continue;
        //return typename eval::matrix_type( eval::matrix_type::Zero() );

        CHECK( it->isConnectedTo0() ) << "invalid face with id=" << it->id();
        CHECK( it->element(0).gm() ) << "invalid geometric transformation assocated to face id="
                                     <<  it->id() << " and element id " << it->element(0).id();
        gm_ptrtype gm = it->element( 0 ).gm();

        //DDLOG(INFO) << "[integrator] evaluate(faces), gm is cached: " << gm->isCached() << "\n";
        for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            __integrators.push_back( im( __f ) );

            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f][__p]->size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( gm, ppts[__f].find( __p )->second ) );
            }
        }

        uint16_type __face_id_in_elt_0 = it->pos_first();

        // get the geometric mapping associated with element 0
        gmc_ptrtype __c0( new gmc_type( gm, it->element( 0 ), __geopc, __face_id_in_elt_0 ) );

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef boost::shared_ptr<eval_expr_type> eval_expr_ptrtype;
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
        eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );
        expr->init( im() );

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
        typedef typename expression_type::template tensor<map2_gmc_type> eval2_expr_type;
        typedef boost::shared_ptr<eval2_expr_type> eval2_expr_ptrtype;
        eval2_expr_ptrtype expr2;

        // true if connected to another element, false otherwise
        bool isConnectedTo1 = it->isConnectedTo1();

        // get the geometric mapping associated with element 1
        gmc_ptrtype __c1;

        //value_type res = 0;
        //value_type res1 = 0;
        if ( isConnectedTo1 )
        {
            uint16_type __face_id_in_elt_1 = it->pos_second();

            __c1 = gmc_ptrtype( new gmc_type( gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );

            map2_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                  fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );

            expr2 = eval2_expr_ptrtype( new eval2_expr_type( expression(), mapgmc ) );
            expr2->init( im() );
        }

        //
        // start the real intensive job:
        // -# iterate over all elements to integrate over
        // -# construct the associated geometric mapping with the reference element
        // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
        // -# assemble the local contribution in the global representation of the bilinear form
        //
        for ( ; it != en; ++it )
        {
            if ( it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << it->id() << " is a ghost face";
                continue;
            }
            // if is a interprocess faces, only integrate in one process
            if ( it->isInterProcessDomain() && it->partition1() > it->partition2() )
                continue;

            if ( it->isConnectedTo1() )
            {
                FEELPP_ASSERT( it->isOnBoundary() == false   )
                    ( it->id() ).error( "face on boundary but connected on both sides" );
                uint16_type __face_id_in_elt_0 = it->pos_first();
                uint16_type __face_id_in_elt_1 = it->pos_second();

                __c0->update( it->element( 0 ), __face_id_in_elt_0 );
                __c1->update( it->element( 1 ), __face_id_in_elt_1 );

#if 0
                std::cout << "face " << it->id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << "  elt 0 : " << it->element( 0 ).id() << "\n"
                          << "  elt 0 G: " << it->element( 0 ).G() << "\n"
                          << "  node elt 0 0 :" << it->element( 0 ).point( it->element( 0 ).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                          << "  node elt 0 1 :" << it->element( 0 ).point( it->element( 0 ).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                          << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                          << "  real nodes 0: " << __c0->xReal() << "\n";
                std::cout << "face " << it->id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << " elt 1 : " << it->element( 1 ).id() << "\n"
                          << "  elt 1 G: " << it->element( 1 ).G() << "\n"
                          << "  node elt 1 0 :" << it->element( 1 ).point( it->element( 1 ).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                          << "  node elt 1 1 :" << it->element( 1 ).point( it->element( 1 ).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
                          << " ref nodes 1 :" << __c1->xRefs() << "\n"
                          << " real nodes 1:" << __c1->xReal() << "\n";
#endif

                __typeof__( im( __face_id_in_elt_0 ) ) im_face ( im( __face_id_in_elt_0 ) );
                //std::cout << "pts = " << im_face.points() << "\n";
                map2_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                      fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );

                expr2->update( mapgmc, __face_id_in_elt_0 );
                const gmc_type& gmc = *__c0;

                __integrators[__face_id_in_elt_0].update( gmc );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    {
                        res( c1,c2 ) += __integrators[__face_id_in_elt_0]( *expr2, c1, c2 );
                    }
            }

            else
            {
                //LOG_IF( !it->isConnectedTo0(), WARN ) << "integration invalid boundary face";
                if ( !it->isConnectedTo0() || it->pos_first() == invalid_uint16_type_value )
                    continue;
                uint16_type __face_id_in_elt_0 = it->pos_first();
                __c0->update( it->element( 0 ), __face_id_in_elt_0 );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                expr->update( mapgmc, __face_id_in_elt_0 );
                //expr->update( mapgmc );
                const gmc_type& gmc = *__c0;

                __integrators[__face_id_in_elt_0].update( gmc );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    {
                        res( c1,c2 ) += __integrators[__face_id_in_elt_0]( *expr, c1, c2 );
                    }
            } // !isConnectedTo1
        } // for loop on face
    }
    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    DLOG(INFO) << "integrating over faces done in " << __timer.elapsed() << "s\n";
    return res;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluate( mpl::int_<MESH_POINTS> ) const
{
    DLOG(INFO)  << "integrating over "
             << std::distance( this->beginElement(), this->endElement() )  << " points\n";

    // first loop on the points, then retrieve the elements to which they belong
    // and evaluate the integrand expression and accumulate it


}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename P0hType>
typename P0hType::element_type
Integrator<Elements, Im, Expr, Im2>::broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const
{
    DLOG(INFO) << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << " elements\n";
    boost::timer __timer;

    //
    // some typedefs
    //
    typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_element_type;
    typedef the_element_type element_type;
    typedef typename the_element_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::POINT, the_element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    //typedef typename eval::gm_type gmc_type;
    //typedef boost::shared_ptr<gmc_type> gmc_ptrtype;

    //typedef typename eval_expr_type::value_type value_type;
    //typedef typename Im::value_type value_type;

    auto p0 = P0h->element( "p0" );
    // set to 0 first
    p0.zero();

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();


        // make sure that we have elements to iterate over (return 0
        // otherwise)
        if ( it == en )
            continue;
            //return p0;

        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        gm_ptrtype gm = it->gm();
        //DDLOG(INFO) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";
        typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( gm,
                                                                           this->im().points() ) );


        it = this->beginElement();
        // wait for all the guys
#ifdef FEELPP_HAS_MPI
        auto const& worldComm = const_cast<MeshBase*>( it->mesh() )->worldComm();

#if 0
        if ( worldComm.size() > 1 )
        {
            worldComm.barrier();
        }
#endif
#endif



        gmc_ptrtype __c( new gmc_type( gm, *it, __geopc ) );
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        eval_expr_type expr( expression(), mapgmc );
        typedef typename eval_expr_type::shape shape;

        //value_type res1 = 0;
        for ( ; it != en; ++it )
        {
            boost::timer ti;
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            expr.update( mapgmc );
            const gmc_type& gmc = *__c;

            M_im.update( gmc );


            for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
            {
                size_type i= P0h->dof()->localToGlobal( it->id(), 0, c1 ).index();
                 double v = M_im( expr, c1, 0 );
                p0.set( i, v );
            }
        }
    }
    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    DLOG(INFO) << "integrating over elements done in " << __timer.elapsed() << "s\n";

    return p0;
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename P0hType>
typename P0hType::element_type
Integrator<Elements, Im, Expr, Im2>::broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const
{
    DLOG(INFO) << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << "faces\n";
    boost::timer __timer;

    //
    // some typedefs
    //
    typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
    typedef typename boost::remove_const<const_t>::type the_face_element_type;
    typedef typename the_face_element_type::super2::entity_type the_element_type;
    typedef typename the_element_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, the_element_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    //typedef typename eval_expr_type::value_type value_type;
    //typedef typename Im::value_type value_type;

    //BOOST_MPL_ASSERT_MSG( the_element_type::nDim > 1, INVALID_DIM, (mpl::int_<the_element_type::nDim>, mpl::int_<the_face_element_type::nDim>, mpl::identity<the_face_element_type>, mpl::identity<the_element_type> ) );;

    QuadMapped<im_type> qm;
    typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename eval::gmpc_type pc_type;
    typedef typename eval::gmpc_ptrtype pc_ptrtype;
    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( im().nFaces() );
    typedef typename im_type::face_quadrature_type face_im_type;

    std::vector<im_face_type> __integrators;

    auto p0 = P0h->element( "p0" );
    // set to 0 first
    p0.zero();

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        // make sure that we have elements to iterate over (return 0
        // otherwise)
        if ( it == en )
            continue;
        //return p0;

        gm_ptrtype gm = it->element( 0 ).gm();

        //DDLOG(INFO) << "[integrator] evaluate(faces), gm is cached: " << gm->isCached() << "\n";
        for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            __integrators.push_back( im( __f ) );

            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f][__p]->size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( gm, ppts[__f].find( __p )->second ) );
            }
        }


        uint16_type __face_id_in_elt_0 = it->pos_first();

        // get the geometric mapping associated with element 0
        gmc_ptrtype __c0( new gmc_type( gm, it->element( 0 ), __geopc, __face_id_in_elt_0 ) );

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef boost::shared_ptr<eval_expr_type> eval_expr_ptrtype;
        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
        eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );
        expr->init( im() );

        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
        typedef typename expression_type::template tensor<map2_gmc_type> eval2_expr_type;
        typedef boost::shared_ptr<eval2_expr_type> eval2_expr_ptrtype;
        eval2_expr_ptrtype expr2;

        // true if connected to another element, false otherwise
        bool isConnectedTo1 = it->isConnectedTo1();

        // get the geometric mapping associated with element 1
        gmc_ptrtype __c1;

        //value_type res = 0;
        //value_type res1 = 0;
        if ( isConnectedTo1 )
        {
            uint16_type __face_id_in_elt_1 = it->pos_second();

            __c1 = gmc_ptrtype( new gmc_type( gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );

            map2_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                  fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );

            expr2 = eval2_expr_ptrtype( new eval2_expr_type( expression(), mapgmc ) );
            expr2->init( im() );
        }

        //
        // start the real intensive job:
        // -# iterate over all elements to integrate over
        // -# construct the associated geometric mapping with the reference element
        // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
        // -# assemble the local contribution in the global representation of the bilinear form
        //
        for ( ; it != en; ++it )
        {

            if ( it->isConnectedTo1() )
            {
                FEELPP_ASSERT( it->isOnBoundary() == false   )
                    ( it->id() ).error( "face on boundary but connected on both sides" );
                uint16_type __face_id_in_elt_0 = it->pos_first();
                uint16_type __face_id_in_elt_1 = it->pos_second();

                __c0->update( it->element( 0 ), __face_id_in_elt_0 );
                __c1->update( it->element( 1 ), __face_id_in_elt_1 );

#if 0
                std::cout << "face " << it->id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << "  elt 0 : " << it->element( 0 ).id() << "\n"
                          << "  elt 0 G: " << it->element( 0 ).G() << "\n"
                          << "  node elt 0 0 :" << it->element( 0 ).point( it->element( 0 ).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                          << "  node elt 0 1 :" << it->element( 0 ).point( it->element( 0 ).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                          << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                          << "  real nodes 0: " << __c0->xReal() << "\n";
                std::cout << "face " << it->id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << " elt 1 : " << it->element( 1 ).id() << "\n"
                          << "  elt 1 G: " << it->element( 1 ).G() << "\n"
                          << "  node elt 1 0 :" << it->element( 1 ).point( it->element( 1 ).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                          << "  node elt 1 1 :" << it->element( 1 ).point( it->element( 1 ).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
                          << " ref nodes 1 :" << __c1->xRefs() << "\n"
                          << " real nodes 1:" << __c1->xReal() << "\n";
#endif

                __typeof__( im( __face_id_in_elt_0 ) ) im_face ( im( __face_id_in_elt_0 ) );
                //std::cout << "pts = " << im_face.points() << "\n";
                map2_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ),
                                      fusion::make_pair<vf::detail::gmc<1> >( __c1 ) );

                expr2->update( mapgmc, __face_id_in_elt_0 );
                const gmc_type& gmc = *__c0;

                __integrators[__face_id_in_elt_0].update( gmc );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                {
                    size_type i0 = P0h->dof()->localToGlobal( it->element( 0 ), 0, c1 ).index();
                    size_type i1 =  P0h->dof()->localToGlobal( it->element( 1 ), 0, c1 ).index();
                    double v = __integrators[__face_id_in_elt_0]( *expr2, c1, 0 );
                    p0.add( i0, v );
                    p0.add( i1, v );
                }
            }

            else
            {
                uint16_type __face_id_in_elt_0 = it->pos_first();
                __c0->update( it->element( 0 ), __face_id_in_elt_0 );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                expr->update( mapgmc, __face_id_in_elt_0 );
                //expr->update( mapgmc );
                const gmc_type& gmc = *__c0;

                __integrators[__face_id_in_elt_0].update( gmc );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                {
                    size_type i0 = P0h->dof()->localToGlobal( it->element( 0 ), 0, c1 ).index();
                    double v = __integrators[__face_id_in_elt_0]( *expr, c1, 0 );
                    p0.add( i0, v );
                }
            } // !isConnectedTo1
        } // for loop on face
    }
    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    DLOG(INFO) << "integrating over faces done in " << __timer.elapsed() << "s\n";
    return p0;
}
/// \endcond

#if 0
/**
 * integrate an expression \c expr over a set of convexes \c elts
 * using the integration rule \c im .
 */
template<typename Elts, typename Im, typename ExprT>
Expr<Integrator<Elts, Im, ExprT, Im> >
integrate( Elts const& elts,
           Im const& im,
           ExprT const& expr,
           GeomapStrategyType gt = GeomapStrategyType::GEOMAP_HO )
{
    typedef Integrator<Elts, Im, ExprT, Im> expr_t;
    return Expr<expr_t>( expr_t( elts, im, expr, gt, im ) );
}
#endif

//Macro which get the good integration order
# define VF_VALUE_OF_IM(O)                                              \
    boost::mpl::if_< boost::mpl::bool_< O::imIsPoly > ,                 \
                     typename boost::mpl::if_< boost::mpl::greater< boost::mpl::int_<O::imorder>, boost::mpl::int_<19> > , \
                                               boost::mpl::int_<19>,    \
                                               boost::mpl::int_<O::imorder> >::type >::type , \
                                                                               boost::mpl::int_<10> >::type::value \
/**/

/**
 * integrate an expression \c expr over a set of convexes \c elts
 * using the integration rule \c im .
 */
template<typename Elts, typename Im, typename ExprT, typename Im2 = Im>
Expr<Integrator<typename Feel::detail::quadptlocrangetype<Elts>::type, Im, ExprT, Im2> >
integrate_impl( Elts const& elts,
                Im const& im,
                ExprT const& expr,
                GeomapStrategyType const& gt,
                Im2 const& im2,
                bool use_tbb,
                bool use_harts,
                int grainsize,
                std::string const& partitioner,
                boost::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<Elts>::type, Im, ExprT > > quadptloc
                = boost::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<Elts>::type, Im, ExprT > >() )
{

    typedef typename Feel::detail::quadptlocrangetype< Elts >::type range_type;
    typedef Integrator<range_type, Im, ExprT, Im2> expr_t;

    typedef typename boost::unwrap_reference< typename boost::tuples::template element<1,range_type>::type >::type element_iterator;
    static const uint16_type geoOrder = element_iterator::value_type::nOrder;
    LOG_IF(WARNING, gt != GeomapStrategyType::GEOMAP_HO && geoOrder == 1 ) << "you use a non standard geomap : ";
    return Expr<expr_t>( expr_t( elts, im, expr, gt, im2, use_tbb, use_harts, grainsize, partitioner, quadptloc ) );
}


/// \cond DETAIL
namespace detail
{
template<typename Args>
struct integrate_type
{
    typedef typename clean2_type<Args,tag::expr,Expr<Cst<double> > >::type _expr_type;
    typedef typename Feel::detail::quadptlocrangetype<typename clean_type<Args,tag::range>::type>::type _range_type;
    typedef typename boost::unwrap_reference< typename boost::tuples::template element<1, _range_type>::type >::type _element_iterator;
    static const uint16_type geoOrder = _element_iterator::value_type::nOrder;

    //typedef _Q< ExpressionOrder<_range_type,_expr_type>::value > the_quad_type;
    typedef typename clean2_type<Args,tag::quad, _Q< ExpressionOrder<_range_type,_expr_type>::value > >::type _quad_type;
    typedef typename clean2_type<Args,tag::quad1, _Q< ExpressionOrder<_range_type,_expr_type>::value_1 > >::type _quad1_type;
    typedef Expr<Integrator<_range_type, _quad_type, _expr_type, _quad1_type> > expr_type;

    typedef boost::shared_ptr<QuadPtLocalization<_range_type,_quad_type,_expr_type > > _quadptloc_ptrtype;
};
} // detail

/// \endcond



} // vf


} // feel



#endif /* FEELPP_INTEGRATORS_HPP */
