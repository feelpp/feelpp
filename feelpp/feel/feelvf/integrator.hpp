/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#ifndef FEELPP_VF_INTEGRATORS_HPP
#define FEELPP_VF_INTEGRATORS_HPP 1

#include <cxxabi.h>
#include <typeinfo>

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
#include <feel/feelvf/matvec.hpp>
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

    static const Feel::size_type context = Expr::context|vm::JACOBIAN;
    static const bool is_terminal = false;

    //static const uint16_type imorder = 0;
    //static const bool imIsPoly = true;

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
    template<typename Func>
    static const bool has_test_basis = Expr::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = Expr::template has_trial_basis<Func>;
    using test_basis = typename Expr::test_basis;
    using trial_basis = typename Expr::trial_basis;

    static const Feel::size_type iDim = boost::tuples::template element<0, Elements>::type::value;
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
        using range_entity_type = std::remove_cv_t<std::remove_reference_t<typename boost::unwrap_reference<typename element_iterator::value_type>::type>>;
        using the_element_type = typename range_entity_type::super2::template Element<range_entity_type>::type;

        using im_type = im_t<the_element_type,expression_value_type>;
        using im2_type = im_t<the_element_type,expression_value_type>;

        typedef the_element_type element_type;
        typedef typename the_element_type::gm_type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;
        typedef typename the_element_type::gm1_type gm1_type;
        typedef std::shared_ptr<gm1_type> gm1_ptrtype;
        static const size_type gmc_context_v = expression_type::context|vm::JACOBIAN;
        typedef typename gm_type::template Context< the_element_type> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm1_type::template Context<the_element_type> gmc1_type;
        typedef std::shared_ptr<gmc1_type> gmc1_ptrtype;
#if 0
        typedef typename gm_type::template precompute<im_type::numPoints>::type gmpc_type;
        typedef typename gm_type::template precompute<im_type::numPoints>::ptrtype gmpc_ptrtype;
#else
        typedef typename gm_type::PreCompute gmpc_type;
        typedef std::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef typename gm1_type::PreCompute gmpc1_type;
        typedef std::shared_ptr<gmpc1_type> gmpc1_ptrtype;
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

    using shape = typename eval::shape;
    typedef typename eval::im_type im_type;
    typedef typename eval::im2_type im2_type;
    typedef typename im_type::face_quadrature_type im_face_type;
    typedef typename im2_type::face_quadrature_type im2_face_type;
    //typedef typename eval::value_type value_type;
    typedef typename eval::matrix_type matrix_type;
    typedef typename eval::matrix_type value_type;
    typedef typename eval::matrix_type evaluate_type;
    using index_type = typename eval::element_type::index_type;
    using size_type = typename eval::element_type::size_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Integrator( Elements const& elts,
                Im const& __im,
                expression_type const& __expr,
                GeomapStrategyType gt,
                Im2 const& __im2,
                bool use_tbb, bool use_harts, int grainsize, std::string const& partitioner,
                std::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > qpl )
        :
        M_elts(),
        M_eltbegin( elts.begin() ),
        M_eltend( elts.end() ),
        M_wc( elts.worldCommPtr() ),
        M_im( __im ),
        M_im2( __im2 ),
        M_expr( __expr ),
        M_gt( gt ),
        M_use_tbb( use_tbb ),
        M_use_harts( use_harts ),
        M_grainsize( grainsize ),
        M_partitioner( partitioner ),
        M_QPL( qpl )
    {
        M_elts.push_back( elts );
        DLOG(INFO) << "im quad order : " << M_im.order();
        DLOG(INFO) << "im quad1 order : " << M_im2.order();
        DLOG(INFO) << "im : " << M_im.points() << " w:" << M_im.weights();
        DLOG(INFO) << "Integrator constructor from expression\n";
    }

    Integrator( std::list<Elements> const& elts,
                Im const& __im,
                expression_type const& __expr,
                GeomapStrategyType gt,
                Im2 const& __im2,
                bool use_tbb, bool use_harts, int grainsize, std::string const& partitioner,
                std::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > qpl )
        :
        M_elts( elts ),
        M_wc( Environment::worldCommPtr() ),
        M_im( __im ),
        M_im2( __im2 ),
        M_expr( __expr ),
        M_gt( gt ),
        M_use_tbb( use_tbb ),
        M_use_harts( use_harts ),
        M_grainsize( grainsize ),
        M_partitioner( partitioner ),
        M_QPL( qpl )
    {
        DLOG(INFO) << "Integrator constructor from expression\n";
        if ( !elts.empty() )
        {
            M_eltbegin = elts.begin()->template get<1>();
            M_eltend = elts.begin()->template get<2>();
            M_wc = elts.begin()->worldCommPtr();
        }
    }

    Integrator( Integrator const& __vfi )
        :
        M_elts( __vfi.M_elts) ,
        M_eltbegin( __vfi.M_eltbegin ),
        M_eltend( __vfi.M_eltend ),
        M_wc( __vfi.M_wc ),
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
#if 0
        typedef _Q< ExpressionOrder<Elements,expr_type>::value > quad_type;
        typedef _Q< ExpressionOrder<Elements,expr_type>::value_1 > quad1_type;
#else
        using expr_order_t = ExpressionOrder<Elements,expr_type>;
        using quad_type = im_t<typename expr_order_t::the_element_type,typename expr_type::value_type>;
        using quad1_type = im_t<typename expr_order_t::the_element_type,typename expr_type::value_type>;
#endif
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
            using expr_order_t = ExpressionOrder<Elements,e_type>;
            using quad_type = im_t<typename expr_order_t::the_element_type,typename e_type::value_type>;
            using quad1_type = im_t<typename expr_order_t::the_element_type,typename e_type::value_type>;
            //typedef _Q< ExpressionOrder<Elements,e_type>::value > quad_type;
            //typedef _Q< ExpressionOrder<Elements,e_type>::value_1 > quad1_type;
            typedef std::shared_ptr<QuadPtLocalization<Elements,quad_type,expr_type > > quadptloc_ptrtype;
            quad_type quad( expr_order_t::value(new_expr) );
            quad1_type quad1( expr_order_t::value_1(new_expr) );

#if 0
            auto the_ims = vf::detail::integrate_im_type<Elements,expr_type>::im( quad,quad1,expr );
            auto const& the_im = the_ims.first;
            auto const& the_im1 = the_ims.second;
#endif


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

    worldcomm_t& worldComm() { return *M_wc; }
    worldcomm_t const& worldComm() const { return *M_wc; }
    worldcomm_ptr_t& worldCommPtr() { return M_wc; }
    worldcomm_ptr_t const& worldCommPtr() const { return M_wc; }

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
    void assemble( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __form ) const;

    template<typename Elem1, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& __v,
                   FormType& __form ) const;

    //typename expression_type::template tensor<Geo_t>::value_type
    template<typename P0hType>
    typename P0hType::element_type  broken( std::shared_ptr<P0hType>& P0h ) const
    {
        auto p0 = broken( P0h, mpl::int_<iDim>() );
        return p0;
    }
    /**
     * evaluate the integral for each entry of the vector \c v
     */
    template<typename T,int M, int N=1>
    decltype(auto)
    evaluate( std::vector<Eigen::Matrix<T,M,N>> const& v,
              bool parallel=true ) const;

    matrix_type
    evaluate( bool parallel=true ) const
     {
         typename eval::matrix_type loc = this->evaluateImpl();

        if ( !parallel )
            return loc;

        else // parallel
        {
            typename eval::matrix_type glo( loc );
            // maybe better to create anoter worldcomm which split the mesh worldcomm
            // with only partition that contains at least one element (Vincent C.)
            // and thus argument worldComm can be remove
            // auto const& worldcomm = const_cast<MeshBase<>*>( this->beginElement()->mesh() )->worldComm();

            if ( M_wc->localSize() > 1 )
            {
                mpi::all_reduce( M_wc->localComm(),
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
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef std::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
        //typedef vf::detail::FormContextBase<map_gmc_type,im_type> fcb_type;
        typedef form_context_type fcb_type;
        typedef std::shared_ptr<fcb_type> fcb_ptrtype;


        Context( FormType& _form,
                 ExprType const& _expr,
                 IMType const& _im,
                 EltType const& _elt )
            :
            M_geopc( new gmpc_type( _form.gm(), _im.points() ) ),
            M_c( new gmc_type( _form.gm(), _elt, M_geopc, invalid_uint16_type_value, _expr.dynamicContext()  ) ),
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
        typedef std::shared_ptr<gm_type> gm_ptrtype;
        typedef typename eval::gmc_type gmc_type;
        typedef typename eval::gmpc_type gmpc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef std::shared_ptr<gmpc_type> gmpc_ptrtype;
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
            M_c( new gmc_type( M_gm, _elt, M_geopc, invalid_uint16_type_value, _expr.dynamicContext() ) ),
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
            M_c( new gmc_type( M_gm, c.M_c->element(), M_geopc, invalid_uint16_type_value, c.M_expr.dynamicContext() ) ),
            M_expr( c.M_expr ),
            M_im( c.M_im ),
            M_ret( eval::matrix_type::Zero() )
        {}

        ContextEvaluate( ContextEvaluate const& c )
            :
            M_gm( new gm_type( *c.M_gm ) ),
            //M_geopc( new gmpc_type( M_gm, c.M_im.points() ) ),
            M_geopc( c.M_geopc ),
            M_c( new gmc_type( M_gm, c.M_c->element(), M_geopc; invalid_uint16_type_value, c.M_expr.dynamicContext() ) ),
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
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool hasRelation/*, mpl::false_*/  ) const;

    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<false> /**/, bool hasRelation ) const;

    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<true> /**/, bool hasRelation ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<false> /**/, bool hasRelation ) const;

    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
        {
            DLOG(INFO) << "Intergrator::assembleWithRelationDifferentMeshType mortar case:" << (FE1::is_mortar || FE2::is_mortar);
            assembleWithRelationDifferentMeshType( __form, mpl::int_<MESH_ELEMENTS>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::true_ ) const;

    template<typename FE,typename VectorType,typename ElemContType>
    void assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleWithRelationDifferentMeshType( vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;

    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
        {
            DLOG(INFO) << "Intergrator::assembleInCaseOfInterpolate :" << (FE1::is_mortar || FE2::is_mortar);
            assembleInCaseOfInterpolate( __form, mpl::int_<MESH_ELEMENTS>(), mpl::bool_<FE1::is_mortar||FE2::is_mortar>() );
        }

    // bilinear interpolation case (no relation)
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
        {
            DLOG(INFO) << "Intergrator::assembleInCaseOfInterpolate mortar case:" << (FE1::is_mortar || FE2::is_mortar );
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
    typename P0hType::element_type  broken( std::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const;
    template<typename P0hType>
    typename P0hType::element_type  broken( std::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const;

    template <int iDimDummy=iDim,std::enable_if_t< iDimDummy == MESH_ELEMENTS , bool> = true>
    typename eval::matrix_type evaluateImpl() const;
    template <int iDimDummy=iDim,std::enable_if_t< iDimDummy == MESH_FACES /*|| ( iDimDummy == MESH_EDGES && eval::gm_type::nDim == 2)*/ , bool> = true>
    typename eval::matrix_type evaluateImpl() const;
    template <int iDimDummy=iDim,std::enable_if_t< iDimDummy == MESH_POINTS , bool> = true>
    typename eval::matrix_type evaluateImpl() const;

private:

    std::list<Elements> M_elts;
    element_iterator M_eltbegin;
    element_iterator M_eltend;
    worldcomm_ptr_t M_wc;
    mutable im_type M_im;
    mutable im2_type M_im2;
    expression_type   M_expr;
    GeomapStrategyType M_gt;
    bool M_use_tbb;
    bool M_use_harts;
    int M_grainsize;
    std::string M_partitioner;
    mutable std::shared_ptr<QuadPtLocalization<Elements, Im, Expr > > M_QPL;
    //     mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, double, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, double> > > M_profile_local_assembly;

    //     mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, double, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, double> > > M_profile_global_assembly;
};

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename Elem1, typename Elem2, typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( std::shared_ptr<Elem1> const& __u,
        std::shared_ptr<Elem2> const& __v,
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

    // specifiy matrix (form2) is in assembly state
    __form.matrixPtr()->setIsClosed( false );

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
    const bool test_related_to_trial = __v->mesh()->isRelatedTo( __u->mesh() ) &&
        ( __u->mesh()->isRelatedTo( boost::unwrap_ref(*it).mesh() ) || __v->mesh()->isRelatedTo( boost::unwrap_ref(*it).mesh() ) );
    DLOG(INFO) << "[integrator::assemble bilinear form] with_relation_mesh (same_mesh: " << same_mesh_type::value
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
Integrator<Elements, Im, Expr, Im2>::assemble( std::shared_ptr<Elem1> const& __v,
        FormType& __form ) const
{
#if 0
    details::GlobalVectorAssembler<iDim, self_type, Elem1, FormType>( *this,
            __v,
            __form );
#endif

    typedef typename boost::is_same<typename eval::gmc_type::element_type,typename Elem1::mesh_type::element_type>::type same_mesh_type;

    // specifiy vector (form1) is in assembly state
    __form.vectorPtr()->setIsClosed( false );

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


    const bool test_related_to_range = __v->mesh()->isRelatedTo( boost::unwrap_ref(*it).mesh() );

    //if ( dynamic_cast<void*>( const_cast<MeshBase<>*>( it->mesh() ) ) == dynamic_cast<void*>( __v->mesh().get() ) )
    if ( test_related_to_range )
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<same_mesh_type::value>(), test_related_to_range/*true*/ );

    else
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<false>(), test_related_to_range/*false*/ );

    //assemble( __form, mpl::int_<iDim>(), mpl::bool_<true>() );
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/, bool /*hasRelation*/ ) const
{
    tic();
    DLOG(INFO) << "[integrator::assemble FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true>\n";

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
            typedef typename eval::gmc_type gmc_type;
            typedef std::shared_ptr<gmc_type> gmc_ptrtype;
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
            typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
            typedef form_context_type fcb_type;
            using focb_ptrtype = std::shared_ptr<fcb_type>;

            typedef typename eval::gm1_type gm1_type;
            typedef typename eval::gmc1_type gmc1_type;
            typedef std::shared_ptr<gmc1_type> gmc1_ptrtype;
            typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
            typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type> form1_context_type;
            typedef form1_context_type fcb1_type;
            using focb1_ptrtype = std::shared_ptr<fcb1_type>;

            // mortar context
            static const bool has_mortar_test = FormType::test_space_type::is_mortar;
            static const bool has_mortar_trial = FormType::trial_space_type::is_mortar;
            static const int mortarTag = (has_mortar_test && has_mortar_trial)? 3 : ( (has_mortar_test)? 1 : ( (has_mortar_trial)? 2 : 0 ) );

            typedef typename FormType::template Context<map_gmc_type, expression_type, im_type,map_gmc_type,map_gmc_type,mortarTag> form_mortar_context_type;
            typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type,map_gmc1_type,map_gmc1_type,mortarTag> form1_mortar_context_type;
            using mortar_focb_ptrtype = std::shared_ptr<form_mortar_context_type>;
            using mortar_focb1_ptrtype = std::shared_ptr<form1_mortar_context_type>;

            bool useGeomapHO = (M_gt == GeomapStrategyType::GEOMAP_HO) || (M_gt == GeomapStrategyType::GEOMAP_OPT );
            bool useGeomapO1 = (M_gt == GeomapStrategyType::GEOMAP_O1) || (M_gt == GeomapStrategyType::GEOMAP_OPT );

            gmc_ptrtype __c;
            focb_ptrtype formc;
            gmc1_ptrtype __c1;
            focb1_ptrtype formc1;
            mortar_focb_ptrtype formcm;
            mortar_focb1_ptrtype formc1m;
            size_type nElt_ho = 0, nElt_o1 = 0;

            for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
            {
                element_iterator it = lit->template get<1>();
                element_iterator en = lit->template get<2>();
                DLOG(INFO) << "integrating over " << std::distance( it, en )  << " elements\n";

                // check that we have elements to iterate over
                if ( it == en )
                    continue;

                auto const& eltInit = boost::unwrap_ref( *it );

                bool rangeAndTestUseSameMesh = eltInit.mesh()->isSameMesh( __form.testSpace()->mesh() );

                size_type idEltTestInit = eltInit.id();
                if ( !rangeAndTestUseSameMesh )
                {
                    if ( eltInit.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                        idEltTestInit = eltInit.mesh()->subMeshToMesh( idEltTestInit );
                    else if ( __form.testSpace()->mesh()->isSubMeshFrom( eltInit.mesh() ) )
                        idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltTestInit );
                }

                auto const& eltTestInit = rangeAndTestUseSameMesh? eltInit :  __form.testSpace()->mesh()->element( idEltTestInit );

                if ( useGeomapHO && !__c )
                {
                    // Precompute some data in the reference element for geometric mapping and reference finite element
                    auto __geopc = __form.gm()->preCompute( this->im().points() );
                    __c = __form.gm()->template context<eval::gmc_context_v>( eltTestInit, __geopc, this->expression().dynamicContext() );
                    auto mapgmc = vf::mapgmc(__c);
                    formc = std::make_shared<form_context_type>( __form, mapgmc, mapgmc, mapgmc,
                                                                 this->expression(), this->im() );
                    if constexpr ( mortarTag > 0 )
                         formcm = std::make_shared<form_mortar_context_type>( __form, mapgmc, mapgmc, mapgmc,
                                                                              this->expression(), this->im() );
                }

                if ( useGeomapO1 && !__c1 )
                {
                    auto __geopc1 = __form.gm1()->preCompute( this->im2().points() );
                    __c1 = __form.gm1()->template context<eval::gmc_context_v>( eltTestInit, __geopc1, this->expression().dynamicContext() );
                    auto mapgmc1 = vf::mapgmc(__c1);
                    formc1 = std::make_shared<form1_context_type>( __form, mapgmc1, mapgmc1, mapgmc1,
                                                                   this->expression(), this->im2() );
                    if constexpr ( mortarTag > 0 )
                          formc1m = std::make_shared<form1_mortar_context_type>( __form, mapgmc1, mapgmc1, mapgmc1,
                                                                                 this->expression(), this->im2() );
                }

                //int nelt = std::distance( this->beginElement(), this->endElement() );

                //
                // start the real intensive job:
                // -# iterate over all elements to integrate over
                // -# construct the associated geometric mapping with the reference element
                // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
                // -# assemble the local contribution in the global representation of the bilinear form
                //
                const bool updateCtxAndIntegrate = !(__form.testSpace()->mesh()->isCartesian()  &&
                                                     ( !hasPOINT(__c->dynamicContext() ) && !hasINTERPOLANT(__c->dynamicContext()) ) );

                for ( ; it != en; ++it )
                {
                    auto const& eltCur = boost::unwrap_ref( *it );

                    size_type idElt = eltCur.id();
                    if ( !rangeAndTestUseSameMesh )
                    {
                        if ( eltCur.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                            idElt = eltCur.mesh()->subMeshToMesh( idElt );
                        else if ( __form.testSpace()->mesh()->isSubMeshFrom( eltCur.mesh() ) )
                            idElt = __form.testSpace()->mesh()->meshToSubMesh( idElt );
                    }

                    if ( useGeomapHO )
                    {
                        if ( formc->isZero( idElt ) )
                            continue;
                    }
                    else
                    {
                        if ( formc1->isZero( idElt ) )
                            continue;
                    }

                    auto const& eltTest = rangeAndTestUseSameMesh? eltCur : __form.testSpace()->mesh()->element( idElt );

                    // 0 : ho , 1 -> o1, 2 -> mortar ho, 3 -> mortar o1
                    int currentAssemblyType = 0;
                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                        currentAssemblyType = 0;
                        break;
                    case GeomapStrategyType::GEOMAP_O1:
                        currentAssemblyType = 1;
                        break;
                    case GeomapStrategyType::GEOMAP_OPT:
                        if ( eltCur.isOnBoundary() )
                            currentAssemblyType = 0;
                        else
                            currentAssemblyType = 1;
                        break;
                    }

                    bool useMortarAssembly = false;
                    if constexpr ( mortarTag > 0 )
                    {
                        if ( currentAssemblyType == 0 )
                            useMortarAssembly = ( has_mortar_test && eltTest.isOnBoundary() ) || ( has_mortar_trial && formc->trialElementIsOnBoundary( idElt ) );
                        else
                            useMortarAssembly = ( has_mortar_test && eltTest.isOnBoundary() ) || ( has_mortar_trial && formc1->trialElementIsOnBoundary( idElt ) );
                        if ( useMortarAssembly )
                            currentAssemblyType += 2;
                    }


                    switch ( currentAssemblyType )
                    {
                    default:
                    case 0:
                        //Feel::cout << "update ctx and integrate HO" << std::endl;
                        if ( updateCtxAndIntegrate || ( nElt_ho == 0 ) )
                        {
                            __c->template update<eval::gmc_context_v>( eltTest );
                            auto mapgmc = vf::mapgmc( __c );
                            formc->update( mapgmc, mapgmc, mapgmc );
                            formc->integrate();
                        }
                        else
                        {
                            __c->setElement( eltTest );
                        }
                        formc->assemble();
                        ++nElt_ho;
                        break;
                    case 1:
                        if ( updateCtxAndIntegrate || ( nElt_o1 == 0 )  )
                        {
                            __c1->template update<eval::gmc_context_v>( eltTest );
                            auto mapgmc1 = vf::mapgmc( __c1 );
                            formc1->update( mapgmc1, mapgmc1, mapgmc1 );
                            formc1->integrate();
                        }
                        else
                        {
                            __c1->setElement( eltTest );
                        }
                        formc1->assemble();
                        ++nElt_o1;
                        break;
                    case 2: // mortar ho
                        if constexpr ( mortarTag > 0 )
                        {
                            __c->template update<eval::gmc_context_v>( eltTest );
                            auto mapgmc = vf::mapgmc( __c );
                            formcm->update( mapgmc, mapgmc, mapgmc );
                            formcm->integrate();
                            formcm->assemble();
                        }
                        break;
                    case 3: //mortar o1
                        if constexpr ( mortarTag > 0 )
                        {
                            __c1->template update<eval::gmc_context_v>( eltTest );
                            auto mapgmc = vf::mapgmc( __c1 );
                            formc1m->update( mapgmc, mapgmc, mapgmc );
                            formc1m->integrate();
                            formc1m->assemble();
                        }
                        break;
                    }


                } // end loop on elements
            } // end loop on list of elements

            toc("Integrator::assemble form MESH_ELEMENTS", FLAGS_v>1);
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
template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( std::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt ,std::shared_ptr<GmcExprType> const& gmcExpr,mpl::int_<0> /**/ )
{
    if constexpr ( std::is_same_v<GmcType,GmcExprType> )
        return gmcExpr;
    else
    {
        throw std::logic_error("GmcExprType must be of the same type as GmcType");
        return nullptr;
    }
}
template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType( std::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                       size_type idElt,std::shared_ptr<GmcExprType> const& _expr,mpl::int_<1> /**/ )
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
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;

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
    gmc_ptrtype gmc = gm->template context<gmc_v>( eltInit.element0(), __geopc, eltInit.pos_first() /*face_id_in_elt_0*/, _expr->dynamicContext() );

    return gmc;
}
template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType( std::shared_ptr<SpaceType> const& /*space*/,
                                        std::shared_ptr<GmcType> /*gmc*/, std::shared_ptr<GmcExprType> /*gmcExpr*/,
                                        size_type /*idElt*/, mpl::int_<0> /**/ )
{
    // nothing to do!
}
template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
void
updateGmcWithRelationDifferentMeshType( std::shared_ptr<SpaceType> const& space,
                                        std::shared_ptr<GmcType> gmc, std::shared_ptr<GmcExprType> const& gmcExpr,
                                        size_type idElt, mpl::int_<1> /**/ )
{
#if 0
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;

    auto const& theface = space->mesh()->face( idElt );
    bool findPermutation=false;
    for ( permutation_type __p( permutation_type::IDENTITY );
          __p < permutation_type( permutation_type::N_PERMUTATIONS ) && !findPermutation; ++__p )
    {
        // update only xReal in gmc
        gmc->template update<vm::POINT>( theface.element0(), theface.pos_first(), __p, false );

        bool check=true;
        for ( uint16_type i=0;i<gmc->nPoints() && check;++i )
            for (uint16_type d=0;d<GmcType::NDim;++d)
                check = check && ( std::abs(gmc->xReal(i)[d]-gmcExpr->xReal(i)[d])<1e-8 );

        // if check compute full gmc context with the good permutation
        if (check) { gmc->template update<gmc_v>( theface.element0(), theface.pos_first(), __p ); findPermutation=true; }
    }
    CHECK(findPermutation) << "the permutation of quad point was not found\n";
#else
    auto const& theface = space->mesh()->face( idElt );
    auto [found,perm] = gmc->template updateFromMatchingNodes<gmc_v>( theface.element0(), theface.pos_first(), gmcExpr );
    CHECK(found) << "the permutation of quad point was not found";
#endif
}


template<typename SpaceTestType, typename SpaceTrialType, typename RangeType>
struct ManageRelationDifferentMeshType
{
    using index_type = typename SpaceTestType::index_type;

    using element_mesh_test_type = typename MeshTraits<typename SpaceTestType::mesh_type>::element_type;
    using element_mesh_trial_type = typename MeshTraits<typename SpaceTrialType::mesh_type>::element_type;
    using element_wrapper_test_type = typename MeshTraits<typename SpaceTestType::mesh_type>::elements_reference_wrapper_type::value_type;
    using element_wrapper_trial_type = typename MeshTraits<typename SpaceTrialType::mesh_type>::elements_reference_wrapper_type::value_type;
    static const int range_idim_type = boost::tuples::template element<0,range_t<RangeType> >::type::value;

    ManageRelationDifferentMeshType( std::shared_ptr<SpaceTestType> const& spaceTest, std::shared_ptr<SpaceTrialType> const& spaceTrial, RangeType const& range )
        :
        M_spaceTest( spaceTest ),
        M_spaceTrial( spaceTrial ),
        M_hasMeshSupportPartialTest( spaceTest->dof()->hasMeshSupport() && spaceTest->dof()->meshSupport()->isPartialSupport() ),
        M_hasMeshSupportPartialTrial( spaceTrial->dof()->hasMeshSupport() && spaceTrial->dof()->meshSupport()->isPartialSupport() ),
        M_rangeIsSubMeshFromTest( false ), M_rangeIsSubMeshFromTrial( false ), M_testIsSubMeshFromRange( false ), M_trialIsSubMeshFromRange( false ), M_rangeAndTestUseSameMesh( false ), M_rangeAndTrialUseSameMesh( false ),
        M_rangeEntityInit( nullptr ), M_eltIndexInitTest( invalid_v<index_type> ),  M_eltIndexInitTrial( invalid_v<index_type> )
        {
            if ( M_hasMeshSupportPartialTest )
                M_meshSupportTest = spaceTest->dof()->meshSupport();
            if ( M_hasMeshSupportPartialTrial )
                M_meshSupportTrial = spaceTrial->dof()->meshSupport();

            bool doCheckRelation = false;
            for ( auto const& entityWrap : range )
            {
                auto const& entity = unwrap_ref( entityWrap );
                if ( !doCheckRelation )
                {
                    M_rangeIsSubMeshFromTest = entity.mesh()->isSubMeshFrom( M_spaceTest->mesh() );
                    M_rangeIsSubMeshFromTrial = entity.mesh()->isSubMeshFrom( M_spaceTrial->mesh() );
                    M_testIsSubMeshFromRange = M_spaceTest->mesh()->isSubMeshFrom( entity.mesh() );
                    M_trialIsSubMeshFromRange = M_spaceTrial->mesh()->isSubMeshFrom( entity.mesh() );
                    M_rangeAndTestUseSameMesh = entity.mesh()->isSameMesh( M_spaceTest->mesh() );
                    M_rangeAndTrialUseSameMesh = entity.mesh()->isSameMesh( M_spaceTrial->mesh() );

                    doCheckRelation = true;
                }

                auto eltsInit = this->eltsRelatedToRange( entity );
                if ( eltsInit.empty() )
                    continue;
                M_faceConnectionIdInit = std::get<0>( eltsInit.front() ); // only for range of faces
                M_eltIndexInitTest = unwrap_ref( std::get<1>( eltsInit.front() ) ).id();
                M_eltIndexInitTrial = unwrap_ref( std::get<2>( eltsInit.front() ) ).id();
                M_rangeEntityInit = &entity;

                break;
            }

        }

    std::vector<std::tuple<uint16_type, element_wrapper_test_type, element_wrapper_trial_type>>
    eltsRelatedToRange( entity_range_t<RangeType> const& entity )
        {
            std::vector<std::tuple<uint16_type,element_wrapper_test_type, element_wrapper_trial_type>> res;
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
                {
                    if ( !entity.isConnectedTo0() )
                        return res;
                    index_type idEntityRangeInit = entity.id();

                    std::set<uint16_type> entityFaceConnectionAvailable;
                    if ( entity.isConnectedTo0() )
                        entityFaceConnectionAvailable.insert( 0 );
                    if ( entity.isConnectedTo1() )
                        entityFaceConnectionAvailable.insert( 1 );
                    if ( entityFaceConnectionAvailable.empty() )
                        return res;

                    std::map<uint16_type,element_wrapper_test_type> eltsRelatedTest;
                    if ( M_rangeAndTestUseSameMesh )
                    {
                        if constexpr ( std::is_same_v< element_mesh_test_type,std::decay_t<decltype(entity.element(0))> > )
                            {
                                for ( uint16_type faceConnectionId : entityFaceConnectionAvailable )
                                    eltsRelatedTest.emplace( faceConnectionId, boost::cref( entity.element( faceConnectionId ) ) );
                            }
                        else
                            CHECK( false ) << "something wrong";
                    }
                    else if ( M_testIsSubMeshFromRange )
                    {
                        index_type eltIdTest = M_spaceTest->mesh()->meshToSubMesh( idEntityRangeInit );
                        if( eltIdTest == invalid_v<index_type> )
                            return res;
                        eltsRelatedTest.emplace( invalid_v<uint16_type>, boost::cref( M_spaceTest->mesh()->element( eltIdTest) ) );
                    }
                    else if ( M_rangeIsSubMeshFromTest )
                    {
                        CHECK( false ) << "TODO";
                        index_type idEltTestInit = entity.mesh()->subMeshToMesh( idEntityRangeInit );
                    }

                    for ( uint16_type faceConnectionId : std::vector<uint16_type>({invalid_v<uint16_type>,0,1}) )
                        {
                            auto itFindConnection = eltsRelatedTest.find( faceConnectionId );
                            if ( itFindConnection == eltsRelatedTest.end() )
                                continue;
                            bool insertElt = true;
                            if ( M_hasMeshSupportPartialTest )
                            {
                                index_type eltIdTest = unwrap_ref( itFindConnection->second ).id();
                                if ( !M_meshSupportTest->hasElement( eltIdTest ) )
                                    insertElt = false;
                            }
                            if ( !insertElt )
                                eltsRelatedTest.erase( faceConnectionId );
                        }


                    std::map<uint16_type,element_wrapper_trial_type> eltsRelatedTrial;
                    if ( M_rangeAndTrialUseSameMesh )
                    {
                        if constexpr ( std::is_same_v< element_mesh_trial_type,std::decay_t<decltype(entity.element(0))> > )
                            {
                                for ( uint16_type faceConnectionId : entityFaceConnectionAvailable )
                                    eltsRelatedTrial.emplace( faceConnectionId, boost::cref( entity.element( faceConnectionId ) ) );
                            }
                        else
                            CHECK( false ) << "something wrong";
                    }
                    else if ( M_trialIsSubMeshFromRange )
                    {
                        index_type eltIdTrial = M_spaceTrial->mesh()->meshToSubMesh( idEntityRangeInit );
                        if( eltIdTrial == invalid_v<index_type> )
                            return res;
                        eltsRelatedTrial.emplace( invalid_v<uint16_type>, boost::cref( M_spaceTrial->mesh()->element( eltIdTrial ) ) );
                    }
                    else if ( M_rangeIsSubMeshFromTrial )
                    {
                        CHECK( false ) << "TODO";
                        index_type idEltTrialInit = entity.mesh()->subMeshToMesh( idEntityRangeInit );
                    }

                    for ( uint16_type faceConnectionId : std::vector<uint16_type>({invalid_v<uint16_type>,0,1}) )
                        {
                            auto itFindConnection = eltsRelatedTrial.find( faceConnectionId );
                            if ( itFindConnection == eltsRelatedTrial.end() )
                                continue;
                            bool insertEltId = true;
                            if ( M_hasMeshSupportPartialTrial )
                            {
                                index_type eltIdTrial = unwrap_ref( itFindConnection->second ).id();
                                if ( !M_meshSupportTrial->hasElement( eltIdTrial ) )
                                    insertEltId = false;
                            }

                            if ( !insertEltId )
                                eltsRelatedTrial.erase( faceConnectionId );
                        }

                    if ( eltsRelatedTest.empty() || eltsRelatedTrial.empty() )
                        return res;
                    bool eltTestHasNoFaceConnection = eltsRelatedTest.find( invalid_v<uint16_type> ) != eltsRelatedTest.end();
                    bool eltTrialHasNoFaceConnection = eltsRelatedTrial.find( invalid_v<uint16_type> ) != eltsRelatedTrial.end();
                    CHECK( eltTestHasNoFaceConnection || eltTrialHasNoFaceConnection ) << "at least one space should not have a connection";
                    if ( eltTestHasNoFaceConnection && eltTrialHasNoFaceConnection )
                    {
                        if ( entity.isGhostFace() ) // test/trial are not connected, do the integration on one proc from the range face property
                            return res;
                        auto const& eltTest = eltsRelatedTest.find( invalid_v<uint16_type> )->second;
                        auto const& eltTrial = eltsRelatedTrial.find( invalid_v<uint16_type> )->second;
                        for ( uint16_type faceConnectionId : entityFaceConnectionAvailable ) // really integrate both side here?
                            res.push_back( std::make_tuple( faceConnectionId, eltTest,eltTrial ) );
                    }
                    else if ( eltTestHasNoFaceConnection )
                    {
                        // check if the face is not at interprocess dof, else need to select one side
                        if ( entityFaceConnectionAvailable.size() == 2 )
                        {
                            bool isGhostFace = entity.isGhostFace();
                            if ( M_rangeAndTrialUseSameMesh )
                            {
                                if ( M_hasMeshSupportPartialTrial )
                                    isGhostFace = M_meshSupportTrial->isGhostFace( entity );
                            }
                            else if ( M_rangeIsSubMeshFromTrial )
                            {
                                CHECK( false ) << "TODO";
                            }
                            if ( isGhostFace )
                                return res;
                        }

                        for ( uint16_type faceConnectionId : entityFaceConnectionAvailable )
                        {
                            auto itFindTrialConnection = eltsRelatedTrial.find( faceConnectionId );
                            if ( itFindTrialConnection != eltsRelatedTrial.end() )
                            {
                                auto const& eltTest = eltsRelatedTest.find( invalid_v<uint16_type> )->second;
                                auto const& eltTrial = itFindTrialConnection->second;
                                res.push_back( std::make_tuple( faceConnectionId, eltTest, eltTrial ) );
                            }
                        }
                    }
                    else if ( eltTrialHasNoFaceConnection )
                    {
                        // check if the face is not at interprocess dof, else need to select one side
                        if ( entityFaceConnectionAvailable.size() == 2 /*&& eltTest.size() == 2*/ )
                        {
                            bool isGhostFace = entity.isGhostFace();
                            if ( M_rangeAndTestUseSameMesh )
                            {
                                if ( M_hasMeshSupportPartialTest )
                                    isGhostFace =  M_meshSupportTest->isGhostFace( entity );
                            }
                            else if ( M_rangeIsSubMeshFromTest )
                            {
                                CHECK( false ) << "TODO";
                            }
                            if ( isGhostFace )
                                return res;
                        }

                        for ( uint16_type faceConnectionId : entityFaceConnectionAvailable )
                        {
                            auto itFindTestConnection = eltsRelatedTest.find( faceConnectionId );
                            if ( itFindTestConnection != eltsRelatedTest.end() )
                            {
                                auto const& eltTest = itFindTestConnection->second;
                                auto const& eltTrial =  eltsRelatedTrial.find( invalid_v<uint16_type> )->second;
                                res.push_back( std::make_tuple( faceConnectionId, eltTest, eltTrial ) );
                            }
                        }

                    }

                } // if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value )
            return res;
        }

    template <size_type gmc_ctx_expr_v,size_type gmc_ctx_test_v,size_type gmc_ctx_trial_v, typename GeoPcExprType, typename ExprType, typename ImTestType, typename ImTrialType, typename GmcExprType, typename GmcTestType, typename GmcTrialType>
    void initGmc( GeoPcExprType & __geopc, ExprType const& expr, ImTestType const& imTest, ImTrialType const& imTrial,
                  std::shared_ptr<GmcExprType>  & gmcExpr, std::shared_ptr<GmcTestType> & gmcFormTest, std::shared_ptr<GmcTrialType> & gmcFormTrial  )
        {
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
            {
                auto const& faceInit = *M_rangeEntityInit;
                uint16_type faceIdInElt = (M_faceConnectionIdInit==0)? faceInit.idInElement0() : faceInit.idInElement1() ;
                gmcExpr = faceInit.element( 0 ).gm()->template context<gmc_ctx_expr_v>( faceInit.element( 0 ), __geopc, faceIdInElt, expr.dynamicContext() );

                if ( M_rangeAndTestUseSameMesh )
                {
                    if constexpr ( std::is_same_v<GmcExprType,GmcTestType> )
                        {
                            gmcFormTest = gmcExpr;
                        }
                    else
                        CHECK( false ) << "not allowed";
                }
                else if ( M_testIsSubMeshFromRange )
                {
                    if constexpr ( GmcTestType::subEntityCoDim == 0 )
                    {
                        typedef typename GmcTestType::gm_type::precompute_type geopc_test_type;
                        auto geopcTest = std::make_shared<geopc_test_type>( M_spaceTest->gm(), imTest.points() );
                        CHECK( M_spaceTest->mesh()->hasElement( M_eltIndexInitTest ) ) << " mesh doesnt have an element id " << M_eltIndexInitTest;
                        auto const& eltInit = M_spaceTest->mesh()->element( M_eltIndexInitTest );
                        gmcFormTest = M_spaceTest->gm()->template context<gmc_ctx_test_v>( eltInit, geopcTest );
                    }
                    else
                        CHECK( false ) << "not allowed";
                }
                else if ( M_rangeIsSubMeshFromTest )
                {
                    CHECK( false ) << "TODO";
                }


                if ( M_rangeAndTrialUseSameMesh )
                {
                    if constexpr ( std::is_same_v<GmcExprType,GmcTrialType> )
                        {
                            gmcFormTrial = gmcExpr;
                        }
                    else
                        CHECK( false ) << "not allowed";
                }
                else if ( M_trialIsSubMeshFromRange )
                {
                    if constexpr ( GmcTrialType::subEntityCoDim == 0 )
                    {
                        typedef typename GmcTrialType::gm_type::precompute_type geopc_trial_type;
                        auto geopcTrial = std::make_shared<geopc_trial_type>( M_spaceTrial->gm(), imTrial.points() );
                        CHECK( M_spaceTrial->mesh()->hasElement( M_eltIndexInitTrial ) ) << " mesh doesnt have an element id " << M_eltIndexInitTrial;
                        auto const& eltInit = M_spaceTrial->mesh()->element( M_eltIndexInitTrial );
                        gmcFormTrial = M_spaceTrial->gm()->template context<gmc_ctx_trial_v>( eltInit, geopcTrial );
                    }
                    else
                        CHECK( false ) << "not allowed";
                }
                else if ( M_rangeIsSubMeshFromTrial )
                {
                    CHECK( false ) << "TODO";
                }
            }
        }


    template <size_type gmc_ctx_expr_v,size_type gmc_ctx_test_v,size_type gmc_ctx_trial_v, typename GmcExprType, typename GmcTestType, typename GmcTrialType>
    void updateGmc( entity_range_t<RangeType> const& entity, std::tuple<uint16_type,element_wrapper_test_type, element_wrapper_trial_type> const& eltsTrialTestRelated,
                    std::shared_ptr<GmcExprType> gmcExpr, std::shared_ptr<GmcTestType> gmcFormTest, std::shared_ptr<GmcTrialType> gmcFormTrial )
        {
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
            {
                uint16_type faceConnection =  std::get<0>( eltsTrialTestRelated );
                uint16_type faceIdInElement = faceConnection==0? entity.idInElement0():entity.idInElement1();
                if ( M_testIsSubMeshFromRange )
                {
                    if constexpr ( GmcTestType::subEntityCoDim == 0 )
                    {
                        auto const& theelt = unwrap_ref( std::get<1>( eltsTrialTestRelated ) );
                        CHECK( gmcFormTest ) << "no gmc defined";
                        gmcFormTest->template update<gmc_ctx_test_v>(theelt);
                    }
                    else
                        CHECK( false ) << "not allowed";
                }
                // TODO optimize if test/trial gmc are the same
                if ( M_trialIsSubMeshFromRange )
                {
                    if constexpr ( GmcTrialType::subEntityCoDim == 0 )
                    {
                        auto const& theelt = unwrap_ref( std::get<2>( eltsTrialTestRelated ) );
                        gmcFormTrial->template update<gmc_ctx_trial_v>(theelt);
                    }
                    else
                        CHECK( false ) << "not allowed";
                }

                // update gmc of the expr on the face
                if ( M_testIsSubMeshFromRange )
                {
                    auto [found,perm] = gmcExpr->template updateFromMatchingNodes<gmc_ctx_expr_v>( entity.element( faceConnection ), faceIdInElement, gmcFormTest );
                    CHECK(found) << "the permutation of quad point was not found";
                }
                else if ( M_trialIsSubMeshFromRange )
                {
                    auto [found,perm] = gmcExpr->template updateFromMatchingNodes<gmc_ctx_expr_v>( entity.element( faceConnection ), faceIdInElement, gmcFormTrial );
                    CHECK(found) << "the permutation of quad point was not found";
                }
                else
                {
                    CHECK( false ) << "not allowed";
                }

            }

        }

    bool hasFoundEntityToInit() const { return M_rangeEntityInit != nullptr; }
    auto const& entityToInit() const { return *M_rangeEntityInit; }

    index_type eltIndexInitTest() const { return M_eltIndexInitTest; }
    index_type eltIndexInitTrial() const { return M_eltIndexInitTrial; }

private :
    std::shared_ptr<SpaceTestType> M_spaceTest;
    std::shared_ptr<SpaceTrialType> M_spaceTrial;
    bool M_hasMeshSupportPartialTest, M_hasMeshSupportPartialTrial;
    bool M_rangeIsSubMeshFromTest, M_rangeIsSubMeshFromTrial, M_testIsSubMeshFromRange, M_trialIsSubMeshFromRange, M_rangeAndTestUseSameMesh, M_rangeAndTrialUseSameMesh;
    typename SpaceTestType::dof_type::mesh_support_ptrtype M_meshSupportTest;
    typename SpaceTrialType::dof_type::mesh_support_ptrtype M_meshSupportTrial;

    entity_range_t<RangeType> const * M_rangeEntityInit;
    uint16_type M_faceConnectionIdInit;
    index_type M_eltIndexInitTest, M_eltIndexInitTrial;
};

template<typename SpaceTestType, typename RangeType>
struct ManageRelationDifferentMeshTypeLinearForm
{
    using index_type = typename SpaceTestType::index_type;

    using element_mesh_test_type = typename MeshTraits<typename SpaceTestType::mesh_type>::element_type;
    using element_wrapper_test_type = typename MeshTraits<typename SpaceTestType::mesh_type>::elements_reference_wrapper_type::value_type;
    static const int range_idim_type = boost::tuples::template element<0,range_t<RangeType> >::type::value;

    ManageRelationDifferentMeshTypeLinearForm( std::shared_ptr<SpaceTestType> const& spaceTest, RangeType const& range )
        :
        M_spaceTest( spaceTest ),
        M_hasMeshSupportPartialTest( spaceTest->dof()->hasMeshSupport() && spaceTest->dof()->meshSupport()->isPartialSupport() ),
        M_rangeIsSubMeshFromTest( false ), M_testIsSubMeshFromRange( false ), M_rangeAndTestUseSameMesh( false ),
        M_rangeEntityInit( nullptr ), M_eltIndexInitTest( invalid_v<index_type> )
        {
            if ( M_hasMeshSupportPartialTest )
                M_meshSupportTest = spaceTest->dof()->meshSupport();

            bool doCheckRelation = false;
            for ( auto const& entityWrap : range )
            {
                auto const& entity = unwrap_ref( entityWrap );
                if ( !doCheckRelation )
                {
                    M_rangeIsSubMeshFromTest = entity.mesh()->isSubMeshFrom( M_spaceTest->mesh() );
                    M_testIsSubMeshFromRange = M_spaceTest->mesh()->isSubMeshFrom( entity.mesh() );
                    M_rangeAndTestUseSameMesh = entity.mesh()->isSameMesh( M_spaceTest->mesh() );
                    doCheckRelation = true;
                }

                auto eltsInit = this->eltsRelatedToRange( entity );
                if ( eltsInit.empty() )
                    continue;
                M_faceConnectionIdInit = std::get<0>( eltsInit.front() ); // only for range of faces
                M_eltIndexInitTest = unwrap_ref( std::get<1>( eltsInit.front() ) ).id();
                M_rangeEntityInit = &entity;

                break;
            }
        }

    std::vector<std::tuple<uint16_type, element_wrapper_test_type>>
    eltsRelatedToRange( entity_range_t<RangeType> const& entity )
        {
            std::vector<std::tuple<uint16_type,element_wrapper_test_type>> res;
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
                {
                    if ( !entity.isConnectedTo0() )
                        return res;
                    index_type idEntityRangeInit = entity.id();

                    std::set<uint16_type> entityFaceConnectionAvailable;
                    if ( entity.isConnectedTo0() )
                        entityFaceConnectionAvailable.insert( 0 );
                    if ( entity.isConnectedTo1() )
                        entityFaceConnectionAvailable.insert( 1 );
                    if ( entityFaceConnectionAvailable.empty() )
                        return res;

                    std::map<uint16_type,element_wrapper_test_type> eltsRelatedTest;
                    if ( M_rangeAndTestUseSameMesh )
                    {
                        if constexpr ( std::is_same_v< element_mesh_test_type,std::decay_t<decltype(entity.element(0))> > )
                            {
                                for ( uint16_type faceConnectionId : entityFaceConnectionAvailable )
                                    eltsRelatedTest.emplace( faceConnectionId, boost::cref( entity.element( faceConnectionId ) ) );
                            }
                        else
                            CHECK( false ) << "something wrong";
                    }
                    else if ( M_testIsSubMeshFromRange )
                    {
                        index_type eltIdTest = M_spaceTest->mesh()->meshToSubMesh( idEntityRangeInit );
                        if( eltIdTest == invalid_v<index_type> )
                            return res;
                        eltsRelatedTest.emplace( invalid_v<uint16_type>, boost::cref( M_spaceTest->mesh()->element( eltIdTest) ) );
                    }
                    else if ( M_rangeIsSubMeshFromTest )
                    {
                        CHECK( false ) << "TODO";
                        index_type idEltTestInit = entity.mesh()->subMeshToMesh( idEntityRangeInit );
                    }

                    for ( uint16_type faceConnectionId : std::vector<uint16_type>({invalid_v<uint16_type>,0,1}) )
                    {
                        auto itFindConnection = eltsRelatedTest.find( faceConnectionId );
                        if ( itFindConnection == eltsRelatedTest.end() )
                            continue;
                        bool insertElt = true;
                        if ( M_hasMeshSupportPartialTest )
                        {
                            index_type eltIdTest = unwrap_ref( itFindConnection->second ).id();
                            if ( !M_meshSupportTest->hasElement( eltIdTest ) )
                                insertElt = false;
                        }
                        if ( !insertElt )
                            eltsRelatedTest.erase( faceConnectionId );
                    }

                    if ( eltsRelatedTest.empty() )
                        return res;
                    bool eltTestHasNoFaceConnection = eltsRelatedTest.find( invalid_v<uint16_type> ) != eltsRelatedTest.end();
                    CHECK( eltTestHasNoFaceConnection ) << "at least one space should not have a connection";
                    if ( eltTestHasNoFaceConnection )
                    {
                        if ( entity.isGhostFace() ) // test/trial are not connected, do the integration on one proc from the range face property
                            return res;
                        auto const& eltTest = eltsRelatedTest.find( invalid_v<uint16_type> )->second;
                        for ( uint16_type faceConnectionId : entityFaceConnectionAvailable ) // really integrate both side here?
                            res.push_back( std::make_tuple( faceConnectionId, eltTest ) );
                    }

                }
            return res;
        }

    template <size_type gmc_ctx_expr_v,size_type gmc_ctx_test_v, typename GeoPcExprType, typename ExprType, typename ImTestType, typename GmcExprType, typename GmcTestType>
    void initGmc( GeoPcExprType & __geopc, ExprType const& expr, ImTestType const& imTest,
                  std::shared_ptr<GmcExprType>  & gmcExpr, std::shared_ptr<GmcTestType> & gmcFormTest  )
        {
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
            {
                auto const& faceInit = *M_rangeEntityInit;
                uint16_type faceIdInElt = (M_faceConnectionIdInit==0)? faceInit.idInElement0() : faceInit.idInElement1() ;
                gmcExpr = faceInit.element( 0 ).gm()->template context<gmc_ctx_expr_v>( faceInit.element( 0 ), __geopc, faceIdInElt, expr.dynamicContext() );

                if ( M_rangeAndTestUseSameMesh )
                {
                    if constexpr ( std::is_same_v<GmcExprType,GmcTestType> )
                        {
                            gmcFormTest = gmcExpr;
                        }
                    else
                        CHECK( false ) << "not allowed";
                }
                else if ( M_testIsSubMeshFromRange )
                {
                    typedef typename GmcTestType::gm_type::precompute_type geopc_test_type;
                    auto geopcTest = std::make_shared<geopc_test_type>( M_spaceTest->gm(), imTest.points() );
                    CHECK( M_spaceTest->mesh()->hasElement( M_eltIndexInitTest ) ) << " mesh doesnt have an element id " << M_eltIndexInitTest;
                    auto const& eltInit = M_spaceTest->mesh()->element( M_eltIndexInitTest );
                    gmcFormTest =  M_spaceTest->gm()->template context<gmc_ctx_test_v>( eltInit, geopcTest );
                }
                else if ( M_rangeIsSubMeshFromTest )
                {
                    CHECK( false ) << "TODO";
                }
            }
        }


    template < size_type gmc_ctx_expr_v,size_type gmc_ctx_test_v, typename GmcExprType, typename GmcTestType>
    void updateGmc( entity_range_t<RangeType> const& entity, std::tuple<uint16_type,element_wrapper_test_type> const& eltsTestRelated,
                    std::shared_ptr<GmcExprType> gmcExpr, std::shared_ptr<GmcTestType> gmcFormTest )
        {
            if constexpr ( range_idim_type == mpl::size_t<MESH_FACES>::value ) // case range of faces
            {
                uint16_type faceConnection =  std::get<0>( eltsTestRelated );
                uint16_type faceIdInElement = faceConnection==0? entity.idInElement0():entity.idInElement1();
                if ( M_testIsSubMeshFromRange )
                {
                    auto const& theelt = unwrap_ref( std::get<1>( eltsTestRelated ) );
                    CHECK( gmcFormTest ) << "no gmc defined";
                    gmcFormTest->template update<gmc_ctx_test_v>(theelt);
                }

                // update gmc of the expr on the face
                if ( M_testIsSubMeshFromRange )
                {
                    auto [found,perm] = gmcExpr->template updateFromMatchingNodes<gmc_ctx_expr_v>( entity.element( faceConnection ), faceIdInElement, gmcFormTest );
                    CHECK(found) << "the permutation of quad point was not found";
                }
                else
                {
                    CHECK( false ) << "not allowed";
                }

            }

        }

    bool hasFoundEntityToInit() const { return M_rangeEntityInit != nullptr; }
    auto const& entityToInit() const { return *M_rangeEntityInit; }

    index_type eltIndexInitTest() const { return M_eltIndexInitTest; }

private :
    std::shared_ptr<SpaceTestType> M_spaceTest;
    bool M_hasMeshSupportPartialTest;
    bool M_rangeIsSubMeshFromTest, M_testIsSubMeshFromRange, M_rangeAndTestUseSameMesh;
    typename SpaceTestType::dof_type::mesh_support_ptrtype M_meshSupportTest;

    entity_range_t<RangeType> const * M_rangeEntityInit;
    uint16_type M_faceConnectionIdInit;
    index_type M_eltIndexInitTest;
};


template<typename SpaceTestType, typename SpaceTrialType, typename RangeType>
auto
manageRelationDifferentMeshType( std::shared_ptr<SpaceTestType> const& spaceTest, std::shared_ptr<SpaceTrialType> const& spaceTrial, RangeType const& range )
{
    return ManageRelationDifferentMeshType<SpaceTestType,SpaceTrialType,RangeType>( spaceTest, spaceTrial, range );
}
template<typename SpaceTestType, typename RangeType>
auto
manageRelationDifferentMeshType( std::shared_ptr<SpaceTestType> const& spaceTest, RangeType const& range )
{
    return ManageRelationDifferentMeshTypeLinearForm<SpaceTestType,RangeType>( spaceTest, range );
}


template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( std::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,std::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::true_ )
{
    return gmcExpr;
}

template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( std::shared_ptr<SpaceType> const& /*space*/,typename GmcType::gm_ptrtype const& /*gm*/, ImType const& im,
                                        size_type /*idElt*/ ,std::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/, mpl::false_ )
{
    CHECK( false ) << "not allowed\n";
    return std::shared_ptr<GmcType>();
}

template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( std::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt ,std::shared_ptr<GmcExprType> gmcExpr,mpl::int_<0> /**/ )
{
    return buildGmcWithRelationDifferentMeshType2<gmc_v,SpaceType,ImType,GmcType,GmcExprType>( space,gm,im, idElt, gmcExpr,
                                                                                         mpl::int_<0>(),
                                                                                         mpl::bool_<boost::is_same<GmcType,GmcExprType>::type::value>() );
}
template<size_type gmc_v,typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
std::shared_ptr<GmcType>
buildGmcWithRelationDifferentMeshType2( std::shared_ptr<SpaceType> const& space,typename GmcType::gm_ptrtype const& gm, ImType const& im,
                                        size_type idElt,std::shared_ptr<GmcExprType> /**/ ,mpl::int_<1> /**/ )
{
    typedef GmcType gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename gmc_type::gm_type::precompute_type geopc_type;
    typedef typename gmc_type::gm_type::precompute_ptrtype geopc_ptrtype;

    geopc_ptrtype geopc( std::make_shared<geopc_type>( gm, im.points() ) );
    CHECK( space->mesh()->hasElement( idElt ) ) << " mesh doesnt have an element id " << idElt;
    auto const& eltInit = space->mesh()->element( idElt );
    return gm->template context<gmc_v>( eltInit, geopc );
}

template<size_type gmc_ctx_v,size_type gmc_ctx_expr_v,typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
typename FaceType::size_type
updateGmcWithRelationDifferentMeshType2( FaceType const& theface, std::shared_ptr<SpaceType> const& /*space*/,
                                         std::shared_ptr<GmcType> /*gmc*/, std::shared_ptr<GmcExprType> /*gmcExpr*/,
                                         size_type idElt, mpl::int_<0> /**/ )
{
    // nothing to do!
    return invalid_v<typename FaceType::size_type>;
}
template<typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
typename FaceType::size_type
updateGmcWithRelationDifferentMeshType21( FaceType const& theface, std::shared_ptr<SpaceType> const& /*space*/,
                                          std::shared_ptr<GmcType> /*gmc*/, std::shared_ptr<GmcExprType> /*gmcExpr*/,
                                          size_type /*idElt*/, mpl::int_<0> /**/ )
{
    // nothing to do!
    return invalid_v<typename FaceType::size_type>;
}

template<size_type gmc_ctx_v,size_type gmc_ctx_expr_v,typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
typename FaceType::size_type
updateGmcWithRelationDifferentMeshType2( FaceType const& theface, std::shared_ptr<SpaceType> const& space,
                                         std::shared_ptr<GmcType>& gmc, std::shared_ptr<GmcExprType>& gmcExpr,
                                         size_type idElt, mpl::int_<1> /**/ )
{
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;

    CHECK( space->mesh()->hasElement( idElt ) ) << " mesh doesnt have an element id " << idElt;
    auto const& theelt = space->mesh()->element( idElt );
    gmc->template update<gmc_ctx_v>(theelt);
    auto [found,perm] = gmcExpr->template updateFromMatchingNodes<gmc_ctx_expr_v>( theface.element0(), theface.idInElement0(), gmc );
    CHECK(found) << "the permutation of quad point was not found\n";
    return idElt;
}

template<typename FaceType, typename SpaceType,typename ImType,typename GmcType,typename GmcExprType>
size_type
updateGmcWithRelationDifferentMeshType21( FaceType const& theface, std::shared_ptr<SpaceType> const& space,
                                          std::shared_ptr<GmcType> gmc, std::shared_ptr<GmcExprType> gmcExpr,
                                          size_type idElt, mpl::int_<1> /**/ )
{
    typedef typename QuadMapped<ImType>::permutation_type permutation_type;

    CHECK( space->mesh()->hasElement( idElt ) ) << " mesh doesnt have an element id " << idElt;
    auto const& theelt = space->mesh()->element( idElt );
    gmc->update(theelt);
    auto [found,perm] = gmcExpr->updateFromMatchingNodes( theface.element1(), theface.idInElement1(), gmc );
    //bool found = gmcExpr->updateFromNeighborMatchingFace( theface.element0(), theface.idInElement0(), gmc );
    CHECK(found) << "the permutation of quad point was not found\n";
    return idElt;
}

} // namespace detail

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                           mpl::int_<MESH_ELEMENTS> /**/, mpl::false_ ) const
{
    DLOG(INFO) << "[integrator::assembleWithRelationDifferentMeshType] vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS>\n";

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    constexpr uint16_type nDimRange = gm_expr_type::nDim;
    constexpr size_type gmc_context_expr_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_expr_type::template Context<typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<typename eval::element_type> gmc1_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef std::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
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
    constexpr uint16_type nDimTest = gm_formTest_type::nDim;
    constexpr uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    constexpr size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc1_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef std::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    constexpr uint16_type nDimTrial = gm_formTrial_type::nDim;
    constexpr uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    constexpr size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc1_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef std::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //static const bool gmTestIsGmExpr = boost::is_same<gm_expr_type,gm_formTest_type>::type::value;
    //static const bool gmTrialIsGmExpr = boost::is_same<gm_expr_type,gm_formTrial_type>::type::value;


    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im1_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im1_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;
    using im1_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;

    im_range_type imRange( M_im.order() );
    im1_range_type im1Range( M_im2.order() );
    im_formtest_type imTest( M_im.order() );
    im_formtrial_type imTrial( M_im.order() );

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

            auto const& eltInit = boost::unwrap_ref( *elt_it );

            const bool rangeIsSubMeshFromTest = eltInit.mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
            const bool rangeIsSubMeshFromTrial = eltInit.mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
            const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( eltInit.mesh() );
            const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( eltInit.mesh() );

            const size_type idEltRangeInit = eltInit.id();
            size_type idEltTestInit  = idEltRangeInit;
            if ( rangeIsSubMeshFromTest )
                idEltTestInit = eltInit.mesh()->subMeshToMesh( idEltRangeInit );
            else if ( testIsSubMeshFromRange )
                idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            size_type idEltTrialInit = idEltRangeInit;
            if ( rangeIsSubMeshFromTrial )
                idEltTrialInit = eltInit.mesh()->subMeshToMesh( idEltRangeInit );
            else if ( trialIsSubMeshFromRange )
                idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            //-----------------------------------------------//
            pc_expr_ptrtype geopcExpr( new pc_expr_type( eltInit.gm(), imRange.points() ) );
            pc1_expr_ptrtype geopc1Expr( new pc1_expr_type( eltInit.gm1(), im1Range.points() ) );
            gmc_expr_ptrtype gmcExpr = eltInit.gm()->template context<gmc_context_expr_v>( eltInit, geopcExpr );
            gmc1_expr_ptrtype gmc1Expr = eltInit.gm1()->template context<gmc_context_expr_v>( eltInit, geopc1Expr );
            //-----------------------------------------------//
            auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType< gmc_context_formTest_v, typename FormType::space_1_type,im_formtest_type,
                                                                              gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                                idEltTestInit, gmcExpr, mpl::int_<gmTestRangeRelation>() );
            auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType< gmc_context_formTrial_v, typename FormType::space_2_type,im_formtrial_type,
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
                    auto const& eltCur = boost::unwrap_ref( *elt_it );

                    const size_type idEltRange = eltCur.id();
                    size_type idEltTest = idEltRange;
                    if ( rangeIsSubMeshFromTest )
                        idEltTest = eltCur.mesh()->subMeshToMesh( idEltRange );
                    else if ( testIsSubMeshFromRange )
                        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
                    size_type idEltTrial = idEltRange;
                    if ( rangeIsSubMeshFromTrial )
                        idEltTrial = eltCur.mesh()->subMeshToMesh( idEltRange );
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
                        gmcExpr->template update<gmc_context_expr_v>(eltCur);
                        detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTest_v,typename FormType::space_1_type,im_formtest_type,
                                                                       gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                        idEltTest, mpl::int_<gmTestRangeRelation>() );
                        detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTrial_v,typename FormType::space_2_type,im_formtrial_type,
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
                            gmcExpr->template update<gmc_context_expr_v>(eltCur);
                            detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTest_v,typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                            idEltTest, mpl::int_<gmTestRangeRelation>() );
                            detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTrial_v,typename FormType::space_2_type,im_formtrial_type,
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
    DLOG(INFO) << "[integrator::assembleWithRelationDifferentMeshType(ELEMENTS)] mortar case\n";

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const size_type gmc_context_expr_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_expr_type::template Context<typename eval::element_type> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<typename eval::element_type> gmc1_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef std::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
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
    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc1_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef std::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type gmTrialRangeRelation = ( nDimTrial > nDimRange )? nDimTrial-nDimRange : nDimRange-nDimTrial;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc1_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef std::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //static const bool gmTestIsGmExpr = boost::is_same<gm_expr_type,gm_formTest_type>::type::value;
    //static const bool gmTrialIsGmExpr = boost::is_same<gm_expr_type,gm_formTrial_type>::type::value;

    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im1_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im1_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;
    using im1_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;

    im_range_type imRange( M_im.order() );
    im1_range_type im1Range( M_im2.order() );
    im_formtest_type imTest( M_im.order() );
    im_formtrial_type imTrial( M_im.order() );


    // mortar context
    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const bool has_mortar_trial = FormType::trial_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test && has_mortar_trial)? 3 : ( (has_mortar_test)? 1 : ( (has_mortar_trial)? 2 : 0 ) );
    BOOST_MPL_ASSERT_MSG( mortarTag < 3,TODO_CASE3, (mpl::int_<mortarTag>) );

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type,mortarTag> form_mortar_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, im1_range_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, im1_range_type, map_gmc1_expr_type, map_gmc1_formTrial_type,mortarTag> form1_mortar_context_type;
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

            auto const& eltInit = boost::unwrap_ref( *elt_it );

            const bool rangeIsSubMeshFromTest = eltInit.mesh()->isSubMeshFrom( __form.testSpace()->mesh() );
            const bool rangeIsSubMeshFromTrial = eltInit.mesh()->isSubMeshFrom( __form.trialSpace()->mesh() );
            const bool testIsSubMeshFromRange = __form.testSpace()->mesh()->isSubMeshFrom( eltInit.mesh() );
            const bool trialIsSubMeshFromRange = __form.trialSpace()->mesh()->isSubMeshFrom( eltInit.mesh() );

            const size_type idEltRangeInit = eltInit.id();
            size_type idEltTestInit  = idEltRangeInit;
            if ( rangeIsSubMeshFromTest )
                idEltTestInit = eltInit.mesh()->subMeshToMesh( idEltRangeInit );
            else if ( testIsSubMeshFromRange )
                idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            size_type idEltTrialInit = idEltRangeInit;
            if ( rangeIsSubMeshFromTrial )
                idEltTrialInit = eltInit.mesh()->subMeshToMesh( idEltRangeInit );
            else if ( trialIsSubMeshFromRange )
                idEltTrialInit = __form.trialSpace()->mesh()->meshToSubMesh( idEltRangeInit );
            //-----------------------------------------------//
            pc_expr_ptrtype geopcExpr( new pc_expr_type( eltInit.gm(), imRange.points() ) );
            pc1_expr_ptrtype geopc1Expr( new pc1_expr_type( eltInit.gm1(), im1Range.points() ) );
            gmc_expr_ptrtype gmcExpr = eltInit.gm()->template context<gmc_context_expr_v>( eltInit, geopcExpr, this->expression().dynamicContext() );
            gmc1_expr_ptrtype gmc1Expr =  eltInit.gm1()->template context<gmc_context_expr_v>( eltInit, geopc1Expr, this->expression().dynamicContext() );
            //gmc_expr_ptrtype gmcExpr( new gmc_expr_type( eltInit.gm(),eltInit, geopcExpr, invalid_uint16_type_value, this->expression().dynamicContext() ) );
            //gmc1_expr_ptrtype gmc1Expr( new gmc1_expr_type( eltInit.gm1(),eltInit, geopc1Expr, invalid_uint16_type_value, this->expression().dynamicContext() ) );
            //-----------------------------------------------//
            auto gmcFormTest = detail::buildGmcWithRelationDifferentMeshType< gmc_context_formTest_v, typename FormType::space_1_type,im_formtest_type,
                                                                              gmc_formTest_type,gmc_expr_type>( __form.testSpace(), __form.testSpace()->gm(), imTest,
                                                                                                                idEltTestInit, gmcExpr, mpl::int_<gmTestRangeRelation>() );
            auto gmcFormTrial = detail::buildGmcWithRelationDifferentMeshType< gmc_context_formTrial_v, typename FormType::space_2_type,im_formtrial_type,
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
                    auto const& eltCur = boost::unwrap_ref( *elt_it );
                    const size_type idEltRange = eltCur.id();
                    size_type idEltTest = idEltRange;
                    if ( rangeIsSubMeshFromTest )
                        idEltTest = eltCur.mesh()->subMeshToMesh( idEltRange );
                    else if ( testIsSubMeshFromRange )
                        idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltRange );
                    size_type idEltTrial = idEltRange;
                    if ( rangeIsSubMeshFromTrial )
                        idEltTrial = eltCur.mesh()->subMeshToMesh( idEltRange );
                    else if ( trialIsSubMeshFromRange )
                        idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltRange );

#if 0
                    if ( formc->isZero( eltTest.element0() /*idElt*/ ) )
                        continue;
#endif
                    DLOG(INFO) << "elt "  << eltCur.id() << " bdy = " << eltCur.isOnBoundary();
                    switch ( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        gmcExpr->template update<gmc_context_expr_v>(eltCur);
                        detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTest_v,typename FormType::space_1_type,im_formtest_type,
                                                                       gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                        idEltTest, mpl::int_<gmTestRangeRelation>() );
                        detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTrial_v,typename FormType::space_2_type,im_formtrial_type,
                                                                       gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                         idEltTrial, mpl::int_<gmTrialRangeRelation>() );
                        if ( ( has_mortar_test && __form.testSpace()->mesh()->element( gmcFormTest->id() ).isOnBoundary() ) ||
                             ( has_mortar_trial && __form.trialSpace()->mesh()->element( gmcFormTrial->id() ).isOnBoundary() ) )
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
                            gmcExpr->template update<gmc_context_expr_v>(eltCur);
                            detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTest_v,typename FormType::space_1_type,im_formtest_type,
                                                                           gmc_formTest_type,gmc_expr_type>(__form.testSpace(), gmcFormTest, gmcExpr,
                                                                                                            idEltTest, mpl::int_<gmTestRangeRelation>() );
                            detail::updateGmcWithRelationDifferentMeshType<gmc_context_formTrial_v,typename FormType::space_2_type,im_formtrial_type,
                                                                           gmc_formTrial_type,gmc_expr_type>(__form.trialSpace(), gmcFormTrial, gmcExpr,
                                                                                                             idEltTrial, mpl::int_<gmTrialRangeRelation>() );

                            if ( ( has_mortar_test && __form.testSpace()->mesh()->element( gmcFormTest->id() ).isOnBoundary() ) ||
                                 ( has_mortar_trial && __form.trialSpace()->mesh()->element( gmcFormTrial->id() ).isOnBoundary() ) )
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
    DLOG(INFO) << "[integrator::assembleInCaseOfInterpolate] vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS>\n";

    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_expr_type::template Context<typename eval::element_type> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type> gmc_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type> gmc_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
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

    auto const& eltInit = boost::unwrap_ref( *elt_it );

    //-----------------------------------------------//
    im_range_type imRange( M_im.order() );
    pc_expr_ptrtype geopcExpr( new pc_expr_type( eltInit.gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr = eltInit.gm()->template context<gmc_context_expr_v>( eltInit, geopcExpr );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
    //-----------------------------------------------//
    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    gmc_formTest_ptrtype gmcFormTest = __form.gm()->template context<gmc_context_formTest_v>( __form.testSpace()->mesh()->element( 0 ), geopcFormTest );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
    //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial = __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial );
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
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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

                    gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                    gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );

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
    DLOG(INFO) << "[integrator::assembleInCaseOfInterpolate] (elements) mortar case";
    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_expr_type::template Context<typename eval::element_type> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    // typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type> gmc_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type> gmc_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    // mortar context
    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const bool has_mortar_trial = FormType::trial_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test && has_mortar_trial)? 3 : ( (has_mortar_test)? 1 : ( (has_mortar_trial)? 2 : 0 ) );
    BOOST_MPL_ASSERT_MSG( mortarTag < 3,TODO_CASE3, (mpl::int_<mortarTag>) );

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_range_type, map_gmc_expr_type, map_gmc_formTrial_type, mortarTag> form_mortar_context_type;
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

    auto const& eltInit = boost::unwrap_ref( *elt_it );

    //-----------------------------------------------//
    im_range_type imRange( M_im.order() );
    pc_expr_ptrtype geopcExpr( new pc_expr_type( eltInit.gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr = eltInit.gm()->template context<gmc_context_expr_v>( eltInit, geopcExpr );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );
    //-----------------------------------------------//
    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), __form.template testFiniteElement<false>()->points() ) );
    VLOG(2) << "pts non mortar : " << __form.template testFiniteElement<false>()->points();
    gmc_formTest_ptrtype gmcFormTest = __form.gm()->template context<gmc_context_formTest_v>( __form.testSpace()->mesh()->beginElement()->second, geopcFormTest );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
    // mortar
    pc_formTest_ptrtype geopcFormTestMortar( new pc_formTest_type( __form.gm(), __form.template testFiniteElement<FormType::test_space_type::is_mortar>()->points() ) );
    VLOG(2) << "pts mortar : " << __form.template testFiniteElement<true>()->points();
    gmc_formTest_ptrtype gmcFormTestMortar = __form.gm()->template context<gmc_context_formTest_v>(  __form.testSpace()->mesh()->beginElement()->second, geopcFormTestMortar );
    map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );
    //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial = __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrial );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
    //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrialMortar( new pc_formTrial_type( __form.gmTrial(), __form.template trialFiniteElement<FormType::trial_space_type::is_mortar>()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrialMortar = __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrialMortar );
    map_gmc_formTrial_type mapgmcFormTrialMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrialMortar ) );

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               imRange ) );
    mortar_focb_ptrtype formcmTest( new form_mortar_context_type( __form,
                                                                  mapgmcFormTestMortar,
                                                                  mapgmcFormTrial,
                                                                  mapgmcExpr,
                                                                  this->expression(),
                                                                  imRange ) );

    mortar_focb_ptrtype formcmTrial( new form_mortar_context_type( __form,
                                                                   mapgmcFormTest,
                                                                   mapgmcFormTrialMortar,
                                                                   mapgmcExpr,
                                                                   this->expression(),
                                                                   imRange ) );

    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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
                bool isFirstExperience = true, isFirstExperienceMortarTest = true, isFirstExperienceMortarTrial = true;
                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    if ( mortarTag == 1 && eltTest.isOnBoundary() )
                    {
                        geopcFormTestMortar->update( gmcExpr_it->template get<2>() );
                        geopcFormTrial->update( gmcExpr_it->template get<3>() );
                        gmcFormTestMortar->template update<gmc_context_formTest_v>( eltTest,geopcFormTestMortar );
                        gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );
                        map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );
                        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                        formcmTest->updateInCaseOfInterpolate( mapgmcFormTestMortar, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->template get<0>() );
                        formcmTest->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperienceMortarTest );
                        isFirstExperienceMortarTest = false;
                    }
                    else if ( mortarTag == 2 && eltTrial.isOnBoundary() )
                    {
                        geopcFormTest->update( gmcExpr_it->template get<2>() );
                        geopcFormTrialMortar->update( gmcExpr_it->template get<3>() );
                        gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                        gmcFormTrialMortar->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrialMortar );
                        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                        map_gmc_formTrial_type mapgmcFormTrialMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrialMortar ) );
                        formcmTrial->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrialMortar, mapgmcExpr, gmcExpr_it->template get<0>() );
                        formcmTrial->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperienceMortarTrial );
                        isFirstExperienceMortarTrial = false;
                    }
                    else
                    {
                        geopcFormTest->update( gmcExpr_it->template get<2>() );
                        geopcFormTrial->update( gmcExpr_it->template get<3>() );
                        gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                        gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );
                        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                        formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->template get<0>() );
                        formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                        isFirstExperience = false;
                    }
                }
                if ( !isFirstExperience )
                    formc->assembleInCaseOfInterpolate();
                if ( !isFirstExperienceMortarTest )
                    formcmTest->assembleInCaseOfInterpolate();
                if ( !isFirstExperienceMortarTrial )
                    formcmTrial->assembleInCaseOfInterpolate();
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
    using im_range_type = im_t<typename eval::the_element_type,expression_value_type>;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_expr_type::template Context<typename eval::element_type> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    static const size_type gmc_context_form_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_form_type::template Context<geoelement_form_type> gmc_form_type;
    typedef std::shared_ptr<gmc_form_type> gmc_form_ptrtype;
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

    auto const& eltInit = boost::unwrap_ref( *elt_it );

    //-----------------------------------------------//

    im_range_type imRange( M_im.order() );
    pc_expr_ptrtype geopcExpr( new pc_expr_type( eltInit.gm(), imRange.points() ) );
    gmc_expr_ptrtype gmcExpr = eltInit.gm()->template context<gmc_context_expr_v>( eltInit, geopcExpr );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    gmc_form_ptrtype gmcForm =  __form.gm()->template context<gmc_context_form_v>( __form.testSpace()->mesh()->element( 0 ), geopcForm );
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
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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
                gmcForm->template update<gmc_context_form_v>( eltTest,geopcForm );
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
boost::tuple<typename Integrator<Elements, Im, Expr, Im2>::size_type,rank_type,uint16_type>
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
        if ( idEltTest == invalid_v<size_type> )
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

    if ( idEltTest != invalid_v<size_type> )
    {
        size_type idEltTrial = idEltTest;
        if ( __form.trialSpace()->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
            idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltTest );
        else if ( __form.testSpace()->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() ) )
            idEltTrial = __form.testSpace()->mesh()->subMeshToMesh( idEltTest );

        if (idEltTrial != invalid_v<size_type>)
            trialEltIsOk=true;
     }

    }
    else // ghost cell
    {
        if ( !faceRange.isConnectedTo1() || faceRange.element1().isGhostCell() )
            return boost::make_tuple( invalid_v<size_type>, procIdElt0, __face_id_in_elt_0 );
        else
            idEltTest = invalid_v<size_type>; // continue algo
    }

    // if first pass not good, restart with faceRange.element1
    if ( ( idEltTest == invalid_v<size_type> || !trialEltIsOk )  && faceRange.isConnectedTo1() && !faceRange.element1().isGhostCell() )
    {
        __face_id_in_elt_0 = faceRange.pos_second();
        procIdElt0 = faceRange.proc_second();
        idEltTest = faceRange.element( 1 ).id();
        if ( faceRange.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
            idEltTest = faceRange.mesh()->subMeshToMesh( idEltTest );
        else if ( __form.testSpace()->mesh()->isSubMeshFrom( faceRange.mesh() ) )
            idEltTest = __form.testSpace()->mesh()->meshToSubMesh( idEltTest );

        trialEltIsOk = false;
        if ( idEltTest != invalid_v<size_type> )
        {
            size_type idEltTrial = idEltTest;
            if ( __form.trialSpace()->mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                idEltTrial = __form.trialSpace()->mesh()->meshToSubMesh( idEltTest );
            else if ( __form.testSpace()->mesh()->isSubMeshFrom( __form.trialSpace()->mesh() ) )
                idEltTrial = __form.testSpace()->mesh()->subMeshToMesh( idEltTest );

            if (idEltTrial != invalid_v<size_type>)
                {

                trialEltIsOk=true;
                }
        }
    }

    // if test or trial id not find, return invalid value
    if ( idEltTest == invalid_v<size_type> || !trialEltIsOk )
        return boost::make_tuple( invalid_v<size_type>, procIdElt0, __face_id_in_elt_0 );
    else
        return boost::make_tuple( idEltTest,procIdElt0, __face_id_in_elt_0 );
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType,typename FaceRangeType>
boost::tuple<typename Integrator<Elements, Im, Expr, Im2>::size_type,rank_type,uint16_type>
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
        if ( idEltTest == invalid_v<size_type> )
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

        if ( idEltTrial == invalid_v<size_type> )
            res=false;
        else
        {
            auto const& faceTrial = __form.trialSpace()->mesh()->element(idEltTrial ).face( faceIdInElt );
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
    tic();
    //
    // some typedefs
    //
    typedef typename eval::gm_type gm_type;
    typedef typename eval::gm1_type gm1_type;
    //typedef typename FormType::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;
    static const size_type gmc_context_face_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_type::template Context<typename eval::element_type,1> gmc_type;
    typedef typename gm1_type::template Context<typename eval::element_type,1> gmc1_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef std::shared_ptr<gmc1_type> gmc1_ptrtype;
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

    bool hasMeshSupportPartialTest = __form.testSpace()->dof()->hasMeshSupport() && __form.testSpace()->dof()->meshSupport()->isPartialSupport();

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        DLOG(INFO) << "Standard integration over "
                   << std::distance( it, en )  << " faces\n";

        // check that we have elements to iterate over
        if ( it == en )
            continue;
        auto const& faceInit = boost::unwrap_ref( *it );

        if ( faceInit.isConnectedTo0() == false )
            continue;

        // true if range/test/trial are same mesh
        bool useSameMesh = this->useSameMesh( __form,faceInit );

        uint16_type __face_id_in_elt_0 = faceInit.pos_first();
        rank_type procIdElt0 = faceInit.proc_first();
        size_type idEltTestInit = faceInit.element( 0 ).id();
        if ( !useSameMesh )
        {
            bool hasFindEltToInit = false;
            while( !hasFindEltToInit)
            {
                if ( faceInit.mesh()->isSubMeshFrom( __form.testSpace()->mesh() ) )
                {
                    idEltTestInit = faceInit.mesh()->subMeshToMesh( idEltTestInit );
                }
                else if ( __form.testSpace()->mesh()->isSubMeshFrom( faceInit.mesh() ) )
                {
                    idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( idEltTestInit );
                    if ( idEltTestInit == invalid_v<size_type> )
                    {
                        if ( faceInit.isConnectedTo1() )
                        {
                            __face_id_in_elt_0 = faceInit.pos_second();
                            procIdElt0 = faceInit.proc_second();
                            idEltTestInit = __form.testSpace()->mesh()->meshToSubMesh( faceInit.element( 1 ).id() );
                        }
                    }
                }
                if ( idEltTestInit == invalid_v<size_type> )
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
                CHECK( idEltTestInit != invalid_v<size_type> ) << "mesh relation fail : no find a corresponding element\n";

        }

        auto const& elt0TestInit = __form.testSpace()->mesh()->element( idEltTestInit );
        //auto const& faceTestInit = elt0TestInit.face( __face_id_in_elt_0 );

        // get the geometric mapping associated with element 0
        //DLOG(INFO) << "element " << faceInit.element(0)  << "face " << __face_id_in_elt_0 << " permutation " << faceInit.element(0).permutation( __face_id_in_elt_0 ) << "\n";
        gm_ptrtype __gm = elt0TestInit.gm();
        gm1_ptrtype __gm1 = elt0TestInit.gm1();
        //DLOG(INFO) << "[integrator] evaluate(faces), gm is cached: " << __gm->isCached() << "\n";
        gmc_ptrtype __c0 = __gm->template context<gmc_context_face_v>( elt0TestInit, __geopc, __face_id_in_elt_0, this->expression().dynamicContext() );
        gmc1_ptrtype __c01 = __gm1->template context<gmc_context_face_v>( elt0TestInit, __geopc1, __face_id_in_elt_0, this->expression().dynamicContext() );

        //
        // the case where the face is connected only to one element
        //
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
        typedef typename FormType::template Context<map_gmc_type, expression_type, face_im_type> form_context_type;
        typedef typename FormType::template Context<map_gmc1_type, expression_type, face_im2_type> form1_context_type;
        typedef std::shared_ptr<form_context_type> form_context_ptrtype;
        typedef std::shared_ptr<form1_context_type> form1_context_ptrtype;
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
        typedef std::shared_ptr<form2_context_type> form2_context_ptrtype;
        typedef std::shared_ptr<form21_context_type> form21_context_ptrtype;
        form2_context_ptrtype form2;
        form21_context_ptrtype form21;

        bool isInitConnectionTo0=false;
        bool isInitConnectionTo1=false;

        // true if connected to another element, false otherwise
        //if ( faceInit.isConnectedTo1() )
        if ( this->faceIntegratorUseTwoConnections(__form,faceInit,elt0TestInit,__face_id_in_elt_0) )
        {
            uint16_type __face_id_in_elt_1 = faceInit.pos_second();
            rank_type procIdElt1 = faceInit.proc_second();
            size_type idElt1TestInit = faceInit.element( 1 ).id();
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
            //CHECK( idElt1TestInit != invalid_v<size_type> ) << "mesh relation fail : no find a corresponding element\n";
            // get element1
            auto const& elt1TestInit = __form.testSpace()->mesh()->element( idElt1TestInit );

            // init linear/bilinear form for two connections
            __c1 = __gm->template context<gmc_context_face_v>( elt1TestInit, __geopc, __face_id_in_elt_1, this->expression().dynamicContext() );
            __c11 = __gm1->template context<gmc_context_face_v>( elt1TestInit, __geopc1, __face_id_in_elt_1, this->expression().dynamicContext() );

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
            auto const& faceCur = boost::unwrap_ref( *it );

            // get some info about test mesh element0 connected to the current face
            uint16_type __face_id_in_elt_0 = faceCur.pos_first();
            rank_type procIdElt0 = faceCur.proc_first();
            size_type idEltTest = faceCur.element( 0 ).id();
            bool swapElt0Elt1WithSameMesh = false;
            bool tryTwoConnections = true;
            if ( !useSameMesh)
            {
                boost::tie( idEltTest, procIdElt0, __face_id_in_elt_0) = this->testElt0IdFromFaceRange(__form,faceCur);
                if ( idEltTest == invalid_v<size_type> )
                    continue;
            }
            else if ( faceCur.isConnectedTo1() )
            {
                if ( !hasMeshSupportPartialTest )
                {
                    if ( faceCur.element0().isGhostCell() )
                    {
                        __face_id_in_elt_0 = faceCur.pos_second();
                        procIdElt0 = faceCur.proc_second();
                        idEltTest = faceCur.element( 1 ).id();
                        swapElt0Elt1WithSameMesh = true;
                    }
                }
                else
                {
                    bool partialSupportHasConnection0 = __form.testSpace()->dof()->meshSupport()->hasElement( idEltTest );
                    size_type idEltOtherConnection = faceCur.element( 1 ).id();
                    bool partialSupportHasConnection1 = __form.testSpace()->dof()->meshSupport()->hasElement( idEltOtherConnection );

                    if ( ( partialSupportHasConnection0 && !partialSupportHasConnection1 ) ||
                         ( !partialSupportHasConnection0 && partialSupportHasConnection1 ) )
                        tryTwoConnections = false;

                    if ( tryTwoConnections )
                    {
                        if ( faceCur.element0().isGhostCell() )
                        {
                            __face_id_in_elt_0 = faceCur.pos_second();
                            procIdElt0 = faceCur.proc_second();
                            idEltTest = faceCur.element( 1 ).id();
                            swapElt0Elt1WithSameMesh = true;
                        }
                    }
                    else
                    {
                        if ( partialSupportHasConnection0 )
                        {
                            if ( faceCur.element0().isGhostCell() )
                                continue;
                        }
                        else
                        {
                            if ( faceCur.element1().isGhostCell() )
                                continue;
                            __face_id_in_elt_0 = faceCur.pos_second();
                            procIdElt0 = faceCur.proc_second();
                            idEltTest = idEltOtherConnection;//faceCur.element( 1 ).id();
                            //swapElt0Elt1WithSameMesh = true;
                        }
                    }
                }
            }

            // element0 (from test mesh) used in integration
            auto const& elt0Test = __form.testSpace()->mesh()->element( idEltTest );
            CHECK( !elt0Test.isGhostCell() ) << "elt0 can't be a ghost element";
            //auto const& faceTest = elt0Test.face( __face_id_in_elt_0 );


            //if ( faceCur.isConnectedTo1() )
            if ( tryTwoConnections && this->faceIntegratorUseTwoConnections(__form,faceCur,elt0Test,__face_id_in_elt_0) )
            {
                if ( faceCur.isGhostFace() )
                {
                    if ( hasMeshSupportPartialTest )
                    {
                        if ( __form.testSpace()->dof()->meshSupport()->isGhostFace( faceCur ) )
                            continue;
                    }
                    else
                    {
                        LOG(WARNING) << "face id : " << faceCur.id() << " is a ghost face" << faceCur.G();
                        continue;
                    }
                }
                // // if is a interprocess faces, only integrate in one process
                // if ( faceCur.isInterProcessDomain() && faceCur.partition1() > faceCur.partition2() )
                //     continue;


                // get some info about test mesh element1 connected to the current face
                uint16_type __face_id_in_elt_1 = faceCur.pos_second();
                rank_type procIdElt1 = faceCur.proc_second();
                size_type idElt1Test = faceCur.element( 1 ).id();
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
                    __face_id_in_elt_1 = faceCur.pos_first();
                     procIdElt1 = faceCur.proc_first();
                     idElt1Test = faceCur.element( 0 ).id();
                }
                CHECK( idElt1Test != invalid_v<size_type> ) << "mesh relation fail : no find a corresponding element\n";

                // element1 (from test mesh) used in integration
                auto const& elt1Test = __form.testSpace()->mesh()->element( idElt1Test );

                if ( !isInitConnectionTo1 )
                {
                    // init linear/bilinear form for element1
                    __c1 = __gm->template context<gmc_context_face_v>( elt1Test, __geopc, __face_id_in_elt_1, this->expression().dynamicContext() );
                    __c11 = __gm1->template context<gmc_context_face_v>( elt1Test, __geopc1, __face_id_in_elt_1, this->expression().dynamicContext() );
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
                    FEELPP_ASSERT( faceCur.isOnBoundary() == false  )
                        ( faceCur.id() ).error( "face on boundary but connected on both sides" );
                    //ti0.restart();
                    __c0->template update<gmc_context_face_v>( elt0Test, __face_id_in_elt_0 );
                    bool found_permutation = __c1->template updateFromNeighborMatchingFace<gmc_context_face_v>( elt1Test, __face_id_in_elt_1, __c0 );
                    CHECK(found_permutation) << "the permutation of quadrature points were not found\n";
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
                    FEELPP_ASSERT( faceCur.isOnBoundary() == false  )
                        ( faceCur.id() ).error( "face on boundary but connected on both sides" );
                    //ti0.restart();
                    __c01->template update<gmc_context_face_v>( elt0Test, __face_id_in_elt_0 );
                    bool found_permutation = __c11->template updateFromNeighborMatchingFace<gmc_context_face_v>( elt1Test, __face_id_in_elt_1, __c01 );
                    CHECK(found_permutation) << "the permutation of quadrature points was not found\n";

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
                if ( !isInitConnectionTo0 )
                {
                    form = form_context_ptrtype( new form_context_type( __form, mapgmc, mapgmc, mapgmc, expression(), face_ims[__face_id_in_elt_0], this->im() ) );
                    form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                    isInitConnectionTo0=true;
                }

                //ti0.restart();
                __c0->template update<gmc_context_face_v>( elt0Test,__face_id_in_elt_0 );
                //t0 += ti0.elapsed();

                //ti1.restart();
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                form->update( mapgmc, mapgmc, mapgmc, face_ims[__face_id_in_elt_0] );
                //t1 += ti1.elapsed();

                //ti2.restart();
                form->integrate();
                //t2 += ti2.elapsed();

                //ti3.restart();
                form->assemble();
                //t3 += ti3.elapsed();
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
    toc("integrating over faces", FLAGS_v>1);
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
                                                                           mpl::int_<MESH_FACES> /**/ /*, mpl::false_*/ ) const
{

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const size_type gmc_context_expr_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_expr_type::template Context<typename eval::element_type,1> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<typename eval::element_type,1> gmc1_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef std::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
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
    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest == nDimRange )? 1 : nDimTest - eval::range_entity_type::nDim;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type contextTest = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTest_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT|vm::JACOBIAN> >::type::value;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<geoelement_formTest_type,gmTestRangeRelation> gmc1_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef std::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;

    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::gm1_2_type gm1_formTrial_type;
    static const uint16_type nDimTrial = gm_formTrial_type::nDim;
    static const uint16_type gmTrialRangeRelation =  ( nDimTrial == nDimRange )? 1 : nDimTrial - eval::range_entity_type::nDim;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type contextTrial = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTrial_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                    mpl::int_<expression_type::context|vm::POINT|vm::JACOBIAN> >::type::value;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc_formTrial_type;
    typedef typename gm1_formTrial_type::template Context<geoelement_formTrial_type,gmTrialRangeRelation> gmc1_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef std::shared_ptr<gmc1_formTrial_type> gmc1_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm1_formTrial_type::precompute_type pc1_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef typename gm1_formTrial_type::precompute_ptrtype pc1_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTrial_ptrtype> > map_gmc1_formTrial_type;

    //-----------------------------------------------------//

    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im1_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im1_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;
    using im_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;
    using im1_formtrial_type = im_t<geoelement_formTrial_type, expression_value_type>;

    im_range_type imRange( M_im.order() );
    im1_range_type im1Range( M_im2.order() );
    im_formtest_type imTest( M_im.order() );
    im_formtrial_type imTrial( M_im.order() );

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
    // mortar context
    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const bool has_mortar_trial = FormType::trial_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test && has_mortar_trial)? 3 : ( (has_mortar_test)? 1 : ( (has_mortar_trial)? 2 : 0 ) );
    BOOST_MPL_ASSERT_MSG( mortarTag < 3,TODO_CASE3, (mpl::int_<mortarTag>) );

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type, map_gmc1_formTrial_type> form1_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type,mortarTag> form_mortar_context_type;

    typedef std::shared_ptr<form_context_type> form_context_ptrtype;
    typedef std::shared_ptr<form1_context_type> form1_context_ptrtype;
    typedef std::shared_ptr<form_mortar_context_type> form_mortar_context_ptrtype;

    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype>,fusion::pair<vf::detail::gmc<1>, gmc_expr_ptrtype> > map2_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype>,fusion::pair<vf::detail::gmc<1>, gmc_formTrial_ptrtype> > map2_gmc_formTrial_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype>,fusion::pair<vf::detail::gmc<1>, gmc_formTest_ptrtype> > map2_gmc_formTest_type;
    typedef typename FormType::template Context<map2_gmc_formTest_type, expression_type, face_im_type, map2_gmc_expr_type, map2_gmc_formTrial_type> form2_context_type;
    typedef std::shared_ptr<form2_context_type> form2_context_ptrtype;

    std::vector<gmc_expr_ptrtype> gmcExpr( 2 );
    std::vector<gmc_formTest_ptrtype> gmcFormTest( 2 );
    std::vector<gmc_formTrial_ptrtype> gmcFormTrial( 2 );
    form_context_ptrtype form;
    form2_context_ptrtype form2;
    form_mortar_context_ptrtype formcm;

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto elt_it = lit->template get<1>();
        auto const elt_en = lit->template get<2>();
        DLOG(INFO) << "face/element integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

        auto mrdmt = Feel::vf::detail::manageRelationDifferentMeshType( __form.testSpace(),__form.trialSpace(), *lit );
        if ( !mrdmt.hasFoundEntityToInit() )
            continue;
        auto const& faceInit = mrdmt.entityToInit();

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
                    __geopc[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( faceInit.element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
                    __geopc1[__f][__p] = pc1_expr_ptrtype(  new pc1_expr_type( faceInit.element( 0 ).gm1(), im1Range.fpoints(__f, __p.value() ) ) );
                }
            }
            hasInitGeoPc=true;
        }


        mrdmt.template initGmc<gmc_context_expr_v,contextTest,contextTrial>( __geopc, this->expression(), imTest, imTrial, gmcExpr[0], gmcFormTest[0], gmcFormTrial[0] );

        for ( ; elt_it != elt_en; ++elt_it )
        {
            auto const& faceCur = boost::unwrap_ref( *elt_it );

            auto eltsTrialTestRelated = mrdmt.eltsRelatedToRange( faceCur );
            if ( eltsTrialTestRelated.empty() )
                continue;

            if ( eltsTrialTestRelated.size() == 1 ) // integrate on boundary faces (one side)
            {
                mrdmt.template updateGmc<gmc_context_expr_v,contextTest,contextTrial>( faceCur, eltsTrialTestRelated[0], gmcExpr[0], gmcFormTest[0], gmcFormTrial[0] );

                uint16_type faceIdInElt = (std::get<0>( eltsTrialTestRelated[0] ) == 0)? faceCur.idInElement0():faceCur.idInElement1();
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr[0] ) );
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest[0] ) );
                map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial[0] ) );
                auto const& eltTest = unwrap_ref( std::get<1>( eltsTrialTestRelated[0] ) );
                auto const& eltTrial = unwrap_ref( std::get<2>( eltsTrialTestRelated[0] ) );
                index_type test_elt_0 = eltTest.id();
                index_type trial_elt_0 = eltTrial.id();

                bool useMortarAssembly = false;
                if constexpr ( mortarTag > 0 )
                {
                    useMortarAssembly = ( has_mortar_test && eltTest.isOnBoundary() ) || ( has_mortar_trial && eltTrial.isOnBoundary() );
                }

                if ( !useMortarAssembly )
                {
                    if ( !form )
                    {
                        form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                            this->expression(), face_ims[faceIdInElt], imRange,imTest,imTrial ) );
                    }
                    form->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr, face_ims[faceIdInElt] );
                    form->integrate();
                    form->assemble( std::make_pair(test_elt_0, trial_elt_0) );
                }
                else
                {
                    if ( !formcm )
                    {
                        formcm = form_mortar_context_ptrtype( new form_mortar_context_type( __form, mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,
                                                                                     this->expression(), face_ims[faceIdInElt], imRange,imTest,imTrial ) );
                    }
                    formcm->update( mapgmcFormTest,mapgmcFormTrial,mapgmcExpr, face_ims[faceIdInElt] );
                    formcm->integrate();
                    formcm->assemble( std::make_pair(test_elt_0, trial_elt_0) );
                }

            }
            else if ( eltsTrialTestRelated.size() == 2 ) // integrate on internal faces (on both sides)
            {
                if ( !gmcExpr[1] )
                    mrdmt.template initGmc<gmc_context_expr_v,contextTest,contextTrial>( __geopc, this->expression(), imTest, imTrial, gmcExpr[1], gmcFormTest[1], gmcFormTrial[1] );

                mrdmt.template updateGmc<gmc_context_expr_v,contextTest,contextTrial>( faceCur, eltsTrialTestRelated[0], gmcExpr[0], gmcFormTest[0], gmcFormTrial[0] );
                mrdmt.template updateGmc<gmc_context_expr_v,contextTest,contextTrial>( faceCur, eltsTrialTestRelated[1], gmcExpr[1], gmcFormTest[1], gmcFormTrial[1] );

                uint16_type faceIdInElt = (std::get<0>( eltsTrialTestRelated[0] ) == 0)? faceCur.idInElement0():faceCur.idInElement1();
                auto mapgmctest2 = mapgmc( gmcFormTest[0], gmcFormTest[1] );
                auto mapgmctrial2 = mapgmc( gmcFormTrial[0], gmcFormTrial[1] );
                auto mapgmcexpr2 = mapgmc( gmcExpr[0], gmcExpr[1] );
                index_type test_elt_0 = unwrap_ref( std::get<1>( eltsTrialTestRelated[0] ) ).id();
                index_type trial_elt_0 = unwrap_ref( std::get<2>( eltsTrialTestRelated[0] ) ).id();
                index_type test_elt_1 = unwrap_ref( std::get<1>( eltsTrialTestRelated[1] ) ).id();
                index_type trial_elt_1 = unwrap_ref( std::get<2>( eltsTrialTestRelated[1] ) ).id();

                if ( !form2 )
                {
                    form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmctest2, mapgmctrial2, mapgmcexpr2,
                                                                           this->expression(), face_ims[faceIdInElt],
                                                                           imRange, imTest, imTrial,
                                                                           mpl::int_<2>() ) );
                }

                form2->update( mapgmctest2,mapgmctrial2,mapgmcexpr2, face_ims[faceIdInElt] );
                form2->integrate();
                form2->assemble( std::make_pair(test_elt_0, trial_elt_0),
                                 std::make_pair(test_elt_1, trial_elt_1) );

            }
        }
    }

}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleWithRelationDifferentMeshType(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
{
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename eval::gm1_type gm1_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_expr_type::template Context<typename eval::element_type,1> gmc_expr_type;
    typedef typename gm1_expr_type::template Context<typename eval::element_type,1> gmc1_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef std::shared_ptr<gmc1_expr_type> gmc1_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm1_expr_type::precompute_type pc1_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef typename gm1_expr_type::precompute_ptrtype pc1_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_expr_ptrtype> > map_gmc1_expr_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype>,fusion::pair<vf::detail::gmc<1>, gmc_expr_ptrtype> > map2_gmc_expr_type;

    // typedef on linearform :
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    // test space
    typedef typename FormType::gm_type gm_formTest_type;
    typedef typename FormType::gm1_type gm1_formTest_type;
    typedef typename FormType::mesh_test_element_type geoelement_formTest_type;
    static const size_type contextTest = mpl::if_< boost::is_same< mpl::int_< gm_expr_type::nDim >,mpl::int_< gm_formTest_type::nDim > >,
                                                   mpl::int_<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT>,
                                                   mpl::int_<expression_type::context|vm::POINT|vm::JACOBIAN> >::type::value;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type> gmc_formTest_type;
    typedef typename gm1_formTest_type::template Context<geoelement_formTest_type> gmc1_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef std::shared_ptr<gmc1_formTest_type> gmc1_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm1_formTest_type::precompute_type pc1_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef typename gm1_formTest_type::precompute_ptrtype pc1_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_formTest_ptrtype> > map_gmc1_formTest_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype>,fusion::pair<vf::detail::gmc<1>, gmc_formTest_ptrtype> > map2_gmc_formTest_type;

    //-----------------------------------------------------//

    static const uint16_type nDimTest = gm_formTest_type::nDim;
    static const uint16_type nDimRange = gm_expr_type::nDim;
    static const uint16_type gmTestRangeRelation = ( nDimTest > nDimRange )? nDimTest-nDimRange : nDimRange-nDimTest;

    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im1_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    using im_formtest_type = im_t<geoelement_formTest_type, expression_value_type>;


    im_range_type imRange( M_im.order() );
    im1_range_type im1Range( M_im2.order() );
    im_formtest_type imTest( M_im.order() );

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
    // mortar context
    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test)? 1 : 0;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTest_type,mortarTag> form_mortar_context_type;
    typedef typename FormType::template Context<map_gmc1_formTest_type, expression_type, face_im1_type, map_gmc1_expr_type> form1_context_type;
    typedef std::shared_ptr<form_context_type> form_context_ptrtype;
    typedef std::shared_ptr<form_mortar_context_type> form_mortar_context_ptrtype;
    typedef std::shared_ptr<form1_context_type> form1_context_ptrtype;

    typedef typename FormType::template Context<map2_gmc_formTest_type, expression_type, face_im_type, map2_gmc_expr_type> form2_context_type;
    typedef typename FormType::template Context<map2_gmc_formTest_type, expression_type, face_im_type, map2_gmc_expr_type, map2_gmc_expr_type, mortarTag> form2_mortar_context_type;
    typedef std::shared_ptr<form2_context_type> form2_context_ptrtype;
    typedef std::shared_ptr<form2_mortar_context_type> form2_mortar_context_ptrtype;

    std::vector<gmc_expr_ptrtype> gmcExpr( 2 );
    std::vector<gmc_formTest_ptrtype> gmcFormTest( 2 );
    form_context_ptrtype form;
    form_mortar_context_ptrtype formMortar;

    form2_context_ptrtype form2;

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto elt_it = lit->template get<1>();
        auto const elt_en = lit->template get<2>();
        DLOG(INFO) << "integrating over " << std::distance( elt_it, elt_en )  << " elements\n";

        auto mrdmt = Feel::vf::detail::manageRelationDifferentMeshType( __form.testSpace(), *lit );
        if ( !mrdmt.hasFoundEntityToInit() )
            continue;
        auto const& faceInit = mrdmt.entityToInit();

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
                    __geopc[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( faceInit.element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
                    __geopc1[__f][__p] = pc1_expr_ptrtype(  new pc1_expr_type( faceInit.element( 0 ).gm1(), im1Range.fpoints(__f, __p.value() ) ) );
                }
            }
            hasInitGeoPc=true;
        }

        mrdmt.template initGmc<gmc_context_expr_v,contextTest>( __geopc, this->expression(), imTest,gmcExpr[0], gmcFormTest[0] );

        for ( ; elt_it != elt_en; ++elt_it )
        {

            auto const& faceCur = boost::unwrap_ref( *elt_it );

            auto eltsTestRelated = mrdmt.eltsRelatedToRange( faceCur );
            if ( eltsTestRelated.empty() )
                continue;

            if ( eltsTestRelated.size() == 1 ) // integrate on boundary faces (one side)
            {
                mrdmt.template updateGmc<gmc_context_expr_v,contextTest>( faceCur, eltsTestRelated[0], gmcExpr[0], gmcFormTest[0] );
                uint16_type faceIdInElt = (std::get<0>( eltsTestRelated[0] ) == 0)? faceCur.idInElement0():faceCur.idInElement1();
                map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr[0] ) );
                map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest[0] ) );
                index_type test_elt_0 = unwrap_ref( std::get<1>( eltsTestRelated[0] ) ).id();

                bool useMortarAssembly = has_mortar_test && unwrap_ref( std::get<1>( eltsTestRelated[0] ) ).isOnBoundary();
                if ( useMortarAssembly )
                {
                    if ( !formMortar )
                        formMortar = form_mortar_context_ptrtype( new form_mortar_context_type( __form, mapgmcFormTest, mapgmcExpr,
                                                                                                this->expression(), face_ims[faceIdInElt], imRange,imTest ) );
                    formMortar->update( mapgmcFormTest,mapgmcFormTest,mapgmcExpr, face_ims[faceIdInElt] );
                    formMortar->integrate();
                    formMortar->assemble( test_elt_0 );
                }
                else
                {
                    if ( !form )
                        form = form_context_ptrtype( new form_context_type( __form, mapgmcFormTest, mapgmcExpr,
                                                                            this->expression(), face_ims[faceIdInElt], imRange,imTest ) );

                    form->update( mapgmcFormTest,mapgmcFormTest,mapgmcExpr, face_ims[faceIdInElt] );
                    form->integrate();
                    form->assemble( test_elt_0 );
                }
            }
            else if ( eltsTestRelated.size() == 2 ) // integrate on internal faces (on both sides)
            {
                if ( !gmcExpr[1] )
                    mrdmt.template initGmc<gmc_context_expr_v,contextTest>( __geopc, this->expression(), imTest, gmcExpr[1], gmcFormTest[1] );

                mrdmt.template updateGmc<gmc_context_expr_v,contextTest>( faceCur, eltsTestRelated[0], gmcExpr[0], gmcFormTest[0] );
                mrdmt.template updateGmc<gmc_context_expr_v,contextTest>( faceCur, eltsTestRelated[1], gmcExpr[1], gmcFormTest[1] );

                uint16_type faceIdInElt = (std::get<0>( eltsTestRelated[0] ) == 0)? faceCur.idInElement0():faceCur.idInElement1();
                auto mapgmctest2 = mapgmc( gmcFormTest[0], gmcFormTest[1] );
                auto mapgmcexpr2 = mapgmc( gmcExpr[0], gmcExpr[1] );
                index_type test_elt_0 = unwrap_ref( std::get<1>( eltsTestRelated[0] ) ).id();
                index_type test_elt_1 = unwrap_ref( std::get<1>( eltsTestRelated[1] ) ).id();

                if ( !form2 )
                {
                    form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmctest2, mapgmcexpr2,
                                                                           this->expression(), face_ims[faceIdInElt],
                                                                           imRange, imTest,
                                                                           mpl::int_<2>() ) );
                }

                form2->update( mapgmctest2, mapgmctest2, mapgmcexpr2, face_ims[faceIdInElt] );
                form2->integrate();
                form2->assemble( test_elt_0, test_elt_1 );
            }
        }

    } // for( auto lit = ... )

}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                 mpl::int_<MESH_FACES> /**/,
                                                                 mpl::false_ ) const
{
    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_expr_type::template Context<typename eval::element_type,1> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    /*typedef typename FormType::gm_type gm_form_type;
      typedef typename FormType::mesh_element_type geoelement_form_type;
      typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
      typedef std::shared_ptr<gmc_form_type> gmc_form_ptrtype;
      typedef typename gm_form_type::precompute_type pc_form_type;
      typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
      typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;*/

    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type> gmc_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type> gmc_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
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

    auto const& faceInit = boost::unwrap_ref( *elt_it );


    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_range_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange( M_im.order() );
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( faceInit.element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );

        }
    }


    uint16_type __face_id_in_elt_0 = faceInit.pos_first();

    gmc_expr_ptrtype gmcExpr = faceInit.element( 0 ).gm()->template context<gmc_context_expr_v>( faceInit.element( 0 ),__geopcExpr,__face_id_in_elt_0, this->expression().dynamicContext() );
    // gmc_expr_ptrtype gmcExpr( new gmc_expr_type( faceInit.element( 0 ).gm(),
    //                                              faceInit.element( 0 ),
    //                                              __geopcExpr,
    //                                              __face_id_in_elt_0, this->expression().dynamicContext() ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(),  __form.testSpace()->fe()->points() ) );
    //gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element( 0 ), geopcFormTest, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_formTest_ptrtype gmcFormTest = __form.gm()->template context<gmc_context_formTest_v>( __form.testSpace()->mesh()->element( 0 ), geopcFormTest, this->expression().dynamicContext() );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );

    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    //gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_formTrial_ptrtype gmcFormTrial =  __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->element( 0 ), geopcFormTrial, this->expression().dynamicContext() );
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
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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

                    gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                    gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );
                    //std::cout << "eltTest.id() " << eltTest.id() << " eltTest.G() " << eltTest.G()
                    //          << " eltTrial.id() "<< eltTrial.id() << " eltTrial.G() " << eltTrial.G() << std::endl;
                    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();
                    DCHECK( __face_id_in_elt_0 != invalid_uint16_type_value ) << "Invalid face id for element " <<  gmcExpr_it->template get<1>()->id();
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
#if 0
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
#else
        CHECK( false ) << "TODO fix compilation";
#endif
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
    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_expr_type::template Context<typename eval::element_type,1> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef vf::detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    /*typedef typename FormType::gm_type gm_form_type;
      typedef typename FormType::mesh_element_type geoelement_form_type;
      typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
      typedef std::shared_ptr<gmc_form_type> gmc_form_ptrtype;
      typedef typename gm_form_type::precompute_type pc_form_type;
      typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
      typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;*/

    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    static const size_type gmc_context_formTest_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTest_type::template Context<geoelement_formTest_type> gmc_formTest_type;
    typedef std::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    static const size_type gmc_context_formTrial_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_formTrial_type::template Context<geoelement_formTrial_type> gmc_formTrial_type;
    typedef std::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    // typedef on formcontext
    typedef typename im_range_type::face_quadrature_type face_im_type;

    // mortar context
    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const bool has_mortar_trial = FormType::trial_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test && has_mortar_trial)? 3 : ( (has_mortar_test)? 1 : ( (has_mortar_trial)? 2 : 0 ) );
    BOOST_MPL_ASSERT_MSG( mortarTag < 3,TODO_CASE_TEST_TRIAL_MORTAR, (mpl::int_<mortarTag>) );

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type, map_gmc_expr_type, map_gmc_formTrial_type, mortarTag> form_mortar_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;
    typedef std::shared_ptr<form_mortar_context_type> form_mortar_context_ptrtype;

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

    auto const& faceInit = boost::unwrap_ref( *elt_it );


    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_range_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange( M_im.order() );
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( faceInit.element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );

        }
    }


    uint16_type __face_id_in_elt_0 = faceInit.pos_first();

    gmc_expr_ptrtype gmcExpr = faceInit.element( 0 ).gm()->template context<gmc_context_expr_v>( faceInit.element( 0 ),__geopcExpr,__face_id_in_elt_0, this->expression().dynamicContext() );
    // gmc_expr_ptrtype gmcExpr( new gmc_expr_type( faceInit.element( 0 ).gm(),
    //                                              faceInit.element( 0 ),
    //                                              __geopcExpr,
    //                                              __face_id_in_elt_0, this->expression().dynamicContext() ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(),  __form.testSpace()->fe()->points() ) );
    //gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->beginElement()->second, geopcFormTest, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_formTest_ptrtype gmcFormTest =  __form.gm()->template context<gmc_context_formTest_v>( __form.testSpace()->mesh()->beginElement()->second, geopcFormTest, this->expression().dynamicContext() );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial = __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrial, this->expression().dynamicContext() );
    //gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrial, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );

    // mortar
    pc_formTest_ptrtype geopcFormTestMortar( new pc_formTest_type( __form.gm(), __form.template testFiniteElement<FormType::test_space_type::is_mortar/*true*/>()->points() ) );
    //gmc_formTest_ptrtype gmcFormTestMortar( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->beginElement()->second, geopcFormTestMortar, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_formTest_ptrtype gmcFormTestMortar = __form.gm()->template context<gmc_context_formTest_v>( __form.testSpace()->mesh()->beginElement()->second, geopcFormTestMortar, this->expression().dynamicContext() );
    map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );
    pc_formTrial_ptrtype geopcFormTrialMortar( new pc_formTrial_type( __form.gmTrial(), __form.template trialFiniteElement<FormType::trial_space_type::is_mortar/*true*/>()->points() ) );
    //gmc_formTrial_ptrtype gmcFormTrialMortar( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrialMortar, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_formTrial_ptrtype gmcFormTrialMortar = __form.gmTrial()->template context<gmc_context_formTrial_v>( __form.trialSpace()->mesh()->beginElement()->second, geopcFormTrialMortar, this->expression().dynamicContext() );
    map_gmc_formTrial_type mapgmcFormTrialMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrialMortar ) );

    //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               imRange ) );
    form_mortar_context_ptrtype formcMortarTest( new form_mortar_context_type( __form,
                                                                               mapgmcFormTestMortar,
                                                                               mapgmcFormTrial,
                                                                               mapgmcExpr,
                                                                               this->expression(),
                                                                               face_ims[__face_id_in_elt_0],
                                                                               imRange ) );
    form_mortar_context_ptrtype formcMortarTrial( new form_mortar_context_type( __form,
                                                                                mapgmcFormTest,
                                                                                mapgmcFormTrialMortar,
                                                                                mapgmcExpr,
                                                                                this->expression(),
                                                                                face_ims[__face_id_in_elt_0],
                                                                                imRange ) );
    //-----------------------------------------------//

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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
                //std::cout << "eltTest.id() " << eltTest.id() << " eltTest.G() " << eltTest.G()
                //          << " eltTrial.id() "<< eltTrial.id() << " eltTrial.G() " << eltTrial.G() << std::endl;

                auto vec_gmcExpr = M_QPL->getUsableDataInFormContext( themapQuad,ptRefTest,ptRefTrial );
                auto gmcExpr_it = vec_gmcExpr.begin();
                auto const gmcExpr_en = vec_gmcExpr.end();
                bool isFirstExperience = true, isFirstExperienceMortarTest = true, isFirstExperienceMortarTrial = true;

                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();

                    if ( mortarTag == 1 && eltTest.isOnBoundary() )
                    {
                        geopcFormTestMortar->update( gmcExpr_it->template get<2>() );
                        geopcFormTrial->update( gmcExpr_it->template get<3>() );
                        gmcFormTestMortar->template update<gmc_context_formTest_v>( eltTest,geopcFormTestMortar );
                        gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );
                        map_gmc_formTest_type mapgmcFormTestMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTestMortar ) );
                        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                        formcMortarTest->updateInCaseOfInterpolate( mapgmcFormTestMortar, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );
                        formcMortarTest->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperienceMortarTest );
                        isFirstExperienceMortarTest = false;
                    }
                    else if ( mortarTag == 2 && eltTrial.isOnBoundary() )
                    {
                        geopcFormTest->update( gmcExpr_it->template get<2>() );
                        geopcFormTrialMortar->update( gmcExpr_it->template get<3>() );
                        gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                        gmcFormTrialMortar->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrialMortar );
                        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                        map_gmc_formTrial_type mapgmcFormTrialMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrialMortar ) );
                        formcMortarTrial->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrialMortar, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );
                        formcMortarTrial->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperienceMortarTrial );
                        isFirstExperienceMortarTrial = false;
                    }
                    else
                    {
                        geopcFormTest->update( gmcExpr_it->template get<2>() );
                        geopcFormTrial->update( gmcExpr_it->template get<3>() );
                        gmcFormTest->template update<gmc_context_formTest_v>( eltTest,geopcFormTest );
                        gmcFormTrial->template update<gmc_context_formTrial_v>( eltTrial,geopcFormTrial );
                        map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTest ) );
                        map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<vf::detail::gmc<0> >( gmcFormTrial ) );
                        formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );
                        formc->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperience );
                        isFirstExperience = false;
                    }
                }

                if ( !isFirstExperience )
                    formc->assembleInCaseOfInterpolate();
                if ( !isFirstExperienceMortarTest )
                    formcMortarTest->assembleInCaseOfInterpolate();
                if ( !isFirstExperienceMortarTrial )
                    formcMortarTrial->assembleInCaseOfInterpolate();
            }

        }
    } // if (!M_QPL->hasPrecompute())
    else
    {
#if 0
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
#else
        CHECK( false ) << "TODO fix compilation";
#endif
    } //else

    delete formc;
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate(vf::detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
{
    using im_range_type = im_t<typename eval::the_element_type, expression_value_type>;
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    static const size_type gmc_context_expr_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_expr_type::template Context< typename eval::element_type,1> gmc_expr_type;
    typedef std::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (test):
    typedef vf::detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    static const size_type gmc_context_form_v = expression_type::context|vm::POINT|vm::JACOBIAN;
    typedef typename gm_form_type::template Context<geoelement_form_type> gmc_form_type;
    typedef std::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;

    // typedef on formcontext
    typedef typename im_range_type::face_quadrature_type face_im_type;

    static const bool has_mortar_test = FormType::test_space_type::is_mortar;
    static const int mortarTag = (has_mortar_test)? 1 : 0;

    typedef typename FormType::template Context<map_gmc_form_type, expression_type, face_im_type,map_gmc_expr_type> form_context_type;
    typedef typename FormType::template Context<map_gmc_form_type, expression_type, face_im_type,map_gmc_expr_type,map_gmc_form_type,mortarTag> form_mortar_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;
    typedef std::shared_ptr<form_mortar_context_type> form_mortar_context_ptrtype;


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

    auto const& faceInit = boost::unwrap_ref( *elt_it );

    //-----------------------------------------------//

    //QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_range_type>::permutation_type permutation_type;
    //typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );

    im_range_type imRange( M_im.order() );
    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( imRange.nFaces() );
    std::vector<face_im_type> face_ims( imRange.nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
    {
        face_ims[__f] = imRange.face( __f );

        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
            __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( faceInit.element( 0 ).gm(), imRange.fpoints(__f, __p.value() ) ) );
        }
    }


    uint16_type __face_id_in_elt_0 = faceInit.pos_first();

    gmc_expr_ptrtype gmcExpr = faceInit.element( 0 ).gm()->template context<gmc_context_expr_v>( faceInit.element( 0 ),__geopcExpr, __face_id_in_elt_0, this->expression().dynamicContext() );
    // gmc_expr_ptrtype gmcExpr( new gmc_expr_type( faceInit.element( 0 ).gm(),
    //                                              faceInit.element( 0 ),
    //                                              __geopcExpr,
    //                                              __face_id_in_elt_0, this->expression().dynamicContext() ) );

    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//
    // test form context
    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), __form.testSpace()->fe()->points() ) );
    //gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->beginElement()->second, geopcForm, invalid_uint16_type_value, this->expression().dynamicContext() ) );
    gmc_form_ptrtype gmcForm = __form.gm()->template context<gmc_context_form_v>( __form.testSpace()->mesh()->beginElement()->second, geopcForm, this->expression().dynamicContext() );
    map_gmc_form_type mapgmcForm( fusion::make_pair<vf::detail::gmc<0> >( gmcForm ) );

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcForm, mapgmcForm,
                                               mapgmcExpr, this->expression(),
                                               face_ims[__face_id_in_elt_0], imRange ) );

    //-----------------------------------------------//
    // test form context : mortar
    form_mortar_context_ptrtype formcMortar;
    pc_form_ptrtype geopcFormMortar;
    gmc_form_ptrtype gmcFormMortar;
    if ( has_mortar_test )
    {
        // mortar
        geopcFormMortar.reset( new pc_form_type( __form.gm(), __form.template testFiniteElement<has_mortar_test>()->points() ) );
        //gmcFormMortar.reset( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->beginElement()->second, geopcFormMortar, invalid_uint16_type_value, this->expression().dynamicContext() ) );
        gmcFormMortar = __form.gm()->template context<gmc_context_form_v>( __form.testSpace()->mesh()->beginElement()->second, geopcFormMortar, this->expression().dynamicContext() );
        map_gmc_form_type mapgmcFormMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormMortar ) );

        formcMortar.reset( new form_mortar_context_type( __form,
                                                         mapgmcFormMortar, mapgmcFormMortar,
                                                         mapgmcExpr, this->expression(),
                                                         face_ims[__face_id_in_elt_0], imRange ) );
    }
    //-----------------------------------------------//

    auto meshTest = __form.testSpace()->mesh();

    if (!M_QPL)
    {
        M_QPL.reset( new QuadPtLocalization<Elements, Im, Expr >( M_elts, M_im ) );
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
            bool isFirstExperience = true, isFirstExperienceMortar = true;
            if ( has_mortar_test && eltTest.isOnBoundary() )
            {
                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcFormMortar->update( gmcExpr_it->template get<2>() );
                    gmcFormMortar->template update<gmc_context_form_v>( eltTest,geopcFormMortar );
                    map_gmc_form_type mapgmcFormMortar( fusion::make_pair<vf::detail::gmc<0> >( gmcFormMortar ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<vf::detail::gmc<0> >( gmcExpr_it->template get<1>() ) );
                    __face_id_in_elt_0 = gmcExpr_it->template get<1>()->faceId();
                    formcMortar->updateInCaseOfInterpolate( mapgmcFormMortar, mapgmcFormMortar, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->template get<0>() );

                    formcMortar->integrateInCaseOfInterpolate( gmcExpr_it->template get<0>(),isFirstExperienceMortar );
                    isFirstExperienceMortar = false;
                }
                formcMortar->assemble();
            }
            else
            {
                for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it )
                {
                    geopcForm->update( gmcExpr_it->template get<2>() );
                    gmcForm->template update<gmc_context_form_v>( eltTest,geopcForm );
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
    }
    else
    {
#if 0
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
#else
        CHECK( false ) << "TODO fix compilation";
#endif
    } // else

    delete formc;
}

namespace detail_integrator
{
template <typename ExprLambdaType, typename ExprXType, typename ExprYType, typename ExprZType>
auto
generateLambdaExpr( ExprLambdaType const& exprLambda, ExprXType const& exprX, ExprYType const& exprY, ExprZType const& exprZ, mpl::int_<1> ) { return exprLambda(exprX); }
template <typename ExprLambdaType, typename ExprXType, typename ExprYType, typename ExprZType>
auto
generateLambdaExpr( ExprLambdaType const& exprLambda, ExprXType const& exprX, ExprYType const& exprY, ExprZType const& exprZ, mpl::int_<2> ) { return exprLambda(vec(exprX,exprY)); }
template <typename ExprLambdaType, typename ExprXType, typename ExprYType, typename ExprZType>
auto
generateLambdaExpr( ExprLambdaType const& exprLambda, ExprXType const& exprX, ExprYType const& exprY, ExprZType const& exprZ, mpl::int_<3> ) { return exprLambda(vec(exprX,exprY,exprZ)); }
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename T, int M,int N>
decltype(auto)
Integrator<Elements, Im, Expr, Im2>::evaluate( std::vector<Eigen::Matrix<T, M,N>> const& v, bool parallel ) const
{
    DLOG(INFO)  << "integrating over "
                << std::distance( this->beginElement(), this->endElement() )  << " elements\n";

    using m_t = Eigen::Matrix<T, eval::shape::M,eval::shape::N>;

    std::vector<m_t> res( v.size() );
    for( auto& e : res ) e = m_t::Zero();

    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        // make sure that we have elements to iterate over (return 0
        // otherwise)
        if ( it == en )
            continue;

        auto const& eltInit = boost::unwrap_ref( *it );
        auto gm = eltInit.gm();
        auto gm1 = eltInit.gm1();

        auto geopc = gm->preCompute( this->im().points() );
        auto geopc1 = gm1->preCompute( this->im2().points() );
        auto const& worldComm = const_cast<MeshBase<>*>( eltInit.mesh() )->worldComm();
        static const size_type thegmc_context_v = context|vm::JACOBIAN;
        auto ctx = gm->template context<thegmc_context_v>( eltInit, geopc );
        auto ctx1 = gm1->template context<thegmc_context_v>( eltInit, geopc1 );
        double x=100,y=101,z=102;
        //auto expr_= expression()(vec(cst_ref(x),cst_ref(y), cst_ref(z)));
        static const int inputDataDim = M;
        auto expr_= detail_integrator::generateLambdaExpr( expression(), cst_ref(x), cst_ref(y), cst_ref(z), mpl::int_<inputDataDim>() );
        auto expr_evaluator = expr_.evaluator( mapgmc(ctx) );
        auto expr_evaluator1 = expr_.evaluator( mapgmc(ctx1) );

        for ( ; it != en; ++it )
        {
            auto const& eltCur = boost::unwrap_ref( *it );
            switch ( M_gt )
            {
            default:
            case  GeomapStrategyType::GEOMAP_HO :
            {
                ctx->template update<thegmc_context_v>( eltCur );
                expr_evaluator.update( mapgmc(ctx) );
                M_im.update( *ctx );

                int i = 0;
                for( auto const& e : v )
                {
                    x=e(0,0);
                    if ( inputDataDim > 1 )
                        y=e(1,0);
                    if ( inputDataDim > 2 )
                        z=e(2,0);
                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        {
                            res[i]( (int)c1,(int)c2 ) += M_im( expr_evaluator, c1, c2 );
                        }
                    ++i;
                }
            }
            break;

            case GeomapStrategyType::GEOMAP_O1:
            {
                ctx1->template update<thegmc_context_v>( eltCur );
                expr_evaluator1.update( mapgmc(ctx1) );
                M_im2.update( *ctx1 );

                int i = 0;
                for( auto const& e : v )
                {
                    x=e(0,0);
                    if ( inputDataDim > 1 )
                        y=e(1,0);
                    if ( inputDataDim > 2 )
                        z=e(2,0);

                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        {
                            res[i]( (int)c1,(int)c2 ) += M_im2( expr_evaluator1, c1, c2 );
                        }
                    ++i;
                }

            }
            break;

            case GeomapStrategyType::GEOMAP_OPT:
            {
                //DDLOG(INFO) << "geomap opt" << "\n";
                if ( eltCur.isOnBoundary() )
                {
                    ctx->template update<thegmc_context_v>( eltCur );
                    expr_evaluator.update( mapgmc(ctx) );
                    M_im.update( *ctx );

                    int i = 0;
                    for( auto const& e : v )
                    {
                        x=e(0,0);
                        if ( inputDataDim > 1 )
                            y=e(1,0);
                        if ( inputDataDim > 2 )
                            z=e(2,0);

                        for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                            for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )

                            {
                                res[i]( (int)c1,(int)c2 ) += M_im( expr_evaluator, c1, c2 );
                            }
                        ++i;
                    }
                }

                else
                {
                    ctx1->template update<thegmc_context_v>( eltCur );
                    expr_evaluator1.update( mapgmc(ctx1) );
                    M_im2.update( *ctx1 );

                    int i = 0;
                    for( auto const& e : v )
                    {
                        x=e(0,0);
                        if ( inputDataDim > 1 )
                            y=e(1,0);
                        if ( inputDataDim > 2 )
                            z=e(2,0);

                        for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                            for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            {
                                res[i]( (int)c1,(int)c2 ) += M_im2( expr_evaluator1, c1, c2 );
                            }
                        ++i;
                    }
                }
            }

            //break;
            }
        }
    }


    if ( !parallel || ( M_wc->localSize() == 1 ) )
        return res;
    else
    {
        std::vector<m_t> resGlob( res.size(), m_t::Zero() );
        mpi::all_reduce( M_wc->localComm(),
                         res,
                         resGlob,
                         [] ( std::vector<m_t> const& x, std::vector<m_t> const& y )
                         {
                             CHECK( x.size() == y.size() ) << "invalid size " << x.size() << " vs " << y.size();
                             std::vector<m_t> r( x.size() );
                             for ( int k=0;k<x.size();++k )
                                 r[k] = x[k] + y[k];
                             return r;
                         } );
        return resGlob;
    }

}




template<typename Elements, typename Im, typename Expr, typename Im2>
template <int iDimDummy,std::enable_if_t< iDimDummy == MESH_ELEMENTS , bool> >
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluateImpl() const
{
    DLOG(INFO)  << "integrating over "
                << std::distance( this->beginElement(), this->endElement() )  << " elements\n";
    tic();

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
            _v.push_back( boost::cref( _it->second ) );

        tbb::blocked_range<decltype( _v.begin() )> r( _v.begin(), _v.end(), M_grainsize );
        context_type thecontext( this->expression(), this->im(), it->second );

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

            toc("integrating over elements", FLAGS_v>1);
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

                toc("integrating over elements", FLAGS_v>1);
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
                    typedef typename boost::unwrap_reference<typename boost::remove_const<const_t>::type>::type the_element_type;
                    typedef the_element_type element_type;
                    typedef typename the_element_type::gm_type gm_type;
                    typedef std::shared_ptr<gm_type> gm_ptrtype;
                    typedef typename eval::gmc_type gmc_type;
                    typedef std::shared_ptr<gmc_type> gmc_ptrtype;

                    typedef typename the_element_type::gm1_type gm1_type;
                    typedef std::shared_ptr<gm1_type> gm1_ptrtype;
                    typedef typename eval::gmc1_type gmc1_type;
                    typedef std::shared_ptr<gmc1_type> gmc1_ptrtype;

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
                        auto const& eltInit = boost::unwrap_ref( *it );

                        //
                        // Precompute some data in the reference element for
                        // geometric mapping and reference finite element
                        //
                        // warning this is not efficient here, we want to use the geometric mapping
                        // from the elements in order to take advantage of the cache if possible
                        // this change hsa been made in order to circumvent a bug which is not yet found
                        //#warning INEFFICIENT CODE HERE : TO DEBUG
                        //gm_ptrtype gm( new gm_type) ;//it->gm();
                        gm_ptrtype gm( eltInit.gm() );
                        //std::cout << "0.5" << std::endl;
                        gm1_ptrtype gm1( new gm1_type ); //it->gm1();
                        //std::cout << "0.6:  " << gm1.use_count() << " " << gm.use_count() << std::endl;
                        //DDLOG(INFO) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";
                        typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( gm,
                                                                                           this->im().points() ) );
                        typename eval::gmpc1_ptrtype __geopc1( new typename eval::gmpc1_type( gm1,
                                                                                              this->im2().points() ) );

                        // possibly high order
                        gmc_ptrtype __c = gm->template context<eval::gmc_context_v>( eltInit, __geopc, this->expression().dynamicContext() );
                        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
                        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
                        eval_expr_type expr( expression(), mapgmc );
                        typedef typename eval_expr_type::shape shape;

                        // order 1
                        gmc1_ptrtype __c1 = gm1->template context<eval::gmc_context_v>( eltInit, __geopc1, this->expression().dynamicContext() );
                        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
                        map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                        typedef typename expression_type::template tensor<map_gmc1_type> eval_expr1_type;
                        eval_expr1_type expr1( expression(), mapgmc1 );

                        //value_type res1 = 0;
                        for ( ; it != en; ++it )
                        {
                            auto const& eltCur = boost::unwrap_ref( *it );
                            switch ( M_gt )
                            {
                            default:
                            case  GeomapStrategyType::GEOMAP_HO :
                            {
                                __c->template update<eval::gmc_context_v>( eltCur );
                                expr.update( mapgmc );
                                M_im.update( *__c );

                                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                                    {
                                        res( c1,c2 ) += M_im( expr, c1, c2 );
                                    }
                            }
                            break;

                            case GeomapStrategyType::GEOMAP_O1:
                            {
                                //DDLOG(INFO) << "geomap o1" << "\n";
                                __c1->template update<eval::gmc_context_v>( eltCur );
                                expr1.update( mapgmc1 );
                                M_im2.update( *__c1 );

                                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                                    {
                                        res( c1,c2 ) += M_im2( expr1, c1, c2 );
                                    }
                            }
                            break;

                            case GeomapStrategyType::GEOMAP_OPT:
                            {
                                //DDLOG(INFO) << "geomap opt" << "\n";
                                if ( eltCur.isOnBoundary() )
                                {
                                    //DDLOG(INFO) << "boundary element using ho" << "\n";
                                    __c->template update<eval::gmc_context_v>( eltCur );
                                    expr.update( mapgmc );
                                    M_im.update( *__c );

                                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                                        {
                                            res( c1,c2 ) += M_im( expr, c1, c2 );
                                        }
                                }

                                else
                                {
                                    //DDLOG(INFO) << "interior element using order 1" << "\n";
                                    __c1->template update<eval::gmc_context_v>( eltCur );
                                    expr1.update( mapgmc1 );
                                    M_im2.update( *__c1 );

                                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                                        {
                                            res( c1,c2 ) += M_im2( expr1, c1, c2 );
                                        }
                                }
                            }
                            break;
                            }
                        }
                    }
#if defined(FEELPP_HAS_HARTS)
                    perf_mng.stop("total") ;
                    auto const& worldComm = const_cast<MeshBase<>*>( eltInit.mesh() )->worldComm();
                    std::cout << Environment::worldComm().rank() <<  " Total: " << perf_mng.getValueInSeconds("total") << std::endl;
#endif

                    toc("integrating over elements", FLAGS_v>1);
                    return res;
                }
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template <int iDimDummy,std::enable_if_t< iDimDummy == MESH_FACES /*|| ( iDimDummy == MESH_EDGES && Integrator<Elements, Im, Expr, Im2>::eval::gm_type::nDim == 2)*/ , bool> >
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluateImpl() const
{
    DLOG(INFO) << "integrating over "
               << std::distance( this->beginElement(), this->endElement() )  << "faces\n";
    tic();
    //
    // some typedefs
    //
    typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
    typedef typename boost::unwrap_reference<typename boost::remove_const<const_t>::type>::type the_face_element_type;
    typedef typename the_face_element_type::super2::entity_type the_element_type;
    typedef typename the_element_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    static const size_type gmc_context_face_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
    typedef typename gm_type::template Context<the_element_type,1> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;

    typename eval::matrix_type res( eval::matrix_type::Zero() );
    //typename eval::matrix_type res0( eval::matrix_type::Zero() );
    //typename eval::matrix_type res1( eval::matrix_type::Zero() );


    element_iterator face_it, face_en;
    bool findFaceForInit = false;
    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len && !findFaceForInit; ++lit )
    {
        face_it = lit->template get<1>();
        face_en = lit->template get<2>();

        // check that we have elements to iterate over
        if ( face_it == face_en )
            continue;

        auto const& faceInit = boost::unwrap_ref( *face_it );
        if ( !faceInit.isConnectedTo0() )
            continue;

        findFaceForInit = true;
    }
    if ( !findFaceForInit )
        return res;

    auto const& faceInit = boost::unwrap_ref( *face_it );


    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename eval::gmpc_type pc_type;
    typedef typename eval::gmpc_ptrtype pc_ptrtype;
    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( im().nFaces() );
    typedef typename im_type::face_quadrature_type face_im_type;

    CHECK( faceInit.isConnectedTo0() ) << "invalid face with id=" << faceInit.id();
    CHECK( faceInit.element(0).gm() ) << "invalid geometric transformation associated to face id="
                                      <<  faceInit.id() << " and element id " << faceInit.element(0).id();

    gm_ptrtype gm = faceInit.element( 0 ).gm();

    //DDLOG(INFO) << "[integrator] evaluate(faces), gm is cached: " << gm->isCached() << "\n";
    std::vector<im_face_type> __integrators;
    __integrators.reserve( this->im().nFaces() );
    for ( uint16_type __f = 0; __f < this->im().nFaces(); ++__f )
    {
        __integrators.push_back( im( __f ) );
        for ( permutation_type __p( permutation_type::IDENTITY );
              __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = pc_ptrtype(  new pc_type( gm, this->im().fpoints(__f, __p.value() ) ) );
        }
    }

    //uint16_type __face_id_in_elt_0 = faceInit.pos_first();

    // get the geometric mapping associated with element 0
    gmc_ptrtype __c0 = gm->template context<gmc_context_face_v>( faceInit.element( 0 ), __geopc, /*__face_id_in_elt_0*/faceInit.pos_first(), this->expression().dynamicContext() );

    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

    typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
    typedef std::shared_ptr<eval_expr_type> eval_expr_ptrtype;
    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
    eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );

    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
    typedef typename expression_type::template tensor<map2_gmc_type> eval2_expr_type;
    typedef std::shared_ptr<eval2_expr_type> eval2_expr_ptrtype;
    eval2_expr_ptrtype expr2;

    // true if connected to another element, false otherwise
    //bool isConnectedTo1 = faceInit.isConnectedTo1();

    // get the geometric mapping associated with element 1
    gmc_ptrtype __c1;


    for( auto lit = M_elts.begin(), len = M_elts.end(); lit != len; ++lit )
    {
        auto it = lit->template get<1>();
        auto en = lit->template get<2>();

        //
        // start the real intensive job:
        // -# iterate over all elements to integrate over
        // -# construct the associated geometric mapping with the reference element
        // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
        // -# assemble the local contribution in the global representation of the bilinear form
        //
        for ( ; it != en; ++it )
        {
            auto const& faceCur = boost::unwrap_ref( *it );
            if ( faceCur.isConnectedTo1() )
            {
                if ( faceCur.isGhostFace() )
                    continue;

                DCHECK( !faceCur.isOnBoundary() ) << "face id " << faceCur.id() << " on boundary but connected on both sides";
                uint16_type __face_id_in_elt_0 = faceCur.pos_first();
                uint16_type __face_id_in_elt_1 = faceCur.pos_second();

                if ( !__c1 )
                {
                    __c1 = gm->template context<gmc_context_face_v>( faceCur.element( 1 ), __geopc, __face_id_in_elt_1, this->expression().dynamicContext() );
                    map2_gmc_type mapgmc = Feel::vf::mapgmc(__c0,__c1);
                    expr2 = eval2_expr_ptrtype( new eval2_expr_type( expression(), mapgmc ) );
                }

                __c0->template update<gmc_context_face_v>( faceCur.element( 0 ), __face_id_in_elt_0 );
                bool found_permutation = __c1->template updateFromNeighborMatchingFace<gmc_context_face_v>( faceCur.element( 1 ), __face_id_in_elt_1, __c0 );
                CHECK(found_permutation) << "the permutation of quadrature points were not found\n";

#if 0
                std::cout << "face " << faceCur.id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << "  elt 0 : " << faceCur.element( 0 ).id() << "\n"
                          << "  elt 0 G: " << faceCur.element( 0 ).G() << "\n"
                          << "  node elt 0 0 :" << faceCur.element( 0 ).point( faceCur.element( 0 ).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                          << "  node elt 0 1 :" << faceCur.element( 0 ).point( faceCur.element( 0 ).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                          << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                          << "  real nodes 0: " << __c0->xReal() << "\n";
                std::cout << "face " << faceCur.id() << "\n"
                          << " id in elt = " << __face_id_in_elt_1 << "\n"
                          << " elt 1 : " << faceCur.element( 1 ).id() << "\n"
                          << "  elt 1 G: " << faceCur.element( 1 ).G() << "\n"
                          << "  node elt 1 0 :" << faceCur.element( 1 ).point( faceCur.element( 1 ).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                          << "  node elt 1 1 :" << faceCur.element( 1 ).point( faceCur.element( 1 ).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
                          << " ref nodes 1 :" << __c1->xRefs() << "\n"
                          << " real nodes 1:" << __c1->xReal() << "\n";
#endif

                map2_gmc_type mapgmc = Feel::vf::mapgmc(__c0,__c1);
                expr2->update( mapgmc );
                __integrators[__face_id_in_elt_0].update( *__c0 );

                for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    {
                        res( c1,c2 ) += __integrators[__face_id_in_elt_0]( *expr2, c1, c2 );
                    }
            }

            else
            {
                if ( !faceCur.isConnectedTo0() )
                    continue;

                uint16_type __face_id_in_elt_0 = faceCur.pos_first();
                __c0->template update<gmc_context_face_v>( faceCur.element( 0 ), __face_id_in_elt_0 );
                map_gmc_type mapgmc = Feel::vf::mapgmc(__c0);
                expr->update( mapgmc );
                __integrators[__face_id_in_elt_0].update( *__c0 );

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
    toc("integrating over faces", FLAGS_v>1);
    return res;
}

 template<typename Elements, typename Im, typename Expr, typename Im2>
 template <int iDimDummy,std::enable_if_t< iDimDummy == MESH_POINTS , bool> >
 typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
 Integrator<Elements, Im, Expr, Im2>::evaluateImpl() const
 {
     DLOG(INFO)  << "integrating over "
                 << std::distance( this->beginElement(), this->endElement() )  << " points\n";

     // first loop on the points, then retrieve the elements to which they belong
     // and evaluate the integrand expression and accumulate it


 }

 template<typename Elements, typename Im, typename Expr, typename Im2>
     template<typename P0hType>
     typename P0hType::element_type
     Integrator<Elements, Im, Expr, Im2>::broken( std::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const
 {
     DLOG(INFO) << "integrating over "
                << std::distance( this->beginElement(), this->endElement() )  << " elements\n";

     tic();

     //
     // some typedefs
     //
     using the_element_type = typename eval::the_element_type;

     typedef typename the_element_type::gm_type gm_type;
     typedef std::shared_ptr<gm_type> gm_ptrtype;
     static const size_type gmc_context_elt_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::POINT;

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

         auto const& eltInit = boost::unwrap_ref( *it );

         //
         // Precompute some data in the reference element for
         // geometric mapping and reference finite element
         //
         gm_ptrtype gm = eltInit.gm();
         //DDLOG(INFO) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";

         //it = this->beginElement();
         // wait for all the guys
#ifdef FEELPP_HAS_MPI
         auto const& worldComm = const_cast<MeshBase<>*>( eltInit.mesh() )->worldComm();

#if 0
         if ( worldComm.size() > 1 )
         {
             worldComm.barrier();
         }
#endif
#endif

         auto geopc = gm->preCompute( this->im().points() );
         auto ctx = gm->template context<gmc_context_elt_v>( eltInit, geopc, this->expression().dynamicContext() );
         auto expr_evaluator = this->expression().evaluator( vf::mapgmc(ctx) );

         for ( ; it != en; ++it )
         {
             auto const& eltCur = unwrap_ref( *it );
             ctx->template update<gmc_context_elt_v>( eltCur );
             expr_evaluator.update( vf::mapgmc(ctx) );

             M_im.update( *ctx );

             for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
             {
                 size_type i= P0h->dof()->localToGlobal( eltCur.id(), 0, c1 ).index();
                 double v = M_im( expr_evaluator, c1, 0 );
                 p0.set( i, v );
             }
         }
     }
     toc("integrating [broken] over elements", FLAGS_v>1);

     return p0;
 }
 template<typename Elements, typename Im, typename Expr, typename Im2>
     template<typename P0hType>
     typename P0hType::element_type
     Integrator<Elements, Im, Expr, Im2>::broken( std::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const
 {
     DLOG(INFO) << "integrating over "
                << std::distance( this->beginElement(), this->endElement() )  << "faces\n";
     tic();

     //
     // some typedefs
     //
     typedef typename boost::remove_reference<typename element_iterator::reference>::type const_t;
     typedef typename boost::unwrap_reference<typename boost::remove_const<const_t>::type>::type the_face_element_type;
     typedef typename the_face_element_type::super2::entity_type the_element_type;
     typedef typename the_element_type::gm_type gm_type;
     typedef std::shared_ptr<gm_type> gm_ptrtype;
     static const size_type gmc_context_face_v = expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT;
     typedef typename gm_type::template Context<the_element_type> gmc_type;
     typedef std::shared_ptr<gmc_type> gmc_ptrtype;
     //typedef typename eval_expr_type::value_type value_type;
     //typedef typename Im::value_type value_type;

     //BOOST_MPL_ASSERT_MSG( the_element_type::nDim > 1, INVALID_DIM, (mpl::int_<the_element_type::nDim>, mpl::int_<the_face_element_type::nDim>, mpl::identity<thise_face_element_type>, mpl::identity<the_element_type> ) );;

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
         auto const& faceInit = boost::unwrap_ref( *it );


         gm_ptrtype gm = faceInit.element( 0 ).gm();

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


         uint16_type __face_id_in_elt_0 = faceInit.pos_first();

         // get the geometric mapping associated with element 0
         auto __c0 = gm->template context<gmc_context_face_v>( faceInit.element( 0 ), __geopc, __face_id_in_elt_0, this->expression().dynamicContext() );

         typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

         typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
         typedef std::shared_ptr<eval_expr_type> eval_expr_ptrtype;
         map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
         eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );
         expr->init( im() );

         typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype>, fusion::pair<vf::detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
         typedef typename expression_type::template tensor<map2_gmc_type> eval2_expr_type;
         typedef std::shared_ptr<eval2_expr_type> eval2_expr_ptrtype;
         eval2_expr_ptrtype expr2;

         // true if connected to another element, false otherwise
         bool isConnectedTo1 = faceInit.isConnectedTo1();

         // get the geometric mapping associated with element 1
         gmc_ptrtype __c1;

         //value_type res = 0;
         //value_type res1 = 0;
         if ( isConnectedTo1 )
         {
             uint16_type __face_id_in_elt_1 = faceInit.pos_second();

             __c1 = gm->template context<gmc_context_face_v>( faceInit.element( 1 ), __geopc, __face_id_in_elt_1, this->expression().dynamicContext() );

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

             auto const& faceCur = boost::unwrap_ref( *it );

             if ( faceCur.isConnectedTo1() )
             {
                 FEELPP_ASSERT( faceCur.isOnBoundary() == false   )
                     ( faceCur.id() ).error( "face on boundary but connected on both sides" );
                 uint16_type __face_id_in_elt_0 = faceCur.pos_first();
                 uint16_type __face_id_in_elt_1 = faceCur.pos_second();

                 __c0->template update<gmc_context_face_v>( faceCur.element( 0 ), __face_id_in_elt_0 );
                 __c1->template update<gmc_context_face_v>( faceCur.element( 1 ), __face_id_in_elt_1 );

#if 0
                 std::cout << "face " << faceCur.id() << "\n"
                           << " id in elt = " << __face_id_in_elt_1 << "\n"
                           << "  elt 0 : " << faceCur.element( 0 ).id() << "\n"
                           << "  elt 0 G: " << faceCur.element( 0 ).G() << "\n"
                           << "  node elt 0 0 :" << faceCur.element( 0 ).point( faceCur.element( 0 ).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                           << "  node elt 0 1 :" << faceCur.element( 0 ).point( faceCur.element( 0 ).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                           << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                           << "  real nodes 0: " << __c0->xReal() << "\n";
                 std::cout << "face " << faceCur.id() << "\n"
                           << " id in elt = " << __face_id_in_elt_1 << "\n"
                           << " elt 1 : " << faceCur.element( 1 ).id() << "\n"
                           << "  elt 1 G: " << faceCur.element( 1 ).G() << "\n"
                           << "  node elt 1 0 :" << faceCur.element( 1 ).point( faceCur.element( 1 ).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                           << "  node elt 1 1 :" << faceCur.element( 1 ).point( faceCur.element( 1 ).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
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
                     size_type i0 = P0h->dof()->localToGlobal( faceCur.element( 0 ), 0, c1 ).index();
                     size_type i1 =  P0h->dof()->localToGlobal( faceCur.element( 1 ), 0, c1 ).index();
                     double v = __integrators[__face_id_in_elt_0]( *expr2, c1, 0 );
                     p0.add( i0, v );
                     p0.add( i1, v );
                 }
             }

             else
             {
                 uint16_type __face_id_in_elt_0 = faceCur.pos_first();
                 __c0->template update<gmc_context_face_v>( faceCur.element( 0 ), __face_id_in_elt_0 );
                 map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c0 ) );
                 expr->update( mapgmc, __face_id_in_elt_0 );
                 //expr->update( mapgmc );
                 const gmc_type& gmc = *__c0;

                 __integrators[__face_id_in_elt_0].update( gmc );

                 for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                 {
                     size_type i0 = P0h->dof()->localToGlobal( faceCur.element( 0 ), 0, c1 ).index();
                     double v = __integrators[__face_id_in_elt_0]( *expr, c1, 0 );
                     p0.add( i0, v );
                 }
             } // !isConnectedTo1
         } // for loop on face
     }
     //std::cout << "res=" << res << "\n";
     //std::cout << "res1=" << res1 << "\n";
     toc("integrating [broken] over faces", FLAGS_v>1);
     return p0;
 }
 /// \endcond

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
                     std::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<Elts>::type, Im, ExprT > > quadptloc
                     = std::shared_ptr<QuadPtLocalization<typename Feel::detail::quadptlocrangetype<Elts>::type, Im, ExprT > >() )
 {

     typedef typename Feel::detail::quadptlocrangetype< Elts >::type range_type;
     typedef Integrator<range_type, Im, ExprT, Im2> expr_t;

     typedef typename boost::tuples::template element<1,range_type>::type element_iterator;
     static const uint16_type geoOrder = boost::unwrap_reference<typename element_iterator::value_type>::type::nOrder;
     LOG_IF(WARNING, gt != GeomapStrategyType::GEOMAP_HO && geoOrder == 1 ) << "you use a non standard geomap : ";
     return Expr<expr_t>( expr_t( elts, im, expr, gt, im2, use_tbb, use_harts, grainsize, partitioner, quadptloc ) );
 }


 /// \cond DETAIL
 namespace detail
 {
 template<typename _range_type, typename _expr_type, typename __quad_type, typename __quad1_type>
 struct integrate_im_type
 {
     using expr_order_t = ExpressionOrder<_range_type,_expr_type>;
     //using _value_type = typename _expr_type::value_type;
     using im_type = im_t<typename expr_order_t::the_element_type, typename _expr_type::value_type>;
     using _quad_type = typename mpl::if_<mpl::or_<std::is_integral<__quad_type>,
                                                   std::is_base_of<_QBase,__quad_type>>,
                                          mpl::identity<im_type>,
                                          mpl::identity<std::remove_const_t<__quad_type>> >::type::type;
     using _quad1_type = typename mpl::if_<mpl::or_<std::is_integral<__quad1_type>,
                                                    std::is_base_of<_QBase,__quad1_type>>,
                                           mpl::identity<im_type>, mpl::identity<std::remove_const_t<__quad1_type>> >::type::type;

     template <typename QuadType,typename Quad1Type>
     static
     std::pair<_quad_type,_quad1_type>
     im( QuadType const& thequad, Quad1Type const& thequad1, _expr_type const& expr, std::enable_if_t< std::is_integral<QuadType>::value && std::is_integral<Quad1Type>::value >* = nullptr )
         {
             quad_order_type exprOrder = expr_order_t::value( expr );
             quad_order_type exprOrder_1 = expr_order_t::value_1( expr );
             if ( thequad == quad_order_from_expression && thequad1 == quad_order_from_expression )
                 return std::make_pair( Feel::im<_quad_type>( exprOrder ), Feel::im<_quad1_type>( exprOrder_1 ) );
             else if ( thequad == quad_order_from_expression )
                 return std::make_pair( Feel::im<_quad_type>( exprOrder ), Feel::im<_quad1_type>( thequad1 ) );
             else if ( thequad1 == quad_order_from_expression )
                 return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad ) );
             else
                 return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad1 ) );
         }
     template <typename QuadType,typename Quad1Type>
     static
     std::pair<_quad_type,_quad1_type>
     im( QuadType const& thequad, Quad1Type const& thequad1, _expr_type const& expr, std::enable_if_t< std::is_integral<QuadType>::value && !std::is_integral<Quad1Type>::value >* = nullptr )
         {
             quad_order_type exprOrder = expr_order_t::value( expr );
             if ( thequad == quad_order_from_expression )
                 return std::make_pair( Feel::im<_quad_type>( exprOrder ), Feel::im<_quad1_type>( thequad1 ) );
             else
                 return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad1 ) );
         }
     template <typename QuadType,typename Quad1Type>
     static
     std::pair<_quad_type,_quad1_type>
     im( QuadType const& thequad, Quad1Type const& thequad1, _expr_type const& expr, std::enable_if_t< !std::is_integral<QuadType>::value && std::is_integral<Quad1Type>::value >* = nullptr )
         {
             if ( thequad1 == quad_order_from_expression )
                 return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad ) );
             else
                 return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad1 ) );
         }

     template <typename QuadType,typename Quad1Type>
     static
     std::pair<_quad_type,_quad1_type>
     im( QuadType const& thequad, Quad1Type const& thequad1, _expr_type const& expr, std::enable_if_t< !std::is_integral<QuadType>::value && !std::is_integral<Quad1Type>::value >* = nullptr )
         {
             return std::make_pair( Feel::im<_quad_type>( thequad ), Feel::im<_quad1_type>( thequad1 ) );
         }
 };
#if 0
 template<typename Args>
 struct integrate_type
 {
     typedef clean2_type<Args,tag::expr,Expr<Cst<double> > > _expr_type;
     typedef typename Feel::detail::quadptlocrangetype<clean_type<Args,tag::range>>::type _range_type;
     typedef typename boost::tuples::template element<1, _range_type>::type _element_iterator;
     static const uint16_type geoOrder = boost::unwrap_reference<typename _element_iterator::value_type>::type::nOrder;
     using _element_type = typename boost::unwrap_reference<typename _element_iterator::value_type>::type;
#if 0
     //typedef _Q< ExpressionOrder<_range_type,_expr_type>::value > the_quad_type;
     static const uint16_type exprOrder = ExpressionOrder<_range_type,_expr_type>::value;
     static const uint16_type exprOrder_1 = ExpressionOrder<_range_type,_expr_type>::value_1;
     /**
      * @return expression order
      */
     static constexpr uint16_type expressionOrder() { return exprOrder ; }
     /**
      * @return expression order in the context of an order 1 geometry
      */
     static constexpr uint16_type expressionOrderG1() { return exprOrder_1 ; }
#endif
     using expr_order_t = ExpressionOrder<_range_type,_expr_type>;
     //using _value_type = typename _expr_type::value_type;
     using im_default_type = im_t<typename expr_order_t::the_element_type, typename _expr_type::value_type>;

     typedef clean2_type<Args,tag::quad, im_default_type> __quad_type;
     typedef clean2_type<Args,tag::quad1, im_default_type > __quad1_type;
     using _im_type = integrate_im_type<_range_type,_expr_type,__quad_type,__quad1_type>;
     using _quad_type = typename _im_type::_quad_type;
     using _quad1_type = typename _im_type::_quad1_type;

     typedef Expr<Integrator<_range_type, _quad_type, _expr_type, _quad1_type> > expr_type;

     typedef std::shared_ptr<QuadPtLocalization<_range_type,_quad_type,_expr_type > > _quadptloc_ptrtype;

 };
 #endif
 } // detail

 /// \endcond
} // vf


} // feel



#endif /* FEELPP_INTEGRATORS_HPP */
