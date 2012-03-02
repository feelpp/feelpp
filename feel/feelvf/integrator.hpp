/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-20

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-20
 */
#ifndef __Integrators_H
#define __Integrators_H 1
#include <cxxabi.h>
#include <typeinfo>
#include <boost/timer.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelvf/block.hpp>


#include <Eigen/Eigen>


//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/sparse_traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_sparse.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/quadmapped.hpp>
#include <feel/feelvf/formcontextbase.hpp>
#include <feel/feelvf/bilinearform.hpp>
#include <feel/feelvf/linearform.hpp>
#include <feel/feeldiscr/quadptlocalization.hpp>
#if defined( HAVE_GOOGLE_PROFILER_H )
#include <google/profiler.h>
#endif

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
class Integrator
{
public:

    /** @name Constants
     */
    //@{

    static const size_type context = Expr::context;

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
        typedef typename gm_type::template Context<expression_type::context, the_element_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm1_type::template Context<expression_type::context, the_element_type> gmc1_type;
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
        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef typename eval_expr_type::shape shape;

        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
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
#endif
        typedef Eigen::Matrix<expression_value_type, shape::M, shape::N> matrix_type;
        static value_type zero( mpl::bool_<false> ) { return value_type::Zero(); }
        static value_type zero( mpl::bool_<true> ) { return 0; }
        static value_type zero()  { return zero( boost::is_scalar<value_type>() ); }

    };

    typedef typename eval::im_type im_type;
    typedef typename eval::im2_type im2_type;
    typedef typename im_type::face_quadrature_type im_face_type;
    typedef typename im2_type::face_quadrature_type im2_face_type;
    //typedef typename eval::value_type value_type;
    typedef typename eval::matrix_type matrix_type;
    typedef typename eval::matrix_type value_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Integrator( Elements const& elts, Im const& /*__im*/, expression_type const& __expr, GeomapStrategyType gt, Im2 const& /*__im2*/, bool use_tbb, int grainsize, std::string const& partitioner )
        :
        M_eltbegin( elts.template get<1>() ),
        M_eltend( elts.template get<2>() ),
        M_im( ),
        M_im2( ),
        M_expr( __expr ),
        M_gt( gt ),
        M_use_tbb( use_tbb ),
        M_grainsize( grainsize ),
        M_partitioner( partitioner )
    {
        Debug( 5065 ) << "Integrator constructor from expression\n";
    }

    Integrator( Integrator const& __vfi )
        :
        M_eltbegin( __vfi.M_eltbegin ),
        M_eltend( __vfi.M_eltend ),
        M_im( __vfi.M_im ),
        M_im2( __vfi.M_im2 ),
        M_expr( __vfi.M_expr ),
        M_gt( __vfi.M_gt ),
        M_use_tbb( __vfi.M_use_tbb ),
        M_grainsize( __vfi.M_grainsize ),
        M_partitioner( __vfi.M_partitioner )
    {
        Debug( 5065 ) << "Integrator copy constructor\n";
    }

    virtual ~Integrator() {}

    //@}

    /** @name Operator overloads
     */
    //@{


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
    im_type const& im() const { return M_im; }

    /**
     * get the integration method
     *
     *
     * @return the integration method
     */
    im2_type const& im2() const { return M_im2; }

    /**
     * get the integration method on face f
     *
     *
     * @return the integration method on face f
     */
    im_face_type  im( uint16_type f ) const { return M_im.face( f ); }

    /**
     * get the integration method on face f
     *
     *
     * @return the integration method on face f
     */
    im2_face_type  im2( uint16_type f ) const { return M_im2.face( f ); }

    /**
     * get the variational expression
     *
     *
     * @return the variational expression
     */
    expression_type const& expression() const { return M_expr; }


    /**
     * get the geometric mapping integrator type
     *
     * @return the geometric mapping integrator type
     */
    GeomapStrategyType geomapIntegratorType() const { return M_gt; }

    /**
     * iterator that points at the beginning of the container that
     * holds the data that will apply the integration upon
     */
    element_iterator beginElement() const { return M_eltbegin; }

    /**
     * iterator that points at the end of the container that
     * holds the data that will apply the integration upon
     */
    element_iterator endElement() const { return M_eltend; }

    /**
     * tell whether the expression is symmetric or not
     *
     *
     * @return true if symmetric, false otherwise
     */
    bool isSymmetric() const { return M_expr.isSymmetric(); }

    //@}

    /** @name  Mutators
     */
    //@{

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
    //typename expression_type::template tensor<Geo_t>::value_type
    matrix_type
    evaluate() const
    {
        return evaluate( mpl::int_<iDim>() );
    }

    matrix_type
    evaluateAndSum() const
    {
        typename eval::matrix_type loc =  evaluate( mpl::int_<iDim>() );
        typename eval::matrix_type glo( loc );
#if defined( HAVE_MPI )
        if ( M_comm.size() > 1 )
            {
                MPI_Allreduce( loc.data().begin(),
                               glo.data().begin(),
                               loc.size1()*loc.size2(),
                               MPI_DOUBLE,
                               MPI_SUM,
                               M_comm );
            }
#endif // HAVE_MPI
        return glo;
    }

#if defined( HAVE_TBB )
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
        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
        //typedef detail::FormContextBase<map_gmc_type,im_type> fcb_type;
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
                                            fusion::make_pair<detail::gmc<0> >( M_c ),
                                            fusion::make_pair<detail::gmc<0> >( M_c ),
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
                auto  mapgmc = fusion::make_pair<detail::gmc<0> >( M_c );
                for( auto _elt = r.begin(); _elt != r.end(); ++_elt )
                {
                    M_c->update( *_elt );
                    M_formc->update(mapgmc, mapgmc );
                    M_formc->integrate();
                    lock.acquire( m );
                    M_formc->assemble();
                    lock.release();
                }
#else
                for( auto _elt = r.begin(); _elt != r.end(); ++_elt )
                {
                    M_c->update( *_elt );
                    M_formc->update( fusion::make_pair<detail::gmc<0> >( M_c ),
                                     fusion::make_pair<detail::gmc<0> >( M_c ) );
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
        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
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
            M_expr( _expr, map_gmc_type( fusion::make_pair<detail::gmc<0> >( M_c ) ) ),
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
                for( auto _elt = r.begin(); _elt != r.end(); ++_elt )
                {
                    M_c->update( *_elt );
                    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( M_c ) );

                    M_expr.update( mapgmc );
                    M_im.update( *M_c );

#if 1
                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        {
                            M_ret(c1,c2) += M_im( M_expr, c1, c2 );
                        }
#endif
                }
#else
#if 0
                    for(int i = 0; i < 10000; ++i )
                        for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                            for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            {
                                M_ret(c1,c2) += i*i;//M_im( M_expr, c1, c2 );
                            }
#endif
#endif
            }
        void join( ContextEvaluate const& other )
            {
                M_ret += other.M_ret;
            }

        value_type result() const { return M_ret; }

        gm_ptrtype M_gm;
        gmpc_ptrtype M_geopc;
        gmc_ptrtype M_c;
        eval_expr_type M_expr;
        im_type M_im;
        value_type M_ret;
    };
#endif // HAVE_TBB
    //@}

private:
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/ ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<false> /**/ ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<true> /**/ ) const;
    template<typename FormType>
    void assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<false> /**/ ) const;

    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate( detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    template<typename FE1,typename FE2,typename ElemContType>
    void assembleInCaseOfInterpolate( detail::BilinearForm<FE1,FE2,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleInCaseOfInterpolate( detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const;
    template<typename FE,typename VectorType,typename ElemContType>
    void assembleInCaseOfInterpolate( detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const;


    template<typename P0hType>
    typename P0hType::element_type  broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const;
    template<typename P0hType>
    typename P0hType::element_type  broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const;

    typename eval::matrix_type evaluate( mpl::int_<MESH_ELEMENTS> ) const;
    typename eval::matrix_type evaluate( mpl::int_<MESH_FACES> ) const;

private:

    mpi::communicator M_comm;
    element_iterator M_eltbegin;
    element_iterator M_eltend;
    mutable im_type M_im;
    mutable im2_type M_im2;
    expression_type const&  M_expr;
    GeomapStrategyType M_gt;
    bool M_use_tbb;
    int M_grainsize;
    std::string M_partitioner;
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

    element_iterator it = this->beginElement();
    element_iterator en = this->endElement();

    if ( it == en )
        return;

    if ( dynamic_cast<void*>(const_cast<MeshBase*>(it->mesh())) == dynamic_cast<void*>(__u->mesh().get()) &&
         dynamic_cast<void*>(const_cast<MeshBase*>(it->mesh())) == dynamic_cast<void*>(__v->mesh().get()) )
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<same_mesh_type::value>() );
    else
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<false>() );

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

    element_iterator it = this->beginElement();
    element_iterator en = this->endElement();

    if ( it == en )
        return;

    if ( dynamic_cast<void*>(const_cast<MeshBase*>(it->mesh())) == dynamic_cast<void*>(__v->mesh().get()) )
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<same_mesh_type::value>() );
    else
        assemble( __form, mpl::int_<iDim>(), mpl::bool_<false>() );

    //assemble( __form, mpl::int_<iDim>(), mpl::bool_<true>() );
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<true> /**/ ) const
{
    Debug( 5065 ) << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << " elements\n";
    boost::timer __timer;

#if defined(HAVE_TBB)
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

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
    typedef typename FormType::template Context<map_gmc_type, expression_type, im_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_type, expression_type, im2_type> form1_context_type;
    //typedef detail::FormContextBase<map_gmc_type,im_type> fcb_type;
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



    element_iterator it = this->beginElement();
    element_iterator en = this->endElement();

    // check that we have elements to iterate over
    if ( it == en )
        return;

    gmc_ptrtype __c( new gmc_type( __form.gm(), *it, __geopc ) );
    gmc1_ptrtype __c1( new gmc1_type( __form.gm1(), *it, __geopc1 ) );


    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
    map_gmc1_type mapgmc1( fusion::make_pair<detail::gmc<0> >( __c1 ) );

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
    double t0 = 0, t1 = 0,t2 = 0,t3 = 0;
    //
    // start the real intensive job:
    // -# iterate over all elements to integrate over
    // -# construct the associated geometric mapping with the reference element
    // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
    // -# assemble the local contribution in the global representation of the bilinear form
    //
    for ( ; it != en; ++it )
    {
        switch( M_gt )
        {
        default:
        case GeomapStrategyType::GEOMAP_HO:
        {
            //ti0.restart();
            __c->update( *it );
            //t0+=ti0.elapsed();
#if 0
            std::cout << "Element: " << it->id() << "\n"
                      << " o - points : " << it->G() << "\n"
                      << " o - quadrature :\n"
                      << "     ref : " << this->im().points() << "\n"
                      << "     real : " << __c->xReal() << "\n";
#endif
            //ti1.restart();
            map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
            formc->update( mapgmc,mapgmc,mapgmc );
            //Debug( 5065 )  << "update gmc : " << ti1.elapsed() << "\n";
            //t1+=ti1.elapsed();

            //ti2.restart();
            formc->integrate();
            //Debug( 5065 )  << "integrate : " << ti2.elapsed() << "\n";
            //t2+=ti2.elapsed();

            //ti3.restart();
            formc->assemble();
            //Debug( 5065 )  << "assemble : " << ti3.elapsed() << "\n";
            //t3+=ti3.elapsed();
        }
        break;
        case GeomapStrategyType::GEOMAP_O1:
        {
            //ti0.restart();
            __c1->update( *it );
            //t0+=ti0.elapsed();
#if 0
            Debug( 5065 ) << "Element: " << it->id() << "\n"
                          << " o - points : " << it->G() << "\n"
                          << " o - quadrature :\n"
                          << "     ref : " << this->im().points() << "\n"
                          << "     real : " << __c->xReal() << "\n";
#endif
            //ti1.restart();
            map_gmc1_type mapgmc1( fusion::make_pair<detail::gmc<0> >( __c1 ) );
            formc1->update( mapgmc1,mapgmc1,mapgmc1 );
            //Debug( 5065 )  << "update gmc : " << ti1.elapsed() << "\n";
            //t1+=ti1.elapsed();

            //ti2.restart();
            formc1->integrate();
            //Debug( 5065 )  << "integrate : " << ti2.elapsed() << "\n";
            //t2+=ti2.elapsed();

            //ti3.restart();
            formc1->assemble();
            //Debug( 5065 )  << "assemble : " << ti3.elapsed() << "\n";
            //t3+=ti3.elapsed();
        }
        break;
        case GeomapStrategyType::GEOMAP_OPT:
        {
            if ( it->isOnBoundary() )
            {
                //ti0.restart();
                __c->update( *it );
                //t0+=ti0.elapsed();
#if 0
                Debug( 5065 ) << "Element: " << it->id() << "\n"
                              << " o - points : " << it->G() << "\n"
                              << " o - quadrature :\n"
                              << "     ref : " << this->im().points() << "\n"
                              << "     real : " << __c->xReal() << "\n";
#endif
                //ti1.restart();
                map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
                formc->update( mapgmc,mapgmc,mapgmc );
                //Debug( 5065 )  << "update gmc : " << ti1.elapsed() << "\n";
                //t1+=ti1.elapsed();

                //ti2.restart();
                formc->integrate();
                //Debug( 5065 )  << "integrate : " << ti2.elapsed() << "\n";
                //t2+=ti2.elapsed();

                //ti3.restart();
                formc->assemble();
                //Debug( 5065 )  << "assemble : " << ti3.elapsed() << "\n";
                //t3+=ti3.elapsed();
            }
            else
            {
                //ti0.restart();
                __c1->update( *it );
                //t0+=ti0.elapsed();
#if 0
                Debug( 5065 ) << "Element: " << it->id() << "\n"
                              << " o - points : " << it->G() << "\n"
                              << " o - quadrature :\n"
                              << "     ref : " << this->im().points() << "\n"
                              << "     real : " << __c->xReal() << "\n";
#endif
                //ti1.restart();
                map_gmc1_type mapgmc1( fusion::make_pair<detail::gmc<0> >( __c1 ) );
                formc1->update( mapgmc1,mapgmc1,mapgmc1 );
                //Debug( 5065 )  << "update gmc : " << ti1.elapsed() << "\n";
                //t1+=ti1.elapsed();

                //ti2.restart();
                formc1->integrate();
                //Debug( 5065 )  << "integrate : " << ti2.elapsed() << "\n";
                //t2+=ti2.elapsed();

                //ti3.restart();
                formc1->assemble();
                //Debug( 5065 )  << "assemble : " << ti3.elapsed() << "\n";
                //t3+=ti3.elapsed();
            }
        }
        break;
        }
        } // end loop on elements

#if 0
    Debug( 5065 ) << "[elements] Overall geometric mapping update time : " << (t0+t1+t2) << " per element:" << (t0+t1+t2)/std::distance( this->beginElement(), this->endElement() ) << "\n";
    Debug( 5065 ) << "[elements] Overall geometric mapping update time : " << t0 << "\n";
    Debug( 5065 ) << "[elements] Overall form update time : " << t1 << "\n";
    Debug( 5065 ) << "[elements] Overall local assembly time : " << t2 << "\n";
    Debug( 5065 ) << "[elements] Overall global assembly time : " << t3 << "\n";
#endif
    Debug( 5065 ) << "integrating over elements done in " << __timer.elapsed() << "s\n";

    delete formc;
    delete formc1;
    }
#if defined( HAVE_TBB )
    else
    {
        element_iterator it = this->beginElement();
        element_iterator en = this->endElement();

        if ( it == en )
            return;

        std::vector<boost::reference_wrapper<const typename eval::element_type> > _v;
        for( auto _it = it; _it != en; ++_it )
            _v.push_back(boost::cref(*_it));
        //tbb::blocked_range<decltype(_v.begin())> r( _v.begin(), _v.end(), M_grainsize );
        tbb::blocked_range<decltype(_v.begin())> r( _v.begin(), _v.end(), std::distance( it, en ) );
        Context<FormType,expression_type, im_type, typename eval::the_element_type> thecontext (__form,
                                                                                                this->expression(),
                                                                                                this->im(),
                                                                                                *it);
        tbb::parallel_for( r,  thecontext);
    }
#endif // HAVE_TBB
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_ELEMENTS> /**/, mpl::bool_<false> /**/ ) const
{
    assembleInCaseOfInterpolate(__form,mpl::int_<MESH_ELEMENTS>());
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate( detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                  mpl::int_<MESH_ELEMENTS> /**/ ) const
{
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    // typedef on form (trial and test):
    typedef detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    // test
    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, im_type,map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    auto elt_it = this->beginElement();
    auto elt_en = this->endElement();

    // check that we have elements to iterate over
    if ( elt_it == elt_en )
        return;

   //-----------------------------------------------//
    pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), this->im().points() ) );
    gmc_expr_ptrtype gmcExpr( new gmc_expr_type(elt_it->gm(),*elt_it, geopcExpr ) );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr ) );
    //-----------------------------------------------//
    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(), __form.testSpace()->fe()->points()/* this->im().points()*/ ) );
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element(0), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<detail::gmc<0> >( gmcFormTest ) );
   //-----------------------------------------------//
    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() /*this->im().points()*/ ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element(0), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<detail::gmc<0> >( gmcFormTrial ) );
   //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               this->im() ) );

   //-----------------------------------------------//

    QuadPtLocalization<Elements, Im, Expr > QPL(this->beginElement(),this->endElement(), this->im() );

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    QPL.update( meshTest,meshTrial );

   //-----------------------------------------------//

    auto res_it = QPL.result().begin();
    auto res_en = QPL.result().end();
    for ( ; res_it != res_en ; ++res_it)
        {
            auto idEltTest = res_it->get<0>();
            auto map = res_it->get<1>();
            auto map_it = map.begin();
            auto map_en = map.end();
            for ( ; map_it != map_en ; ++map_it)
                {
                    auto idEltTrial = map_it->first;
                    auto eltTrial = meshTrial->element( idEltTrial );
                    auto eltTest = meshTest->element( idEltTest );

                    auto ptRefTest = map_it->second.get<1>();
                    auto ptRefTrial = map_it->second.get<2>();
                    auto themapQuad = map_it->second.get<0>();

                    auto vec_gmcExpr = QPL.getUsableDataInFormContext(themapQuad,ptRefTest,ptRefTrial);

                    auto gmcExpr_it = vec_gmcExpr.begin();
                    auto gmcExpr_en = vec_gmcExpr.end();
                    bool isFirstExperience = true;
                    for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it)
                        {

                            geopcFormTest->update(gmcExpr_it->get<2>());
                            geopcFormTrial->update(gmcExpr_it->get<3>());

                            gmcFormTest->update(eltTest,geopcFormTest);
                            gmcFormTrial->update(eltTrial,geopcFormTrial);

                            map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<detail::gmc<0> >( gmcFormTest ) );
                            map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<detail::gmc<0> >( gmcFormTrial ) );
                            map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr_it->get<1>() ) );

                            formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr, gmcExpr_it->get<0>() );

                            formc->integrateInCaseOfInterpolate( gmcExpr_it->get<0>(),isFirstExperience );
                            isFirstExperience = false;
                        }
                    formc->assembleInCaseOfInterpolate();

                }
        }


} // assembleInCaseOfInterpolate

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate( detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_ELEMENTS> /**/ ) const
{

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;

    // typedef on formcontext
    typedef typename FormType::template Context<map_gmc_form_type, expression_type, im_type,map_gmc_expr_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    auto elt_it = this->beginElement();
    auto elt_en = this->endElement();

    // check that we have elements to iterate over
    if ( elt_it == elt_en )
        return;

   //-----------------------------------------------//

    pc_expr_ptrtype geopcExpr( new pc_expr_type( elt_it->gm(), this->im().points() ) );
    gmc_expr_ptrtype gmcExpr( new gmc_expr_type(elt_it->gm(),*elt_it, geopcExpr ) );
    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), this->im().points() ) );
    gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->element(0), geopcForm ) );
    map_gmc_form_type mapgmcForm( fusion::make_pair<detail::gmc<0> >( gmcForm ) );

   //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcForm,
                                               mapgmcForm,
                                               mapgmcExpr,
                                               this->expression(),
                                               this->im() ) );

   //-----------------------------------------------//

    QuadPtLocalization<Elements, Im, Expr > QPL(this->beginElement(),this->endElement(), this->im() );

    auto meshTest = __form.testSpace()->mesh();

    QPL.update( meshTest );

   //-----------------------------------------------//

    auto res_it = QPL.resultLinear().begin();
    auto res_en = QPL.resultLinear().end();
    for ( ; res_it != res_en ; ++res_it)
        {

            auto idEltTest = res_it->get<0>();
            auto eltTest = meshTest->element( idEltTest );

            auto ptRefTest = res_it->get<2>();
            auto themapQuad = res_it->get<1>();

            auto vec_gmcExpr = QPL.getUsableDataInFormContext(themapQuad,ptRefTest);
            auto gmcExpr_it = vec_gmcExpr.begin();
            auto gmcExpr_en = vec_gmcExpr.end();
            bool isFirstExperience = true;
            for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it)
                {
                    geopcForm->update(gmcExpr_it->get<2>());
                    gmcForm->update(eltTest,geopcForm);
                    map_gmc_form_type mapgmcForm( fusion::make_pair<detail::gmc<0> >( gmcForm ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr_it->get<1>() ) );
                    formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,gmcExpr_it->get<0>() );

                    formc->integrateInCaseOfInterpolate( gmcExpr_it->get<0>(),isFirstExperience );
                    isFirstExperience = false;
                }

            formc->assemble();
        }

}



template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<true> /**/ ) const
{
    Debug( 5065 ) << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << " faces\n";
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

    QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;
    typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );


    std::vector<std::map<permutation_type, pc_ptrtype> > __geopc( im().nFaces() );
    std::vector<std::map<permutation_type, pc1_ptrtype> > __geopc1( im().nFaces() );
    typedef typename im_type::face_quadrature_type face_im_type;
    typedef typename im2_type::face_quadrature_type face_im2_type;

    //__typeof__(im(__face_id_in_elt_0 ) ) im_face ( im(__face_id_in_elt_0 ) );
    std::vector<face_im_type> face_ims( im().nFaces() );
    std::vector<face_im2_type> face_ims2( im2().nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            face_ims[__f] = this->im( __f );
            face_ims2[__f] = this->im2( __f );

            for( permutation_type __p( permutation_type::IDENTITY );
                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( __form.gm(), ppts[__f].find(__p)->second ) );
                __geopc1[__f][__p] = pc1_ptrtype(  new pc1_type( __form.gm1(), ppts[__f].find(__p)->second ) );
            }
        }

    element_iterator it = beginElement();
    element_iterator en = endElement();

    // check that we have elements to iterate over
    if ( it == en )
        return;

    uint16_type __face_id_in_elt_0 = it->pos_first();

    // get the geometric mapping associated with element 0
    //Debug( 5065 ) << "element " << it->element(0)  << "face " << __face_id_in_elt_0 << " permutation " << it->element(0).permutation( __face_id_in_elt_0 ) << "\n";
    gm_ptrtype __gm = it->element(0).gm();
    gm1_ptrtype __gm1 = it->element(0).gm1();
    //Debug( 5065 ) << "[integrator] evaluate(faces), gm is cached: " << __gm->isCached() << "\n";
    gmc_ptrtype __c0( new gmc_type( __gm, it->element( 0 ), __geopc, __face_id_in_elt_0 ) );
    gmc1_ptrtype __c01( new gmc1_type( __gm1, it->element( 0 ), __geopc1, __face_id_in_elt_0 ) );



    //
    // the case where the face is connected only to one element
    //
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
    typedef typename FormType::template Context<map_gmc_type, expression_type, face_im_type> form_context_type;
    typedef typename FormType::template Context<map_gmc1_type, expression_type, face_im2_type> form1_context_type;
    typedef boost::shared_ptr<form_context_type> form_context_ptrtype;
    typedef boost::shared_ptr<form1_context_type> form1_context_ptrtype;
    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
    map_gmc1_type mapgmc1( fusion::make_pair<detail::gmc<0> >( __c01 ) );
    form_context_ptrtype form;
    form1_context_ptrtype form1;

    //
    // the case where the face is connected only to two elements
    //
    // get the geometric mapping associated with element 1
    gmc_ptrtype __c1;
    gmc1_ptrtype __c11;

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype>, fusion::pair<detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc1_ptrtype>, fusion::pair<detail::gmc<1>, gmc1_ptrtype> > map21_gmc_type;
    typedef typename FormType::template Context<map2_gmc_type, expression_type, face_im_type> form2_context_type;
    typedef typename FormType::template Context<map21_gmc_type, expression_type, face_im2_type> form21_context_type;
    typedef boost::shared_ptr<form2_context_type> form2_context_ptrtype;
    typedef boost::shared_ptr<form21_context_type> form21_context_ptrtype;
    form2_context_ptrtype form2;
    form21_context_ptrtype form21;

    bool isInitConnectionTo0=false;
    bool isInitConnectionTo1=false;
    // true if connected to another element, false otherwise
    if ( it->isConnectedTo1() )
    {
        uint16_type __face_id_in_elt_1 = it->pos_second();

        __c1 = gmc_ptrtype( new gmc_type( __gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );
        __c11 = gmc1_ptrtype( new gmc1_type( __gm1, it->element( 1 ), __geopc1, __face_id_in_elt_1 ) );
        map2_gmc_type mapgmc2( fusion::make_pair<detail::gmc<0> >( __c0 ),
                               fusion::make_pair<detail::gmc<1> >( __c1 ) );
        map21_gmc_type mapgmc21( fusion::make_pair<detail::gmc<0> >( __c01 ),
                                 fusion::make_pair<detail::gmc<1> >( __c11 ) );

        form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
        form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
        isInitConnectionTo1=true;
    }
    else
    {
        form = form_context_ptrtype( new form_context_type( __form, mapgmc, mapgmc, mapgmc, expression(), face_ims[__face_id_in_elt_0], this->im() ) );
        form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
        isInitConnectionTo0=true;
    }

    boost::timer ti0,ti1, ti2, ti3;
    double t0 = 0, t1 = 0,t2 = 0,t3 = 0;
    Debug( 5065 ) << "[Integrator::faces/forms] starting...\n";
    //
    // start the real intensive job:
    // -# iterate over all elements to integrate over
    // -# construct the associated geometric mapping with the reference element
    // -# loop over quadrature loop and assemble the local matrix associated with the bilinear form
    // -# assemble the local contribution in the global representation of the bilinear form
    //
    for ( ; it != en; ++it )
        {
            if ( it->isConnectedTo1())
                {
                    if (!isInitConnectionTo1)
                        {
                            uint16_type __face_id_in_elt_1 = it->pos_second();

                            __c1 = gmc_ptrtype( new gmc_type( __gm, it->element( 1 ), __geopc, __face_id_in_elt_1 ) );
                            __c11 = gmc1_ptrtype( new gmc1_type( __gm1, it->element( 1 ), __geopc1, __face_id_in_elt_1 ) );
                            map2_gmc_type mapgmc2( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                                   fusion::make_pair<detail::gmc<1> >( __c1 ) );
                            map21_gmc_type mapgmc21( fusion::make_pair<detail::gmc<0> >( __c01 ),
                                                     fusion::make_pair<detail::gmc<1> >( __c11 ) );

                            form2 = form2_context_ptrtype( new form2_context_type( __form, mapgmc2, mapgmc2, mapgmc2, expression(), face_ims[__face_id_in_elt_0], this->im(), mpl::int_<2>() ) );
                            form21 = form21_context_ptrtype( new form21_context_type( __form, mapgmc21, mapgmc21, mapgmc21, expression(), face_ims2[__face_id_in_elt_0], this->im2(), mpl::int_<2>() ) );
                            isInitConnectionTo1=true;
                        }
                    switch( M_gt )
                    {
                    default:
                    case GeomapStrategyType::GEOMAP_HO:
                    {
                        FEELPP_ASSERT( it->isOnBoundary() == false  )
                            ( it->id() ).error( "face on boundary but connected on both sides");
                        //ti0.restart();
                        // get the id of the face in each adjacent element
                        uint16_type __face_id_in_elt_0 = it->pos_first();
                        uint16_type __face_id_in_elt_1 = it->pos_second();

                        __c0->update( it->element(0), __face_id_in_elt_0 );
                        __c1->update( it->element(1), __face_id_in_elt_1 );
                        //t0 += ti0.elapsed();

                        //ti1.restart();
                        map2_gmc_type mapgmc2 = map2_gmc_type( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                                               fusion::make_pair<detail::gmc<1> >( __c1 ) );
                        form2->update( mapgmc2, mapgmc2, mapgmc2, face_ims[__face_id_in_elt_0], mpl::int_<2>() );
                        //t1 += ti1.elapsed();

                        //ti2.restart();
                        form2->integrate( );
                        //t2 += ti2.elapsed();

                        //ti3.restart();
                        form2->assemble( it->element(0).id(), it->element(1).id() );
                        //t3 += ti3.elapsed();
                    }
                    break;
                    case GeomapStrategyType::GEOMAP_O1:
                    case GeomapStrategyType::GEOMAP_OPT:
                    {
                        FEELPP_ASSERT( it->isOnBoundary() == false  )
                            ( it->id() ).error( "face on boundary but connected on both sides");
                        //ti0.restart();
                        // get the id of the face in each adjacent element
                        uint16_type __face_id_in_elt_0 = it->pos_first();
                        uint16_type __face_id_in_elt_1 = it->pos_second();

                        __c01->update( it->element(0), __face_id_in_elt_0 );
                        __c11->update( it->element(1), __face_id_in_elt_1 );
                        //t0 += ti0.elapsed();

                        //ti1.restart();
                        map21_gmc_type mapgmc21 = map21_gmc_type( fusion::make_pair<detail::gmc<0> >( __c01 ),
                                                                  fusion::make_pair<detail::gmc<1> >( __c11 ) );
                        form21->update( mapgmc21, mapgmc21, mapgmc21, face_ims2[__face_id_in_elt_0], mpl::int_<2>() );
                        //t1 += ti1.elapsed();

                        //ti2.restart();
                        form21->integrate( );
                        //t2 += ti2.elapsed();

                        //ti3.restart();
                        form21->assemble( it->element(0).id(), it->element(1).id() );
                        //t3 += ti3.elapsed();
                    }
                    break;
                    }
                }
            else
            {
#if 1
                uint16_type __face_id_in_elt_0 = it->pos_first();

                if (!isInitConnectionTo0)
                    {
                        form = form_context_ptrtype( new form_context_type( __form, mapgmc, mapgmc, mapgmc, expression(), face_ims[__face_id_in_elt_0], this->im() ) );
                        form1 = form1_context_ptrtype( new form1_context_type( __form, mapgmc1, mapgmc1, mapgmc1, expression(), face_ims2[__face_id_in_elt_0], this->im2() ) );
                        isInitConnectionTo0=true;
                    }

                //ti0.restart();
                __c0->update( it->element(0),__face_id_in_elt_0 );
                //t0 += ti0.elapsed();

                FEELPP_ASSERT( __face_id_in_elt_0 == __c0->faceId() )
                    ( __face_id_in_elt_0 )
                    ( __c0->faceId() ).warn ( "invalid face id" );

                //ti1.restart();
                map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
                form->update( mapgmc, mapgmc, mapgmc, face_ims[__face_id_in_elt_0] );
                //t1 += ti1.elapsed();

                //ti2.restart();
                form->integrate( );
                //t2 += ti2.elapsed();

                //ti3.restart();
                form->assemble( it->element(0).id() );
                //t3 += ti3.elapsed();
#else
                map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
                form->update( mapgmc, mapgmc, mapgmc, face_ims[__face_id_in_elt_0] );
                form->integrate( );
#endif
            } // end loop on elements

        }
#if 0
    Debug( 5065 ) << "[faces] Overall integration time : " << (t0+t1+t2+t3) << " per element:" << (t0+t1+t2+t3)/std::distance( this->beginElement(), this->endElement() ) << "for " << std::distance( this->beginElement(), this->endElement() ) << "elements\n";
    Debug( 5065 ) << "[faces] Overall geometric mapping update time : " << t0 << "\n";
    Debug( 5065 ) << "[faces] Overall form update time : " << t1 << "\n";
    Debug( 5065 ) << "[faces] Overall local assembly time : " << t2 << "\n";
    Debug( 5065 ) << "[faces] Overall global assembly time : " << t3 << "\n";
#endif
    Debug( 5065 ) << "integrating over faces done in " << __timer.elapsed() << "s\n";
    //std::cout << "integrating over faces done in " << __timer.elapsed() << "s\n";
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FormType>
void
Integrator<Elements, Im, Expr, Im2>::assemble( FormType& __form, mpl::int_<MESH_FACES> /**/, mpl::bool_<false> /**/ ) const
{
    assembleInCaseOfInterpolate(__form,mpl::int_<MESH_FACES>());
}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE1,typename FE2,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate( detail::BilinearForm<FE1,FE2,ElemContType>& __form,
                                                                  mpl::int_<MESH_FACES> /**/ ) const
{

    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (trial and test):
    typedef detail::BilinearForm<FE1,FE2,ElemContType> FormType;
    /*typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;*/

    typedef typename FormType::gm_1_type gm_formTest_type;
    typedef typename FormType::mesh_element_1_type geoelement_formTest_type;
    typedef typename gm_formTest_type::template Context<expression_type::context|vm::POINT,geoelement_formTest_type> gmc_formTest_type;
    typedef boost::shared_ptr<gmc_formTest_type> gmc_formTest_ptrtype;
    typedef typename gm_formTest_type::precompute_type pc_formTest_type;
    typedef typename gm_formTest_type::precompute_ptrtype pc_formTest_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_formTest_ptrtype> > map_gmc_formTest_type;
    // trial
    typedef typename FormType::gm_2_type gm_formTrial_type;
    typedef typename FormType::mesh_element_2_type geoelement_formTrial_type;
    typedef typename gm_formTrial_type::template Context<expression_type::context|vm::POINT,geoelement_formTrial_type> gmc_formTrial_type;
    typedef boost::shared_ptr<gmc_formTrial_type> gmc_formTrial_ptrtype;
    typedef typename gm_formTrial_type::precompute_type pc_formTrial_type;
    typedef typename gm_formTrial_type::precompute_ptrtype pc_formTrial_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_formTrial_ptrtype> > map_gmc_formTrial_type;



    // typedef on formcontext
    typedef typename im_type::face_quadrature_type face_im_type;

    typedef typename FormType::template Context<map_gmc_formTest_type, expression_type, face_im_type,map_gmc_expr_type, map_gmc_formTrial_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    auto elt_it = this->beginElement();
    auto elt_en = this->endElement();

    // check that we have elements to iterate over
    if ( elt_it == elt_en )
        return;

   //-----------------------------------------------//

    QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;
    typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );

    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( im().nFaces() );
    std::vector<face_im_type> face_ims( im().nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            face_ims[__f] = this->im( __f );

            for( permutation_type __p( permutation_type::IDENTITY );
                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element(0).gm(), ppts[__f].find(__p)->second ) );
            }
        }


    uint16_type __face_id_in_elt_0 = elt_it->pos_first();

    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->element(0).gm(),
                                                 elt_it->element( 0 ),
                                                 __geopcExpr,
                                                 __face_id_in_elt_0 ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_formTest_ptrtype geopcFormTest( new pc_formTest_type( __form.gm(),  __form.testSpace()->fe()->points() ) );
    gmc_formTest_ptrtype gmcFormTest( new gmc_formTest_type( __form.gm(), __form.testSpace()->mesh()->element(0), geopcFormTest ) );
    map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<detail::gmc<0> >( gmcFormTest ) );

    pc_formTrial_ptrtype geopcFormTrial( new pc_formTrial_type( __form.gmTrial(), __form.trialSpace()->fe()->points() ) );
    gmc_formTrial_ptrtype gmcFormTrial( new gmc_formTrial_type( __form.gmTrial(), __form.trialSpace()->mesh()->element(0), geopcFormTrial ) );
    map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<detail::gmc<0> >( gmcFormTrial ) );

   //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcFormTest,
                                               mapgmcFormTrial,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               this->im() ) );

   //-----------------------------------------------//

    QuadPtLocalization<Elements, Im, Expr > QPL(this->beginElement(),this->endElement()/*, this->im()*/ );

    auto meshTrial = __form.trialSpace()->mesh();
    auto meshTest = __form.testSpace()->mesh();

    QPL.update( meshTest,meshTrial );

   //-----------------------------------------------//

    auto res_it = QPL.result().begin();
    auto res_en = QPL.result().end();
    for ( ; res_it != res_en ; ++res_it)
        {
            auto idEltTest = res_it->get<0>();
            auto map = res_it->get<1>();
            auto map_it = map.begin();
            auto map_en = map.end();
            for ( ; map_it != map_en ; ++map_it)
                {
                    auto idEltTrial = map_it->first;
                    auto eltTrial = meshTrial->element( idEltTrial );
                    auto eltTest = meshTest->element( idEltTest );

                    auto ptRefTest = map_it->second.get<1>();
                    auto ptRefTrial = map_it->second.get<2>();
                    auto themapQuad = map_it->second.get<0>();

                    auto vec_gmcExpr = QPL.getUsableDataInFormContext(themapQuad,ptRefTest,ptRefTrial);
                    auto gmcExpr_it = vec_gmcExpr.begin();
                    auto gmcExpr_en = vec_gmcExpr.end();
                    bool isFirstExperience = true;
                    for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it)
                        {
                            geopcFormTest->update(gmcExpr_it->get<2>());
                            geopcFormTrial->update(gmcExpr_it->get<3>());

                            gmcFormTest->update(eltTest,geopcFormTest);
                            gmcFormTrial->update(eltTrial,geopcFormTrial);

                            map_gmc_formTest_type mapgmcFormTest( fusion::make_pair<detail::gmc<0> >( gmcFormTest ) );
                            map_gmc_formTrial_type mapgmcFormTrial( fusion::make_pair<detail::gmc<0> >( gmcFormTrial ) );
                            map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr_it->get<1>() ) );
                            __face_id_in_elt_0 = gmcExpr_it->get<1>()->faceId();
                            formc->updateInCaseOfInterpolate( mapgmcFormTest, mapgmcFormTrial, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->get<0>() );

                            formc->integrateInCaseOfInterpolate( gmcExpr_it->get<0>(),isFirstExperience );
                            isFirstExperience = false;
                        }

                    formc->assembleInCaseOfInterpolate();
                }

        }



}

template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename FE,typename VectorType,typename ElemContType>
void
Integrator<Elements, Im, Expr, Im2>::assembleInCaseOfInterpolate( detail::LinearForm<FE,VectorType,ElemContType>& __form, mpl::int_<MESH_FACES> /**/ ) const
{
    // typedef on integral mesh (expr) :
    typedef typename eval::gm_type gm_expr_type;
    typedef typename gm_expr_type::template Context<expression_type::context|vm::JACOBIAN|vm::KB|vm::NORMAL|vm::POINT, typename eval::element_type> gmc_expr_type;
    typedef boost::shared_ptr<gmc_expr_type> gmc_expr_ptrtype;
    typedef typename gm_expr_type::precompute_type pc_expr_type;
    typedef typename gm_expr_type::precompute_ptrtype pc_expr_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_expr_ptrtype> > map_gmc_expr_type;

    //typedef on form (test):
    typedef detail::LinearForm<FE,VectorType,ElemContType> FormType;
    typedef typename FormType::gm_type gm_form_type;
    typedef typename FormType::mesh_test_element_type geoelement_form_type;
    typedef typename gm_form_type::template Context<expression_type::context|vm::POINT,geoelement_form_type> gmc_form_type;
    typedef boost::shared_ptr<gmc_form_type> gmc_form_ptrtype;
    typedef typename gm_form_type::precompute_type pc_form_type;
    typedef typename gm_form_type::precompute_ptrtype pc_form_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_form_ptrtype> > map_gmc_form_type;

    // typedef on formcontext
    typedef typename im_type::face_quadrature_type face_im_type;

    typedef typename FormType::template Context<map_gmc_form_type, expression_type, face_im_type,map_gmc_expr_type> form_context_type;
    typedef form_context_type fcb_type;
    typedef fcb_type* focb_ptrtype;

    //-----------------------------------------------//

    auto elt_it = this->beginElement();
    auto elt_en = this->endElement();

    // check that we have elements to iterate over
    if ( elt_it == elt_en )
        return;

   //-----------------------------------------------//

    QuadMapped<im_type> qm;
    typedef typename QuadMapped<im_type>::permutation_type permutation_type;
    typename QuadMapped<im_type>::permutation_points_type ppts( qm( im() ) );

    std::vector<std::map<permutation_type, pc_expr_ptrtype> > __geopcExpr( im().nFaces() );
    std::vector<face_im_type> face_ims( im().nFaces() );

    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            face_ims[__f] = this->im( __f );

            for( permutation_type __p( permutation_type::IDENTITY );
                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f].find(__p)->second.size2() != 0 ).warn( "invalid quadrature type" );
                __geopcExpr[__f][__p] = pc_expr_ptrtype(  new pc_expr_type( elt_it->element(0).gm(), ppts[__f].find(__p)->second ) );
            }
        }


    uint16_type __face_id_in_elt_0 = elt_it->pos_first();

    gmc_expr_ptrtype gmcExpr( new gmc_expr_type( elt_it->element(0).gm(),
                                                 elt_it->element( 0 ),
                                                 __geopcExpr,
                                                 __face_id_in_elt_0 ) );


    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr ) );

    //-----------------------------------------------//

    pc_form_ptrtype geopcForm( new pc_form_type( __form.gm(), this->im().points() ) );
    gmc_form_ptrtype gmcForm( new gmc_form_type( __form.gm(), __form.testSpace()->mesh()->element(0), geopcForm ) );
    map_gmc_form_type mapgmcForm( fusion::make_pair<detail::gmc<0> >( gmcForm ) );

   //-----------------------------------------------//

    focb_ptrtype formc( new form_context_type( __form,
                                               mapgmcForm,
                                               mapgmcForm,
                                               mapgmcExpr,
                                               this->expression(),
                                               face_ims[__face_id_in_elt_0],
                                               this->im() ) );

   //-----------------------------------------------//

    QuadPtLocalization<Elements, Im, Expr > QPL(this->beginElement(),this->endElement() /*, this->im()*/ );

    auto meshTest = __form.testSpace()->mesh();

    QPL.update( meshTest );

   //-----------------------------------------------//

    auto res_it = QPL.resultLinear().begin();
    auto res_en = QPL.resultLinear().end();
    for ( ; res_it != res_en ; ++res_it)
        {

            auto idEltTest = res_it->get<0>();
            auto eltTest = meshTest->element( idEltTest );

            auto ptRefTest = res_it->get<2>();
            auto themapQuad = res_it->get<1>();
            //geopcForm->update(ptRefTest);
            //gmcForm->update(eltTest,geopcForm);
            //std::cout <<  "\ngmcbeginForm " << gmcForm->xReal();



            auto vec_gmcExpr = QPL.getUsableDataInFormContext(themapQuad,ptRefTest);
            auto gmcExpr_it = vec_gmcExpr.begin();
            auto gmcExpr_en = vec_gmcExpr.end();
            bool isFirstExperience = true;
            //std::cout << "\n start \n";
            for ( ; gmcExpr_it != gmcExpr_en ; ++gmcExpr_it)
                {
                    geopcForm->update(gmcExpr_it->get<2>());
                    gmcForm->update(eltTest,geopcForm);

                    //std::cout << "\ngmcExpr " << gmcExpr_it->get<1>()->xReal()
                    //          << "\ngmcForm " << gmcForm->xReal();

                    map_gmc_form_type mapgmcForm( fusion::make_pair<detail::gmc<0> >( gmcForm ) );
                    map_gmc_expr_type mapgmcExpr( fusion::make_pair<detail::gmc<0> >( gmcExpr_it->get<1>() ) );

                    __face_id_in_elt_0 = gmcExpr_it->get<1>()->faceId();
                    formc->updateInCaseOfInterpolate( mapgmcForm, mapgmcForm, mapgmcExpr,face_ims[__face_id_in_elt_0],gmcExpr_it->get<0>() );

                    formc->integrateInCaseOfInterpolate( gmcExpr_it->get<0>(),isFirstExperience );
                    isFirstExperience = false;
                }

            formc->assemble();
        }


}




template<typename Elements, typename Im, typename Expr, typename Im2>
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluate( mpl::int_<MESH_ELEMENTS> ) const
{
    Debug(5065)  << "integrating over "
                  << std::distance( this->beginElement(), this->endElement() )  << " elements\n";
    boost::timer __timer;

#if defined(HAVE_TBB)
    //std::cout << "Integrator Uses TBB: " << M_use_tbb << "\n";
    if ( !M_use_tbb )
#else
    if ( 1 )
#endif
    {
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

        element_iterator it = this->beginElement();
        element_iterator en = this->endElement();

        // make sure that we have elements to iterate over (return 0
        // otherwise)
        if ( it == en )
            return typename eval::matrix_type(eval::matrix_type::Zero());

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
        gm1_ptrtype gm1( new gm1_type);//it->gm1();
        //std::cout << "0.6:  " << gm1.use_count() << " " << gm.use_count() << std::endl;
        //Debug(5065) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";
        typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( gm,
                                                                           this->im().points() ) );
        //std::cout << "1" << std::endl;
        typename eval::gmpc1_ptrtype __geopc1( new typename eval::gmpc1_type( gm1,
                                                                              this->im().points() ) );

        //std::cout << "2" << std::endl;
        it = this->beginElement();
        // wait for all the guys

#ifdef HAVE_MPI
        if ( M_comm.size() > 1 )
        {
            M_comm.barrier();
        }
#endif


        // possibly high order
        gmc_ptrtype __c( new gmc_type( gm, *it, __geopc ) );
        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
        //std::cout << "3" << std::endl;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        eval_expr_type expr( expression(), mapgmc );
        typedef typename eval_expr_type::shape shape;
        //std::cout << "4" << std::endl;

        // order 1
        gmc1_ptrtype __c1( new gmc1_type( gm1, *it, __geopc1 ) );
        typedef fusion::map<fusion::pair<detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;
        map_gmc1_type mapgmc1( fusion::make_pair<detail::gmc<0> >( __c1 ) );
        //std::cout << "5" << std::endl;
        typedef typename expression_type::template tensor<map_gmc1_type> eval_expr1_type;
        eval_expr1_type expr1( expression(), mapgmc1 );

        //std::cout << "6" << std::endl;
        typename eval::matrix_type res( eval::matrix_type::Zero() );


        //value_type res1 = 0;
        for ( ; it != en; ++it )
        {
            switch( M_gt )
            {
            case  GeomapStrategyType::GEOMAP_HO :
            {
                __c->update( *it );
                map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
                expr.update( mapgmc );
                const gmc_type& gmc = *__c;

                M_im.update( gmc );


                for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        res(c1,c2) += M_im( expr, c1, c2 );
                    }
            }
            break;
            case GeomapStrategyType::GEOMAP_O1:
            {
                //Debug(5065) << "geomap o1" << "\n";
                __c1->update( *it );
                map_gmc1_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c1 ) );
                expr1.update( mapgmc );
                const gmc1_type& gmc = *__c1;

                M_im.update( gmc );


                for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        res(c1,c2) += M_im( expr1, c1, c2 );
                    }
            }
            break;
            case GeomapStrategyType::GEOMAP_OPT:
            {
                //Debug(5065) << "geomap opt" << "\n";
                if ( it->isOnBoundary() )
                {
                    //Debug(5065) << "boundary element using ho" << "\n";
                    __c->update( *it );
                    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
                    expr.update( mapgmc );
                    const gmc_type& gmc = *__c;

                    M_im.update( gmc );


                    for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        {
                            res(c1,c2) += M_im( expr, c1, c2 );
                        }
                }
                else
                {
                    //Debug(5065) << "interior element using order 1" << "\n";
                    __c1->update( *it );
                    map_gmc1_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c1 ) );
                    expr1.update( mapgmc );
                    const gmc1_type& gmc = *__c1;

                    M_im.update( gmc );

                    for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        {
                            res(c1,c2) += M_im( expr1, c1, c2 );
                        }
                }
            }
            //break;
            }
        }
        Debug( 5065 ) << "integrating over elements done in " << __timer.elapsed() << "s\n";
        return res;
    }
    else
    {
#if defined(HAVE_TBB)
        element_iterator it = this->beginElement();
        element_iterator en = this->endElement();
        typedef ContextEvaluate<expression_type, im_type, typename eval::the_element_type> context_type;
        if ( it == en )
            return eval::zero();

        std::vector<boost::reference_wrapper<const typename eval::element_type> > _v;
        for( auto _it = it; _it != en; ++_it )
            _v.push_back(boost::cref(*_it));
        tbb::blocked_range<decltype(_v.begin())> r( _v.begin(), _v.end(), M_grainsize );
        context_type thecontext( this->expression(), this->im(), *it );
        if ( M_partitioner == "auto" )
            tbb::parallel_reduce( r,  thecontext );
        else if ( M_partitioner == "simple" )
            tbb::parallel_reduce( r,  thecontext, tbb::simple_partitioner() );
        //else if ( M_partitioner == "affinity" )
        //tbb::parallel_reduce( r,  thecontext, tbb::affinity_partitioner() );
        return thecontext.result();
#endif // HAVE_TBB
    }

}
template<typename Elements, typename Im, typename Expr, typename Im2>
typename Integrator<Elements, Im, Expr, Im2>::eval::matrix_type
Integrator<Elements, Im, Expr, Im2>::evaluate( mpl::int_<MESH_FACES> ) const
{
    Debug( 5065 ) << "integrating over "
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

    element_iterator it = beginElement();
    element_iterator en = endElement();

    // make sure that we have elements to iterate over (return 0
    // otherwise)
    if ( it == en )
        return typename eval::matrix_type(eval::matrix_type::Zero());

    gm_ptrtype gm = it->element(0).gm();
    //Debug(5065) << "[integrator] evaluate(faces), gm is cached: " << gm->isCached() << "\n";
    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            __integrators.push_back( im(__f) );
            for( permutation_type __p( permutation_type::IDENTITY );
                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f][__p]->size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( gm, ppts[__f].find(__p)->second ) );
            }
        }

    uint16_type __face_id_in_elt_0 = it->pos_first();

    // get the geometric mapping associated with element 0
    gmc_ptrtype __c0( new gmc_type( gm, it->element( 0 ), __geopc, __face_id_in_elt_0 ) );

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

    typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
    typedef boost::shared_ptr<eval_expr_type> eval_expr_ptrtype;
    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
    eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );
    expr->init( im() );

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype>, fusion::pair<detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
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

            map2_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                  fusion::make_pair<detail::gmc<1> >( __c1 ) );

            expr2 = eval2_expr_ptrtype( new eval2_expr_type( expression(), mapgmc ) );
            expr2->init( im() );
        }

    typename eval::matrix_type res( eval::matrix_type::Zero() );
    typename eval::matrix_type res0(eval::matrix_type::Zero() );
    typename eval::matrix_type res1(eval::matrix_type::Zero() );

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
                        ( it->id() ).error( "face on boundary but connected on both sides");
                    uint16_type __face_id_in_elt_0 = it->pos_first();
                    uint16_type __face_id_in_elt_1 = it->pos_second();

                    __c0->update( it->element(0), __face_id_in_elt_0 );
                    __c1->update( it->element(1), __face_id_in_elt_1 );

#if 0
                    std::cout << "face " << it->id() << "\n"
                              << " id in elt = " << __face_id_in_elt_1 << "\n"
                              << "  elt 0 : " << it->element(0).id() << "\n"
                              << "  elt 0 G: " << it->element(0).G() << "\n"
                              << "  node elt 0 0 :" << it->element(0).point( it->element(0).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                              << "  node elt 0 1 :" << it->element(0).point( it->element(0).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                              << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                              << "  real nodes 0: " << __c0->xReal() << "\n";
                    std::cout << "face " << it->id() << "\n"
                              << " id in elt = " << __face_id_in_elt_1 << "\n"
                              << " elt 1 : " << it->element(1).id() << "\n"
                              << "  elt 1 G: " << it->element(1).G() << "\n"
                              << "  node elt 1 0 :" << it->element(1).point( it->element(1).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                              << "  node elt 1 1 :" << it->element(1).point( it->element(1).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
                              << " ref nodes 1 :" << __c1->xRefs() << "\n"
                              << " real nodes 1:" << __c1->xReal() << "\n";
#endif

                    __typeof__(im(__face_id_in_elt_0 ) ) im_face ( im(__face_id_in_elt_0 ) );
                    //std::cout << "pts = " << im_face.points() << "\n";
                    map2_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                          fusion::make_pair<detail::gmc<1> >( __c1 ) );

                    expr2->update( mapgmc, __face_id_in_elt_0 );
                    const gmc_type& gmc = *__c0;

                    __integrators[__face_id_in_elt_0].update( gmc );
                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            {
                                res(c1,c2) += __integrators[__face_id_in_elt_0]( *expr2, c1, c2 );
                            }
                }
            else
                {
                    uint16_type __face_id_in_elt_0 = it->pos_first();
                    __c0->update( it->element(0), __face_id_in_elt_0 );
                    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
                    expr->update( mapgmc, __face_id_in_elt_0 );
                    //expr->update( mapgmc );
                    const gmc_type& gmc = *__c0;

                    __integrators[__face_id_in_elt_0].update( gmc );

                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                        for( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                            {
                                res(c1,c2) += __integrators[__face_id_in_elt_0]( *expr, c1, c2 );
                            }
                } // !isConnectedTo1
        } // for loop on face

    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    Debug( 5065 ) << "integrating over faces done in " << __timer.elapsed() << "s\n";
    return res;
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename P0hType>
typename P0hType::element_type
Integrator<Elements, Im, Expr, Im2>::broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_ELEMENTS> ) const
{
    Debug( 5065 ) << "integrating over "
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

    element_iterator it = this->beginElement();
    element_iterator en = this->endElement();

    auto p0 = P0h->element( "p0" );
    // set to 0 first
    p0.zero();
    // make sure that we have elements to iterate over (return 0
    // otherwise)
    if ( it == en )
        return p0;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    gm_ptrtype gm = it->gm();
    //Debug(5065) << "[integrator] evaluate(elements), gm is cached: " << gm->isCached() << "\n";
    typename eval::gmpc_ptrtype __geopc( new typename eval::gmpc_type( gm,
                                                                       this->im().points() ) );


    it = this->beginElement();
    // wait for all the guys
#ifdef HAVE_MPI
        if ( M_comm.size() > 1 )
            {
                M_comm.barrier();
            }
#endif



    gmc_ptrtype __c( new gmc_type( gm, *it, __geopc ) );
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );

    typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
    eval_expr_type expr( expression(), mapgmc );
    typedef typename eval_expr_type::shape shape;

    //value_type res1 = 0;
    for ( ; it != en; ++it )
        {
            boost::timer ti;
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
            expr.update( mapgmc );
            const gmc_type& gmc = *__c;

            M_im.update( gmc );


            for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
            {
                size_type i;
                boost::tie( i, boost::tuples::ignore, boost::tuples::ignore) = P0h->dof()->localToGlobal( it->id(), 0, c1 );
                double v = M_im( expr, c1, 0 );
                p0.set( i, v );
            }
        }
    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    Debug( 5065 ) << "integrating over elements done in " << __timer.elapsed() << "s\n";

    return p0;
}
template<typename Elements, typename Im, typename Expr, typename Im2>
template<typename P0hType>
typename P0hType::element_type
Integrator<Elements, Im, Expr, Im2>::broken( boost::shared_ptr<P0hType>& P0h, mpl::int_<MESH_FACES> ) const
{
        Debug( 5065 ) << "integrating over "
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

    element_iterator it = beginElement();
    element_iterator en = endElement();

    auto p0 = P0h->element( "p0" );
    // set to 0 first
    p0.zero();
    // make sure that we have elements to iterate over (return 0
    // otherwise)
    if ( it == en )
        return p0;

    gm_ptrtype gm = it->element(0).gm();
    //Debug(5065) << "[integrator] evaluate(faces), gm is cached: " << gm->isCached() << "\n";
    for ( uint16_type __f = 0; __f < im().nFaces(); ++__f )
        {
            __integrators.push_back( im(__f) );
            for( permutation_type __p( permutation_type::IDENTITY );
                 __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                //FEELPP_ASSERT( ppts[__f][__p]->size2() != 0 ).warn( "invalid quadrature type" );
                __geopc[__f][__p] = pc_ptrtype(  new pc_type( gm, ppts[__f].find(__p)->second ) );
            }
        }


    uint16_type __face_id_in_elt_0 = it->pos_first();

    // get the geometric mapping associated with element 0
    gmc_ptrtype __c0( new gmc_type( gm, it->element( 0 ), __geopc, __face_id_in_elt_0 ) );

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

    typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
    typedef boost::shared_ptr<eval_expr_type> eval_expr_ptrtype;
    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
    eval_expr_ptrtype expr( new eval_expr_type( expression(), mapgmc ) );
    expr->init( im() );

    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype>, fusion::pair<detail::gmc<1>, gmc_ptrtype> > map2_gmc_type;
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

            map2_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                  fusion::make_pair<detail::gmc<1> >( __c1 ) );

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
                        ( it->id() ).error( "face on boundary but connected on both sides");
                    uint16_type __face_id_in_elt_0 = it->pos_first();
                    uint16_type __face_id_in_elt_1 = it->pos_second();

                    __c0->update( it->element(0), __face_id_in_elt_0 );
                    __c1->update( it->element(1), __face_id_in_elt_1 );

#if 0
                    std::cout << "face " << it->id() << "\n"
                              << " id in elt = " << __face_id_in_elt_1 << "\n"
                              << "  elt 0 : " << it->element(0).id() << "\n"
                              << "  elt 0 G: " << it->element(0).G() << "\n"
                              << "  node elt 0 0 :" << it->element(0).point( it->element(0).fToP( __face_id_in_elt_0, 0 ) ).node() << "\n"
                              << "  node elt 0 1 :" << it->element(0).point( it->element(0).fToP( __face_id_in_elt_0, 1 ) ).node() << "\n"
                              << "  ref nodes 0 :" << __c0->xRefs() << "\n"
                              << "  real nodes 0: " << __c0->xReal() << "\n";
                    std::cout << "face " << it->id() << "\n"
                              << " id in elt = " << __face_id_in_elt_1 << "\n"
                              << " elt 1 : " << it->element(1).id() << "\n"
                              << "  elt 1 G: " << it->element(1).G() << "\n"
                              << "  node elt 1 0 :" << it->element(1).point( it->element(1).fToP( __face_id_in_elt_1, 1 ) ).node() << "\n"
                              << "  node elt 1 1 :" << it->element(1).point( it->element(1).fToP( __face_id_in_elt_1, 0 ) ).node() << "\n"
                              << " ref nodes 1 :" << __c1->xRefs() << "\n"
                              << " real nodes 1:" << __c1->xReal() << "\n";
#endif

                    __typeof__(im(__face_id_in_elt_0 ) ) im_face ( im(__face_id_in_elt_0 ) );
                    //std::cout << "pts = " << im_face.points() << "\n";
                    map2_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ),
                                          fusion::make_pair<detail::gmc<1> >( __c1 ) );

                    expr2->update( mapgmc, __face_id_in_elt_0 );
                    const gmc_type& gmc = *__c0;

                    __integrators[__face_id_in_elt_0].update( gmc );
                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        size_type i0;
                        boost::tie( i0, boost::tuples::ignore, boost::tuples::ignore) = P0h->dof()->localToGlobal( it->element(0), 0, c1 );
                        size_type i1;
                        boost::tie( i1, boost::tuples::ignore, boost::tuples::ignore) = P0h->dof()->localToGlobal( it->element(1), 0, c1 );
                        double v = __integrators[__face_id_in_elt_0]( *expr2, c1, 0 );
                        p0.add( i0, v );
                        p0.add( i1, v );
                    }
                }
            else
                {
                    uint16_type __face_id_in_elt_0 = it->pos_first();
                    __c0->update( it->element(0), __face_id_in_elt_0 );
                    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c0 ) );
                    expr->update( mapgmc, __face_id_in_elt_0 );
                    //expr->update( mapgmc );
                    const gmc_type& gmc = *__c0;

                    __integrators[__face_id_in_elt_0].update( gmc );

                    for( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        size_type i0;
                        boost::tie( i0, boost::tuples::ignore, boost::tuples::ignore) = P0h->dof()->localToGlobal( it->element(0), 0, c1 );
                        double v = __integrators[__face_id_in_elt_0]( *expr, c1, 0 );
                        p0.add( i0, v );
                    }
                } // !isConnectedTo1
        } // for loop on face

    //std::cout << "res=" << res << "\n";
    //std::cout << "res1=" << res1 << "\n";
    Debug( 5065 ) << "integrating over faces done in " << __timer.elapsed() << "s\n";
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

#if 0
/**
 * integrate an expression \c expr over a set of convexes \c elts
 * using an automatic integration rule .
 */
template<typename Elts, typename ExprT>
Expr<Integrator<Elts, _Q< ExpressionOrder<Elts,ExprT>::value >, ExprT, _Q< ExpressionOrder<Elts,ExprT>::value > > >
integrate( Elts const& elts,
           ExprT const& expr,
           GeomapStrategyType gt = GeomapStrategyType::GEOMAP_HO )
{

    Debug(5065) << "[integrate] order to integrate = " << ExpressionOrder<Elts,ExprT>::value << "\n";
    _Q< ExpressionOrder<Elts,ExprT>::value > quad;
    return integrate_impl( elts, quad, expr, gt, quad );

}
#endif

/**
 * integrate an expression \c expr over a set of convexes \c elts
 * using the integration rule \c im .
 */
template<typename Elts, typename Im, typename ExprT, typename Im2 = Im>
Expr<Integrator<Elts, Im, ExprT, Im2> >
integrate_impl( Elts const& elts,
                Im const& im,
                ExprT const& expr,
                GeomapStrategyType const& gt,
                Im2 const& im2,
                bool use_tbb,
                int grainsize,
                std::string const& partitioner )
{
    typedef Integrator<Elts, Im, ExprT, Im2> expr_t;
    return Expr<expr_t>( expr_t( elts, im, expr, gt, im2, use_tbb, grainsize, partitioner ) );
}


/// \cond DETAIL
namespace detail
{
template<typename TheArgs, typename Tag>
struct clean_type
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<TheArgs, Tag>::type
                >::type
            >::type
        >::type type;
};
template<typename TheArgs, typename Tag, typename Default>
struct clean2_type
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<TheArgs, Tag, Default>::type
                >::type
            >::type
        >::type type;
};

template<typename Args>
struct integrate_type
{
    typedef typename clean_type<Args,tag::expr>::type _expr_type;
    typedef typename clean_type<Args,tag::range>::type _range_type;
    //typedef _Q< ExpressionOrder<_range_type,_expr_type>::value > the_quad_type;
    typedef typename clean2_type<Args,tag::quad, _Q< ExpressionOrder<_range_type,_expr_type>::value > >::type _quad_type;
    typedef typename clean2_type<Args,tag::quad1, _Q< ExpressionOrder<_range_type,_expr_type>::value_1 > >::type _quad1_type;
    typedef Expr<Integrator<_range_type, _quad_type, _expr_type, _quad1_type> > expr_type;

};
}
/// \endcond

/**
 *
 * \brief projection/interpolation of an expresion onto a noal functionspace
 *
 * \arg space the function space to project onto
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 * \arg sum sum the multiple nodal  contributions  if applicable (false by default)
 */
BOOST_PARAMETER_FUNCTION(
    (typename detail::integrate_type<Args>::expr_type), // return type
    integrate,    // 2. function name

    tag,           // 3. namespace of tag types

    (required
     (range, *  )
     (expr,   *)
        ) // 4. one required parameter, and

    (optional
     (quad,   *, typename detail::integrate_type<Args>::_quad_type() )
     (geomap, *, GeomapStrategyType::GEOMAP_OPT )
     (quad1,   *, typename detail::integrate_type<Args>::_quad1_type() )
     (use_tbb,   (bool), false )
     (grainsize,   (int), 100 )
     (partitioner,   *, "auto" )
     (verbose,   (bool), false )
     )
    )
{

    auto ret =  integrate_impl( range, quad, expr, geomap, quad1, use_tbb, grainsize, partitioner );
    if ( verbose )
    {
        std::cout << " -- integrate: size(range) = " << std::distance(ret.expression().beginElement(),
                                                                      ret.expression().endElement() ) << "\n";
        std::cout << " -- integrate: quad = " << ret.expression().im().nPoints() << "\n";
        std::cout << " -- integrate: quad1 = " << ret.expression().im2().nPoints() << "\n";
        //std::cout << " -- integrate: geomap = " << geomap << "\n";
    }

    return ret;
}


} // vf
} // feel
#include <feel/feelvf/integratoron.hpp>


#endif /* __Integrator_H */
