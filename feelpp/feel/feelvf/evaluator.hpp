/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-05-27

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file evaluator.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-05-27
 */
#ifndef __FEELPP_EVALUATORS_H
#define __FEELPP_EVALUATORS_H 1

#include <boost/timer.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
enum EvaluatorType
{
    EVAL_NODAL = 0             /**< Nodal projection */
};
namespace details
{
template<typename ElementT, typename NodeT>
class EvaluatorData: public boost::tuple<ElementT,NodeT>
{
public:
    using super = boost::tuple<ElementT,NodeT>;
    using element_type = ElementT;
    using node_type = NodeT;

    EvaluatorData( element_type e, node_type n ) :
        super( e, n )
        {}

    //! access to element i
    typename element_type::Scalar operator()( int i ) const { return this->template get<0>()(i); }

    //! return the elements
    element_type const& data() const { return this->template get<0>(); }

    node_type const& nodes() const { return this->template get<1>(); }
};
/**
 * \class Evaluator
 * \brief work class to evaluate expressions at sets of points
 *
 * @author Christophe Prud'homme
 */
template<EvaluatorType iDim, typename IteratorRange, typename Pset, typename ExprT>
class Evaluator
{
public:


    /** @name Typedefs
     */
    //@{

    static const size_type context = ExprT::context|vm::POINT;

    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef typename boost::unwrap_reference<typename iterator_type::value_type>::type mesh_element_fromiterator_type;
    typedef typename boost::remove_const< typename boost::remove_reference< mesh_element_fromiterator_type >::type >::type mesh_element_type;
    typedef IteratorRange range_iterator;
    typedef typename mpl::if_<mpl::bool_<mesh_element_type::is_simplex>,
                              mpl::identity<typename Pset::template Apply<mesh_element_type::nRealDim, value_type, Simplex>::type >,
                              mpl::identity<typename Pset::template Apply<mesh_element_type::nRealDim, value_type, Hypercube>::type >
                              >::type::type pointset_type;
    typedef Eigen::Tensor<value_type,4> element_type;
    using node_type = Eigen::Tensor<value_type,3>;
    //typedef Eigen::Matrix<value_type,mesh_element_type::nRealDim,Eigen::Dynamic> node_type;
    using eval_element_type = EvaluatorData<element_type,node_type>;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Evaluator( IteratorRange const& r,
               Pset const& pset,
               expression_type const& __expr,
               GeomapStrategyType geomap_strategy )
        :
        M_range( r ),
        //M_pset( pset.template get<value_type>( typename mesh_element_type::convex_type{} ) ),
        M_pset( pset.template getGeoEntity<value_type,mesh_element_type>() ),
        M_expr( __expr ),
        M_geomap_strategy( geomap_strategy )
    {
        DVLOG(2) << "Evaluator constructor from expression\n";
    }


    Evaluator( Evaluator const& __vfi )
        :
        M_range( __vfi.M_range ),
        M_pset( __vfi.M_pset ),
        M_expr( __vfi.M_expr ),
        M_geomap_strategy( __vfi.M_geomap_strategy )
    {
        DVLOG(2) << "Evaluator copy constructor\n";
    }

    virtual ~Evaluator() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    eval_element_type operator()() const
    {
        return operator()( idim_type() );
    }

    //@}

    /** @name Accessors
     */
    //@{

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

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{
    //! polynomial order
    constexpr uint16_type polynomialOrder() const { return 1; }

    //! expression is polynomial?
    constexpr bool isPolynomial() const { return true; }

    //@}

private:

    eval_element_type operator()( mpl::size_t<MESH_ELEMENTS> ) const;
    eval_element_type operator()( mpl::size_t<MESH_FACES> ) const;

private:

    range_iterator M_range;
    pointset_type M_pset;
    expression_type const&  M_expr;
    GeomapStrategyType M_geomap_strategy;
};

template<EvaluatorType iDim, typename Iterator, typename Pset, typename ExprT>
typename Evaluator<iDim, Iterator, Pset, ExprT>::eval_element_type
Evaluator<iDim, Iterator, Pset, ExprT>::operator()( mpl::size_t<MESH_ELEMENTS> ) const
{
    boost::timer __timer;

    typedef typename mesh_element_type::gm_type gm_type;
    typedef typename gm_type::template Context<mesh_element_type> gm_context_type;
    typedef typename mesh_element_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<mesh_element_type> gm1_context_type;


    typedef std::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef std::shared_ptr<gm1_context_type> gm1_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm1_context_ptrtype> > map_gmc1_type;
    //typedef typename expression_type::template tensor<map_gmc_type,fusion::map<fusion::pair<vf::detail::gmc<0>,std::shared_ptr<fecontext_type> > > > t_expr_type;
    //typedef decltype( basis_type::isomorphism( M_expr ) ) the_expression_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename iso_expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::value_type value_type;

    // we should manipulate the same type of functions on the left and
    // on the right
    //BOOST_STATIC_ASSERT(( boost::is_same<return_value_type, typename functionspace_type::return_value_type>::value ));

    typedef typename t_expr_type::shape shape;


    auto it = boost::get<1>( M_range );
    auto en = boost::get<2>( M_range );

    int npoints = M_pset.points().size2();
    element_type __v( std::distance( it, en ), npoints, shape::M, shape::N );
    node_type __p( std::distance(it,en), npoints, mesh_element_type::nRealDim );
    //node_type __p( std::distance( it, en ), npoints );
    __v.setZero();
    __p.setZero();



    // return if no elements
    if ( it == en )
        return eval_element_type( __v, __p );

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    auto const& initElt = boost::unwrap_ref( *it );
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( initElt.gm(),
                                                                                         M_pset.points() ) );
    typename gm1_type::precompute_ptrtype __geopc1( new typename gm1_type::precompute_type( initElt.gm1(),
                                                                                            M_pset.points() ) );



    gm_context_ptrtype __c = initElt.gm()->template context<context>( initElt,__geopc );
    gm1_context_ptrtype __c1 =  initElt.gm1()->template context<context>( initElt,__geopc1 );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type tensor_expr( M_expr, mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( M_expr, mapgmc1 );

    for ( int e = 0; it!=en ; ++it, ++e )
    {
        auto const& curElt =  boost::unwrap_ref( *it );
        switch ( M_geomap_strategy )
        {
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->template update<context>( curElt );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                {
                    __p(e,p,c1) = __c->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                    {
                        __v(e, p, c1, c2 ) = tensor_expr.evalq( c1, c2, p );
                    }
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->template update<context>( curElt );
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            tensor_expr1.update( mapgmc1 );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                {
                    __p(e,p,c1) = __c1->xReal(p)[c1];
                }
                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                    {
                        __v( e, p, c1, c2 ) = tensor_expr1.evalq( c1, c2, p );
                    }
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_OPT:
        {
            if ( curElt.isOnBoundary() )
            {
                // HO if on boundary
                __c->template update<context>( curElt );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                tensor_expr.update( mapgmc );

                for ( uint16_type p = 0; p < npoints; ++p )
                {
                    for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                    {
                        __p(e,p,c1) = __c->xReal(p)[c1];
                    }

                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    {
                        for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                        {
                            __v( e, p, c1, c2 ) = tensor_expr.evalq( c1, c2, p );
                        }
                    }
                }
            }

            else
            {
                __c1->template update<context>( curElt );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                tensor_expr1.update( mapgmc1 );


                for ( uint16_type p = 0; p < npoints; ++p )
                {
                    for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                    {
                        __p(e,p,c1) = __c1->xReal(p)[c1];
                    }

                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    {
                        for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                        {
                            __v( e, p, c1, c2 ) = tensor_expr1.evalq( c1, c2, p );
                        }
                    }
                }
            }
        }
        break;
        }
    }
    return eval_element_type( __v, __p );
}

template<EvaluatorType iDim, typename Iterator, typename Pset, typename ExprT>
typename Evaluator<iDim, Iterator, Pset, ExprT>::eval_element_type
Evaluator<iDim, Iterator, Pset, ExprT>::operator()( mpl::size_t<MESH_FACES> ) const
{

    boost::timer __timer;

    VLOG(2) << "evaluator(MESH_FACES) " << "\n";
    //
    // a few typedefs
    //

    // mesh element
    typedef typename mesh_element_type::entity_type geoelement_type;
    //typedef typename geoelement_type::face_type face_type;
    typedef mesh_element_type face_type;

    // geometric mapping context
    typedef typename geoelement_type::gm_type gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;

    typedef typename gm_type::template Context<geoelement_type,1> gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm1_type::template Context<geoelement_type,1> gmc1_type;
    typedef std::shared_ptr<gmc1_type> gmc1_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc1_ptrtype> > map_gmc1_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    //typedef decltype( basis_type::isomorphism( M_expr ) ) the_expression_type;
    typedef expression_type the_expression_type;
    typedef typename boost::remove_reference<typename boost::remove_const<the_expression_type>::type >::type iso_expression_type;
    typedef typename iso_expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename iso_expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //

    auto __face_it = boost::get<1>( M_range );
    auto __face_en = boost::get<2>( M_range );

    int npoints = M_pset.fpoints(0,1).size2();
    //element_type __v( M_pset.fpoints(0,1).size2()*std::distance( __face_it, __face_en )*shape::M );
    element_type __v( std::distance( __face_it, __face_en ), npoints, shape::M, shape::N );
    //node_type __p( mesh_element_type::nRealDim, M_pset.fpoints(0,1).size2()*std::distance( __face_it, __face_en ) );
    node_type __p( std::distance( __face_it, __face_en ), npoints, mesh_element_type::nRealDim );

    __v.setZero();
    __p.setZero();
    VLOG(2) << "pset: " << M_pset.fpoints(0,1);
    VLOG(2) << "Checking trivial result...";

    if ( __face_it == __face_en )
        return eval_element_type( __v, __p );

    gm_ptrtype __gm( new gm_type );
    gm1_ptrtype __gm1( new gm1_type );



    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename geoelement_type::permutation_type permutation_type;
    typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
    typedef typename gm_type::precompute_type geopc_type;
    typedef typename gm1_type::precompute_ptrtype geopc1_ptrtype;
    typedef typename gm1_type::precompute_type geopc1_type;

    std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( M_pset.nFaces() );
    std::vector<std::map<permutation_type, geopc1_ptrtype> > __geopc1( M_pset.nFaces() );

    VLOG(2) << "computing geopc...";
    for ( uint16_type __f = 0; __f < M_pset.nFaces(); ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
                __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, M_pset.fpoints(__f, __p.value() ) ) );
            __geopc1[__f][__p] = geopc1_ptrtype(  new geopc1_type( __gm1, M_pset.fpoints(__f, __p.value() ) ) );
            DVLOG(2) << "pset " << __f << " : " << M_pset.fpoints(__f, __p.value() );
            CHECK( __geopc[__f][__p]->nPoints()  ) << "invalid number of points for geopc";
            CHECK( __geopc1[__f][__p]->nPoints() ) << "invalid number of points for geopc1";
        }
    }

    auto const& initFace = boost::unwrap_ref( *__face_it );
    uint16_type __face_id = initFace.pos_first();
    gmc_ptrtype __c = __gm->template context<context>( initFace.element( 0 ), __geopc, __face_id );
    gmc1_ptrtype __c1 = __gm1->template context<context>( initFace.element( 0 ), __geopc1, __face_id );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type expr( M_expr, mapgmc );
    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
    t_expr1_type expr1( M_expr, mapgmc1 );




    size_type nbFaceDof = invalid_v<size_type>;

    for ( int e = 0; __face_it != __face_en; ++__face_it, ++e )
    {
        auto const& curFace = boost::unwrap_ref( *__face_it );
        FEELPP_ASSERT( curFace.isOnBoundary() && !curFace.isConnectedTo1() )
        ( curFace.marker() )
        ( curFace.isOnBoundary() )
        ( curFace.ad_first() )
        ( curFace.pos_first() )
        ( curFace.ad_second() )
        ( curFace.pos_second() )
        ( curFace.id() ).warn( "inconsistent data face" );
        DVLOG(2) << "[evaluator] FACE_ID = " << curFace.id()
                      << " element id= " << curFace.ad_first()
                      << " pos in elt= " << curFace.pos_first()
                      << " marker: " << curFace.marker() << "\n";
        DVLOG(2) << "[evaluator] FACE_ID = " << curFace.id() << " real pts=" << curFace.G() << "\n";

        uint16_type __face_id = curFace.pos_first();


        switch ( M_geomap_strategy )
        {
        default:
        case GeomapStrategyType::GEOMAP_OPT:
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->template update<context>( curFace.element( 0 ), __face_id );
            DVLOG(2) << "[evaluator::GEOMAP_HO|GEOMAP_OPT] FACE_ID = " << curFace.id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "[evaluator::GEOMAP_HO|GEOMAP_OPT] FACE_ID = " << curFace.id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            expr.update( mapgmc );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nRealDim; ++c1 )
                {
                    __p(e,p,c1) = __c->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                    {
                        __v( e,p,c1,c2) = expr.evalq( c1, c2, p );
                    }
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->template update<context>( curFace.element( 0 ), __face_id );
            DVLOG(2) << "[evaluator::GEOMAP_O1] FACE_ID = " << curFace.id() << "  ref pts=" << __c1->xRefs() << "\n";
            DVLOG(2) << "[evaluator::GEOMAP_O1] FACE_ID = " << curFace.id() << " real pts=" << __c1->xReal() << "\n";

            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

            expr1.update( mapgmc1 );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nRealDim; ++c1 )
                {
                    __p(e,p,c1) = __c1->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                    {
                        __v( e,p,c1,c2) = expr1.evalq( c1, c2, p );
                    }
                }

            }
        }
        break;
        }

    } // face_it


    return eval_element_type( __v, __p );
}
}
/// \endcond

/// \cond DETAIL
namespace detail
{
//template<typename Args>
template<typename ArgExprType,typename ArgPsetType,typename ArgRangeType>
struct evaluate
{
    using _expr_type = std::decay_t<ArgExprType>;
    using _pset_type = std::decay_t<ArgPsetType>;
    using _range_type = std::decay_t<ArgRangeType>;
    // typedef clean_type<Args,tag::expr> _expr_type;
    // typedef clean_type<Args,tag::pset> _pset_type;
    // typedef clean_type<Args,tag::range> _range_type;
    typedef details::Evaluator<EVAL_NODAL, _range_type, _pset_type, Expr<_expr_type> > eval_t;
    typedef typename eval_t::mesh_element_type mesh_element_type;
    typedef typename eval_t::eval_element_type element_type;
    static const uint16_type nDim = mesh_element_type::nDim;
    static const uint16_type nRealDim = mesh_element_type::nRealDim;
};
}
/// \endcond


template<typename IteratorRange, typename Pset, typename ExprT>
typename details::Evaluator<EVAL_NODAL, IteratorRange, Pset, Expr<ExprT> >::eval_element_type
evaluate_impl( IteratorRange const& range_it,
               Pset const& pset,
               Expr<ExprT> const& __expr,
               GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef details::Evaluator<EVAL_NODAL, IteratorRange, Pset, Expr<ExprT> > proj_t;
    proj_t p( range_it, pset, __expr, geomap );
    return p();
}


/**
 *
 * \brief evaluate an expression over a range of element at a set of points defined in the reference element
 *
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg pset set of points (e.g. quadrature points) in the reference elements to be transformed in the real elements
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 */
#if 0
BOOST_PARAMETER_FUNCTION(
    ( typename vf::detail::evaluate<Args>::element_type ), // return type
    evaluate,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( pset, * )
      ( expr, * )
    ) // 4. one required parameter, and

    ( optional
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
    )
#endif

template <typename ... Ts>
auto evaluate( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && pset = args.get(_pset);
    auto && expr = args.get(_expr);
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT );

    LOG(INFO) << "evaluate expression..." << std::endl;
    return evaluate_impl( range, pset, expr, geomap );
    LOG(INFO) << "evaluate expression done." << std::endl;
}

template<typename Dimensions>
std::array<int,2> indices2Array( const int& index, const Dimensions& dimensions, const int options )
{
    std::array<int,2> out;

    if ( options&Eigen::RowMajor )
    {
        div_t res = std::div( index, (int) (dimensions[3]*dimensions[2]*dimensions[1]) );
        out[0] = res.quot;
        res = std::div( res.rem, (int) (dimensions[3]*dimensions[2]) );
        out[1] = res.quot;
    }
    else
    {
        div_t res = std::div( index, (int) (dimensions[0]*dimensions[1]*dimensions[2]) );
        res = std::div( res.rem, (int) (dimensions[0]*dimensions[1]) );
        res = std::div( res.rem, (int) dimensions[0] );
        out[1] = res.quot;
        out[0] = res.rem;
    }

    return out;
}

/**
 * \brief data returned by normLinf
 */
template<int Dim>
struct normLinfData : public boost::tuple<double, Eigen::Matrix<double, Dim,1> >
{
    typedef boost::tuple<double, Eigen::Matrix<double, Dim,1> > super;

    normLinfData( double v, Eigen::Matrix<double, Dim,1> const& x  ) : super( v, x ) {}
    normLinfData() = default;
    normLinfData( normLinfData const&  ) = default;
    normLinfData( normLinfData &&  ) = default;
    normLinfData& operator=( normLinfData const&  ) = default;
    normLinfData& operator=( normLinfData &&  ) = default;

    /**
     * \return the maximum absolute value
     */
    double value() const { return this->template get<0>(); }

    /**
     * \return the maximum absolute value
     */
    double operator()() const { return this->template get<0>(); }

    /**
     * \return the point at which the expression is maximal
     */
    Eigen::Matrix<double, Dim,1> const& arg() const { return this->template get<1>(); }

    /**
     * Serialization for minmaxData
     */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->template get<0>();
            ar & this->template get<1>();
        }
};
/**
 *
 * \brief evaluate an expression over a range of element at a set of points defined in the reference element
 *
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg pset set of points (e.g. quadrature points) in the reference elements to be transformed in the real elements
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 */
#if 0
    BOOST_PARAMETER_FUNCTION(
    ( normLinfData<vf::detail::evaluate<Args>::nRealDim> ), // return type
    normLinf,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( pset, * )
      ( expr, * )
    ) // 4. one required parameter, and

    ( optional
      ( worldcomm,       (worldcomm_ptr_t), Environment::worldCommPtr() )
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
    )
)
#endif

template <typename ... Ts>
auto normLinf( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && pset = args.get(_pset);
    auto && expr = args.get(_expr);
    worldcomm_ptr_t worldcomm =args.get_else(_worldcomm,Environment::worldCommPtr() );
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT );

    // using _expr_type = std::decay_t<decltype(expr)>;
    // using _pset_type = std::decay_t<decltype(pset)>;
    // using _range_type = std::decay_t<decltype(range)>;
    //typedef typename details::Evaluator<EVAL_NODAL, _range_type, _pset_type, Expr<_expr_type>>::eval_element_type eval_element_type;
    using evaluate_helper_type = vf::detail::evaluate<decltype(expr),decltype(pset),decltype(range)>;
    using eval_element_type = typename evaluate_helper_type::eval_t;
    typedef typename eval_element_type::element_type element_type;

    int proc_number = worldcomm->globalRank();
    int world_size = worldcomm->size();
    constexpr int nRealDim = evaluate_helper_type::nRealDim;


    LOG(INFO) << "evaluate expression..." << std::endl;

    auto e = evaluate_impl( range, pset, expr, geomap );

    int index = 0;
    double maxe = 1e-30;
    Eigen::Matrix<double, nRealDim,1> n;

    if ( e.data().size() )
    {
        Eigen::Tensor<Eigen::DenseIndex, 0> arg_max =  e.data().abs().argmax() ;
        index = arg_max(0);
        maxe = std::abs( e.data()(index) );
    }

    if ( e.nodes().size() )
    {
        auto dimensions = e.data().dimensions();
        auto ids = indices2Array( index, dimensions, element_type::Options );

        int elem_id=ids[0];
        int pt_id=ids[1];

        for (int i=0; i<nRealDim; i++ )
            n[i]=e.nodes()( elem_id, pt_id ,i );
    }

    LOG(INFO) << "proc "<<proc_number<<" index at which function (size: " << e.data().size() << ") is maximal: "<< index << " coord = \n"<< n <<"\n";

    std::vector<normLinfData<nRealDim>> D_world( world_size );
    normLinfData<nRealDim> D( maxe, n );

    mpi::all_gather( worldcomm->globalComm(),
                     D,
                     D_world );

    auto it_max = std::max_element( D_world.begin() , D_world.end(),
                                    []( auto const& d1, auto const& d2 )
                                    { return d1.value() < d2.value(); } );

    int position = it_max - D_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global max = "<<it_max->value()<<" at position "<<position<<" with coord : \n "<<it_max->arg()<<"\n";

    // some extra check
    if (VLOG_IS_ON(2))
    {
        int index2=0;
        double maxe2 = 0;
        for( int i = 0; i < e.data().size(); ++i )
        {
            if ( math::abs(e.data()(i) ) > maxe2 )
            {
                maxe2 = math::abs(e.data()(i));
                index2 = i;
            }
        }
        LOG_ASSERT( index2 == index ) << " index2 = " << index2 <<  " and index  = " << index << "\n";
    }

    LOG(INFO) << "evaluate expression done." << std::endl;
    return normLinfData<nRealDim>( it_max->value(), it_max->arg() );

}


/**
 * \brief data returned by minmax
 */
template<int Dim>
struct minmaxData : public boost::tuple<double,double, Eigen::Matrix<double, Dim,2> >
{
    typedef boost::tuple<double,double, Eigen::Matrix<double, Dim,2> > super;
    minmaxData( super const& s ) : super( s ) {}
    minmaxData() = default;
    minmaxData( minmaxData const& ) = default;
    minmaxData( minmaxData && ) = default;
    minmaxData& operator=( minmaxData const& ) = default;
    minmaxData& operator=( minmaxData && ) = default;

    /**
     * \return the minimum absolute value
     */
    double min() const { return this->template get<0>(); }

    /**
     * \return the maximum absolute value
     */
    double max() const { return this->template get<1>(); }


    /**
     * \return the point at which the expression is maximal
     */
    Eigen::Matrix<double, Dim,1> argmin() const { return this->template get<2>().col(0); }

    /**
     * \return the point at which the expression is maximal
     */
    Eigen::Matrix<double, Dim,1> argmax() const { return this->template get<2>().col(1); }

    /**
     * coordinates of the points where  min and max are attained
     */
    Eigen::Matrix<double, Dim,2> const& coords() const { return this->template get<2>(); }

    /**
     * Serialization for minmaxData
     */
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->template get<0>();
            ar & this->template get<1>();
            ar & this->template get<2>();
        }

};


/**
 *
 * \brief evaluate an expression over a range of element at a set of points defined in the reference element
 *
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg pset set of points (e.g. quadrature points) in the reference elements to be transformed in the real elements
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 */
#if 0
BOOST_PARAMETER_FUNCTION(
    ( minmaxData<vf::detail::evaluate<Args>::nRealDim> ), // return type
    minmax,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( pset, * )
      ( expr, * )
    ) // 4. one required parameter, and

    ( optional
      ( worldcomm,       (worldcomm_ptr_t), Environment::worldCommPtr() )
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
    )
)
#endif

template <typename ... Ts>
auto minmax( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && pset = args.get(_pset);
    auto && expr = args.get(_expr);
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT );
    worldcomm_ptr_t worldcomm = args.get_else(_worldcomm,Environment::worldCommPtr() );

    // typedef detail::clean_type<Args,tag::expr> _expr_type;
    // typedef detail::clean_type<Args,tag::pset> _pset_type;
    // typedef detail::clean_type<Args,tag::range> _range_type;
    // typedef typename details::Evaluator<EVAL_NODAL, _range_type, _pset_type, Expr<_expr_type>>::eval_element_type eval_element_type;

    using evaluate_helper_type = vf::detail::evaluate<decltype(expr),decltype(pset),decltype(range)>;
    using eval_element_type = typename evaluate_helper_type::eval_t;
    typedef typename eval_element_type::element_type element_type;

    constexpr int nRealDim = evaluate_helper_type::nRealDim;
    int proc_number = worldcomm->globalRank();

    LOG(INFO) << "evaluate minmax(expression)..." << std::endl;
    auto e = evaluate_impl( range, pset, expr, geomap );

    int indexmin = 0;
    int indexmax = 0;
    double mine = std::numeric_limits<double>::max();
    double maxe = std::numeric_limits<double>::lowest();

    if ( e.data().size() )
    {
        Eigen::Tensor<Eigen::DenseIndex, 0> arg_max = e.data().argmax() ;
        Eigen::Tensor<Eigen::DenseIndex, 0> arg_min = e.data().argmin() ;

        indexmin = arg_min(0);
        indexmax = arg_max(0);
        mine = e.data()(indexmin);
        maxe = e.data()(indexmax);
    }

    Eigen::Matrix<double, nRealDim,2> n = Eigen::Matrix<double, nRealDim,2>::Zero();

    if ( e.nodes().size() )
    {
        auto dimensions = e.data().dimensions();
        const int options = element_type::Options;

        auto ids = indices2Array( indexmin, dimensions, options );
        int elem_id=ids[0];
        int pt_id=ids[1];
        for (int i=0; i<nRealDim; i++ )
            n(i,0)=e.nodes()( elem_id, pt_id ,i );

        ids = indices2Array( indexmax, dimensions, options );
        elem_id=ids[0];
        pt_id=ids[1];
        for (int i=0; i<nRealDim; i++ )
            n(i,1)=e.nodes()( elem_id, pt_id ,i );
    }

    LOG(INFO) << "proc "<<proc_number <<" index at which function (size: " << e.data().size() << ") is minimal: " << indexmin << " coord = \n"<<n.col(0)<<"\n";
    LOG(INFO) << "proc "<<proc_number <<" index at which function (size: " << e.data().size() << ") is maximal: " << indexmax << " coord = \n"<<n.col(1)<<"\n";

    int world_size = worldcomm->size();
    minmaxData<nRealDim> D( boost::make_tuple( mine, maxe, n) );
    std::vector<minmaxData<nRealDim>> D_world( world_size );
    mpi::all_gather( worldcomm->globalComm(), D, D_world );

    auto it_min = std::min_element( D_world.begin() , D_world.end(),
                                    []( auto const& d1, auto const& d2 )
                                    { return d1.min() < d2.min(); } );
    int positionmin = it_min - D_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global min = "<<it_min->min()<<" at position "<<positionmin<<" with coord : \n "<< it_min->argmin()<<"\n";

    auto it_max = std::max_element( D_world.begin() , D_world.end(),
                                    []( auto const& d1, auto const& d2 )
                                    { return d1.max() < d2.max(); } );
    int positionmax = it_max - D_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global max = "<<it_max->max()<<" at position "<<positionmax<<" with coord : \n "<< it_max->argmax()<<"\n";

    Eigen::Matrix<double, nRealDim,2> coords;
    coords.col(0) = it_min->argmin();
    coords.col(1) = it_max->argmax();

    LOG(INFO) << "evaluate minmax(expression) done." << std::endl;

    return minmaxData<nRealDim> (boost::make_tuple( it_min->min(), it_max->max(), coords ));
}

#if 0
/// \cond DETAIL
namespace detail{
template <typename Args>
struct maxPerCellData
{
    typedef clean_type<Args,tag::expr> _expr_type;
    typedef typename _expr_type::value_type value_type;

    typedef Eigen::Tensor<value_type,3> element_type;
}; // maxpercelldata

} //detail
#endif

#if 0
BOOST_PARAMETER_FUNCTION(
    ( typename vf::detail::maxPerCellData<Args>::element_type ), // return type
    maxPerCell,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( pset, * )
      ( expr, * )
      ) // 4. one required parameter, and

    ( optional
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
      )
#endif

template <typename ... Ts>
auto maxPerCell( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && pset = args.get(_pset);
    auto && expr = args.get(_expr);
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT );

    using _expr_type = std::decay_t<decltype(expr)>;
    typedef typename _expr_type::value_type value_type;

    auto e = evaluate_impl( range, pset, expr, geomap );
    std::array<std::ptrdiff_t,1> reduction_axis;
    reduction_axis[0] = 2;

    Eigen::Tensor<value_type,4> data = e.data();

    Eigen::Tensor<value_type,3> reduced_data = data.maximum( reduction_axis );

    return reduced_data;
}

} // vf
} // feel


#endif /* __FEELPP_EVALUATORS_H */
