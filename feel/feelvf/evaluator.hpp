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


    static const uint16_type imorder = 1;
    static const bool imIsPoly = true;

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef typename iterator_type::value_type mesh_element_type;
    typedef IteratorRange range_iterator;
    typedef typename mpl::if_<mpl::bool_<mesh_element_type::is_simplex>,
                              mpl::identity<typename Pset::template apply<mesh_element_type::nRealDim, value_type, Simplex>::type >,
                              mpl::identity<typename Pset::template apply<mesh_element_type::nRealDim, value_type, Hypercube>::type >
                              >::type::type pointset_type;
    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> element_type;
    typedef Eigen::Matrix<value_type,mesh_element_type::nRealDim,Eigen::Dynamic> node_type;
    typedef boost::tuple<element_type,node_type> eval_element_type;
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
        M_pset(),
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
    typedef typename gm_type::template Context<context, mesh_element_type> gm_context_type;
    typedef typename mesh_element_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<context, mesh_element_type> gm1_context_type;


    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef boost::shared_ptr<gm1_context_type> gm1_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm1_context_ptrtype> > map_gmc1_type;
    //typedef typename expression_type::template tensor<map_gmc_type,fusion::map<fusion::pair<vf::detail::gmc<0>,boost::shared_ptr<fecontext_type> > > > t_expr_type;
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


    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    int npoints = M_pset.points().size2();
    element_type __v( M_pset.points().size2()*std::distance( it, en )*shape::M );
    node_type __p( mesh_element_type::nDim, M_pset.points().size2()*std::distance( it, en ) );
    __v.setZero();
    __p.setZero();



    // return if no elements
    if ( it == en )
        return boost::make_tuple( __v, __p );

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( it->gm(),
                                                                                         M_pset.points() ) );
    typename gm1_type::precompute_ptrtype __geopc1( new typename gm1_type::precompute_type( it->gm1(),
                                                                                            M_pset.points() ) );



    gm_context_ptrtype __c( new gm_context_type( it->gm(),*it,__geopc ) );
    gm1_context_ptrtype __c1( new gm1_context_type( it->gm1(),*it,__geopc1 ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type tensor_expr( M_expr, mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( M_expr, mapgmc1 );

    for ( int e = 0; it!=en ; ++it, ++e )
    {
        switch ( M_geomap_strategy )
        {
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                {
                    __p(c1, e*npoints+p) = __c->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    __v( e*npoints*shape::M+shape::M*p+c1) = tensor_expr.evalq( c1, 0, p );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( *it );
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            tensor_expr1.update( mapgmc1 );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                {
                    __p(c1, e*npoints+p) = __c1->xReal(p)[c1];
                }
                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    __v( e*npoints*shape::M+shape::M*p+c1) = tensor_expr1.evalq( c1, 0, p );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_OPT:
        {
            if ( it->isOnBoundary() )
            {
                // HO if on boundary
                __c->update( *it );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                tensor_expr.update( mapgmc );

                for ( uint16_type p = 0; p < npoints; ++p )
                {
                    for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                    {
                        __p(c1, e*npoints+p) = __c->xReal(p)[c1];
                    }
                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    {
                        __v( e*npoints*shape::M+shape::M*p+c1) = tensor_expr.evalq( c1, 0, p );
                    }
                }
            }

            else
            {
                __c1->update( *it );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                tensor_expr1.update( mapgmc1 );


                for ( uint16_type p = 0; p < npoints; ++p )
                {
                    for ( uint16_type c1 = 0; c1 < mesh_element_type::nDim; ++c1 )
                    {
                        __p(c1, e*npoints+p) = __c1->xReal(p)[c1];
                    }
                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    {
                        __v( e*npoints*shape::M+shape::M*p+c1) = tensor_expr1.evalq( c1, 0, p );
                    }
                }
            }
        }
        break;
        }
    }
    return boost::make_tuple( __v, __p );
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
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef boost::shared_ptr<gm1_type> gm1_ptrtype;

    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm1_type::template Context<context, geoelement_type> gmc1_type;
    typedef boost::shared_ptr<gmc1_type> gmc1_ptrtype;
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

    iterator_type __face_it, __face_en;
    boost::tie( boost::tuples::ignore, __face_it, __face_en ) = M_range;

    int npoints = M_pset.fpoints(0,1).size2();
    element_type __v( M_pset.fpoints(0,1).size2()*std::distance( __face_it, __face_en )*shape::M );
    node_type __p( mesh_element_type::nRealDim, M_pset.fpoints(0,1).size2()*std::distance( __face_it, __face_en ) );
    __v.setZero();
    __p.setZero();
    VLOG(2) << "pset: " << M_pset.fpoints(0,1);
    VLOG(2) << "Checking trivial result...";

    if ( __face_it == __face_en )
        return boost::make_tuple( __v, __p );

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

    uint16_type __face_id = __face_it->pos_first();
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );
    gmc1_ptrtype __c1( new gmc1_type( __gm1, __face_it->element( 0 ), __geopc1, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type expr( M_expr, mapgmc );
    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
    t_expr1_type expr1( M_expr, mapgmc1 );




    size_type nbFaceDof = invalid_size_type_value;

    for ( int e = 0; __face_it != __face_en; ++__face_it, ++e )
    {
        FEELPP_ASSERT( __face_it->isOnBoundary() && !__face_it->isConnectedTo1() )
        ( __face_it->marker() )
        ( __face_it->isOnBoundary() )
        ( __face_it->ad_first() )
        ( __face_it->pos_first() )
        ( __face_it->ad_second() )
        ( __face_it->pos_second() )
        ( __face_it->id() ).warn( "inconsistent data face" );
        DVLOG(2) << "[evaluator] FACE_ID = " << __face_it->id()
                      << " element id= " << __face_it->ad_first()
                      << " pos in elt= " << __face_it->pos_first()
                      << " marker: " << __face_it->marker() << "\n";
        DVLOG(2) << "[evaluator] FACE_ID = " << __face_it->id() << " real pts=" << __face_it->G() << "\n";

        uint16_type __face_id = __face_it->pos_first();


        switch ( M_geomap_strategy )
        {
        default:
        case GeomapStrategyType::GEOMAP_OPT:
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( __face_it->element( 0 ), __face_id );
            DVLOG(2) << "[evaluator::GEOMAP_HO|GEOMAP_OPT] FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "[evaluator::GEOMAP_HO|GEOMAP_OPT] FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            expr.update( mapgmc );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nRealDim; ++c1 )
                {
                    __p(c1, e*npoints+p) = __c->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    __v( e*npoints*shape::M+shape::M*p+c1) = expr.evalq( c1, 0, p );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( __face_it->element( 0 ), __face_id );
            DVLOG(2) << "[evaluator::GEOMAP_O1] FACE_ID = " << __face_it->id() << "  ref pts=" << __c1->xRefs() << "\n";
            DVLOG(2) << "[evaluator::GEOMAP_O1] FACE_ID = " << __face_it->id() << " real pts=" << __c1->xReal() << "\n";

            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

            expr1.update( mapgmc1 );

            for ( uint16_type p = 0; p < npoints; ++p )
            {
                for ( uint16_type c1 = 0; c1 < mesh_element_type::nRealDim; ++c1 )
                {
                    __p(c1, e*npoints+p) = __c1->xReal(p)[c1];
                }

                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                {
                    __v( e*npoints*shape::M+shape::M*p+c1) = expr1.evalq( c1, 0, p );
                }

            }
        }
        break;
        }

    } // face_it


    return boost::make_tuple( __v, __p );
}
}
/// \endcond

/// \cond DETAIL
namespace detail
{
template<typename Args>
struct evaluate
{
    typedef typename clean_type<Args,tag::expr>::type _expr_type;
    typedef typename clean_type<Args,tag::pset>::type _pset_type;
    typedef typename clean_type<Args,tag::range>::type _range_type;
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
)
{
    LOG(INFO) << "evaluate expression..." << std::endl;
    return evaluate_impl( range, pset, expr, geomap );
    LOG(INFO) << "evaluate expression done." << std::endl;
}

/**
 * \brief data returned by normLinf
 */
template<int Dim>
struct normLinfData : public boost::tuple<double, Eigen::Matrix<double, Dim,1> >
{
    typedef boost::tuple<double, Eigen::Matrix<double, Dim,1> > super;
    normLinfData( double v, Eigen::Matrix<double, Dim,1> const& x  ) : super( v, x ) {}

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
      ( worldcomm,       (WorldComm), Environment::worldComm() )
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
    )
)
{

    int proc_number = worldcomm.globalRank();

    LOG(INFO) << "evaluate expression..." << std::endl;

    auto e = evaluate_impl( range, pset, expr, geomap );
    int index;
    double maxe = e.template get<0>().array().abs().maxCoeff(&index);

    Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,1> n = e.template get<1>().col(index);
    LOG(INFO) << "proc "<<proc_number<<" index at which function (size: " << e.template get<0>().array().size() << ") is maximal: "<< index << " coord = \n"<<n<<"\n";

    int world_size = worldcomm.size();
    std::vector<double> maxe_world( world_size );
    mpi::all_gather( worldcomm.globalComm(),
                     maxe,
                     maxe_world );

    std::vector< Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,1> > n_world( world_size );
    mpi::all_gather( worldcomm.globalComm(),
                     n,
                     n_world );

    auto it_max = std::max_element( maxe_world.begin() , maxe_world.end() );
    int position = it_max - maxe_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global max = "<<*it_max<<" at position "<<position<<" with coord : \n "<<n<<"\n";

    // some extra check
    if (VLOG_IS_ON(2))
    {
        int index2=0;
        double maxe2 = 0;
        for( int i = 0; i < e.template get<0>().size(); ++i )
        {
            if ( math::abs(e.template get<0>()(i)) > maxe2 )
            {
                maxe2 = math::abs(e.template get<0>()(i));
                index2 = i;
            }
        }
        LOG_ASSERT( index2 == index ) << " index2 = " << index2 <<  " and index  = " << index << "\n";
    }
    LOG(INFO) << "evaluate expression done." << std::endl;
    Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,1> x = n_world[position];
    return normLinfData<vf::detail::evaluate<Args>::nRealDim> (*it_max, x);
}


/**
 * \brief data returned by minmax
 */
template<int Dim>
struct minmaxData : public boost::tuple<double,double, Eigen::Matrix<double, Dim,2> >
{
    typedef boost::tuple<double,double, Eigen::Matrix<double, Dim,2> > super;
    minmaxData( super const& s ) : super( s ) {}
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
      ( worldcomm,       (WorldComm), Environment::worldComm() )
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
    )
)
{

    int proc_number = worldcomm.globalRank();

    LOG(INFO) << "evaluate minmax(expression)..." << std::endl;

    auto e = evaluate_impl( range, pset, expr, geomap );
    int indexmin;
    int indexmax;
    double mine = e.template get<0>().array().minCoeff(&indexmin);
    double maxe = e.template get<0>().array().maxCoeff(&indexmax);

    Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,2> n;
    n.col(0) = e.template get<1>().col(indexmin);
    n.col(1) = e.template get<1>().col(indexmax);
    LOG(INFO) << "proc "<<proc_number<<" index at which function (size: " << e.template get<0>().array().size() << ") is minimal: "<< indexmin << " coord = \n"<<n.col(0)<<"\n";
    LOG(INFO) << "proc "<<proc_number<<" index at which function (size: " << e.template get<0>().array().size() << ") is maximal: "<< indexmax << " coord = \n"<<n.col(1)<<"\n";

    int world_size = worldcomm.size();
    std::vector<double> mine_world( world_size );
    std::vector<double> maxe_world( world_size );
    mpi::all_gather( worldcomm.globalComm(),mine,mine_world );
    mpi::all_gather( worldcomm.globalComm(),maxe,maxe_world );

    std::vector< Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,2> > n_world( world_size );
    mpi::all_gather( worldcomm.globalComm(),n,n_world );

    auto it_min = std::min_element( mine_world.begin() , mine_world.end() );
    int positionmin = it_min - mine_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global min = "<<*it_min<<" at position "<<positionmin<<" with coord : \n "<<n_world[positionmin]<<"\n";
    auto it_max = std::max_element( maxe_world.begin() , maxe_world.end() );
    int positionmax = it_max - maxe_world.begin();
    LOG(INFO)<<"proc "<<proc_number<<" : global max = "<<*it_max<<" at position "<<positionmax<<" with coord : \n "<<n_world[positionmax]<<"\n";

    Eigen::Matrix<double, vf::detail::evaluate<Args>::nRealDim,2> coords;
    coords.col(0) = n_world[positionmin].col(0);
    coords.col(1) = n_world[positionmax].col(1);

    LOG(INFO) << "evaluate minmax(expression) done." << std::endl;

    return minmaxData<vf::detail::evaluate<Args>::nRealDim> (boost::make_tuple( *it_min, *it_max, coords ));
}

} // vf
} // feel


#endif /* __FEELPP_EVALUATORS_H */
