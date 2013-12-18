/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-05-31

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Université Joseph Fourier (Grenoble I)

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
   \file projectors.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-05-31
 */
#ifndef __Projectors_H
#define __Projectors_H 1

#include <boost/timer.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace vf
{
namespace details
{
/**
 * \class Projector
 * \brief base class for projectors
 *
 * @author Christophe Prud'homme
 * @see Projector3D, Projector2D, ProjectorOn
 */
template<ProjectorType iDim, typename FunctionSpaceType, typename IteratorRange, typename ExprT>
class Projector
{
public:


    /** @name Typedefs
     */
    //@{

    static const size_type context = ExprT::context|vm::POINT;

    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::basis_type basis_type;
    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;


    static const uint16_type imorder = functionspace_type::basis_type::nOrder;
    static const bool imIsPoly = true;

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Projector( functionspace_ptrtype const& __functionspace,
               IteratorRange const& r,
               expression_type const& __expr,
               GeomapStrategyType geomap_strategy )
        :
        M_functionspace( __functionspace ),
        M_range( r ),
        M_expr( __expr ),
        M_geomap_strategy( geomap_strategy )
    {
        DVLOG(2) << "Projector constructor from expression\n";
    }


    Projector( Projector const& __vfi )
        :
        M_functionspace( __vfi.M_functionspace ),
        M_range( __vfi.M_range ),
        M_expr( __vfi.M_expr ),
        M_geomap_strategy( __vfi.M_geomap_strategy )
    {
        DVLOG(2) << "Projector copy constructor\n";
    }

    virtual ~Projector() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    element_type operator()( const bool sum = false ) const
    {
        return operator()( sum, idim_type() );
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

    element_type operator()( const bool sum, mpl::size_t<MESH_ELEMENTS> ) const;
    element_type operator()( const bool sum, mpl::size_t<MESH_FACES> ) const;
    element_type operator()( const bool sum, mpl::size_t<MESH_POINTS> ) const;

private:

    functionspace_ptrtype const& M_functionspace;
    range_iterator M_range;
    expression_type const&  M_expr;
    GeomapStrategyType M_geomap_strategy;
};

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( const bool sum, mpl::size_t<MESH_ELEMENTS> ) const
{
    boost::timer __timer;

    element_type __v( M_functionspace );
    FEELPP_ASSERT( __v.size() == M_functionspace->dof()->nDof() )( __v.size() )( M_functionspace->dof()->nDof() ).warn( "invalid size" );
    __v.setZero();

    typedef typename functionspace_type::fe_type fe_type;
    fe_type* __fe = __v.functionSpace()->fe().get();

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type element_type;

    typedef typename element_type::gm_type gm_type;
    typedef typename gm_type::template Context<context, element_type> gm_context_type;
    typedef typename element_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<context, element_type> gm1_context_type;


    typedef typename fe_type::template Context<context, fe_type, gm_type, element_type> fecontext_type;
    typedef typename fe_type::template Context<context, fe_type, gm1_type, element_type> fecontext1_type;

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

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( __v.functionSpace()->gm(),
            __fe->points() ) );
    typename gm1_type::precompute_ptrtype __geopc1( new typename gm1_type::precompute_type( __v.mesh()->gm1(),
            __fe->points() ) );



    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = M_range;

    // return if no elements
    if ( it == en )
        return __v;

    gm_context_ptrtype __c( new gm_context_type( __v.functionSpace()->gm(),*it,__geopc ) );
    gm1_context_ptrtype __c1( new gm1_context_type( __v.mesh()->gm1(),*it,__geopc1 ) );

    typedef typename t_expr_type::shape shape;
    static const bool is_rank_ok = ( shape::M == FunctionSpaceType::nComponents1 &&
                                     shape::N == FunctionSpaceType::nComponents2 );

    BOOST_MPL_ASSERT_MSG( mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          ( mpl::int_<shape::M>, mpl::int_<FunctionSpaceType::nComponents>, shape ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type tensor_expr( basis_type::isomorphism( M_expr ), mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( basis_type::isomorphism( M_expr ), mapgmc1 );

    std::vector<bool> points_done( __v.functionSpace()->dof()->nLocalDof()/__v.nComponents );
    std::fill( points_done.begin(), points_done.end(),false );

    for ( ; it!=en ; ++it )
    {
        switch ( M_geomap_strategy )
        {
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {
                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                {
                    if ( sum )
                        __v.plus_assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );

                    else
                        __v.assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( *it );
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            tensor_expr1.update( mapgmc1 );

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {
                for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                    //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                {
                    if ( sum )
                        __v.plus_assign( it->id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );

                    else
                        __v.assign( it->id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );
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

                for ( uint16_type __j = 0; __j < ndofv; ++__j )
                {
                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                        //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        if ( sum )
                            __v.plus_assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );

                        else
                            __v.assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                    }
                }
            }

            else
            {
                __c1->update( *it );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                tensor_expr1.update( mapgmc1 );

                for ( uint16_type __j = 0; __j < ndofv; ++__j )
                {
                    for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                        //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        if ( sum )
                            __v.plus_assign( it->id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );

                        else
                            __v.assign( it->id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );
                    }
                }
            }
        }
        break;
        }

        //if P0 continuous finish loop here
        if (fe_type::isLagrangeP0Continuous )
        {
            break;
        }
    }


    return __v;
}

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( const bool sum, mpl::size_t<MESH_FACES> ) const
{
    boost::timer __timer;

    element_type __v( M_functionspace );
    __v.setZero();

    DVLOG(2) << "call project(MESH_FACES) " << "\n";
    //
    // a few typedefs
    //

    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename geoelement_type::face_type face_type;

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


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    typedef typename fe_type::template Context< context, fe_type, gm1_type, geoelement_type> fecontext1_type;
    typedef boost::shared_ptr<fecontext1_type> fecontext1_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

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
    DVLOG(2)  << "assembling Dirichlet conditions\n";

    dof_type const* __dof = __v.functionSpace()->dof().get();

    fe_type const* __fe = __v.functionSpace()->fe().get();

    iterator_type __face_it, __face_en;
    boost::tie( boost::tuples::ignore, __face_it, __face_en ) = M_range;

    if ( __face_it == __face_en )
        return __v;

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

    DVLOG(2)  << "[integratoron] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
    std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );
    std::vector<std::map<permutation_type, geopc1_ptrtype> > __geopc1( geoelement_type::numTopologicalFaces );

    for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
                __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
            __geopc1[__f][__p] = geopc1_ptrtype(  new geopc1_type( __gm1, __fe->points( __f ) ) );
            //DVLOG(2) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
            FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points for geopc" );
            FEELPP_ASSERT( __geopc1[__f][__p]->nPoints() ).error( "invalid number of points for geopc1" );
        }
    }

    uint16_type __face_id = __face_it->pos_first();
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );
    gmc1_ptrtype __c1( new gmc1_type( __gm1, __face_it->element( 0 ), __geopc1, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type expr( basis_type::isomorphism( M_expr ), mapgmc );
    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
    t_expr1_type expr1( basis_type::isomorphism( M_expr ), mapgmc1 );




    size_type nbFaceDof = invalid_size_type_value;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    DVLOG(2)  << "[projector::operator(MESH_FACES)] nbFaceDof = " << nbFaceDof << "\n";

    std::vector<int> dofs;
    std::vector<value_type> values;

    for ( ; __face_it != __face_en; ++__face_it )
    {
        FEELPP_ASSERT( __face_it->isOnBoundary() && !__face_it->isConnectedTo1() )
        ( __face_it->marker() )
        ( __face_it->isOnBoundary() )
        ( __face_it->ad_first() )
        ( __face_it->pos_first() )
        ( __face_it->ad_second() )
        ( __face_it->pos_second() )
        ( __face_it->id() ).warn( "inconsistent data face" );
        DVLOG(2) << "[projector] FACE_ID = " << __face_it->id()
                      << " element id= " << __face_it->ad_first()
                      << " pos in elt= " << __face_it->pos_first()
                      << " marker: " << __face_it->marker() << "\n";
        DVLOG(2) << "[projector] FACE_ID = " << __face_it->id() << " real pts=" << __face_it->G() << "\n";

        uint16_type __face_id = __face_it->pos_first();

        std::pair<size_type,size_type> range_dof( std::make_pair( __v.start(),
                __v.functionSpace()->nDof() ) );
        DVLOG(2)  << "[projector] dof start = " << range_dof.first << "\n";
        DVLOG(2)  << "[projector] dof range = " << range_dof.second << "\n";

        switch ( M_geomap_strategy )
        {
        default:
        case GeomapStrategyType::GEOMAP_OPT:
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( __face_it->element( 0 ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            expr.update( mapgmc );

            for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                {
                    for ( uint16_type l = 0; l < nbFaceDof; ++l )
                    {
                        typename expression_type::value_type __value = expr.evalq( c1, c2, l );

                        size_type thedof =  __v.start() +
                                            boost::get<0>( __dof->faceLocalToGlobal( __face_it->id(), l, c1 ) );

                        if ( sum )
                            __v( thedof ) +=  __value;

                        else
                            __v( thedof ) =  __value;
                    }
                }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( __face_it->element( 0 ), __face_id );
            DVLOG(2) << "[projector] FACE_ID = " << __face_it->id() << "  ref pts=" << __c1->xRefs() << "\n";
            DVLOG(2) << "[projector] FACE_ID = " << __face_it->id() << " real pts=" << __c1->xReal() << "\n";

            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

            expr1.update( mapgmc1 );


            for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                {
                    for ( uint16_type l = 0; l < nbFaceDof; ++l )
                    {
                        typename expression_type::value_type __value = expr1.evalq( c1, c2, l );

                        size_type thedof =  __v.start() +
                                            boost::get<0>( __dof->faceLocalToGlobal( __face_it->id(), l, c1 ) );

                        if ( sum )
                            __v( thedof ) +=  __value;

                        else
                            __v( thedof ) =  __value;
                    }
                }
        }
        break;
        }

    } // face_it

    return __v;
}

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( const bool sum, mpl::size_t<MESH_POINTS> ) const
{
    boost::timer __timer;

    element_type __v( M_functionspace );
    __v.setZero();

    iterator_type pt_it, pt_en;
    boost::tie( boost::tuples::ignore, pt_it, pt_en ) = M_range;

    if ( pt_it == pt_en )
        return __v;
#if 0
    BOOST_FOREACH( auto dof, M_functionspace->dof()->markerToDof( pt_it->marker() ) );
    {

        // get the first element to which the point/dof belong and then build
        // the proper geomap context in order to evaluate the expression at the
        // point

#if 0
        if ( sum )
            __v.add( dof.second,  );
#endif
    }
#endif
    return __v;
}

} // detail
/// \endcond


/**
 * \brief nodal projection of \p __expr onto the the subspace of \p __functionspace described by the range \p range_it
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr over the range \p __range_it
 */
template<typename FunctionSpaceType, typename IteratorRange, typename ExprT>
typename FunctionSpaceType::element_type
project( boost::shared_ptr<FunctionSpaceType> const& __functionspace,
         IteratorRange const& range_it,
         Expr<ExprT> const& __expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr, geomap );
    return p();
}

template<typename FunctionSpaceType, typename IteratorRange, typename ExprT>
typename FunctionSpaceType::element_type
project_impl( boost::shared_ptr<FunctionSpaceType> const& __functionspace,
              IteratorRange const& range_it,
              Expr<ExprT> const& __expr,
              GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr, geomap );
    return p();
}

/**
 * \brief nodal projection of \p __expr onto the space \p __functionspace
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr
 */
template<typename FunctionSpaceType, typename ExprT>
typename FunctionSpaceType::element_type
project( boost::shared_ptr<FunctionSpaceType> const& __functionspace, Expr<ExprT> const& __expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    return project( __functionspace, elements( __functionspace->mesh() ), __expr, geomap );
}

/**
 * \brief nodal projection of \p __expr onto the the subspace of \p __functionspace described by the range \p range_it
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr over the range \p __range_it
 */
template<typename FunctionSpaceType, typename IteratorRange, typename ExprT>
typename FunctionSpaceType::element_type
sum( boost::shared_ptr<FunctionSpaceType> const& __functionspace,
     IteratorRange const& range_it,
     Expr<ExprT> const& __expr,
     GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr, geomap );
    return p( true );
}

/**
 * \brief nodal projection of \p __expr onto the space \p __functionspace
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr
 */
template<typename FunctionSpaceType, typename ExprT>
typename FunctionSpaceType::element_type
sum( boost::shared_ptr<FunctionSpaceType> const& __functionspace, Expr<ExprT> const& __expr )
{
    return sum( __functionspace, elements( __functionspace->mesh() ), __expr );
}

/**
 * \brief nodal projection of \p __expr onto the space \p __functionspace
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr
 */
template<typename FunctionSpaceType, typename ExprT>
typename FunctionSpaceType::element_type
project( FunctionSpaceType const& __functionspace, Expr<ExprT> const& __expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef __typeof__( __functionspace->mesh()->elementsRange() ) IteratorRange;
    typedef details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, __functionspace->mesh()->elementsRange(), __expr, geomap );
    return p();
}

/**
 * \brief nodal projection of \p __expr onto the the subspace of \p __functionspace described by the range \p range_it
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr over the range \p __range_it
 */
template<typename FunctionSpaceType, typename IteratorRange, typename ExprT>
typename FunctionSpaceType::element_type
project( FunctionSpaceType const& __functionspace,
         IteratorRange const& range_it,
         Expr<ExprT> const& __expr,
         GeomapStrategyType geomap = GeomapStrategyType::GEOMAP_HO )
{
    typedef details::Projector<NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr, geomap );
    return p();
}

/// \cond DETAIL
namespace detail
{
template<typename S>
struct space_ptr
{
    typedef typename S::element_type type;
};

template<typename S>
struct space_value
{
    typedef S type;
};

template<typename Args>
struct project
{
    typedef typename clean_type<Args,tag::space>::type the_space_type;
    typedef typename mpl::if_<is_shared_ptr<the_space_type>,
                              mpl::identity<space_ptr<the_space_type> >,
                              mpl::identity<space_value<the_space_type> > >::type::type space_type;
    typedef typename space_type::type _space_type;
    typedef boost::shared_ptr<_space_type> _space_ptrtype;
    typedef typename _space_type::element_type element_type;
    //typedef typename clean_type<Args,tag::expr>::type _expr_type;
    //typedef typename clean_type<Args,tag::range>::type _range_type;
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
    ( typename vf::detail::project<Args>::element_type ), // return type
    project,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( space, *( boost::is_convertible<mpl::_,boost::shared_ptr<Feel::FunctionSpaceBase> > ) )
      ( expr, * )
    ) // 4. one required parameter, and

    ( optional
      ( range,          *, elements( space->mesh() )  )
      ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
      ( accumulate,     *( boost::is_integral<mpl::_> ), false )
    )
)
{
#if 0
    typedef typename vf::detail::project<Args>::_space_type _space_type;
    typedef typename vf::detail::project<Args>::_range_type _range_type;
    typedef typename vf::detail::project<Args>::_expr_type _expr_type;

    typedef details::Projector<NODAL, _space_type, _range_type, Expr<_expr_type> > proj_t;
    proj_t p( space, range, expr,geomap );
    return p( sum );
#else
    //    if ( accumulate  )
    //return sum( space, range, expr, geomap );
    return project_impl( space, range, expr, geomap );
#endif
}

} // vf
} // feel


#endif /* __Projectors_H */
