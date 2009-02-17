/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-05-31

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
   \file projectors.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-05-31
 */
#ifndef __Projectors_H
#define __Projectors_H 1

#include <boost/timer.hpp>
#include <boost/signal.hpp>

namespace Life
{
namespace vf
{
/// \cond detail
enum ProjectorType
{
    PROJ_NODAL = 0,             /**< Nodal projection */
    PROJ_L2 = 1,                /**< L2 projection */
    PROJ_H1 = 2,                /**< H1 projection */
};
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
    typedef ExprT expression_type;
    typedef typename expression_type::value_type value_type;

    typedef typename boost::tuples::template element<0, IteratorRange>::type idim_type;
    typedef typename boost::tuples::template element<1, IteratorRange>::type iterator_type;
    typedef IteratorRange range_iterator;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Projector( functionspace_ptrtype const& __functionspace,
               IteratorRange const& r,
               expression_type const& __expr )
        :
        _M_functionspace( __functionspace ),
        _M_range( r ),
        _M_expr( __expr )
        {
            Debug( 5065 ) << "Projector constructor from expression\n";
        }


    Projector( Projector const& __vfi )
        :
        _M_functionspace( __vfi._M_functionspace ),
        _M_range( __vfi._M_range ),
        _M_expr( __vfi._M_expr )
        {
            Debug( 5065 ) << "Projector copy constructor\n";
        }

    virtual ~Projector() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    element_type operator()() const
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
    expression_type const& expression() const { return _M_expr; }

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

    element_type operator()( mpl::size_t<MESH_ELEMENTS> ) const;
    element_type operator()( mpl::size_t<MESH_FACES> ) const;

private:

    functionspace_ptrtype const& _M_functionspace;
    range_iterator _M_range;
    expression_type const&  _M_expr;
};

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( mpl::size_t<MESH_ELEMENTS> ) const
{
    boost::timer __timer;

    element_type __v( _M_functionspace );
    __v.clear();

    typedef typename functionspace_type::fe_type fe_type;
    fe_type* __fe = __v.functionSpace()->fe().get();

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type element_type;

    typedef typename functionspace_type::gm_type gm_type;
    typedef typename gm_type::template Context<context, element_type> gm_context_type;

    typedef typename fe_type::template Context<context, fe_type, gm_type, element_type> fecontext_type;

    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    //typedef typename expression_type::template tensor<map_gmc_type,fusion::map<fusion::pair<detail::gmc<0>,boost::shared_ptr<fecontext_type> > > > t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
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



    const uint16_type ndofv = functionspace_type::fe_type::nDof;
    LIFE_ASSERT( __v.size() == _M_functionspace->dof()->nDof() )( __v.size() )( _M_functionspace->dof()->nDof() ).warn( "invalid size" );
    __v.resize( _M_functionspace->dof()->nDof() );
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = _M_range;

    gm_context_ptrtype __c( new gm_context_type( __v.functionSpace()->gm(),
                                                 *it,
                                                 __geopc ) );
    typedef typename t_expr_type::shape shape;
    static const bool is_rank_ok = ( shape::M == FunctionSpaceType::nComponents1 &&
                                     shape::N == FunctionSpaceType::nComponents2 );

    BOOST_MPL_ASSERT_MSG( mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          (mpl::int_<shape::M>, mpl::int_<FunctionSpaceType::nComponents>, shape ) );

    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
    t_expr_type tensor_expr( _M_expr, mapgmc );

    std::vector<bool> points_done( __v.functionSpace()->dof()->nLocalDof()/__v.nComponents );
    std::fill( points_done.begin(), points_done.end(),false );
    for ( ; it!=en ; ++it )
    {
        __c->update( *it, __geopc );
        map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
        tensor_expr.update( mapgmc );

        for ( uint16_type __j = 0; __j < ndofv;++__j )
            {
                for ( uint16_type c1 = 0; c1 < shape::M;++c1 )
                    //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        __v.assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                    }
            }
    }
    return __v;
}

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( mpl::size_t<MESH_FACES> ) const
{
    boost::timer __timer;

    element_type __v( _M_functionspace );
    __v.zero();

#if 0
    if ( !boost::is_same<Elem1, typename Elem::functionspace_type>::value ||
         !boost::is_same<Elem2, typename Elem::functionspace_type>::value )
        return;
#endif
    Debug( 5066 ) << "call project() " << "\n";
    //
    // a few typedefs
    //

    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename geoelement_type::face_type face_type;

    // geometric mapping context
    typedef typename element_type::functionspace_type::mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gmc_ptrtype> > map_gmc_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    Debug(5066)  << "assembling Dirichlet conditions\n";

    dof_type const* __dof = __v.functionSpace()->dof().get();

    fe_type const* __fe = __v.functionSpace()->fe().get();

    iterator_type __face_it, __face_en;
    boost::tie( boost::tuples::ignore, __face_it, __face_en ) = _M_range;
    if ( __face_it == __face_en )
        return __v;

    gm_ptrtype __gm( new gm_type );



    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
    typedef typename gm_type::precompute_type geopc_type;
    Debug(5066)  << "[integratoron] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
    std::vector<geopc_ptrtype> __geopc( geoelement_type::numTopologicalFaces );
    for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
        {
            __geopc[__f] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
            //Debug(5066) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
            LIFE_ASSERT( __geopc[__f]->nPoints() ).error( "invalid number of points" );
        }

    uint16_type __face_id = __face_it->pos_first();
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element(0), __geopc[__face_id], __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
    t_expr_type expr( _M_expr, mapgmc );




    size_type nbFaceDof = invalid_size_type_value;
    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );
    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    Debug(5066)  << "nbFaceDof = " << nbFaceDof << "\n";
    //const size_type nbFaceDof = __fe->boundaryFE()->points().size2();

    std::vector<int> dofs;
    std::vector<value_type> values;

#if 0
    typedef typename functionspace_type::fe_type fe_type;
    fe_type* __fe = __v.functionSpace()->fe().get();

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type element_type;

    typedef typename functionspace_type::gm_type gm_type;
    typedef typename gm_type::template Context<vm::POINT, element_type> gm_context_type;

    typedef typename fe_type::template Context<vm::POINT, fe_type, gm_type, element_type> fecontext_type;

    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef fusion::map<fusion::pair<detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef typename expression_type::template tensor<map_gmc_type,fusion::map<fusion::pair<detail::gmc<0>,boost::shared_ptr<fecontext_type> > > > t_expr_type;
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



    const uint16_type ndofv = functionspace_type::fe_type::nDof;
    LIFE_ASSERT( __v.size() == _M_functionspace->dof()->nDof() )( __v.size() )( _M_functionspace->dof()->nDof() ).warn( "invalid size" );
    __v.resize( _M_functionspace->dof()->nDof() );
    iterator_type it, en;
    boost::tie( boost::tuples::ignore, it, en ) = _M_range;

    gm_context_ptrtype __c( new gm_context_type( __v.functionSpace()->gm(),
                                                 *it,
                                                 __geopc ) );
    typedef typename t_expr_type::shape shape;
    static const bool is_rank_ok = ( shape::M == FunctionSpaceType::nComponents1 &&
                                     shape::N == FunctionSpaceType::nComponents2 );

    BOOST_MPL_ASSERT_MSG( mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          (mpl::int_<shape::M>, mpl::int_<FunctionSpaceType::nComponents>, shape ) );
#endif
    //element_iterator __face_it = this->beginElement();
    for ( ;
          __face_it != __face_en;
          ++__face_it )
        {
            LIFE_ASSERT( __face_it->isOnBoundary() && !__face_it->isConnectedTo1() )
                ( __face_it->marker() )
                ( __face_it->isOnBoundary() )
                ( __face_it->ad_first() )
                ( __face_it->pos_first() )
                ( __face_it->ad_second() )
                ( __face_it->pos_second() )
                ( __face_it->id() ).warn( "inconsistent data face" );
            Debug(5066) << "[projector] FACE_ID = " << __face_it->id()
                       << " element id= " << __face_it->ad_first()
                       << " pos in elt= " << __face_it->pos_first()
                       << " marker: " << __face_it->marker() << "\n";
            Debug(5066) << "[projector] FACE_ID = " << __face_it->id() << " real pts=" << __face_it->G() << "\n";

            uint16_type __face_id = __face_it->pos_first();
            __c->update( __face_it->element(0), __geopc[__face_id], __face_id );

            Debug(5066) << "[projector] FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
            Debug(5066) << "[projector] FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );

            expr.update( mapgmc );

            std::pair<size_type,size_type> range_dof( std::make_pair( __v.start(),
                                                                      __v.functionSpace()->nDof()) );
            Debug(5066)  << "[projector] dof start = " << range_dof.first << "\n";
            Debug(5066)  << "[projector] dof range = " << range_dof.second << "\n";

            for ( uint16_type c1 = 0; c1 < shape::M;++c1 )
                for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        for ( uint16_type l = 0; l < nbFaceDof; ++l )
                            {
                                typename expression_type::value_type __value = expr.evalq( c1, c2, l );

                                size_type thedof =  __v.start() +
                                    boost::get<0>(__dof->faceLocalToGlobal( __face_it->id(), l, c1 ));
                                __v( thedof ) =  __value;
                            }
                    }
        } // face_it
#if 0
    for ( ; it!=en ; ++it )
    {
        __c->update( *it, __geopc );

        map_gmc_type mapgmc( fusion::make_pair<detail::gmc<0> >( __c ) );
        t_expr_type tensor_expr( _M_expr, mapgmc );

        for ( uint16_type __j = 0; __j < ndofv;++__j )
            {
                for ( uint16_type c1 = 0; c1 < shape::M;++c1 )
                    //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                        {
                            __v.assign( it->id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                        }
            }
    }
#endif
    return __v;
}

}
/// \endcond

/**
 * \brief nodal projection of \p __expr onto the space \p __functionspace
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr
 */
template<typename FunctionSpaceType, typename ExprT>
typename FunctionSpaceType::element_type
project( boost::shared_ptr<FunctionSpaceType> const& __functionspace, Expr<ExprT> const& __expr )
{
    typedef __typeof__( __functionspace->mesh()->elementsRange()) IteratorRange;
    typedef details::Projector<PROJ_NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, __functionspace->mesh()->elementsRange(), __expr );
    return p();
}

/**
 * \brief nodal projection of \p __expr onto the the subspace of \p __functionspace described by the range \p range_it
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr over the range \p __range_it
 */
template<typename FunctionSpaceType, typename IteratorRange, typename ExprT>
typename FunctionSpaceType::element_type
project( boost::shared_ptr<FunctionSpaceType> const& __functionspace,
         IteratorRange const& range_it,
         Expr<ExprT> const& __expr )
{
    typedef details::Projector<PROJ_NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr );
    return p();
}

/**
 * \brief nodal projection of \p __expr onto the space \p __functionspace
 *
 * \return the element of the space \p __functionspace resulting from the nodal projection of \p __expr
 */
template<typename FunctionSpaceType, typename ExprT>
typename FunctionSpaceType::element_type
project( FunctionSpaceType const& __functionspace, Expr<ExprT> const& __expr )
{
    typedef __typeof__( __functionspace->mesh()->elementsRange()) IteratorRange;
    typedef details::Projector<PROJ_NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, __functionspace->mesh()->elementsRange(), __expr );
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
         Expr<ExprT> const& __expr )
{
    typedef details::Projector<PROJ_NODAL, FunctionSpaceType, IteratorRange, Expr<ExprT> > proj_t;
    proj_t p( __functionspace, range_it, __expr );
    return p();
}

} // vf
} // life


#endif /* __Projectors_H */

