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

#include <feel/feelvf/detail/clean.hpp>
#include <feel/feelvf/expr.hpp>

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

    typedef IteratorRange range_iterator;
    typedef typename mpl::if_< boost::is_std_list<range_iterator>,
                               mpl::identity<range_iterator>,
                               mpl::identity<std::list<range_iterator> > >::type::type::value_type range_iterator_type;
    typedef typename boost::tuples::template element<0, range_iterator_type>::type idim_type;
    typedef typename boost::tuples::template element<1, range_iterator_type>::type iterator_type;


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
    __v.on( _range=M_range, _expr=M_expr, _geomap=M_geomap_strategy, _accumulate=sum );
    return __v;
}

template<ProjectorType iDim, typename FunctionSpaceType, typename Iterator, typename ExprT>
typename Projector<iDim, FunctionSpaceType, Iterator, ExprT>::element_type
Projector<iDim, FunctionSpaceType, Iterator, ExprT>::operator()( const bool sum, mpl::size_t<MESH_FACES> ) const
{
    boost::timer __timer;

    element_type __v( M_functionspace );
    __v.setZero();
    __v.on( _range=M_range, _expr=M_expr, _geomap=M_geomap_strategy, _accumulate=sum );
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
