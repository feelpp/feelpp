/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-15

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Universite Joseph Fourier (Grenoble I)

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
   \file integratoron.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-15
 */
#ifndef __INTEGRATORON_HPP
#define __INTEGRATORON_HPP 1

#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <feel/feelalg/enums.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
template<typename T>
struct access_value
{
};
#if defined( FEELPP_HAS_QD_REAL )
template<>
struct access_value<dd_real>
{
    access_value( dd_real val, int /*n*/ )
    {
        v = val;
    }
    dd_real operator()() const
    {
        return v;
    }
    dd_real v;
};
template<>
struct access_value<qd_real>
{
    access_value( qd_real val, int /*n*/ )
    {
        v = val;
    }
    qd_real operator()() const
    {
        return v;
    }
    qd_real v;
};
#endif /*FEELPP_HAS_QD_REAL*/

#if defined(FEELPP_HAS_MPFR)
template<>
struct access_value<mp_type>
{
    access_value( mp_type val, int /*n*/ )
    {
        v = val;
    }
    mp_type operator()() const
    {
        return v;
    }
    mp_type v;
};
#endif /* FEELPP_HAS_MPFR */

template<>
struct access_value<double>
{
    access_value( double val, int /*n*/ )
    {
        v = val;
    }
    double operator()() const
    {
        return v;
    }
    double v;
};
template<>
struct access_value<int>
{
    access_value( int val, int /*n*/ )
    {
        v = val;
    }
    double operator()() const
    {
        return v;
    }
    double v;
};
template<>
struct access_value<node_type>
{
    access_value( node_type vec, int n )
    {
        v = vec[n];
    }
    double operator()() const
    {
        return v;
    }
    double v;
};
/*!
  \class IntegratorOnExpr
  \brief Handle Dirichlet condition



  @author Christophe Prud'homme
  @see
*/
template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr >
class IntegratorOnExpr
{
public:


    /** @name Typedefs
     */
    //@{
    static const size_type context = OnExpr::context|vm::POINT;
    static const size_type is_terminal = false;

    static const uint16_type imorder = OnExpr::imorder;
    static const bool imIsPoly = OnExpr::imIsPoly;

    typedef typename boost::tuples::template element<1, ElementRange>::type element_iterator;

    typedef Elem element_type;
    typedef RhsElem rhs_element_type;
    typedef typename element_type::value_type value_type;
    typedef typename element_type::return_type return_type;
    typedef boost::function<return_type ( node_type const& )> bc_type;
    typedef OnExpr expression_type;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = true;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = boost::is_same<Func,typename element_type::functionspace_type::basis_type>::value;
    };

    static const uint16_type nComponents = element_type::nComponents;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    IntegratorOnExpr( ElementRange const& __elts,
                      element_type const& __u,
                      rhs_element_type const& __rhs,
                      expression_type const& __expr,
                      size_type __on )
        :
        _M_eltbegin( __elts.template get<1>() ),
        _M_eltend( __elts.template get<2>() ),
        _M_u( __u ),
        _M_rhs( __rhs ),
        _M_expr( __expr ),
        _M_on_strategy( __on )
    {
    }
    IntegratorOnExpr( IntegratorOnExpr const& ioe )
        :
        _M_eltbegin( ioe._M_eltbegin ),
        _M_eltend( ioe._M_eltend ),
        _M_u( ioe._M_u ),
        _M_rhs( ioe._M_rhs ),
        _M_expr( ioe._M_expr ),
        _M_on_strategy( ioe._M_on_strategy )
    {
    }

    ~IntegratorOnExpr() {}

    //@}

    /** @name Accessors
    */
    //@{


    /**
     * iterator that points at the beginning of the container that
     * holds the data that will apply the Dirichlet condition upon
     */
    element_iterator beginElement() const
    {
        return _M_eltbegin;
    }

    /**
     * iterator that points at the end of the container that
     * holds the data that will apply the Dirichlet condition upon
     */
    element_iterator endElement() const
    {
        return _M_eltend;
    }


    //@}
    /** @name  Methods
     */
    //@{

    /**
     * assembly routine for Dirichlet condition
     *
     */
    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __u,
                   boost::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        typedef typename Elem::functionspace_type functionspace_type;
        Debug( 5066 ) << "[IntegratorOn::assemble()] is_same: "
                      << mpl::bool_<boost::is_same<functionspace_type,Elem1>::value>::value << "\n";
        assemble( __u, __v, __f, mpl::bool_<boost::is_same<functionspace_type,Elem1>::value>() );
    }
    //@}
private:

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& /*__u*/,
                   boost::shared_ptr<Elem2> const& /*__v*/,
                   FormType& /*__f*/, mpl::bool_<false> ) const {}

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( boost::shared_ptr<Elem1> const& __u,
                   boost::shared_ptr<Elem2> const& __v,
                   FormType& __f, mpl::bool_<true> ) const;

private:

    element_iterator _M_eltbegin;
    element_iterator _M_eltend;

    element_type const& _M_u;
    mutable rhs_element_type _M_rhs;
    expression_type _M_expr;
    Context _M_on_strategy;
};

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::assemble( boost::shared_ptr<Elem1> const& /*__u*/,
        boost::shared_ptr<Elem2> const& /*__v*/,
        FormType& __form,
        mpl::bool_<true> ) const
{
#if 0

    if ( !boost::is_same<Elem1, typename Elem::functionspace_type>::value ||
            !boost::is_same<Elem2, typename Elem::functionspace_type>::value )
        return;

#endif
    Debug( 5066 ) << "call on::assemble() " << "\n";
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
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    Debug( 5066 )  << "assembling Dirichlet conditions\n";
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;

    element_iterator __face_it = this->beginElement();
    element_iterator __face_en = this->endElement();
    if ( __face_it != __face_en )
    {

    dof_type const* __dof = _M_u.functionSpace()->dof().get();

    fe_type const* __fe = _M_u.functionSpace()->fe().get();

    gm_ptrtype __gm( new gm_type );


    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    typedef typename geoelement_type::permutation_type permutation_type;
    typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
    typedef typename gm_type::precompute_type geopc_type;
    Debug( 5066 )  << "[integratoron] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
    std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );

    for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
    {
        for ( permutation_type __p( permutation_type::IDENTITY );
                __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
        {
            __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
            //Debug(5066) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
            FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points" );
        }
    }

    uint16_type __face_id = __face_it->pos_first();
    gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    //t_expr_type expr( _M_expr, mapgmc );


    Debug( 5066 )  << "face_type::numVertices = " << face_type::numVertices << ", fe_type::nDofPerVertex = " << fe_type::nDofPerVertex << "\n"
                   << "face_type::numEdges = " << face_type::numEdges << ", fe_type::nDofPerEdge = " << fe_type::nDofPerEdge << "\n"
                   << "face_type::numFaces = " << face_type::numFaces << ", fe_type::nDofPerFace = " << fe_type::nDofPerFace << "\n";

    size_type nbFaceDof = invalid_size_type_value;

    if ( !fe_type::is_modal )
        nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                      face_type::numEdges * fe_type::nDofPerEdge +
                      face_type::numFaces * fe_type::nDofPerFace );

    else
        nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

    Debug( 5066 )  << "nbFaceDof = " << nbFaceDof << "\n";
    //const size_type nbFaceDof = __fe->boundaryFE()->points().size2();

    for ( ;
            __face_it != this->endElement();
            ++__face_it )
    {
        FEELPP_ASSERT( __face_it->isOnBoundary() && !__face_it->isConnectedTo1() )
        ( __face_it->marker() )
        ( __face_it->isOnBoundary() )
        ( __face_it->ad_first() )
        ( __face_it->pos_first() )
        ( __face_it->ad_second() )
        ( __face_it->pos_second() )
        ( __face_it->id() ).warn( "inconsistent data face" );
        Debug( 5066 ) << "FACE_ID = " << __face_it->id()
                      << " element id= " << __face_it->ad_first()
                      << " pos in elt= " << __face_it->pos_first()
                      << " marker: " << __face_it->marker() << "\n";
        Debug( 5066 ) << "FACE_ID = " << __face_it->id() << " face pts=" << __face_it->G() << "\n";

        uint16_type __face_id = __face_it->pos_first();
        __c->update( __face_it->element( 0 ), __face_id );

        Debug( 5066 ) << "FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
        Debug( 5066 ) << "FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

        t_expr_type expr( _M_expr, mapgmc );
        expr.update( mapgmc );

        std::pair<size_type,size_type> range_dof( std::make_pair( _M_u.start(),
                _M_u.functionSpace()->nDof() ) );
        Debug( 5066 )  << "[integratoron] dof start = " << range_dof.first << "\n";
        Debug( 5066 )  << "[integratoron] dof range = " << range_dof.second << "\n";

        for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
            for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
            {
                for ( uint16_type l = 0; l < nbFaceDof; ++l )
                {
                    Debug( 5066 ) << "[integratoronexpr] local dof=" << l
                                  << " |comp1=" << c1 << " comp 2= " << c2 << " | pt = " <<  __c->xReal( l ) << "\n";
                    typename expression_type::value_type __value = expr.evalq( c1, c2, l );
                    Debug( 5066 ) << "[integratoronexpr] value=" << __value << "\n";

                    // global Dof
                    size_type thedof =  _M_u.start() +
                                        boost::get<0>( __dof->faceLocalToGlobal( __face_it->id(), l, c1 ) );

                    //size_type thedof_nproc = __dof->dofNProc( thedof );
                    if ( std::find( dofs.begin(),
                                    dofs.end(),
                                    thedof ) != dofs.end() )
                        continue;

                    if ( _M_on_strategy.test( ON_ELIMINATION ) )
                    {
                        Debug( 5066 ) << "Eliminating row " << thedof << " using value : " << __value << "\n";

                        // this can be quite expensive depending on the
                        // matrix storage format.
                        //__form.diagonalize( thedof, range_dof, _M_rhs, __value, thedof_nproc );
#if !defined(FEELPP_ENABLE_MPI_MODE)
                        dofs.push_back( thedof );
                        values.push_back( __value );
#else

                        // only the real dof ( not the ghosts )
                        if ( __form.testSpace()->mapOn().dofGlobalClusterIsOnProc( __form.testSpace()->mapOn().mapGlobalProcessToGlobalCluster( thedof ) ) )
                        {
                            dofs.push_back( thedof );
                            values.push_back( __value );
                        }

#endif


                        //_M_rhs.set( thedof, __value );
                    }

                    else if (  _M_on_strategy.test( ON_PENALISATION ) &&
                               !_M_on_strategy.test( ON_ELIMINATION ) )
                    {
                        __form.set( thedof, thedof, 1.0*1e30 );
                        _M_rhs->set( thedof, __value*1e30 );
                    }
                } // loop on space components

            } // loop on face dof
    }

    } // __face_it != __face_en

    __form.zeroRows( dofs, values, *_M_rhs, _M_on_strategy );
}

#if 1

namespace detail
{
template<typename Args>
struct integratoron_type
{
    typedef typename clean_type<Args,tag::range>::type _range_type;
    typedef typename clean_type<Args,tag::rhs>::type _rhs_type;
    typedef typename clean_type<Args,tag::element>::type _element_type;
    typedef typename clean_type<Args,tag::expr>::type _expr_type;


    typedef IntegratorOnExpr<_range_type, _element_type, _rhs_type,
            typename mpl::if_<boost::is_arithmetic<_expr_type>,
            mpl::identity<Expr<Cst<_expr_type> > >,
            mpl::identity<_expr_type> >::type::type> type;
    typedef Expr<type> expr_type;
};

}
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
    ( typename vf::detail::integratoron_type<Args>::expr_type ), // return type
    on,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( element, *  )
      ( rhs, *  )
      ( expr,   * )
    ) // 4. one required parameter, and

    ( optional
      ( type,   ( size_type ), ON_ELIMINATION|ON_ELIMINATION_KEEP_DIAGONAL )
      ( verbose,   ( bool ), false )
    )
)
{
    typename vf::detail::integratoron_type<Args>::type ion( range, element, rhs, expr, type );
    return typename vf::detail::integratoron_type<Args>::expr_type( ion );
}

#else
namespace detail
{

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
Expr<IntegratorOnExpr<ElementRange, Elem,RhsElem,
     typename mpl::if_<boost::is_arithmetic<OnExpr>,
     mpl::identity<Expr<Cst<OnExpr> > >,
     mpl::identity<OnExpr> >::type::type> >
     on( ElementRange const& __r,
         Elem const& __u,
         RhsElem&  __rhs,
         OnExpr const& __e,
         size_type __on,
         mpl::bool_<false> )
{
    typedef typename mpl::if_<boost::is_arithmetic<OnExpr>,
            mpl::identity<Expr<Cst<OnExpr> > >,
            mpl::identity<OnExpr> >::type::type expr_type;
    expr_type expr( __e );
    typedef IntegratorOnExpr<ElementRange, Elem,RhsElem, expr_type> expr_t;
    return Expr<expr_t>( expr_t( __r, __u, __rhs, expr, __on ) );
    //return on( __r, __u, __rhs, expr, __on, mpl::bool_<boost::is_arithmetic<OnExpr>::value>() );
}
template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
Expr<IntegratorOnExpr<ElementRange, Elem, typename RhsElem::value_type,
     typename mpl::if_<boost::is_arithmetic<OnExpr>,
     mpl::identity<Expr<Cst<OnExpr> > >,
     mpl::identity<OnExpr> >::type::type> >
     on( ElementRange const& __r,
         Elem const& __u,
         RhsElem&  __rhs,
         OnExpr const& __e,
         size_type __on,
         mpl::bool_<true> )
{
    return on( __r, __u, __rhs, __e, __on, mpl::bool_<false>() );
}

} // detail namespace
/// \endcond

/**
 * Apply Dirichlet boundary condition over of a set of \p Dof.
 *
 * \param __r tuple of iterators over geometric entities (faces)
 * \param __u element of a functional space whose \p Dof will be constrained
 * \param __rhs right hand side of the system which shall be modified
 * \param __e expression which gives the value of the constraint
 * \param  __on type of strategy to impose the constraint
 *
 */
template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
Expr<IntegratorOnExpr<ElementRange, Elem,
     typename mpl::if_<mpl::or_<is_shared_ptr<RhsElem>, boost::is_pointer<RhsElem> >,
     mpl::identity<typename RhsElem::value_type>,
     mpl::identity<RhsElem> >::type::type,
     typename mpl::if_<boost::is_arithmetic<OnExpr>,
     mpl::identity<Expr<Cst<OnExpr> > >,
     mpl::identity<OnExpr> >::type::type> >
     on( ElementRange const& __r,
         Elem const& __u,
         RhsElem&  __rhs,
         OnExpr const& __e,
         size_type __on = ON_ELIMINATION|ON_ELIMINATION_KEEP_DIAGONAL )
{
    returnvf::detail::on( __r, __u, __rhs, __e, __on,  mpl::or_<is_shared_ptr<RhsElem>, boost::is_pointer<RhsElem> >() );

}

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
Expr<IntegratorOnExpr<ElementRange, Elem,
     typename mpl::if_<mpl::or_<is_shared_ptr<RhsElem>, boost::is_pointer<RhsElem> >,
     mpl::identity<typename RhsElem::value_type>,
     mpl::identity<RhsElem> >::type::type,
     typename mpl::if_<boost::is_arithmetic<OnExpr>,
     mpl::identity<Expr<Cst<OnExpr> > >,
     mpl::identity<OnExpr> >::type::type> >
     on( ElementRange const& __r,
         boost::shared_ptr<Elem>  __u,
         boost::shared_ptr<RhsElem>  __rhs,
         OnExpr const& __e,
         size_type __on = ON_ELIMINATION|ON_ELIMINATION_KEEP_DIAGONAL )
{
    returnvf::detail::on( __r, __u, __rhs, __e, __on,  mpl::or_<is_shared_ptr<RhsElem>, boost::is_pointer<RhsElem> >() );

}
#endif // 0

} // vf
} // feel
#endif
