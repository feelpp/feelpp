/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-03

  Copyright (C) 2013 Feel++ Consortium

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
   \file on.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-03
 */
#ifndef FEELPP_ELEMENTON_HPP
#define FEELPP_ELEMENTON_HPP 1

#include <boost/timer.hpp>
#include <boost/foreach.hpp>

#include <feel/feelalg/enums.hpp>

namespace Feel
{
namespace vf
{
///\cond detail
/*!
  \class ElementOnExpr
  \brief Handle Dirichlet condition

  @author Christophe Prud'homme
  @see
*/
template<typename ElementRange, typename Elem, typename OnExpr >
class ElementOnExpr
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

    ElementOnExpr( ElementRange const& __elts,
                   element_type& __u,
                   expression_type const& __expr )
        :
        M_eltbegin( __elts.template get<1>() ),
        M_eltend( __elts.template get<2>() ),
        M_u( __u ),
        M_expr( __expr )
    {
    }
    ElementOnExpr( ElementOnExpr const& ioe )
        :
        M_eltbegin( ioe.M_eltbegin ),
        M_eltend( ioe.M_eltend ),
        M_u( ioe.M_u ),
        M_expr( ioe.M_expr )
    {
    }

    ~ElementOnExpr() {}

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
        return M_eltbegin;
    }

    /**
     * iterator that points at the end of the container that
     * holds the data that will apply the Dirichlet condition upon
     */
    element_iterator endElement() const
    {
        return M_eltend;
    }


    //@}
    /** @name  Methods
     */
    //@{
    void apply();
    //@}

private:

    element_iterator M_eltbegin;
    element_iterator M_eltend;

    element_type& M_u;
    expression_type M_expr;
};

template<typename ElementRange, typename Elem, typename OnExpr>
void
ElementOnExpr<ElementRange, Elem, OnExpr>::apply()
{
    DVLOG(2) << "call on::assemble() " << "\n";
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
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;
    element_iterator __face_it = this->beginElement();
    element_iterator __face_en = this->endElement();
    if ( __face_it != __face_en )
    {
        // get the first face properly connected
        for( ; __face_it != __face_en; ++__face_it )
            if ( __face_it->isConnectedTo0() )
                break;


        dof_type const* __dof = M_u.functionSpace()->dof().get();

        fe_type const* __fe = M_u.functionSpace()->fe().get();

        gm_ptrtype __gm( new gm_type );


        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        typedef typename geoelement_type::permutation_type permutation_type;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        DVLOG(2)  << "[elementon] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
        std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );

        for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
                //DVLOG(2) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
                FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points" );
            }
        }

        uint16_type __face_id = __face_it->pos_first();
        gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
        //t_expr_type expr( M_expr, mapgmc );


        DVLOG(2)  << "face_type::numVertices = " << face_type::numVertices << ", fe_type::nDofPerVertex = " << fe_type::nDofPerVertex << "\n"
                  << "face_type::numEdges = " << face_type::numEdges << ", fe_type::nDofPerEdge = " << fe_type::nDofPerEdge << "\n"
                  << "face_type::numFaces = " << face_type::numFaces << ", fe_type::nDofPerFace = " << fe_type::nDofPerFace << "\n";

        size_type nbFaceDof = invalid_size_type_value;

        if ( !fe_type::is_modal )
            nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                          face_type::numEdges * fe_type::nDofPerEdge +
                          face_type::numFaces * fe_type::nDofPerFace );

        else
            nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

        DVLOG(2)  << "nbFaceDof = " << nbFaceDof << "\n";
        //const size_type nbFaceDof = __fe->boundaryFE()->points().size2();

        for ( ;
              __face_it != this->endElement();
              ++__face_it )
        {
            if ( !__face_it->isConnectedTo0() )
            {
                LOG( WARNING ) << "face not connected" << *__face_it;

                continue;
            }
            // do not process the face if it is a ghost face: belonging to two
            // processes and being in a process id greater than the one
            // corresponding face
            if ( __face_it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << __face_it->id() << " is a ghost face";
                continue;
            }

            DVLOG(2) << "FACE_ID = " << __face_it->id()
                     << " element id= " << __face_it->ad_first()
                     << " pos in elt= " << __face_it->pos_first()
                     << " marker: " << __face_it->marker() << "\n";
            DVLOG(2) << "FACE_ID = " << __face_it->id() << " face pts=" << __face_it->G() << "\n";

            uint16_type __face_id = __face_it->pos_first();
            __c->update( __face_it->element( 0 ), __face_id );

            DVLOG(2) << "FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            t_expr_type expr( M_expr, mapgmc );
            expr.update( mapgmc );

            std::pair<size_type,size_type> range_dof( std::make_pair( M_u.start(),
                                                                      M_u.functionSpace()->nDof() ) );
            DVLOG(2)  << "[elementon] dof start = " << range_dof.first << "\n";
            DVLOG(2)  << "[elementon] dof range = " << range_dof.second << "\n";

            for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                {
                    for ( uint16_type l = 0; l < nbFaceDof; ++l )
                    {
                        DVLOG(2) << "[elementonexpr] local dof=" << l
                                 << " |comp1=" << c1 << " comp 2= " << c2 << " | pt = " <<  __c->xReal( l ) << "\n";
                        typename expression_type::value_type __value = expr.evalq( c1, c2, l );
                        DVLOG(2) << "[elementonexpr] value=" << __value << "\n";

                        // global Dof
                        size_type thedof =  M_u.start() +
                            boost::get<0>( __dof->faceLocalToGlobal( __face_it->id(), l, c1 ) );

                        //size_type thedof_nproc = __dof->dofNProc( thedof );
                        if ( std::find( dofs.begin(),
                                        dofs.end(),
                                        thedof ) != dofs.end() )
                            continue;

                        M_u( thedof ) = __value;
                    } // loop on space components

                } // loop on face dof
        }

    } // __face_it != __face_en
}


namespace detail
{

template<typename Args>
struct elementon_type
{
    typedef typename clean_type<Args,tag::range>::type _range_type;
    typedef typename clean_type<Args,tag::element>::type _element_type;
    typedef typename clean_type<Args,tag::expr>::type _expr_type;
    typedef ElementOnExpr<_range_type, _element_type,
                      typename mpl::if_<boost::is_arithmetic<_expr_type>,
                                        mpl::identity<Expr<Cst<_expr_type> > >,
                                        mpl::identity<_expr_type> >::type::type> type;
};



}
/**
 *
 * \brief projection/interpolation of an expresion onto a nodal functionspace
 *
 * \arg space the function space to project onto
 * \arg range the range of mesh elements to apply the projection (the remaining parts are set to 0)
 * \arg expr the expression to project
 * \arg geomap the type of geomap to use (make sense only using high order meshes)
 */
BOOST_PARAMETER_FUNCTION(
    ( typename vf::detail::elementon_type<Args>::type ), // return type
    on2,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( element, *  )
      ( expr,   * )
        ) // 4. one required parameter, and

    ( optional
      ( prefix,   ( std::string ), "" )
      ( verbose,   ( bool ), boption(_prefix=prefix,_name="on.verbose") )
        )
    )
{
    typename vf::detail::elementon_type<Args>::type ion( range,element,expr );
    if ( verbose )
    {
        LOG(INFO) << "set Dof over : "<< nelements(range) << " faces";
    }
    return ion;
}


} // vf
} // feel
#endif
