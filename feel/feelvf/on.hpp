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
#ifndef FEELPP_INTEGRATORON_HPP
#define FEELPP_INTEGRATORON_HPP 1

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
    typedef value_type evaluate_type;
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
        M_elts(),
        M_eltbegin( __elts.template get<1>() ),
        M_eltend( __elts.template get<2>() ),
        M_u( __u ),
        M_rhs( __rhs ),
        M_expr( __expr ),
        M_on_strategy( __on )
    {
        M_elts.push_back( __elts );
    }
    IntegratorOnExpr( std::list<ElementRange> const& __elts,
                      element_type const& __u,
                      rhs_element_type const& __rhs,
                      expression_type const& __expr,
                      size_type __on )
        :
        M_elts( __elts ),
        M_u( __u ),
        M_rhs( __rhs ),
        M_expr( __expr ),
        M_on_strategy( __on )
    {
        if ( __elts.size() )
        {
            M_eltbegin = __elts.begin()->template get<1>();
            M_eltend = __elts.begin()->template get<2>();
        }
    }
    IntegratorOnExpr( IntegratorOnExpr const& ioe )
        :
        M_elts( ioe.M_elts ),
        M_eltbegin( ioe.M_eltbegin ),
        M_eltend( ioe.M_eltend ),
        M_u( ioe.M_u ),
        M_rhs( ioe.M_rhs ),
        M_expr( ioe.M_expr ),
        M_on_strategy( ioe.M_on_strategy )
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
        DVLOG(2) << "[IntegratorOn::assemble()] is_same: "
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

    std::list<ElementRange> M_elts;
    element_iterator M_eltbegin;
    element_iterator M_eltend;

    element_type const& M_u;
    mutable rhs_element_type M_rhs;
    expression_type M_expr;
    Context M_on_strategy;
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
    DVLOG(2) << "call on::assemble() " << "\n";
    //
    // a few typedefs
    //

    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename element_type::functionspace_type::mesh_type::shape_type geoshape_type;
    typedef typename geoelement_type::face_type face_type;

    typedef typename element_type::functionspace_type::fe_type fe_type;

    // geometric mapping context
    typedef typename element_type::functionspace_type::mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    //typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef typename mpl::if_< mpl::or_<is_hdiv_conforming<fe_type>, is_hcurl_conforming<fe_type> >,
                               typename gm_type::template Context<context|vm::JACOBIAN|vm::KB|vm::TANGENT|vm::NORMAL, geoelement_type>,
                               typename gm_type::template Context<context, geoelement_type> >::type gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

    static const uint16_type nDim = geoshape_type::nDim;

    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::shape shape;

    // make sure that the form is close, ie the associated matrix is assembled
    __form.matrix().close();
    // make sure that the right hand side is closed, ie the associated vector is assembled
    M_rhs->close();

    //
    // start
    //
    DVLOG(2)  << "assembling Dirichlet conditions\n";
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;
    element_iterator __face_it = this->beginElement();
    element_iterator __face_en = this->endElement();

    bool findAFace = false;
    for( auto& lit : M_elts )
    {
        __face_it = lit.template get<1>();
        __face_en = lit.template get<2>();
        if ( __face_it != __face_en )
        {
            findAFace=true;
            break;
        }
    }
    if ( findAFace )
    {
        // get the first face properly connected
        bool findAFaceToInit=false;
        for( auto& lit : M_elts )
        {
            __face_it = lit.template get<1>();
            __face_en = lit.template get<2>();
            for( ; __face_it != __face_en; ++__face_it )
            {
                if ( boost::unwrap_ref(*__face_it).isConnectedTo0() )
                {
                    findAFaceToInit=true;
                    break;
                }
            }
            if ( findAFaceToInit ) break;
        }
        CHECK( findAFaceToInit ) << "not find a face to init\n";

        auto const& faceForInit = boost::unwrap_ref( *__face_it );

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
        DVLOG(2)  << "[integratoron] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
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

        uint16_type __face_id = faceForInit.pos_first();
        gmc_ptrtype __c( new gmc_type( __gm, faceForInit.element( 0 ), __geopc, __face_id ) );

        map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
        t_expr_type expr( M_expr, mapgmc );


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

        auto IhLoc = __fe->faceLocalInterpolant();
        for( auto& lit : M_elts )
        {
        __face_it = lit.template get<1>();
        __face_en = lit.template get<2>();
        for ( ;
              __face_it != __face_en;//this->endElement();
              ++__face_it )
        {
            auto const& theface = boost::unwrap_ref( *__face_it );

            if ( !theface.isConnectedTo0() )
            {
                LOG( WARNING ) << "face not connected" << theface;

                continue;
            }
            // do not process the face if it is a ghost face: belonging to two
            // processes and being in a process id greater than the one
            // corresponding face
            if ( theface.isGhostFace() )
            {
                LOG(WARNING) << "face id : " << theface.id() << " is a ghost face";
                continue;
            }

            DVLOG(2) << "FACE_ID = " << theface.id()
                     << " element id= " << theface.ad_first()
                     << " pos in elt= " << theface.pos_first()
                     << " marker: " << theface.marker() << "\n";
            DVLOG(2) << "FACE_ID = " << theface.id() << " face pts=" << theface.G() << "\n";

            uint16_type __face_id = theface.pos_first();
            __c->update( theface.element( 0 ), __face_id );

            DVLOG(2) << "FACE_ID = " << theface.id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "FACE_ID = " << theface.id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );

            t_expr_type expr( M_expr, mapgmc );
            expr.update( mapgmc );

            std::pair<size_type,size_type> range_dof( std::make_pair( M_u.start(),
                                                                      M_u.functionSpace()->nDof() ) );
            DVLOG(2)  << "[integratoron] dof start = " << range_dof.first << "\n";
            DVLOG(2)  << "[integratoron] dof range = " << range_dof.second << "\n";

            //use interpolant
            __fe->faceInterpolate( expr, IhLoc );

            for( auto const& ldof : M_u.functionSpace()->dof()->faceLocalDof( theface.id() ) )
                {
                    size_type thedof = M_u.start()+ ldof.index(); // global dof
                    DCHECK( ldof.localDofInFace() < IhLoc.size() ) 
                        << "Invalid local dof index in face for face Interpolant "
                        << ldof.localDofInFace() << ">=" << IhLoc.size();
                    double __value = ldof.sign()*IhLoc( ldof.localDofInFace() );

                    if ( std::find( dofs.begin(),
                                    dofs.end(),
                                    thedof ) != dofs.end() )
                        continue;

                    if ( M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                        {
                            DVLOG(2) << "Eliminating row " << thedof << " using value : " << __value << "\n";

                            // this can be quite expensive depending on the
                            // matrix storage format.
                            //__form.diagonalize( thedof, range_dof, M_rhs, __value, thedof_nproc );

                            // only the real dof ( not the ghosts )
                            //if ( __form.testSpace()->mapOn().dofGlobalClusterIsOnProc( __form.testSpace()->mapOn().mapGlobalProcessToGlobalCluster( thedof ) ) )
                            {
                                dofs.push_back( thedof );
                                values.push_back(  __value );
                            }

                            //M_rhs.set( thedof, __value );
                        }

                    else if (  M_on_strategy.test( ContextOn::PENALISATION ) &&
                               !M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                        {
                            __form.set( thedof, thedof, 1.0*1e30 );
                            M_rhs->set( thedof, __value*1e30 );
                        }
                }
        }// __face_it != __face_en
        } // for( auto& lit : M_elts )
    }// findAFace

    if ( __form.rowStartInMatrix()!=0)
    {
        auto const thedofshift = __form.rowStartInMatrix();
        for (auto& itd : dofs)
            itd+=thedofshift;
    }

    auto x = M_rhs->clone();
    CHECK( values.size() == dofs.size() ) << "Invalid dofs/values size: " << dofs.size() << "/" << values.size();
    x->setVector( dofs.data(), dofs.size(), values.data() );
    x->close();
    __form.zeroRows( dofs, *x, *M_rhs, M_on_strategy );
    x.reset();
}


namespace detail
{
template<typename T >
struct v_ptr1
{
    typedef T type;
};
template<typename T >
struct v_ptr2
{
    typedef typename T::vector_ptrtype type;
};

template<typename Args>
struct integratoron_type
{
    typedef typename clean_type<Args,tag::range>::type _range_base_type;
    typedef typename clean_type<Args,tag::rhs>::type _rhs_type;
    typedef typename clean_type<Args,tag::element>::type _element_type;
    typedef typename clean_type<Args,tag::expr>::type _expr_type;

    typedef typename mpl::if_< boost::is_std_list<_range_base_type>,
                               mpl::identity<_range_base_type>,
                               mpl::identity<std::list<_range_base_type> > >::type::type::value_type _range_type;

#if 1
    typedef typename mpl::if_<Feel::detail::is_vector_ptr<_rhs_type>,
                              mpl::identity<v_ptr1<_rhs_type> >,
                              mpl::identity<v_ptr2<_rhs_type> > >::type::type::type the_rhs_type;
#else
    typedef _rhs_type the_rhs_type;
#endif
typedef IntegratorOnExpr<_range_type, _element_type, the_rhs_type,
            typename mpl::if_<boost::is_arithmetic<_expr_type>,
            mpl::identity<Expr<Cst<_expr_type> > >,
            mpl::identity<_expr_type> >::type::type> type;
    typedef Expr<type> expr_type;
};

template<typename V>
typename V::vector_ptrtype
getRhsVector( V const&  v, mpl::false_ )
{
    return v.vectorPtr();
}

template<typename V>
V
getRhsVector( V const&  v, mpl::true_ )
{
    return v;
}

template<typename V>
typename mpl::if_<Feel::detail::is_vector_ptr<V>,
                  mpl::identity<v_ptr1<V> >,
                  mpl::identity<v_ptr2<V> > >::type::type::type
getRhsVector( V const&  v )
{
    return getRhsVector( v, Feel::detail::is_vector_ptr<V>() );
}


}
///\endcond detail
/**
 *
 * \brief projection/interpolation of an expresion onto a nodal functionspace
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
      ( prefix,   ( std::string ), "" )
      ( type,   ( std::string ), soption(_prefix=prefix,_name="on.type") )
      ( verbose,   ( bool ), boption(_prefix=prefix,_name="on.verbose") )
        )
    )
{
    typename vf::detail::integratoron_type<Args>::type ion( range,
                                                            element,
                                                            Feel::vf::detail::getRhsVector(rhs),
                                                            expr,
                                                            size_type(ContextOnMap[type]) );
    if ( verbose )
    {
        LOG(INFO) << "Dirichlet condition over : "<< nelements(range) << " faces";
        switch( ContextOnMap[type] )
        {
        case ContextOn::ELIMINATION:
            LOG(INFO) << "treatment of Dirichlet condition: " << type << " (elimination, unsymmetric)";
            break;
        // case size_type(ContextOn::ELIMINATION)|size_type(ContextOn::ELIMINATION_KEEP_DIAGONAL):
        //     LOG(INFO) << "treatment of Dirichlet condition: " << type << " (elimination and keep diagonal, unsymmetric)";
        //     break;
        case ContextOn::SYMMETRIC:
            LOG(INFO) << "treatment of Dirichlet condition: " << type << " (elimination, symmetric, more expensive than unsymmetric treatment)";
            break;
        // case size_type(ContextOn::ELIMINATION_SYMMETRIC)|size_type(ContextOn::ELIMINATION_KEEP_DIAGONAL):
        //     LOG(INFO) << "treatment of Dirichlet condition: " << type << " (elimination and keep diagonal, symmetric, more expensive than unsymmetric treatment)";
        //     break;
        case ContextOn::PENALISATION:
            LOG(INFO) << "treatment of Dirichlet condition: " << type << " (penalisation, symmetric, very big value on diagonal)";
            break;
        default:
            break;
        }
    }
    //typename vf::detail::integratoron_type<Args>::type ion( range, element, rhs, expr, type );
    return typename vf::detail::integratoron_type<Args>::expr_type( ion );
}


} // vf
} // feel
#endif
