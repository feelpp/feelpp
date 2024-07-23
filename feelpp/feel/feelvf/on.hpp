/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-15

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
#ifndef FEELPP_VF_ON_H
#define FEELPP_VF_ON_H

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/vectorcondensed.hpp>
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

    //static inline const uint16_type imorder = OnExpr::imorder;
    //static const bool imIsPoly = OnExpr::imIsPoly;

    using on_type =  typename ElementRange::idim_t;
    using element_iterator =  typename ElementRange::iterator_t;

    typedef Elem element_type;
    typedef RhsElem rhs_element_type;
    typedef typename element_type::value_type value_type;
    typedef typename element_type::return_type return_type;
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

    template<typename Func>
    static const bool has_test_basis = false;
    template<typename Func>
    static const bool has_trial_basis = boost::is_same<Func,typename element_type::functionspace_type::basis_type>::value;
    using test_basis = std::nullptr_t;
    using trial_basis = typename element_type::functionspace_type::basis_type;


    static inline const uint16_type nComponents = element_type::nComponents;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    IntegratorOnExpr() = delete;

    IntegratorOnExpr( ElementRange const& __elts,
                      element_type const& __u,
                      rhs_element_type const& __rhs,
                      expression_type const& __expr,
                      size_type __on,
                      double value_on_diag )
        :
        M_elts(),
        M_eltbegin( __elts.begin() ),
        M_eltend( __elts.end() ),
        M_u( __u ),
        M_rhs( __rhs ),
        M_expr( __expr ),
        M_on_strategy( __on ),
        M_value_on_diagonal( value_on_diag )
    {
        M_elts.push_back( __elts );
    }

    IntegratorOnExpr( ElementRange const& __elts,
                      element_type && __u,
                      rhs_element_type const& __rhs,
                      expression_type const& __expr,
                      size_type __on,
                      double value_on_diag )
        :
        M_elts(),
        M_eltbegin( __elts.begin() ),
        M_eltend( __elts.end() ),
        M_uFromRValue( std::make_shared<element_type>( std::forward<element_type>(__u) ) ),
        M_u( *M_uFromRValue ),
        M_rhs( __rhs ),
        M_expr( __expr ),
        M_on_strategy( __on ),
        M_value_on_diagonal( value_on_diag )
    {
        M_elts.push_back( __elts );
    }

    IntegratorOnExpr( std::list<ElementRange> const& __elts,
                      element_type const& __u,
                      rhs_element_type const& __rhs,
                      expression_type const& __expr,
                      size_type __on,
                      double value_on_diag )
        :
        M_elts( __elts ),
        M_u( __u ),
        M_rhs( __rhs ),
        M_expr( __expr ),
        M_on_strategy( __on ),
        M_value_on_diagonal( value_on_diag )
    {
        if ( __elts.size() )
        {
            M_eltbegin = __elts.begin()->begin();
            M_eltend = __elts.begin()->end();
        }
    }
    IntegratorOnExpr( IntegratorOnExpr const& ioe ) = default;
    IntegratorOnExpr( IntegratorOnExpr && ioe ) = default;
    ~IntegratorOnExpr() = default;

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
    void assemble( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const
    {
        typedef typename Elem::functionspace_type functionspace_type;
        static constexpr bool is_same_space = boost::is_same<functionspace_type,Elem1>::value;
        static constexpr bool is_comp_space = boost::is_same<functionspace_type,typename Elem1::component_functionspace_type>::value;
        VLOG(2) << fmt::format("[IntegratorOn::assemble()] is_same: {} is_comp: {}", is_same_space, is_comp_space);
        if constexpr ( ( is_same_space || is_comp_space ) && ( on_type::value == MESH_ELEMENTS ) )
            onElements( __u, __v, __f );
        if constexpr ( ( is_same_space || is_comp_space ) && ( on_type::value == MESH_FACES ) )
            onFaces( __u, __v, __f );
        if constexpr ( ( is_same_space || is_comp_space ) && ( on_type::value == MESH_EDGES ) )
            onEdges( __u, __v, __f );
        if constexpr ( ( is_same_space || is_comp_space ) && ( on_type::value == MESH_POINTS ) )
            onPoints( __u, __v, __f );
    }
    //@}
private:

    template<typename Elem1, typename Elem2, typename FormType>
    void assemble( std::shared_ptr<Elem1> const& /*__u*/,
                   std::shared_ptr<Elem2> const& /*__v*/,
                   FormType& /*__f*/, mpl::bool_<false>, on_type ) const {}

    template<typename Elem1, typename Elem2, typename FormType>
    void onElements( std::shared_ptr<Elem1> const& __u,
                     std::shared_ptr<Elem2> const& __v,
                     FormType& __f ) const;

    template<typename Elem1, typename Elem2, typename FormType>
    void onFaces( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const;

    template<typename Elem1, typename Elem2, typename FormType>
    void onEdges( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const;

    template<typename Elem1, typename Elem2, typename FormType>
    void onPoints( std::shared_ptr<Elem1> const& __u,
                   std::shared_ptr<Elem2> const& __v,
                   FormType& __f ) const;

private:

    std::list<ElementRange> M_elts;
    element_iterator M_eltbegin;
    element_iterator M_eltend;

    std::shared_ptr<element_type> M_uFromRValue;
    element_type const& M_u;
    mutable rhs_element_type M_rhs;
    expression_type M_expr;
    Context M_on_strategy;
    double M_value_on_diagonal { 1.0 };
};

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onElements( std::shared_ptr<Elem1> const& /*__u*/,
                                                                    std::shared_ptr<Elem2> const& /*__v*/,
                                                                    FormType& __form ) const
{
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onElements()");
    
    typedef typename Elem::functionspace_type functionspace_type;
    static constexpr bool is_same_space = boost::is_same<functionspace_type,Elem1>::value;
    static constexpr bool is_comp_space = Elem1::is_vectorial && Elem1::is_product && boost::is_same<functionspace_type,typename Elem1::component_functionspace_type>::value;

    if (  M_on_strategy.test( ContextOn::PENALISATION ) )
    {
        // make sure that the form is close, ie the associated matrix is assembled
        __form.matrix().close();
        // make sure that the right hand side is closed, ie the associated vector is assembled
        M_rhs->close();
    }

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto elt_it = this->beginElement();
    auto elt_en = this->endElement();

    bool findAElt = false;
    for( auto& lit : M_elts )
    {
        elt_it = lit.begin();
        elt_en = lit.end();
        if ( elt_it != elt_en )
        {
            findAElt=true;
            break;
        }
    }
    if ( findAElt )
    {
        auto const& eltForInit = boost::unwrap_ref( *elt_it );

        auto const dof = M_u.functionSpace()->dof();
        auto const fe = M_u.functionSpace()->fe();
        auto const mesh = M_u.functionSpace()->mesh();
        auto gm = mesh->gm();
        auto geopc = gm->preCompute( fe->points() );
        auto ctx = gm->template context<context>( eltForInit, geopc );
        auto expr_evaluator = M_expr.evaluator( mapgmc(ctx) );
        auto IhLoc = fe->localInterpolant();

        int compDofShift = (is_comp_space)? ((int)M_u.component()) : 0;
        auto const& trialDofIdToContainerId = __form.dofIdToContainerIdTrial();

        for( auto& lit : M_elts )
        {
            elt_it = lit.begin();
            elt_en = lit.end();
            DVLOG(2) << "element nb: " << std::distance(elt_it,elt_en);
            for ( ; elt_it != elt_en; ++elt_it )
            {
                auto const& curElt = unwrap_ref( *elt_it );
                ctx->template update<context>( curElt, geopc );
                expr_evaluator.update( mapgmc( ctx ) );
                fe->interpolate( expr_evaluator, IhLoc );

                auto const& s = dof->localToGlobalSigns( curElt.id() );
                for( auto const& ldof : dof->localDof( curElt.id() ) )
                {
                    //size_type index = ldof.second.index();
                    index_type thedof = (is_comp_space)? compDofShift+Elem1::nComponents*ldof.second.index() : ldof.second.index();
                    thedof = trialDofIdToContainerId[ thedof ];
                    if ( std::find( dofs.begin(),dofs.end(),thedof ) != dofs.end() )
                        continue;

                    uint16_type thelocdof = ldof.first.localDof();
                    double __value = s(thelocdof)*IhLoc( thelocdof );

                    if ( M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                    {
                        DVLOG(2) << "Eliminating row " << thedof << " using value : " << __value << "\n";
                        dofs.push_back( thedof );
                        values.push_back(  __value );
                    }

                    else if (  M_on_strategy.test( ContextOn::PENALISATION ) &&
                               !M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                    {
                        __form.set( thedof, thedof, 1.0*1e30 );
                        M_rhs->set( thedof, __value*1e30 );
                    }
                }
            }
        }
    }

    auto x = M_rhs->clone();
    CHECK( values.size() == dofs.size() ) << "Invalid dofs/values size: " << dofs.size() << "/" << values.size();
    x->setVector( dofs.data(), dofs.size(), values.data() );
    x->close();
    __form.zeroRows( dofs, *x, *M_rhs, M_on_strategy, M_value_on_diagonal );
    x.reset();
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onElements() done.");
}

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onFaces( std::shared_ptr<Elem1> const& /*__u*/,
                                                                 std::shared_ptr<Elem2> const& /*__v*/,
                                                                 FormType& __form ) const
{
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onFaces() ");
#if 0

    if ( !boost::is_same<Elem1, typename Elem::functionspace_type>::value ||
         !boost::is_same<Elem2, typename Elem::functionspace_type>::value )
        return;

#endif
    typedef typename Elem::functionspace_type functionspace_type;
    static constexpr bool is_same_space = boost::is_same<functionspace_type,Elem1>::value;
    static constexpr bool is_comp_space = Elem1::is_vectorial && Elem1::is_product && boost::is_same<functionspace_type,typename Elem1::component_functionspace_type>::value;
    
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
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    static const size_type gmc_v = is_hdiv_conforming_v<fe_type> || is_hcurl_conforming_v<fe_type> ? context|vm::JACOBIAN|vm::KB|vm::TANGENT|vm::NORMAL : context;

    constexpr uint16_type nDim = geoshape_type::nDim;

    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    if (  M_on_strategy.test( ContextOn::PENALISATION ) )
    {
        // make sure that the form is close, ie the associated matrix is assembled
        __form.matrix().close();
        // make sure that the right hand side is closed, ie the associated vector is assembled
        M_rhs->close();
    }

    //
    // start
    //
    LOG(INFO)  << fmt::format("[IntegratorOnExpr<>::onFaces()] start find a face");

    std::vector<int> dofs;
    std::vector<value_type> values;
    element_iterator __face_it = this->beginElement();
    element_iterator __face_en = this->endElement();

    bool findAFace = false;
    for( auto& lit : M_elts )
    {
        __face_it = lit.begin();
        __face_en = lit.end();
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
            __face_it = lit.begin();
            __face_en = lit.end();
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

        auto __gm = faceForInit.element( 0 ).gm();

        auto __geopc = __gm->preComputeOnFaces( __gm, [&__fe]( int f ) { return __fe->points(f); } );

        uint16_type __face_id = faceForInit.pos_first();
        auto ctx = __gm->template context<gmc_v>( faceForInit.element( 0 ), __geopc, __face_id, M_expr.dynamicContext() );

        auto expr_evaluator = M_expr.evaluator( vf::mapgmc(ctx) );

        DVLOG(2)  << "face_type::numVertices = " << face_type::numVertices << ", fe_type::nDofPerVertex = " << fe_type::nDofPerVertex << "\n"
                  << "face_type::numEdges = " << face_type::numEdges << ", fe_type::nDofPerEdge = " << fe_type::nDofPerEdge << "\n"
                  << "face_type::numFaces = " << face_type::numFaces << ", fe_type::nDofPerFace = " << fe_type::nDofPerFace << "\n";
#if 0
        index_type nbFaceDof = invalid_v<index_type>;
        if ( !fe_type::is_modal )
            nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                          face_type::numEdges * fe_type::nDofPerEdge +
                          face_type::numFaces * fe_type::nDofPerFace );
        else
            nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;
        DVLOG(2)  << "nbFaceDof = " << nbFaceDof << "\n";
        //const index_type nbFaceDof = __fe->boundaryFE()->points().size2();
#endif
        int compDofShift = (is_comp_space)? ((int)M_u.component()) : 0;
        auto const& trialDofIdToContainerId = __form.dofIdToContainerIdTrial();

        bool hasMeshSupportPartial = __dof->hasMeshSupport() && __dof->meshSupport()->isPartialSupport();


        auto IhLoc = __fe->faceLocalInterpolant();
        for( auto& lit : M_elts )
        {
        __face_it = lit.begin();
        __face_en = lit.end();
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
            if ( hasMeshSupportPartial )
            {
                if ( __dof->meshSupport()->isGhostFace( theface ) )
                    continue;
            }
            else if ( theface.isGhostFace() )
            {
                LOG(WARNING) << "face id : " << theface.id() << " is a ghost face";
                continue;
            }

            __face_id = theface.pos_first();
            uint16_type faceConnectionId = 0;
            if ( hasMeshSupportPartial )
            {
                auto const& elt0 = theface.element( 0 );
                if ( !__dof->meshSupport()->hasElement( elt0.id() ) || elt0.isGhostCell() )
                {
                    __face_id = theface.pos_second();
                    faceConnectionId = 1;
                }
                DCHECK( __dof->meshSupport()->hasElement( theface.element( faceConnectionId ).id() ) ) << "element not present in partial support";
                DCHECK( !theface.element( faceConnectionId ).isGhostCell() ) << "face connection to a ghost element is forbidden";
            }
            else if ( theface.element( 0 ).isGhostCell() )
            {
                DCHECK( theface.isConnectedTo1() ) << "invalid face, no other connection";
                __face_id = theface.pos_second();
                faceConnectionId = 1;
            }
            ctx->template update<gmc_v>( theface.element( faceConnectionId ), __face_id );

            DVLOG(2) << "FACE_ID = " << theface.id()
                     << " element id= " << ((faceConnectionId==0)? theface.ad_first() : theface.ad_second())
                     << " pos in elt= " << ((faceConnectionId==0)? theface.pos_first() : theface.pos_second());
            DVLOG(2) << "FACE_ID = " << theface.id() << " face pts=" << theface.G() << "\n";
            DVLOG(2) << "FACE_ID = " << theface.id() << "  ref pts=" << ctx->xRefs() << "\n";
            DVLOG(2) << "FACE_ID = " << theface.id() << " real pts=" << ctx->xReal() << "\n";

            expr_evaluator.update( mapgmc( ctx ) );

#if 0
            std::pair<index_type,index_type> range_dof( std::make_pair( M_u.start(),
                                                                      M_u.functionSpace()->nDof() ) );
            DVLOG(2)  << "[integratoron] dof start = " << range_dof.first << "\n";
            DVLOG(2)  << "[integratoron] dof range = " << range_dof.second << "\n";
#endif
            //use interpolant
            __fe->faceInterpolate( expr_evaluator, IhLoc );

            for( auto const& ldof : M_u.functionSpace()->dof()->faceLocalDof( theface.id() ) )
                {
                    index_type thedof = (is_comp_space)? compDofShift+Elem1::nComponents*ldof.index() : ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    DCHECK( ldof.localDofInFace() < IhLoc.size() ) 
                        << "Invalid local dof index in face for face Interpolant "
                        << ldof.localDofInFace() << ">=" << IhLoc.size();
                    double __value = ldof.sign()*IhLoc( ldof.localDofInFace() );
                    DVLOG(3) << " on " << theface.id() << " thedof "<< thedof << " = " << __value
                             << " start=" << M_u.start() << " ldof=" << ldof.index() << "\n";
                    if ( std::find( dofs.begin(),
                                    dofs.end(),
                                    thedof ) != dofs.end() )
                        continue;

                    if ( M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                        {
                            DVLOG(2) << fmt::format("Eliminating row {} using value {}, ldof.index={}", thedof, __value, ldof.index());

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
                            // __form.set( thedof, thedof, 1.0*1e30 );
                            // M_rhs->set( thedof, __value*1e30 );
                        }
                }
        }// __face_it != __face_en
        } // for( auto& lit : M_elts )
    }// findAFace

    auto x = M_rhs->clone();
    CHECK( values.size() == dofs.size() ) << "Invalid dofs/values size: " << dofs.size() << "/" << values.size();
    x->setVector( dofs.data(), dofs.size(), values.data() );
    x->close();

    LOG(INFO) << fmt::format("IntegratorOnExpr<>::onFaces() zeroRows set {} rhs values and set {} dofs", values.size(), dofs.size());
    __form.zeroRows( dofs, *x, *M_rhs, M_on_strategy, M_value_on_diagonal );
    x.reset();
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onFaces() done.");
}


template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onEdges( std::shared_ptr<Elem1> const& /*__u*/,
                                                                  std::shared_ptr<Elem2> const& /*__v*/,
                                                                  FormType& __form ) const
{


    typedef typename Elem::functionspace_type functionspace_type;
    static constexpr bool is_same_space = boost::is_same<functionspace_type,Elem1>::value;
    static constexpr bool is_comp_space = Elem1::is_vectorial && Elem1::is_product && boost::is_same<functionspace_type,typename Elem1::component_functionspace_type>::value;
    VLOG(2) << "call on::assemble(edges): " << is_comp_space<< "\n";

    if (  M_on_strategy.test( ContextOn::PENALISATION ) )
    {
        // make sure that the form is close, ie the associated matrix is assembled
        __form.matrix().close();
        // make sure that the right hand side is closed, ie the associated vector is assembled
        M_rhs->close();
    }

    //
    // start
    //
    DVLOG(2)  << "assembling Dirichlet conditions\n";

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto edge_it = this->beginElement();
    auto edge_en = this->endElement();

    bool findAEdge = false;
    for( auto& lit : M_elts )
    {
        edge_it = lit.begin();
        edge_en = lit.end();
        if ( edge_it != edge_en )
        {
            findAEdge=true;
            break;
        }
    }
    if ( findAEdge )
    {
        auto const& edgeForInit = boost::unwrap_ref( *edge_it );

        auto const* __dof = M_u.functionSpace()->dof().get();
        auto const* __fe = M_u.functionSpace()->fe().get();
        auto const* mesh = M_u.functionSpace()->mesh().get();
        auto gm = mesh->gm();
        index_type eid = edgeForInit.elements().begin()->first;
        uint16_type edgeid_in_element = edgeForInit.elements().begin()->second;
        auto const& elt = mesh->element( eid );
#if 0
        //auto geopc = gm->preComputeOnEdges([&__fe]( int f ){ return __fe->edgePoints(f); } );
        auto geopc = gm->preCompute( __fe->edgePoints(edgeid_in_element) );
        // TODO: create context for edge associated to element
        auto ctx =  gm->template context<context>( elt, geopc);
#else
        typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
        typedef typename geoelement_type::template PermutationSubEntity<2> permutation_type;
        typedef typename geoelement_type::gm_type::precompute_ptrtype geopc_ptrtype;
        std::vector<std::map<permutation_type, geopc_ptrtype> > geopc( geoelement_type::numEdges );
        for ( uint16_type __f = 0; __f < geoelement_type::numEdges; ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY ); __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                geopc[__f][__p] = gm->preCompute( __fe->edgePoints(__f) );
            }
        }
        auto ctx = gm->template context<context,2>( elt, geopc, edgeid_in_element );
#endif

        auto expr_evaluator = M_expr.evaluator( mapgmc(ctx) );
        auto IhLoc = __fe->edgeLocalInterpolant();

        int compDofShift = (is_comp_space)? ((int)M_u.component()) : 0;
        auto const& trialDofIdToContainerId = __form.dofIdToContainerIdTrial();

        for( auto& lit : M_elts )
        {
            edge_it = lit.begin();
            edge_en = lit.end();
            DVLOG(2) << "edge nb: " << std::distance(edge_it,edge_en);
            for ( ;
                  edge_it != edge_en;//this->endElement();
                  ++edge_it )
            {
                auto const& theedge = boost::unwrap_ref( *edge_it );
                // do not process the edge if it is a ghost edge: belonging to two
                // processes and being in a process id greater than the one
                // corresponding edge
                if ( theedge.isGhostCell() )
                {
                    LOG(WARNING) << "edge id : " << theedge.id() << " is a ghost edge";
                    continue;
                }

                eid = invalid_v<index_type>;
                for ( auto const& eltConnectedToEdge : theedge.elements() )
                {
                    index_type eltIdConnected = eltConnectedToEdge.first;
                    if ( __dof->isElementDone( eltIdConnected ) )
                    {
                        eid = eltIdConnected;
                        edgeid_in_element = eltConnectedToEdge.second;
                        break;
                    }
                }
                if ( eid == invalid_v<index_type> )
                    continue;

                auto const& elt = mesh->element( eid );
                //geopc = gm->preCompute( __fe->edgePoints(edgeid_in_element) );
                ////geopc = gm->preComputeAtEdges( __fe->edgePoints(ptid_in_element) );
                //ctx->update( elt, edgeid_in_element, geopc );
                ctx->template update<context>( elt, edgeid_in_element );
                expr_evaluator.update( mapgmc( ctx ) );
                __fe->edgeInterpolate( expr_evaluator, IhLoc );

                for( auto const& ldof : M_u.functionSpace()->dof()->edgeLocalDof( eid, edgeid_in_element ) )
                {
                    index_type thedof = (is_comp_space)? compDofShift+Elem1::nComponents*ldof.index() : ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];
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
        }// edge_it != edge_en

        } // for( auto& lit : M_elts )

    }// findAEdge

    auto x = M_rhs->clone();
    CHECK( values.size() == dofs.size() ) << "Invalid dofs/values size: " << dofs.size() << "/" << values.size();
    x->setVector( dofs.data(), dofs.size(), values.data() );
    x->close();
    __form.zeroRows( dofs, *x, *M_rhs, M_on_strategy, M_value_on_diagonal );
    x.reset();

}

template<typename ElementRange, typename Elem, typename RhsElem, typename OnExpr>
template<typename Elem1, typename Elem2, typename FormType>
void
IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onPoints( std::shared_ptr<Elem1> const& /*__u*/,
                                                                  std::shared_ptr<Elem2> const& /*__v*/,
                                                                  FormType& __form ) const
{
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onPoints()");
#if 0

    if ( !boost::is_same<Elem1, typename Elem::functionspace_type>::value ||
         !boost::is_same<Elem2, typename Elem::functionspace_type>::value )
        return;

#endif
    typedef typename Elem::functionspace_type functionspace_type;
    static constexpr bool is_product = functionspace_type::is_product;
    static constexpr bool is_same_space = boost::is_same<functionspace_type,Elem1>::value;
    static constexpr bool is_comp_space = Elem1::is_vectorial && Elem1::is_product && boost::is_same<functionspace_type,typename Elem1::component_functionspace_type>::value;
    VLOG(2) << "call on::assemble(MESH_POINTS): " << is_comp_space<< "\n";

    VLOG(2) << "On::assemble on Mesh Points";
    // TODO : check that we do not use hdiv hcurl or other type of elements
    const size_type context = OnExpr::context|vm::POINT;
    VLOG(2)  << "assembling Dirichlet conditions\n";
    auto mesh = M_u.functionSpace()->mesh().get();
    auto const* dof = M_u.functionSpace()->dof().get();
    auto const* __fe = M_u.functionSpace()->fe().get();
    auto gm = mesh->gm();

    if (  M_on_strategy.test( ContextOn::PENALISATION ) )
    {
        // make sure that the form is close, ie the associated matrix is assembled
        __form.matrix().close();
        // make sure that the right hand side is closed, ie the associated vector is assembled
        M_rhs->close();
    }

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto pt_it = this->beginElement();
    auto pt_en = this->endElement();

    bool findAPt = false;
    for( auto& lit : M_elts )
    {
        pt_it = lit.begin();
        pt_en = lit.end();
        if ( pt_it != pt_en )
        {
            findAPt=true;
            break;
        }
    }
    if ( findAPt )
    {
        // get the first pt properly connected
        bool findAPtToInit=false;
        for( auto& lit : M_elts )
        {
            pt_it = lit.begin();
            pt_en = lit.end();
            for( ; pt_it != pt_en; ++pt_it )
            {
                if ( boost::unwrap_ref( *pt_it).elements().size() )
                {
                    findAPtToInit=true;
                    break;
                }
            }
            if ( findAPtToInit ) break;
        }
        CHECK( findAPtToInit ) << "a point to initialize the Dirichlet constraint\n";

        auto const& thept = boost::unwrap_ref( *pt_it );

        index_type eid = thept.elements().begin()->first;
        uint16_type ptid_in_element = thept.elements().begin()->second;

        auto const& elt = mesh->element( eid );
        //auto geopc = gm->preCompute( __fe->vertexPoints(ptid_in_element) );
        //auto ctx = gm->template context<context>( elt, geopc );
        typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
        typedef typename geoelement_type::template PermutationSubEntity<geoelement_type::nDim> permutation_type;
        typedef typename geoelement_type::gm_type::precompute_ptrtype geopc_ptrtype;
        std::vector<std::map<permutation_type, geopc_ptrtype> > geopc( geoelement_type::numVertices );
        for ( uint16_type __f = 0; __f < geoelement_type::numVertices; ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY ); __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                geopc[__f][__p] = gm->preCompute( __fe->vertexPoints(__f) );
            }
        }
        auto ctx = gm->template context<context,geoelement_type::nDim>( elt, geopc, ptid_in_element );

        auto expr_evaluator = M_expr.evaluator( mapgmc(ctx) );
        auto IhLoc = __fe->vertexLocalInterpolant();

        int compDofShift = (is_comp_space)? ((int)M_u.component()) : 0;
        auto const& trialDofIdToContainerId = __form.dofIdToContainerIdTrial();

        for( auto& lit : M_elts )
        {
            pt_it = lit.begin();
            pt_en = lit.end();
            DVLOG(2) << "point " << boost::unwrap_ref( *pt_it ).id() << " nb: " << std::distance(pt_it,pt_en);

            if ( pt_it == pt_en )
                continue;

            for ( ;
                  pt_it != pt_en;
                  ++pt_it )
            {
                auto const& thept = boost::unwrap_ref( *pt_it );

                eid = invalid_v<index_type>;
                for ( auto const& eltConnectedToPoint : thept.elements() )
                {
                    index_type eltIdConnected = eltConnectedToPoint.first;
                    if ( dof->isElementDone( eltIdConnected ) )
                    {
                        eid = eltIdConnected;
                        ptid_in_element = eltConnectedToPoint.second;
                        break;
                    }
                }
                if ( eid == invalid_v<index_type> )
                    continue;

                auto const& elt = mesh->element( eid );
                //geopc = gm->preCompute( __fe->vertexPoints(ptid_in_element) );
                //ctx->update( elt, ptid_in_element, geopc, mpl::int_<0>() );
                ctx->template update<context>( elt, ptid_in_element );
                expr_evaluator.update( mapgmc( ctx ) );
                __fe->vertexInterpolate( expr_evaluator, IhLoc );

                for( int c = 0; c < (is_product?nComponents:1); ++c )
                {
                    index_type index = dof->localToGlobal( eid, ptid_in_element, c ).index();
                    index_type thedof = (is_comp_space)? compDofShift+Elem1::nComponents*index : index;
                    thedof = trialDofIdToContainerId[ thedof ];

                    double __value = IhLoc( c );

                    if ( std::find( dofs.begin(),
                                    dofs.end(),
                                    thedof ) != dofs.end() )
                        continue;

                    if ( M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                    {
                        DVLOG(3) << "Eliminating row " << thedof << " using value : " << __value << "\n";
                        dofs.push_back( thedof );
                        values.push_back(  __value );
                    }
                    else if (  M_on_strategy.test( ContextOn::PENALISATION ) &&
                               !M_on_strategy.test( ContextOn::ELIMINATION|ContextOn::SYMMETRIC ) )
                    {
                        __form.set( thedof, thedof, 1.0*1e30 );
                        M_rhs->set( thedof, __value*1e30 );
                    }
                }
            }// pt_it != pt_en
        } // for( auto& lit : M_elts )
    }// findAFace

    auto x = M_rhs->clone();
    CHECK( values.size() == dofs.size() ) << "Invalid dofs/values size: " << dofs.size() << "/" << values.size();
    x->setVector( dofs.data(), dofs.size(), values.data() );
    x->close();
    __form.zeroRows( dofs, *x, *M_rhs, M_on_strategy, M_value_on_diagonal );
    x.reset();
    LOG(INFO) << fmt::format("IntegratorOnExpr<ElementRange, Elem, RhsElem,  OnExpr>::onPoints() done");
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


template<typename ArgRangeType,typename ArgRhsType,typename ArgElementType,typename ArgExprType>
struct integratoron_type
{
    using _range_base_type = std::decay_t<ArgRangeType>;
    using _rhs_type = std::decay_t<ArgRhsType>;
    using _element_type = std::decay_t<ArgElementType>;
    using _expr_type = std::decay_t<ArgExprType>;

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
auto
getRhsVector( V const&  v )
{
    LOG(INFO) << fmt::format("call getRhsVector(V const& v) : is_vector_ptr: {} is_vectorblock: {}, is_vectorcondensed: {}",Feel::detail::is_vector_ptr_v<V>, Feel::is_vectorblock_v<V>, Feel::is_vectorcondensed_v<V>);

    if constexpr ( Feel::detail::is_vector_ptr_v<V> )
        return v;
    if constexpr ( Feel::is_vectorblock_v<V> )
        return v.getVector();
    else if constexpr ( Feel::is_vectorcondensed_v<V> )
        return v.getVector();
    else
        return v.vectorPtr();
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
template <typename ... Ts>
auto on( Ts && ... v )
{
    LOG(INFO) << fmt::format("call on(Ts&&... v) : {}",sizeof...(Ts));
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    //auto && element = args.get(_element);
    // NOTE : the line below allows to defined element as an lvalue or rvalue (as given when called the function) : TODO : add api in NApp
    auto && element = std::move(args.template getArgument<na::element>()).value();
    auto && rhs = args.get(_rhs);
    auto && expr = args.get(_expr);
    po::variables_map const& vm = args.get_else( _vm, Environment::vm() );
    std::string const& prefix = args.get_else(_prefix,"");
    std::string const& type = args.get_else_invocable(_type,[&prefix,&vm](){ return soption(_prefix=prefix,_name="on.type",_vm=vm); } );
    bool verbose = args.get_else_invocable(_verbose,[&prefix,&vm](){ return boption(_prefix=prefix,_name="on.verbose",_vm=vm); } );
    double value_on_diagonal = args.get_else_invocable(_value_on_diagonal,[&prefix,&vm](){ return doption(_prefix=prefix,_name="on.value_on_diagonal",_vm=vm); } );

    using integratoron_helper_type = vf::detail::integratoron_type<decltype( range ), decltype( Feel::vf::detail::getRhsVector( rhs ) ), decltype( element ), decltype( expr )>;

    typename integratoron_helper_type::type ion( range,
                                                 std::forward<decltype(element)>( element ),
                                                 Feel::vf::detail::getRhsVector(rhs),
                                                 expr,
                                                 size_type(ContextOnMap[type]),
                                                 value_on_diagonal );
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
    auto r = typename integratoron_helper_type::expr_type( std::move(ion) );
    LOG(INFO) << fmt::format("call on(Ts&&... v) : {} done",sizeof...(Ts));
    //typename vf::detail::integratoron_type<Args>::type ion( range, element, rhs, expr, type );
    return r;
}


} // vf
} // feel
#endif
