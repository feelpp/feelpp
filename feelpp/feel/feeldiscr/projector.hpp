/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
       Date: 2011-04-25

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file projector.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2011-04-25
 */
#ifndef FEELPP_DISCR_PROJECTOR_H
#define FEELPP_DISCR_PROJECTOR_H

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelmesh/intersect.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
//#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/twovalued.hpp>
#include <feel/feelvf/trace.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelmesh/ranges.hpp>

namespace Feel
{
    template<class DomainSpace, class DualImageSpace> class Projector;

#if 0
namespace detail
{
template<typename Args>
struct projector_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::element_type domain_type;
    typedef typename vf::detail::clean_type<Args,tag::imageSpace>::element_type image_type;
    typedef std::shared_ptr<Projector<domain_type,image_type> > return_type;
};

template<typename Args>
struct lift_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::element_type domain_type;
    typedef std::shared_ptr<Projector<domain_type,domain_type> > lift_return_type;
};

} // detail
#endif
/**
 * \class Projector
 * \brief Projection made easy
 *
 * @author Vincent Doyeux
 * @see OperatorLinear
 */
template<class DomainSpace, class DualImageSpace>
class Projector : public OperatorLinear< functionspace_type<DomainSpace>,functionspace_type<DualImageSpace> >
{
    typedef OperatorLinear< functionspace_type<DomainSpace>,functionspace_type<DualImageSpace> > super;
    typedef Projector<DomainSpace, DualImageSpace> self_type;

public :

    /** @name Typedefs
     */
    //@{
    typedef super ol_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::dual_image_space_type  dual_image_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename super::dual_image_space_ptrtype  dual_image_space_ptrtype;
    typedef typename domain_space_type::element_type domain_element_type;
    typedef typename dual_image_space_type::element_type dual_image_element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;

    typedef Range<typename dual_image_space_type::mesh_type, MESH_ELEMENTS> range_elements_type;
    typedef Range<typename dual_image_space_type::mesh_type, MESH_FACES> range_faces_type;
#if 0
    template<typename Args,typename IntEltsDefault>
    struct integrate_type
    {
        typedef vf::detail::clean_type<Args,tag::expr> _expr_type;
        typedef vf::detail::clean2_type<Args,tag::range,IntEltsDefault> _range_type;
        typedef typename boost::tuples::template element<1, _range_type>::type _element_iterator;
        static inline const uint16_type geoOrder = boost::unwrap_reference<typename _element_iterator::value_type>::type::nOrder;

        using expr_order_t = ExpressionOrder<_range_type,_expr_type>;
        using im_default_type = im_t<typename expr_order_t::the_element_type, typename _expr_type::value_type>;
        typedef vf::detail::clean2_type<Args,tag::quad, im_default_type> __quad_type;
        typedef vf::detail::clean2_type<Args,tag::quad1, im_default_type > __quad1_type;
        using _im_type = vf::detail::integrate_im_type<_range_type,_expr_type,__quad_type,__quad1_type>;
        using _quad_type = typename _im_type::_quad_type;
        using _quad1_type = typename _im_type::_quad1_type;

        
        //typedef vf::detail::clean2_type<Args,tag::quad, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value > > _quad_type;
        //typedef vf::detail::clean2_type<Args,tag::quad1, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value_1 > > _quad1_type;
    };
#endif
    //@}
    /** @name Constructors, destructor
     */
    //@{

    Projector( domain_space_ptrtype     domainSpace,
               dual_image_space_ptrtype dualImageSpace,
               backend_ptrtype abackend = Backend<double>::build( soption( _name="backend" ) ),
               ProjectorType proj_type=L2,
               double epsilon = 0.01,
               double gamma = 20,
               DirichletType dirichlet_type = WEAK
             )
        :
        ol_type( domainSpace, dualImageSpace, abackend, false ),
        //M_backend( abackend ),
        M_epsilon( epsilon ),
        M_gamma( gamma ),
        M_proj_type( proj_type ),
        M_dir( dirichlet_type )

    {
        if ( M_proj_type == CIP && domainSpace->dof()->worldComm().localSize() > 1 )
            CHECK( domainSpace->dof()->buildDofTableMPIExtended() && dualImageSpace->dof()->buildDofTableMPIExtended() ) << "functionspaces must have extended dof table";

        size_type thepattern = ( M_proj_type == CIP )? size_type(Pattern::EXTENDED) : size_type(Pattern::COUPLED);
        M_matrixFull = this->backend()->newMatrix( _trial=this->dualImageSpace(), _test=this->domainSpace(),_pattern=thepattern );

        this->matPtr() = M_matrixFull;

        if ( M_proj_type == LIFT )
            M_matrixCst = this->backend()->newMatrix( _trial=this->dualImageSpace(), _test=this->domainSpace(),_pattern=thepattern );
        else
            M_matrixCst = M_matrixFull; // same pointer

        M_ie = this->backend()->newVector( this->domainSpace() );

        // Init ranges
        M_rangeElements = self_type::rangeElements( this->domainSpace(), this->dualImageSpace() );
        M_rangeBoundaryFaces = self_type::rangeBoundaryFaces( this->domainSpace(), this->dualImageSpace() );

        initMatrix<domain_element_type>();
    }

    ~Projector() override {}

    //~Projector() {}
    //@}

    /** @name  Methods
     */
    //@{

    static range_elements_type rangeElements( domain_space_ptrtype const& domainSpace, dual_image_space_ptrtype const& imageSpace )
        {
            return rangeElements( domainSpace,imageSpace, mpl::bool_<boost::is_same<dual_image_element_type,domain_element_type>::value>() );
        }

    static range_elements_type rangeElements( domain_space_ptrtype const& domainSpace, dual_image_space_ptrtype const& imageSpace, mpl::true_ )
    {
        if ( !imageSpace->mesh()->isSameMesh( domainSpace->mesh() ) )
            return elements( imageSpace->mesh() );

        bool hasMeshSupportPartialDomain = domainSpace->dof()->hasMeshSupport()
            && domainSpace->dof()->meshSupport()->isPartialSupport();
        bool hasMeshSupportPartialImage = imageSpace->dof()->hasMeshSupport()
            && imageSpace->dof()->meshSupport()->isPartialSupport();

        range_elements_type rangeElts;
        if ( hasMeshSupportPartialDomain && hasMeshSupportPartialImage )
            rangeElts = intersect( domainSpace->dof()->meshSupport()->rangeElements(), imageSpace->dof()->meshSupport()->rangeElements() );
        else if ( hasMeshSupportPartialDomain )
            rangeElts = domainSpace->dof()->meshSupport()->rangeElements();
        else if ( hasMeshSupportPartialImage )
            rangeElts = imageSpace->dof()->meshSupport()->rangeElements();
        else
            rangeElts = elements( imageSpace->mesh() );

        return rangeElts;
    }
    static range_elements_type rangeElements( domain_space_ptrtype const& domainSpace, dual_image_space_ptrtype const& imageSpace, mpl::false_ )
        {
            return elements( imageSpace->mesh() );
        }

    static range_faces_type rangeBoundaryFaces( domain_space_ptrtype const& domainSpace, dual_image_space_ptrtype const& imageSpace )
    {
        range_faces_type rangeFaces;

        bool hasMeshSupportPartialDomain = domainSpace->dof()->hasMeshSupport() 
            && domainSpace->dof()->meshSupport()->isPartialSupport();
        bool hasMeshSupportPartialImage = imageSpace->dof()->hasMeshSupport() 
            && imageSpace->dof()->meshSupport()->isPartialSupport();
        if ( hasMeshSupportPartialDomain && hasMeshSupportPartialImage )
            rangeFaces = intersect( domainSpace->dof()->meshSupport()->rangeBoundaryFaces(), imageSpace->dof()->meshSupport()->rangeBoundaryFaces() );
        else if ( hasMeshSupportPartialDomain )
            rangeFaces = domainSpace->dof()->meshSupport()->rangeBoundaryFaces();
        else if ( hasMeshSupportPartialImage )
            rangeFaces = imageSpace->dof()->meshSupport()->rangeBoundaryFaces();
        else
            rangeFaces = boundaryfaces( imageSpace->mesh() );

        return rangeFaces;
    }

    template<typename Range, typename Expr, typename Elem>
    void applyOn( Range range, Expr const& expr, Elem const& de )
    {
        applyOn( range, expr, de, mpl::int_<Range::entities()>() );
    }

    template<typename RangeT, typename ExprT>
    dual_image_element_type
    projectL2( RangeT range, ExprT const& expr )
    {
        auto uDomain = this->domainSpace()->element();
        auto sol = this->dualImageSpace()->element();

        M_ie->zero();

        form1( _test=this->domainSpace(), _vector=M_ie ) +=
            integrate( _range=range, _expr=inner(expr,id( uDomain ) ) );

        this->backend()->solve( _matrix=M_matrixFull, _solution=sol, _rhs=M_ie );

        return sol;
    }

    template <typename ... Ts>
    dual_image_element_type project( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && expr = args.get(_expr);
        auto && range = args.get_else_invocable(_range, [this]() { return self_type::rangeElements(this->domainSpace(), this->dualImageSpace()); } );
        auto && quad = args.get_else(_quad,quad_order_from_expression );
        auto && quad1 = args.get_else(_quad1,quad_order_from_expression );
        GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT);
        auto && grad_expr = args.get_else(_grad_expr, vf::zero<domain_space_type::nComponents,domain_space_type::nDim>() );
        auto && div_expr = args.get_else(_div_expr,cst(0.) );
        auto && curl_expr = args.get_else(_curl_expr,vf::zero< mpl::if_<mpl::equal_to<mpl::int_<domain_space_type::nComponents>, mpl::int_<1> >,
                                          mpl::int_<1>,
                                          typename mpl::if_<mpl::equal_to<mpl::int_<domain_space_type::nDim>, mpl::int_<3> >,
                                          mpl::int_<3>, mpl::int_<1> >::type >::type::value, 1>() );



        using namespace vf;

        using _integrate_helper_type = Feel::detail::integrate_type<decltype(expr),decltype(range),decltype(quad),decltype(quad1)>;
        auto the_ims = _integrate_helper_type::_im_type::im( quad,quad1,expr );
        auto const& the_im = the_ims.first;
        auto const& the_im1 = the_ims.second;


#if 0
        auto the_ims = integrate_type<Args,range_elements_type>::_im_type::im( quad,quad1,expr );
        auto const& the_im = the_ims.first;
        auto const& the_im1 = the_ims.second;
#endif
        //typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad)>::type >::type thequad_type;
        //typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad1)>::type >::type thequad1_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(range)>::type >::type therange_type;
        using element_iterator = typename therange_type::iterator_t;
        using element_t = typename therange_type::element_t;
        constexpr uint16_type geoOrder = element_t::nOrder;
        constexpr uint16_type nOrderImageSpace = dual_image_space_type::basis_type::nOrder;
        constexpr uint16_type quadOrderId = nOrderImageSpace*geoOrder;
        constexpr uint16_type quadOrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1)*geoOrder:0;
        constexpr uint16_type quad1OrderId = nOrderImageSpace;
        constexpr uint16_type quad1OrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1):0;

        auto sol = this->dualImageSpace()->element();

        this->setRHSAndBC<domain_element_type>( 
                sol, 
                expr, range, 
                the_im,the_im1,//quad, quad1, 
                geomap,
                grad_expr, div_expr, curl_expr
                );

        this->backend()->solve( _matrix=M_matrixFull, _solution=sol, _rhs=M_ie );

        return sol;
    }

    template<typename RhsExpr>
    dual_image_element_type
    operator()( RhsExpr const& rhs_expr )
    {
        return this->project( _expr=rhs_expr );
    }
    template<typename RhsExpr>
    void
    operator()( dual_image_element_type& de, RhsExpr const& rhs_expr )
    {
        de = this->project( _expr=rhs_expr );
    }

    dual_image_element_type
    operator()( domain_element_type const& de )
    {
        auto ie = this->dualImageSpace()->element();
        this->apply( de, ie );
        return ie ;
    }

    void
    operator()( domain_element_type const& de,dual_image_element_type &ie )
    {
        this->apply( de, ie );
    }

    template<typename Range, typename Expr>
    dual_image_element_type
    operator()( Range const& range ,Expr const& expr )
    {
        return this->project( _expr=expr, _range=range );
    }

    void
    apply( domain_element_type const& de,
           dual_image_element_type& ie ) override
    {
        this->backend()->solve( _matrix=M_matrixFull, _solution=ie, _rhs=this->backend()->toBackendVectorPtr( de ) );
    }

    template<typename RhsExpr>
    void
    apply( RhsExpr const& rhs_expr,
           dual_image_element_type& ie )
    {
        ie = this->project( _expr=rhs_expr );
    }


    template< typename Expr>
    dual_image_element_type derivate( Expr const& expr )
    {
        return this->derivateImpl<domain_element_type>( expr );
    }

    double epsilon() const { return M_epsilon; }
    double gamma() const { return M_gamma; }
    ProjectorType projectorType() const { return M_proj_type; }

    //@}


private :

    template<typename T>
        void initMatrix( typename std::enable_if<is_tensor2_field<T>::value>::type* = nullptr )
    {
        auto uDomain = this->domainSpace()->element();
        auto uImage = this->dualImageSpace()->element();


        auto a = form2 ( _trial=this->dualImageSpace(),
                         _test=this->domainSpace(),
                         _matrix=M_matrixCst );

        a = integrate( _range=this->M_rangeElements,
                       _expr=inner( idt( uImage ), id( uDomain ) ) );
    }
    template<typename T>
        void initMatrix( typename std::enable_if<mpl::not_<is_tensor2_field<T>>::value>::type* = nullptr )
    {
        using namespace vf;
        auto uDomain = this->domainSpace()->element();
        auto uImage = this->dualImageSpace()->element();

        auto a = form2 ( _trial=this->dualImageSpace(),
                         _test=this->domainSpace(),
                         _matrix=M_matrixCst );

        if ( M_proj_type != LIFT )
        {
            a = integrate( _range=this->M_rangeElements,
                           _expr=inner( idt( uImage ), id( uDomain ) ) );
        }
        switch ( M_proj_type )
        {
        case L2:
        {
        }
        break;

        case H1:
        {
            a += integrate( _range=this->M_rangeElements,
                            _expr=inner( gradt( uImage ), grad( uDomain ) ) );
        }
        break;

        case DIFF:
        {
            a += integrate( _range=this->M_rangeElements,
                            _expr=M_epsilon*inner( gradt( uImage ), grad( uDomain ) ) );

            //weak boundary conditions
            a += integrate( _range=this->M_rangeBoundaryFaces,
                            _expr= M_epsilon*( -trans( id( uDomain ) )*gradt( uImage )*vf::N() ) );
            a += integrate( _range=this->M_rangeBoundaryFaces,
                            _expr= M_epsilon*( -trans( idt( uImage ) )* grad( uDomain )*vf::N() ) );
            a += integrate( _range=this->M_rangeBoundaryFaces,
                            _expr= M_gamma * trans( idt( uImage ) ) /*trial*/
                            *id( uDomain ) / vf::hFace()   /*test*/
                            );
        }
        break;

        case HDIV:
        {
            a += integrate( _range=this->M_rangeElements,
                            _expr=divt( uImage )*div( uDomain ) );
        }
        break;

        case HCURL:
        {
            a += integrate( _range=this->M_rangeElements,
                           // only for 2D, need to specialize this for 3D
                            _expr=curlzt( uImage )*curlz( uDomain ) );
        }
        break;

        case LIFT:
        {
            a = integrate( _range=this->M_rangeElements,
                           _expr= inner( gradt( uImage ), grad( uDomain ) ) );
        }
        break;

        case CIP:
        {
            a += integrate( _range=internalfaces( this->dualImageSpace()->mesh() ),
                            _expr=M_gamma*hFace()*hFace()*
                            /**/  inner( jumpt( gradt(uImage) ),jump( grad(uDomain) ) ) );
        }
        break;

        case NODAL:
            break;

        }

        M_matrixCst->close();
    }

    template<typename T, 
        typename ExprT, typename RangeT, 
        typename QuadT, typename Quad1T, 
        typename GradExprT, typename DivExprT, typename CurlExprT>
    void setRHSAndBC( 
            dual_image_element_type const& sol, 
            ExprT expr, RangeT range,
            QuadT quad, Quad1T quad1,
            GeomapStrategyType geomap,
            GradExprT grad_expr, DivExprT div_expr, CurlExprT curl_expr,
            typename std::enable_if<is_tensor2_field<T>::value>::type* = nullptr )
    {
        if( M_proj_type != L2 )
            throw std::logic_error( "Only L2 projection supported for rank-2 tensors" );

        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad)>::type >::type thequad_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad1)>::type >::type thequad1_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(range)>::type >::type therange_type;
        typedef typename boost::tuples::template element<1, therange_type>::type element_iterator;
        constexpr uint16_type geoOrder = boost::unwrap_reference<typename element_iterator::value_type>::type::nOrder;
        constexpr uint16_type nOrderImageSpace = dual_image_space_type::basis_type::nOrder;
        constexpr quad_order_type quadOrder = quad.order();//thequad_type::CompileTimeOrder;
        constexpr quad_order_type quadOrderId = nOrderImageSpace*geoOrder;
        constexpr quad_order_type quad1Order = quad1.order();//thequad1_type::CompileTimeOrder;
        constexpr quad_order_type quad1OrderId = nOrderImageSpace;

        //auto uImage = this->dualImageSpace()->element();
        auto uDomain = this->domainSpace()->element();
        M_ie->zero();

        form1( _test=this->dualImageSpace(), _vector=M_ie ) += integrate( 
            _range=range, _expr=inner(expr,id( uDomain ) ),
                _quad=quadOrder+quadOrderId,
                _quad1=quad1Order+quad1OrderId,
                _geomap=geomap 
                );
    }

    template<typename T, 
        typename ExprT, typename RangeT, 
        typename QuadT, typename Quad1T, 
        typename GradExprT, typename DivExprT, typename CurlExprT>
    void setRHSAndBC( 
            dual_image_element_type const& sol, 
            ExprT expr, RangeT range,
            QuadT quad, Quad1T quad1,
            GeomapStrategyType geomap,
            GradExprT grad_expr, DivExprT div_expr, CurlExprT curl_expr,
            typename std::enable_if<mpl::not_<is_tensor2_field<T>>::value>::type* = nullptr )
    {
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad)>::type >::type thequad_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad1)>::type >::type thequad1_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(range)>::type >::type therange_type;
        using element_iterator = typename therange_type::iterator_t;
        using element_t = typename therange_type::element_t;
        constexpr uint16_type geoOrder = element_t::nOrder;
        constexpr uint16_type nOrderImageSpace = dual_image_space_type::basis_type::nOrder;
        const quad_order_type quadOrder = quad.order();//thequad_type::CompileTimeOrder;
        constexpr quad_order_type quadOrderId = nOrderImageSpace*geoOrder;
        constexpr quad_order_type quadOrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1)*geoOrder:0;
        const quad_order_type quad1Order = quad1.order();//thequad1_type::CompileTimeOrder;
        constexpr quad_order_type quad1OrderId = nOrderImageSpace;
        constexpr quad_order_type quad1OrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1):0;

        auto uImage = this->dualImageSpace()->element();
        auto uDomain = this->domainSpace()->element();
        M_ie->zero();

        if ( M_proj_type != LIFT )
        {
            //typedef typename integrate_type<Args,decltype( elements( this->dualImageSpace()->mesh() ) )>::_quad_type myquad;
            form1( _test=this->domainSpace(), _vector=M_ie ) +=
                integrate( _range=range, _expr=inner(expr,id( uDomain ) ),
                           //integrate( _range=range, _expr=trans(expr)*id( uImage ),
                           _quad=quadOrder+quadOrderId,
                           _quad1=quad1Order+quad1OrderId,
                           _geomap=geomap );

            switch( M_proj_type )
            {
            case H1:
                form1( _test=this->domainSpace(), _vector=M_ie ) +=
                    integrate( _range=range, _expr=trace(grad_expr*trans(grad( uDomain )) ),
                               _quad=quadOrder+quadOrderGrad,
                               _quad1=quad1Order+quad1OrderGrad,
                               _geomap=geomap );
                break;
            case HDIV:
                form1( _test=this->dualImageSpace(), _vector=M_ie ) +=
                    integrate( _range=range, _expr=div_expr*div( uDomain ),
                               _quad=quadOrder+quadOrderGrad,//quad,
                               _quad1=quad1Order+quad1OrderGrad,//quad1,
                               _geomap=geomap );
                break;
            case HCURL:
                form1( _test=this->domainSpace(), _vector=M_ie ) +=
                    integrate( _range=range, _expr=trans(curl_expr)*curl( uDomain ),
                               _quad=quadOrder+quadOrderGrad,
                               _quad1=quad1Order+quad1OrderGrad,
                               _geomap=geomap );
                break;
            case L2:
            default:
                break;
            }
        }

        else if ( ( M_proj_type == LIFT ) && ( M_dir == WEAK ) )
        {
            form1( _test=this->domainSpace(), _vector=M_ie ) +=
                integrate( _range=range,
                           _expr=inner( expr, -grad( uImage )*vf::N() +
                                        M_gamma / vf::hFace() *id( uImage ) ),
                           _quad=quadOrder+quadOrderId,
                           _quad1=quad1Order+quad1OrderId,
                           _geomap=geomap );
        }

        //weak boundary conditions
        if ( M_proj_type == DIFF )
        {
            form1( _test=this->domainSpace(), _vector=M_ie ) +=
                integrate( _range=this->M_rangeBoundaryFaces,
                           _expr=inner( expr, M_epsilon*( -grad( uImage )*vf::N() ) +
                                                          M_gamma / vf::hFace() *id( uImage ) ),
                           _quad=quadOrder+quadOrderId,
                           _quad1=quad1Order+quad1OrderId,
                           _geomap=geomap );
        }


        if ( M_proj_type == LIFT )
        {
            M_matrixFull->zero();
            M_matrixFull->addMatrix( 1., M_matrixCst );
            auto bilinearForm = form2( _trial=this->dualImageSpace(), _test=this->domainSpace(), _matrix=M_matrixFull );

            if ( M_dir == WEAK )
            {
                bilinearForm +=
                    integrate( _range=range, _expr=
                               -inner( id( uDomain ), gradt( sol )*vf::N() )
                               -inner( idt( sol ), grad( uDomain )*vf::N() )
                               + M_gamma * inner( idt( sol ), id( uDomain ) ) / vf::hFace(),
                               _quad=quadOrder+quadOrderId,
                               _quad1=quad1Order+quad1OrderId,
                               _geomap=geomap );
            }
            else if ( M_dir == STRONG )
            {
                this->applyOn( range, expr, sol );
            }
        }
    }

    template<typename T, typename ExprT,
        typename std::enable_if<mpl::not_<is_tensor2_field<T>>::value, int>::type = 0>
    dual_image_element_type derivateImpl( ExprT const& expr )
    {
        auto uDomain = this->domainSpace()->element();
        auto uImage = this->dualImageSpace()->element();
        M_ie->zero();
        form1( _test=this->domainSpace(), _vector=M_ie ) +=
            integrate(_range = elements( this->dualImageSpace()->mesh() ),
                      _expr = - trace( expr * trans( grad( uDomain ) ) ) );

        form1(_test=this->domainSpace(), _vector=M_ie ) +=
            integrate(_range = boundaryfaces( this->dualImageSpace()->mesh() ),
                      _expr =  trans( id( uDomain ) ) * expr * vf::N() );

        this->backend()->solve( _matrix=M_matrixFull, _solution=uImage, _rhs=M_ie );
        return uImage;
    }

    template<typename T, typename ExprT,
        typename std::enable_if<is_tensor2_field<T>::value, int>::type = 0>
    dual_image_element_type derivateImpl( ExprT const& expr )
    {
        static_assert( mpl::not_<is_tensor2_field<T>>::value, 
                "Derivation in projector class is not yet supported for tensor2 fields." );
    }

    template<typename Range, typename Expr, typename Elem >
    void applyOn( Range range, Expr const& expr, Elem const& ie, mpl::int_<MESH_ELEMENTS> ){}

    template<typename Range, typename Expr, typename Elem>
    void applyOn( Range range, Expr const& expr, Elem const& ie, mpl::int_<MESH_FACES> )
    {
        form2 ( _trial=this->dualImageSpace(),
                _test=this->domainSpace(),
                _matrix=M_matrixFull ) +=  on( _range=range , _element=ie, _rhs=M_ie, _expr=expr );
    }

    //backend_ptrtype M_backend;
    const double M_epsilon;
    const double M_gamma;
    const ProjectorType M_proj_type;
    DirichletType M_dir;
    matrix_ptrtype M_matrixCst;
    matrix_ptrtype M_matrixFull;
    //domain_element_type M_de;
    vector_ptrtype M_ie;
    range_elements_type M_rangeElements;
    range_faces_type M_rangeBoundaryFaces;

};//Projector

template<class DomainSpace, class DualImageSpace>
using projector_type =  Projector<DomainSpace,DualImageSpace>;

template<class DomainSpace, class DualImageSpace>
using projector_ptrtype =  std::shared_ptr<Projector<DomainSpace,DualImageSpace>>;



/**
 * this function returns a \c Projector \c shared_ptr with
 *
 * \param domainSpace
 * \param imageSpace
 * \param backend
 */

template<typename TDomainSpace, typename TDualImageSpace>
std::shared_ptr< Projector<TDomainSpace, TDualImageSpace> >
projector( std::shared_ptr<TDomainSpace> const& domainspace,
           std::shared_ptr<TDualImageSpace> const& imagespace,
           typename Projector<TDomainSpace, TDualImageSpace>::backend_ptrtype const& abackend = Backend<double>::build( soption( _name="backend" ) ),
           ProjectorType proj_type=L2, double epsilon=0.01, double gamma = 20, DirichletType dirichlet_type = WEAK)
{
    typedef Projector<TDomainSpace, TDualImageSpace > Proj_type;
    return std::make_shared<Proj_type>( domainspace, imagespace, abackend, proj_type, epsilon, gamma, dirichlet_type ) ;
}


template <typename ... Ts>
auto opProjection( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    auto && imageSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_imageSpace);
    ProjectorType type = args.get_else(_type,L2);
    double penaldir = args.get_else(_penaldir,20.);
    auto && backend = args.get_else(_backend, Backend<double>::build( soption( _name="backend" ) ) );

    return projector( domainSpace,imageSpace, backend, type, 0.01, penaldir );
}

template <typename ... Ts>
auto opLift( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.template get<NA::constraint::is_convertible<std::shared_ptr<FunctionSpaceBase>>::apply>(_domainSpace);
    DirichletType type = args.get_else(_type,WEAK);
    double penaldir = args.get_else(_penaldir,20.);
    auto && backend = args.get_else_invocable(_backend, [](){ return Backend<double>::build( soption( _name="backend" ) ); } );

    return projector( domainSpace, domainSpace, backend, ProjectorType::LIFT, 0.01 , penaldir, type );
}

} //namespace Feel


#endif
