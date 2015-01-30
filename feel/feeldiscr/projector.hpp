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
#ifndef _PROJECTOR_HPP_
#define _PROJECTOR_HPP_

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/backend.hpp>
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

namespace Feel
{
    template<class DomainSpace, class DualImageSpace> class Projector;

namespace detail
{
template<typename Args>
struct projector_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::type::element_type domain_type;
    typedef typename vf::detail::clean_type<Args,tag::imageSpace>::type::element_type image_type;
    typedef boost::shared_ptr<Projector<domain_type,image_type> > return_type;
};

template<typename Args>
struct lift_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::type::element_type domain_type;
    typedef boost::shared_ptr<Projector<domain_type,domain_type> > lift_return_type;
};

} // detail
/**
 * \class Projector
 * \brief Projection made easy
 *
 * @author Vincent Doyeux
 * @see OperatorLinear
 */
template<class DomainSpace, class DualImageSpace>
class Projector : public OperatorLinear<DomainSpace, DualImageSpace>
{
    typedef Projector<DomainSpace,DualImageSpace> super;

public :

    /** @name Typedefs
     */
    //@{


    typedef OperatorLinear<DomainSpace, DualImageSpace> ol_type;

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
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;

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
        M_matrixFull = this->backend()->newMatrix( _trial=this->domainSpace(), _test=this->dualImageSpace() );

        this->matPtr() = M_matrixFull;

        if ( M_proj_type == LIFT )
            M_matrixCst = this->backend()->newMatrix( _trial=this->domainSpace(), _test=this->dualImageSpace() ) ;
        else
            M_matrixCst = M_matrixFull; // same pointer

        ie = this->backend()->newVector( this->dualImageSpace() );

        initMatrix();
    }

    //~Projector() {}
    //@}

    /** @name  Methods
     */
    //@{
    template<typename Args,typename IntEltsDefault>
    struct integrate_type
    {
        typedef typename vf::detail::clean_type<Args,tag::expr>::type _expr_type;
        typedef typename vf::detail::clean2_type<Args,tag::range,IntEltsDefault>::type _range_type;
        typedef typename boost::tuples::template element<1, _range_type>::type _element_iterator;
        static const uint16_type geoOrder = _element_iterator::value_type::nOrder;

        typedef typename vf::detail::clean2_type<Args,tag::quad, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value > >::type _quad_type;
        typedef typename vf::detail::clean2_type<Args,tag::quad1, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value_1 > >::type _quad1_type;
    };

    template<typename Range, typename Expr>
    void applyOn( Range range, Expr expr)
    {
        typedef typename boost::tuples::template element<0, Range>::type idim_type;
        applyOn( range, expr, mpl::int_<idim_type::value>() );
    }


    BOOST_PARAMETER_MEMBER_FUNCTION( ( domain_element_type ),
                                     project,
                                     tag,
                                     ( required
                                       ( expr,   * )
                                     )
                                     ( optional
                                       ( range,   *, elements( this->dualImageSpace()->mesh() )  )
                                       ( quad,   *, ( typename integrate_type<Args,decltype( elements( this->dualImageSpace()->mesh() ) )>::_quad_type() ) )
                                       ( quad1,   *, ( typename integrate_type<Args,decltype( elements( this->dualImageSpace()->mesh() ) )>::_quad1_type() ) )
                                       ( geomap, *, (integrate_type<Args,decltype( elements( this->dualImageSpace()->mesh() ) )>::geoOrder > 1 )?
                                         GeomapStrategyType::GEOMAP_OPT:GeomapStrategyType::GEOMAP_HO )
                                       (grad_expr, *, ( vf::zero<domain_space_type::nComponents,domain_space_type::nDim>() ))
                                       (div_expr, *, cst(0.) )
                                       (curl_expr, *,  ( vf::zero<  mpl::if_<mpl::equal_to<mpl::int_<domain_space_type::nComponents>, mpl::int_<1> >,
                                                         mpl::int_<1>, mpl::int_<domain_space_type::nDim> >::type::value, 1>() ) )
                                       )
                                   )
    {
        using namespace vf;

        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad)>::type >::type thequad_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(quad1)>::type >::type thequad1_type;
        typedef typename boost::remove_reference<typename boost::remove_const< decltype(range)>::type >::type therange_type;
        typedef typename boost::tuples::template element<1, therange_type>::type element_iterator;
        static const uint16_type geoOrder = element_iterator::value_type::nOrder;
        static const uint16_type nOrderImageSpace = dual_image_space_type::basis_type::nOrder;
        static const uint16_type quadOrderId = nOrderImageSpace*geoOrder;
        static const uint16_type quadOrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1)*geoOrder:0;
        static const uint16_type quad1OrderId = nOrderImageSpace;
        static const uint16_type quad1OrderGrad = (nOrderImageSpace>0)?(nOrderImageSpace-1):0;

        auto sol = this->domainSpace()->element();
        auto uImage = this->dualImageSpace()->element();

        ie->zero();

        if ( M_proj_type != LIFT )
        {
            //typedef typename integrate_type<Args,decltype( elements( this->dualImageSpace()->mesh() ) )>::_quad_type myquad;
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=range, _expr=inner(expr,id( this->dualImageSpace()->element() ) ),
                //integrate( _range=range, _expr=trans(expr)*id( uImage ),
                           _quad=_Q<thequad_type::order+quadOrderId>(),
                           _quad1=_Q<thequad1_type::order+quad1OrderId>(),
                           _geomap=geomap );

            switch( M_proj_type )
                {
                case H1:
                    form1( _test=this->dualImageSpace(), _vector=ie ) +=
                        integrate( _range=range, _expr=trace(grad_expr*trans(grad( uImage )) ),
                                   _quad=_Q<thequad_type::order+quadOrderGrad>(),
                                   _quad1=_Q<thequad1_type::order+quad1OrderGrad>(),
                                   _geomap=geomap );
                    break;
                case HDIV:
                            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                                integrate( _range=range, _expr=div_expr*div( uImage ),
                                           _quad=_Q<thequad_type::order+quadOrderGrad>(),//quad,
                                           _quad1=_Q<thequad1_type::order+quad1OrderGrad>(),//quad1,
                                           _geomap=geomap );
                            break;
                case HCURL:
                            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                                integrate( _range=range, _expr=trans(curl_expr)*curl( uImage ),
                                           _quad=_Q<thequad_type::order+quadOrderGrad>(),
                                           _quad1=_Q<thequad1_type::order+quad1OrderGrad>(),
                                           _geomap=geomap );
                            break;
                case L2:
                default:
                    break;
                }
        }
        else if ( ( M_proj_type == LIFT ) && ( M_dir == WEAK ) )
        {
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=range,
                           _expr=inner( expr, -grad( uImage )*vf::N() +
                                        M_gamma / vf::hFace() *id( uImage ) ),
                           _quad=_Q<thequad_type::order+quadOrderId>(),
                           _quad1=_Q<thequad1_type::order+quad1OrderId>(),
                           _geomap=geomap );
        }

        //weak boundary conditions
        if ( M_proj_type == DIFF )
        {
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=boundaryfaces( this->dualImageSpace()->mesh() ),
                           _expr=inner( expr, M_epsilon*( -grad( uImage )*vf::N() +
                                                          M_gamma / vf::hFace() *id( uImage ) ) ),
                           _quad=_Q<thequad_type::order+quadOrderId>(),
                           _quad1=_Q<thequad1_type::order+quad1OrderId>(),
                           _geomap=geomap );
        }


        if ( M_proj_type == LIFT )
        {
            M_matrixFull->zero();
            M_matrixFull->addMatrix( 1., M_matrixCst );
            auto bilinearForm = form2( _trial=this->domainSpace(), _test=this->dualImageSpace(), _matrix=M_matrixFull );

            if ( M_dir == WEAK )
            {
                bilinearForm +=
                    integrate( _range=range, _expr=
                               ( -trans( id( this->dualImageSpace()->element() ) )*gradt( this->domainSpace()->element() )*vf::N()
                                 -trans( idt( this->domainSpace()->element() ) )* grad( this->dualImageSpace()->element() )*vf::N()
                                 + M_gamma * trans( idt( this->domainSpace()->element() ) ) /*trial*/
                                 *id( this->dualImageSpace()->element() ) / vf::hFace()   /*test*/
                                 ),
                               _quad=_Q<thequad_type::order+quadOrderId>(),
                               _quad1=_Q<thequad1_type::order+quad1OrderId>(),
                               _geomap=geomap );
            }
            else if ( M_dir == STRONG )
            {
                this->applyOn(range, expr);
            }
        }

        this->backend()->solve( _matrix=M_matrixFull, _solution=sol, _rhs=ie );

        return sol;
    }

    template<typename RhsExpr>
    domain_element_type
    operator()( RhsExpr const& rhs_expr )
    {
        return this->project( rhs_expr );
    }

    template<typename RhsExpr>
    void
    operator()( domain_element_type& de, RhsExpr const& rhs_expr )
    {
        de = this->project( rhs_expr );
    }

    domain_element_type
    operator()( image_element_type const& ie )
    {
        domain_element_type de = this->domainSpace()->element();
        this->backend()->solve( M_matrixFull, de, ie );
        return de ;
    }

    void
    operator()( domain_element_type &de, image_element_type const& ie )
    {
        this->backend()->solve( M_matrixFull, de, ie );
    }

    template<typename Range, typename Expr>
    domain_element_type
    operator()( Range const& range ,Expr const& expr )
    {
        return this->project( expr, range );
    }

		using OperatorLinear<DomainSpace, DualImageSpace>::apply;
    void
    apply( domain_element_type& de,
           image_element_type const& ie )
    {
        this->backend()->solve( M_matrixFull, de, ie );
    }

    template<typename RhsExpr>
    void
    apply( domain_element_type& de,
           RhsExpr const& rhs_expr )
    {
        de=this->project( rhs_expr );
    }


    template< typename Expr>
    domain_element_type derivate( Expr expr )
    {
        de = this->domainSpace()->element();
        ie = this->backend()->newVector( this->dualImageSpace() );

        form1( _test=this->dualImageSpace(), _vector=ie, _init=true );

        form1( _test=this->dualImageSpace(), _vector=ie ) +=
            integrate(_range = elements( this->dualImageSpace()->mesh() ),
                      _expr = - trace( expr * trans( grad( this->dualImageSpace()->element() ) ) ) );

        form1(_test=this->dualImageSpace(), _vector=ie) +=
            integrate(_range = boundaryfaces( this->dualImageSpace()->mesh() ),
                      _expr =  trans( id( this->dualImageSpace()->element() ) ) * expr * vf::N() );

        ie->close();

        this->backend()->solve( M_matrixFull, de, ie );
        return de;
    }


    //@}


private :

    void initMatrix()
    {
        using namespace vf;
        auto uDomain = this->domainSpace()->element();
        auto uImage = this->dualImageSpace()->element();

        auto a = form2 ( _trial=this->domainSpace(),
                         _test=this->dualImageSpace(),
                         _matrix=M_matrixCst );

        if ( M_proj_type != LIFT )
        {
            a = integrate( _range=elements( this->dualImageSpace()->mesh() ),
                           _expr=inner( idt( uDomain ), id( uImage ) ) );
        }
        switch ( M_proj_type )
        {
        case L2:
        {
        }
        break;

        case H1:
        {
            a += integrate( _range=elements( this->dualImageSpace()->mesh() ),
                            _expr=inner( gradt( uDomain ), grad( uImage ) ) );
        }
        break;

        case DIFF:
        {
            a += integrate( _range=elements( this->dualImageSpace()->mesh() ),
                            _expr=M_epsilon*inner( gradt( uDomain ), grad( uImage ) ) );

            //weak boundary conditions
            a += integrate( _range=boundaryfaces( this->dualImageSpace()->mesh() ),
                            _expr= M_epsilon*( -trans( id( uImage ) )*gradt( uDomain )*vf::N() ) );
            a += integrate( _range=boundaryfaces( this->dualImageSpace()->mesh() ),
                            _expr= M_epsilon*( -trans( idt( uDomain ) )* grad( uImage )*vf::N() ) );
            a += integrate( _range=boundaryfaces( this->dualImageSpace()->mesh() ),
                            _expr= M_epsilon*( M_gamma * trans( idt( uDomain ) ) /*trial*/
                                        *id( uImage ) / vf::hFace()   /*test*/
                                        ) );
        }
        break;

        case HDIV:
        {
            a += integrate( _range=elements( this->dualImageSpace()->mesh() ),
                            _expr=divt( uDomain )*div( uImage ) );
        }
        break;

        case HCURL:
        {
            a += integrate( _range=elements( this->dualImageSpace()->mesh() ),
                           // only for 2D, need to specialize this for 3D
                           _expr=curlzt( uDomain )*curlz( uImage ) );
        }
        break;

        case LIFT:
        {
            a = integrate( _range=elements( this->dualImageSpace()->mesh() ),
                           _expr= inner( gradt( uDomain ), grad( uImage ) ) );
        }
        break;

        case CIP:
        {
            a += integrate( _range=internalfaces( this->dualImageSpace()->mesh() ),
                            _expr=M_gamma*hFace()*hFace()*
                            /**/  inner( jumpt( gradt(uDomain) ),jump( grad(uImage) ) ) );
        }
        break;

        case NODAL:
            break;
        }

        M_matrixCst->close();
    }

    template<typename Range, typename Expr>
    void applyOn( Range range, Expr expr, mpl::int_<MESH_ELEMENTS> ){}

    template<typename Range, typename Expr>
    void applyOn( Range range, Expr expr, mpl::int_<MESH_FACES> )
    {
        form2 ( _trial=this->domainSpace(),
                _test=this->dualImageSpace(),
                _matrix=M_matrixFull ) +=  on( _range=range , _element=de, _rhs=ie, _expr=expr );
    }

    //backend_ptrtype M_backend;
    const double M_epsilon;
    const double M_gamma;
    const ProjectorType M_proj_type;
    DirichletType M_dir;
    matrix_ptrtype M_matrixCst;
    matrix_ptrtype M_matrixFull;
    domain_element_type de;
    vector_ptrtype ie;

};//Projector

#if 1
/**
 * this function returns a \c Projector \c shared_ptr with
 *
 * \param domainSpace
 * \param imageSpace
 * \param backend
 */

template<typename TDomainSpace, typename TDualImageSpace>
boost::shared_ptr< Projector<TDomainSpace, TDualImageSpace> >
projector( boost::shared_ptr<TDomainSpace> const& domainspace,
           boost::shared_ptr<TDualImageSpace> const& imagespace,
           typename Projector<TDomainSpace, TDualImageSpace>::backend_ptrtype const& abackend = Backend<double>::build( soption( _name="backend" ) ),
           ProjectorType proj_type=L2, double epsilon=0.01, double gamma = 20, DirichletType dirichlet_type = WEAK)
{
    typedef Projector<TDomainSpace, TDualImageSpace > Proj_type;
    boost::shared_ptr<Proj_type> proj( new Proj_type( domainspace, imagespace, abackend, proj_type, epsilon, gamma, dirichlet_type ) );
    return proj;
}

BOOST_PARAMETER_FUNCTION( ( typename Feel::detail::projector_args<Args>::return_type ),
                          opProjection,
                          tag,
                          ( required
                            ( domainSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                            ( imageSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                          )
                          ( optional
                            ( type, (ProjectorType), L2 )
                            ( penaldir, *( boost::is_arithmetic<mpl::_> ), 20. )
                            ( backend, *, Backend<double>::build( soption( _name="backend" ) ) )
                          ) )
{
    return projector( domainSpace,imageSpace, backend, type, 0.01, penaldir );
}

BOOST_PARAMETER_FUNCTION( ( typename Feel::detail::lift_args<Args>::lift_return_type ),
                          opLift,
                          tag,
                          ( required
                            ( domainSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                          )
                          ( optional
                            ( type, (DirichletType), WEAK )
                            ( penaldir, *( boost::is_arithmetic<mpl::_> ), 20. )
                            ( backend, *, Backend<double>::build( soption( _name="backend" ) ) )
                            ) )
{
    return projector( domainSpace, domainSpace, backend, ProjectorType::LIFT, 0.01 , penaldir, type );
}


#endif



} //namespace Feel


#endif
