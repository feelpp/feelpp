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
#include <feel/feelvf/vf.hpp>

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
               backend_ptrtype abackend = Feel::backend(_rebuild=true),
               ProjectorType proj_type=L2,
               double epsilon = 0.01,
               double gamma = 20,
               DirichletType dirichlet_type = WEAK
             )
        :
        ol_type( domainSpace, dualImageSpace, abackend ),
        M_backend( abackend ),
        M_epsilon( epsilon ),
        M_gamma( gamma ),
        M_proj_type( proj_type ),
        M_dir( dirichlet_type )

    {
        M_matrix = M_backend->newMatrix( _trial=domainSpace, _test=dualImageSpace ) ;
        initMatrix();
    }

    ~Projector() {}
    //@}

    /** @name  Methods
     */
    //@{
    template<typename Args,typename IntEltsDefault>
    struct integrate_type
    {
        typedef typename vf::detail::clean_type<Args,tag::expr>::type _expr_type;
        typedef typename vf::detail::clean2_type<Args,tag::range,IntEltsDefault>::type _range_type;
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
                                       ( geomap, *, GeomapStrategyType::GEOMAP_OPT )
                                       (grad_expr, *, ( zero<domain_space_type::nComponents,domain_space_type::nDim>() ))
                                       (div_expr, *, cst(0.) )
                                       (curl_expr, *,  ( zero<  mpl::if_<mpl::equal_to<mpl::int_<domain_space_type::nComponents>, mpl::int_<1> >, mpl::int_<1>, mpl::int_<domain_space_type::nDim> >::type::value, 1>() ) )
                                       )
                                   )
    {
        using namespace vf;

        auto sol = this->domainSpace()->element();

        ie = M_backend->newVector( this->dualImageSpace() );

        form1( _test=this->dualImageSpace(), _vector=ie, _init=true );

        if ( (M_proj_type != LIFT) )
        {
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=range, _expr=expr*id( this->dualImageSpace()->element() ),
                           _quad=quad, _quad1=quad1, _geomap=geomap );

            switch( M_proj_type )
                {
                case H1:
                    form1( _test=this->dualImageSpace(), _vector=ie ) +=
                        integrate( _range=range, _expr=trace(grad_expr*trans(grad( this->dualImageSpace()->element() )) ),
                                   _quad=quad, _quad1=quad1, _geomap=geomap );
                    break;
                case HDIV:
                            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                                integrate( _range=range, _expr=div_expr*div( this->dualImageSpace()->element() ),
                                           _quad=quad, _quad1=quad1, _geomap=geomap );
                            break;
                case HCURL:
                            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                                integrate( _range=range, _expr=trans(curl_expr)*curl( this->dualImageSpace()->element() ),
                                           _quad=quad, _quad1=quad1, _geomap=geomap );
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
                           _expr=expr*( -grad( this->dualImageSpace()->element() )*vf::N() +
                                        M_gamma / vf::hFace() *id( this->dualImageSpace()->element() ) ),
                           _quad=quad, _quad1=quad1, _geomap=geomap );
        }

        //weak boundary conditions
        if ( M_proj_type == DIFF )
        {
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=boundaryfaces( this->dualImageSpace()->mesh() ),
                           _expr=expr*M_epsilon*( -grad( this->domainSpace()->element() )*vf::N() +
                                                  M_gamma / vf::hFace() *id( this->dualImageSpace()->element() ) ),
                           _quad=quad, _quad1=quad1, _geomap=geomap );
        }

        ie->close();

        M_matrixFull = M_backend->newMatrix( _trial=this->domainSpace(), _test=this->dualImageSpace() );
        auto bilinearForm = form2( _trial=this->domainSpace(), _test=this->dualImageSpace(), _matrix=M_matrixFull );

        if ( ( M_proj_type == LIFT ) && ( M_dir == WEAK ) )
        {
            bilinearForm +=
                integrate( _range=range, _expr=
                           ( -trans( id( this->dualImageSpace()->element() ) )*gradt( this->domainSpace()->element() )*vf::N()
                             -trans( idt( this->domainSpace()->element() ) )* grad( this->dualImageSpace()->element() )*vf::N()
                             + M_gamma * trans( idt( this->domainSpace()->element() ) ) /*trial*/
                             *id( this->dualImageSpace()->element() ) / vf::hFace()   /*test*/
                             ), _quad=quad, _quad1=quad1, _geomap=geomap );
        }

        M_matrixFull->close();
        M_matrixFull->addMatrix( 1., M_matrix );

        if ( ( M_proj_type == LIFT ) && ( M_dir == STRONG )  )
            this->applyOn(range, expr);

        M_backend->solve( M_matrixFull, sol, ie );

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
        M_backend->solve( M_matrix, de, ie );
        return de ;
    }

    void
    operator()( domain_element_type &de, image_element_type const& ie )
    {
        M_backend->solve( M_matrix, de, ie );
    }

    template<typename Range, typename Expr>
    domain_element_type
    operator()( Range const& range ,Expr const& expr )
    {
        return this->project( expr, range );
    }

    void
    apply( domain_element_type& de,
           image_element_type const& ie )
    {
        M_backend->solve( M_matrix, de, ie );
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
        ie = M_backend->newVector( this->dualImageSpace() );

        form1( _test=this->dualImageSpace(), _vector=ie, _init=true );

        form1( _test=this->dualImageSpace(), _vector=ie ) +=
            integrate(_range = elements( this->dualImageSpace()->mesh() ),
                      _expr = - trace( expr * trans( grad( this->dualImageSpace()->element() ) ) ) );

        form1(_test=this->dualImageSpace(), _vector=ie) +=
            integrate(_range = boundaryfaces( this->dualImageSpace()->mesh() ),
                      _expr =  trans( id( this->dualImageSpace()->element() ) ) * expr * vf::N() );

        ie->close();

        M_backend->solve( M_matrix, de, ie );
        return de;
    }


    //@}


private :

    void initMatrix()
    {
        using namespace vf;
        auto a = form2 ( _trial=this->domainSpace(),
                         _test=this->dualImageSpace(),
                         _matrix=M_matrix,
                         _init=true );

        switch ( M_proj_type )
        {
        case L2:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                         );
        }
        break;

        case H1:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                           +
                           trace( gradt( this->domainSpace()->element() )
                                  * trans( grad( this->dualImageSpace()->element() ) ) )
                           );
        }
        break;

        case DIFF:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                           +
                           M_epsilon *
                           trace( gradt( this->domainSpace()->element() )
                                  * trans( grad( this->dualImageSpace()->element() ) ) )
                           );
            //weak boundary conditions
            a += integrate( boundaryfaces( this->dualImageSpace()->mesh() ),
                            M_epsilon*( -trans( id( this->dualImageSpace()->element() ) )*gradt( this->domainSpace()->element() )*vf::N() ) );
            a += integrate( boundaryfaces( this->dualImageSpace()->mesh() ),
                            M_epsilon*( -trans( idt( this->domainSpace()->element() ) )* grad( this->dualImageSpace()->element() )*vf::N() ) );
            a += integrate( boundaryfaces( this->dualImageSpace()->mesh() ),
                            M_epsilon*( M_gamma * trans( idt( this->domainSpace()->element() ) ) /*trial*/
                                        *id( this->dualImageSpace()->element() ) / vf::hFace()   /*test*/
                                        ) );
        }
        break;

        case HDIV:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                           +
                           ( divt( this->domainSpace()->element() ) *
                             div( this->dualImageSpace()->element() ) )
                           );
        }
        break;

        case HCURL:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                           +
                           // only for 2D, need to specialize this for 3D
                           curlzt( this->domainSpace()->element() )
                           * curlz( this->dualImageSpace()->element() )
                         );
        }
        break;

        case LIFT:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trace( gradt( this->domainSpace()->element() )
                                  * trans( grad( this->dualImageSpace()->element() ) ) )
                           );

        }
        break;

        case CIP:
        {
            a = integrate( elements( this->dualImageSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->dualImageSpace()->element() ) /*test*/
                           );

            a += integrate( internalfaces( this->dualImageSpace()->mesh() ),
                            M_gamma * hFace() * hFace()
                           * trans(jumpt( gradt(this->domainSpace()->element()) ))
                           * jump( grad(this->dualImageSpace()->element()) )
                           );
        }
        break;

        case NODAL:
            break;
        }

        M_matrix->close();
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

    backend_ptrtype M_backend;
    const double M_epsilon;
    const double M_gamma;
    const ProjectorType M_proj_type;
    DirichletType M_dir;
    matrix_ptrtype M_matrix;
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
           typename Projector<TDomainSpace, TDualImageSpace>::backend_ptrtype const& abackend = Feel::backend(_rebuild=true),
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
                            ( backend, *, Feel::backend(_rebuild=true) )
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
                            ( backend, *, backend(_rebuild=true) )
                            ) )
{
    return projector( domainSpace, domainSpace, backend, ProjectorType::LIFT, 0.01 , penaldir, type );
}


#endif



} //namespace Feel


#endif
