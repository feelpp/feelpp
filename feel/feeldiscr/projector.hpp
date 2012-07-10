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
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
template<class DomainSpace, class DualImageSpace> class Projector;
    enum ProjectorType {L2=0, H1=1, DIFF=2, HDIV=3, HCURL=4, LIFT=5};
namespace detail
{
template<typename Args>
struct projector_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::type::value_type domain_type;
    typedef typename vf::detail::clean_type<Args,tag::imageSpace>::type::value_type image_type;
    //typedef typename vf::detail::clean_type<Args,tag::type>::type type;
    //typedef typename vf::detail::clean_type<Args,tag::backend>::backend_type backend_type;

    typedef boost::shared_ptr<Projector<domain_type,image_type> > return_type;
    typedef Projector<domain_type,image_type> projector_type;
    typedef boost::shared_ptr<Projector<domain_type,domain_type> > lift_return_type;
};

template<typename Args>
struct lift_args
{
    typedef typename vf::detail::clean_type<Args,tag::domainSpace>::type::value_type domain_type;
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

    // typedef Operator<DomainSpace, DualImageSpace> super_type;
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
               backend_ptrtype backend = Backend<double>::build( BACKEND_PETSC ),
               ProjectorType proj_type=L2,
               double epsilon = 0.01,
               double gamma = 20,
               WeakDirichlet dirichlet_type = WEAK
             )
        :
        ol_type( domainSpace, dualImageSpace, backend ),
        M_backend( backend ),
        M_epsilon( epsilon ),
        M_gamma( gamma ),
        M_proj_type( proj_type ),
        M_dir( dirichlet_type ),
        M_matrix( M_backend->newMatrix( _trial=domainSpace, _test=dualImageSpace ) )
    {
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
        //typedef _Q< ExpressionOrder<_range_type,_expr_type>::value > the_quad_type;
        typedef typename vf::detail::clean2_type<Args,tag::quad, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value > >::type _quad_type;
        typedef typename vf::detail::clean2_type<Args,tag::quad1, _Q< vf::ExpressionOrder<_range_type,_expr_type>::value_1 > >::type _quad1_type;
    };

    BOOST_PARAMETER_MEMBER_FUNCTION( ( domain_element_type ),
                                     project,
                                     tag,
                                     ( required
                                       ( expr,   * ) )
                                     ( optional
                                       ( range,   *,( M_proj_type != LIFT ) ? elements( this->domainSpace()->mesh() ) : boundaryfaces(this->domainSpace()->mesh()) )
                                       ( quad,   *, ( typename integrate_type<Args,decltype( elements( this->domainSpace()->mesh() ) )>::_quad_type() ) )
                                       ( quad1,   *, ( typename integrate_type<Args,decltype( elements( this->domainSpace()->mesh() ) )>::_quad1_type() ) )
                                       ( geomap, *, GeomapStrategyType::GEOMAP_OPT )
                                     ) )
    {
        using namespace vf;
        domain_element_type de = this->domainSpace()->element();

        auto ie = M_backend->newVector( this->dualImageSpace() );
        form1( _test=this->dualImageSpace(), _vector=ie, _init=true );

        if ( M_proj_type != LIFT )
        {
            form1( _test=this->dualImageSpace(), _vector=ie ) +=
                integrate( _range=range, _expr=expr * id( this->dualImageSpace()->element() ),
                           _quad=quad, _quad1=quad1, _geomap=geomap );
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
                integrate( _range=boundaryfaces( this->domainSpace()->mesh() ),
                           _expr=expr*M_epsilon*( -grad( this->domainSpace()->element() )*vf::N() +
                                                  M_gamma / vf::hFace() *id( this->dualImageSpace()->element() ) ),
                           _quad=quad );
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
                             ) );
        }

        M_matrixFull->close();
        M_matrixFull->addMatrix( 1., M_matrix );

        if ( ( M_proj_type == LIFT ) && ( M_dir == STRONG ) )
        {
            form2 ( _trial=this->domainSpace(),
                    _test=this->dualImageSpace(),
                    _matrix=M_matrixFull ) +=  on( range , de, ie, expr );
        }

        M_backend->solve( M_matrixFull, de, ie );

        return de;
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
            a = integrate( elements( this->domainSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->domainSpace()->element() ) /*test*/
                         );
        }
        break;

        case H1:
        {
            a = integrate( elements( this->domainSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->domainSpace()->element() ) /*test*/
                           +
                           trace( gradt( this->domainSpace()->element() )
                                  * trans( grad( this->domainSpace()->element() ) ) )
                         );
        }
        break;

        case DIFF:
        {
            a = integrate( elements( this->domainSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->domainSpace()->element() ) /*test*/
                           +
                           M_epsilon *
                           trace( gradt( this->domainSpace()->element() )
                                  * trans( grad( this->domainSpace()->element() ) ) )
                         );
            //weak boundary conditions
            a += integrate( boundaryfaces( this->domainSpace()->mesh() ),
                            M_epsilon*( -trans( id( this->domainSpace()->element() ) )*gradt( this->domainSpace()->element() )*vf::N() ) );
            a += integrate( boundaryfaces( this->domainSpace()->mesh() ),
                            M_epsilon*( -trans( idt( this->domainSpace()->element() ) )* grad( this->domainSpace()->element() )*vf::N() ) );
            a += integrate( boundaryfaces( this->domainSpace()->mesh() ),
                            M_epsilon*( M_gamma * trans( idt( this->domainSpace()->element() ) ) /*trial*/
                                        *id( this->domainSpace()->element() ) / vf::hFace()   /*test*/
                                      ) );
        }
        break;

        case HDIV:
        {
            a = integrate( elements( this->domainSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->domainSpace()->element() ) /*test*/
                           +
                           ( divt( this->domainSpace()->element() ) *
                             div( this->domainSpace()->element() ) )
                         );
        }
        break;

        case HCURL:
        {
            a = integrate( elements( this->domainSpace()->mesh() ),
                           trans( idt( this->domainSpace()->element() ) ) /*trial*/
                           *id( this->domainSpace()->element() ) /*test*/
                           +
                           // only for 2D, need to specialize this for 3D
                           curlzt( this->domainSpace()->element() )
                           * curlz( this->domainSpace()->element() )
                         );
        }
        break;

        case LIFT:
            {
                a = integrate( elements( this->domainSpace()->mesh() ),
                               trace( gradt( this->domainSpace()->element() )
                                      * trans( grad( this->domainSpace()->element() ) ) )
                               );

        }
        break;
        }

        M_matrix->close();
    }

    backend_ptrtype M_backend;
    const double M_epsilon;
    const double M_gamma;
    ProjectorType M_proj_type;
    WeakDirichlet M_dir;
    matrix_ptrtype M_matrix;
    matrix_ptrtype M_matrixFull;

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
           typename Projector<TDomainSpace, TDualImageSpace>::backend_ptrtype const& backend = Backend<double>::build( BACKEND_PETSC ),
           ProjectorType proj_type=L2, double epsilon=0.01, double gamma = 20, WeakDirichlet dirichlet_type = WEAK )
{
    typedef Projector<TDomainSpace, TDualImageSpace> Proj_type;
    boost::shared_ptr<Proj_type> proj( new Proj_type( domainspace, imagespace, backend, proj_type, epsilon, gamma, dirichlet_type ) );
    return proj;
}

BOOST_PARAMETER_FUNCTION( ( typename detail::projector_args<Args>::return_type ),
                          opProjection,
                          tag,
                          ( required
                            ( domainSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                            ( imageSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                          )
                          ( optional
                            ( type, *, L2 )
                            ( backend, *, Backend<double>::build( BACKEND_PETSC ) )
                          ) )
{
    return projector( domainSpace,imageSpace, backend, type );
}

BOOST_PARAMETER_FUNCTION( ( typename detail::lift_args<Args>::lift_return_type ),
                          opLift,
                          tag,
                          ( required
                            ( domainSpace,   *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
                          )
                          ( optional
                            ( type, *, WeakDirichlet::WEAK )
                            ( backend, *, Backend<double>::build( BACKEND_PETSC ) )
                            ) )
{
    return projector( domainSpace, domainSpace, backend, ProjectorType::LIFT, 0.01 , 20, type );
}


#endif



} //namespace Feel


#endif
