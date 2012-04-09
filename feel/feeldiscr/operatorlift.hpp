/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
       Date: 2011-08-07

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
   \file operatorlift.cpp
   \author Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
   \date 2011-08-07
 */
#ifndef _OPERATORLIFT_HPP_
#define _OPERATORLIFT_HPP_

#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelvf/vf.hpp>
#include<iostream>
namespace Feel
{

/* use weak or strong Dirichlet condition */
// enum WeakDirichlet{ strong=0, weak=1 };
/**
 * \class OperatorLift
 * \brief OperatorLift made easy
 *
 * @author Abdoulaye Samake
 * @see OperatorLinear
 */
template<class DomainSpace>
class OperatorLift : public OperatorLinear<DomainSpace, DomainSpace>
{
    typedef OperatorLift<DomainSpace> super;

public :

    /** @name Typedefs
     */
    //@{

    typedef OperatorLinear<DomainSpace, DomainSpace> ol_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::element_type domain_element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DomainSpace> image_element_type;

    //@}
    /** @name Constructors, destructor
     */
    //@{

    OperatorLift( domain_space_ptrtype domainSpace,
                  backend_ptrtype backend = Backend<double>::build( BACKEND_PETSC ),
                  double gamma = 20,
                  WeakDirichlet dirichlet_type = WEAK
                )
        :
        ol_type( domainSpace, domainSpace, backend ),
        M_backend( backend ),
        M_gamma( gamma ),
        M_dir( dirichlet_type ),
        M_matrix( M_backend->newMatrix( domainSpace, domainSpace ) )
    {
        initMatrix();
    }

    ~OperatorLift() {}
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

    BOOST_PARAMETER_MEMBER_FUNCTION( ( domain_element_type ),
                                     lift,
                                     tag,
                                     ( required
                                       ( range,  * )
                                       ( expr,   * )
                                     )
                                     ( optional
                                       ( quad,   *, ( typename integrate_type<Args,decltype( elements( this->domainSpace()->mesh() ) )>::_quad_type() ) )
                                       ( quad1,   *, ( typename integrate_type<Args,decltype( elements( this->domainSpace()->mesh() ) )>::_quad1_type() ) )
                                       ( geomap, *, GeomapStrategyType::GEOMAP_OPT )
                                     ) )

    {
        using namespace vf;

        domain_element_type de = this->domainSpace()->element();

        auto ie = M_backend->newVector( this->domainSpace() );

        //weak dirichlet boundary conditions

        form1( _test=this->domainSpace(), _vector=ie, _init=true );

        if ( M_dir == WEAK )
        {

            form1( _test=this->domainSpace(), _vector=ie ) +=
                integrate( _range=range,
                           _expr=expr*( -grad( this->domainSpace()->element() )*vf::N() +
                                        M_gamma / vf::hFace() *id( this->domainSpace()->element() ) ),
                           _quad=quad );
        }

        M_matrixFull = M_backend->newMatrix( this->domainSpace(), this->domainSpace() );

        form2( _test=this->domainSpace(), _trial=this->domainSpace(), _matrix=M_matrixFull, _init=true );

        if ( M_dir == WEAK )
        {

            form2 ( _trial=this->domainSpace(),
                    _test=this->domainSpace(),
                    _matrix=M_matrixFull ) +=
                        integrate( _range=range, _expr=
                                       ( -trans( id( this->domainSpace()->element() ) )*gradt( this->domainSpace()->element() )*vf::N()
                                         -trans( idt( this->domainSpace()->element() ) )* grad( this->domainSpace()->element() )*vf::N()
                                         + M_gamma * trans( idt( this->domainSpace()->element() ) ) /*trial*/
                                         *id( this->domainSpace()->element() ) / vf::hFace()   /*test*/
                                       ) );

        }

        M_matrixFull->close();
        M_matrixFull->addMatrix( 1., M_matrix );

        if ( M_dir == STRONG )
        {
            form2 ( _trial=this->domainSpace(),
                    _test=this->domainSpace(),
                    _matrix=M_matrixFull ) +=  on( range , de, ie, expr );
        }

        M_backend->solve( M_matrixFull, de, ie );

        return de;
    }

    template<typename Elts, typename RhsExpr>
    domain_element_type
    operator()( Elts elts ,RhsExpr const& rhs_expr )
    {
        return this->lift( elts, rhs_expr );
    }



    template<typename Elts, typename RhsExpr>
    void
    operator()( Elts elts, domain_element_type& de, RhsExpr const& rhs_expr )
    {
        de = this->lift( elts, rhs_expr );
    }

    //@}


private :

    void initMatrix()
    {
        using namespace vf;


        form2 ( _trial=this->domainSpace(),
                _test=this->domainSpace(),
                _matrix=M_matrix,
                _init=true ) =
                    integrate( elements( this->domainSpace()->mesh() ),
                               trace( gradt( this->domainSpace()->element() )
                                      * trans( grad( this->domainSpace()->element() ) ) )
                             );

        M_matrix->close();

    }

    backend_ptrtype M_backend;
    matrix_ptrtype M_matrix;
    matrix_ptrtype M_matrixFull;
    const double M_gamma;
    WeakDirichlet M_dir;

};//OperatorLift

/**
 * this function returns a \c OperatorLift \c shared_ptr with
 *
 * \param domainSpace
 * \param backend
 */

template<typename TDomainSpace>
boost::shared_ptr< OperatorLift<TDomainSpace> >
operatorLift( boost::shared_ptr<TDomainSpace> const& domainspace,
              typename OperatorLift<TDomainSpace>::backend_ptrtype const& backend = Backend<double>::build( BACKEND_PETSC ),
              double gamma = 20, WeakDirichlet dirichlet_type = WEAK )
{
    typedef OperatorLift<TDomainSpace> Proj_type;
    boost::shared_ptr<Proj_type> Lift( new Proj_type( domainspace, backend, gamma, dirichlet_type ) );
    return Lift;
}

} //namespace Feel


#endif
