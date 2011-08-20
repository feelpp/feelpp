/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
       Date: 2011-08-18

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
   \file operatortrace.cpp
   \author Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
   \date 2011-08-18
 */
#ifndef _OPERATORTRACE_HPP_
#define _OPERATORTRACE_HPP_

#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelvf/vf.hpp>
#include<iostream>
namespace Feel
{

/* use weak or strong Dirichlet condition */
    // enum WeakDirichlet{ strong=0, weak=1 };
/**
 * \class OperatorTrace
 * \brief OperatorTrace made easy
 *
 * @author Abdoulaye Samake
 * @see OperatorLinear
 */
template<class TraceSpace>
class OperatorTrace : public OperatorLinear<TraceSpace, TraceSpace>
{
    typedef OperatorTrace<TraceSpace> super;

public :

    /** @name Typedefs
     */
    //@{

    typedef OperatorLinear<TraceSpace, TraceSpace> ol_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename domain_space_type::element_type trace_element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<TraceSpace> image_element_type;

    //@}
    /** @name Constructors, destructor
     */
    //@{

    OperatorTrace(domain_space_ptrtype traceSpace,
                  backend_ptrtype backend = Backend<double>::build(BACKEND_PETSC))
        :
        ol_type(traceSpace, traceSpace, backend),
        M_backend(backend)
    {
    }

    ~OperatorTrace() {}
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

    BOOST_PARAMETER_MEMBER_FUNCTION((trace_element_type),
                                    trace,
                                    tag,
                                    (required
                                     (range,  *)
                                     (expr,   *)
                                     ) )

    {
        using namespace vf;

        trace_element_type te = this->domainSpace()->element();

        te = vf::project(_space=this->domainSpace(), _range=range, _expr=expr );

        return te;
    }

    //@}


private :

    backend_ptrtype M_backend;

};//OperatorTrace

/**
 * this function returns a \c OperatorTrace \c shared_ptr with
 *
 * \param traceSpace
 * \param backend
 */

template<typename TTraceSpace>
boost::shared_ptr< OperatorTrace<TTraceSpace> >
operatorTrace( boost::shared_ptr<TTraceSpace> const& tracespace,
               typename OperatorTrace<TTraceSpace>::backend_ptrtype const& backend = Backend<double>::build(BACKEND_PETSC) )
{
    typedef OperatorTrace<TTraceSpace> Trace_type;
    boost::shared_ptr<Trace_type> Trace( new Trace_type(tracespace, backend ) );
    return Trace;
}

} //namespace Feel


#endif
