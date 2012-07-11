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
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{

/**
 * \class OperatorTrace
 * \brief OperatorTrace made easy
 *
 * @author Abdoulaye Samake
 * @see OperatorLinear
 */
template<class fs_type>
class OperatorTrace
{
    typedef OperatorTrace<fs_type> super;

public :

    /** @name Typedefs
     */
    //@{

    typedef fs_type functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::trace_functionspace_type trace_functionspace_type;
    typedef typename trace_functionspace_type::element_type trace_element_type;

    //@}
    /** @name Constructors, destructor
     */
    //@{

    OperatorTrace( functionspace_ptrtype domainSpace )
        :
        M_domainSpace( domainSpace )
    {
    }

    ~OperatorTrace() {}
    //@}

    /** @name  Methods
     */
    //@{

    BOOST_PARAMETER_MEMBER_FUNCTION( ( trace_element_type ),
                                     trace,
                                     tag,
                                     ( required
                                       ( expr,   * )
                                     )
                                     ( optional
                                       ( range,  *, boundaryfaces( M_domainSpace->mesh() ) )
                                     ) )

    {
        using namespace vf;


        auto Th = M_domainSpace->trace( range ) ;

        trace_element_type te = Th->element();

        te = vf::project( _space=Th, _range=elements( Th->mesh() ), _expr=expr );

        return te;
    }

    //@}


private :

    functionspace_ptrtype M_domainSpace;
};//OperatorTrace

/**
 * this function returns a \c OperatorTrace \c shared_ptr with
 *
 * \param traceSpace
 *
 */

template<typename self_type>
boost::shared_ptr< OperatorTrace<self_type> >
operatorTrace( boost::shared_ptr<self_type> const& domainspace )
{
    typedef OperatorTrace<self_type> Trace_type;
    boost::shared_ptr<Trace_type> trace( new Trace_type( domainspace ) );
    return trace;
}

} //namespace Feel


#endif
