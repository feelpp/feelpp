/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Winkelmann
              Vincent Doyeux
              Thibaut Metivet
   Date     : Fri Mar 30 15:18:07 2018

   Copyright (C) 2014-2018 Feel++ Consortium

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
#ifndef FEELPP_REINIT_FMS_IMPL_HPP
#define FEELPP_REINIT_FMS_IMPL_HPP 1

#define FEELPP_INSTANTIATE_FMS 1

#include <feel/feells/reinit_fms.hpp>

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/operators.hpp>

//#define FM_EXPORT 1

#if defined(FM_EXPORT)
#include <feel/feelfilters/exporter.hpp>
#endif

#define REINITIALIZERFMS_CLASS_TEMPLATE_DECLARATIONS \
    template<typename FunctionSpaceType, typename periodicity_type> \
    /**/ 
#define REINITIALIZERFMS_CLASS_TEMPLATE_TYPE \
    ReinitializerFMS<FunctionSpaceType, periodicity_type> \
    /**/ 

namespace Feel
{

REINITIALIZERFMS_CLASS_TEMPLATE_DECLARATIONS
REINITIALIZERFMS_CLASS_TEMPLATE_TYPE::ReinitializerFMS( 
        functionspace_ptrtype const& __functionspace,
        periodicity_type __periodicity)
:
    super_type( __functionspace, __periodicity )
{}

REINITIALIZERFMS_CLASS_TEMPLATE_DECLARATIONS
void 
REINITIALIZERFMS_CLASS_TEMPLATE_TYPE::processDof( size_type idOnProc, value_type val, heap_data_type const& opt_data  )
{
    this->setDofDistance( idOnProc, val );
    this->setDofStatus( idOnProc, super_type::DONE );
}


} // Feel

#undef REINITIALIZERFMS_CLASS_TEMPLATE_TYPE 
#undef REINITIALIZERFMS_CLASS_TEMPLATE_DECLARATIONS 

#endif
