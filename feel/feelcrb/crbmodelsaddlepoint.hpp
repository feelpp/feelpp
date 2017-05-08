/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-09

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
 * @file   crbmodelsaddlepoint.hpp
 * @author wahl
 * @date   Tue Nov 25 16:24:55 2014
 */


#ifndef __CRBModelSaddlePoint_H
#define __CRBModelSaddlePoint_H 1

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbmodel.hpp>

namespace Feel
{
template<typename ModelType>
class CRBModelSaddlePoint :
        public CRBModel<ModelType>
{
    typedef CRBModel<ModelType> super;
public :
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;


    CRBModelSaddlePoint( CRBModelMode mode = CRBModelMode::PFEM, int level=0, bool doInit = true ) :
        super ( mode, level, doInit )
    {}

    CRBModelSaddlePoint( model_ptrtype const& model , CRBModelMode mode = CRBModelMode::PFEM, bool doInit = true ) :
        super ( model, mode, doInit )
    {}

protected :


}; // class CRBModelSaddlepoint


} // namespace Feel

#endif // __CRBModelSaddlePoint_H
