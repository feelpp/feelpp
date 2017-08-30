/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-17

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file fluidmechanics.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelpde/preconditionerblockns.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::FluidMechanics( std::string const& prefix,
                                                    bool buildMesh,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    :
    super_type( prefix, buildMesh, worldComm, subPrefix, rootRepository )
{
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::self_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::New( std::string const& prefix, bool buildMesh,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type>( prefix, buildMesh, worldComm, subPrefix, rootRepository );

}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::init( bool buildModelAlgebraicFactory )
{
    super_type::init( buildModelAlgebraicFactory, this->shared_from_this() );


}



//---------------------------------------------------------------------------------------------------------//


} // namespace FeelModels

} // namespace Feel
