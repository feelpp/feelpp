/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
 Date: 2020-04-29

 Copyright (C) 2018 Université de Strasbourg
 Copyright (C) 2020 Inria

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
 \file levelsetspacemanager.hpp
 \author Thibaut Metivet <thibaut.metivet@inria.fr>
 \date 2020-04-29
 */

#include <feel/feelmodels/levelset/levelsetspacemanager.hpp>

namespace Feel {
namespace FeelModels {

#define LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS \
    template< typename ConvexType, typename BasisType, typename PeriodicityType, typename BasisPnType > \
        /**/
#define LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE \
    LevelSetSpaceManager<ConvexType, BasisType, PeriodicityType, BasisPnType> \
        /**/

LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::LevelSetSpaceManager( 
        mesh_ptrtype const& mesh,
        std::string const& prefix,
        std::string const& rootRepository
        ) :
    M_prefix( prefix ),
    M_rootRepository( rootRepository ),
    M_mesh( mesh ),
    M_worldsComm( worldscomm_ptrtype(1,mesh->worldCommPtr()) ),
    M_buildExtendedDofTable( false ),
    M_functionSpaceCreated( false )
{
    M_rangeMeshElements = elements( M_mesh );
}

LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::setPeriodicity( periodicity_type const& p ) 
{
    if( M_functionSpaceCreated )
        Feel::cout << "WARNING !! Setting periodicity after spaces creation ! Need to re-build them !" << std::endl;
    M_periodicity = p;
}

LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::createFunctionSpaceDefault()
{
    if( !M_spaceVectorial )
    {
        M_spaceVectorial = space_vectorial_type::New( 
                _mesh=this->mesh(), 
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity()
                );
    }
    if( !M_spaceScalar )
    {
        if( M_buildExtendedDofTable )
        {
            std::vector<bool> extendedDT( 1, M_buildExtendedDofTable );
            M_spaceScalar = space_scalar_type::New( 
                    _mesh=this->mesh(), 
                    _worldscomm=this->worldsComm(),
                    _extended_doftable=extendedDT,
                    _periodicity=this->periodicity()
                    );
        }
        else
        {
            M_spaceScalar = M_spaceVectorial->compSpace();
        }
    }
    if( !M_spaceMarkers )
    {
        M_spaceMarkers = space_markers_type::New( 
                _mesh=this->mesh(), 
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity(),
                _extended_doftable=std::vector<bool>(1, true)
                );
    }

    M_functionSpaceCreated = true;
}

LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::createFunctionSpaceIsoPN()
{
    if( !M_spaceScalarPN )
    {
        M_spaceScalarPN = space_scalar_PN_type::New( 
                _mesh=this->mesh(),
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity()
                );
    }
    if( !M_spaceVectorialPN )
    {
        M_spaceVectorialPN = space_vectorial_PN_type::New( 
                _mesh=this->mesh(),
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity()
                );
    }
    if( !M_meshIsoPN )
    {
        M_opLagrangeP1isoPN = lagrangeP1(
                _space=this->M_spaceScalarPN
                );
        M_meshIsoPN = M_opLagrangeP1isoPN->mesh();
        M_rangeMeshIsoPNElements = elements( M_meshIsoPN );
    }
    if( !M_spaceVectorialIsoPN )
    {
        M_spaceVectorialIsoPN = space_vectorial_type::New( 
                _mesh=this->meshIsoPN(),
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity()
                );
    }
    if( !M_spaceScalarIsoPN )
    {
        if( M_buildExtendedDofTable )
        {
            std::vector<bool> extendedDT( 1, M_buildExtendedDofTable );
            M_spaceScalarIsoPN = space_scalar_type::New(
                    _mesh=this->meshIsoPN(),
                    _worldscomm=this->worldsComm(),
                    _extended_doftable=extendedDT,
                    _periodicity=this->periodicity()
                    );
        }
        else
        {
            M_spaceScalarIsoPN = M_spaceVectorialIsoPN->compSpace();
        }
    }
    if( !M_spaceMarkersIsoPN )
    {
        M_spaceMarkersIsoPN = space_markers_type::New( 
                _mesh=this->meshIsoPN(),
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity(),
                _extended_doftable=std::vector<bool>(1, true)
                );
    }

    if( !M_opInterpolationScalarFromPN )
    {
        M_opInterpolationScalarFromPN = opInterpolation(
                _domainSpace = M_spaceScalarPN,
                _imageSpace = M_spaceScalarIsoPN,
                _type = InterpolationNonConforme(false)
                );
    }
    if( !M_opInterpolationScalarToPN )
    {
        M_opInterpolationScalarToPN = opInterpolation(
                _domainSpace = M_spaceScalarIsoPN,
                _imageSpace = M_spaceScalarPN,
                _type = InterpolationNonConforme(false)
                );
    }
    if( !M_opInterpolationVectorialFromPN )
    {
        M_opInterpolationVectorialFromPN = opInterpolation(
                _domainSpace = M_spaceVectorialPN,
                _imageSpace = M_spaceVectorialIsoPN,
                _type = InterpolationNonConforme(false)
                );
    }
    if( !M_opInterpolationVectorialToPN )
    {
        M_opInterpolationVectorialToPN = opInterpolation(
                _domainSpace = M_spaceVectorialIsoPN,
                _imageSpace = M_spaceVectorialPN,
                _type = InterpolationNonConforme(false)
                );
    }

    M_functionSpaceCreated = true;
}

LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::createFunctionSpaceHovisu()
{
    CHECK( M_spaceScalar ) << "M_spaceScalar must be built to build Hovisu function spaces\n";
    CHECK( M_spaceVectorial ) << "M_spaceVectorial must be built to build Hovisu function spaces\n";

    if( !M_meshHovisu )
    {
        M_opLagrangeP1Hovisu = lagrangeP1(
                _space=this->M_spaceScalar,
                _path=this->M_rootRepository,
                _prefix=prefixvm( this->M_prefix, "hovisu" ),
                //_rebuild=!this->doRestart(),
                _parallel=false
                );
        M_meshHovisu = M_opLagrangeP1Hovisu->mesh();
        M_rangeMeshHovisuElements = elements( M_meshHovisu );
    }
    if( !M_spaceVectorialHovisu )
    {
        M_spaceVectorialHovisu = space_vectorial_hovisu_type::New( 
                _mesh=this->meshHovisu(),
                _worldscomm=this->worldsComm()
                );
    }
    if( !M_spaceScalarHovisu )
    {
        M_spaceScalarHovisu = M_spaceVectorialHovisu->compSpace();
    }

    if( !M_opInterpolationScalarToHovisu )
    {
        M_opInterpolationScalarToHovisu = opInterpolation(
                _domainSpace = M_spaceScalar,
                _imageSpace = M_spaceScalarHovisu,
                _type = InterpolationNonConforme(false, true, false, 15)
                );
    }
    if( !M_opInterpolationVectorialToHovisu )
    {
        M_opInterpolationVectorialToHovisu = opInterpolation(
                _domainSpace = M_spaceVectorial,
                _imageSpace = M_spaceVectorialHovisu,
                _type = InterpolationNonConforme(false, true, false, 15)
                );
    }

    M_functionSpaceCreated = true;
}


LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE::createFunctionSpaceTensor2Symm()
{
    if( !M_spaceTensor2Symm )
    {
        M_spaceTensor2Symm = space_tensor2symm_type::New( 
                _mesh=this->mesh(), 
                _worldscomm=this->worldsComm(),
                _periodicity=this->periodicity()
                );
    }

    M_functionSpaceCreated = true;
}

#undef LEVELSETSPACEMANAGER_CLASS_TEMPLATE_DECLARATIONS
#undef LEVELSETSPACEMANAGER_CLASS_TEMPLATE_TYPE

} // namespace FeelModels
} // namespace Feel

