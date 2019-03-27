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
 \file meshale.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#include <feel/feelmodels/modelmesh/meshale.hpp>

namespace Feel
{
namespace FeelModels
{


template< class Convex >
MeshALE<Convex>::MeshALE(mesh_ptrtype mesh_moving,
                         std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                         ModelBaseRepository const& modelRep )
    :
    super_type( prefixvm(prefix,"alemesh"),worldcomm,"",modelRep ),
    M_referenceMesh( mesh_moving->createP1mesh() ),
    M_movingMesh(mesh_moving),
    M_isOnReferenceMesh( true ), M_isOnMovingMesh( true ),
    M_aleFactory( ale_map_type::build(M_referenceMesh ,this->prefix(), this->worldCommPtr(), modelRep ) ),
    M_Xhref( M_aleFactory->functionSpace() ),
    M_Xhmove( ale_map_functionspace_type::New(_mesh=M_movingMesh ) ),
    M_identity_ale(new ale_map_element_type( M_Xhmove) ),
    M_dispP1ToHO_ref(new ale_map_element_ref_type( M_Xhref) ),
    M_displacementOnMovingBoundary_HO_ref(new ale_map_element_type( M_Xhmove) ),
    M_displacementOnMovingBoundary_P1_ref(new ale_map_element_ref_type( M_Xhref) ),
    M_displacement(new ale_map_element_type( M_Xhmove) ),
    M_displacement_ref(new ale_map_element_ref_type( M_Xhref) ),
    M_meshVelocity(new ale_map_element_type( M_Xhmove) ),
    M_drm( new DofRelationshipMap_type(M_Xhref,M_Xhmove ) ),
    M_cpt_export(0),
    M_isARestart(boption(_name="ts.restart")),
    M_restartPath(soption(_name="ts.restart.path"))
    //M_doExport(option(_name="export",_prefix=prefixvm(prefix,"alemesh")).template as<bool>())
{

    this->log(prefixvm(this->prefix(),"MeshALE"),"constructor", "start");

    // update M_identity_ale
    this->updateIdentityMap();

    // compute dist between P1(ref) to Ho mesh
    if ( mesh_type::nOrder == mesh_ref_type::nOrder )
    {
        M_dispP1ToHO_ref->on(_range=elements(M_referenceMesh),
                             _expr=vf::P() );
        for (size_type i=0;i<M_dispP1ToHO_ref->nLocalDof();++i)
            (*M_dispP1ToHO_ref)(i) = (*M_identity_ale)(M_drm->dofRelMap()[i]) - (*M_dispP1ToHO_ref)(i);
    }

    this->initTimeStep();

    this->log(prefixvm(this->prefix(),"MeshALE"),"constructor", "finish");
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::init()
{
    this->log(prefixvm(this->prefix(),"MeshALE"),"init", "start");

    if (!M_isARestart)
    {
        *M_displacement_ref = *M_dispP1ToHO_ref;

        M_bdf_ale_identity->start(*M_identity_ale);
        M_bdf_ale_velocity->start(*M_meshVelocity);
        M_bdf_ale_displacement_ref->start(*M_displacement_ref);
    }
    else
    {
        M_bdf_ale_identity->restart();
        M_bdf_ale_velocity->restart();
        M_bdf_ale_displacement_ref->restart();

        // get displacement from t0
        *M_displacement_ref = M_bdf_ale_displacement_ref->unknown(0);

        // transfert displacement on the mobile mesh
        for (size_type i=0;i<M_displacement->nLocalDof();++i)
            (*M_displacement)(M_drm->dofRelMap()[i])=(*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i);

        //for (uint i=0;i<M_displacement->nDof();++i)
        //    (*M_displacement)(M_drm->dofRelMap()[i]) -= M_bdf_ale_displacement_ref->unknown(0)(i);

        //move the mesh
        M_mesh_mover.apply(M_movingMesh, *M_displacement );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;

        this->log(prefixvm(this->prefix(),"MeshALE"),"init", "restart on moving mesh");

        // rebuild dof point (necessary in updateIdentityMap with extended dof table)
        M_Xhmove->rebuildDofPoints();
        // update M_identity_ale
        this->updateIdentityMap();


#if 1
        *M_meshVelocity = M_bdf_ale_velocity->unknown(0);
#else
        // not good
        *M_meshVelocity = *M_identity_ale;
        M_meshVelocity->scale( M_bdf_ale_identity->polyDerivCoefficient(0) );
        M_meshVelocity->add(-1.0,M_bdf_ale_identity->polyDeriv());
#endif


    }


    for ( auto const& faceWrap : markedfaces(M_movingMesh,this->aleFactory()->flagSet("moving") ) )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto facedof = M_displacementOnMovingBoundary_HO_ref->functionSpace()->dof()->faceLocalDof( face.id() );
        for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
            M_dofsOnMovingBoundary_HO.insert( it->index() );
    }

    M_aleFactory->init();

    this->log(prefixvm(this->prefix(),"MeshALE"),"init", "finish");
}

template< class Convex >
void
MeshALE<Convex>::initTimeStep()
{
    // create time schemes
    double timeInitial = doption(_name="ts.time-initial");
    double timestep = doption(_name="ts.time-step");
    double timeFinal = doption(_name="ts.time-final");
    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(this->rootRepository())/fs::path( prefixvm(this->prefix(),prefixvm(this->subPrefix(),"ts")) );

    M_bdf_ale_identity = bdf( _space=M_Xhmove,
                              _name="identity"+suffixName,
                              _prefix=this->prefix(),
                              _initial_time=timeInitial,
                              _final_time=timeFinal,
                              _time_step=timestep,
                              _restart=M_isARestart,
                              _restart_path=M_restartPath,
                              _restart_at_last_save=boption(_name="ts.restart.at-last-save"),
                              _save=boption(_name="ts.save"),_freq=ioption(_name="ts.save.freq"),
                              _format=myFileFormat
                              );
    M_bdf_ale_identity->setPathSave( ( saveTsDir/"identity" ).string() );

    M_bdf_ale_velocity = bdf( _space=M_Xhmove,
                              _name="velocity"+suffixName,
                              _prefix=this->prefix(),
                              _initial_time=timeInitial,
                              _final_time=timeFinal,
                              _time_step=timestep,
                              _restart=M_isARestart,
                              _restart_path=M_restartPath,
                              _restart_at_last_save=boption(_name="ts.restart.at-last-save"),
                              _save=boption(_name="ts.save"),_freq=ioption(_name="ts.save.freq"),
                              _format=myFileFormat
                              );
    M_bdf_ale_velocity->setPathSave( ( saveTsDir/"velocity" ).string() );

    M_bdf_ale_displacement_ref = bdf( _space=M_Xhref,
                                      _name="displacement_ref"+suffixName,
                                      _prefix=this->prefix(),
                                      _initial_time=timeInitial,
                                      _final_time=timeFinal,
                                      _time_step=timestep,
                                      _restart=M_isARestart,
                                      _restart_path=M_restartPath,
                                      _restart_at_last_save=boption(_name="ts.restart.at-last-save"),
                                      _save=boption(_name="ts.save"),_freq=ioption(_name="ts.save.freq"),
                                      _format=myFileFormat
                                      );
    M_bdf_ale_displacement_ref->setPathSave( ( saveTsDir/"displacement_ref" ).string() );
}

template< class Convex >
std::shared_ptr<std::ostringstream>
MeshALE<Convex>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n||==============================================||"
           << "\n||----------------Info : MeshALE----------------||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix();
    *_ostr << this->aleFactory()->getInfo()->str()
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::addBoundaryFlags(std::string __type, std::string __marker)
{
    CHECK( this->referenceMesh()->hasMarker(__marker) ) << " marker " << __marker << " is not define in reference mesh\n";
    M_aleFactory->addBoundaryFlags( __type, this->referenceMesh()->markerName(__marker) );
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::mesh_ref_ptrtype
MeshALE<Convex>::referenceMesh()
{
    return M_referenceMesh;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::mesh_ptrtype
MeshALE<Convex>::movingMesh()
{
    return M_movingMesh;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_functionspace_ptrtype
MeshALE<Convex>::functionSpace()
{
    return M_Xhmove;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_functionspace_ref_ptrtype
MeshALE<Convex>::functionSpaceInRef()
{
    return M_Xhref;
}

//------------------------------------------------------------------------------------------------//
template< class Convex >
typename MeshALE<Convex>::ale_map_ptrtype const&
MeshALE<Convex>::aleFactory() const
{
    return M_aleFactory;
}
//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ref_type const&
MeshALE<Convex>::mapInRef() const
{
    return M_aleFactory->map();
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_type
MeshALE<Convex>::map()
{
    auto mapRef= M_aleFactory->map();

    ale_map_element_type mapMove(M_Xhmove);

    for (size_type i=0;i<mapMove.nLocalDof();++i)
        (mapMove)(M_drm->dofRelMap()[i])=(mapRef)(i);

    return mapMove;

}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::identityALE()
{
    return M_identity_ale;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::dispP1ToHO_ref()
{
    return M_dispP1ToHO_ref;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::displacementOnMovingBoundary()
{
    return M_displacementOnMovingBoundary_HO_ref;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::displacementOnMovingBoundaryInRef()
{
    return M_displacementOnMovingBoundary_P1_ref;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::displacement()
{
    return M_displacement;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::displacementInRef()
{
    return M_displacement_ref;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::velocity() { return M_meshVelocity;}

template< class Convex >
typename MeshALE<Convex>::ale_map_element_ptrtype const&
MeshALE<Convex>::velocity() const { return M_meshVelocity; }

//------------------------------------------------------------------------------------------------//

template< class Convex >
typename MeshALE<Convex>::DofRelationshipMap_ptrtype
MeshALE<Convex>::dofRelationShipMap()
{
    return M_drm;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::updateIdentityMap()
{
    M_identity_ale->on(_range=elements( M_identity_ale->mesh() ),
                       _expr=vf::P() );
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::revertReferenceMesh( bool updateMeshMeasures )
{
    if ( !this->isOnReferenceMesh() )
    {

#if 1
        auto temporaryDisp = *M_displacement;
#else
        auto temporaryDisp = M_Xhmove->element();
        for (size_type i=0;i<temporaryDisp.nLocalDof();++i)
            temporaryDisp(M_drm->dofRelMap()[i])=(*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i);
#endif
        temporaryDisp.scale(-1.);
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply(M_movingMesh, temporaryDisp );
        M_mesh_mover.setUpdateMeshMeasures( true );
        M_isOnReferenceMesh = true;
        M_isOnMovingMesh = false;
    }
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::revertMovingMesh( bool updateMeshMeasures )
{
    if ( !this->isOnMovingMesh() )
    {
#if 1
        auto temporaryDisp = *M_displacement;
#else
        auto temporaryDisp = M_Xhmove->element();
        for (size_type i=0;i<temporaryDisp.nLocalDof();++i)
            temporaryDisp(M_drm->dofRelMap()[i])=(*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i);
#endif
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply(M_movingMesh, temporaryDisp );
        M_mesh_mover.setUpdateMeshMeasures( true );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;
    }
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::updateBdf()
{
    M_bdf_ale_identity->shiftRight( *M_identity_ale );
    M_bdf_ale_velocity->shiftRight( *M_meshVelocity );
    M_bdf_ale_displacement_ref->shiftRight( *M_displacement_ref );

    M_bdf_ale_identity->next();
    M_bdf_ale_velocity->next();
    M_bdf_ale_displacement_ref->next();
}
//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::exportResults(double time)
{
    std::string exportRepository = (fs::path(this->rootRepository()) / prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports"))).string();
    if (!M_exporter_ref)
    {
        //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter_ref = exporter( _mesh=this->referenceMesh(),
                                   //ame=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                                   _name=prefixvm(this->prefix(),"exportMeshALE_ref"),
                                   _geo=geoExportType,
                                   _worldcomm=this->worldComm(),
                                   _path=exportRepository );
    }

    M_exporter_ref->step( time )->add( prefixvm(this->prefix(),"ref_displacement"), *M_displacement_ref );
    M_exporter_ref->step( time )->add( prefixvm(this->prefix(),"ref_map"),M_aleFactory->map());
    M_exporter_ref->step( time )->add( prefixvm(this->prefix(),"ref_dispOnMovingBoundary"),*this->displacementOnMovingBoundaryInRef());
    M_exporter_ref->save();

    if ( mesh_type::nOrder == 1 )
    {
        if (!M_exporter)
        {
            //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
            std::string geoExportType="static";//change_coords_only, change, static
            M_exporter = exporter( _mesh=this->movingMesh(),
                                   //ame=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                                   _name=prefixvm(this->prefix(),"exportMeshALE_ho"),
                                   _geo=geoExportType,
                                   _worldcomm=this->worldComm(),
                                   _path=exportRepository );
        }

        //M_exporter->step( time )->setMesh( M_movingMesh );
        M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_displacement"), *M_displacement );
        M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_meshvelocity"), *M_meshVelocity );
        M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_identity"), *M_identity_ale );
        //M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_identitypolyderiv"), M_bdf_ale_identity->polyDeriv() );
        M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_dispOnMovingBoundary"), *M_displacementOnMovingBoundary_HO_ref );
        M_exporter->save();
    }
}

template< class Convex >
void
MeshALE<Convex>::updateImpl()
{
    //---------------------------------------------------------------------------------------------//
    // transform disp from ref_ho -> ref_p1
    for (size_type i=0;i<M_dispP1ToHO_ref->nLocalDof();++i)
        (*M_displacementOnMovingBoundary_P1_ref)(i) = (*M_displacementOnMovingBoundary_HO_ref)(M_drm->dofRelMap()[i]) + (*M_dispP1ToHO_ref)(i);

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    boost::mpi::timer mpiTimer;
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", "start generateMap");
    //---------------------------------------------------------------------------------------------//
    //generate the ale map
    M_aleFactory->generateMap( *M_displacementOnMovingBoundary_P1_ref, M_bdf_ale_displacement_ref->unknown(0) );
    //---------------------------------------------------------------------------------------------//
    std::ostringstream ostrTime; ostrTime << "boost::mpi::timer :" << mpiTimer.elapsed() << "s";
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", "finish generateMap in "+ostrTime.str());
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//

    //---------------------------------------------------------------------------------------------//
    //usefull for implicit fsi coupling (reset when we are in ptfixe cycle)
    auto resetDisp = *M_displacement;
    resetDisp.scale(-1.);
    M_mesh_mover.apply(M_movingMesh, resetDisp );

    //---------------------------------------------------------------------------------------------//
    // get displacement on ref mesh
    *M_displacement_ref = M_aleFactory->displacement();

    //---------------------------------------------------------------------------------------------//
    // transfert to the move mesh
    for (size_type i=0;i<M_displacement->nLocalDof();++i)
        (*M_displacement)(M_drm->dofRelMap()[i])=(*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i);/////
#if 0
    //---------------------------------------------------------------------------------------------//
    // update to real move
    auto const& dispRefPrevious = M_bdf_ale_displacement_ref->unknown(0);
    for (size_type i=0;i<M_displacement->nLocalDof();++i)
        (*M_displacement)(M_drm->dofRelMap()[i]) -= dispRefPrevious(i) - (*M_dispP1ToHO_ref)(i);//////
#endif
    //---------------------------------------------------------------------------------------------//
    //move the mesh
    M_mesh_mover.apply(M_movingMesh, *M_displacement );

    M_isOnReferenceMesh = false;
    M_isOnMovingMesh = true;

    // rebuild dof point (necessary in updateIdentityMap with extended dof table)
    M_Xhmove->rebuildDofPoints();

    // up identity
    this->updateIdentityMap();

    // compute mesh velocity
    *M_meshVelocity = *M_identity_ale;
    M_meshVelocity->scale( M_bdf_ale_identity->polyDerivCoefficient(0) );
    M_meshVelocity->add(-1.0,M_bdf_ale_identity->polyDeriv());

}



} // namespace FeelModels
} // namespace Feel
