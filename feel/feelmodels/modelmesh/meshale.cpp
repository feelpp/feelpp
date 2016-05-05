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

template <class Convex>
MeshALE<Convex>::MeshALE( mesh_ptrtype mesh_moving,
                          std::string const& prefix, WorldComm const& worldcomm,
                          bool moveGhostEltFromExtendedStencil,
                          std::string const& rootRepository )
    : super_type( prefixvm( prefix, "alemesh" ), worldcomm, "", rootRepository ),
      M_referenceMesh( mesh_moving->createP1mesh() ),
      M_movingMesh( mesh_moving ),
      M_isOnReferenceMesh( true ), M_isOnMovingMesh( true ),
      M_aleFactory( ale_map_type::build( M_referenceMesh, this->prefix(), this->worldComm(), moveGhostEltFromExtendedStencil ) ),
      M_Xhref( M_aleFactory->functionSpace() ),
      M_Xhmove( ale_map_functionspace_type::New( _mesh = M_movingMesh,
                                                 _worldscomm = std::vector<WorldComm>( 1, worldcomm ),
                                                 _extended_doftable = std::vector<bool>( 1, moveGhostEltFromExtendedStencil ) ) ),
      M_identity_ale( new ale_map_element_type( M_Xhmove ) ),
      M_dispP1ToHO_ref( new ale_map_element_ref_type( M_Xhref ) ),
      M_displacementOnMovingBoundary_HO_ref( new ale_map_element_type( M_Xhmove ) ),
      M_displacementOnMovingBoundary_P1_ref( new ale_map_element_ref_type( M_Xhref ) ),
      M_displacement( new ale_map_element_type( M_Xhmove ) ),
      M_displacement_ref( new ale_map_element_ref_type( M_Xhref ) ),
      M_meshVelocity( new ale_map_element_type( M_Xhmove ) ),
      M_drm( new DofRelationshipMap_type( M_Xhref, M_Xhmove ) ),
      M_cpt_export( 0 ),
      M_isARestart( boption( _name = "ts.restart" ) ),
      M_restartPath( soption( _name = "ts.restart.path" ) )
//M_doExport(option(_name="export",_prefix=prefixvm(prefix,"alemesh")).template as<bool>())
{

    this->log( prefixvm( this->prefix(), "MeshALE" ), "constructor", "start" );

    //std::cout << "MESHALE START " << mesh_moving->worldComm().glogodRank() << std::endl;
    EntityProcessType entityProcess = ( moveGhostEltFromExtendedStencil ) ? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    //*M_dispP1ToHO_ref = vf::project( _space=M_Xhref, _range=elements(M_referenceMesh,entityProcess),
    //                                 _expr=vf::P(),
    //                                 _geomap=GeomapStrategyType::GEOMAP_OPT/*HO*/ );
    M_dispP1ToHO_ref->on( _range = elements( M_referenceMesh, entityProcess ),
                          _expr = vf::P() /*,_geomap=GeomapStrategyType::GEOMAP_OPT*/ /*HO*/ );

    // update M_identity_ale
    this->updateIdentityMap();

    // compute dist between P1(ref) to Ho mesh
    for ( size_type i = 0; i < M_dispP1ToHO_ref->nLocalDof(); ++i )
        ( *M_dispP1ToHO_ref )( i ) = ( *M_identity_ale )( M_drm->dofRelMap()[i] ) - ( *M_dispP1ToHO_ref )( i );

    if ( mesh_type::nOrder == mesh_ref_type::nOrder )
        for ( size_type k = 0; k < M_dispP1ToHO_ref->functionSpace()->nLocalDof(); ++k )
            CHECK( math::abs( M_dispP1ToHO_ref->operator()( k ) ) < 1e-14 ) << "error must be null : " << M_dispP1ToHO_ref->operator()( k );

    // create time schemes

    double timestep = doption( _name = "ts.time-step" );
    std::string myFileFormat = soption( _name = "ts.file-format" ); // without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = ( boost::format( "_rank%1%_%2%" ) % this->worldComm().rank() % this->worldComm().size() ).str();

    M_bdf_ale_identity = bdf( _vm = Environment::vm(), _space = M_Xhmove,
                              _name = prefixvm( prefix, "meshale_identity" + suffixName ),
                              _prefix = this->prefix(),
                              _initial_time = doption( _name = "ts.time-initial" ),
                              _final_time = doption( _name = "ts.time-final" ),
                              _time_step = timestep,
                              _restart = M_isARestart,
                              _restart_path = M_restartPath,
                              _restart_at_last_save = boption( _name = "ts.restart.at-last-save" ),
                              _save = boption( _name = "ts.save" ), _freq = ioption( _name = "ts.save.freq" ) );
    M_bdf_ale_identity->setfileFormat( myFileFormat );
    M_bdf_ale_identity->setPathSave( ( fs::path( this->rootRepository() ) /
                                       fs::path( prefixvm( this->prefix(), ( boost::format( "alemesh.bdf_o_%1%_dt_%2%" ) % timestep % M_bdf_ale_identity->bdfOrder() ).str() ) ) )
                                         .string() );

    M_bdf_ale_velocity = bdf( _vm = Environment::vm(), _space = M_Xhmove,
                              _name = prefixvm( prefix, "meshale_velocity" + suffixName ),
                              _prefix = this->prefix(),
                              _initial_time = doption( _name = "bdf.time-initial" ),
                              _final_time = doption( _name = "ts.time-final" ),
                              _time_step = timestep,
                              _restart = M_isARestart,
                              _restart_path = M_restartPath,
                              _restart_at_last_save = boption( _name = "ts.restart.at-last-save" ),
                              _save = boption( _name = "ts.save" ), _freq = ioption( _name = "ts.save.freq" ) );
    M_bdf_ale_velocity->setfileFormat( myFileFormat );
    M_bdf_ale_velocity->setPathSave( ( fs::path( this->rootRepository() ) /
                                       fs::path( prefixvm( this->prefix(), ( boost::format( "alemesh.bdf_o_%1%_dt_%2%" ) % timestep % M_bdf_ale_velocity->bdfOrder() ).str() ) ) )
                                         .string() );

    M_bdf_ale_displacement_ref = bdf( _vm = Environment::vm(), _space = M_Xhref,
                                      _name = prefixvm( prefix, "meshale_displacement_ref" + suffixName ),
                                      _prefix = this->prefix(),
                                      _initial_time = doption( _name = "ts.time-initial" ),
                                      _final_time = doption( _name = "ts.time-final" ),
                                      _time_step = timestep,
                                      _restart = M_isARestart,
                                      _restart_path = M_restartPath,
                                      _restart_at_last_save = boption( _name = "ts.restart.at-last-save" ),
                                      _save = boption( _name = "ts.save" ), _freq = ioption( _name = "ts.save.freq" ) );
    M_bdf_ale_displacement_ref->setfileFormat( myFileFormat );
    M_bdf_ale_displacement_ref->setPathSave( ( fs::path( this->rootRepository() ) /
                                               fs::path( prefixvm( this->prefix(), ( boost::format( "alemesh.bdf_o_%1%_dt_%2%" ) % timestep % M_bdf_ale_displacement_ref->bdfOrder() ).str() ) ) )
                                                 .string() );

    this->log( prefixvm( this->prefix(), "MeshALE" ), "constructor", "finish" );
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::init()
{
    this->log( prefixvm( this->prefix(), "MeshALE" ), "init", "start" );

    if ( !M_isARestart )
    {
        *M_displacement_ref = *M_dispP1ToHO_ref;

        M_bdf_ale_identity->start( *M_identity_ale );
        M_bdf_ale_velocity->start( *M_meshVelocity );
        M_bdf_ale_displacement_ref->start( *M_displacement_ref );
    }
    else
    {
        M_bdf_ale_identity->restart();
        M_bdf_ale_velocity->restart();
        M_bdf_ale_displacement_ref->restart();

        // get displacement from t0
        *M_displacement_ref = M_bdf_ale_displacement_ref->unknown( 0 );

        // transfert displacement on the mobile mesh
        for ( size_type i = 0; i < M_displacement->nLocalDof(); ++i )
            ( *M_displacement )( M_drm->dofRelMap()[i] ) = ( *M_displacement_ref )( i ) - ( *M_dispP1ToHO_ref )( i );

        //for (uint i=0;i<M_displacement->nDof();++i)
        //    (*M_displacement)(M_drm->dofRelMap()[i]) -= M_bdf_ale_displacement_ref->unknown(0)(i);

        //move the mesh
        M_mesh_mover.apply( M_movingMesh, *M_displacement );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;

        this->log( prefixvm( this->prefix(), "MeshALE" ), "init", "restart on moving mesh" );

#if 0
        //important!!!
        M_displacement->zero();
#endif

#if 0
        *M_identity_ale = vf::project( _space=M_identity_ale->functionSpace(), _range=elements( M_identity_ale->mesh() ),
                                       _expr=vf::P(),
                                       _geomap=GeomapStrategyType::GEOMAP_OPT/*HO*/);
#else

        // rebuild dof point (necessary in updateIdentityMap with extended dof table)
        M_Xhmove->rebuildDofPoints();
        // update M_identity_ale
        this->updateIdentityMap();
#endif

#if 1
        *M_meshVelocity = M_bdf_ale_velocity->unknown( 0 );
#else
        // not good
        *M_meshVelocity = *M_identity_ale;
        M_meshVelocity->scale( M_bdf_ale_identity->polyDerivCoefficient( 0 ) );
        M_meshVelocity->add( -1.0, M_bdf_ale_identity->polyDeriv() );
#endif
    }

    this->log( prefixvm( this->prefix(), "MeshALE" ), "init", "finish" );
}

template <class Convex>
boost::shared_ptr<std::ostringstream>
MeshALE<Convex>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

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

template <class Convex>
void MeshALE<Convex>::addBoundaryFlags( std::string __type, std::string __marker )
{
    CHECK( this->referenceMesh()->hasMarker( __marker ) ) << " marker " << __marker << " is not define in reference mesh\n";
    M_aleFactory->addBoundaryFlags( __type, this->referenceMesh()->markerName( __marker ) );
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::mesh_ref_ptrtype
MeshALE<Convex>::referenceMesh()
{
    return M_referenceMesh;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::mesh_ptrtype
MeshALE<Convex>::movingMesh()
{
    return M_movingMesh;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_functionspace_ptrtype
MeshALE<Convex>::functionSpace()
{
    return M_Xhmove;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_functionspace_ref_ptrtype
MeshALE<Convex>::functionSpaceInRef()
{
    return M_Xhref;
}

//------------------------------------------------------------------------------------------------//
template <class Convex>
typename MeshALE<Convex>::ale_map_ptrtype const&
MeshALE<Convex>::aleFactory() const
{
    return M_aleFactory;
}
//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ref_type const&
MeshALE<Convex>::mapInRef() const
{
    return M_aleFactory->map();
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_type
MeshALE<Convex>::map()
{
    auto mapRef = M_aleFactory->map();

    ale_map_element_type mapMove( M_Xhmove );

    for ( size_type i = 0; i < mapMove.nLocalDof(); ++i )
        ( mapMove )( M_drm->dofRelMap()[i] ) = ( mapRef )( i );

    return mapMove;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::identityALE()
{
    return M_identity_ale;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::dispP1ToHO_ref()
{
    return M_dispP1ToHO_ref;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::displacementOnMovingBoundary()
{
    return M_displacementOnMovingBoundary_HO_ref;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::displacementOnMovingBoundaryInRef()
{
    return M_displacementOnMovingBoundary_P1_ref;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::displacement()
{
    return M_displacement;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ref_ptrtype
MeshALE<Convex>::displacementInRef()
{
    return M_displacement_ref;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ptrtype
MeshALE<Convex>::velocity()
{
    return M_meshVelocity;
}

template <class Convex>
typename MeshALE<Convex>::ale_map_element_ptrtype const&
MeshALE<Convex>::velocity() const
{
    return M_meshVelocity;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
typename MeshALE<Convex>::DofRelationshipMap_ptrtype
MeshALE<Convex>::dofRelationShipMap()
{
    return M_drm;
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::updateIdentityMap()
{
    bool moveGhostEltFromExtendedStencil = M_Xhmove->dof()->buildDofTableMPIExtended() && M_movingMesh->worldComm().localSize() > 1;
    EntityProcessType entityProcess = ( moveGhostEltFromExtendedStencil ) ? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    //*M_identity_ale = vf::project( _space=M_identity_ale->functionSpace(),
    //                               _range=elements( M_identity_ale->mesh(),entityProcess ),
    //                               _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_OPT/*HO*/ );
    M_identity_ale->on( _range = elements( M_identity_ale->mesh(), entityProcess ),
                        _expr = vf::P(),
                        _geomap = ( mesh_type::nOrder > 1 ) ? GeomapStrategyType::GEOMAP_OPT : GeomapStrategyType::GEOMAP_HO );

#if 0
    if ( M_Xhmove->dof()->buildDofTableMPIExtended() && M_movingMesh->worldComm().localSize()>1 )
    {
        auto rangeMove = interprocessfaces(M_movingMesh);
        for ( auto face_it = rangeMove.get<1>(), face_en = rangeMove.get<2>() ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( uint16_type l =0; l < ale_map_functionspace_type::dof_type::fe_type::nLocalDof; ++l )
            {
                int ncdof  = ale_map_functionspace_type::dof_type::is_product?ale_map_functionspace_type::dof_type::nComponents:1;
                for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                {
                    const size_type thedof = M_Xhmove->dof()->localToGlobal(eltOffProc,l,c1).index();
                    const double dofptCoord = M_Xhmove->dof()->dofPoint( thedof ).template get<0>()[c1];
                    M_identity_ale->set(thedof, dofptCoord);
                }
            }
        }
    }
#endif
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::revertReferenceMesh( bool updateMeshMeasures )
{
    if ( !this->isOnReferenceMesh() )
    {

#if 1
        auto temporaryDisp = *M_displacement;
#else
        auto temporaryDisp = M_Xhmove->element();
        for ( size_type i = 0; i < temporaryDisp.nLocalDof(); ++i )
            temporaryDisp( M_drm->dofRelMap()[i] ) = ( *M_displacement_ref )( i ) - ( *M_dispP1ToHO_ref )( i );
#endif
        temporaryDisp.scale( -1. );
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply( M_movingMesh, temporaryDisp );
        M_mesh_mover.setUpdateMeshMeasures( true );
        M_isOnReferenceMesh = true;
        M_isOnMovingMesh = false;
    }
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::revertMovingMesh( bool updateMeshMeasures )
{
    if ( !this->isOnMovingMesh() )
    {
#if 1
        auto temporaryDisp = *M_displacement;
#else
        auto temporaryDisp = M_Xhmove->element();
        for ( size_type i = 0; i < temporaryDisp.nLocalDof(); ++i )
            temporaryDisp( M_drm->dofRelMap()[i] ) = ( *M_displacement_ref )( i ) - ( *M_dispP1ToHO_ref )( i );
#endif
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply( M_movingMesh, temporaryDisp );
        M_mesh_mover.setUpdateMeshMeasures( true );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;
    }
}

//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::updateBdf()
{
    M_bdf_ale_identity->shiftRight( *M_identity_ale );
    M_bdf_ale_velocity->shiftRight( *M_meshVelocity );
    M_bdf_ale_displacement_ref->shiftRight( *M_displacement_ref );

    M_bdf_ale_identity->next();
    M_bdf_ale_velocity->next();
    M_bdf_ale_displacement_ref->next();
#if 0
    //Attention !!!!!
    M_displacement->zero();
#endif
}
//------------------------------------------------------------------------------------------------//

template <class Convex>
void MeshALE<Convex>::exportResults( double time )
{
    if ( !M_exporter_ref )
    {
        //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
        std::string geoExportType = "static"; //change_coords_only, change, static
        M_exporter_ref = exporter( _mesh = this->referenceMesh(),
                                   //ame=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                                   _name = prefixvm( this->prefix(), "exportMeshALE_ref" ),
                                   _geo = geoExportType,
                                   _worldcomm = this->worldComm() );
    }

    M_exporter_ref->step( time )->add( prefixvm( this->prefix(), "ref_displacement" ), *M_displacement_ref );
    M_exporter_ref->step( time )->add( prefixvm( this->prefix(), "ref_map" ), M_aleFactory->map() );
    M_exporter_ref->step( time )->add( prefixvm( this->prefix(), "ref_dispOnMovingBoundary" ), *this->displacementOnMovingBoundaryInRef() );
    M_exporter_ref->save();

    if ( mesh_type::nOrder == 1 )
    {
        if ( !M_exporter )
        {
            //auto const geoExportType = ExporterGeometry::EXPORTER_GEOMETRY_STATIC;
            std::string geoExportType = "static"; //change_coords_only, change, static
            M_exporter = exporter( _mesh = this->movingMesh(),
                                   //ame=prefixvm(this->prefix(), prefixvm(this->subPrefix(),"Export")),
                                   _name = prefixvm( this->prefix(), "exportMeshALE_ho" ),
                                   _geo = geoExportType,
                                   _worldcomm = this->worldComm() );
        }

        //M_exporter->step( time )->setMesh( M_movingMesh );
        M_exporter->step( time )->add( prefixvm( this->prefix(), "moving_displacement" ), *M_displacement );
        M_exporter->step( time )->add( prefixvm( this->prefix(), "moving_meshvelocity" ), *M_meshVelocity );
        M_exporter->step( time )->add( prefixvm( this->prefix(), "moving_identity" ), *M_identity_ale );
        //M_exporter->step( time )->add( prefixvm(this->prefix(),"moving_identitypolyderiv"), M_bdf_ale_identity->polyDeriv() );
        M_exporter->step( time )->add( prefixvm( this->prefix(), "moving_dispOnMovingBoundary" ), *M_displacementOnMovingBoundary_HO_ref );
        M_exporter->save();
    }
}

template <class Convex>
void MeshALE<Convex>::updateImpl( Vector<double> const& dispInput )
{
    //M_displacementOnMovingBoundary_HO_ref->zero();
    *M_displacementOnMovingBoundary_HO_ref = dispInput;

    //---------------------------------------------------------------------------------------------//
    // transform disp from ref_ho -> ref_p1
    for ( size_type i = 0; i < M_dispP1ToHO_ref->nLocalDof(); ++i )
        ( *M_displacementOnMovingBoundary_P1_ref )( i ) = ( *M_displacementOnMovingBoundary_HO_ref )( M_drm->dofRelMap()[i] ) + ( *M_dispP1ToHO_ref )( i );

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    boost::mpi::timer mpiTimer;
    this->log( prefixvm( this->prefix(), "MeshALE" ), "updateImpl", "start generateMap" );
    //---------------------------------------------------------------------------------------------//
    //generate the ale map
    M_aleFactory->generateMap( *M_displacementOnMovingBoundary_P1_ref, M_bdf_ale_displacement_ref->unknown( 0 ) );
    //---------------------------------------------------------------------------------------------//
    std::ostringstream ostrTime;
    ostrTime << "boost::mpi::timer :" << mpiTimer.elapsed() << "s";
    this->log( prefixvm( this->prefix(), "MeshALE" ), "updateImpl", "finish generateMap in " + ostrTime.str() );
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//

    //---------------------------------------------------------------------------------------------//
    //usefull for implicit fsi coupling (reset when we are in ptfixe cycle)
    auto resetDisp = *M_displacement;
    resetDisp.scale( -1. );
    M_mesh_mover.apply( M_movingMesh, resetDisp );

    //---------------------------------------------------------------------------------------------//
    // get displacement on ref mesh
    *M_displacement_ref = M_aleFactory->displacement();

    //---------------------------------------------------------------------------------------------//
    // transfert to the move mesh
    for ( size_type i = 0; i < M_displacement->nLocalDof(); ++i )
        ( *M_displacement )( M_drm->dofRelMap()[i] ) = ( *M_displacement_ref )( i ) - ( *M_dispP1ToHO_ref )( i ); /////
#if 0
    //---------------------------------------------------------------------------------------------//
    // update to real move
    auto const& dispRefPrevious = M_bdf_ale_displacement_ref->unknown(0);
    for (size_type i=0;i<M_displacement->nLocalDof();++i)
        (*M_displacement)(M_drm->dofRelMap()[i]) -= dispRefPrevious(i) - (*M_dispP1ToHO_ref)(i);//////
#endif
    //---------------------------------------------------------------------------------------------//
    //move the mesh
    M_mesh_mover.apply( M_movingMesh, *M_displacement );

    M_isOnReferenceMesh = false;
    M_isOnMovingMesh = true;

    // rebuild dof point (necessary in updateIdentityMap with extended dof table)
    M_Xhmove->rebuildDofPoints();

    // up identity
    this->updateIdentityMap();

    // compute mesh velocity
    *M_meshVelocity = *M_identity_ale;
    M_meshVelocity->scale( M_bdf_ale_identity->polyDerivCoefficient( 0 ) );
    M_meshVelocity->add( -1.0, M_bdf_ale_identity->polyDeriv() );
}

} // namespace FeelModels
} // namespace Feel
