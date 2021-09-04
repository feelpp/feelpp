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
    M_Xhmove( ale_map_functionspace_type::New(_mesh=M_movingMesh ) ),
    M_identity_ale(new ale_map_element_type( M_Xhmove) ),
    M_displacementOnMovingBoundary_HO_ref(new ale_map_element_type( M_Xhmove) ),
    M_displacement(new ale_map_element_type( M_Xhmove) ),
    M_meshVelocity(new ale_map_element_type( M_Xhmove) ),
    M_fieldTmp( new  ale_map_element_type( M_Xhmove) ),
    M_isARestart(boption(_name="ts.restart")),
    M_restartPath(soption(_name="ts.restart.path"))
    //M_doExport(option(_name="export",_prefix=prefixvm(prefix,"alemesh")).template as<bool>())
{
    this->log(prefixvm(this->prefix(),"MeshALE"),"constructor", "start");

    // update M_identity_ale
    this->updateIdentityMap();
#if 0
    // compute dist between P1(ref) to Ho mesh
    if ( mesh_type::nOrder != mesh_ref_type::nOrder )
    {
        M_dispP1ToHO_ref->on(_range=elements(M_referenceMesh),
                             _expr=vf::P() );
        for (size_type i=0;i<M_dispP1ToHO_ref->nLocalDof();++i)
            (*M_dispP1ToHO_ref)(i) = (*M_identity_ale)(M_drm->dofRelMap()[i]) - (*M_dispP1ToHO_ref)(i);
    }
#endif


    this->log(prefixvm(this->prefix(),"MeshALE"),"constructor", "finish");
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::init()
{
    this->log(prefixvm(this->prefix(),"MeshALE"),"init", "start");

    if ( M_computationalDomains.empty() )
        this->setWholeMeshAsComputationalDomain( "default");

    for ( auto & [name,aleRefDomain] : M_computationalDomains )
        aleRefDomain.init( M_isARestart );

    this->initTimeStep();

    if (!M_isARestart)
    {
        //this->updateIdentityMap();
        M_bdf_ale_identity->start(*M_identity_ale);
        M_bdf_ale_velocity->start(*M_meshVelocity);
    }
    else
    {
        M_bdf_ale_identity->restart();
        M_bdf_ale_velocity->restart();

        // transfert displacement on the mobile mesh
        for ( auto & [name,aleRefDomain] : M_computationalDomains )
            aleRefDomain.updateDisplacement( *M_displacement );

#if 0
        //move the mesh
        M_mesh_mover.apply( M_movingMesh, *M_displacement );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;

        this->log(prefixvm(this->prefix(),"MeshALE"),"init", "restart on moving mesh");

        // update M_identity_ale
        this->updateIdentityMap();
#endif
        *M_identity_ale = M_bdf_ale_identity->unknown(0);
        *M_meshVelocity = M_bdf_ale_velocity->unknown(0);
        M_isOnReferenceMesh = true;
        M_isOnMovingMesh = false;
    }

#if 0
    for ( auto const& [name,aleRefDomain] : M_computationalDomains ) // TODO VINCENT
    {
        auto therange = markedfaces(M_movingMesh,aleRefDomain.aleFactory()->flagSet("moving") );
        //M_dofsMultiProcessOnMovingBoundary_HO = M_displacementOnMovingBoundary_HO_ref->functionSpace()->dofs( therange, ComponentType::NO_COMPONENT, true );
        auto dofMultiProcessOnRange = M_displacementOnMovingBoundary_HO_ref->functionSpace()->dofs( therange, ComponentType::NO_COMPONENT, true );
        M_dofsMultiProcessOnMovingBoundary_HO.insert( dofMultiProcessOnRange.begin(), dofMultiProcessOnRange.end() );
    }
#endif

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
    for ( auto const& [name,aleRefDomain] : M_computationalDomains )
        *_ostr << aleRefDomain.aleFactory()->getInfo()->str()
               << "\n||==============================================||"
               << "\n";
    return _ostr;
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::setWholeMeshAsComputationalDomain( std::string const& name )
{
    CHECK( M_computationalDomains.empty() ) << "can't define on whole mesh";
    ComputationalDomain hola( this );
    M_computationalDomains.emplace( std::make_pair( name, std::move( hola ) ) );
}

template< class Convex >
void
MeshALE<Convex>::setComputationalDomain( std::string const& name, range_elements_type const& rangeElt )
{
    ComputationalDomain hola( this,rangeElt );
    M_computationalDomains.emplace( std::make_pair( name, std::move( hola ) ) );
}

template< class Convex >
void
MeshALE<Convex>::addBoundaryFlags(std::string const& bctype, std::string const& marker)
{
    //CHECK( this->referenceMesh()->hasMarker(marker) ) << " marker " << marker << " is not define in reference mesh\n";
    for ( auto & [name,aleRefDomain] : M_computationalDomains )
        aleRefDomain.addBoundaryFlags( bctype, marker );
}

//------------------------------------------------------------------------------------------------//
#if 0
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
#endif
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
        *M_fieldTmp =  *M_displacement;
        M_fieldTmp->scale(-1.);
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply(M_movingMesh, *M_fieldTmp );
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
        M_mesh_mover.setUpdateMeshMeasures( updateMeshMeasures );
        M_mesh_mover.apply(M_movingMesh, *M_displacement );
        M_mesh_mover.setUpdateMeshMeasures( true );
        M_isOnReferenceMesh = false;
        M_isOnMovingMesh = true;
    }
}

//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::updateTimeStep()
{
    M_bdf_ale_identity->next( *M_identity_ale );
    M_bdf_ale_velocity->next( *M_meshVelocity );
    for ( auto & [name,aleRefDomain] : M_computationalDomains )
        aleRefDomain.updateTimeStep();
}
//------------------------------------------------------------------------------------------------//

template< class Convex >
void
MeshALE<Convex>::exportResults(double time)
{
    std::string exportRepository = (fs::path(this->rootRepository()) / prefixvm(this->prefix(), prefixvm(this->subPrefix(),"exports"))).string();
#if 0
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
#endif
    if constexpr ( mesh_type::nOrder <= 2 )
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
MeshALE<Convex>::updateDisplacementImposedForUse()
{
    sync( *M_displacementOnMovingBoundary_HO_ref, "=", M_dofsMultiProcessOnMovingBoundary_HO );
    M_dofsMultiProcessOnMovingBoundary_HO.clear();
}

template< class Convex >
void
MeshALE<Convex>::updateImpl()
{
    this->updateDisplacementImposedForUse();

    //---------------------------------------------------------------------------------------------//
    // transform disp from ref_ho -> ref_p1
#if 0
    for (size_type i=0;i<M_displacementOnMovingBoundary_P1_ref->nLocalDof();++i)
    {
        if constexpr ( mesh_type::nOrder != mesh_ref_type::nOrder )
            (*M_displacementOnMovingBoundary_P1_ref)(i) = (*M_displacementOnMovingBoundary_HO_ref)(M_drm->dofRelMap()[i]) + (*M_dispP1ToHO_ref)(i);
        else
            (*M_displacementOnMovingBoundary_P1_ref)(i) = (*M_displacementOnMovingBoundary_HO_ref)(M_drm->dofRelMap()[i]);
    }
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    boost::mpi::timer mpiTimer;
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", "start generateMap");
    //---------------------------------------------------------------------------------------------//
    //generate the ale map
    M_aleFactory->generateMap( *M_displacementOnMovingBoundary_P1_ref, M_bdf_ale_displacement_ref->unknown(0) );
    //---------------------------------------------------------------------------------------------//
    double tElapsed =  mpiTimer.elapsed();
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", (boost::format("finish generateMap in %1% s")%tElapsed).str() );
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //usefull for implicit fsi coupling (reset when we are in ptfixe cycle)
    //auto resetDisp = *M_displacement;
    //resetDisp.scale(-1.);
    //M_mesh_mover.apply(M_movingMesh, resetDisp );


    //---------------------------------------------------------------------------------------------//
    // get displacement on ref mesh
    *M_displacement_ref = M_aleFactory->displacement();
#else
    for ( auto & [name,aleRefDomain] : M_computationalDomains )
    {
        aleRefDomain.generateMap(*M_displacementOnMovingBoundary_HO_ref);
    }
#endif

    // save previous disp
    if ( M_isOnMovingMesh )
        *M_fieldTmp = *M_displacement;

    // copy displacement imposed (maybe some part are only imposed,no required computation)
    *M_displacement = *M_displacementOnMovingBoundary_HO_ref;

    //---------------------------------------------------------------------------------------------//
    // get displacement on mesh and move it
    if ( M_isOnReferenceMesh )
    {
#if 0
        for (size_type i=0;i<M_displacement->nLocalDof();++i)
        {
            size_type j = M_drm->dofRelMap()[i];
            if constexpr ( mesh_type::nOrder != mesh_ref_type::nOrder )
                             M_displacement->set( j, (*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i) );
            else
                M_displacement->set( j, (*M_displacement_ref)(i) );
        }
#else
        for ( auto & [name,aleRefDomain] : M_computationalDomains )
            aleRefDomain.updateDisplacement( *M_displacement );
#endif
        //move the mesh
        M_mesh_mover.apply( M_movingMesh, *M_displacement );
    }
    else // M_isOnMovingMesh
    {
        auto displacementToApply = M_fieldTmp;
        //*displacementToApply = *M_displacement; // save previous disp
#if 0
        for (size_type i=0;i<M_displacement->nLocalDof();++i)
        {
            size_type j = M_drm->dofRelMap()[i];
            if constexpr ( mesh_type::nOrder != mesh_ref_type::nOrder )
                             M_displacement->set( j, (*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i) );
            else
                M_displacement->set( j, (*M_displacement_ref)(i) );

            displacementToApply->set( j, (*M_displacement)(j) - (*displacementToApply)(j) );
        }
#else
        for ( auto & [name,aleRefDomain] : M_computationalDomains )
            aleRefDomain.updateDisplacement( *M_displacement );
        for (size_type j=0;j<M_displacement->nLocalDof();++j)
            displacementToApply->set( j, (*M_displacement)(j) - (*displacementToApply)(j) );
#endif
        //move the mesh
        M_mesh_mover.apply( M_movingMesh, *displacementToApply );
    }
    M_isOnReferenceMesh = false;
    M_isOnMovingMesh = true;

    //---------------------------------------------------------------------------------------------//

    // rebuild dof point (necessary in updateIdentityMap with extended dof table)
    //M_Xhmove->rebuildDofPoints();

    // up identity
    this->updateIdentityMap();

    // compute mesh velocity
    *M_meshVelocity = *M_identity_ale;
    M_meshVelocity->scale( M_bdf_ale_identity->polyDerivCoefficient(0) );
    M_meshVelocity->add(-1.0,M_bdf_ale_identity->polyDeriv());

}


template< class Convex >
MeshALE<Convex>::ComputationalDomain::ComputationalDomain( self_type const* meshALE )
    :
    M_meshALE( meshALE ),
    M_aleFactory( ale_map_type::build(M_meshALE->referenceMesh(), M_meshALE->prefix(), M_meshALE->worldCommPtr(), M_meshALE->repository()  ) )
{
    this->build();
}

template< class Convex >
MeshALE<Convex>::ComputationalDomain::ComputationalDomain( self_type const* meshALE, range_elements_type const& rangeElt )
    :
    M_meshALE( meshALE )
{
    // transform rangeElt (on moving mesh) into a range of element belongs to reference mesh
    CHECK( M_meshALE->referenceMesh()->isSubMeshFrom( M_meshALE->movingMesh() ) ) << "referenceMesh and movingMesh should be related";
    typename MeshTraits<mesh_ref_type>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<mesh_ref_type>::elements_reference_wrapper_type );
    for ( auto const& eltWrap : rangeElt )
    {
        auto const& elt = unwrap_ref( eltWrap );
        size_type eltRefId = M_meshALE->referenceMesh()->meshToSubMesh( elt.id() );
        if ( eltRefId == invalid_v<size_type> )
            continue;
        auto const& eltRef = M_meshALE->referenceMesh()->element( eltRefId );
        myelts->push_back( boost::cref( eltRef ) );
    }
    myelts->shrink_to_fit();
    range_elements_ref_type rangeEltOnRef = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),myelts->begin(),myelts->end(),myelts );

    // build data
    M_aleFactory = ale_map_type::build(M_meshALE->referenceMesh(), rangeEltOnRef, M_meshALE->prefix(), M_meshALE->worldCommPtr(), M_meshALE->repository()  );
    this->build();
}

template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::build()
{
    M_drm.reset( new DofRelationshipMap_type(/*M_Xhref*/ M_aleFactory->functionSpace(),M_meshALE->functionSpace()/* Xhmove*/ ) );
    auto M_Xhref = M_aleFactory->functionSpace();
    M_dispP1ToHO_ref.reset(new ale_map_element_ref_type( M_Xhref) );
    M_displacementOnMovingBoundary_P1_ref.reset(new ale_map_element_ref_type( M_Xhref) );
    M_displacement_ref.reset(new ale_map_element_ref_type( M_Xhref) );

#if 0 // TODO VINCENT
    // compute dist between P1(ref) to Ho mesh
    if ( mesh_type::nOrder != mesh_ref_type::nOrder )
    {
        M_dispP1ToHO_ref->on(_range=elements(/*M_referenceMesh*/support(M_Xhref)),
                             _expr=vf::P() );
        for (size_type i=0;i<M_dispP1ToHO_ref->nLocalDof();++i)
            (*M_dispP1ToHO_ref)(i) = (*M_identity_ale)(M_drm->dofRelMap()[i]) - (*M_dispP1ToHO_ref)(i);
    }
#endif
}

template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::init( bool M_isARestart )
{
#if 1
    // create time schemes
    double timeInitial = doption(_name="ts.time-initial");
    double timestep = doption(_name="ts.time-step");
    double timeFinal = doption(_name="ts.time-final");
    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%M_meshALE->worldComm().rank()%M_meshALE->worldComm().size() ).str();
    fs::path saveTsDir = fs::path(M_meshALE->rootRepository())/fs::path( prefixvm(M_meshALE->prefix(),prefixvm(M_meshALE->subPrefix(),"ts")) );
    std::string M_restartPath = soption(_name="ts.restart.path");
    auto M_Xhref = M_aleFactory->functionSpace();
    M_bdf_ale_displacement_ref = bdf( _space=M_Xhref,
                                      _name="displacement_ref"+suffixName,
                                      _prefix=M_meshALE->prefix(),
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
#endif

    if (!M_isARestart)
    {
        *M_displacement_ref = *M_dispP1ToHO_ref;
        M_bdf_ale_displacement_ref->start(*M_displacement_ref);
    }
    else
    {
        M_bdf_ale_displacement_ref->restart();

        // get displacement from t0
        *M_displacement_ref = M_bdf_ale_displacement_ref->unknown(0);
    }

    M_aleFactory->init();
}

template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::generateMap( ale_map_element_type const& displacementOnMovingBoundary_HO_ref )
{
    for (size_type i=0;i<M_displacementOnMovingBoundary_P1_ref->nLocalDof();++i)
    {
        if constexpr ( mesh_type::nOrder != mesh_ref_type::nOrder )
            (*M_displacementOnMovingBoundary_P1_ref)(i) = displacementOnMovingBoundary_HO_ref(M_drm->dofRelMap()[i]) + (*M_dispP1ToHO_ref)(i);
        else
            (*M_displacementOnMovingBoundary_P1_ref)(i) = displacementOnMovingBoundary_HO_ref(M_drm->dofRelMap()[i]);
    }

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
#if 0
    boost::mpi::timer mpiTimer;
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", "start generateMap");
#endif
    //---------------------------------------------------------------------------------------------//
    //generate the ale map
    M_aleFactory->generateMap( *M_displacementOnMovingBoundary_P1_ref, M_bdf_ale_displacement_ref->unknown(0) );
    //---------------------------------------------------------------------------------------------//
#if 0
    double tElapsed =  mpiTimer.elapsed();
    this->log(prefixvm(this->prefix(),"MeshALE"),"updateImpl", (boost::format("finish generateMap in %1% s")%tElapsed).str() );
#endif
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//

    // get displacement on ref mesh
    *M_displacement_ref = M_aleFactory->displacement();
}


template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::updateDisplacement( ale_map_element_type & displacement )
{
    for (size_type i=0;i<M_displacement_ref->nLocalDof();++i)
    {
        size_type j = M_drm->dofRelMap()[i];
        if constexpr ( mesh_type::nOrder != mesh_ref_type::nOrder )
            displacement.set( j, (*M_displacement_ref)(i) - (*M_dispP1ToHO_ref)(i) );
        else
            displacement.set( j, (*M_displacement_ref)(i) );
    }
}

template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::updateTimeStep()
{
    M_bdf_ale_displacement_ref->next( *M_displacement_ref );
}

template< class Convex >
void
MeshALE<Convex>::ComputationalDomain::addBoundaryFlags(std::string const& bctype, std::string const& marker)
{
    CHECK( M_meshALE->referenceMesh()->hasMarker(marker) ) << " marker " << marker << " is not define in reference mesh\n";
    M_aleFactory->addBoundaryFlags( bctype, M_meshALE->referenceMesh()->markerName( marker ) );
}

} // namespace FeelModels
} // namespace Feel
