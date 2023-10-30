/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2017-02-09

  Copyright (C) 2017 Feel++ Consortium

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
 \file converteracusimdatabase.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2017-02-09
 */
#ifndef FEELPP_CONVERTERACUSIMDATABASE_HPP
#define FEELPP_CONVERTERACUSIMDATABASE_HPP 1

#if defined( FEELPP_HAS_ACUSIM )

#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <acusim.h>
#include <adb.h>
//#include <h3dpublic_export.h>
#include<sstream>
#include<fstream>
#include <math.h>
//#include <nrutil.h>


namespace Feel
{

template<typename MeshType >
class ConverterAcusimDatabase
{
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Vectorial,Continuous> > > space_vectorial_P1_type;
    typedef std::shared_ptr<space_vectorial_P1_type> space_vectorial_P1_ptrtype;
    typedef typename space_vectorial_P1_type::component_functionspace_type space_scalar_P1_type;
    typedef std::shared_ptr<space_scalar_P1_type> space_scalar_P1_ptrtype;
public :

    ConverterAcusimDatabase()
        :
        M_acusimWorkRepository( "ACUSIM.DIR" ),
        M_buildOnlyScalarSpace( true ),
        M_nNodes( 0 ),
        M_nodeIds( nullptr ),
        M_nOutSteps( 0 ),
        M_outSteps( nullptr ),
        M_outTimes( nullptr )
        {}

    WorldComm const& worldComm() const { return M_feelppDatabase.worldComm(); }

    void setAcusimRepository( std::string const& idir )
        {
            M_acusimRepository = fs::canonical( idir );
        }
    void setAcusimProblemName(std::string const& name )
        {
            M_acusimProblemName = name;
        }
    void setAcusimWorkRepository(std::string const& dir )
        {
            M_acusimWorkRepository = dir;
        }
    void setBuildOnlyScalarSpace( bool b )
        {
            M_buildOnlyScalarSpace = b;
        }

    void setOutputRepository( std::string const& odir )
        {
            M_feelppDatabase.setDbRepository( odir );
        }
    void addField( std::string const& name )
        {
            M_convertedFieldNames.insert( name );
        }

    void init()
        {
            //M_adbHd = adbNew( (char*)"test", (char*)"ACUSIM.DIR", (char*)"/home/u2/chabannes/po_smallmesh/DATA_PO/", 0 );
            M_adbHd = adbNew( (char*)M_acusimProblemName.c_str(), (char*)M_acusimWorkRepository.c_str(), (char*)M_acusimRepository.string().c_str(), 0 );
            CHECK( M_adbHd != NULL ) << "Error opening the data base " << adbGetError();
            // init database loading
            bool success = adbOpenRun(M_adbHd,1);
            CHECK( success ) << "error on opening the acusim database name="<<M_acusimProblemName<<" rep: " << M_acusimRepository;
        }
    void close()
        {
            if ( M_nodeIds )
                delete [] M_nodeIds;
            if ( M_outSteps )
                delete [] M_outSteps;
            if ( M_outTimes )
                delete [] M_outTimes;

            bool success = adbCloseRun( M_adbHd );

            adbFree( M_adbHd );
        }

    void run();

private :
    mesh_ptrtype loadMesh();

    void loadTimeSetInfo();

    void loadFieldInfo();

    template<typename SpaceType>
    void
    generateDofMappingP1( std::shared_ptr<SpaceType> const& space );

    template<typename ElementType>
    void
    applyFieldConversion( int fieldIndex, int timeStepIndex, double * outFieldValues, ElementType & u );

    template<typename ElementType>
    void
    applyFieldConversion( int fieldIndex, int timeStepIndex, double * outFieldValues, ElementType & u, int nComp1, int nComp2, int comp1, int comp2 );

    template<typename SpaceType>
    void
    runSaveFieldInTimeSet( std::shared_ptr<SpaceType> const& space, std::string const& fieldName, int fieldIndex );

    template<typename SpaceType>
    void
    runSaveFieldInTimeSet( std::shared_ptr<SpaceType> const& space, std::string const& fieldName, int fieldIndex, int nComp1, int nComp2 );
private :

    fs::path M_acusimRepository;
    std::string M_acusimProblemName;
    std::string M_acusimWorkRepository;

    bool M_buildOnlyScalarSpace;

    FeelppDatabase<mesh_type> M_feelppDatabase;

    AdbHd M_adbHd;

    std::map<size_type,size_type> M_nodeIdToContainerId;
    std::map<uint16_type,std::vector<size_type>> M_feelDofIdToAcusimDofId;
    int M_nNodes;
    int * M_nodeIds;
    //timeset info
    int M_nOutSteps;
    int * M_outSteps;
    double * M_outTimes;

    std::set<std::string> M_convertedFieldNames;
    std::map<std::string,std::tuple<int,int>> M_fieldNameToFieldInfo;

    space_scalar_P1_ptrtype M_spaceP1Scalar;
    space_vectorial_P1_ptrtype M_spaceP1Vectorial;
};

template <typename MeshType>
void
ConverterAcusimDatabase<MeshType>::run()
{

    this->init();

    bool doMeshConversion = true;
    if ( doMeshConversion )
    {
        auto meshSeq = this->loadMesh();
        M_feelppDatabase.saveMeshParititioned( meshSeq, this->worldComm().size() );
        if ( meshSeq )
        {
            meshSeq->clear();
            meshSeq.reset();
        }
    }

    Environment::worldComm().barrier();
    auto mem3  = Environment::logMemoryUsage("memory usage after update for use");
    std::cout << "[run] resident memory --3-- : " << mem3.memory_usage/1.e9  << "GBytes\n";
    Environment::worldComm().barrier();

    this->loadTimeSetInfo();
    this->loadFieldInfo();

    bool hasScalarSpace = false, hasVectorialSpace = false, hasTensor2Space = false;
    for ( auto const& fieldInfo : M_fieldNameToFieldInfo )
    {
        if ( M_convertedFieldNames.find( fieldInfo.first ) == M_convertedFieldNames.end() )
            continue;
        if ( std::get<1>( fieldInfo.second ) == 1 )
            hasScalarSpace = true;
        else if ( std::get<1>( fieldInfo.second ) == mesh_type::nRealDim )
            hasVectorialSpace = true;
        else if ( std::get<1>( fieldInfo.second ) == (mesh_type::nRealDim*mesh_type::nRealDim) )
            hasTensor2Space = true;

    }
    if ( this->worldComm().isMasterRank() )
        std::cout << "hasScalarSpace : " << hasScalarSpace << "\n"
                  << "hasVectorialSpace : " << hasVectorialSpace << "\n";

    if ( FLAGS_v >= 1 )
    {
        auto mem  = Environment::logMemoryUsage("memory usage after update for use");
        std::cout << "[run] resident memory before loadMesh // : " << mem.memory_usage/1.e9  << "GBytes\n";
        //Environment::worldComm().barrier();
    }

    auto mesh = M_feelppDatabase.loadMesh( /*MESH_UPDATE_FACES_MINIMAL|*/MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED );

    if ( FLAGS_v >= 1 )
    {
        auto mem  = Environment::logMemoryUsage("memory usage after update for use");
        std::cout << "[run] resident memory after loadMesh // : " << mem.memory_usage/1.e9  << "GBytes\n";
        //Environment::worldComm().barrier();
    }

    if ( M_buildOnlyScalarSpace )
    {
        M_spaceP1Scalar = space_scalar_P1_type::New( mesh );
        this->generateDofMappingP1( M_spaceP1Scalar );
    }
    else
    {
        if ( hasVectorialSpace )
        {
            M_spaceP1Vectorial = space_vectorial_P1_type::New( mesh );
            this->generateDofMappingP1( M_spaceP1Vectorial );
            if ( hasScalarSpace )
            {
                M_spaceP1Scalar = M_spaceP1Vectorial->compSpace();
                this->generateDofMappingP1( M_spaceP1Scalar );
            }
        }
        else if ( hasScalarSpace )
        {
            M_spaceP1Scalar = space_scalar_P1_type::New( mesh );
            this->generateDofMappingP1( M_spaceP1Scalar );
        }
    }




    for ( std::string const& fieldName : M_convertedFieldNames )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "start conversion of field : " << fieldName << "\n";

        auto itFindName = M_fieldNameToFieldInfo.find( fieldName );
        CHECK( itFindName != M_fieldNameToFieldInfo.end() ) << "no field " << fieldName << " found in database";
        int fieldIndex = std::get<0>( itFindName->second );
        int outVarDim = std::get<1>( itFindName->second );

        if ( this->worldComm().isMasterRank() )
            std::cout << "outVarDim " << outVarDim << "\n";

        if ( M_buildOnlyScalarSpace )
        {
            if ( outVarDim == 1 )
                this->runSaveFieldInTimeSet( M_spaceP1Scalar,fieldName, fieldIndex, 1, 1 );
            else if ( outVarDim == mesh_type::nRealDim )
                this->runSaveFieldInTimeSet( M_spaceP1Scalar,fieldName, fieldIndex, mesh_type::nRealDim, 1 );
            else if ( outVarDim == (mesh_type::nRealDim*mesh_type::nRealDim) )
                this->runSaveFieldInTimeSet( M_spaceP1Scalar,fieldName, fieldIndex, mesh_type::nRealDim, mesh_type::nRealDim );
        }
        else
        {
            if ( outVarDim == 1 )
                this->runSaveFieldInTimeSet( M_spaceP1Scalar,fieldName, fieldIndex );
            else if ( outVarDim == mesh_type::nRealDim )
                this->runSaveFieldInTimeSet( M_spaceP1Vectorial,fieldName, fieldIndex );
            else if ( outVarDim == (mesh_type::nRealDim*mesh_type::nRealDim) )
                CHECK( false )  << "not implemented";
        }
    }

    this->close();
}

template <typename MeshType>
template<typename SpaceType>
void
ConverterAcusimDatabase<MeshType>::runSaveFieldInTimeSet( std::shared_ptr<SpaceType> const& space, std::string const& fieldName, int fieldIndex )
{
    int outVarDim = space->nComponents;
    double * outFieldValues = new double[M_nNodes*outVarDim];

    auto u = space->element();

    for ( int k=0;k<M_nOutSteps;++k )
    {
        double time = M_outTimes[k];
        if ( this->worldComm().isMasterRank() )
            std::cout << "export timestep " << time << " ("<<k<<"/"<<M_nOutSteps<<")\n";

        this->applyFieldConversion( fieldIndex, k, outFieldValues, u );

        M_feelppDatabase.save( time,fieldName, u );
    }

    delete [] outFieldValues;
}

template <typename MeshType>
template<typename SpaceType>
void
ConverterAcusimDatabase<MeshType>::runSaveFieldInTimeSet( std::shared_ptr<SpaceType> const& space, std::string const& fieldName, int fieldIndex, int nComp1, int nComp2 )
{
    std::vector<std::string> compIdToCompName = { "X","Y","Z" };
    int nComp = nComp1*nComp2;
    double * outFieldValues = new double[M_nNodes*nComp];
    auto u = space->element();

    for ( int k=0;k<M_nOutSteps;++k )
    {
        double time = M_outTimes[k];
        if ( this->worldComm().isMasterRank() )
            std::cout << "export timestep " << time << " ("<<k<<"/"<<M_nOutSteps<<")\n";

        bool success = adbGetReals2( M_adbHd,(char*)"outExtValues", outFieldValues, fieldIndex, k/*timeStepIndex*/, nComp ,M_nNodes);

        for ( int c1=0;c1<nComp1;++c1 )
        {
            for ( int c2=0;c2<nComp2;++c2 )
            {
                this->applyFieldConversion( fieldIndex, k, outFieldValues, u, nComp1, nComp2, c1, c2 );

                std::string fieldNameInDb = fieldName;
                if ( nComp1 > 1 )
                    fieldNameInDb += "_" + compIdToCompName[ c1 ];
                if ( nComp2 > 1)
                    fieldNameInDb += "_"+ compIdToCompName[ c2 ];
                M_feelppDatabase.save( time,fieldNameInDb, u );
            }
        }
    }

    delete [] outFieldValues;
}




template <typename MeshType>
typename ConverterAcusimDatabase<MeshType>::mesh_ptrtype
ConverterAcusimDatabase<MeshType>::loadMesh()
{
    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;

#if 0
	int usrIds[nNodes];
	//usrIds = (int*) malloc(nNodes*sizeof(int));
    success = adbGetInts0( M_adbHd,(char*)"usrIds", usrIds, nNodes , 1);
#else
    bool success = adbGetInt0( M_adbHd, (char*)"nNodes",&M_nNodes);
	CHECK( success ) << "fail to get nNodes\n";
    M_nodeIds = new int[M_nNodes];
    success = adbGetInts0( M_adbHd,(char*)"usrIds", M_nodeIds, M_nNodes , 1);
    for ( int k=0;k<M_nNodes;++k )
        M_nodeIdToContainerId[ M_nodeIds[k] ] = k;
#endif

    mesh_ptrtype mesh;
    if ( !this->worldComm().isMasterRank() )
        return mesh;

    mesh.reset( new mesh_type( Environment::worldCommSeqPtr() ) );
    //-------------------------------------------------------------//
    // load mesh points
    //-------------------------------------------------------------//
	//int nNodes =0;
	//bool success = adbGetInt0  ( M_adbHd, (char*)"nNodes",&nNodes);

    double *crd = new double[M_nNodes*3];
	std::cout << "PO : load nNodes " << M_nNodes << "\n";
    success = adbGetReals0 ( M_adbHd,(char*)"crd", crd, 3 ,M_nNodes);


	// add points in mesh
	node_type coords( 3 );
    int ptid = 0;
    rank_type partId = mesh->worldComm().localRank();
	for( int i=0;i<M_nNodes;i++)
      {
		  coords[0]=crd[3*i];
		  coords[1]=crd[3*i+1];
		  coords[2]=crd[3*i+2];
		  // feel mesh point
		 point_type pt( /*ptid*/M_nodeIds[i], coords, false );
         pt.setProcessIdInPartition( partId );
         pt.setProcessId( partId );
         mesh->addPoint( pt );
	  }
    delete [] crd;

    //-------------------------------------------------------------//
    // load mesh elements
    //-------------------------------------------------------------//
    int eltid = 0;
    element_type e;
    e.setProcessIdInPartition( partId );
    e.setProcessId( partId );
    int markerId = 1;

    // get number of group of element
    int nElms = 0;
    success = adbGetInt0(M_adbHd,(char*)"nElms",&nElms);
    std::cout << "number of domain element " << nElms << "\n";
    int numberOfDomain = nElms; // seems necessary, to understand!!!
    int nElmElems = 0 , nElmElemNodes = 0;
    for (int j=0;j<numberOfDomain/*nElms*/;++j)
    {
        char elmName[1024] ;
        success = adbGetStr1(M_adbHd,(char*)"elmName", elmName,j);
        std::cout << "  Load element group "<< j << " with name : " << elmName << "\n";
        success = adbGetInt1(M_adbHd,(char*)"nElmElems",&nElmElems,j);
        std::cout << "  number of element : " << nElmElems << "\n";
        success = adbGetInt1(M_adbHd,(char*)"nElmElemNodes",&nElmElemNodes,j);
        std::cout << "  number of point in element : " << nElmElemNodes << "\n";
        CHECK( nElmElemNodes == element_type::numPoints ) << "number of point in element wrong " << nElmElemNodes <<" vs " << element_type::numPoints;

        // add marker
        std::string marker = Feel::detail::markerNameFromAcusimName( std::string(elmName), "tet4" );
        std::cout << "  marker element : " << marker << "\n";
        mesh->addMarkerName( marker, markerId, mesh_type::nDim );
        e.setMarker( markerId );


        // load element connectivity
        int * elmCnn = new int[nElmElems*nElmElemNodes];
        success = adbGetInts1(M_adbHd,(char*)"elmCnn",elmCnn,j,nElmElemNodes,nElmElems);
        // add elements in mesh
        for( int i=0;i<nElmElems;++i )
        {
            for ( uint16_type p = 0; p < element_type::numPoints; ++p )
            {
                int ptid = elmCnn[element_type::numPoints*i+p];
                e.setPoint( p, mesh->point( ptid ) );
            }
            mesh->addElement( e, true );
        }
        delete [] elmCnn;
        ++markerId;
    }


    //-------------------------------------------------------------//
    // load mesh faces
    //-------------------------------------------------------------//
    face_type f;
    f.setProcessIdInPartition( partId );
    f.setProcessId( partId );

    //int eltidConnected = 0, faceid = 0, faceidTemp = 0, ptids = 0;
    int faceid = 0;

    // get number of group of element
	int nOsfs = 0;
	success = adbGetInt0(M_adbHd,(char*)"nOsfs",&nOsfs);
    std::cout << "number of domain face " << nOsfs << "\n";
    int numberOfDomainFace = nOsfs; // seems necessary, to understand!!!
    int nOsfSrfs = 0, nOsfSrfNodes = 0;
	for( int j=0;j<numberOfDomainFace/*nOsfs*/;++j )
	{
        char osfName[1024] ;
        success = adbGetStr1(M_adbHd,(char*)"osfName",/*&*/osfName,j);
        std::cout << "  Load face group "<< j << " with name : " << osfName << "\n";
        success = adbGetInt1(M_adbHd,(char*)"nOsfSrfs",&nOsfSrfs,j);
        std::cout << "  number of face : " << nOsfSrfs << "\n";
        success = adbGetInt1(M_adbHd,(char*)"nOsfSrfNodes",&nOsfSrfNodes,j);
        std::cout << "  number of point in face : " << nOsfSrfNodes << "\n";
        CHECK( nOsfSrfNodes == face_type::numPoints ) << "number of point in element wrong " << nOsfSrfNodes <<" vs " << face_type::numPoints;

        std::string marker = Feel::detail::markerNameFromAcusimName( std::string(osfName), "tria3" );
        std::cout << "  marker face : " << marker << "\n";

        // add marker
        mesh->addMarkerName( marker, markerId, mesh_type::nDim - 1 );
        f.setMarker( markerId );


        // load face connectivity
        int * osfSrfCnn = new int[nOsfSrfNodes*nOsfSrfs];
        success = adbGetInts1(M_adbHd,(char*)"osfSrfCnn",osfSrfCnn,j,nOsfSrfNodes,nOsfSrfs);

        // add faces in mesh
        for( int i=0;i<nOsfSrfs;++i)
        {
            f.setId( faceid++ );
            for ( uint16_type p = 0; p < face_type::numPoints; ++p )
            {
                int ptid = osfSrfCnn[face_type::numPoints*i+p];
                f.setPoint( p, mesh->point( ptid ) );
            }
            mesh->addFace( f );
        }
        delete [] osfSrfCnn;
        ++markerId;
	}


    size_type updateCtx = size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED);
    mesh->components().reset();
    mesh->components().set( updateCtx );
    mesh->updateForUse();

    return mesh;
}

template <typename MeshType>
void
ConverterAcusimDatabase<MeshType>::loadTimeSetInfo()
{
    // number of time step
    M_nOutSteps = 0;
    bool success = adbGetInt0(M_adbHd, (char*)"nOutSteps", &M_nOutSteps);
    CHECK( success ) << "error on loading nOutSteps";
    if ( this->worldComm().isMasterRank() )
        std::cout << "nOutSteps : " << M_nOutSteps << "\n";
    // get time step index
    M_outSteps = new int[M_nOutSteps];
    success = adbGetInts0(M_adbHd,(char*)"outSteps",M_outSteps,M_nOutSteps,1);
    CHECK( success ) << "error on loading outSteps";
    // get time step values
    M_outTimes = new double[M_nOutSteps];
    success = adbGetReals0(M_adbHd,(char*)"outTimes",M_outTimes,M_nOutSteps,1);
    CHECK( success ) << "error on loading outTimes";
}

template <typename MeshType>
void
ConverterAcusimDatabase<MeshType>::loadFieldInfo()
{
    M_fieldNameToFieldInfo.clear();
    int nOutVars=0;
    bool success = adbGetInt0(M_adbHd,(char*)"nOutExtVars",&nOutVars);
    CHECK( success ) << "error on loading nOutExtVars";

    for (int i=0;i<nOutVars;i++)
    {
        char outVarName[1024];
        success = adbGetStr1(M_adbHd,(char*)"outExtVarName",outVarName,i);
        CHECK( success ) << "error on loading nOutExtVars";

        int outVarDim = 0;
        success = adbGetInt1(M_adbHd,(char*)"outExtVarDim",&outVarDim,i);
        CHECK( success ) << "error on loading outExtVarDim";

        M_fieldNameToFieldInfo[ std::string(outVarName) ] = std::make_tuple( i,outVarDim );
    }

    if ( this->worldComm().isMasterRank() )
    {
        std::cout << "---------------------------------------------\n"
                  << "Field available in Acusim Database :\n";
        for ( auto const& fieldInfo : M_fieldNameToFieldInfo )
            std::cout << " - " << fieldInfo.first << " [" << std::get<1>( fieldInfo.second ) << "] \n";
        std::cout << "---------------------------------------------\n";
    }
}

template <typename MeshType>
template<typename SpaceType>
void
ConverterAcusimDatabase<MeshType>::generateDofMappingP1( std::shared_ptr<SpaceType> const& space )
{
    auto mesh = space->mesh();
    auto doftable = space->dof();
    uint16_type nComp = doftable->nComponents;

    auto & feelDofIdToAcusimDofId = M_feelDofIdToAcusimDofId[nComp];
    feelDofIdToAcusimDofId.resize( space->nLocalDofWithGhost() );
    const size_type nLocDof = space->dof()->nLocalDof(true);
    for (auto const& eltWrap : elements(mesh) )
    {
        auto const& elt = boost::unwrap_ref( eltWrap );
        const size_type eltId = elt.id();
        for ( uint16_type j=0 ; j<nLocDof ; ++j )
        {
            auto const& thepoint = elt.point(j);
            size_type ptId = thepoint.id();
            auto itFindPt = M_nodeIdToContainerId.find( ptId );
            CHECK( itFindPt != M_nodeIdToContainerId.end() ) << "pt not found";
            size_type indexDataLoaded = itFindPt->second;
            for ( uint16_type comp=0 ; comp<nComp ; ++comp )
            {
                size_type dofIndex = doftable->localToGlobal( eltId, j , comp ).index();
#if 0
                auto const& nodeAtNode = doftable->dofPoint( dofIndex ).template get<0>();
                CHECK( ( std::abs(nodeAtNode[0] - thepoint[0]) < 1e-9 ) &&
                       ( std::abs(nodeAtNode[1] - thepoint[1]) < 1e-9 ) &&
                       ( std::abs(nodeAtNode[2] - thepoint[2]) < 1e-9 ) ) << "points and dof must be same : " << nodeAtNode << " vs " << thepoint;
#endif
                feelDofIdToAcusimDofId[ dofIndex ] = indexDataLoaded*nComp+ comp;
            }
        }
    }
}


template <typename MeshType>
template<typename ElementType>
void
ConverterAcusimDatabase<MeshType>::applyFieldConversion( int fieldIndex, int timeStepIndex, double * outFieldValues, ElementType & u )
{
    auto space = u.functionSpace();
    int outVarDim = space->nComponents;
    auto itFindDofMapping = M_feelDofIdToAcusimDofId.find( outVarDim );
    CHECK( itFindDofMapping != M_feelDofIdToAcusimDofId.end() ) << "mapping dof not found with nComp " << outVarDim;
    auto const& feelDofIdToAcusimDofId = itFindDofMapping->second;
    bool success = adbGetReals2( M_adbHd,(char*)"outExtValues", outFieldValues, fieldIndex,timeStepIndex, outVarDim ,M_nNodes);
    for ( size_type k=0;k<space->nLocalDofWithGhost() ;++k )
        u.set( k, outFieldValues[ feelDofIdToAcusimDofId[k] ] );
}


template <typename MeshType>
template<typename ElementType>
void
ConverterAcusimDatabase<MeshType>::applyFieldConversion( int fieldIndex, int timeStepIndex, double * outFieldValues, ElementType & u, int nComp1, int nComp2, int comp1, int comp2 )
{
    int nComp = nComp1*nComp2;
    auto space = u.functionSpace();
    int nCompInSpace = space->nComponents;
    auto itFindDofMapping = M_feelDofIdToAcusimDofId.find( nCompInSpace );
    CHECK( itFindDofMapping != M_feelDofIdToAcusimDofId.end() ) << "mapping dof not found with nComp " << nCompInSpace;
    auto const& feelDofIdToAcusimDofId = itFindDofMapping->second;

    //bool success = adbGetReals2( M_adbHd,(char*)"outExtValues", outFieldValues, fieldIndex,timeStepIndex, nComp ,M_nNodes);
    for ( size_type k=0;k<space->nLocalDofWithGhost() ;++k )
        u.set( k, outFieldValues[ feelDofIdToAcusimDofId[k]*nComp+comp1*nComp2+comp2  ] );
    //u.set( k, outFieldValues[ feelDofIdToAcusimDofId[k]*nComp+comp1  ] );
}



} // namespace Feel

#endif // defined( FEELPP_HAS_ACUSIM )

#endif // FEELPP_CONVERTERACUSIMDATABASE_HPP
