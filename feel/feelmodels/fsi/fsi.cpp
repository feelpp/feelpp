/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-07-05

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
   \file fsi.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-07-05
 */

#include <feel/feelmodels/fsi/fsi.hpp>
#include <feel/feelmodels/modelmesh/fsimesh.hpp>

namespace Feel
{
namespace FeelModels
{

template< class FluidType, class SolidType >
FSI<FluidType,SolidType>::FSI(std::string prefix,WorldComm const& worldComm )
    :
    super_type( prefix, worldComm ),
    M_meshSize( doption(_name="hsize",_prefix=this->prefix()) ),
    M_tagFileNameMeshGenerated( soption(_name="mesh-save.tag",_prefix=this->prefix()) ),
    M_fsiCouplingType( soption(_name="coupling-type",_prefix=this->prefix()) ),
    M_fsiCouplingBoundaryCondition( soption(_name="coupling-bc",_prefix=this->prefix()) ),
    M_interfaceFSIisConforme( boption(_name="conforming-interface",_prefix=this->prefix()) ),
    M_tolPtFixe( doption(_name="fixpoint.tol",_prefix=this->prefix()) ),
    M_initialTheta( doption(_name="fixpoint.initialtheta",_prefix=this->prefix()) ),
    M_minTheta( doption(_name="fixpoint.min_theta",_prefix=this->prefix()) ),
    M_fixPointMaxIt( ioption(_name="fixpoint.maxit",_prefix=this->prefix()) ),
    M_previousTimeOrder(0),M_currentTimeOrder(1),
    M_reusePrecOptFluid(false),
    M_reusePrecRebuildAtFirstFSIStepOptFluid( boption(_name="fluid.reuse-prec.rebuild-at-first-fsi-step",_prefix=this->prefix()) ),
    M_reuseJacOptFluid(false),
    M_reuseJacRebuildAtFirstNewtonStepOptFluid(false),
    M_reuseJacRebuildAtFirstFSIStepOptFluid( boption(_name="fluid.reuse-jac.rebuild-at-first-fsi-step",_prefix=this->prefix()) ),
    M_reusePrecOptSolid(false),
    M_reusePrecRebuildAtFirstFSIStepOptSolid( boption(_name="solid.reuse-prec.rebuild-at-first-fsi-step",_prefix=this->prefix()) ),
    M_reuseJacOptSolid(false),
    M_reuseJacRebuildAtFirstNewtonStepOptSolid(false),
    M_reuseJacRebuildAtFirstFSIStepOptSolid( boption(_name="solid.reuse-jac.rebuild-at-first-fsi-step",_prefix=this->prefix()) )
{
    this->log("fsi","constructor","start");

    if ( Environment::vm().count( prefixvm(this->prefix(),"fluid-mesh.markers" ) ) )
    {
        std::vector<std::string> mymarkers = Environment::vm()[prefixvm(this->prefix(),"fluid-mesh.markers").c_str()].template as<std::vector<std::string> >();
        for ( std::string const& marker : mymarkers )
            M_markersNameFluid.insert( marker );
    }
    if ( Environment::vm().count( prefixvm(this->prefix(),"solid-mesh.markers" ) ) )
    {
        std::vector<std::string> mymarkers = Environment::vm()[prefixvm(this->prefix(),"solid-mesh.markers").c_str()].template as<std::vector<std::string> >();
        for ( std::string const& marker : mymarkers )
            M_markersNameSolid.insert( marker );
    }

    this->log("fsi","constructor","finish");
}
//---------------------------------------------------------------------------------------------------------//

template <typename FluidType,typename SolidType>
void
FSI<FluidType,SolidType>::createMesh()
{
    this->log("FSI","createMesh","start");

    FSIMesh<typename fluid_type::convex_type> fsimeshTool(this->prefix(),this->worldComm());

    int nPart = this->worldComm().size();//this->nPartitions();

    fs::path meshesdirectories = fs::path(this->appliRepositoryWithoutNumProc()) / fs::path("meshes");
    fs::path gp;
    if (this->hasMshfileStr())
        gp = this->mshfileStr();
    else
        gp = this->geofileStr();

    std::string nameMeshFile = gp.stem().string();
    fs::path mshFileNameFluidPart1, mshFileNameSolidPart1,mshFileNameFluidPartN,mshFileNameSolidPartN;
    fs::path mshFileNameFSI;
    if ( !M_tagFileNameMeshGenerated.empty() )
    {
        mshFileNameFluidPart1 = (boost::format("%1%_fluid_%2%_M%3%_p1.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated ).str();
        mshFileNameSolidPart1 = (boost::format("%1%_solid_%2%_M%3%_p1.msh") %nameMeshFile %solid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated ).str();
        mshFileNameFluidPartN = (boost::format("%1%_fluid_%2%_M%3%_p%4%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %nPart ).str();
        mshFileNameSolidPartN = (boost::format("%1%_solid_%2%_M%3%_p%4%.msh") %nameMeshFile %solid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %nPart ).str();
    }
    else
    {
        mshFileNameFluidPart1 = (boost::format("%1%_fluid_%2%_p1.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() ).str();
        mshFileNameSolidPart1 = (boost::format("%1%_solid_%2%_p1.msh") %nameMeshFile %solid_type::mesh_type::shape_type::name() ).str();
        mshFileNameFluidPartN = (boost::format("%1%_fluid_%2%_p%3%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %nPart ).str();
        mshFileNameSolidPartN = (boost::format("%1%_solid_%2%_p%3%.msh") %nameMeshFile %solid_type::mesh_type::shape_type::name() %nPart ).str();
    }
    M_mshfilepathFluidPart1 = meshesdirectories / mshFileNameFluidPart1;
    M_mshfilepathSolidPart1 = meshesdirectories / mshFileNameSolidPart1;
    M_mshfilepathFluidPartN = meshesdirectories / mshFileNameFluidPartN;
    M_mshfilepathSolidPartN = meshesdirectories / mshFileNameSolidPartN;

    fsimeshTool.setNumberOfPartitions( this->worldComm().localSize() );
    fsimeshTool.setFluidMshPathPart1( M_mshfilepathFluidPart1 );
    fsimeshTool.setSolidMshPathPart1( M_mshfilepathSolidPart1 );
    fsimeshTool.setFluidMshPathPartN( M_mshfilepathFluidPartN );
    fsimeshTool.setSolidMshPathPartN( M_mshfilepathSolidPartN );
    fsimeshTool.setMarkersNameFluidVolume( M_markersNameFluid );
    fsimeshTool.setMarkersNameSolidVolume( M_markersNameSolid );


    if (this->hasMshfileStr())
    {
        fsimeshTool.setMshPathFSI( fs::path(this->mshfileStr()) );
        fsimeshTool.buildFSIMeshFromMsh();
    }
    else if (this->hasGeofileStr())
    {
        fs::path mshFileNameFSI;
        if ( !M_tagFileNameMeshGenerated.empty() )
            mshFileNameFSI = (boost::format("%1%_fsi_%2%_M%3%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated ).str();
        else
            mshFileNameFSI = (boost::format("%1%_fsi_%2%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() ).str();
        fs::path mshPathFSI = meshesdirectories / mshFileNameFSI;

        fsimeshTool.setGeoPathFSI( fs::path(this->geofileStr()) );
        fsimeshTool.setMshPathFSI( mshPathFSI );
        fsimeshTool.setMeshSize( this->meshSize() );
        fsimeshTool.buildFSIMeshFromGeo();
    }

    this->log("FSI","createMesh","finish");
}

//---------------------------------------------------------------------------------------------------------//

namespace detail
{

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh2d( typename FluidType::self_ptrtype const& FM, mpl::bool_<false> /**/ )
{
    auto submeshStruct = createSubmesh( FM->getMeshALE()->referenceMesh(), markedfaces( FM->getMeshALE()->referenceMesh(), FM->markersNameMovingBoundary()/*"Paroi"*/) );
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
      submeshStruct->faces().modify( submeshStruct->faceIterator( itp->id() ) , Feel::detail::UpdateMarker( submeshStruct->markerName("Fixe") ) );

    typedef SubMeshData smd_type;
    typedef boost::shared_ptr<smd_type> smd_ptrtype;
    smd_ptrtype smd( new smd_type(FM->mesh()) );
    for ( auto const& e : elements(submeshStruct) )
      {
        auto const& theface = FM->getMeshALE()->referenceMesh()->face( submeshStruct->subMeshToMesh(e.id()) );
        size_type idElt2 = FM->getMeshALE()->dofRelationShipMap()->geoElementMap()[ theface.element0().id() ].first;
        //std::cout << " e.G() " << e.G() << " other.G() " <<  theface.G() << std::endl;
        auto const& theface2 = FM->mesh()->element(idElt2,e.processId()).face(theface.pos_first());
        smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), theface2.id() ) );
      }
    submeshStruct->setSubMeshData( smd );

    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename FluidType::mesh_type::trace_mesh_ptrtype
createMeshStruct1dFromFluidMesh2d( typename FluidType::self_ptrtype const& FM, mpl::bool_<true> /**/ )
{
    auto submeshStruct = createSubmesh( FM->mesh(), markedfaces( FM->mesh(),FM->markersNameMovingBoundary() ) );
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
      submeshStruct->faces().modify( submeshStruct->faceIterator( itp->id() ) , Feel::detail::UpdateMarker( submeshStruct->markerName("Fixe") ) );

    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, mpl::int_<2> /**/ )
{
  static const bool hasSameOrderGeo = FluidType::mesh_type::nOrder == SolidType::mesh_1dreduced_type::nOrder;
  return createMeshStruct1dFromFluidMesh2d<FluidType,SolidType>(FM, mpl::bool_<hasSameOrderGeo>() );
}

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, mpl::int_<3> /**/ )
{
    CHECK( false ) << "not possible";
    return typename SolidType::mesh_1dreduced_ptrtype();
}

template <typename FluidType,typename SolidType>
typename SolidType::mesh_1dreduced_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM )
{
  return createMeshStruct1dFromFluidMesh<FluidType,SolidType>( FM, mpl::int_<SolidType::nDim>() );
}

} // namespace detail

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::init()
{
    this->log("ToolBoxFSI","init","start");

    // create fsimesh and partitioned meshes if require
    if ( this->hasMshfileStr() || this->hasGeofileStr() )
        this->createMesh();

    // fluid model build
    if ( !M_fluid )
    {
        M_fluid = fluid_ptrtype( new fluid_type("fluid",false,this->worldComm() ) );
        if ( !M_mshfilepathFluidPartN.empty() )
            M_fluid->setMshfileStr(M_mshfilepathFluidPartN.string());
        M_fluid->build();
    }

    // solid model build
    if ( !M_solid )
    {
        M_solid = solid_ptrtype( new solid_type("solid",false,this->worldComm() ) );
        bool doExtractSubmesh = boption(_name="solid-mesh.extract-1d-from-fluid-mesh",_prefix=this->prefix() );
        if ( doExtractSubmesh )
        {
            CHECK( !M_fluid->markersNameMovingBoundary().empty() ) << "no marker moving boundary in fluid model";

            if ( M_fluid->doRestart() )
                M_fluid->getMeshALE()->revertReferenceMesh();
            auto submeshStruct = detail::createMeshStruct1dFromFluidMesh<fluid_type,solid_type>( M_fluid );
            if ( M_fluid->doRestart() )
                M_fluid->getMeshALE()->revertMovingMesh();

            // TODO ( save 1d mesh and reload )
            if ( M_fluid->doRestart() )
                M_solid->build(submeshStruct);
            else
                M_solid->loadMesh(submeshStruct);

        }
        else
        {
            if ( !M_mshfilepathSolidPartN.empty() )
                M_solid->setMshfileStr( M_mshfilepathSolidPartN.string() );
            M_solid->build();
        }
    }

    CHECK( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
           this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
           this->fsiCouplingBoundaryCondition() == "nitsche" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" ) << "invalid fsiCouplingBoundaryCondition : " << this->fsiCouplingBoundaryCondition();
    M_fluid->couplingFSIcondition(this->fsiCouplingBoundaryCondition());
    M_solid->couplingFSIcondition(this->fsiCouplingBoundaryCondition());

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        M_fluid->useFSISemiImplicitScheme(true);
        M_solid->useFSISemiImplicitScheme(true);
    }

    // specific value for robin
    M_solid->muFluidFSI( M_fluid->densityViscosityModel()->cstMu() );
    double gammaNitsche = doption(_name="coupling-robin-robin.gamma",_prefix=this->prefix());
    M_fluid->gammaNitschFSI( gammaNitsche );
    M_solid->gammaNitschFSI( gammaNitsche );
    double gamma0Nitsche = doption(_name="coupling-robin-robin.gamma0",_prefix=this->prefix());
    M_fluid->gamma0NitschFSI( gamma0Nitsche );

    // save if reuse prec option at the begining
    M_reusePrecOptFluid = M_fluid->backend()->reusePrec();
    M_reuseJacOptFluid = M_fluid->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptFluid = M_fluid->backend()->reuseJacRebuildAtFirstNewtonStep();
    M_reusePrecOptSolid = M_solid->backend()->reusePrec();
    M_reuseJacOptSolid = M_solid->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptSolid = M_solid->backend()->reuseJacRebuildAtFirstNewtonStep();
    //-------------------------------------------------------------------------//
    // call init
    M_fluid->init();
    M_solid->init();
    //-------------------------------------------------------------------------//
    if (M_solid->isStandardModel())
    {
        M_fluid->setCouplingFSI_solidIs1dReduced(false);
        M_fluid->setRhoSolidFSI( M_solid->mechanicalProperties()->cstRho() );
    }
    else if ( M_solid->is1dReducedModel() )
    {
        M_fluid->setCouplingFSI_solidIs1dReduced(true);
        M_fluid->setRhoSolidFSI( M_solid->mechanicalProperties()->cstRho()*0.1 ); // 0.1 is thickness (TODO)
    }

    //M_solid->createAdditionalFunctionSpacesFSI();
    //-------------------------------------------------------------------------//
    // init interpolation tool
    M_interpolationFSI.reset(new interpolationFSI_type(M_fluid,M_solid));

    if ( M_solid->isStandardModel() )
        M_interpolationFSI->updateFieldVelocitySolidPreviousPrevious( M_solid->timeStepNewmark()->previousVelocity() );
    else if ( M_solid->is1dReducedModel() )
        M_interpolationFSI->updateFieldVelocitySolid1dReducedPreviousPrevious( M_solid->timeStepNewmark1dReduced()->previousVelocity() );

    //-------------------------------------------------------------------------//
    // init aitken relaxation tool
    //AitkenType aitkenType = AitkenType::AITKEN_METHOD_1;
    std::string aitkenType = "method1";
    if (this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
        this->fsiCouplingBoundaryCondition()=="robin-neumann" || this->fsiCouplingBoundaryCondition()=="robin-neumann-genuine" ||
        this->fsiCouplingBoundaryCondition()=="robin-neumann-generalized" ||
        this->fsiCouplingBoundaryCondition() == "nitsche")
    {
        //aitkenType = AitkenType::FIXED_RELAXATION_METHOD;
        aitkenType = "fixed-relaxation";
        M_initialTheta=0;
        M_fixPointConvergenceFSI.reset(new fixpointconvergenceFSI_type(M_solid) );
    }
    M_aitkenFSI.reset(new aitkenrelaxationFSI_type(M_solid,
                                                   aitkenType,//AITKEN_METHOD_1,
                                                   M_initialTheta,
                                                   M_tolPtFixe,
                                                   M_minTheta));

    //-------------------------------------------------------------------------//
    // build interface operator for generalized robin-neumann
    if (this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized")
    {
        auto fieldInit = M_fluid->meshVelocity2().functionSpace()->elementPtr();
        M_fluid->setCouplingFSI_RNG_evalForm1( fieldInit );
    }
    if (M_solid->isStandardModel() && this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" &&
        boption(_name="coupling-robin-neumann-generalized.use-interface-operator",_prefix=this->prefix() ) )
    {
#if 1
        auto spaceDisp = M_solid->functionSpaceDisplacement();
        auto const& uDisp = M_solid->fieldDisplacement();
        auto mygraph = stencil(_test=spaceDisp,_trial=spaceDisp)->graph();

        //----------------------------------------------------------------------------------//
        // get dofs on markedfaces fsi
        std::set<size_type> dofMarkerFsi;
        for ( auto const& faceMarked : markedfaces(M_solid->mesh(),M_solid->markerNameFSI() ) )
        {
            auto __face_it = faceMarked.template get<1>();
            auto __face_en = faceMarked.template get<2>();
            for( ; __face_it != __face_en; ++__face_it )
            {
                auto const& face = *__face_it;
                for ( uint16_type l = 0; l < M_solid->functionSpaceDisplacement()->dof()->nLocalDofOnFace(true); ++l )
                {
                    for (uint16_type c1=0;c1<solid_type::nDim;c1++)
                    {
                        size_type gdof = boost::get<0>(M_solid->functionSpaceDisplacement()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                        dofMarkerFsi.insert( gdof );
                    }
                }
            }
        }
#if 0 // bug
        std::set<size_type> dofMarkerFsiTest;
        for ( std::string const markName : M_solid->markerNameFSI() )
        {
            //functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 )
            auto setofdof = M_solid->functionSpaceDisplacement()->dof()->markerToDof( markName );
            for( auto it = setofdof.first, en = setofdof.second; it != en; ++ it )
                dofMarkerFsiTest.insert( it->second );
        }
        std::cout << "size " << dofMarkerFsiTest.size() << " vs " << dofMarkerFsi.size() << "\n";
#endif

        //----------------------------------------------------------------------------------//
        // mass matrix on solid
        auto matMass = M_solid->backend()->newMatrix(0,0,0,0,mygraph);
        form2( _trial=spaceDisp, _test=spaceDisp,_matrix=matMass )
            += integrate(
                //_range=markedfaces(M_solid->mesh(),M_solid->markerNameFSI()),
                _range=elements(M_solid->mesh()),
                //_quad=_Q<1>(),
                _expr=/*M_solid->mechanicalProperties()->cstRho()**/inner(idt(uDisp),id(uDisp)) );
        matMass->close();
        auto vecDiag = M_solid->backend()->newVector(spaceDisp);

        auto vecDiagMassLumped = M_solid->backend()->newVector(spaceDisp);
        auto unityVec = M_solid->backend()->newVector(spaceDisp);
        for (int k=0;k<unityVec->map().nLocalDofWithGhost();++k)
            unityVec->set(k,1);
        unityVec->close();
        matMass->multVector( unityVec,vecDiagMassLumped );
        vecDiagMassLumped->close();

        //----------------------------------------------------------------------------------//
        // mass matrix on interface
        auto matMassInterface = M_solid->backend()->newMatrix(0,0,0,0,mygraph);
        form2( _trial=spaceDisp, _test=spaceDisp,_matrix=matMassInterface )
            += integrate( _range=markedfaces(M_solid->mesh(),M_solid->markerNameFSI()),
                          //_quad=_Q<1>(),
                          _expr=inner(idt(uDisp),id(uDisp)) );
        matMassInterface->close();
        auto vecDiagMassInterface = M_solid->backend()->newVector(spaceDisp);
        M_solid->backend()->diag(matMassInterface,vecDiagMassInterface );
        vecDiagMassInterface->close();


        auto unityVecInterface = M_solid->backend()->newVector(spaceDisp);
        for (int k=0;k<unityVecInterface->map().nLocalDofWithGhost();++k)
            if ( dofMarkerFsi.find( k ) != dofMarkerFsi.end() ) 
                unityVecInterface->set(k,1);
        unityVecInterface->close();
        auto vecDiagMassInterfaceLumped = M_solid->backend()->newVector(spaceDisp);
        matMassInterface->multVector( unityVecInterface/*unityVec*/,vecDiagMassInterfaceLumped );
        vecDiagMassInterfaceLumped->close();


        auto bbbb = sum(spaceDisp,elements(spaceDisp->mesh()),cst(1./4.)*meas()*one());
        for ( size_type k : dofMarkerFsi )
        {
            if ( vecDiag->map().dofGlobalProcessIsGhost( k ) ) continue;
            size_type gcdof = vecDiag->map().mapGlobalProcessToGlobalCluster(k);
            auto const& graphRow = matMass->graph()->row(gcdof);
            double sumCol=0,sumCol2=0;
            for ( size_type idCol : graphRow.template get<2>() )
            {
                sumCol += matMass->operator()(k,idCol);
                if ( k==idCol ) std::cout << "diag eval mass " << matMass->operator()(k,idCol) << "\n";

                if ( dofMarkerFsi.find( idCol ) != dofMarkerFsi.end() ) 
                    sumCol2 += matMass->operator()(k,idCol);

            }
            std::cout << "val interface op " << sumCol << " VS " << sumCol2 << " VS " << bbbb(k) << " VS " << vecDiagMassLumped->operator()(k)
                      << " Interface " << matMassInterface->operator()(k,k) << " VS " << vecDiagMassInterfaceLumped->operator()(k)
                      << " in row " << k << "with "<< graphRow.template get<2>().size() <<" dof in row\n";
            //vecDiag->set( k, sumCol );
            //vecDiag->set( k, vecDiagMassLumped->operator()(k)/matMassInterface->operator()(k,k) );///THE NEWS
            vecDiag->set( k, vecDiagMassLumped->operator()(k)/vecDiagMassInterface->operator()(k) );///THE NEWS
            std::cout << "use " << vecDiagMassLumped->operator()(k) << " and " << vecDiagMassInterface->operator()(k)
                      << " = " << vecDiagMassLumped->operator()(k)/vecDiagMassInterface->operator()(k) << "\n";
            //vecDiag->set( k, vecDiagMassLumped->operator()(k)/vecDiagMassInterfaceLumped->operator()(k) );///THE NEWS
        }
        vecDiag->close();
#endif

        /*auto bbbb = sum(spaceDisp,elements(spaceDisp->mesh()),cst(1./4.)*meas()*one());
        auto vecDiag = M_solid->backend()->newVector(spaceDisp);
        *vecDiag = bbb;
         vecDiag.close();*/


        M_interpolationFSI->setRobinNeumannInterfaceOperator( vecDiag );
        M_interpolationFSI->transfertRobinNeumannInterfaceOperatorS2F();
    }
    //-------------------------------------------------------------------------//

    this->log("ToolBoxFSI","init","finish");
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solve()
{
    boost::mpi::timer mytimer;

    if ( this->fsiCouplingBoundaryCondition()=="dirichlet-neumann" )
    {
        this->solveImpl1();
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
              this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
              this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        this->solveImpl2();
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        this->solveImpl3();
    }

    double timeElapsed = mytimer.elapsed();
    if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
        std::cout << "["<<prefixvm(this->prefix(),"FSI") <<"] finish fsi solve in " << timeElapsed << "\n";
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateBackendOptimisation( bool restartFullStepFluid,bool restartFullStepSolid )
{
    // ici on impose la reconstruction du preconditionneur au debut du solveur FSI
    // maybe add an option to do that
    if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid) )
    {
        M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(restartFullStepFluid);
    }

    if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
    {
        M_fluid->backend()->setReusePrec(!restartFullStepFluid);
    }

    if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
    {
        M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(restartFullStepSolid);
    }

    if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
    {
        M_solid->backend()->setReusePrec(!restartFullStepSolid);
    }
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl1()
{
    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();

        M_fluid->setRebuildLinearPartInJacobian(true);M_fluid->setRebuildCstPartInLinearSystem(true);
        M_solid->setRebuildLinearPartInJacobian(true);M_solid->setRebuildCstPartInLinearSystem(true);
        M_fluid->setRebuildCstPartInResidual(true);M_solid->setRebuildCstPartInResidual(true);
    }
#if 0
    // ici on impose la reconstruction du preconditionneur au debut du solveur FSI
    // maybe add an option to do that
    if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid) )
    {
        M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(true);
    }

    if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
    {
        M_fluid->backend()->setReusePrec(false);
    }

    if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
    {
        M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(true);
    }

    if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
    {
        M_solid->backend()->setReusePrec(false);
    }
#else
    this->updateBackendOptimisation(true,true);
#endif
    // predictor disp
    M_solid->predictorDispl();
    // coupling fluid structure
    this->aitkenRelaxTool()->restart();
    //aitkenFSI.setTheta(initialTheta);

    boost::mpi::timer timerCur,timerIter;

    while (!this->aitkenRelaxTool()->isFinished() && this->aitkenRelaxTool()->nIterations()<this->fixPointMaxIt())
    {
        timerIter.restart();
        timerCur.restart();
        //--------------------------------------------------------------//

        this->aitkenRelaxTool()->saveOldSolution();

        //--------------------------------------------------------------//
        // ALE solver
        timerCur.restart();
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluid->updateALEmesh();
            double tALE = timerCur.elapsed();
            this->log("main","update ale","finish in "+(boost::format("%1% s") % tALE).str() );
        }
        else  if (this->fsiCouplingType()=="Semi-Implicit")
        {
            this->interpolationTool()->transfertVelocity();
            double tALE = timerCur.elapsed();
            this->log("main","transfert velocity","finish in "+(boost::format("%1% s") % tALE).str());
        }
        // robin-neumann as fsi condition
        if (this->fsiCouplingBoundaryCondition()=="robin-neumann")
        {
            this->interpolationTool()->transfertStressS2F();
        }

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        M_fluid->solve();
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//

        //--------------------------------------------------------------//
        timerCur.restart();
        // revert ref mesh
        M_fluid->getMeshALE()->revertReferenceMesh();
        // transfert stress
        this->interpolationTool()->transfertStress();
        // revert moving mesh
        M_fluid->getMeshALE()->revertMovingMesh();


        //if (this->fsiCouplingBoundaryCondition()=="robin-robin")
        //    this->interpolationTool()->transfertVelocityF2S(true);

        double t3 = timerCur.elapsed();
        this->log("main","transfert stress","finish in "+(boost::format("%1% s") % t3).str() );
        //--------------------------------------------------------------//

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        M_solid->solve();
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
#if 0
        M_fluid->updateTime(M_fluid->time()+1);
        M_fluid->exportResults();
        M_solid->updateTime(M_solid->time()+1);
        M_solid->exportResults();
#endif

        //--------------------------------------------------------------//
        timerCur.restart();
        //compute and apply aitken relaxation
        //if ( this->fsiCouplingBoundaryCondition()=="dirichlet-neumann" )
        this->aitkenRelaxTool()->applyRelaxation();
        // update velocity and acceleration
        M_solid->updateVelocity();
        // aitken relaxtion
        if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
            this->aitkenRelaxTool()->printInfo();
        this->aitkenRelaxTool()->shiftRight();

        double t4 = timerCur.elapsed();
        this->log("main","apply relax and up vel/acc struct","finish in "+(boost::format("%1% s") % t4).str() );
        //--------------------------------------------------------------//

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluid->setRebuildLinearPartInJacobian(false);M_fluid->setRebuildCstPartInLinearSystem(false);
            M_solid->setRebuildLinearPartInJacobian(false);M_solid->setRebuildCstPartInLinearSystem(false);
            M_fluid->setRebuildCstPartInResidual(false);M_solid->setRebuildCstPartInResidual(false);
        }
        /*else if (this->fsiCouplingType()=="Implicit" &&
         M_reuseJacOptFluid && !M_reuseJacRebuildAtFirstNewtonStepOptFluid &&
         !M_fluid->useLinearJacobianInResidual() )
         {
         M_fluid->setRebuildLinearPartInJacobian(false);
         }*/
#if 0
        if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid ) )
        {
            M_fluid->backend()->setReuseJacRebuildAtFirstNewtonStep(false);
        }
        if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
        {
            M_fluid->backend()->setReusePrec(true);
        }

        if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
        {
            M_solid->backend()->setReuseJacRebuildAtFirstNewtonStep(false);
        }
        if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
        {
            M_solid->backend()->setReusePrec(true);
        }
#else
        this->updateBackendOptimisation(false,false);
#endif
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("main","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
    }

}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl2()
{
    boost::mpi::timer timerIter;

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
#if 0
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();
#endif

        M_fluid->setRebuildLinearPartInJacobian(true);M_fluid->setRebuildCstPartInLinearSystem(true);
        M_solid->setRebuildLinearPartInJacobian(true);M_solid->setRebuildCstPartInLinearSystem(true);
        M_fluid->setRebuildCstPartInResidual(true);M_solid->setRebuildCstPartInResidual(true);
    }

    this->updateBackendOptimisation(true,true);

    /*if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        M_solid->updateSubMeshDispFSIFromPrevious();
     }*/


    int cptFSI=0;
    double residualConvergence=1;

    bool solveStruct = this->fluidModel()->timeStepBDF()->iteration() > 1 || this->fixPointMaxIt()==1;

    while ( ( residualConvergence > this->tolPtFixe() || cptFSI <3 ) && cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        M_fixPointConvergenceFSI->saveOldSolution();

        //--------------------------------------------------------------//
        if (solveStruct)
        {
            M_fluid->getMeshALE()->revertReferenceMesh();
            // transfert stress
            this->interpolationTool()->transfertStress();
            // revert moving mesh
            M_fluid->getMeshALE()->revertMovingMesh();
            if ( ( this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" ||
                   this->fsiCouplingBoundaryCondition()=="nitsche" ) &&
                 M_solid->isStandardModel() )
            {
                bool useExtrap = boption(_prefix=this->prefix(),_name="transfert-velocity-F2S.use-extrapolation");
                this->interpolationTool()->transfertVelocityF2S(cptFSI,useExtrap);
            }
            M_solid->solve();
            M_solid->updateVelocity();
        }
        //--------------------------------------------------------------//
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluid->updateALEmesh();
        }
        else
        {
            this->interpolationTool()->transfertVelocity();
        }
        M_fluid->solve();
        //--------------------------------------------------------------//

        residualConvergence = M_fixPointConvergenceFSI->computeConvergence();

        if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
            std::cout << "["<<prefixvm(this->prefix(),"FSI") <<"] iteration " << cptFSI
                      << " residualConvergence " << std::scientific << residualConvergence << "\n";

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluid->setRebuildLinearPartInJacobian(false);M_fluid->setRebuildCstPartInLinearSystem(false);
            M_fluid->setRebuildCstPartInResidual(false);

            if (solveStruct)
            {
                M_solid->setRebuildLinearPartInJacobian(false);M_solid->setRebuildCstPartInLinearSystem(false);
                M_solid->setRebuildCstPartInResidual(false);
            }
        }

        this->updateBackendOptimisation(false,!solveStruct);


        if (!solveStruct) solveStruct=true;

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("main","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");

        ++cptFSI;
    }

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();
    }


}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl3()
{
    boost::mpi::timer timerIter;

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluid->updateALEmesh();

        M_fluid->setRebuildLinearPartInJacobian(true);M_fluid->setRebuildCstPartInLinearSystem(true);
        M_solid->setRebuildLinearPartInJacobian(true);M_solid->setRebuildCstPartInLinearSystem(true);
        M_fluid->setRebuildCstPartInResidual(true);M_solid->setRebuildCstPartInResidual(true);
    }

    this->updateBackendOptimisation(true,true);

    double manualScalingRNG = doption(_name="coupling-robin-neumann-generalized.manual-scaling",_prefix=this->prefix());

    int cptFSI=0;
    double residualConvergence=1;

    //bool solveStruct = this->fluidModel()->timeStepBDF()->iteration() > 1 || this->fixPointMaxIt()==1;

    while ( ( residualConvergence > this->tolPtFixe() || cptFSI <3 ) && cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        //timerCur.restart();
        M_fixPointConvergenceFSI->saveOldSolution();
        //--------------------------------------------------------------//
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluid->updateALEmesh();
        }
        else
        {
            //bool useExtrap = (this->fluidModel()->timeStepBDF()->iteration() > 1) && cptFSI==0;
            //this->interpolationTool()->transfertVelocity(useExtrap);
        }
        this->interpolationTool()->transfertRobinNeumannGeneralizedS2F( cptFSI, manualScalingRNG );
        M_fluid->solve();
        //--------------------------------------------------------------//
        M_fluid->getMeshALE()->revertReferenceMesh();
        // transfert stress
        this->interpolationTool()->transfertStress();
        // revert moving mesh
        M_fluid->getMeshALE()->revertMovingMesh();
        M_solid->solve();
        M_solid->updateVelocity();
        //--------------------------------------------------------------//

        residualConvergence = M_fixPointConvergenceFSI->computeConvergence();
        if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
            std::cout << "["<<prefixvm(this->prefix(),"FSI") <<"] iteration " << cptFSI
                      << " residualConvergence " << std::scientific << residualConvergence << "\n";

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluid->setRebuildLinearPartInJacobian(false);M_fluid->setRebuildCstPartInLinearSystem(false);
            M_fluid->setRebuildCstPartInResidual(false);
            M_solid->setRebuildLinearPartInJacobian(false);M_solid->setRebuildCstPartInLinearSystem(false);
            M_solid->setRebuildCstPartInResidual(false);
        }

        this->updateBackendOptimisation(false,false);


        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("main","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("main","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");

        ++cptFSI;
    }

}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateTimeStep()
{
    M_previousTimeOrder=M_fluid->timeStepBDF()->timeOrder();

    M_fluid->updateTimeStep();

    if ( M_interpolationFSI )
    {
        if (M_solid->isStandardModel())
            M_interpolationFSI->updateFieldVelocitySolidPreviousPrevious( M_solid->timeStepNewmark()->previousVelocity() );
        else if ( M_solid->is1dReducedModel() )
            M_interpolationFSI->updateFieldVelocitySolid1dReducedPreviousPrevious( M_solid->timeStepNewmark1dReduced()->previousVelocity() );
    }
    M_solid->updateTimeStep();

    M_currentTimeOrder=M_fluid->timeStepBDF()->timeOrder();
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
boost::shared_ptr<std::ostringstream>
FSI<FluidType,SolidType>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << this->fluidModel()->getInfo()->str()
           << this->solidModel()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n||-----------------Info : FSI-------------------||"
           << "\n||==============================================||"
           << "\n   Coupling type : " << this->fsiCouplingType()
           << "\n   Coupling BC : " << this->fsiCouplingBoundaryCondition();
    if (M_interfaceFSIisConforme)
        *_ostr << "\n   Interface property : conformal";
    else
        *_ostr << "\n   Interface property : non conformal";

    *_ostr << "\n   Fix point parameters"
           << "\n     -- method : Aitken"
           << "\n     -- tolerance  : " << M_tolPtFixe
           << "\n     -- initial theta  : " << M_initialTheta
           << "\n     -- min theta  : " << M_minTheta
           << "\n     -- maxit : " << M_fixPointMaxIt;
    *_ostr << "\n||==============================================||"
           << "\n";

    return _ostr;
}


} // namespace FeelModels
} // namespace Feel
