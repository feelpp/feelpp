/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
FSI<FluidType,SolidType>::FSI(std::string const& prefix,WorldComm const& worldComm, std::string const& rootRepository )
    :
    super_type( prefix, worldComm, "", self_type::expandStringFromSpec( rootRepository ) ),
    M_meshSize( doption(_name="hsize",_prefix=this->prefix()) ),
    M_tagFileNameMeshGenerated( soption(_name="mesh-save.tag",_prefix=this->prefix()) ),
    M_fsiCouplingType( soption(_name="coupling-type",_prefix=this->prefix()) ),
    M_fsiCouplingBoundaryCondition( soption(_name="coupling-bc",_prefix=this->prefix()) ),
    M_interfaceFSIisConforme( boption(_name="conforming-interface",_prefix=this->prefix()) ),
    M_fixPointTolerance( doption(_name="fixpoint.tol",_prefix=this->prefix()) ),
    M_fixPointInitialTheta( doption(_name="fixpoint.initialtheta",_prefix=this->prefix()) ),
    M_fixPointMinTheta( doption(_name="fixpoint.min_theta",_prefix=this->prefix()) ),
    M_fixPointMaxIt( ioption(_name="fixpoint.maxit",_prefix=this->prefix()) ),
    M_fixPointMinItConvergence( ioption(_name="fixpoint.minit-convergence",_prefix=this->prefix()) ),
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
    M_reuseJacRebuildAtFirstFSIStepOptSolid( boption(_name="solid.reuse-jac.rebuild-at-first-fsi-step",_prefix=this->prefix()) ),
    M_reusePrecActivatedAfterNbFsiIterationFluid( ioption(_name="fluid.reuse-prec.activated-after-n-fsi-it",_prefix=this->prefix()) ),
    M_reusePrecActivatedAfterNbFsiIterationSolid( ioption(_name="solid.reuse-prec.activated-after-n-fsi-it",_prefix=this->prefix()) ),
    M_reusePrecActivatedToleranceFluid( doption(_name="fluid.reuse-prec.activated-only-if-greater-than-tol",_prefix=this->prefix()) ),
    M_reusePrecActivatedToleranceSolid( doption(_name="solid.reuse-prec.activated-only-if-greater-than-tol",_prefix=this->prefix()) ),
    M_couplingNitscheFamily_gamma( doption(_name="coupling-nitsche-family.gamma",_prefix=this->prefix()) ),
    M_couplingNitscheFamily_gamma0( doption(_name="coupling-nitsche-family.gamma0",_prefix=this->prefix()) ),
    M_couplingNitscheFamily_alpha( doption(_name="coupling-nitsche-family.alpha",_prefix=this->prefix()) ),
    M_couplingRNG_manualScaling( doption(_name="coupling-robin-neumann-generalized.manual-scaling",_prefix=this->prefix()) ),
    M_couplingRNG_useInterfaceOperator( boption(_name="coupling-robin-neumann-generalized.use-interface-operator",_prefix=this->prefix() ) )
{
    this->log("FSI","constructor","start");

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

    this->log("FSI","constructor","finish");
}

template <typename FluidType,typename SolidType>
std::string
FSI<FluidType,SolidType>::expandStringFromSpec( std::string const& expr )
{
    std::string res = expr;
    res = fluid_type::expandStringFromSpec( res );
    res = solid_type::expandStringFromSpec( res );
    return res;
}

//---------------------------------------------------------------------------------------------------------//

template <typename FluidType,typename SolidType>
void
FSI<FluidType,SolidType>::createMesh()
{
    this->log("FSI","createMesh","start");

    FSIMesh<typename fluid_type::convex_type> fsimeshTool(this->prefix(),this->worldComm());

    int nPart = this->worldComm().size();//this->nPartitions();

    fs::path meshesdirectories;
    if ( Environment::vm().count( prefixvm(this->prefix(),"mesh-save.directory" ) ) )
    {
        auto meshesdirgiven = fs::path(soption(_prefix=this->prefix(),_name="mesh-save.directory") );
        if ( meshesdirgiven.is_relative() )
            meshesdirectories = fs::path(Environment::rootRepository()) / meshesdirgiven;
        else
            meshesdirectories = meshesdirgiven;
    }
    else
        meshesdirectories = fs::path(this->repository().rootWithoutNumProc()) / fs::path("meshes");

    fs::path gp;
    if (this->hasMeshFile())
        gp = this->meshFile();
    else
        gp = this->geoFile();

    std::string nameMeshFile = gp.stem().string();
    fs::path mshFileNameFluidPart1, mshFileNameSolidPart1,mshFileNameFluidPartN,mshFileNameSolidPartN;
    fs::path mshFileNameFSI;
#ifdef FEELPP_HAS_HDF5
    std::string meshFileExtension = "json";
#else
    std::string meshFileExtension = "msh";
#endif
    if ( !M_tagFileNameMeshGenerated.empty() )
    {
        mshFileNameFluidPart1 = (boost::format("%1%_fluid_%2%_%3%_p1.%4%") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %meshFileExtension ).str();
        mshFileNameSolidPart1 = (boost::format("%1%_solid_%2%_%3%_p1.%4%") %nameMeshFile %solid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %meshFileExtension ).str();
        mshFileNameFluidPartN = (boost::format("%1%_fluid_%2%_%3%_p%4%.%5%") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %nPart %meshFileExtension ).str();
        mshFileNameSolidPartN = (boost::format("%1%_solid_%2%_%3%_p%4%.%5%") %nameMeshFile %solid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated %nPart %meshFileExtension ).str();
    }
    else
    {
        mshFileNameFluidPart1 = (boost::format("%1%_fluid_%2%_p1.%3%") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %meshFileExtension ).str();
        mshFileNameSolidPart1 = (boost::format("%1%_solid_%2%_p1.%3%") %nameMeshFile %solid_type::mesh_type::shape_type::name() %meshFileExtension ).str();
        mshFileNameFluidPartN = (boost::format("%1%_fluid_%2%_p%3%.%4%") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %nPart %meshFileExtension ).str();
        mshFileNameSolidPartN = (boost::format("%1%_solid_%2%_p%3%.%4%") %nameMeshFile %solid_type::mesh_type::shape_type::name() %nPart %meshFileExtension ).str();
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

    fsimeshTool.setForceRebuild( boption(_prefix=this->prefix(),_name="mesh-save.force-rebuild" ) );

    if (this->hasMeshFile())
    {
        fsimeshTool.setMshPathFSI( fs::path(this->meshFile()) );
        fsimeshTool.buildFSIMeshFromMsh();
    }
    else if (this->hasGeoFile())
    {
        fs::path mshFileNameFSI;
        if ( !M_tagFileNameMeshGenerated.empty() )
            mshFileNameFSI = (boost::format("%1%_fsi_%2%_%3%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() %M_tagFileNameMeshGenerated ).str();
        else
            mshFileNameFSI = (boost::format("%1%_fsi_%2%.msh") %nameMeshFile %fluid_type::mesh_type::shape_type::name() ).str();
        fs::path mshPathFSI = meshesdirectories / mshFileNameFSI;

        fsimeshTool.setGeoPathFSI( fs::path(this->geoFile()) );
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
    auto submeshStruct = createSubmesh( FM->meshALE()->referenceMesh(), markedfaces( FM->meshALE()->referenceMesh(), FM->markersNameMovingBoundary()/*"Paroi"*/) );
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
        submeshStruct->faceIterator( unwrap_ref(*itp).id() )->second.setMarker( submeshStruct->markerName("Fixe") );

    typedef SubMeshData smd_type;
    typedef boost::shared_ptr<smd_type> smd_ptrtype;
    smd_ptrtype smd( new smd_type(FM->mesh()) );
    for ( auto const& ew : elements(submeshStruct) )
    {
        auto const& e = unwrap_ref(ew);
        auto const& theface = FM->meshALE()->referenceMesh()->face( submeshStruct->subMeshToMesh(e.id()) );
        size_type idElt2 = FM->meshALE()->dofRelationShipMap()->geoElementMap().at( theface.element0().id() ).first;
        //std::cout << " e.G() " << e.G() << " other.G() " <<  theface.G() << std::endl;
        auto const& theface2 = FM->mesh()->element(idElt2).face(theface.pos_first());
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
        submeshStruct->faceIterator( unwrap_ref(*itp).id() )->second.setMarker( submeshStruct->markerName("Fixe") );

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
    this->log("FSI","init","start");

    // create fsimesh and partitioned meshes if require
    if ( this->hasMeshFile() || this->hasGeoFile() )
        this->createMesh();

    // fluid model build
    if ( !M_fluidModel )
    {
        M_fluidModel = fluid_ptrtype( new fluid_type("fluid",false,this->worldComm(), "", this->repository() ) );
        if ( !M_mshfilepathFluidPartN.empty() )
            M_fluidModel->setMeshFile(M_mshfilepathFluidPartN.string());
        //M_fluidModel->build();
        M_fluidModel->init();
    }

    // solid model build
    if ( !M_solidModel )
    {
        M_solidModel = solid_ptrtype( new solid_type("solid",false,this->worldComm(), "", this->repository() ) );
        bool doExtractSubmesh = boption(_name="solid-mesh.extract-1d-from-fluid-mesh",_prefix=this->prefix() );
        if ( doExtractSubmesh )
        {
            CHECK( !M_fluidModel->markersNameMovingBoundary().empty() ) << "no marker moving boundary in fluid model";

            if ( M_fluidModel->doRestart() )
                M_fluidModel->meshALE()->revertReferenceMesh();
            auto submeshStruct = detail::createMeshStruct1dFromFluidMesh<fluid_type,solid_type>( M_fluidModel );
            if ( M_fluidModel->doRestart() )
                M_fluidModel->meshALE()->revertMovingMesh();

            // TODO ( save 1d mesh and reload )
            if ( M_fluidModel->doRestart() )
                M_solidModel->build(submeshStruct);
            else
                M_solidModel->loadMesh(submeshStruct);

        }
        else
        {
            if ( !M_mshfilepathSolidPartN.empty() )
                M_solidModel->setMeshFile( M_mshfilepathSolidPartN.string() );
            M_solidModel->build();
        }
        M_solidModel->init();
    }

    CHECK( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
           this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
           this->fsiCouplingBoundaryCondition() == "nitsche" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" ) << "invalid fsiCouplingBoundaryCondition : " << this->fsiCouplingBoundaryCondition();
    M_fluidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());
    M_solidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        M_fluidModel->useFSISemiImplicitScheme(true);
        M_solidModel->useFSISemiImplicitScheme(true);
    }

    // specific value for robin
    M_solidModel->muFluidFSI( M_fluidModel->materialProperties()->cstMu() );
    M_fluidModel->setCouplingFSI_Nitsche_gamma( M_couplingNitscheFamily_gamma );
    M_solidModel->gammaNitschFSI( M_couplingNitscheFamily_gamma );
    M_fluidModel->setCouplingFSI_Nitsche_gamma0( M_couplingNitscheFamily_gamma0 );
    M_fluidModel->setCouplingFSI_Nitsche_alpha( M_couplingNitscheFamily_alpha );

    // save if reuse prec option at the begining
    M_reusePrecOptFluid = M_fluidModel->backend()->reusePrec();
    M_reuseJacOptFluid = M_fluidModel->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptFluid = M_fluidModel->backend()->reuseJacRebuildAtFirstNewtonStep();
    M_reusePrecOptSolid = M_solidModel->backend()->reusePrec();
    M_reuseJacOptSolid = M_solidModel->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptSolid = M_solidModel->backend()->reuseJacRebuildAtFirstNewtonStep();
    //-------------------------------------------------------------------------//
    // call init
    // M_fluidModel->init();
    // M_solidModel->init();
    this->updateTime( this->timeStepBase()->time() );
    //-------------------------------------------------------------------------//
    M_fluidModel->setCouplingFSI_solidIs1dReduced( M_solidModel->is1dReducedModel() );
    //-------------------------------------------------------------------------//
    // init interpolation tool
    M_interpolationFSI.reset(new interpolationFSI_type(M_fluidModel,M_solidModel));
    //-------------------------------------------------------------------------//
    // init aitken relaxation tool
    std::string aitkenType = "method1";
#if 0
    aitkenType = "fixed-relaxation";
    M_fixPointInitialTheta=0;
#endif
    M_fixPointConvergenceFSI.reset(new fixpointconvergenceFSI_type(M_solidModel) );
    M_aitkenFSI.reset(new aitkenrelaxationFSI_type(M_solidModel,
                                                   aitkenType,//AITKEN_METHOD_1,
                                                   M_fixPointInitialTheta,
                                                   M_fixPointTolerance,
                                                   M_fixPointMinTheta));

    //-------------------------------------------------------------------------//
    // build interface operator for generalized robin-neumann
    if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        auto fieldInit = M_fluidModel->meshVelocity2().functionSpace()->elementPtr();
        M_fluidModel->setCouplingFSI_RNG_evalForm1( fieldInit );
    }
    if ( M_solidModel->is1dReducedModel() )
        M_couplingRNG_useInterfaceOperator = false;
    if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        if ( !M_couplingRNG_useInterfaceOperator )
        {
            M_fluidModel->setCouplingFSI_RNG_useInterfaceOperator( false );
            if ( boption(_name="coupling-robin-neumann-generalized.without-interface-operator.precompute-mass-matrix",_prefix=this->prefix() ) )
                M_fluidModel->couplingFSI_RNG_updateForUse();
        }
        else if ( M_solidModel->isStandardModel() && M_couplingRNG_useInterfaceOperator )
    {
        auto fieldInitBis = M_fluidModel->functionSpaceVelocity()->elementPtr();
        M_fluidModel->setCouplingFSI_RNG_evalForm1Bis( fieldInitBis );

#if 1
        auto spaceDisp = M_solidModel->functionSpaceDisplacement();
        auto const& uDisp = M_solidModel->fieldDisplacement();
        auto mygraph = stencil(_test=spaceDisp,_trial=spaceDisp)->graph();

        //----------------------------------------------------------------------------------//
        // get dofs on markedfaces fsi
        std::set<size_type> dofMarkerFsi;
        for ( auto const& faceWrap : markedfaces(M_solidModel->mesh(),M_solidModel->markerNameFSI() ) )
        {
            //auto __face_it = faceMarked.template get<1>();
            //auto __face_en = faceMarked.template get<2>();
            //for( ; __face_it != __face_en; ++__face_it )
            //{
                auto const& face = boost::unwrap_ref( faceWrap );//*__face_it );
                for ( uint16_type l = 0; l < M_solidModel->functionSpaceDisplacement()->dof()->nLocalDofOnFace(true); ++l )
                {
                    for (uint16_type c1=0;c1<solid_type::nDim;c1++)
                    {
                        size_type gdof = boost::get<0>(M_solidModel->functionSpaceDisplacement()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                        dofMarkerFsi.insert( gdof );
                    }
                }
                //}
        }
#if 0 // bug
        std::set<size_type> dofMarkerFsiTest;
        for ( std::string const markName : M_solidModel->markerNameFSI() )
        {
            //functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 )
            auto setofdof = M_solidModel->functionSpaceDisplacement()->dof()->markerToDof( markName );
            for( auto it = setofdof.first, en = setofdof.second; it != en; ++ it )
                dofMarkerFsiTest.insert( it->second );
        }
        std::cout << "size " << dofMarkerFsiTest.size() << " vs " << dofMarkerFsi.size() << "\n";
#endif

        //----------------------------------------------------------------------------------//
        // mass matrix on solid
        auto matMass = M_solidModel->backend()->newMatrix(0,0,0,0,mygraph);
        form2( _trial=spaceDisp, _test=spaceDisp,_matrix=matMass )
            += integrate(
                //_range=markedfaces(M_solidModel->mesh(),M_solidModel->markerNameFSI()),
                _range=elements(M_solidModel->mesh()),
                //_quad=_Q<1>(),
                _expr=/*M_solidModel->mechanicalProperties()->cstRho()**/inner(idt(uDisp),id(uDisp)) );
        matMass->close();
        auto vecDiag = M_solidModel->backend()->newVector(spaceDisp);

        auto vecDiagMassLumped = M_solidModel->backend()->newVector(spaceDisp);
        auto unityVec = M_solidModel->backend()->newVector(spaceDisp);
        for (int k=0;k<unityVec->map().nLocalDofWithGhost();++k)
            unityVec->set(k,1);
        unityVec->close();
        matMass->multVector( unityVec,vecDiagMassLumped );
        vecDiagMassLumped->close();

        //----------------------------------------------------------------------------------//
        // mass matrix on interface
        auto matMassInterface = M_solidModel->backend()->newMatrix(0,0,0,0,mygraph);
        form2( _trial=spaceDisp, _test=spaceDisp,_matrix=matMassInterface )
            += integrate( _range=markedfaces(M_solidModel->mesh(),M_solidModel->markerNameFSI()),
                          //_quad=_Q<1>(),
                          _expr=inner(idt(uDisp),id(uDisp)) );
        matMassInterface->close();
        auto vecDiagMassInterface = M_solidModel->backend()->newVector(spaceDisp);
        M_solidModel->backend()->diag(matMassInterface,vecDiagMassInterface );
        vecDiagMassInterface->close();


        auto unityVecInterface = M_solidModel->backend()->newVector(spaceDisp);
        for (int k=0;k<unityVecInterface->map().nLocalDofWithGhost();++k)
            if ( dofMarkerFsi.find( k ) != dofMarkerFsi.end() ) 
                unityVecInterface->set(k,1);
        unityVecInterface->close();
        auto vecDiagMassInterfaceLumped = M_solidModel->backend()->newVector(spaceDisp);
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

            vecDiag->set( k, vecDiagMassLumped->operator()(k) );

        }
        vecDiag->close();
#endif


        M_interpolationFSI->setRobinNeumannInterfaceOperator( vecDiag );
        M_interpolationFSI->transfertRobinNeumannInterfaceOperatorS2F();
    }
    } // if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    //-------------------------------------------------------------------------//

    M_fluidModel->algebraicFactory()->addFunctionLinearAssembly( boost::bind( &self_type::updateLinearPDE_Fluid,
                                                                              boost::ref( *this ), _1 ) );
    M_fluidModel->algebraicFactory()->addFunctionLinearPostAssembly( boost::bind( &self_type::updateLinearPDEStrongDirichletBC_Fluid,
                                                                                  boost::ref( *this ), _1, _2 ) );

    this->log("FSI","init","finish");
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
FSI<FluidType,SolidType>::updateBackendOptimisation( int iterationFSI, double lastErrorRelative )
{
    bool isFirstFsiIteration = (iterationFSI == 0);

    bool canActivatedReusePrecFluid = ( (iterationFSI) >= M_reusePrecActivatedAfterNbFsiIterationFluid );
    bool canActivatedReusePrecSolid = ( (iterationFSI) >= M_reusePrecActivatedAfterNbFsiIterationSolid );
    if ( iterationFSI > 1 )
    {
        if ( lastErrorRelative > M_reusePrecActivatedToleranceSolid )
            canActivatedReusePrecSolid = false;
        if ( lastErrorRelative > M_reusePrecActivatedToleranceFluid )
            canActivatedReusePrecFluid = false;
    }

    // reuse-prec
    if ( M_reusePrecOptFluid && ( M_previousTimeOrder!=M_currentTimeOrder || M_reusePrecRebuildAtFirstFSIStepOptFluid ) )
    {
        bool doReusePrec = !isFirstFsiIteration && canActivatedReusePrecFluid;
        M_fluidModel->backend()->setReusePrec(doReusePrec);
    }
    if ( M_reusePrecOptSolid && M_reusePrecRebuildAtFirstFSIStepOptSolid )
    {
        bool doReusePrec = !isFirstFsiIteration && canActivatedReusePrecSolid;
        M_solidModel->backend()->setReusePrec(doReusePrec);
    }

    // reuse-jac.RebuildAtFirstNewtonStep
    if ( M_reuseJacOptFluid && (!M_reuseJacRebuildAtFirstNewtonStepOptFluid  && M_reuseJacRebuildAtFirstFSIStepOptFluid) )
    {
        M_fluidModel->backend()->setReuseJacRebuildAtFirstNewtonStep(isFirstFsiIteration);
    }
    if ( M_reuseJacOptSolid && ( !M_reuseJacRebuildAtFirstNewtonStepOptSolid && M_reuseJacRebuildAtFirstFSIStepOptSolid ) )
    {
        M_solidModel->backend()->setReuseJacRebuildAtFirstNewtonStep(isFirstFsiIteration);
    }


    this->log("FSI","updateBackendOptimisation",
              (boost::format("backend fluid reuse prec : %1%")%M_fluidModel->backend()->reusePrec() ).str());
    this->log("FSI","updateBackendOptimisation",
              (boost::format("backend solid reuse prec : %1%")%M_solidModel->backend()->reusePrec() ).str());
    if ( true )
        this->log("FSI","updateBackendOptimisation",
                  (boost::format("backend fluid reuse jac : %1%")%M_fluidModel->backend()->reuseJac() ).str());
    if ( true )
        this->log("FSI","updateBackendOptimisation",
                  (boost::format("backend solid reuse jac : %1%")%M_solidModel->backend()->reuseJac() ).str());
    if ( M_fluidModel->backend()->reuseJac() )
        this->log("FSI","updateBackendOptimisation",
                  (boost::format("backend fluid reuse jac (RebuildAtFirstNewtonStep) : %1%")%M_fluidModel->backend()->reuseJacRebuildAtFirstNewtonStep() ).str());
    if ( M_solidModel->backend()->reuseJac() )
        this->log("FSI","updateBackendOptimisation",
                  (boost::format("backend solid reuse jac (RebuildAtFirstNewtonStep) : %1%")%M_solidModel->backend()->reuseJacRebuildAtFirstNewtonStep() ).str());


}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl1()
{

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluidModel->updateALEmesh();

        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
    }

    // predictor disp
    M_solidModel->predictorDispl();
    M_solidModel->updateVelocity();

    // coupling fluid structure
    bool useAitken = true;//boption(_name="coupling-dirichlet-neumann.use-aitken",_prefix=this->prefix());
    if ( useAitken )
        this->aitkenRelaxTool()->restart();
    //aitkenFSI.setTheta(initialTheta);

    bool hasConverged = false;
    int cptFSI=0;
    double residualRelativeConvergence=1;
    //this->updateBackendOptimisation(true,true);
    boost::mpi::timer timerCur,timerIter;

    while ( (!this->aitkenRelaxTool()->isFinished() || cptFSI < this->fixPointMinItConvergence() ) &&
            cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        timerCur.restart();
        //--------------------------------------------------------------//
        this->aitkenRelaxTool()->saveOldSolution();
        this->updateBackendOptimisation(cptFSI,residualRelativeConvergence);
        //--------------------------------------------------------------//
        // ALE solver
        timerCur.restart();
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluidModel->updateALEmesh();
            double tALE = timerCur.elapsed();
            this->log("FSI","update ale","finish in "+(boost::format("%1% s") % tALE).str() );
        }
        else  if (this->fsiCouplingType()=="Semi-Implicit")
        {
            this->interpolationTool()->transfertVelocity();
            double tALE = timerCur.elapsed();
            this->log("FSI","transfert velocity","finish in "+(boost::format("%1% s") % tALE).str());
        }
        //--------------------------------------------------------------//
        M_fluidModel->solve();
        //--------------------------------------------------------------//
        timerCur.restart();
        // revert ref mesh
        //M_fluidModel->meshALE()->revertReferenceMesh();
        // transfert stress
        this->interpolationTool()->transfertStress();
        // revert moving mesh
        //M_fluidModel->meshALE()->revertMovingMesh();
        double t3 = timerCur.elapsed();
        this->log("FSI","transfert stress","finish in "+(boost::format("%1% s") % t3).str() );
        //--------------------------------------------------------------//
        M_solidModel->solve();
        //--------------------------------------------------------------//
        timerCur.restart();
        //compute and apply aitken relaxation
        this->aitkenRelaxTool()->applyRelaxation();
        // update velocity and acceleration
        M_solidModel->updateVelocity();
        // aitken relaxtion
        if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
            this->aitkenRelaxTool()->printInfo();
        this->aitkenRelaxTool()->shiftRight();
        hasConverged = this->aitkenRelaxTool()->isFinished();
        residualRelativeConvergence = this->aitkenRelaxTool()->residualNorm();

        double t4 = timerCur.elapsed();
        this->log("FSI","apply relax and up vel/acc struct","finish in "+(boost::format("%1% s") % t4).str() );
        //--------------------------------------------------------------//

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluidModel->setRebuildLinearPartInJacobian(false);M_fluidModel->setRebuildCstPartInLinearSystem(false);
            M_solidModel->setRebuildLinearPartInJacobian(false);M_solidModel->setRebuildCstPartInLinearSystem(false);
            M_fluidModel->setRebuildCstPartInResidual(false);M_solidModel->setRebuildCstPartInResidual(false);
        }
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("FSI","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        ++cptFSI;
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
        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
    }
    bool useAitken = boption(_name="coupling-nitsche-family.use-aitken",_prefix=this->prefix());
    if ( useAitken )
        this->aitkenRelaxTool()->restart();

    bool solveStruct = this->fluidModel()->timeStepBDF()->iteration() > 1 || this->fixPointMaxIt()==1;

    bool hasConverged = false;
    int cptFSI=0;
    double residualRelativeConvergence=1;
    while ( ( !hasConverged || cptFSI < this->fixPointMinItConvergence() ) && cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        if ( useAitken )
            this->aitkenRelaxTool()->saveOldSolution();
        else
            M_fixPointConvergenceFSI->saveOldSolution();

        this->updateBackendOptimisation(cptFSI,residualRelativeConvergence);

        //--------------------------------------------------------------//
        if (solveStruct)
        {
            //M_fluidModel->meshALE()->revertReferenceMesh();
            // transfert stress
            this->interpolationTool()->transfertStress();
            // revert moving mesh
            //M_fluidModel->meshALE()->revertMovingMesh();
            if ( ( this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" ||
                   this->fsiCouplingBoundaryCondition()=="nitsche" ) &&
                 M_solidModel->isStandardModel() )
            {
                bool useExtrap = boption(_prefix=this->prefix(),_name="transfert-velocity-F2S.use-extrapolation");
                this->interpolationTool()->transfertVelocityF2S(cptFSI,useExtrap);
            }
            M_solidModel->solve();
            M_solidModel->updateVelocity();

            if ( useAitken )
            {
                this->aitkenRelaxTool()->applyRelaxation();
                // update velocity and acceleration
                M_solidModel->updateVelocity();
                // aitken relaxtion
                if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
                    this->aitkenRelaxTool()->printInfo();
                this->aitkenRelaxTool()->shiftRight();
                hasConverged = this->aitkenRelaxTool()->isFinished();
                residualRelativeConvergence = this->aitkenRelaxTool()->residualNorm();
            }

        }

        //--------------------------------------------------------------//
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluidModel->updateALEmesh();
        }
        else
        {
            this->interpolationTool()->transfertVelocity();
        }
        M_fluidModel->solve();
        //--------------------------------------------------------------//

        if ( !useAitken )
        {
            residualRelativeConvergence = M_fixPointConvergenceFSI->computeConvergence();
            hasConverged = residualRelativeConvergence < this->fixPointTolerance();
            if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
                std::cout << "["<<prefixvm(this->prefix(),"FSI") <<"] iteration " << cptFSI
                          << " residualRelativeConvergence " << std::scientific << residualRelativeConvergence << "\n";
        }

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluidModel->setRebuildLinearPartInJacobian(false);M_fluidModel->setRebuildCstPartInLinearSystem(false);
            M_fluidModel->setRebuildCstPartInResidual(false);

            if (solveStruct)
            {
                M_solidModel->setRebuildLinearPartInJacobian(false);M_solidModel->setRebuildCstPartInLinearSystem(false);
                M_solidModel->setRebuildCstPartInResidual(false);
            }
        }

        if (!solveStruct) solveStruct=true;

        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("FSI","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");

        ++cptFSI;
    }

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->interpolationTool()->transfertDisplacement();
        M_fluidModel->updateALEmesh();
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
        M_fluidModel->updateALEmesh();

        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
    }

    double manualScalingRNG = M_couplingRNG_manualScaling;

    bool useAitken = boption(_name="coupling-robin-neumann-generalized.use-aitken",_prefix=this->prefix());
    if ( useAitken )
        this->aitkenRelaxTool()->restart();

    bool hasConverged = false;
    int cptFSI=0;
    double residualRelativeConvergence=1;
    while ( ( !hasConverged || cptFSI <this->fixPointMinItConvergence() ) && cptFSI < this->fixPointMaxIt() )
    {
        timerIter.restart();
        //timerCur.restart();
        this->updateBackendOptimisation(cptFSI,residualRelativeConvergence);

        if ( useAitken )
            this->aitkenRelaxTool()->saveOldSolution();
        else
            M_fixPointConvergenceFSI->saveOldSolution();


        //--------------------------------------------------------------//
        if (this->fsiCouplingType()=="Implicit")
        {
            this->interpolationTool()->transfertDisplacement();
            M_fluidModel->updateALEmesh();
        }
        else
        {
            //bool useExtrap = (this->fluidModel()->timeStepBDF()->iteration() > 1) && cptFSI==0;
            //this->interpolationTool()->transfertVelocity(useExtrap);
        }
        this->interpolationTool()->transfertRobinNeumannGeneralizedS2F( cptFSI, manualScalingRNG );
        M_fluidModel->solve();

        //--------------------------------------------------------------//
        //M_fluidModel->meshALE()->revertReferenceMesh();
        // transfert stress
        this->interpolationTool()->transfertStress();
        // revert moving mesh
        //M_fluidModel->meshALE()->revertMovingMesh();
        M_solidModel->solve();
        M_solidModel->updateVelocity();
        //--------------------------------------------------------------//

        if ( !useAitken )
        {
            residualRelativeConvergence = M_fixPointConvergenceFSI->computeConvergence();
            hasConverged = residualRelativeConvergence < this->fixPointTolerance();
            if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
                std::cout << "["<<prefixvm(this->prefix(),"FSI") <<"] iteration " << cptFSI
                          << " residualRelativeConvergence " << std::scientific << residualRelativeConvergence << "\n";
        }
        else
        {
            this->aitkenRelaxTool()->applyRelaxation();
            // update velocity and acceleration
            M_solidModel->updateVelocity();
            // aitken relaxtion
            if (this->worldComm().isMasterRank() && this->verboseSolverTimer())
                this->aitkenRelaxTool()->printInfo();
            this->aitkenRelaxTool()->shiftRight();
            hasConverged = this->aitkenRelaxTool()->isFinished();
            residualRelativeConvergence = this->aitkenRelaxTool()->residualNorm();
        }

        if (this->fsiCouplingType()=="Semi-Implicit")
        {
            M_fluidModel->setRebuildLinearPartInJacobian(false);M_fluidModel->setRebuildCstPartInLinearSystem(false);
            M_fluidModel->setRebuildCstPartInResidual(false);
            M_solidModel->setRebuildLinearPartInJacobian(false);M_solidModel->setRebuildCstPartInLinearSystem(false);
            M_solidModel->setRebuildCstPartInResidual(false);
        }
        //--------------------------------------------------------------//
        //--------------------------------------------------------------//
        double timeElapsedIter = timerIter.elapsed();
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");
        this->log("FSI","iteration fsi","finish in "+(boost::format("%1% s") % timeElapsedIter).str());
        this->log("FSI","iteration fsi","--------------------xxxxxxxxxxxxxxxxxxxx--------------------");

        ++cptFSI;
    }

}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateTime(double time)
{
    M_fluidModel->updateTime(time);
    M_solidModel->updateTime(time);
    super_type::updateTime(time);
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateTimeStep()
{
    M_previousTimeOrder=M_fluidModel->timeStepBDF()->timeOrder();

    M_fluidModel->updateTimeStep();
    M_solidModel->updateTimeStep();

    M_currentTimeOrder=M_fluidModel->timeStepBDF()->timeOrder();

    this->updateTime( this->timeStepBase()->time() );
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

    bool useAitken = false;
    if ( this->fsiCouplingBoundaryCondition()  == "dirichlet-neumann" )
    {
        useAitken = true;
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
              this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
              this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        *_ostr << "\n   Nitsche family parameters"
               << "\n     -- gamma : " << M_couplingNitscheFamily_gamma;
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
              this->fsiCouplingBoundaryCondition() == "nitsche" )
            *_ostr << "\n     -- gamma0 : " << M_couplingNitscheFamily_gamma0;
        *_ostr << "\n     -- alpha : " << M_couplingNitscheFamily_alpha;

    }
    else if (this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        *_ostr << "\n   Generalized Robin-Neumann"
               << "\n     -- use interface operator : " << std::boolalpha << M_couplingRNG_useInterfaceOperator
               << "\n     -- manual scaling : " << M_couplingRNG_manualScaling;
    }

    *_ostr << "\n   Solver";
    if ( useAitken )
        *_ostr << "\n     -- method        : fix point with Aitken relaxation";
    else
        *_ostr << "\n     -- method        : fix point";

    *_ostr << "\n     -- tolerance     : " << M_fixPointTolerance
           << "\n     -- maxit         : " << M_fixPointMaxIt;

    if ( useAitken )
        *_ostr << "\n     -- initial theta : " << M_fixPointInitialTheta
               << "\n     -- min theta     : " << M_fixPointMinTheta;

    *_ostr << "\n||==============================================||"
           << "\n";

    return _ostr;
}


} // namespace FeelModels
} // namespace Feel
