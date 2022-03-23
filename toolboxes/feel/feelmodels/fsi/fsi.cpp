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
#include <feel/feelpde/operatorpcd.hpp>

namespace Feel
{
namespace FeelModels
{

template< class FluidType, class SolidType >
FSI<FluidType,SolidType>::FSI( std::string const& prefix,
                               std::string const& keyword,
                               worldcomm_ptr_t const& worldComm,
                               ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, "", modelRep ),
    ModelPhysics<mesh_fluid_type::nRealDim>( "fsi" ),
    ModelBase( prefix, keyword, worldComm, "", modelRep ),
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
    M_coulingRNG_strategyTimeStepCompatibility( soption(_name="coupling-robin-neumann-generalized.strategy-time-step-compatibility",_prefix=this->prefix()) )
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

//---------------------------------------------------------------------------------------------------------//

template <typename FluidType,typename SolidType>
void
FSI<FluidType,SolidType>::createMesh()
{
    this->log("FSI","createMesh","start");
    FSIMesh<typename fluid_type::convex_type> fsimeshTool(this->prefix(),this->worldCommPtr());

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
#if 0
    fs::path gp;
    if (this->hasMeshFile())
        gp = this->meshFile();
    else
        gp = this->geoFile();

    std::string nameMeshFile = gp.stem().string();
#endif
    std::string nameMeshFile = "fsi";
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


    this->modelMesh( this->keyword() ).importConfig().setStraightenMesh( false );
    this->modelMesh( this->keyword() ).importConfig().setupSequentialAndLoadByMasterRankOnly();
    this->modelMesh( this->keyword() ).importConfig().setMeshComponents( MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES );
    if ( !fs::exists( M_mshfilepathFluidPart1 ) || !fs::exists( M_mshfilepathSolidPart1 ) || fsimeshTool.forceRebuild() )
    {
        super_type::super_model_meshes_type::updateForUse<mesh_fluid_type>( this->keyword() );
        if ( this->worldComm().isMasterRank() )
            fsimeshTool.buildSubMesh( this->modelMesh( this->keyword() ).template mesh<mesh_fluid_type>()  );
    }
    if ( nPart > 1 )
        fsimeshTool.buildMeshesPartitioning();

    this->log("FSI","createMesh","finish");
}

//---------------------------------------------------------------------------------------------------------//

namespace detail
{

template <typename FluidType,typename SolidType>
typename SolidType::solid_1dreduced_type::mesh_ptrtype
createMeshStruct1dFromFluidMesh2d( typename FluidType::self_ptrtype const& FM, std::set<std::string> const& markersFSI, mpl::bool_<false> /**/ )
{
    auto submeshStruct = createSubmesh( _mesh=FM->meshMotionTool()->referenceMesh(), _range=markedfaces( FM->meshMotionTool()->referenceMesh(), markersFSI ) );
#if 0
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
        submeshStruct->faceIterator( unwrap_ref(*itp).id() )->second.setMarker( submeshStruct->markerName("Fixe") );
#endif
    typedef SubMeshData<typename FluidType::mesh_type::index_type> smd_type;
    typedef std::shared_ptr<smd_type> smd_ptrtype;
    smd_ptrtype smd( new smd_type(FM->mesh()) );
    for ( auto const& ew : elements(submeshStruct) )
    {
        auto const& e = unwrap_ref(ew);
        auto const& theface = FM->meshMotionTool()->referenceMesh()->face( submeshStruct->subMeshToMesh(e.id()) );
        size_type idElt2 = FM->meshMotionTool()->dofRelationShipMap()->geoElementMap().at( theface.element0().id() ).first;
        //std::cout << " e.G() " << e.G() << " other.G() " <<  theface.G() << std::endl;
        auto const& theface2 = FM->mesh()->element(idElt2).face(theface.pos_first());
        //smd->bm.insert( typename smd_type::bm_type::value_type( e.id(), theface2.id() ) );
        smd->bm.insert( { e.id(), theface2.id() } );
    }
    submeshStruct->setSubMeshData( smd );

    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename FluidType::mesh_type::trace_mesh_ptrtype
createMeshStruct1dFromFluidMesh2d( typename FluidType::self_ptrtype const& FM, std::set<std::string> const& markersFSI, mpl::bool_<true> /**/ )
{
    auto submeshStruct = createSubmesh( _mesh=FM->mesh(), _range=markedfaces( FM->mesh(),markersFSI ) );
#if 0
    auto hola = boundaryfaces(submeshStruct);
    for ( auto itp = hola.template get<1>(),enp = hola.template get<2>() ; itp!=enp ; ++itp )
        submeshStruct->faceIterator( unwrap_ref(*itp).id() )->second.setMarker( submeshStruct->markerName("Fixe") );
#endif
    return submeshStruct;
}

template <typename FluidType,typename SolidType>
typename SolidType::solid_1dreduced_type::mesh_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, std::set<std::string> const& markersFSI, mpl::int_<2> /**/ )
{
    static const bool hasSameOrderGeo = FluidType::mesh_type::nOrder == SolidType::solid_1dreduced_type::mesh_type::nOrder;
    return createMeshStruct1dFromFluidMesh2d<FluidType,SolidType>(FM, markersFSI, mpl::bool_<hasSameOrderGeo>() );
}

template <typename FluidType,typename SolidType>
typename SolidType::solid_1dreduced_type::mesh_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, std::set<std::string> const& markersFSI, mpl::int_<3> /**/ )
{
    CHECK( false ) << "not possible";
    return typename SolidType::solid_1dreduced_type::mesh_ptrtype();
}

template <typename FluidType,typename SolidType>
typename SolidType::solid_1dreduced_type::mesh_ptrtype
createMeshStruct1dFromFluidMesh( typename FluidType::self_ptrtype const& FM, std::set<std::string> const& markersFSI )
{
    return createMeshStruct1dFromFluidMesh<FluidType,SolidType>( FM, markersFSI, mpl::int_<SolidType::nDim>() );
}

} // namespace detail



template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updatePhysics( typename super_physics_type::PhysicsTreeNode & physicsTree, ModelModels const& models )
{
    if ( !M_fluidModel )
        M_fluidModel = std::make_shared<fluid_type>("fluid","fluid",this->worldCommPtr(), "", this->repository() );
    if ( !M_solidModel )
        M_solidModel = std::make_shared<solid_type>("solid","solid",this->worldCommPtr(), "", this->repository() );

    physicsTree.addChild( M_fluidModel, models );
    physicsTree.addChild( M_solidModel, models );

    physicsTree.updateMaterialSupportFromChildren( "merge" );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::init()
{
    this->log("FSI","init","start");

    this->initModelProperties();

    this->initPhysics( this->shared_from_this(), this->modelProperties().models() );

    // physical properties
    if ( !M_materialsProperties )
    {
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    if ( this->modelProperties().jsonData().contains("Meshes") )
        super_type::super_model_meshes_type::setup( this->modelProperties().jsonData().at("Meshes"), {this->keyword()} );

    // create fsimesh and partitioned meshes if require
    if ( !this->modelMesh( this->keyword() ).importConfig().inputFilename().empty() && !this->doRestart() )
        this->createMesh();

    std::set<std::string> markersFSI_fluid;// = { "fsiWall" /*"fsi-wall"*/ }; // this->fluidModel()->markersFSI()
    std::set<std::string> markersFSI_solid;// = { "fsiWall" /*"fsi-wall"*/ }; // this->solidModel()->markerNameFSI()

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFSIData = std::static_pointer_cast<ModelPhysicFSI<mesh_fluid_type::nRealDim>>(physicData);
        markersFSI_fluid.insert( physicFSIData->interfaceFluid().begin(),physicFSIData->interfaceFluid().end() );
        markersFSI_solid.insert( physicFSIData->interfaceSolid().begin(),physicFSIData->interfaceSolid().end() );
    }
    // if ( this->worldComm().isMasterRank() )
    // {
    //     std::cout << "markersFSI_fluid :  " << markersFSI_fluid << std::endl;
    //     std::cout << "markersFSI_solid :  " << markersFSI_solid << std::endl;
    // }


    // fluid model build
    //if ( !M_fluidModel )
    {

        // if ( this->hasModelMesh( M_fluidModel->keyword() ) )
        //     M_heatModel->setModelMeshAsShared( this->modelMesh() );
        //M_fluidModel = std::make_shared<fluid_type>("fluid","fluid",this->worldCommPtr(), "", this->repository() );
        if ( !M_mshfilepathFluidPartN.empty() )
            M_fluidModel->modelMesh( M_fluidModel->keyword() ).importConfig().setupInputMeshFilenameWithoutApplyPartitioning( M_mshfilepathFluidPartN.string() );

        M_fluidModel->setManageParameterValues( false );
        if ( !M_fluidModel->modelPropertiesPtr() )
        {
            M_fluidModel->setModelProperties( this->modelPropertiesPtr() );
            M_fluidModel->setManageParameterValuesOfModelProperties( false );
        }
        M_fluidModel->setMaterialsProperties( M_materialsProperties );

        // temporary fix (else use in dirichle-neunamm bc in residual) TODO !!!!
        M_fluidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());
        M_fluidModel->init( false );
    }

    this->log("FSI","init","fluid init done");


    // up mesh motion tool
    this->fluidModel()->meshMotionTool()->addMarkersInBoundaryCondition( "moving", markersFSI_fluid );
     //this->fluidModel()->meshMotionTool()->setDisplacementImposedOnInitialDomainOverFaces( M_fluidModel->keyword()+"_"+this->keyword(), markersFSI_fluid );

    // revert fluid reference mesh if restart
    if ( this->fluidModel()->doRestart() )
    {
        this->fluidModel()->meshMotionTool()->revertReferenceMesh();
        // need to rebuild this dof point because updated in meshale after restart
        this->fluidModel()->meshMotionTool()->displacement()->functionSpace()->rebuildDofPoints();
        //this->fluidModel()->functionSpaceVelocity()->rebuildDofPoints();
    }



    // solid model build
    //if ( !M_solidModel )
    {
        //M_solidModel = std::make_shared<solid_type>("solid","solid",this->worldCommPtr(), "", this->repository() );
        bool doExtractSubmesh = boption(_name="solid-mesh.extract-1d-from-fluid-mesh",_prefix=this->prefix() );
        if ( doExtractSubmesh )
        {
            //CHECK( !M_fluidModel->markersFSI().empty() ) << "no marker moving boundary in fluid model";

            //if ( M_fluidModel->doRestart() )
            //M_fluidModel->meshMotionTool()->revertReferenceMesh();
            auto submeshStruct = detail::createMeshStruct1dFromFluidMesh<fluid_type,solid_type>( M_fluidModel, markersFSI_fluid );
            //if ( M_fluidModel->doRestart() )
            //M_fluidModel->meshMotionTool()->revertMovingMesh();

            M_solidModel->createsSolid1dReduced();
            // TODO ( save 1d mesh and reload )
            if ( M_fluidModel->doRestart() )
                M_solidModel->solid1dReduced()->setMesh(submeshStruct);
            else
                M_solidModel->solid1dReduced()->setMesh(submeshStruct);

        }
        else
        {
            if ( !M_mshfilepathSolidPartN.empty() )
                M_solidModel->modelMesh( M_solidModel->keyword() ).importConfig().setupInputMeshFilenameWithoutApplyPartitioning( M_mshfilepathSolidPartN.string() );
        }


        M_solidModel->setManageParameterValues( false );
        if ( !M_solidModel->modelPropertiesPtr() )
        {
            M_solidModel->setModelProperties( this->modelPropertiesPtr() );
            M_solidModel->setManageParameterValuesOfModelProperties( false );
        }
        M_solidModel->setMaterialsProperties( M_materialsProperties );

        // temporary fix TODO !!!!
        M_solidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());
        if (this->fsiCouplingType()=="Semi-Implicit")
            M_solidModel->useFSISemiImplicitScheme(true);

        M_solidModel->init();
    }

    this->log("FSI","init","solid init done");

    CHECK( this->fluidModel()->hasMeshMotion() ) << "aiaie";

    CHECK( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
           this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
           this->fsiCouplingBoundaryCondition() == "nitsche" ||
           this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" ) << "invalid fsiCouplingBoundaryCondition : " << this->fsiCouplingBoundaryCondition();
    M_fluidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());
    M_solidModel->couplingFSIcondition(this->fsiCouplingBoundaryCondition());

    M_fluidModel->setApplyMovingMeshBeforeSolve( false );

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        M_fluidModel->useFSISemiImplicitScheme(true);
        M_solidModel->useFSISemiImplicitScheme(true);
    }


    M_rangeFSI_fluid = markedfaces( this->fluidModel()->mesh(),markersFSI_fluid );
    auto submeshfsi_fluid = createSubmesh( _mesh=this->fluidModel()->mesh(),_range=M_rangeFSI_fluid,_view=1 );
    M_spaceNormalStress_fluid = fluid_type::space_normalstress_type::New(_mesh=submeshfsi_fluid );
    M_fieldNormalStressRefMesh_fluid.reset( new typename fluid_type::element_normalstress_type( M_spaceNormalStress_fluid ) );

    M_XhMeshVelocityInterface = space_fluid_meshvelocityonboundary_type::New(_mesh=submeshfsi_fluid );
    // mesh velocity only on moving interface
    M_meshVelocityInterface.reset(new element_fluid_meshvelocityonboundary_type( M_XhMeshVelocityInterface ) );

    M_meshDisplacementOnInterface_fluid = this->fluidModel()->meshMotionTool()->displacement()->functionSpace()->elementPtr();

    this->log("FSI","init","fsi functionspace done");

    this->fluidModel()->updateRangeDistributionByMaterialName( "interface_fsi", M_rangeFSI_fluid );

    M_fluidModel->addUpdateInHousePreconditionerPCD( "fsi",
                                                     std::bind( &self_type::initInHousePreconditionerPCD_fluid, std::ref( *this ), std::placeholders::_1 ),
                                                     std::bind( &self_type::updateInHousePreconditionerPCD_fluid, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 ) );
    M_fluidModel->initAlgebraicFactory();
    this->log("FSI","init","fluid initAlgebraicFactory done");


    if ( M_solidModel->isStandardModel() )
    {
        M_rangeFSI_solid = markedfaces(this->solidModel()->mesh(),markersFSI_solid);

        //M_spaceNormalStressFromFluid_solid = space_solid_normalstressfromfluid_type::New( _mesh=M_solidModel->mesh() );
        auto subfsimesh = createSubmesh( _mesh=this->solidModel()->mesh(),_range=M_rangeFSI_solid, _view=true ) ;
        M_spaceNormalStressFromFluid_solid = space_solid_normalstressfromfluid_type::New( _mesh=subfsimesh );
        M_fieldNormalStressFromFluid_solid.reset(new element_solid_normalstressfromfluid_type( M_spaceNormalStressFromFluid_solid ) );
        if ( this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
             this->fsiCouplingBoundaryCondition() == "robin-robin" ||
             this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
             this->fsiCouplingBoundaryCondition() == "nitsche" )
            M_fieldVelocityInterfaceFromFluid_solid.reset( new typename solid_type::element_vectorial_type( M_solidModel->functionSpaceDisplacement() ) );
    }
    else if ( M_solidModel->is1dReducedModel() )
    {
        // normal stress as source term
        M_spaceNormalStressFromFluid_solid1dReduced = space_solid1dreduced_normalstressfromfluid_vect_type::New(_mesh=M_solidModel->solid1dReduced()->mesh() );
        M_fieldNormalStressFromFluidScalar_solid1dReduced.reset( new element_solid1dreduced_normalstressfromfluid_scal_type( M_spaceNormalStressFromFluid_solid1dReduced->compSpace() ) );
        M_fieldNormalStressFromFluidVectorial_solid1dReduced.reset( new element_solid1dreduced_normalstressfromfluid_vect_type( M_spaceNormalStressFromFluid_solid1dReduced ) );
    }

    this->log("FSI","init","fsi functionspace [solid] done");

    // save if reuse prec option at the begining
    M_reusePrecOptFluid = M_fluidModel->backend()->reusePrec();
    M_reuseJacOptFluid = M_fluidModel->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptFluid = M_fluidModel->backend()->reuseJacRebuildAtFirstNewtonStep();
    M_reusePrecOptSolid = M_solidModel->backend()->reusePrec();
    M_reuseJacOptSolid = M_solidModel->backend()->reuseJac();
    M_reuseJacRebuildAtFirstNewtonStepOptSolid = M_solidModel->backend()->reuseJacRebuildAtFirstNewtonStep();
    //-------------------------------------------------------------------------//
    this->updateTime( this->timeStepBase()->time() );
    //-------------------------------------------------------------------------//
    // init interpolation tool
    this->initInterpolation();
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
    auto XhFluidVelocity = this->fluidModel()->functionSpaceVelocity();
    M_dofsMultiProcessVelocitySpaceOnFSI_fluid = XhFluidVelocity->dofs( M_rangeFSI_fluid, ComponentType::NO_COMPONENT, true );
    if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        this->initCouplingRobinNeumannGeneralized();
    }
    else if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        auto dofsToAdd = XhFluidVelocity->dofs(  M_rangeFSI_fluid );
        XhFluidVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("fluid.velocity",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        //auto dofsMultiProcessToAdd = XhFluidVelocity->dofs( M_rangeFSI_fluid, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("fluid.velocity",MESH_FACES).insert( M_dofsMultiProcessVelocitySpaceOnFSI_fluid/*dofsMultiProcessToAdd*/.begin(), M_dofsMultiProcessVelocitySpaceOnFSI_fluid/*dofsMultiProcessToAdd*/.end() );
    }

    if ( ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
           this->fsiCouplingBoundaryCondition() == "nitsche" ) && this->solidModel()->isStandardModel() )
    {
        M_fieldsGradVelocity_fluid.resize( fluid_type::nDim );
        for (int k = 0;k<fluid_type::nDim ;++k)
            M_fieldsGradVelocity_fluid[k] = M_spaceNormalStress_fluid->elementPtr();
        M_fieldsGradVelocity_solid.resize( fluid_type::nDim );
        for (int k = 0;k<fluid_type::nDim ;++k)
            M_fieldsGradVelocity_solid[k] = M_spaceNormalStressFromFluid_solid->elementPtr();
    }
    //-------------------------------------------------------------------------//

    if ( this->fluidModel()->doRestart() )
    {
        this->fluidModel()->meshMotionTool()->revertMovingMesh();
        this->fluidModel()->meshMotionTool()->displacement()->functionSpace()->rebuildDofPoints();
    }

    //-------------------------------------------------------------------------//

    M_fluidModel->algebraicFactory()->addFunctionLinearAssembly( std::bind( &self_type::updateLinearPDE_Fluid,
                                                                            std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionJacobianAssembly( std::bind( &self_type::updateJacobian_Fluid,
                                                                              std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionResidualAssembly( std::bind( &self_type::updateResidual_Fluid,
                                                                              std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionLinearDofElimination( std::bind( &self_type::updateLinearPDEDofElimination_Fluid,
                                                                                  std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionNewtonInitialGuess( std::bind( &self_type::updateNewtonInitialGuess_Fluid,
                                                                                std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionJacobianDofElimination( std::bind( &self_type::updateJacobianDofElimination_Fluid,
                                                                                    std::ref( *this ), std::placeholders::_1 ) );
    M_fluidModel->algebraicFactory()->addFunctionResidualDofElimination( std::bind( &self_type::updateResidualDofElimination_Fluid,
                                                                                    std::ref( *this ), std::placeholders::_1 ) );
    if ( M_solidModel->isStandardModel() )
    {
        M_solidModel->algebraicFactory()->addFunctionLinearAssembly( std::bind( &self_type::updateLinearPDE_Solid,
                                                                                std::ref( *this ), std::placeholders::_1 ) );
        M_solidModel->algebraicFactory()->addFunctionJacobianAssembly( std::bind( &self_type::updateJacobian_Solid,
                                                                                  std::ref( *this ), std::placeholders::_1 ) );
        M_solidModel->algebraicFactory()->addFunctionResidualAssembly( std::bind( &self_type::updateResidual_Solid,
                                                                                  std::ref( *this ), std::placeholders::_1 ) );
    }
    else if ( M_solidModel->is1dReducedModel() )
    {
        M_solidModel->algebraicFactory()->addFunctionLinearAssembly( std::bind( &self_type::updateLinearPDE_Solid1dReduced,
                                                                                std::ref( *this ), std::placeholders::_1 ) );
    }
    this->log("FSI","init","finish");
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initCouplingRobinNeumannGeneralized()
{
    bool useAlgebraicInnerProductWithLumping = false;

    M_couplingRNG_evalForm1 = M_fluidModel->functionSpaceVelocity()/*meshVelocity2().functionSpace()*/->elementPtr();
    if ( M_fluidModel->doRestart() )
    {
        M_fluidModel->updateNormalStressOnReferenceMesh( "interface_fsi", M_fieldNormalStressRefMesh_fluid );
    }

    if ( this->solidModel()->isStandardModel() )
    {
        auto Vh = this->solidModel()->functionSpaceDisplacement();
        auto mesh = Vh->mesh();

        if ( !this->solidModel()->massMatrixLumped() )
            this->solidModel()->updateMassMatrixLumped();
        auto massMatrixLumped = this->solidModel()->massMatrixLumped();

        this->solidModel()->setUseMassMatrixLumped( boption(_name="coupling-robin-neumann-generalized.use-mass-matrix-lumped-in-solid",_prefix=this->prefix()) );

        M_coulingRNG_operatorDiagonalOnFluid = this->fluidModel()->functionSpaceVelocity()/*meshVelocity2().functionSpace()*/->elementPtr();

        //--------------------------------------------------------
        auto rangeFSI = M_rangeFSI_solid;
        double areaFSI = measure(_range=rangeFSI);

        std::set<size_type> dofsIdOnFSI;
        for ( auto const& faceWrap : rangeFSI )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = Vh->dof()->faceLocalDof( face.id() );
            for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                dofsIdOnFSI.insert( it->index() );
        }

        auto vecDiagMassLumped2 = this->solidModel()->backend()->newVector(Vh);
        massMatrixLumped->diagonal(vecDiagMassLumped2);

        if ( Environment::vm().count( prefixvm(this->prefix(),"coupling-robin-neumann-generalized.use-operator-constant" ) ) )
        {
            // consider the operator B as c*Id and c is the constant given
            M_coulingRNG_operatorDiagonalOnFluid->setConstant( doption(_name="coupling-robin-neumann-generalized.use-operator-constant", _prefix=this->prefix()) );
        }
        else if ( true )
        {
            // compute algebraic counterpart of the inner_product on Sigma with include operator B ( this is the version of the paper)
            auto vecDiagMassMatrixLumpedUBLAS = Vh->element( vecDiagMassLumped2 );
            auto qqqq = Vh->element();
            qqqq = vecDiagMassMatrixLumpedUBLAS;
            //M_opVelocity2dTo2dconf->apply( qqqq, *M_coulingRNG_operatorDiagonalOnFluid );
            auto rrrr = M_XhMeshVelocityInterface->element();
            M_opVelocity2dTo2dconf->apply( qqqq, rrrr );
            M_coulingRNG_operatorDiagonalOnFluid->on( _range=M_rangeFSI_fluid,_expr=idv(rrrr));
            sync( *M_coulingRNG_operatorDiagonalOnFluid, "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );
            useAlgebraicInnerProductWithLumping = true;
        }
#if 0
        else if ( true )
        {
            // compute algebraic counterpart of the operator B (my work V2 but not enought accurate)
            auto const& u = this->solidModel()->fieldDisplacement();
            auto massMatrixRestrictFSI = backend()->newMatrix(_test=Vh,_trial=Vh);
            form2(_test=Vh,_trial=Vh,_matrix=massMatrixRestrictFSI ) =
                integrate(_range=rangeFSI,_expr=inner( idt(u), id(u) ) );
            massMatrixRestrictFSI->close();
            auto vecDiagMassMatrixLumpedRestrictFSI = backend()->newVector(Vh);
            if ( Vh->fe()->nOrder==1)
            {
                auto unityVec = this->solidModel()->backend()->newVector(Vh);
                unityVec->setConstant(1);
                massMatrixRestrictFSI->multVector( unityVec,vecDiagMassMatrixLumpedRestrictFSI );
            }
            else
            {
                auto sumRow = this->solidModel()->backend()->newVector(Vh);
                auto unityVec = this->solidModel()->backend()->newVector(Vh);
                unityVec->setConstant(1);
                massMatrixRestrictFSI->multVector( unityVec,sumRow );
                double sumMatrix = sumRow->sum();
                vecDiagMassMatrixLumpedRestrictFSI = massMatrixRestrictFSI->diagonal();
                double sumDiag = vecDiagMassMatrixLumpedRestrictFSI->sum();
                vecDiagMassMatrixLumpedRestrictFSI->scale( sumMatrix/sumDiag );
            }

            auto vecDiagMassMatrixLumpedRestrictFSI_OP = Vh->element();
            auto vecDiagMassMatrixLumpedUBLAS = Vh->element( vecDiagMassLumped2 );
            for ( size_type k : dofsIdOnFSI )
                vecDiagMassMatrixLumpedRestrictFSI_OP(k) = vecDiagMassMatrixLumpedUBLAS(k)/vecDiagMassMatrixLumpedRestrictFSI->operator()(k);
            sync( vecDiagMassMatrixLumpedRestrictFSI_OP, "=", dofsIdOnFSI );

            //M_opVelocity2dTo2dconf->apply( vecDiagMassMatrixLumpedRestrictFSI_OP, *M_coulingRNG_operatorDiagonalOnFluid );
            auto rrrr = M_XhMeshVelocityInterface->element();
            M_opVelocity2dTo2dconf->apply( vecDiagMassMatrixLumpedRestrictFSI_OP, rrrr );
            M_coulingRNG_operatorDiagonalOnFluid->on( _range=M_rangeFSI_fluid,_expr=idv(rrrr));
            sync( *M_coulingRNG_operatorDiagonalOnFluid, "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );
        }
        else
        {
            // compute algebraic counterpart of the operator B (my work V1 but not enought accurate)
            auto vecDiagMassLumpedRestrictedFSI = Vh->element();
            auto markedDofsOnFSI = Vh->element();
            for ( size_type k : dofsIdOnFSI )
            {
                vecDiagMassLumpedRestrictedFSI.set( k, vecDiagMassLumped2->operator()(k) );
                markedDofsOnFSI.set( k, 1. );
            }
            sync( vecDiagMassLumpedRestrictedFSI, "=", dofsIdOnFSI );
            sync( markedDofsOnFSI, "=", dofsIdOnFSI );

            size_type countActiveDofOnFSILocal = 0;
            for ( size_type k=0;k<Vh->nLocalDofWithGhost();++k )
            {
                if ( markedDofsOnFSI(k) > 0.5 && !Vh->dof()->dofGlobalProcessIsGhost( k ) )
                    ++countActiveDofOnFSILocal;
            }
            size_type countActiveDofOnFSI = countActiveDofOnFSILocal;
            if ( this->solidModel()->worldComm().size() > 1 )
                mpi::all_reduce( this->solidModel()->worldComm(), countActiveDofOnFSILocal, countActiveDofOnFSI, std::plus<size_type>() );


            bool useOperatorProportionalToIdentity = boption(_name="coupling-robin-neumann-generalized.use-operator-proportional-to-identity",_prefix=this->prefix());
            if ( useOperatorProportionalToIdentity )
            {
                // alternative version (less accurate) :  compute a scalar c in order to define operator B = c*Id_{fsi}
                double sumDiagMassLumpedRestrictedFSI = vecDiagMassLumpedRestrictedFSI.sum();
                double myscalingbc = sumDiagMassLumpedRestrictedFSI/(Vh->nComponents*areaFSI);
                for ( size_type k : dofsIdOnFSI )
                    vecDiagMassLumpedRestrictedFSI.set( k,myscalingbc );
                sync( vecDiagMassLumpedRestrictedFSI, "=", dofsIdOnFSI );
            }
            else
            {
                // diagonal entries of the operator B
                vecDiagMassLumpedRestrictedFSI.scale( countActiveDofOnFSI/(Vh->nComponents*areaFSI) );
            }

            //M_opVelocity2dTo2dconf->apply( vecDiagMassLumpedRestrictedFSI, *M_coulingRNG_operatorDiagonalOnFluid );
            auto rrrr = M_XhMeshVelocityInterface->element();
            M_opVelocity2dTo2dconf->apply( vecDiagMassLumpedRestrictedFSI,rrrr );
            M_coulingRNG_operatorDiagonalOnFluid->on( _range=M_rangeFSI_fluid,_expr=idv(rrrr));
            sync( *M_coulingRNG_operatorDiagonalOnFluid, "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );
        }
#endif
    }
    else if ( this->solidModel()->is1dReducedModel() )
    {
#if 1 // NEW

        auto VhSolid1dVelocityVect = this->solidModel()->solid1dReduced()->fieldVelocityVect1dReduced().functionSpace();
        auto qqqq = VhSolid1dVelocityVect->element();
        auto se = this->solidModel()->symbolsExpr();
        for ( auto const& [physicName,physicData] :  this->solidModel()->solid1dReduced()->physicsFromCurrentType() )
        {
            //auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nRealDim>>(physicData);
            for ( std::string const& matName : this->solidModel()->solid1dReduced()->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->solidModel()->solid1dReduced()->materialsProperties()->rangeMeshElementsByMaterial( this->solidModel()->solid1dReduced()->mesh(),matName );
                auto const& densityProp = this->solidModel()->solid1dReduced()->materialsProperties()->density( matName );
                auto densityExpr = expr( densityProp.expr(), se );
                double thickness = this->solidModel()->solid1dReduced()->thickness1dReduced();
                qqqq.on(_range=range,_expr=thickness*densityExpr*one());
            }
        }

        auto rrrr = M_XhMeshVelocityInterface->element();
        M_opVelocity1dToNdconf->apply( qqqq, rrrr );

        auto VhFluid = this->fluidModel()->functionSpaceVelocity();//meshVelocity2().functionSpace();
        M_coulingRNG_operatorDiagonalOnFluid = VhFluid->elementPtr();
        M_coulingRNG_operatorDiagonalOnFluid->on( _range=M_rangeFSI_fluid,_expr=idv(rrrr));
        sync( *M_coulingRNG_operatorDiagonalOnFluid, "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );

#else
        auto VhFluid = this->fluidModel()->functionSpaceVelocity();//meshVelocity2().functionSpace();
        M_coulingRNG_operatorDiagonalOnFluid = VhFluid->elementPtr();
        std::set<size_type> dofsIdOnFSIFluid;
        for ( auto const& faceWrap : M_rangeFSI_fluid )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = VhFluid->dof()->faceLocalDof( face.id() );
            for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                dofsIdOnFSIFluid.insert( it->index() );
        }
        M_coulingRNG_operatorDiagonalOnFluid->on(_range=M_rangeFSI_fluid,_expr=one());
        sync( *M_coulingRNG_operatorDiagonalOnFluid, "=", dofsIdOnFSIFluid );
#endif
    }


    if ( true )
    {
        // create matrix which represent time derivative  bc operator
        auto dmFullFluidSpace = this->fluidModel()->algebraicFactory()->sparsityMatrixGraph()->mapRowPtr();
        M_coulingRNG_vectorTimeDerivative = this->fluidModel()->backend()->newVector( dmFullFluidSpace );
        int nBlock = this->fluidModel()->nBlockMatrixGraph();
        auto VhVelocity = this->fluidModel()->functionSpaceVelocity();
        auto ru = stencilRange<0,0>(M_rangeFSI_fluid);
        auto const& u = M_fluidModel->fieldVelocity();
        if ( useAlgebraicInnerProductWithLumping )
        {
            auto graph = std::make_shared<graph_type>( dmFullFluidSpace,dmFullFluidSpace );
            graph->addMissingZeroEntriesDiagonal();
            graph->close();
            M_coulingRNG_matrixTimeDerivative = this->fluidModel()->backend()->newMatrix( 0,0,0,0,graph );
            auto thediagg = this->fluidModel()->backend()->newVector( dmFullFluidSpace );
            auto thediaggUBLAS_u = VhVelocity->element( thediagg, 0 );
            thediaggUBLAS_u = *M_coulingRNG_operatorDiagonalOnFluid;
            M_coulingRNG_matrixTimeDerivative->setDiagonal( thediagg );
        }
        else
        {
            BlocksBaseGraphCSR myblockGraphTimeDerivative(nBlock,nBlock);
            auto mygraphTimeDerivative = stencil(_test=VhVelocity,_trial=VhVelocity,
                                                 _diag_is_nonzero=false,_close=true,
                                                 _range=stencilRangeMap(ru) )->graph();
            myblockGraphTimeDerivative(0,0) = mygraphTimeDerivative;
            for ( int k=1;k<nBlock;++k )
            {
                auto mapPtr = this->fluidModel()->algebraicBlockVectorSolution()->operator()(k)->mapPtr();
                graph_ptrtype zeroGraph = std::make_shared<graph_type>( mapPtr,mapPtr );
                zeroGraph->zero();
                myblockGraphTimeDerivative(k,k) = zeroGraph;
            }
            M_coulingRNG_matrixTimeDerivative = this->fluidModel()->backend()->newBlockMatrix(_block=myblockGraphTimeDerivative);
            auto myB = this->couplingRNG_operatorExpr( mpl::int_<fluid_type::nDim>() );
            form2(_test=VhVelocity,_trial=VhVelocity,_matrix=M_coulingRNG_matrixTimeDerivative/*,_pattern=size_type(Pattern::DEFAULT)*/ ) +=
                integrate( _range=M_rangeFSI_fluid,
                           _expr=inner(myB*idt(u),id(u)),
                           _geomap=this->geomap() );
            M_coulingRNG_matrixTimeDerivative->close();
        }

        // create matrix which represent stress bc operator
        auto VhStress = M_spaceNormalStress_fluid;//this->fluidModel()->fieldNormalStressRefMeshPtr()->functionSpace();
        M_coulingRNG_vectorStress = this->fluidModel()->backend()->newVector( VhStress );
        auto dofStress = VhStress->dof();
        BlocksBaseGraphCSR myblockGraph(nBlock,1);
        auto mygraph = stencil(_test=VhVelocity,_trial=VhStress,
                               _diag_is_nonzero=false,_close=true,
                               _range=stencilRangeMap(ru) )->graph();
        myblockGraph(0,0) = mygraph;
        for ( int k=1;k<nBlock;++k )
        {
            auto mapPtr = this->fluidModel()->algebraicBlockVectorSolution()->operator()(k)->mapPtr();
            graph_ptrtype zeroGraph = std::make_shared<graph_type>( mapPtr,dofStress );
            zeroGraph->zero();
            myblockGraph(k,0) = zeroGraph;
        }
        M_coulingRNG_matrixStress = this->fluidModel()->backend()->newBlockMatrix(_block=myblockGraph);
        // assembly stress matrix
        form2(_test=VhVelocity,_trial=VhStress,_matrix=M_coulingRNG_matrixStress ) +=
            integrate( _range=M_rangeFSI_fluid,
                       _expr= inner( idt(M_fieldNormalStressRefMesh_fluid),id(u)),
                       _geomap=this->geomap() );
        M_coulingRNG_matrixStress->close();
    }

}

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
FSI<FluidType,SolidType>::initInHousePreconditionerPCD_fluid( operatorpcdbase_fluid_type & opPCDBase ) const
{
    typedef Feel::Alternatives::OperatorPCD<typename fluid_type::space_velocity_type,typename fluid_type::space_pressure_type> op_pcd_type;
    op_pcd_type * opPCD = dynamic_cast<op_pcd_type*>(&opPCDBase);
    CHECK( opPCD ) << "fails to cast OperatorPCD";

    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        for ( auto const& rangeFacesMat : this->fluidModel()->rangeDistributionByMaterialName()->rangeMeshFacesByMaterial( "interface_fsi" ) )
        {
            std::string matName = rangeFacesMat.first;
            auto const& rangeFaces = rangeFacesMat.second;
            opPCD->addRangeDirichletBC( "FSI_FluidDirichlet_" + matName, rangeFaces );
        }
    }
    else
        opPCD->addRangeNeumannBC( "FSI_FluidNeumann", M_rangeFSI_fluid );
}
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateInHousePreconditionerPCD_fluid( operatorpcdbase_fluid_type & opPCDBase, DataUpdateBase & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        typedef Feel::Alternatives::OperatorPCD<typename fluid_type::space_velocity_type,typename fluid_type::space_pressure_type> op_pcd_type;
        op_pcd_type * opPCD = dynamic_cast<op_pcd_type*>(&opPCDBase);
        CHECK( opPCD ) << "fails to cast OperatorPCD";
        for ( auto const& rangeFacesMat : this->fluidModel()->rangeDistributionByMaterialName()->rangeMeshFacesByMaterial( "interface_fsi" ) )
        {
            std::string matName = rangeFacesMat.first;
            auto rhoExpr = this->fluidModel()->materialsProperties()->density( matName ).expr();
            opPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, "FSI_FluidDirichlet_" + matName, idv(this/*M_fluidModel*/->meshVelocity2()) );
        }
    }
}



//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );
    for ( auto [physicName,physicData] : this->physics/*FromCurrentType*/() )
        physicData->updateParameterValues( paramValues );

    this->setParameterValues( paramValues );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::setParameterValues( std::map<std::string,double> const& paramValues )
{
    // for ( auto const& [param,val] : paramValues )
    //     M_currentParameterValues[param] = val;

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_fluidModel->setParameterValues( paramValues );
    M_solidModel->setParameterValues( paramValues );
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl1()
{

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->transfertDisplacementAndApplyMeshMoving();
        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        //M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
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
            this->transfertDisplacementAndApplyMeshMoving();
            double tALE = timerCur.elapsed();
            this->log("FSI","update ale","finish in "+(boost::format("%1% s") % tALE).str() );
        }
        else  if (this->fsiCouplingType()=="Semi-Implicit")
        {
            this->transfertVelocity();
            double tALE = timerCur.elapsed();
            this->log("FSI","transfert velocity","finish in "+(boost::format("%1% s") % tALE).str());
        }
        //--------------------------------------------------------------//
        M_fluidModel->solve();
        //--------------------------------------------------------------//
        timerCur.restart();
        // revert ref mesh
        //M_fluidModel->meshMotionTool()->revertReferenceMesh();
        // transfert stress
        this->transfertStress();
        // revert moving mesh
        //M_fluidModel->meshMotionTool()->revertMovingMesh();
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
            //M_fluidModel->setRebuildCstPartInResidual(false);M_solidModel->setRebuildCstPartInResidual(false);
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
        this->transfertDisplacementAndApplyMeshMoving();
        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        //M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
    }
    bool useAitken = boption(_name="coupling-nitsche-family.use-aitken",_prefix=this->prefix());
    if ( useAitken )
        this->aitkenRelaxTool()->restart();

    bool solveStruct = true;//this->fluidModel()->timeStepBDF()->iteration() > 1 || this->fixPointMaxIt()==1;

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
            //M_fluidModel->meshMotionTool()->revertReferenceMesh();
            // transfert stress
            this->transfertStress();
            // revert moving mesh
            //M_fluidModel->meshMotionTool()->revertMovingMesh();
            if ( ( this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" ||
                   this->fsiCouplingBoundaryCondition()=="nitsche" ) &&
                 M_solidModel->isStandardModel() )
            {
                bool useExtrap = boption(_prefix=this->prefix(),_name="transfert-velocity-F2S.use-extrapolation");
                this->transfertVelocityF2S(cptFSI,useExtrap);
                if ( this->fluidModel()->hasNonNewtonianViscosity() )
                    this->transfertGradVelocityF2S();
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
            this->transfertDisplacementAndApplyMeshMoving();
        }
        else
        {
            this->transfertVelocity();
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
            //M_fluidModel->setRebuildCstPartInResidual(false);

            if (solveStruct)
            {
                M_solidModel->setRebuildLinearPartInJacobian(false);M_solidModel->setRebuildCstPartInLinearSystem(false);
                //M_solidModel->setRebuildCstPartInResidual(false);
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

    if (false && this->fsiCouplingType()=="Semi-Implicit")
    {
        this->transfertDisplacementAndApplyMeshMoving();
    }


}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::solveImpl3()
{
    boost::mpi::timer timerIter;

    if (this->fsiCouplingType()=="Semi-Implicit")
    {
        this->transfertDisplacementAndApplyMeshMoving();
        M_fluidModel->setRebuildLinearPartInJacobian(true);M_fluidModel->setRebuildCstPartInLinearSystem(true);
        M_solidModel->setRebuildLinearPartInJacobian(true);M_solidModel->setRebuildCstPartInLinearSystem(true);
        //M_fluidModel->setRebuildCstPartInResidual(true);M_solidModel->setRebuildCstPartInResidual(true);
    }

#if 0
    // predictor disp : WARNING explicit is instable
    M_solidModel->predictorDispl();
    M_solidModel->updateVelocity();
#endif

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
        if ( false )
            this->updateBackendOptimisation(cptFSI,residualRelativeConvergence);

        if ( useAitken )
            this->aitkenRelaxTool()->saveOldSolution();
        else
            M_fixPointConvergenceFSI->saveOldSolution();


        //--------------------------------------------------------------//
        if (this->fsiCouplingType()=="Implicit")
        {
            this->transfertDisplacementAndApplyMeshMoving();
        }
        else
        {
            //bool useExtrap = (this->fluidModel()->timeStepBDF()->iteration() > 1) && cptFSI==0;
            //this->transfertVelocity(useExtrap);
        }
        this->transfertRobinNeumannGeneralizedS2F( cptFSI );
        M_fluidModel->solve();

        //--------------------------------------------------------------//
        //M_fluidModel->meshMotionTool()->revertReferenceMesh();
        // transfert stress
        this->transfertStress();
        // revert moving mesh
        //M_fluidModel->meshMotionTool()->revertMovingMesh();
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
            //M_fluidModel->setRebuildCstPartInResidual(false);
            M_solidModel->setRebuildLinearPartInJacobian(false);M_solidModel->setRebuildCstPartInLinearSystem(false);
            //M_solidModel->setRebuildCstPartInResidual(false);
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
FSI<FluidType,SolidType>::startTimeStep()
{
    // need to transfert some quantities because theta scheme require a residual evaluation at previous
    // this evaluation is done when call M_solidModel->startTimeStep()
    if ( M_solidModel->timeStepping() == "Theta" )
    {
        this->transfertStress();
        if ( ( this->fsiCouplingBoundaryCondition()=="robin-robin" || this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" ||
               this->fsiCouplingBoundaryCondition()=="nitsche" ) &&
             M_solidModel->isStandardModel() )
        {
            bool useExtrap = boption(_prefix=this->prefix(),_name="transfert-velocity-F2S.use-extrapolation");
            this->transfertVelocityF2S(0,useExtrap);
        }
    }


    M_fluidModel->startTimeStep();
    M_solidModel->startTimeStep();
    this->updateTime( M_fluidModel->currentTime() );
    this->updateParameterValues();
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateTimeStep()
{
    M_previousTimeOrder=M_fluidModel->timeStepBDF()->timeOrder();

    M_fluidModel->updateTimeStep();
    M_solidModel->updateTimeStep();

    M_currentTimeOrder=M_fluidModel->timeStepBDF()->timeOrder();

    this->updateTime( M_fluidModel->currentTime() );
    this->updateParameterValues();
}

//---------------------------------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
std::shared_ptr<std::ostringstream>
FSI<FluidType,SolidType>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
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
        *_ostr << "\n   Generalized Robin-Neumann";
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
#endif
    return _ostr;
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateInformationObject( nl::json & p ) const
{
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    p["Toolbox Fluid"] = M_fluidModel->journalSection().to_string();
    p["Toolbox Solid"] = M_solidModel->journalSection().to_string();

    nl::json subPt;
    subPt["type"] = this->fsiCouplingType();
    subPt["BC"] = this->fsiCouplingBoundaryCondition();
    subPt["Interface property"] = M_interfaceFSIisConforme? "conformal":"non conformal";

    bool useAitken = false;
    if ( this->fsiCouplingBoundaryCondition()  == "dirichlet-neumann" )
    {
        useAitken = true;
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
              this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
              this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        nl::json subPt2;
        subPt2["gamma"] = M_couplingNitscheFamily_gamma;
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
              this->fsiCouplingBoundaryCondition() == "nitsche" )
            subPt2["gamma0"] = M_couplingNitscheFamily_gamma0;
        subPt2["alpha"] = M_couplingNitscheFamily_alpha;
        subPt["Nitsche family parameters"] = subPt2;
    }
    else if (this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
    }

    p["Coupling"] = std::move( subPt );

    subPt.clear();
    subPt["method"] = useAitken? "fix point with Aitken relaxation" : "fix point";
    subPt["tolerance"] = M_fixPointTolerance;
    subPt["maxit"] = M_fixPointMaxIt;
    if ( useAitken )
    {
        subPt["initial theta"] = M_fixPointInitialTheta;
        subPt["min theta"] = M_fixPointMinTheta;
    }
    p["Numerical Solver"] = std::move( subPt );
}

template< class FluidType, class SolidType >
tabulate_informations_ptr_t
FSI<FluidType,SolidType>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( jsonInfo.contains( "Coupling" ) )
    {
        auto tabInfoCoupling = TabulateInformationsSections::New( tabInfoProp );
        auto const& jsonInfoCoupling = jsonInfo.at("Coupling");
        Feel::Table tabInfoCouplingOther;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoCouplingOther, jsonInfoCoupling, tabInfoProp );
        tabInfoCouplingOther.format().setShowAllBorders( false ).setColumnSeparator(":").setHasRowSeparator( false );
        tabInfoCoupling->add( "", TabulateInformations::New( tabInfoCouplingOther, tabInfoProp ) );
        if ( jsonInfoCoupling.contains( "Nitsche family parameters" ) )
        {
            Feel::Table tabInfoCouplingNitsche;
            TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoCouplingNitsche, jsonInfoCoupling.at( "Nitsche family parameters" ), tabInfoProp );
            tabInfoCouplingNitsche.format().setShowAllBorders( false ).setColumnSeparator(":").setHasRowSeparator( false );
            tabInfoCoupling->add( "Nitsche family parameters", TabulateInformations::New( tabInfoCouplingNitsche, tabInfoProp ) );
        }
        tabInfo->add( "Coupling", tabInfoCoupling );
    }

    // Numerical Solver
    if ( jsonInfo.contains( "Numerical Solver" ) )
    {
        Feel::Table tabInfoNumSolver;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoNumSolver, jsonInfo.at("Numerical Solver"), tabInfoProp );
        tabInfoNumSolver.format().setShowAllBorders( false ).setColumnSeparator(":").setHasRowSeparator( false );
        tabInfo->add( "Numerical Solver",  TabulateInformations::New( tabInfoNumSolver, tabInfoProp ) );
    }

    // generate sub toolboxes info
    if ( M_fluidModel && jsonInfo.contains( "Toolbox Fluid" ) )
    {
        nl::json::json_pointer jsonPointerFluid( jsonInfo.at( "Toolbox Fluid" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerFluid ) )
        {
            auto tabInfos_fluid = M_fluidModel->tabulateInformations( JournalManager::journalData().at( jsonPointerFluid ), tabInfoProp );
            tabInfo->add( "Toolbox Fluid", tabInfos_fluid );
        }
    }
    if ( M_solidModel && jsonInfo.contains( "Toolbox Solid" ) )
    {
        nl::json::json_pointer jsonPointerSolid( jsonInfo.at( "Toolbox Solid" ).template get<std::string>() );
        if ( JournalManager::journalData().contains( jsonPointerSolid ) )
        {
            auto tabInfos_solid = M_solidModel->tabulateInformations( JournalManager::journalData().at( jsonPointerSolid ), tabInfoProp );
            tabInfo->add( "Toolbox Solid", tabInfos_solid );
        }
    }

    return tabInfo;
}

} // namespace FeelModels
} // namespace Feel
