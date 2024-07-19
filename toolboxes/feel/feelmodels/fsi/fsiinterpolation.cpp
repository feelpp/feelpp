
namespace Feel
{
namespace FeelModels
{
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initInterpolation()
{


    boost::mpi::timer btime;double thet;
    if (this->solidModel()->isStandardModel())
    {
        //---------------------------------------------//
        this->initDispInterpolation();
        std::ostringstream ostr1;ostr1<<btime.elapsed()<<"s";
        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initDispInterpolation","done in "+ostr1.str(),
                                                   this->worldComm(),this->verboseAllProc());
        btime.restart();
        //---------------------------------------------//
        this->initStressInterpolation();
        std::ostringstream ostr2;ostr2<<btime.elapsed()<<"s";
        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStressInterpolation","done in "+ostr2.str(),
                                                   this->worldComm(),this->verboseAllProc());
        btime.restart();
        //---------------------------------------------//
        if (this->fsiCouplingType()=="Semi-Implicit" ||
            this->fsiCouplingBoundaryCondition()=="robin-neumann-generalized" ||
            this->fsiCouplingBoundaryCondition()=="robin-neumann" || this->fsiCouplingBoundaryCondition()=="robin-robin" ||
            this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" || this->fsiCouplingBoundaryCondition()=="nitsche" )
        {
            this->initVelocityInterpolation();
            std::ostringstream ostr3;ostr3<<btime.elapsed()<<"s";
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initVelocityInterpolation","done in "+ostr3.str(),
                                                       this->worldComm(),this->verboseAllProc());
        }
        btime.restart();
        //---------------------------------------------//
#if 0
        if (this->fsiCouplingBoundaryCondition()=="robin-neumann")
        {
            this->initStressInterpolationS2F();
            std::ostringstream ostr4;ostr4<<btime.elapsed()<<"s";
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStressInterpolationS2F","done in "+ostr4.str(),
                                                       this->worldComm(),this->verboseAllProc());

        }
#endif
        //---------------------------------------------//
        if ( /*this->fsiCouplingBoundaryCondition()=="robin-neumann" ||*/ this->fsiCouplingBoundaryCondition()=="robin-robin" ||
             this->fsiCouplingBoundaryCondition()=="robin-robin-genuine" || this->fsiCouplingBoundaryCondition()=="nitsche" )
        {
            this->initVelocityInterpolationF2S();
        }
    }
    else if ( this->solidModel()->is1dReducedModel() )
    {
        //---------------------------------------------//
        this->initDisp1dToNdInterpolation();
        thet = btime.elapsed();
        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initDisp1dToNdInterpolation",
                                                   (boost::format("done in %1%s")% thet).str(),
                                                   this->worldComm(),this->verboseAllProc());
        btime.restart();
        //---------------------------------------------//
        this->initStress1dToNdInterpolation();
        thet = btime.elapsed();
        if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initStress1dToNdInterpolation",
                                                   (boost::format("done in %1%s")% thet).str(),
                                                   this->worldComm(),this->verboseAllProc());
        btime.restart();
        //---------------------------------------------//
        if (this->fsiCouplingType()=="Semi-Implicit" || this->fsiCouplingBoundaryCondition()=="robin-neumann-generalized")
        {
            this->initVelocity1dToNdInterpolation();
            thet = btime.elapsed();
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","initVelocity1dToNdInterpolation",
                                                       (boost::format("done in %1%s")% thet).str(),
                                                       this->worldComm(),this->verboseAllProc());
        }
        //---------------------------------------------//
    }

}

//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//



template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initDispInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank() )
            std::cout << "initDispInterpolation() CONFORME"  << std::endl;
        M_opDisp2dTo2dconf = opInterpolation(_domainSpace=this->solidModel()->functionSpaceDisplacement(),
                                             _imageSpace=this->fluidModel()->meshMotionTool()->displacement()->functionSpace(),
                                             _range=M_rangeFSI_fluid,
                                             _type=InterpolationConforme(),
                                             _backend=this->fluidModel()->backend() );
    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initDispInterpolation() NONCONFORME" << std::endl;
        M_opDisp2dTo2dnonconf = opInterpolation(_domainSpace=this->solidModel()->functionSpaceDisplacement(),
                                                _imageSpace=this->fluidModel()->meshMotionTool()->displacement()->functionSpace(),
                                                _range=M_rangeFSI_fluid,
                                                _type=InterpolationNonConforme(),
                                                _backend=this->fluidModel()->backend() );
    }
}


//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initDisp1dToNdInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initDisp1dToNdInterpolation() CONFORME"  << std::endl;
        M_opDisp1dToNdconf = opInterpolation(_domainSpace=this->solidModel()->solid1dReduced()->fieldDisplacementVect1dReduced().functionSpace(),
                                             _imageSpace=this->fluidModel()->meshMotionTool()->displacement()->functionSpace(),//M_fluid->meshDisplacementOnInterface().functionSpace(),
                                             _range=M_rangeFSI_fluid,
                                             _type=InterpolationConforme(),
                                             _backend=this->fluidModel()->backend() );
    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initDisp1dToNdInterpolation() NONCONFORME" << std::endl;
        M_opDisp1dToNdnonconf = opInterpolation(_domainSpace=this->solidModel()->solid1dReduced()->fieldDisplacementVect1dReduced().functionSpace(),
                                                _imageSpace=this->fluidModel()->meshMotionTool()->displacement()->functionSpace(),//M_fluid->meshDisplacementOnInterface().functionSpace(),
                                                _range=M_rangeFSI_fluid,
                                                _type=InterpolationNonConforme(),
                                                _backend=this->fluidModel()->backend() );
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initStressInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initStressInterpolation() CONFORME" << std::endl;
        M_opStress2dTo2dconf = opInterpolation(_domainSpace=M_spaceNormalStress_fluid,
                                               _imageSpace=M_spaceNormalStressFromFluid_solid,
                                               _range=elements(M_spaceNormalStressFromFluid_solid->mesh()),
                                               _type=InterpolationConforme(),
                                               _backend=this->fluidModel()->backend() );
    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initStressInterpolation() NONCONFORME" << std::endl;
        M_opStress2dTo2dnonconf = opInterpolation(_domainSpace=M_spaceNormalStress_fluid,
                                                  _imageSpace=M_spaceNormalStressFromFluid_solid,
                                                  _range=elements(M_spaceNormalStressFromFluid_solid->mesh()),
                                                  _type=InterpolationNonConforme(true,true,true,15),
                                                  _backend=this->fluidModel()->backend() );
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initStress1dToNdInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initStress1dToNdInterpolation() CONFORME" << std::endl;
        M_opStress1dToNdconf = opInterpolation(_domainSpace=M_spaceNormalStress_fluid,
                                               _imageSpace=M_spaceNormalStressFromFluid_solid1dReduced,
                                               _range=elements(M_solidModel->solid1dReduced()->mesh()),
                                               _type=InterpolationConforme(),
                                               _backend=M_fluidModel->backend() );

    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initStress1dToNdInterpolation() NONCONFORME" << std::endl;
        M_opStress1dToNdnonconf = opInterpolation(_domainSpace=M_spaceNormalStress_fluid,
                                                  _imageSpace=M_spaceNormalStressFromFluid_solid1dReduced,
                                                  _range=elements(M_solidModel->solid1dReduced()->mesh()),
                                                  _type=InterpolationNonConforme(),
                                                  _backend=M_fluidModel->backend() );
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initVelocityInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocityInterpolation() CONFORME" << std::endl;
        M_opVelocity2dTo2dconf = opInterpolation(_domainSpace=this->solidModel()->fieldVelocity().functionSpace(),
                                                 _imageSpace=this->/*fluidModel()->*/meshVelocity2().functionSpace(),
                                                 //_range=M_rangeFSI_fluid,
                                                 _type=InterpolationConforme(),
                                                 _backend=this->fluidModel()->backend() );

    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocityInterpolation() NONCONFORME" << std::endl;
        M_opVelocity2dTo2dnonconf = opInterpolation(_domainSpace=this->solidModel()->fieldVelocity().functionSpace(),
                                                    _imageSpace=this->/*fluidModel()->*/meshVelocity2().functionSpace(),
                                                    //_range=M_rangeFSI_fluid,
                                                    _type=InterpolationNonConforme(),
                                                    _backend=this->fluidModel()->backend() );
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initVelocity1dToNdInterpolation()
{
    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocity1dToNdInterpolation() CONFORME" << std::endl;
        M_opVelocity1dToNdconf = opInterpolation(_domainSpace=this->solidModel()->solid1dReduced()->fieldVelocityVect1dReduced().functionSpace(),
                                                 _imageSpace=this->/*fluidModel()->*/meshVelocity2().functionSpace(),
                                                 //_range=M_rangeFSI_fluid,
                                                 _type=InterpolationConforme(),
                                                 _backend=this->fluidModel()->backend() );
    }
    else
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocity1dToNdInterpolation() NONCONFORME" << std::endl;
        M_opVelocity1dToNdnonconf = opInterpolation(_domainSpace=this->solidModel()->solid1dReduced()->fieldVelocityVect1dReduced().functionSpace(),
                                                    _imageSpace=this->/*fluidModel()->*/meshVelocity2().functionSpace(),
                                                    //_range=M_rangeFSI_fluid,
                                                    _type=InterpolationNonConforme(),
                                                    _backend=this->fluidModel()->backend() );
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initVelocityInterpolationF2S()
{

    if (M_interfaceFSIisConforme)
    {
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocityInterpolationF2S() CONFORME" << std::endl;
        M_opVelocity2dTo2dconfF2S = opInterpolation(_domainSpace=this->fluidModel()->functionSpaceVelocity(),
                                                    _imageSpace=this->fieldVelocityInterfaceFromFluidPtr_solid()->functionSpace(),
                                                    _range=M_rangeFSI_solid,
                                                    _type=InterpolationConforme(),
                                                    _backend=this->fluidModel()->backend() );

    }
    else
    {
        CHECK( false ) << "TODO";
#if 0
        if (this->verbose() && this->fluidModel()->worldComm().isMasterRank())
            std::cout << "initVelocityInterpolation() NONCONFORME" << std::endl;
        M_opVelocity2dTo2dnonconfF2S = opInterpolation(_domainSpace=this->fluidModel()->meshVelocity2().functionSpace(),
                                                       _imageSpace=this->solidModel()->fieldVelocity()->functionSpace(),
                                                       _range=M_rangeFSI_solid,
                                                       _type=InterpolationNonConforme(),
                                                       _backend=this->fluidModel()->backend() );
#endif
    }
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::initStressInterpolationS2F()
{
#if 0
    M_opStress2dTo2dconfS2F = opInterpolation(_domainSpace=this->solidModel()->fieldNormalStressFromStructPtr()->functionSpace(),
                                              _imageSpace=this->fluidModel()->normalStressFromStruct()->functionSpace(),
                                              _range=M_rangeFSI_fluid,
                                              _type=InterpolationConforme(),
                                              _backend=this->fluidModel()->backend() );
#else
    CHECK( false ) << "not implemented";
#endif
}

//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertDisplacement()
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertDisplacement", "start",
                                               this->worldComm(),this->verboseAllProc());

    if (M_solidModel->isStandardModel())
    {
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opDisp2dTo2dconf ) << "interpolation operator not build";
            M_opDisp2dTo2dconf->apply(M_solidModel->fieldDisplacement(),
                                      *M_meshDisplacementOnInterface_fluid/*  *(M_fluidModel->meshDisplacementOnInterface()) */ );
        }
        else
        {
            CHECK( M_opDisp2dTo2dnonconf ) << "interpolation operator not build";
            M_opDisp2dTo2dnonconf->apply(M_solidModel->fieldDisplacement(),
                                         *M_meshDisplacementOnInterface_fluid /**(M_fluidModel->meshDisplacementOnInterface())*/ );
        }
    }
    else if ( M_solidModel->is1dReducedModel() )
    {
        M_solidModel->solid1dReduced()->updateInterfaceDispFrom1dDisp();
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opDisp1dToNdconf ) << "interpolation operator not build";
            M_opDisp1dToNdconf->apply(M_solidModel->solid1dReduced()->fieldDisplacementVect1dReduced(),
                                      *M_meshDisplacementOnInterface_fluid /**(M_fluidModel->meshDisplacementOnInterface())*/ );
        }
        else
        {
            CHECK( M_opDisp1dToNdnonconf ) << "interpolation operator not build";
            M_opDisp1dToNdnonconf->apply(M_solidModel->solid1dReduced()->fieldDisplacementVect1dReduced(),
                                         *M_meshDisplacementOnInterface_fluid /**(M_fluidModel->meshDisplacementOnInterface())*/ );
        }
    }
    else
        CHECK( false ) << "something wrong";

    this->fluidModel()->meshMotionTool()->updateDisplacementImposed( idv(M_meshDisplacementOnInterface_fluid), M_rangeFSI_fluid );

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertDisplacement", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertDisplacementAndApplyMeshMoving()
{
    this->transfertDisplacement();
    M_fluidModel->updateALEmesh();
    // up mesh velocity on interface from mesh velocity
    M_meshVelocityInterface->on( _range=M_rangeFSI_fluid,
                                 _expr=vf::idv(M_fluidModel->meshMotionTool()->velocity()),
                                 _geomap=this->geomap() );
    //sync( *M_meshVelocityInterface, "=", M_dofsVelocityInterfaceOnMovingBoundary);
}

//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertStress()
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStress", "start",
                                               this->worldComm(),this->verboseAllProc());
    M_fieldNormalStressRefMesh_fluid->zero();
    M_fluidModel->updateNormalStressOnReferenceMesh( "interface_fsi", M_fieldNormalStressRefMesh_fluid );

    if (M_solidModel->isStandardModel())
    {
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opStress2dTo2dconf ) << "interpolation operator not build";
            M_opStress2dTo2dconf->apply( *M_fieldNormalStressRefMesh_fluid, *M_fieldNormalStressFromFluid_solid );
        }
        else
        {
            CHECK( M_opStress2dTo2dnonconf ) << "interpolation operator not build";
            M_opStress2dTo2dnonconf->apply( *M_fieldNormalStressRefMesh_fluid, *M_fieldNormalStressFromFluid_solid );
        }
    }
    else if ( M_solidModel->is1dReducedModel() )
    {
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opStress1dToNdconf )  << "interpolation operator not build";
            M_opStress1dToNdconf->apply(*M_fieldNormalStressRefMesh_fluid, *M_fieldNormalStressFromFluidVectorial_solid1dReduced );
        }
        else
        {
            CHECK( M_opStress1dToNdnonconf )  << "interpolation operator not build";
            M_opStress1dToNdnonconf->apply(*M_fieldNormalStressRefMesh_fluid, *M_fieldNormalStressFromFluidVectorial_solid1dReduced );
        }
#if 0
        M_solidModel->updateInterfaceScalStressDispFromVectStress();
#else
        // TODO : create generic method for this
        M_fieldNormalStressFromFluidScalar_solid1dReduced->on( _range=elements(M_solidModel->solid1dReduced()->mesh()),
                                                               _expr=-inner(idv(M_fieldNormalStressFromFluidVectorial_solid1dReduced),oneY()) );
#endif
    }

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStress", "finish",
                                               this->worldComm(),this->verboseAllProc());

}

//-----------------------------------------------------------------------------------//

#if 0
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannInterfaceOperatorS2F()
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannInterfaceOperatorS2F", "start",
                                               this->worldComm(),this->verboseAllProc());

    if (!M_solidModel->isStandardModel()) return;

    auto fieldInterpolated = this->fluidModel()->functionSpaceVelocity()->elementPtr();//M_fluid->meshVelocity2().functionSpace()->elementPtr();
    auto fieldToTransfert = M_solidModel->functionSpaceDisplacement()->elementPtr(); //fieldVelocityPtr()->functionSpace()->elementPtr();
    CHECK( fieldToTransfert->map().nLocalDofWithGhost() == this->robinNeumannInterfaceOperator()->map().nLocalDofWithGhost() ) << "invalid compatibility size";
    *fieldToTransfert = *this->robinNeumannInterfaceOperator();
    if (M_interfaceFSIisConforme)
    {
        M_opVelocityBis2dTo2dconf/*auto opI*/ = opInterpolation(_domainSpace=this->solidModel()->functionSpaceDisplacement(),
                                   _imageSpace=this->fluidModel()->functionSpaceVelocity(),
                                                                _range=M_rangeFSI_fluid,
                                   _type=InterpolationConforme(),
                                   _backend=this->fluidModel()->backend() );
        M_opVelocityBis2dTo2dconf/*opI*/->apply( *fieldToTransfert, *fieldInterpolated );
    }
    else
    {
        M_opVelocityBis2dTo2dnonconf/*auto opI*/ = opInterpolation(_domainSpace=this->solidModel()->functionSpaceDisplacement(),
                                   _imageSpace=this->fluidModel()->functionSpaceVelocity(),
                                                                   _range=M_rangeFSI_fluid,
                                   _type=InterpolationNonConforme(),
                                   _backend=this->fluidModel()->backend() );
        M_opVelocityBis2dTo2dnonconf/*opI*/->apply( *fieldToTransfert, *fieldInterpolated );
    }

    M_fluid->setCouplingFSI_RNG_interfaceOperator( fieldInterpolated );
    M_fluid->setCouplingFSI_RNG_useInterfaceOperator( true );
    M_fluid->couplingFSI_RNG_updateForUse();

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannInterfaceOperatorS2F", "finish",
                                               this->worldComm(),this->verboseAllProc());
}
#endif

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertVelocity(bool useExtrap)
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertVelocity", "start",
                                               this->worldComm(),this->verboseAllProc());

    if (M_solidModel->isStandardModel())
    {
        typename solid_type::element_displacement_ptrtype fieldToTransfert;
        if ( !useExtrap )
        {
            fieldToTransfert = M_solidModel->fieldVelocityPtr();
        }
        else
        {
            fieldToTransfert = M_solidModel->fieldVelocity().functionSpace()->elementPtr();
            fieldToTransfert->add(  2.0, this->solidModel()->timeStepNewmark()->previousVelocity() );
            fieldToTransfert->add( -1.0, this->solidModel()->timeStepNewmark()->previousVelocity(1) );
        }

        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dconf->apply( *fieldToTransfert, this/*M_fluidModel*/->meshVelocity2() );
        }
        else
        {
            CHECK( M_opVelocity2dTo2dnonconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dnonconf->apply( *fieldToTransfert, this/*M_fluidModel*/->meshVelocity2() );
        }
    }
    else if ( M_solidModel->is1dReducedModel() )
    {

        typename solid_type::solid_1dreduced_type::element_displacement_ptrtype/*element_vect_1dreduced_ptrtype*/ fieldToTransfert;
        if( !useExtrap )
        {
            M_solidModel->solid1dReduced()->updateInterfaceVelocityFrom1dVelocity();
            fieldToTransfert = M_solidModel->solid1dReduced()->fieldVelocityVect1dReducedPtr();
        }
        else
        {
            if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertVelocity", "use extrapolation (1dReduced)",
                                                       this->worldComm(),this->verboseAllProc());
            //auto fieldExtrapolated = M_solidModel->fieldVelocityScal1dReduced().functionSpace()->elementPtr();
            auto fieldExtrapolated = this->solidModel()->solid1dReduced()->timeStepNewmark()->previousVelocity().functionSpace()->elementPtr();
            fieldExtrapolated->add(  2.0, this->solidModel()->solid1dReduced()->timeStepNewmark()->previousVelocity() );
            fieldExtrapolated->add( -1.0, this->solidModel()->solid1dReduced()->timeStepNewmark()->previousVelocity(1) );
            fieldToTransfert = M_solidModel->solid1dReduced()->extendVelocity1dReducedVectorial( *fieldExtrapolated );
        }

        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity1dToNdconf ) << "interpolation operator not build";
            M_opVelocity1dToNdconf->apply(*fieldToTransfert/*M_solid->fieldVelocityVect1dReduced()*/,
                                          this/*M_fluidModel*/->meshVelocity2() );
        }
        else
        {
            CHECK( M_opVelocity1dToNdnonconf ) << "interpolation operator not build";
            M_opVelocity1dToNdnonconf->apply(*fieldToTransfert/*M_solid->fieldVelocityVect1dReduced()*/,
                                             this/*M_fluidModel*/->meshVelocity2() );
        }
    }

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertVelocity", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

#if 0
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannGeneralizedS2F( int iterationFSI )
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "start",
                                               this->worldComm(),this->verboseAllProc());

    if ( this->solidModel()->isStandardModel())
    {
        auto fieldToTransfert = M_solidModel->fieldVelocity().functionSpace()->elementPtr();
        double dt = M_solidModel->timeStepNewmark()->timeStep();
        double gamma = M_solidModel->timeStepNewmark()->gamma();
        double beta = M_solidModel->timeStepNewmark()->beta();
        // time derivative acceleration in solid
        if ( true/*false*/ )
            fieldToTransfert->add( -1, M_solidModel->timeStepNewmark()->currentAcceleration());
        else
        {
            if ( iterationFSI == 0 )
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i+1) );
            }
            else
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i) );
            }
        }

        for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
            fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                   this->solidModel()->timeStepNewmark()->previousVelocity(i) );

        // transfer solid to fluid
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dconf->apply( *fieldToTransfert, *this->couplingRNG_evalForm1() );
        }
        else
        {
            CHECK( M_opVelocity2dTo2dnonconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dnonconf->apply( *fieldToTransfert, *this->couplingRNG_evalForm1() );
        }
        // time derivative acceleration in fluid
#if 0
        auto UPolyDeriv = this->fluidModel()->timeStepBDF()->polyDeriv();
        auto uPolyDeriv = UPolyDeriv.template element<0>();
        this->couplingRNG_evalForm1()->add( -1.0, uPolyDeriv );
        M_couplingRNG_coeffForm2 = this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0);
#else
        M_couplingRNG_coeffForm2 = this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0);
#endif

    }
    else if ( this->solidModel()->is1dReducedModel() )
    {
        typename solid_type::element_vect_1dreduced_ptrtype fieldToTransfert;
        auto timeStepNewmark1dReduced = this->solidModel()->solid1dReduced()->timeStepNewmark();
        auto fieldExtrapolated2 = timeStepNewmark1dReduced->previousVelocity().functionSpace()->elementPtr();
        double dt = timeStepNewmark1dReduced->timeStep();
        double gamma = timeStepNewmark1dReduced->gamma();
        double beta = timeStepNewmark1dReduced->beta();
        double scaleTimeDisc = M_solidModel->mechanicalProperties()->cstRho()*M_solidModel->solid1dReduced()->thickness1dReduced();
        // time derivative acceleration in solid
        if ( true )
            fieldExtrapolated2->add( -1, timeStepNewmark1dReduced->currentAcceleration());
        else
        {
            if ( iterationFSI == 0 )
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         timeStepNewmark1dReduced->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             timeStepNewmark1dReduced->previousVelocity(i+1) );
            }
            else
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         timeStepNewmark1dReduced->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             timeStepNewmark1dReduced->previousVelocity(i) );
            }
        }
        // time derivative acceleration in fluid
        for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
            fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                     timeStepNewmark1dReduced->previousVelocity(i) );
        fieldExtrapolated2->scale( scaleTimeDisc );
        fieldToTransfert = M_solidModel->extendVelocity1dReducedVectorial( *fieldExtrapolated2 );
        // transfer solid to fluid
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity1dToNdconf ) << "interpolation operator not build";
            M_opVelocity1dToNdconf->apply(*fieldToTransfert,
                                          *this->couplingRNG_evalForm1() );
        }
        else
        {
            CHECK( M_opVelocity1dToNdnonconf ) << "interpolation operator not build";
            M_opVelocity1dToNdnonconf->apply(*fieldToTransfert,
                                             *this->couplingRNG_evalForm1() );
        }
        // time derivative acceleration in fluid
        //auto UPolyDeriv = this->fluidModel()->timeStepBDF()->polyDeriv();
        //auto uPolyDeriv = UPolyDeriv.template element<0>();
        //this->couplingRNG_evalForm1()->add( -scaleTimeDisc, uPolyDeriv );

        M_couplingRNG_coeffForm2 = scaleTimeDisc*this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0);
    }

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "finish",
                                               this->worldComm(),this->verboseAllProc());
}
#elif 0
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannGeneralizedS2F( int iterationFSI )
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "start",
                                               this->worldComm(),this->verboseAllProc());
    bool useOriginalMethod = true;
    if (M_solidModel->isStandardModel())
    {
        auto fieldToTransfert = M_solidModel->fieldVelocity().functionSpace()->elementPtr();
        //typename solid_type::element_displacement_ptrtype fieldToTransfert;
        double dt = M_solidModel->timeStepNewmark()->timeStep();
        double gamma = M_solidModel->timeStepNewmark()->gamma();
        double beta = M_solidModel->timeStepNewmark()->beta();

#if 0
        if ( useOriginalMethod || (iterationFSI == 0) )
            fieldToTransfert->add( -1./(dt*gamma) , M_solidModel->timeStepNewmark()->previousVelocity() );
            //fieldToTransfert->add( (1./(dt*gamma))*( (gamma/beta)-1. ) -1./(beta*dt) , M_solidModel->timeStepNewmark()->previousVelocity() );
        //else
        //fieldToTransfert->add( (1./(dt*gamma))*( (gamma/beta)-1. ) -1./(beta*dt) , M_solidModel->timeStepNewmark()->currentVelocity() );

        if ( useOriginalMethod || (iterationFSI == 0) )
        {
            //fieldToTransfert->add( (1./gamma)*(gamma/(2*beta) - 1) - (1./(2*beta) -1), M_solidModel->timeStepNewmark()->previousAcceleration() );
            fieldToTransfert->add( -(1.-gamma)/gamma, M_solidModel->timeStepNewmark()->previousAcceleration() );
            fieldToTransfert->add( -1.0, M_solidModel->timeStepNewmark()->currentAcceleration() );
        }
        else
        {
            //fieldToTransfert->add( (1./gamma)*(gamma/(2*beta) - 1) - (1./(2*beta) -1) + 1.0, M_solidModel->timeStepNewmark()->currentAcceleration() );
            //fieldToTransfert->add( (1./gamma)*(gamma/(2*beta) - 1) - (1./(2*beta) -1) - 1.0, M_solidModel->timeStepNewmark()->currentAcceleration() );
            fieldToTransfert->add(-(1./(dt*gamma)), M_solidModel->timeStepNewmark()->currentVelocity() );
        }
#else
        fieldToTransfert->add( -1./(dt*gamma) , M_solidModel->timeStepNewmark()->previousVelocity() );
        fieldToTransfert->add( -(1.-gamma)/gamma, M_solidModel->timeStepNewmark()->previousAcceleration() );
        if ( true )
            fieldToTransfert->add( -1.0, M_solidModel->timeStepNewmark()->currentAcceleration() );
        else
        {
            if ( iterationFSI == 0 )
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i+1) );
            }
            else
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i) );
            }
        }
#endif

        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dconf->apply( *fieldToTransfert, *this->couplingRNG_evalForm1() );
        }
        else
        {
            CHECK( M_opVelocity2dTo2dnonconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dnonconf->apply( *fieldToTransfert, *this->couplingRNG_evalForm1() );
        }

        M_couplingRNG_coeffForm2 = (1./(dt*gamma));
    }
    else if ( M_solidModel->is1dReducedModel() )
    {
        typename solid_type::element_vect_1dreduced_ptrtype fieldToTransfert;
        auto fieldExtrapolated2 = this->solidModel()->timeStepNewmark1dReduced()->previousVelocity().functionSpace()->elementPtr();
        double dt = M_solidModel->timeStepNewmark1dReduced()->timeStep();
        double gamma = M_solidModel->timeStepNewmark1dReduced()->gamma();
        double beta = M_solidModel->timeStepNewmark1dReduced()->beta();
        double scaleTimeDisc = M_solidModel->mechanicalProperties()->cstRho()*M_solidModel->thickness1dReduced();

#if 0
        if ( useOriginalMethod || (iterationFSI == 0)  )
            fieldExtrapolated2->add( (1./(dt*gamma))*( (gamma/beta)-1. ) -1./(beta*dt) , M_solidModel->timeStepNewmark1dReduced()->previousVelocity() );
        else
        {
            //fieldExtrapolated2->add( (1./(dt*gamma))*( (gamma/beta)-1. ) -1./(beta*dt) , M_solidModel->timeStepNewmark1dReduced()->currentVelocity() );
        }


        if ( useOriginalMethod || (iterationFSI == 0)  )
        {
            fieldExtrapolated2->add( (1./gamma)*(gamma/(2*beta) - 1) - (1./(2*beta) -1), M_solidModel->timeStepNewmark1dReduced()->previousAcceleration() );
            fieldExtrapolated2->add( -1.0, M_solidModel->timeStepNewmark1dReduced()->currentAcceleration() );
        }
        else
        {
            //fieldExtrapolated2->add( (1./gamma)*(gamma/(2*beta) - 1) - (1./(2*beta) -1) + 1.0, M_solidModel->timeStepNewmark1dReduced()->currentAcceleration() );
            fieldExtrapolated2->add(-(1./(dt*gamma)), M_solidModel->timeStepNewmark1dReduced()->currentVelocity() );
        }
        fieldExtrapolated2->scale( scaleTimeDisc );
#else
        fieldExtrapolated2->add( -1./(dt*gamma) , M_solidModel->timeStepNewmark1dReduced()->previousVelocity() );
        fieldExtrapolated2->add( -(1.-gamma)/gamma, M_solidModel->timeStepNewmark1dReduced()->previousAcceleration() );
        if ( true )
            fieldExtrapolated2->add( -1.0, M_solidModel->timeStepNewmark1dReduced()->currentAcceleration() );
        else
        {
            if ( iterationFSI == 0 )
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         this->solidModel()->timeStepNewmark1dReduced()->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             this->solidModel()->timeStepNewmark1dReduced()->previousVelocity(i+1) );
            }
            else
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         this->solidModel()->timeStepNewmark1dReduced()->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             this->solidModel()->timeStepNewmark1dReduced()->previousVelocity(i) );
            }
        }
        fieldExtrapolated2->scale( scaleTimeDisc );
#endif
        fieldToTransfert = M_solidModel->extendVelocity1dReducedVectorial( *fieldExtrapolated2 );

        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity1dToNdconf ) << "interpolation operator not build";
            M_opVelocity1dToNdconf->apply(*fieldToTransfert,
                                          *this->couplingRNG_evalForm1() );
        }
        else
        {
            CHECK( M_opVelocity1dToNdnonconf ) << "interpolation operator not build";
            M_opVelocity1dToNdnonconf->apply(*fieldToTransfert,
                                             *this->couplingRNG_evalForm1() );
        }

        M_couplingRNG_coeffForm2 = scaleTimeDisc*(1./(dt*gamma));
    }

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "finish",
                                               this->worldComm(),this->verboseAllProc());

}
#else

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannGeneralizedS2F( int iterationFSI )
{
    if ( this->solidModel()->timeStepping() == "Newmark" )
        this->transfertRobinNeumannGeneralizedS2F_BdfNewmark( iterationFSI );
    else
        this->transfertRobinNeumannGeneralizedS2F_BdfBdf( iterationFSI );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannGeneralizedS2F_BdfBdf( int iterationFSI )
{
   if ( this->solidModel()->isStandardModel())
    {
        auto fieldToTransfert = M_solidModel->fieldVelocity().functionSpace()->elementPtr();
        double dt = this->solidModel()->timeStepBdfDisplacement()->timeStep();
        if ( iterationFSI == 0 )
        {
            fieldToTransfert->add( -2, this->solidModel()->timeStepBdfVelocity/*timeStepBdfDisplacement*/()->unknown(0) );
            fieldToTransfert->add( 1, this->solidModel()->timeStepBdfVelocity/*timeStepBdfDisplacement*/()->unknown(1) );
            fieldToTransfert->scale( 1./dt );
        }
        else
        {
            //fieldToTransfert->add( -2, this->solidModel()->fieldVelocity() );
            //fieldToTransfert->add( 1, this->solidModel()->timeStepBdfVelocity()->unknown(0) );
            fieldToTransfert->add( -1, this->solidModel()->fieldVelocity() );
            fieldToTransfert->scale( 1./dt );
            //for ( uint8_type i = 0; i < this->timeOrder(); ++i )
            //M_polyDeriv->add( this->polyDerivCoefficient( i+1 ), *M_unknowns[i] );
        }

        M_couplingRNG_coeffForm2 = 1./dt;

        // transfer solid to fluid
        auto fieldInterpolatedOnInterface = M_XhMeshVelocityInterface->element();
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dconf->apply( *fieldToTransfert, fieldInterpolatedOnInterface );
        }
        else
        {
            CHECK( M_opVelocity2dTo2dnonconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dnonconf->apply( *fieldToTransfert, fieldInterpolatedOnInterface );
        }
        this->couplingRNG_evalForm1()->on( _range=M_rangeFSI_fluid,_expr=idv(fieldInterpolatedOnInterface));
        sync( *this->couplingRNG_evalForm1(), "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );
    }
   else if ( this->solidModel()->is1dReducedModel() )
   {
       CHECK( false ) << "TODO";
   }
}


template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertRobinNeumannGeneralizedS2F_BdfNewmark( int iterationFSI )
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "start",
                                               this->worldComm(),this->verboseAllProc());

    bool useMeanBdfNewmark = (M_coulingRNG_strategyTimeStepCompatibility == "mean");
    if ( this->solidModel()->isStandardModel())
    {
        auto fieldToTransfert = M_solidModel->fieldVelocity().functionSpace()->elementPtr();
        double dt = M_solidModel->timeStepNewmark()->timeStep();
        double gamma = M_solidModel->timeStepNewmark()->gamma();
        double beta = M_solidModel->timeStepNewmark()->beta();
        // time derivative acceleration in solid (newmark)
        fieldToTransfert->add( -1, M_solidModel->timeStepNewmark()->currentAcceleration());
        // time derivative acceleration in fluid (newamrk)
        fieldToTransfert->add( -1./(dt*gamma) , this->solidModel()->timeStepNewmark()->previousVelocity() );
        fieldToTransfert->add( -(1.-gamma)/gamma, this->solidModel()->timeStepNewmark()->previousAcceleration() );

        if ( useMeanBdfNewmark )
        {
            // time derivative acceleration in solid (bdf)
            if ( iterationFSI == 0 )
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i+1) );
            }
            else
            {
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                       this->solidModel()->timeStepNewmark()->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldToTransfert->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                           this->solidModel()->timeStepNewmark()->previousVelocity(i) );
            }

            // time derivative acceleration in fluid (bdf)
            for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                fieldToTransfert->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                       this->solidModel()->timeStepNewmark()->previousVelocity(i) );

            fieldToTransfert->scale( 0.5 );
        }

        // transfer solid to fluid
        auto fieldInterpolatedOnInterface = M_XhMeshVelocityInterface->element();
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dconf->apply( *fieldToTransfert, fieldInterpolatedOnInterface );
        }
        else
        {
            CHECK( M_opVelocity2dTo2dnonconf ) << "interpolation operator not build";
            M_opVelocity2dTo2dnonconf->apply( *fieldToTransfert, fieldInterpolatedOnInterface );
        }
        this->couplingRNG_evalForm1()->on( _range=M_rangeFSI_fluid,_expr=idv(fieldInterpolatedOnInterface));
        sync( *this->couplingRNG_evalForm1(), "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );

        // time derivative acceleration in fluid
#if 0
        auto UPolyDeriv = this->fluidModel()->timeStepBDF()->polyDeriv();
        auto uPolyDeriv = UPolyDeriv.template element<0>();
        this->couplingRNG_evalForm1()->add( -1.0, uPolyDeriv );
        M_couplingRNG_coeffForm2 = this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0);
#else
        //M_couplingRNG_coeffForm2 = scaleTimeDisc*this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0);
        if ( useMeanBdfNewmark )
            M_couplingRNG_coeffForm2 = 0.5*(this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0) + (1./(dt*gamma)));
        else
            M_couplingRNG_coeffForm2 = (1./(dt*gamma));
#endif

    }
    else if ( this->solidModel()->is1dReducedModel() )
    {
        typename solid_type::solid_1dreduced_type::element_displacement_ptrtype fieldToTransfert;
        auto timeStepNewmark1dReduced = this->solidModel()->solid1dReduced()->timeStepNewmark();
        auto fieldExtrapolated2 = timeStepNewmark1dReduced->previousVelocity().functionSpace()->elementPtr();
        double dt = timeStepNewmark1dReduced->timeStep();
        double gamma = timeStepNewmark1dReduced->gamma();
        double beta = timeStepNewmark1dReduced->beta();

        //double rhoValue = 1.1;//0; // M_solidModel->mechanicalProperties()->cstRho()
        //CHECK( false ) << "TODO : rho";
        //double scaleTimeDisc = rhoValue*M_solidModel->solid1dReduced()->thickness1dReduced();
        double scaleTimeDisc = 1;
        // time derivative acceleration in solid (newmark)
        fieldExtrapolated2->add( -1, timeStepNewmark1dReduced->currentAcceleration());
        // time derivative acceleration in fluid (newamrk)
        fieldExtrapolated2->add( -1./(dt*gamma) , timeStepNewmark1dReduced->previousVelocity() );
        fieldExtrapolated2->add( -(1.-gamma)/gamma, timeStepNewmark1dReduced->previousAcceleration() );

        if ( useMeanBdfNewmark )
        {
            // time derivative acceleration in solid (bdf)
            if ( iterationFSI == 0 )
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         timeStepNewmark1dReduced->previousVelocity(0) );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             timeStepNewmark1dReduced->previousVelocity(i+1) );
            }
            else
            {
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( 0 ),
                                         timeStepNewmark1dReduced->currentVelocity() );
                for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                    fieldExtrapolated2->add( this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                             timeStepNewmark1dReduced->previousVelocity(i) );
            }

            // time derivative acceleration in fluid (bdf)
            for ( uint8_type i = 0; i < this->fluidModel()->timeStepBDF()->timeOrder(); ++i )
                fieldExtrapolated2->add( -this->fluidModel()->timeStepBDF()->polyDerivCoefficient( i+1 ),
                                         timeStepNewmark1dReduced->previousVelocity(i) );
            fieldExtrapolated2->scale( 0.5 );
        }

        fieldExtrapolated2->scale( scaleTimeDisc );

        fieldToTransfert = M_solidModel->solid1dReduced()->extendVelocity1dReducedVectorial( *fieldExtrapolated2 );
        // transfer solid to fluid
        auto fieldInterpolatedOnInterface = M_XhMeshVelocityInterface->element();
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity1dToNdconf ) << "interpolation operator not build";
            M_opVelocity1dToNdconf->apply(*fieldToTransfert,fieldInterpolatedOnInterface );
        }
        else
        {
            CHECK( M_opVelocity1dToNdnonconf ) << "interpolation operator not build";
            M_opVelocity1dToNdnonconf->apply(*fieldToTransfert,fieldInterpolatedOnInterface );
        }
        this->couplingRNG_evalForm1()->on( _range=M_rangeFSI_fluid,_expr=idv(fieldInterpolatedOnInterface));
        sync( *this->couplingRNG_evalForm1(), "=", M_dofsMultiProcessVelocitySpaceOnFSI_fluid );
        // time derivative acceleration in fluid
        //auto UPolyDeriv = this->fluidModel()->timeStepBDF()->polyDeriv();
        //auto uPolyDeriv = UPolyDeriv.template element<0>();
        //this->couplingRNG_evalForm1()->add( -scaleTimeDisc, uPolyDeriv );

        if ( useMeanBdfNewmark )
            M_couplingRNG_coeffForm2 = scaleTimeDisc*(0.5*this->fluidModel()->timeStepBDF()->polyDerivCoefficient(0)+0.5*(1./(dt*gamma)));
        else
            M_couplingRNG_coeffForm2 = scaleTimeDisc*(1./(dt*gamma));
    }

    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertRobinNeumannGeneralizedS2F", "finish",
                                               this->worldComm(),this->verboseAllProc());
}
#endif
//-----------------------------------------------------------------------------------//

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertStressS2F()
{
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStressS2F", "start",
                                               this->worldComm(),this->verboseAllProc());
#if 0
    this->solidModel()->updateNormalStressFromStruct();
    CHECK( M_opStress2dTo2dconfS2F ) << "interpolation operator not build";
    M_opStress2dTo2dconfS2F->apply(*this->solidModel()->fieldNormalStressFromStructPtr(),
                                   *this->fluidModel()->normalStressFromStruct());
#else
    CHECK( false ) << "not implemented";
#endif
    if (this->verbose()) Feel::FeelModels::Log("InterpolationFSI","transfertStressS2F", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertVelocityF2S( int iterationFSI, bool _useExtrapolation )
{
    bool useExtrapolation = ( iterationFSI == 0) && _useExtrapolation && (this->fluidModel()->timeStepBDF()->iteration() > 2);
    if ( useExtrapolation )
    {
        if ( true )
        {
            // bdf extrapolation
            auto const& velExtrap = this->fluidModel()->timeStepBDF()->poly();
            if (M_interfaceFSIisConforme)
            {
                CHECK( M_opVelocity2dTo2dconfF2S ) << "interpolation operator not build";
                M_opVelocity2dTo2dconfF2S->apply( velExtrap,*this->fieldVelocityInterfaceFromFluidPtr_solid() );
            }
            else
            {
                CHECK(false) << "TODO\n";
            }
        }
        else
        {
            // extrap in Explicit strategies for incompressible fluid-structure interaction problems: Nitche ....
#if 0
            auto velExtrap = this->fluidModel()->functionSpaceVelocity()->element();
            velExtrap.add(  2.0, this->fluidModel()->getSolution()->template element<0>() );
            velExtrap.add( -1.0, this->fluidModel()->timeStepBDF()->unknown(0).template element<0>() );
#else
            auto velExtrap = this->fluidModel()->functionSpaceVelocity()->element();
            velExtrap.add(  2.0, this->fluidModel()->timeStepBDF()->unknown(0) );
            velExtrap.add( -1.0, this->fluidModel()->timeStepBDF()->unknown(1) );
#endif
            if (M_interfaceFSIisConforme)
            {
                CHECK( M_opVelocity2dTo2dconfF2S ) << "interpolation operator not build";
                M_opVelocity2dTo2dconfF2S->apply( velExtrap,*this->fieldVelocityInterfaceFromFluidPtr_solid() );
            }
            else
            {
                CHECK(false) << "TODO\n";
            }
        }
    }
    else // no extrap : take current solution
    {
        if (M_interfaceFSIisConforme)
        {
            CHECK( M_opVelocity2dTo2dconfF2S ) << "interpolation operator not build";
#if 0
            M_opVelocity2dTo2dconfF2S->apply( this->fluidModel()->timeStepBDF()->unknown(0),
                                              *this->fieldVelocityInterfaceFromFluidPtr_solid() );
#else
            M_opVelocity2dTo2dconfF2S->apply( this->fluidModel()->fieldVelocity(),
                                              *this->fieldVelocityInterfaceFromFluidPtr_solid() );
#endif
        }
        else
        {
            CHECK(false) << "TODO\n";
        }
    }
}



template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::transfertGradVelocityF2S()
{
    this->log("FSI","transfertGradVelocityF2S","start");

    bool meshIsOnRefAtBegin = this->fluidModel()->meshMotionTool()->isOnReferenceMesh();
    if ( !meshIsOnRefAtBegin )
        this->fluidModel()->meshMotionTool()->revertReferenceMesh( false );

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();
    auto Fa = Id+gradv( this->fluidModel()->meshMotionTool()->displacement());
    auto fsirange = M_rangeFSI_fluid;
    for ( int k = 0;k<fluid_type::nDim ;++k )
    {
        Component comp = ( k==0 )? Component::X : ( k == 1 )? Component::Y : Component::Z;
        //auto uComp = this->fluidModel()->fieldVelocityPressurePtr()->template element<0>().comp(comp);
        auto const& uComp = this->fluidModel()->fieldVelocity().comp(comp);
        M_fieldsGradVelocity_fluid[k]->on(_range=fsirange,_expr=trans(inv(Fa))*trans(gradv(uComp)));
    }

    if ( !meshIsOnRefAtBegin )
        this->fluidModel()->meshMotionTool()->revertMovingMesh( false );

    if ( M_interfaceFSIisConforme )
    {
        CHECK( M_opStress2dTo2dconf ) << "interpolation operator not build";
        for ( int k = 0;k<fluid_type::nDim ;++k )
            M_opStress2dTo2dconf->apply( *M_fieldsGradVelocity_fluid[k], *M_fieldsGradVelocity_solid[k] );
    }
    else
    {
        CHECK( M_opStress2dTo2dnonconf ) << "interpolation operator not build";
        for ( int k = 0;k<fluid_type::nDim ;++k )
            M_opStress2dTo2dnonconf->apply( *M_fieldsGradVelocity_fluid[k], *M_fieldsGradVelocity_solid[k] );
    }

    this->log("FSI","transfertGradVelocityF2S","finish");
}



} // namespace FeelModels
} // namespace Feel
