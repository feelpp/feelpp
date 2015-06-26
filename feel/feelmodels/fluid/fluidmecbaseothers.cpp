/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels/fluid/fluidmecbase.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/stdmathfunctors.hpp>
#include <feel/feelvf/symm.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/val.hpp>
#include <feel/feelvf/pow.hpp>
#include <feel/feelvf/det.hpp>
#include <feel/feelvf/inv.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/geometricdata.hpp>
//#include <fsi/fsicore/variousfunctions.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::useExtendedDofTable() const
{
    if ( this->worldComm().localSize() == 1 ) return false;
    bool useExtendedDofTable=false;
    for ( bool hasExt : M_Xh->extendedDofTableComposite() )
        useExtendedDofTable = useExtendedDofTable || hasExt;
    return useExtendedDofTable;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::restartExporters()
{
    if (this->doRestart() && this->restartPath().empty() )
    {
        if (!M_isHOVisu)
        {
            // if restart and same directory, update the exporter for new value, else nothing (create a new exporter)
            if ( true )//nOrderGeo == 1 && this->application()->vm()["exporter.format"].as< std::string >() == "ensight")
            {
                if ( M_exporter->doExport() ) M_exporter->restart(this->timeInitial());
            }
            else
            {
                //if ( M_exporter_gmsh->doExport() ) M_exporter_gmsh->restart(this->timeInitial());
            }
        }
        else
        {
#if defined(FEELPP_HAS_VTK)
            if ( M_exporter_ho->doExport() ) M_exporter_ho->restart(this->timeInitial());
#endif
        }
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::string ALEmode;if (M_isMoveDomain) ALEmode="Yes"; else ALEmode="No";
    std::string StateTemporal;if (this->isStationary()) StateTemporal="Stationary"; else StateTemporal="Transient";

    std::string ResartMode;
    if (this->doRestart()) ResartMode = "Yes";
    else ResartMode = "No";

    std::string stabAll_str;
    if (this->doStabConvectionEnergy()) stabAll_str=(stabAll_str.empty())?"convection energy":stabAll_str+" - convection energy";
    if (this->doCIPStabConvection()) stabAll_str=(stabAll_str.empty())?"CIP Convection":stabAll_str+" - CIP Convection";
    if (this->doCIPStabDivergence()) stabAll_str=(stabAll_str.empty())?"CIP Divergence":stabAll_str+" - CIP Divergence";
    if (this->doCIPStabPressure()) stabAll_str=(stabAll_str.empty())?"CIP Pressure":stabAll_str+" - CIP Pressure";
    if (this->doStabDivDiv()) stabAll_str=(stabAll_str.empty())?"DivDiv":stabAll_str+" - DivDiv";
    //if (this->doCstPressureStab()) stabAll_str=(stabAll_str.empty())?"CstPressure":stabAll_str+" - CstPressure";
    if (stabAll_str.empty()) stabAll_str="OFF";

    std::string hovisuMode,myexporterType;
    int myexporterFreq = 1;
    if (M_isHOVisu)
    {
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        if ( hovisuSpaceUsed == "velocity" || hovisuSpaceUsed == "pressure" )
            hovisuMode = "ON (with OperatorLagrangeP1 on "+hovisuSpaceUsed+")";
        else if ( hovisuSpaceUsed == "p1" )
            hovisuMode = "ON (with createP1mesh)";
#if defined(FEELPP_HAS_VTK)
        myexporterType=M_exporter_ho->type();
        myexporterFreq = M_exporter_ho->freq();
#endif
    }
    else
    {
        hovisuMode = "OFF";
        myexporterType = M_exporter->type();
        myexporterFreq = M_exporter->freq();
    }
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : FluidMechanics---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Appli Repository : " << this->appliRepository()
           << "\n   Physical Model"
           << "\n     -- pde name  : " << M_pdeType
           << "\n     -- stress tensor law  : " << this->densityViscosityModel()->dynamicViscosityLaw()
           << "\n     -- time mode : " << StateTemporal
           << "\n     -- ale mode  : " << ALEmode
           << "\n   Physical Parameters"
           << "\n     -- rho : " << this->densityViscosityModel()->cstRho()
           << "\n     -- mu  : " << this->densityViscosityModel()->cstMu()
           << "\n     -- nu  : " << this->densityViscosityModel()->cstNu();
    *_ostr << this->densityViscosityModel()->getInfo()->str();
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoALEMeshBC();
    if ( M_fluidOutletsBCType.size()>0 )
    {
        for ( auto itBC = M_fluidOutletsBCType.begin(), enBC = M_fluidOutletsBCType.end() ; itBC!=enBC ; ++itBC )
        {
            *_ostr << "\n       -- fluid outlets (" << itBC->first << " with " << this->nFluidOutlet() << " outlets) : ";
            int cptMark = 0;
            for ( auto itMark = itBC->second.begin(), enMark = itBC->second.end() ; itMark!= enMark ; ++itMark, ++cptMark)
            {
                if ( cptMark > 0) *_ostr << " , ";
                *_ostr << *itMark;
            }
        }
    }
    *_ostr << "\n   Space Discretization"
           << "\n     -- msh file name   : " << this->mshfileStr()
           << "\n     -- nb elt in mesh  : " << M_mesh->numGlobalElements()//numElements()
        // << "\n     -- nb elt in mesh  : " << M_mesh->numElements()
        // << "\n     -- nb face in mesh : " << M_mesh->numFaces()
           << "\n     -- geometry order  : " << nOrderGeo
           << "\n     -- velocity order  : " << nOrderVelocity
           << "\n     -- pressure order  : " << nOrderPressure
           << "\n     -- nb dof (u+p)    : " << M_Xh->nDof() << " (" << M_Xh->nLocalDof() << ")"
           << "\n     -- nb dof (u)      : " << M_Xh->template functionSpace<0>()->nDof()
           << "\n     -- nb dof (p)      : " << M_Xh->template functionSpace<1>()->nDof()
           << "\n     -- stabilisation   : " << stabAll_str;
    if ( this->definePressureCst() )
    {
        *_ostr << "\n     -- define cst pressure  : " << this->definePressureCstMethod();
        if ( this->definePressureCstMethod() == "penalisation" )
            *_ostr << " ( beta=" << this->definePressureCstPenalisationBeta() << ")";
    }

    *_ostr << "\n   Time Discretization"
           << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
           << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
           << "\n     -- time step    : " << this->timeStepBase()->timeStep()
           << "\n     -- order        : " << this->timeStepBDF()->timeOrder()
           << "\n     -- restart mode : " << ResartMode
           << "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
    if ( this->timeStepBase()->saveFreq() )
        *_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
               << "\n     -- file format save : " << this->timeStepBase()->fileFormat();
    *_ostr << "\n   Exporter"
           << "\n     -- type            : " << myexporterType
           << "\n     -- high order visu : " << hovisuMode
           << "\n     -- freq save : " << myexporterFreq
           << "\n   Processors"
           << "\n     -- number of proc environnement : " << Environment::worldComm().globalSize()
           << "\n     -- environement rank : " << Environment::worldComm().rank()
           << "\n     -- global rank : " << this->worldComm().globalRank()
           << "\n     -- local rank : " << this->worldComm().localRank()
        //<< "\n   Matrix Index Start"
        //  << "\n     -- rowstart : " << this->rowStartInMatrix()
        //  << "\n     -- colstart : " << this->colStartInMatrix()
           << "\n   Numerical Solver"
           << "\n     -- solver : " << M_pdeSolver;
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
        *_ostr << this->getMeshALE()->getInfo()->str();
#endif
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

//---------------------------------------------------------------------------------------------------------//
#if 0
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::printInfo() const
{
    if ( this->verboseAllProc() ) std::cout << this->getInfo()->str();
    else if (this->worldComm().isMasterRank() )
        std::cout << this->getInfo()->str();
}

//---------------------------------------------------------------------------------------------------------//

void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::saveInfo() const
{
    std::string nameFile = prefixvm(this->prefix(),"FluidMechanics.info");
    this->saveInfo( nameFile );
    /*if (this->worldComm().isMasterRank() )
     {
     std::string nameFile = prefixvm(this->prefix(),"FluidMechanics.info");
     std::ofstream file(nameFile.c_str(), std::ios::out);
     file << this->getInfo()->str();
     file.close();
     }*/
}
#endif
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeType(std::string __type)
{
    // if pde change -> force to rebuild all algebraic data at next solve
    if ( __type != M_pdeType )
        this->setNeedToRebuildCstPart(true);

    if ( __type == "Stokes" )
    {
        M_pdeType="Stokes";
        M_pdeSolver="LinearSystem";
    }
    else if ( __type == "Oseen" )
    {
        M_pdeType="Oseen";
        M_pdeSolver="LinearSystem";
    }
    else if ( __type == "Navier-Stokes" )
    {
        M_pdeType="Navier-Stokes";
        M_pdeSolver="Newton";
    }
    else
        CHECK( false ) << "invalid pdeType "<< __type << "\n";
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeType() const
{
    return M_pdeType;
}
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeSolver(std::string __type)
{
    // if solver change -> force to rebuild all algebraic data at next solve
    if ( __type != M_pdeSolver )
        this->setNeedToRebuildCstPart(true);

    if ( __type == "LinearSystem" )
        M_pdeSolver="LinearSystem";
    else if ( __type == "PtFixe" )
        M_pdeSolver="PtFixe";
    else if ( __type == "Newton" )
        M_pdeSolver="Newton";
    else
        CHECK( false ) << "invalid pdeSolver " << __type << "\n";
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeSolver() const
{
    return M_pdeSolver;
}

//---------------------------------------------------------------------------------------------------------//
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::stressTensorLawType(std::string __type)
{
    // if viscosity model change -> force to rebuild all algebraic data at next solve
    if ( __type != this->densityViscosityModel()->dynamicViscosityLaw() )
        this->setNeedToRebuildCstPart(true);

    this->densityViscosityModel()->setDynamicViscosityLaw( __type );
}

//---------------------------------------------------------------------------------------------------------//
#if 0
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::stressTensorLawType() const
{
    return M_stressTensorLaw;
}
#endif
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::setDoExport(bool b)
{
    if (!M_isHOVisu)
    {
        CHECK( M_exporter )  << "export not init\n";
        M_exporter->setDoExport(b);
    }
    else
    {
#if defined(FEELPP_HAS_VTK)
        CHECK( M_exporter )  << "hoexport not init\n";
        M_exporter_ho->setDoExport(b);
#endif
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("FluidMechanics","exportResults", (boost::format("start at time %1%")%time).str() );
    this->timerTool("PostProcessing").start();

    if (this->isMoveDomain() && (M_doExportMeshALE || M_doExportAll ) )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        this->getMeshALE()->exportResults( time );
#endif
    }

    if ( false )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        //ExporterGeometry::EXPORTER_GEOMETRY_STATIC
        std::string geoExportType="static";//change_coords_only, change, static

        if ( !M_exporterFluidOutlet )
            M_exporterFluidOutlet = exporter( _mesh=M_fluidOutletWindkesselMesh,
                                              _name="ExportFluidOutet",
                                              _geo=geoExportType,
                                              _worldcomm=M_Xh->worldComm(),
                                              _path=this->exporterPath() );

        M_exporterFluidOutlet->step( time )->add( prefixvm(this->prefix(),"fluidoutlet-disp"),
                                                  prefixvm(this->prefix(),prefixvm(this->subPrefix(),"fluidoutlet-disp")),
                                                  *M_fluidOutletWindkesselMeshDisp );
        M_exporterFluidOutlet->save();
#endif
    }

    if (!M_isHOVisu)
    {
        this->exportResultsImpl( time );
    }
    else
    {
        this->exportResultsImplHO( time );
    }

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("FluidMechanics","exportResults", "finish" );

} // FluidMechanics::exportResult


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    //if (this->worldComm().globalSize()==1) this->updateVorticity(mpl::int_<nDim>());

    //if ( true )//nOrderGeo == 1 && this->application()->vm()["exporter.format"].as< std::string >() == "ensight")
    //{
    if ( !M_exporter->doExport() ) return;


#if defined( FEELPP_MODELS_HAS_MESHALE )
    // because write geofile at each step ( TODO fix !!! )
    if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
        this->getMeshALE()->revertReferenceMesh();
#endif
    //M_exporter->step( time )->setMesh( M_mesh );
    M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity")),
                                   M_Solution->template element<0>() );
    this->log("FluidMechanics","exportResults", "velocity done" );

    M_exporter->step( time )->add( prefixvm(this->prefix(),"pressure"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure")),
                                   M_Solution->template element<1>() );
    if (M_doExportVorticity /*|| M_doExportAll*/)
    {
        this->updateVorticity(mpl::int_<nDim>());
        M_exporter->step( time )->add( prefixvm(this->prefix(),"vorticity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"vorticity")),
                                       *M_vorticity );
    }
    if (M_doExportNormalStress || M_doExportAll)
    {
        this->updateNormalStressOnCurrentMesh();
        M_exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")),
                                       this->fieldNormalStress() );
    }
    if (M_doExportWallShearStress || M_doExportAll)
    {
        this->updateWallShearStress();
        M_exporter->step( time )->add( prefixvm(this->prefix(),"wallshearstress"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")),
                                       this->fieldWallShearStress() );
    }
    if ( M_doExportViscosity || M_doExportAll )
    {
        if ( !M_XhNormalBoundaryStress ) this->createFunctionSpacesNormalStress();
        auto uCur = M_Solution->template element<0>();
        auto pCur = M_Solution->template element<1>();
        auto myViscosity = Feel::vf::FeelModels::fluidMecViscosity<2*nOrderVelocity>(uCur,pCur,*this->densityViscosityModel());
        auto viscosityField = M_XhNormalBoundaryStress->compSpace()->element(myViscosity);
        M_exporter->step( time )->add( prefixvm(this->prefix(),"viscosity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"viscosity")),
                                       viscosityField );
    }
    if ( this->isMoveDomain() )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        auto drm = M_meshALE->dofRelationShipMap();
        auto thedisp = M_meshALE->functionSpace()->element();

        for (size_type i=0;i<thedisp.nLocalDof();++i)
            thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i);

        M_exporter->step( time )->add( prefixvm(this->prefix(),"displacement"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacement")),
                                       thedisp );

        if (M_doExportMeshDisplacementOnInterface || M_doExportAll)
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"displacementOnInterface"),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacementOnInterface")),
                                           this->meshDisplacementOnInterface() );
        }
#endif
    }

    //----------------------//
    M_exporter->save();
    this->log("FluidMechanics","exportResults", "save done" );

    //----------------------//
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
        this->getMeshALE()->revertMovingMesh();
#endif
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResultsImplHO( double time )
{
#if defined(FEELPP_HAS_VTK)
    if ( !M_exporter_ho->doExport() ) return;

    // because write geofile at each step ( TODO fix !!! )
    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //    this->getMeshALE()->revertReferenceMesh();

    M_opIvelocity->apply(M_Solution->template element<0>(),*M_velocityVisuHO);
    M_opIpressure->apply(M_Solution->template element<1>(),*M_pressureVisuHO);
    //M_exporter_ho->step( time )->setMesh( M_velocityVisuHO->mesh() );
    M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity_ho")), *M_velocityVisuHO );
    M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure_ho")), *M_pressureVisuHO );

    if (M_isMoveDomain)
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        auto drm = M_meshALE->dofRelationShipMap();
        auto thedisp = M_meshALE->functionSpace()->element();
        for (size_type i=0;i<thedisp.nLocalDof();++i)
        {
#if 0
            thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i) - (*M_meshALE->dispP1ToHO_ref())(i) ;
            //thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i);
#else
            thedisp(i) = M_meshALE->displacement()->operator()(i);
#endif
        }

        M_opImeshdisp->apply( thedisp , *M_meshdispVisuHO);

        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshdisp_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshdisp_ho")), *M_meshdispVisuHO );

        if ( M_doExportAll )
        {
            auto themeshvelocityOnInterfaceVisuHO = M_XhVectorialVisuHO->element();
            auto hola2 = vf::project(_space=this->functionSpaceVelocity(),_range=elements(M_Xh->mesh()),_expr=vf::idv(this->meshVelocity2Ptr()));
            //auto hola2 = vf::project(_space=M_Xh->functionSpace<0>(),_range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
            //                         _expr=vf::idv(this->meshVelocity2Ptr()));
            M_opIvelocity->apply( hola2/* *this->meshVelocity2Ptr()*/,themeshvelocityOnInterfaceVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshVelocityOnInterface_ho"),
                                              prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshVelocityOnInterface_ho")),
                                              themeshvelocityOnInterfaceVisuHO );
            auto themeshvelocityOnVolumeHO = M_XhVectorialVisuHO->element();
            auto hola = vf::project(_space=this->functionSpaceVelocity(),_range=elements(M_Xh->mesh()),_expr=vf::idv(this->meshVelocity()));
            M_opIvelocity->apply( hola/* *this->meshVelocity2Ptr()*/,themeshvelocityOnVolumeHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshVelocityOnVolume_ho"),
                                              prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshVelocityOnVolume_ho")),
                                              themeshvelocityOnVolumeHO );
            //--------------------------------------------------------------//
            auto thedisplacementInRef = M_meshALE->functionSpace()->element();
            for (size_type i=0;i<thedisplacementInRef.nLocalDof();++i)
                thedisplacementInRef(drm->dofRelMap()[i])=M_meshALE->displacementInRef()->operator()(i) - M_meshALE->dispP1ToHO_ref()->operator()(i) ;
            //thedisplacementInRef(drm->dofRelMap()[i])=M_meshALE->aleFactory()->displacement().operator()(i) - M_meshALE->dispP1ToHO_ref()->operator()(i) ;
            auto thedisplacementInRefVisuHO = M_XhVectorialVisuHO->element();
            M_opImeshdisp->apply( thedisplacementInRef, thedisplacementInRefVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshale_dispinref_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshale_dispinref_ho")),
                                              thedisplacementInRefVisuHO );
            //--------------------------------------------------------------//
            auto thedispP1ToHO = M_meshALE->functionSpace()->element();
            for (size_type i=0;i<thedispP1ToHO.nLocalDof();++i)
                thedispP1ToHO(drm->dofRelMap()[i])=M_meshALE->dispP1ToHO_ref()->operator()(i);
            auto thedispP1ToHOVisuHO = M_XhVectorialVisuHO->element();
            M_opImeshdisp->apply( thedispP1ToHO , thedispP1ToHOVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshdisp_P1toHO_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshdisp_P1toHO_ho")),
                                              thedispP1ToHOVisuHO );
            //--------------------------------------------------------------//
            auto theidentityALEVisuHO = M_XhVectorialVisuHO->element();
            M_opImeshdisp->apply( *M_meshALE->identityALE() , theidentityALEVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshale_identity_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshale_identity_ho")),
                                              theidentityALEVisuHO );
            //--------------------------------------------------------------//
            auto thedisplacementOnMovingBoundaryVisuHO = M_XhVectorialVisuHO->element();
            M_opImeshdisp->apply( *M_meshALE->displacementOnMovingBoundary(),thedisplacementOnMovingBoundaryVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshale_disponboundary_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshale_disponboundary_ho")),
                                              thedisplacementOnMovingBoundaryVisuHO );
            //--------------------------------------------------------------//
            auto thedisplacementOnMovingBoundaryInRef = M_meshALE->functionSpace()->element();
            for (size_type i=0;i<thedispP1ToHO.nLocalDof();++i)
                thedisplacementOnMovingBoundaryInRef(drm->dofRelMap()[i])=M_meshALE->displacementOnMovingBoundaryInRef()->operator()(i) - M_meshALE->dispP1ToHO_ref()->operator()(i) ;
            auto thedisplacementOnMovingBoundaryInRefVisuHO = M_XhVectorialVisuHO->element();
            M_opImeshdisp->apply( thedisplacementOnMovingBoundaryInRef,thedisplacementOnMovingBoundaryInRefVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshale_disponboundaryinref_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshale_disponboundaryinref_ho")),
                                              thedisplacementOnMovingBoundaryInRefVisuHO );
            //--------------------------------------------------------------//
        }
#endif // HAS_MESHALE
    }

    if (M_doExportNormalStress)
    {
        this->updateNormalStressOnCurrentMesh();
        M_opIstress->apply( this->fieldNormalStress(),*M_normalStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"normalstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")), *M_normalStressVisuHO );
    }
    if (M_doExportWallShearStress)
    {
        this->updateWallShearStress();
        M_opIstress->apply( this->fieldWallShearStress(),*M_fieldWallShearStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"wallshearstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")), *M_fieldWallShearStressVisuHO );
    }

    M_exporter_ho->save();

    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //   this->getMeshALE()->revertMovingMesh();

#endif
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("FluidMechanics","solve", "start" );
    this->timerTool("Solve").start();

    // copy velocity/pressure in algebraic vector solution (maybe velocity/pressure has been changed externaly)
    this->updateBlockVectorSolution();

    if ( this->startBySolveStokesStationary() && !this->isStationary() &&
         !this->hasSolveStokesStationaryAtKickOff() && !this->doRestart() )
    {
        this->log("FluidMechanics","solve", "start by solve stokes stationary" );

        std::string saveStressTensorLawType = this->densityViscosityModel()->dynamicViscosityLaw();//stressTensorLawType();
        std::string savePdeType = this->pdeType();
        // prepare Stokes-stationary config
        this->stressTensorLawType( "newtonian" );
        this->pdeType( "Stokes" );
        this->setStationary( true );
        // possibility to config a specific time which appear in bc
        double timeUsed = doption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary.time-value-used-in-bc");
        this->updateTime( timeUsed );

        this->solve();
        this->hasSolveStokesStationaryAtKickOff( true );

        if ( boption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary.do-export") )
            this->exportResults(this->timeInitial());
        // revert parameters
        this->stressTensorLawType( saveStressTensorLawType );
        this->pdeType( savePdeType );
        this->setStationary( false );

        this->initTimeStep();
    }

    if ( this->startBySolveNewtonian() && this->densityViscosityModel()->dynamicViscosityLaw() != "newtonian" &&
         !this->hasSolveNewtonianAtKickOff() && !this->doRestart() )
    {
        this->log("FluidMechanics","solve", "start by solve newtonian" );

        std::string saveStressTensorLawType = this->densityViscosityModel()->dynamicViscosityLaw();//stressTensorLawType();
        this->stressTensorLawType("newtonian");
        this->solve();
        this->hasSolveNewtonianAtKickOff( true );
        this->stressTensorLawType( saveStressTensorLawType );
    }


    //--------------------------------------------------
    // solution vector
    //auto Uvec = this->backend()->newBlockVector(_block=M_blockVectorSolution/*M_blockMonolitic*/);
    //--------------------------------------------------
    // run solver
    if (M_pdeSolver=="LinearSystem")
    {
        M_algebraicFactory->linearSolver(this->blockVectorSolution().vector()/*Uvec*/);
    }
    else if ( M_pdeSolver == "PtFixe")
    {
        M_algebraicFactory->AlgoPtFixe(this->blockVectorSolution().vector()/*Uvec*/);
    }
    else if ( M_pdeSolver == "Newton")
    {
        //Uvec->close(); //????????
        M_algebraicFactory->AlgoNewton2(this->blockVectorSolution().vector()/*Uvec*/);
    }

    //Uvec->close();

    //--------------------------------------------------
    //get solution compute with the numerical method
    //*M_Solution = *Uvec;
    //M_blockVectorSolution.localize(Uvec);
    M_blockVectorSolution.localize();

    // update windkessel solution ( todo put fluidOutletWindkesselPressureDistal and Proximal in blockVectorSolution )
    int cptBlock=1;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        ++cptBlock;
    if (this->hasMarkerDirichletBClm())
        ++cptBlock;
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        for (int k=0;k<M_nFluidOutlet;++k)
        {
            if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
            {
                M_fluidOutletWindkesselPressureDistal[k] = M_blockVectorSolution(cptBlock)->operator()(0);
                M_fluidOutletWindkesselPressureProximal[k]= M_blockVectorSolution(cptBlock)->operator()(1);
            }
            ++cptBlock;
        }
    }

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }
    this->log("FluidMechanics","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( sparse_matrix_ptrtype const& mat,
                                                                     vector_ptrtype const& vecSol ) const
{
    if ( this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasInHousePreconditioners( "blockns" ) )
    {
        this->updateInHousePreconditionerPCD( mat,vecSol );
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateLinearPDE(const vector_ptrtype& X,sparse_matrix_ptrtype& A , vector_ptrtype& F,bool _buildCstPart,
                                                        sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                                        bool _doClose, bool _doBCStrongDirichlet) const
{
    if ( M_pdeType == "Stokes" || M_pdeType == "Oseen" )
        updateOseen(A,F,_buildCstPart,A_extended,_BuildExtendedPart,_doClose,_doBCStrongDirichlet);
    else if ( M_pdeType == "Navier-Stokes")
        updatePtFixe(X,A,F,_buildCstPart,_doClose,_doBCStrongDirichlet);
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateVorticity(mpl::int_<2> /***/)
{

    if (!M_Xh_vorticity) this->createFunctionSpacesVorticity();

    //BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,2),
    *M_vorticity = vf::project(_space=M_vorticity->functionSpace(),
                               _range=elements(M_vorticity->mesh()),
                               _expr=vf::abs(vf::dxv(M_Solution->template element<0>().template comp<ComponentType::Y>())
                                             -vf::dyv(M_Solution->template element<0>().template comp<ComponentType::X>())),
                               _geomap=this->geomap() );
    //std::cout << "\nWARNING UPDATE VORTICITY\n"
    //);
}
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateVorticity(mpl::int_<3> /***/)
{
#if 0
    if (!M_Xh_vorticity) this->createFunctionSpacesVorticity();

    BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,3),
                *M_vorticity = vf::project(_space=M_vorticity->functionSpace(),
                                           _range=elements(M_vorticity->mesh()),
                                           _expr=vf::curlv(M_Solution->element<0>()),
                                           _geomap=this->geomap() ),
                std::cout << "\nWARNING UPDATE VORTICITY\n"
                );
#endif
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateBdf()
{
    this->log("FluidMechanics","updateBdf", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // windkessel outlet
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel")
    {
        bool doWriteOnDisk = this->worldComm().isMasterRank();
        // warning, in the implicit case, not all process have info
        // only process which have a face on outlet
        if (this->fluidOutletWindkesselCoupling() == "implicit" )
            doWriteOnDisk = (bool)(M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0);

        for (int k=0;k<M_nFluidOutlet;++k)
        {
            if ( doWriteOnDisk )
            {
                std::string nameFile = this->appliRepository() + "/" + prefixvm(this->prefix(),(boost::format("bloodFlowOutlet.windkessel%1%.data") %k ).str());
                std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
                file << M_bdf_fluid->iteration() << " " << this->time() << " " << M_fluidOutletWindkesselPressureDistal[k] << " " << M_fluidOutletWindkesselPressureProximal[k] << "\n";
                file.close();
            }

            const int sizeOld = M_fluidOutletWindkesselPressureDistal_old[k].size();
            for (int l=0;l<sizeOld-1;++l)
                M_fluidOutletWindkesselPressureDistal_old[k][sizeOld-l-1] = M_fluidOutletWindkesselPressureDistal_old[k][sizeOld-l-2];
            M_fluidOutletWindkesselPressureDistal_old[k][0] = M_fluidOutletWindkesselPressureDistal[k];
        }
    }


    int previousTimeOrder = this->timeStepBDF()->timeOrder();

    M_bdf_fluid->next( *M_Solution );

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain())
        M_meshALE->updateBdf();
#endif
    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf_fluid->time() );

    // maybe rebuild cst jacobian or linear
    if ( M_algebraicFactory &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (this->pdeSolver() == "Newton" && !this->rebuildLinearPartInJacobian() )
        {
            this->log("FluidMechanics","updateBdf", "do rebuildCstJacobian" );
            M_algebraicFactory->rebuildCstJacobian(M_Solution);
        }
        else if (this->pdeSolver() == "LinearSystem" && !this->rebuildCstPartInLinearSystem())
        {
            this->log("FluidMechanics","updateBdf", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(M_Solution);
        }
    }

    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("FluidMechanics","updateBdf", "finish" );
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressOnCurrentMesh( std::list<std::string> const& listMarkers )
{
    this->log("FluidMechanics","updateNormalStressOnCurrentMesh", "start" );

    // current solution
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    if( !M_fieldNormalStress )
        this->createFunctionSpacesNormalStress();

    auto const Id = eye<nDim,nDim>();
    // deformations tensor
    auto defv = sym(gradv(u));
    // stress tensor
    auto Sigmav = -idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv;

    M_fieldNormalStress->zero();
    if ( listMarkers.empty() )
        M_fieldNormalStress->on(_range=boundaryfaces(this->mesh()),
                                       _expr=Sigmav*N(),
                                       _geomap=this->geomap() );
    else
        M_fieldNormalStress->on(_range=markedfaces(this->mesh(),listMarkers),
                                _expr=Sigmav*N(),
                                _geomap=this->geomap() );

    this->log("FluidMechanics","updateNormalStressOnCurrentMesh", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressOnReferenceMesh( std::list<std::string> const& listMarkers )
{
    if (this->useFSISemiImplicitScheme())
        this->updateNormalStressOnReferenceMeshOptSI( listMarkers );
    else
        this->updateNormalStressOnReferenceMeshStandard( listMarkers );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressOnReferenceMeshStandard( std::list<std::string> const& listMarkers )
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    using namespace Feel::vf;

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshStandard", "start" );

    // current solution
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();

    //Tenseur des deformations
    auto defv = sym(gradv(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Tenseur des contraintes (trial)
    auto Sigmav = -idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv;

    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());

#if 0
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = Fa(0,0);
    auto Fa12 = Fa(0,1);
    auto Fa21 = Fa(1,0);
    auto Fa22 = Fa(1,1);
    //sans le determinant devant car il s annule avec un terme apres
    auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = Fa(0,0);
    auto Fa12 = Fa(0,1);
    auto Fa13 = Fa(0,2);
    auto Fa21 = Fa(1,0);
    auto Fa22 = Fa(1,1);
    auto Fa23 = Fa(1,2);
    auto Fa31 = Fa(2,0);
    auto Fa32 = Fa(2,1);
    auto Fa33 = Fa(2,2);
    //sans le determinant devant car il s annule avec un terme apres
    auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
                           Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
                           Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
                           );
#endif
#else
    auto InvFa = det(Fa)*inv(Fa);
#endif
    this->getMeshALE()->revertReferenceMesh();

    M_fieldNormalStressRefMesh->zero();
    if ( listMarkers.empty() )
        M_fieldNormalStressRefMesh->on(_range=boundaryfaces(this->mesh()),
                                       _expr=val(Sigmav*trans(InvFa)*N()),
                                       _geomap=this->geomap() );
    else
        M_fieldNormalStressRefMesh->on(_range=markedfaces(this->mesh(),listMarkers/*this->markersNameMovingBoundary()*/),
                                       _expr=val(Sigmav*trans(InvFa)*N()),
                                       _geomap=this->geomap() );

    this->getMeshALE()->revertMovingMesh();

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshStandard", "finish" );
#endif
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressOnReferenceMeshOptPrecompute( std::list<std::string> const& listMarkers )
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    using namespace Feel::vf;

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshOptPrecompute", "start" );

    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());
#if 0
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = Fa(0,0);
    auto Fa12 = Fa(0,1);
    auto Fa21 = Fa(1,0);
    auto Fa22 = Fa(1,1);
    //sans le determinant devant car il s annule avec un terme apres
    auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = Fa(0,0);
    auto Fa12 = Fa(0,1);
    auto Fa13 = Fa(0,2);
    auto Fa21 = Fa(1,0);
    auto Fa22 = Fa(1,1);
    auto Fa23 = Fa(1,2);
    auto Fa31 = Fa(2,0);
    auto Fa32 = Fa(2,1);
    auto Fa33 = Fa(2,2);
    //sans le determinant devant car il s annule avec un terme apres
    auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
                           Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
                           Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
                           );
#endif
#else
    auto InvFa = det(Fa)*inv(Fa);
#endif



    if (!M_saveALEPartNormalStress) M_saveALEPartNormalStress = M_XhMeshALEmapDisc->elementPtr();

    this->getMeshALE()->revertReferenceMesh();


    M_saveALEPartNormalStress->zero();
    if ( listMarkers.empty() )
        M_saveALEPartNormalStress->on(_range=boundaryfaces(this->mesh()),
                                      _expr=val(trans(InvFa)*N()),
                                      _geomap=this->geomap() );
    else
        M_saveALEPartNormalStress->on(_range=markedfaces(this->mesh(),listMarkers/*this->markersNameMovingBoundary()*/),
                                      _expr=val(trans(InvFa)*N()),
                                      _geomap=this->geomap() );
#if 0
    *M_saveALEPartNormalStress = vf::project(_space=M_XhMeshALEmapDisc,
                                             _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                             _expr=val(trans(InvFa)*N()),
                                             _geomap=this->geomap() );
#endif

    this->getMeshALE()->revertMovingMesh();

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshOptPrecompute", "finish" );
#endif
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressOnReferenceMeshOptSI( std::list<std::string> const& listMarkers )
{
#if defined( FEELPP_MODELS_HAS_MESHALE )
    using namespace Feel::vf;

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshOptSI", "start" );

    // current solution
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();

    //Tenseur des deformations
    auto defv = sym(gradv(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Tenseur des contraintes (trial)
    auto Sigmav = -idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv;

    this->getMeshALE()->revertReferenceMesh();

    if ( M_saveALEPartNormalStress )
    {
        if ( listMarkers.empty() )
            M_fieldNormalStressRefMesh->on(_range=boundaryfaces(this->mesh()),
                                           _expr=val(Sigmav*idv(M_saveALEPartNormalStress)),
                                           _geomap=this->geomap() );
        else
            M_fieldNormalStressRefMesh->on(_range=markedfaces(this->mesh(),listMarkers/*this->markersNameMovingBoundary()*/),
                                           _expr=val(Sigmav*idv(M_saveALEPartNormalStress)),
                                           _geomap=this->geomap() );
    }
    else
    {
        if ( listMarkers.empty() )
            M_fieldNormalStressRefMesh->on(_range=boundaryfaces(this->mesh()),
                                      _expr=Sigmav*N(),
                                      _geomap=this->geomap() );
        else
            M_fieldNormalStressRefMesh->on(_range=markedfaces(this->mesh(),listMarkers/*this->markersNameMovingBoundary()*/),
                                           _expr=Sigmav*N(),
                                           _geomap=this->geomap() );
    }

    this->getMeshALE()->revertMovingMesh();

    this->log("FluidMechanics","updateNormalStressOnReferenceMeshOptSI", "finish" );
#endif
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateWallShearStress()
{
    using namespace Feel::vf;

    this->log("FluidMechanics","updateWallShearStress", "start" );

    if ( !M_fieldWallShearStress ) this->createFunctionSpacesNormalStress();

    // current solution
    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();
    auto p = solFluid->template element<1>();

    //Tenseur des deformations
    auto defv = sym(gradv(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Tenseur des contraintes (trial)
    auto Sigmav = (-idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv);

    M_fieldWallShearStress->on(_range=boundaryfaces(this->mesh()),
                               _expr=Sigmav*vf::N() - (trans(Sigmav*vf::N())*vf::N())*vf::N(),
                               _geomap=this->geomap() );

    this->log("FluidMechanics","updateWallShearStress", "finish" );
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeDiff(const vector_ptrtype& X1,const vector_ptrtype& X2)
{
    using namespace Feel::vf;

    auto U1 = M_Xh->element();U1 = *X1;
    auto u1 = U1.template element<0>();
    auto p1 = U1.template element<1>();

    auto U2 = M_Xh->element();U2 = *X2;
    auto u2 = U2.template element<0>();
    auto p2 = U2.template element<1>();

    double err_u = integrate(_range=elements(M_Xh->mesh()),
                             _expr= trans(idv(u1)-idv(u2))*(idv(u1)-idv(u2)),
                             _geomap=this->geomap() ).evaluate()(0,0);
    double err_p = integrate(_range=elements(M_Xh->mesh()),
                             _expr= (idv(p1)-idv(p2))*(idv(p1)-idv(p2)),
                             _geomap=this->geomap() ).evaluate()(0,0);

    return std::sqrt(err_u+err_p);
}

//---------------------------------------------------------------------------------------------------------//

#if defined( FEELPP_MODELS_HAS_MESHALE )

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateALEmesh()
{
    this->log("FluidMechanics","updateALEmesh", "start");

    //-------------------------------------------------------------------//
    // compute ALE map
    //std::vector< mesh_ale_type::ale_map_element_type> polyBoundarySet = { *M_meshDisplacementOnInterface };
    M_meshALE->update(*M_meshDisplacementOnInterface/*polyBoundarySet*/);

    //-------------------------------------------------------------------//
    // up mesh velocity on interface from mesh velocity
#if 1
    // no mpi comm
    /**M_meshVelocityInterface = vf::project( _space=M_meshVelocityInterface->functionSpace(),
     _range=markedfaces(this->mesh(),this->markersNameMovingBoundary().front()),
     _expr=vf::idv(M_meshALE->velocity()),
     _geomap=this->geomap() );*/
    M_meshVelocityInterface->on(//_range=boundaryelements(this->mesh()),
        _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
        //_range=boundaryfaces(this->mesh()),
        //_range=elements(this->mesh()),
        _expr=vf::idv(M_meshALE->velocity()),
        _geomap=this->geomap() );
#else
    // with mpi comm
    auto vectVelInterface = backend()->newVector( M_meshVelocityInterface->functionSpace() );
    modifVec(markedfaces(this->mesh(),this->markersNameMovingBoundary()),
             *M_meshVelocityInterface,
             vectVelInterface,
             vf::idv(M_meshALE->velocity()) );

    vectVelInterface->close();
    *M_meshVelocityInterface = *vectVelInterface;
#endif

    //-------------------------------------------------------------------//
    // semi implicit optimisation
    if (this->useFSISemiImplicitScheme())
    {
        //this->getMeshALE()->revertReferenceMesh();
        this->updateNormalStressOnReferenceMeshOptPrecompute(this->markersNameMovingBoundary());
        //this->getMeshALE()->revertMovingMesh();
    }

    //-------------------------------------------------------------------//
    // move winkessel submesh
    if ( this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        // revert on reference mesh
        M_fluidOutletWindkesselMeshDisp->scale(-1);
        M_fluidOutletWindkesselMeshMover.apply( M_fluidOutletWindkesselMesh, *M_fluidOutletWindkesselMeshDisp );
        // interpolate disp
        M_fluidOutletWindkesselOpMeshDisp->apply( *M_meshALE->displacement(), *M_fluidOutletWindkesselMeshDisp );
        // apply disp
        M_fluidOutletWindkesselMeshMover.apply( M_fluidOutletWindkesselMesh, *M_fluidOutletWindkesselMeshDisp );

        if (this->verbose())
        {
            double normWind = M_fluidOutletWindkesselMeshDisp->l2Norm();
            double normDisp = M_meshALE->displacement()->l2Norm();
            double normDispImposed = M_meshDisplacementOnInterface->l2Norm();
            this->log("FluidMechanics","updateALEmesh",
                      (boost::format( "normWind %1% normDisp %2% normDispImposed %3% ") %normWind %normDisp %normDispImposed ).str() );
        }
    }

    //-------------------------------------------------------------------//
    //ho visu
#if defined(FEELPP_HAS_VTK)
    if (M_isHOVisu && false) // useless now ( this export mesh stay fix!)
    {
        auto drm = M_meshALE->dofRelationShipMap();
        //revert ref
        M_meshdispVisuHO->scale(-1);
        M_meshmover_visu_ho.apply(M_XhVectorialVisuHO->mesh(),*M_meshdispVisuHO);
        // get disp from ref
        auto thedisp = M_meshALE->functionSpace()->element();
        for (uint i=0;i<thedisp.nLocalDof();++i)
            thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i) - (*M_meshALE->dispP1ToHO_ref())(i) ;
        //transfert disp on ho mesh
        M_opImeshdisp->apply( thedisp , *M_meshdispVisuHO);
        // move ho mesh
        M_meshmover_visu_ho.apply(M_XhVectorialVisuHO->mesh(),*M_meshdispVisuHO);
    }
#endif


    this->log("FluidMechanics","updateALEmesh", "finish");
}
#endif

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FluidMechanicsBase< ConvexType,BasisVelocityType,
                    BasisPressureType,BasisDVType,UsePeriodicity >::savePressureAtPoints(const std::list<boost::tuple<std::string,typename mesh_type::node_type> > & __listPt,
                                                                                         bool extrapolate)
{
    auto it = __listPt.begin();
    auto en = __listPt.end();

    for ( ; it!=en ; ++it)
    {
        auto nameFile = prefixvm(this->prefix(),"eval-pressure-")+boost::get<0>(*it)+".data";
        if (M_nameFilesPressureAtPoints.find(nameFile)==M_nameFilesPressureAtPoints.end())
        {
            if (!this->doRestart() && this->worldComm().globalRank()==this->worldComm().masterRank())
            {
                //ATTENTION EFFACER LES DONNEES REDONTANTE EN CAS DE RESTART
                std::ofstream newfilePres(nameFile.c_str(), std::ios::out | std::ios::trunc);
                newfilePres.close();
                M_nameFilesPressureAtPoints.insert(nameFile);
            }
        }
        auto const pressureOnPt = this->fieldVelocityPressurePtr()->template element<1>()(boost::get<1>(*it),extrapolate);

        if (this->worldComm().globalRank()==this->worldComm().masterRank())
        {
            std::ofstream filePres(nameFile.c_str(), std::ios::out | std::ios::app);
            filePres << this->time();
            //for (uint16_type i = 0; i<mesh_type::nRealDim;++i)
            uint16_type c = 0; //comp
            filePres << " " << pressureOnPt(c,0,0);
            filePres << "\n";
            filePres.close();
        }
    }

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
Eigen::Matrix<typename FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::value_type,
              FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nDim,1>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeForce(std::string markerName)
{
    using namespace Feel::vf;

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();
    auto p = solFluid->template element<1>();
    //Tenseur des deformations
    auto defv = sym(gradv(u));
    // Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Tenseur des contraintes
    auto Sigmav = val(-idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv);

    return integrate(_range=markedfaces(M_mesh,markerName),
                     _expr= Sigmav*N(),
                     _geomap=this->geomap() ).evaluate();

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeFlowRate(std::string marker)
{
    using namespace Feel::vf;

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();

    double res = integrate(_range=markedfaces(this->mesh(),marker),
                           _expr= -trans(idv(u))*N(),
                           _geomap=this->geomap() ).evaluate()(0,0);

    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMeanPressure()
{
    using namespace Feel::vf;


    auto solFluid = this->fieldVelocityPressurePtr();
    auto p = solFluid->template element<1>();

    double area = integrate(_range=elements(this->mesh()),
                            _expr=cst(1.) ).evaluate()(0,0);
    double res = (1./area)*integrate(_range=elements(this->mesh()),
                                     _expr= idv(p),
                                     _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMeanDivergence()
{
    using namespace Feel::vf;

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();

    double area = integrate(_range=elements(this->mesh()),
                            _expr=cst(1.) ).evaluate()(0,0);
    double res = (1./area)*integrate(_range=elements(this->mesh()),
                                     _expr= divv(u),
                                     _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeNormL2Divergence()
{
    using namespace Feel::vf;

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();

    double res = math::sqrt(integrate(_range=elements(this->mesh()),
                                      _expr= pow(divv(u),2),
                                      _geomap=this->geomap() ).evaluate()(0,0));
    return res;
}

//---------------------------------------------------------------------------------------------------------//

#if defined( FEELPP_MODELS_HAS_MESHALE )

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::vector<double>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeAveragedPreassure( std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                                                  std::vector<op_interp_pressure_ptrtype> const& opInterp,
                                                                  bool computeOnRefMesh,
                                                                  std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp )
{
    int nbSlice = setMeshSlices.size();
    std::vector<double> res(nbSlice);

    auto solFluid = this->fieldVelocityPressurePtr();
    auto p = solFluid->template element<1>();

    bool computeOnRefMesh2 = this->isMoveDomain() && computeOnRefMesh;
    if ( !computeOnRefMesh2 )
    {
        for (uint16_type i = 0 ; i<nbSlice ; ++i)
        {
            auto meshSlice = setMeshSlices[i];
            double area = integrate(_range=elements(meshSlice),
                                    _expr=cst(1.) ).evaluate()(0,0);
            if ( opInterp[i] )
            {
                auto pInterp = opInterp[i]->operator()( p );
                res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                              _expr=idv(pInterp) ).evaluate()(0,0));
            }
            else
            {
                res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                              _expr=idv(p) ).evaluate()(0,0));
            }
        }
    }
    else
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        this->getMeshALE()->revertReferenceMesh();
        auto const Id = eye<nDim,nDim>();
        for (uint16_type i = 0 ; i<nbSlice ; ++i)
        {
            auto meshSlice = setMeshSlices[i];
            double area = integrate(_range=elements(meshSlice),
                                    _expr=cst(1.) ).evaluate()(0,0);
            if ( opInterp[i] )
            {
                auto pInterp = opInterp[i]->operator()( p );
                auto dispInterp = opInterpMeshDisp[i]->operator()( *M_meshALE->displacement() );
                // Deformation tensor
                //double holl=integrate(_range=elements(meshSlice),_expr=inner(gradv(dispInterp),Id2) ).evaluate()(0,0);
                //double holl=integrate(_range=elements(meshSlice),_expr=inner(gradv(*M_meshALE->displacement()),Id) ).evaluate()(0,0);
                auto Fa = Id+gradv(dispInterp);
                auto detFa = det(Fa);

                res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                              _expr=idv(pInterp)*detFa ).evaluate()(0,0));
            }
            else
            {
                // Deformation tensor
                auto Fa = Id+gradv(*M_meshALE->displacement());
                auto detFa = det(Fa);
                res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                              _expr=idv(p)*detFa ).evaluate()(0,0));
            }
        }
        this->getMeshALE()->revertMovingMesh();
#endif
    }


    return res;

}
FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::vector<double>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeFlowRate(std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                                        std::vector<op_interp_velocity_ptrtype> const& opInterp,
                                                        bool computeOnRefMesh,
                                                        std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp )
{
    int nbSlice = setMeshSlices.size();
    std::vector<double> res(nbSlice);
    std::vector<double> res2(nbSlice);

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();
    auto dirVelocity = oneX();//vec(cst(1.0),cst(0.));

    bool computeOnRefMesh2 = this->isMoveDomain() && computeOnRefMesh;
    //if ( !computeOnRefMesh2 )
    {
        for (uint16_type i = 0 ; i<nbSlice ; ++i)
        {
            auto meshSlice = setMeshSlices[i];
            res2[i] = integrate(_range=elements(meshSlice),
                                _expr=inner(idv(u),dirVelocity) ).evaluate()(0,0);
        }
    }
    //else
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        this->getMeshALE()->revertReferenceMesh();

        auto const Id = eye<nDim,nDim>();
        for (uint16_type i = 0 ; i<nbSlice ; ++i)
        {
            auto meshSlice = setMeshSlices[i];
            if ( opInterp[i] )
            {
                auto uInterp = opInterp[i]->operator()( u );
                auto dispInterp = opInterpMeshDisp[i]->operator()( *M_meshALE->displacement() );
                //auto Fa = Id+gradv(dispInterp);
                auto Fa = Id+gradv(*M_meshALE->displacement());
                auto detFa = det(Fa);
                res[i] = integrate(_range=elements(meshSlice),
                                   _expr=inner(idv(uInterp),dirVelocity)*detFa ).evaluate()(0,0);
            }
            else
            {
                // Deformation tensor
                auto Fa = Id+gradv(*M_meshALE->displacement());
                auto detFa = det(Fa);
                res[i] = integrate(_range=elements(meshSlice),
                                   _expr=inner(idv(u),dirVelocity)*detFa ).evaluate()(0,0);
            }
        }

        this->getMeshALE()->revertMovingMesh();
#endif
    }

    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        std::cout <<"res[i] " << res[i] <<  " res2[i] " <<  res2[i] << " jac " << std::abs(res[i] - res2[i] ) << "\n";
        CHECK( std::abs(res[i] - res2[i] ) < 1e-5 ) << "res[i] " << res[i] <<  " res2[i] " <<  res2[i] << " jac " << std::abs(res[i] - res2[i] ) << "\n";
    }

    return res;
}

#endif // HAS_MESHALE

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::block_pattern_type
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::blockPattern() const
{
    size_type pat_uu = size_type(Pattern::COUPLED);
    size_type pat_pp = size_type(Pattern::ZERO);

    bool doStabPSPG = boption(_prefix=this->prefix(),_name="stabilisation-pspg");
    bool doStabGLS = boption(_prefix=this->prefix(),_name="stabilisation-gls");
    if ( doStabPSPG || doStabGLS ||
         (this->definePressureCst() && this->definePressureCstMethod() == "penalisation") ||
         (this->isMoveDomain() && ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "nitsche" ) )
         )
        pat_pp = size_type(Pattern::COUPLED);

    if ( (this->doCIPStabConvection() || this->doCIPStabDivergence() ) && !this->applyCIPStabOnlyOnBoundaryFaces() )
        pat_uu = size_type(Pattern::EXTENDED);
    if ( this->doCIPStabPressure() )
        pat_pp = size_type(Pattern::EXTENDED);

    return vf::Blocks<2,2,size_type>() << pat_uu                       << size_type(Pattern::COUPLED)
                                       << size_type(Pattern::COUPLED)  << pat_pp;
} // blockPattern

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
int
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        ++nBlock;
    if (this->hasMarkerDirichletBClm())
        ++nBlock;
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
        nBlock += this->nFluidOutlet();
    return nBlock;
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("FluidMechanics","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;
    // first field (velocity+pressure)
    auto ru = stencilRange<0,0>(marked3faces(this->mesh(),1));
    auto rp = stencilRange<1,1>(marked3faces(this->mesh(),1));
    myblockGraph(indexBlock,indexBlock) =stencil(_test=this->functionSpace(),_trial=this->functionSpace(),
                                                 _pattern_block=this->blockPattern(),
                                                 //_diag_is_nonzero=false,_close=false,
                                                 _diag_is_nonzero=(nBlock==1),_close=(nBlock==1),
                                                 _range_extended=stencilRangeMap(ru,rp) )->graph();
    ++indexBlock;

    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        BlocksStencilPattern patCouplingLM(1,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingLM(0,1) = size_type(Pattern::COUPLED);

        myblockGraph(indexBlock,0) = stencil(_test=M_XhMeanPressureLM,_trial=this->functionSpace(),
                                             _pattern_block=patCouplingLM,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=this->functionSpace(),_trial=M_XhMeanPressureLM,
                                             _pattern_block=patCouplingLM.transpose(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    if (this->hasMarkerDirichletBClm())
    {
        BlocksStencilPattern patCouplingLM(1,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingLM(0,0) = size_type(Pattern::COUPLED);

        myblockGraph(indexBlock,0) = stencil(_test=this->XhDirichletLM(),_trial=this->functionSpace(),
                                             _pattern_block=patCouplingLM,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=this->functionSpace(),_trial=this->XhDirichletLM(),
                                             _pattern_block=patCouplingLM.transpose(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        BlocksStencilPattern patCouplingFirstCol(2,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingFirstCol(0,0) = size_type(Pattern::COUPLED);
        patCouplingFirstCol(1,0) = size_type(Pattern::COUPLED);

        BlocksStencilPattern patCouplingFirstRow(2,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingFirstRow(0,1) = size_type(Pattern::COUPLED);

        BlocksStencilPattern patCouplingDiag(2,2,size_type(Pattern::COUPLED));
        patCouplingDiag(0,1) = size_type(Pattern::ZERO);

        for (int k=0;k<this->nFluidOutlet();++k)
        {
            // first column
            auto rangeFirstCol1 = stencilRange<0,0>( markedelements(this->fluidOutletWindkesselMesh(),this->fluidOutletMarkerName(k)) );
            auto rangeFirstCol2 = stencilRange<1,0>( markedelements(this->fluidOutletWindkesselMesh(),this->fluidOutletMarkerName(k)) );
            myblockGraph(indexBlock+k,0) = stencil(_test=M_fluidOutletWindkesselSpace,_trial=this->functionSpace(),
                                                   _pattern_block=patCouplingFirstCol,
                                                   _range=stencilRangeMap(rangeFirstCol1,rangeFirstCol2),
                                                   _diag_is_nonzero=false,_close=false)->graph();

            // first line
#if 0
            myblockGraph(0,indexBlock+k) = stencil(_test=this->functionSpace(),_trial=M_fluidOutletWindkesselSpace,
                                                   _pattern_block=patCouplingFirstRow,
                                                   _diag_is_nonzero=false,_close=false)->graph();
#else
            myblockGraph(0,indexBlock+k) = myblockGraph(indexBlock+k,0)->transpose();
#endif

            // on diag
            myblockGraph(indexBlock+k,indexBlock+k) = stencil(_test=M_fluidOutletWindkesselSpace,_trial=M_fluidOutletWindkesselSpace,
                                                              //_pattern=(size_type)Pattern::COUPLED,
                                                              _pattern_block=patCouplingDiag,
                                                              _diag_is_nonzero=false,_close=false)->graph();
        }
        indexBlock += this->nFluidOutlet();
    }

    myblockGraph.close();
    this->log("FluidMechanics","buildBlockMatrixGraph", "finish" );

    return myblockGraph;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::graph_ptrtype
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    auto blockGraph = this->buildBlockMatrixGraph();
    if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
        return blockGraph(0,0);
    else
        return graph_ptrtype( new graph_type( blockGraph ) );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::indexsplit_ptrtype
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildIndexSplit() const
{
    // WARNING : this function is not use now

    //indexsplit_type theres(space_fluid_type::nSpaces+1);

    indexsplit_ptrtype res = this->functionSpace()->dofIndexSplit();

    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {
        int procMasterWindkessel=0;bool findProc=false;
        for (int proc=0;proc<this->worldComm().globalSize() && !findProc;++proc)
            if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhostOnProc(proc) > 0 )
            {
                procMasterWindkessel=proc;
                findProc=true;
            }
        CHECK(findProc)<< "not find the proc\n";

        if (this->worldComm().globalRank() > procMasterWindkessel )
        {
            for ( int k=0;k<=1;++k)
            {
                const size_type splitSize = (*res)[k].size();
                for ( size_type i=0; i<splitSize; ++i)
                    (*res)[k][i]+=2*this->nFluidOutlet();
            }
        }

        if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
        {
            const size_type velocitySplitSize = (*res)[0].size();
            (*res)[0].resize( velocitySplitSize + 2*this->nFluidOutlet(), 0 );
            size_type cpt = velocitySplitSize;
            const size_type startSplitW = (*res)[0][0]+this->functionSpace()->nLocalDofWithoutGhost();
            for (int k=0;k<this->nFluidOutlet();++k)
            {
                (*res)[0][cpt]=startSplitW+2*k;// res[0][cpt-1]+1;
                (*res)[0][cpt+1]=startSplitW+2*k+1;//res[0][cpt]+1;
                cpt+=2;
            }
        }
    }

    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->functionSpace()->nLocalDofWithGhost();
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        res += M_XhMeanPressureLM->nLocalDofWithGhost();
    if (this->hasMarkerDirichletBClm())
        res += this->XhDirichletLM()->nLocalDofWithGhost();
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
        res += 2*this->nFluidOutlet();
    return res;
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateBlockVectorSolution()
{
    // copy velocity/pressure in block
    auto & vecAlgebraic = M_blockVectorSolution.vector();
    auto const& fieldVelPres = this->fieldVelocityPressure();
    for (int k=0;k< this->functionSpace()->nLocalDofWithGhost() ;++k)
        vecAlgebraic->set(k, fieldVelPres(k) );

    // do nothing for others block (fields define only in blockVectorSolution)
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {}
    if (this->hasMarkerDirichletBClm())
    {}
    if (this->hasFluidOutlet() && this->fluidOutletType()=="windkessel" && this->fluidOutletWindkesselCoupling() == "implicit" )
    {}

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::couplingFSI_RNG_updateForUse()
{
    if ( M_couplingFSI_RNG_matrix ) return;

    this->log("FluidMechanics","couplingFSI_RNG_updateForUse", "start" );

    //-----------------------------------------------------------------//
    // build matrix
    //-----------------------------------------------------------------//
#if 0 // not work if we use a different pattern matrix with matrix
    auto ru = stencilRange<0,0>(markedfaces(this->mesh(),this->markersNameMovingBoundary()));
    auto myblockpattern = vf::Blocks<2,2,size_type>() << size_type(Pattern::COUPLED) << size_type(Pattern::ZERO)
                                                      << size_type(Pattern::ZERO)  << size_type(Pattern::ZERO);
    auto mygraph = stencil(_test=this->functionSpace(),_trial=this->functionSpace(),
                           _pattern_block=myblockpattern,
                           //_diag_is_nonzero=false,_close=false,
                           _diag_is_nonzero=false,_close=true,
                           _range=stencilRangeMap(ru) )->graph();
    M_couplingFSI_RNG_matrix = this->backend()->newMatrix(0,0,0,0,mygraph);
#else
    M_couplingFSI_RNG_matrix = this->backend()->newMatrix(0,0,0,0,this->algebraicFactory()->sparsityMatrixGraph());
#endif


    //-----------------------------------------------------------------//
    // assembly : first version
    //-----------------------------------------------------------------//
#if 0 /////////////////////////
    this->getMeshALE()->revertReferenceMesh();
    form2( _test=Xh,_trial=Xh,_matrix=M_couplingFSI_RNG_matrix,
           _pattern=size_type(Pattern::COUPLED) ) +=
        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                   _expr=this->couplingFSI_RNG_coeffForm2()*inner(idt(u),id(u)),
                   _geomap=this->geomap() );
    M_couplingFSI_RNG_matrix->close();
    this->getMeshALE()->revertMovingMesh();


    // scale with diagonal operator
    auto myoperator = this->couplingFSI_RNG_interfaceOperator();
    CHECK( /*M_couplingFSI_RNG_matrix->mapRow()*/this->functionSpaceVelocity()->nLocalDofWithGhost() == myoperator->map().nLocalDofWithGhost() ) << "invalid compatibility size";

    std::set<size_type> dofMarkerFsi,dofMarkerFsiP1;
#if 0 // bug!!!!
    for ( std::string const markName : this->markersNameMovingBoundary() )
    {
        //functionSpace()->dof()->faceLocalToGlobal( __face_it->id(), l, c1 )
        auto setofdof = this->functionSpaceVelocity()->dof()->markerToDof( markName );
        for( auto it = setofdof.first, en = setofdof.second; it != en; ++ it )
            dofMarkerFsi.insert( it->second );
    }
#endif
    for ( auto const& faceMarked : markedfaces(this->mesh(),this->markersNameMovingBoundary() ) )
    {
        auto __face_it = faceMarked.template get<1>();
        auto __face_en = faceMarked.template get<2>();
        for( ; __face_it != __face_en; ++__face_it )
            //for ( auto const& face : faceMarked )
        {
            auto const& face = *__face_it;
            for ( uint16_type l = 0; l < mesh_type::face_type::numVertices; ++l )
            {
                for (uint16_type c1=0;c1<nDim;c1++)
                {
                    size_type gdof = boost::get<0>(this->functionSpaceVelocity()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    dofMarkerFsiP1.insert( gdof );
                }
            }
            for ( uint16_type l = 0; l < this->functionSpaceVelocity()->dof()->nLocalDofOnFace(true); ++l )
            {
                for (uint16_type c1=0;c1<nDim;c1++)
                {
                    size_type gdof = boost::get<0>(this->functionSpaceVelocity()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    dofMarkerFsi.insert( gdof );
                }
            }

        }
    }
    //std::cout << "nLocalDofOnFace(true)" << this->functionSpaceVelocity()->dof()->nLocalDofOnFace(true) << "\n";
    //std::cout << "nLocalDofOnFace(false)" << this->functionSpaceVelocity()->dof()->nLocalDofOnFace(false) << "\n";
    //std::cout << "mesh_type::element_type::numVertices" << mesh_type::face_type::numVertices << "\n";
    std::cout << "size dofMarkerFsi " << dofMarkerFsi.size() << " size dofMarkerFsiP1 " << dofMarkerFsiP1.size() <<"\n";

    //for ( size_type k = 0 ; k<myoperator->map().nLocalDofWithGhost() ; ++k )
    for ( size_type k : dofMarkerFsi )
    {
        if ( myoperator->map().dofGlobalProcessIsGhost( k ) ) continue;
        size_type gcdof = myoperator->map().mapGlobalProcessToGlobalCluster(k);
#if 0
        double scalDiag = myoperator->operator()(k);
        double valDiag = M_couplingFSI_RNG_matrix->operator()( gcdof,gcdof);
        M_couplingFSI_RNG_matrix->set(k,k,valDiag*scalDiag);
#else
        auto const& graphRow = M_couplingFSI_RNG_matrix->graph()->row(gcdof);
        for ( size_type idCol : graphRow.template get<2>() )
        {
            if ( dofMarkerFsiP1.find(idCol) == dofMarkerFsiP1.end() ) continue; // Only P1 dof
            double scalDiag = myoperator->operator()(idCol);
            double valEntryMat = M_couplingFSI_RNG_matrix->operator()( gcdof,idCol);
            M_couplingFSI_RNG_matrix->set(k,idCol,valEntryMat*scalDiag);
        }
#endif
    }
    M_couplingFSI_RNG_matrix->close();

#endif ///////////


    //-----------------------------------------------------------------//
    // assembly : second version
    //-----------------------------------------------------------------//
    auto myoperator = this->couplingFSI_RNG_interfaceOperator();
    //M_couplingFSI_RNG_matrix->zero();
#if 0
    // CELUI QUI MARCHAIT A PEU PRES
    for ( size_type k : dofMarkerFsi )
        //for ( size_type k : dofMarkerFsiP1 )
    {
        double scalDiag = myoperator->operator()(k);
        M_couplingFSI_RNG_matrix->set(k,k,this->couplingFSI_RNG_coeffForm2()*scalDiag);
    }
#else
    std::vector<bool> dofdone(this->functionSpaceVelocity()->nLocalDofWithGhost(),false);
    for ( auto const& faceMarked : markedfaces(this->mesh(),this->markersNameMovingBoundary() ) )
    {
        auto __face_it = faceMarked.template get<1>();
        auto __face_en = faceMarked.template get<2>();
        for( ; __face_it != __face_en; ++__face_it )
        {
            auto const& face = *__face_it;
            for ( uint16_type l = 0; l < mesh_type::face_type::numVertices; ++l ) // only P1
                //for ( uint16_type l = 0; l < this->functionSpaceVelocity()->dof()->nLocalDofOnFace(true); ++l )
            {
                for (uint16_type c1=0;c1<nDim;c1++)
                {
                    size_type gdofVelFluid = boost::get<0>(this->functionSpaceVelocity()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    if ( dofdone[gdofVelFluid] ) continue;
                    double scalDiag = myoperator->operator()(gdofVelFluid);
                    M_couplingFSI_RNG_matrix->add(gdofVelFluid,gdofVelFluid,this->couplingFSI_RNG_coeffForm2()*scalDiag);
                    dofdone[gdofVelFluid] = true;
                }
            }
        }
    }
#endif
    M_couplingFSI_RNG_matrix->close();

    this->log("FluidMechanics","couplingFSI_RNG_updateForUse", "finish" );

}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::couplingFSI_RNG_updateLinearPDE( vector_ptrtype& F) const
{

#if 0
    this->getMeshALE()->revertReferenceMesh();
    form1( _test=Xh, _vector=F,
           _rowstart=rowStartInVector ) +=
        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                   _expr= -inner(idv(this->couplingFSI_RNG_evalForm1()),id(u)),
                   _geomap=this->geomap() );
#elif 0
    // this one has almost worked!!!!!!!!!!!
    // on dirait que le dofdone ralentit conv
    auto myoperator = this->couplingFSI_RNG_interfaceOperator();
    std::vector<bool> dofdone(this->functionSpaceVelocity()->nLocalDofWithGhost(),false);
    for ( auto const& faceMarked : markedfaces(this->mesh(),this->markersNameMovingBoundary() ) )
    {
        auto __face_it = faceMarked.template get<1>();
        auto __face_en = faceMarked.template get<2>();
        for( ; __face_it != __face_en; ++__face_it )
        {
            auto const& face = *__face_it;
            for ( uint16_type l = 0; l < mesh_type::face_type::numVertices; ++l ) // only P1
            {
                for (uint16_type c1=0;c1<nDim;c1++)
                {
                    size_type gdofVelFluid = boost::get<0>(this->functionSpaceVelocity()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    //if ( dofdone[gdofVelFluid] ) continue;
                    size_type gdofVelInterf = boost::get<0>(this->couplingFSI_RNG_evalForm1()->functionSpace()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    double thevalue = this->couplingFSI_RNG_evalForm1()->operator()(gdofVelInterf);//*(*myoperator)(gdofVelFluid);
                    F->add( gdofVelFluid, -thevalue );
                    //dofdone[gdofVelFluid] = true;
                }
            }
        }
    }
#else
    auto myoperator = this->couplingFSI_RNG_interfaceOperator();
    std::vector<bool> dofdone(this->functionSpaceVelocity()->nLocalDofWithGhost(),false);
    for ( auto const& faceMarked : markedfaces(this->mesh(),this->markersNameMovingBoundary() ) )
    {
        auto __face_it = faceMarked.template get<1>();
        auto __face_en = faceMarked.template get<2>();
        for( ; __face_it != __face_en; ++__face_it )
        {
            auto const& face = *__face_it;
            //for ( uint16_type l = 0; l < this->functionSpaceVelocity()->dof()->nLocalDofOnFace(true); ++l )
            for ( uint16_type l = 0; l < mesh_type::face_type::numVertices; ++l ) // only P1
            {
                for (uint16_type c1=0;c1<nDim;c1++)
                {
                    size_type gdofVelFluid = boost::get<0>(this->functionSpaceVelocity()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    if ( dofdone[gdofVelFluid] ) continue;
                    //double thevalue = this->couplingFSI_RNG_evalForm1Bis()->operator()(gdofVelFluid);//*(*myoperator)(gdofVelFluid);
                    double thevalue = this->couplingFSI_RNG_evalForm1Bis()->operator()(gdofVelFluid)*(*myoperator)(gdofVelFluid);

                    //size_type gdofVelInterf = boost::get<0>(this->couplingFSI_RNG_evalForm1()->functionSpace()->dof()->faceLocalToGlobal(face.id(), l, c1 ));
                    //double thevalue = this->couplingFSI_RNG_evalForm1()->operator()(gdofVelInterf);//*(*myoperator)(gdofVelFluid);

                    F->add( gdofVelFluid, -thevalue );
                    dofdone[gdofVelFluid] = true;
                }
            }
        }
    }
#endif

}






//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updatePtFixe(const vector_ptrtype& Xold, sparse_matrix_ptrtype& A , vector_ptrtype& F,
                                                     bool _buildCstPart,
                                                     bool _doClose, bool _doBCStrongDirichlet ) const
{}






} // namespace FeelModels
} // namespace Feel


