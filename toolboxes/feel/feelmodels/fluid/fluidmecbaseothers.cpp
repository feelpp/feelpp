/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

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
#include <feel/feelvf/mean.hpp>
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
boost::shared_ptr<std::ostringstream>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::string ALEmode;if (M_isMoveDomain) ALEmode="Yes"; else ALEmode="No";
    std::string StateTemporal;if (this->isStationary()) StateTemporal="Stationary"; else StateTemporal="Transient";

    std::string ResartMode;
    if (this->doRestart()) ResartMode = "Yes";
    else ResartMode = "No";

    std::string stabAll_str;
    if (M_stabilizationGLS) stabAll_str=(stabAll_str.empty())? this->stabilizationGLSType():stabAll_str+" - "+this->stabilizationGLSType();
    if (this->doStabConvectionEnergy()) stabAll_str=(stabAll_str.empty())?"convection energy":stabAll_str+" - convection energy";
    if (this->doCIPStabConvection()) stabAll_str=(stabAll_str.empty())?"CIP Convection":stabAll_str+" - CIP Convection";
    if (this->doCIPStabDivergence()) stabAll_str=(stabAll_str.empty())?"CIP Divergence":stabAll_str+" - CIP Divergence";
    if (this->doCIPStabPressure()) stabAll_str=(stabAll_str.empty())?"CIP Pressure":stabAll_str+" - CIP Pressure";
    if (this->doStabDivDiv()) stabAll_str=(stabAll_str.empty())?"DivDiv":stabAll_str+" - DivDiv";
    //if (this->doCstPressureStab()) stabAll_str=(stabAll_str.empty())?"CstPressure":stabAll_str+" - CstPressure";
    if (stabAll_str.empty()) stabAll_str="OFF";

    std::string hovisuMode,myexporterType;
    int myexporterFreq = 1;
    if ( M_exporter_ho )
    {
        std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
        if ( hovisuSpaceUsed == "velocity" || hovisuSpaceUsed == "pressure" )
            hovisuMode = "ON (with OperatorLagrangeP1 on "+hovisuSpaceUsed+")";
        else if ( hovisuSpaceUsed == "p1" )
            hovisuMode = "ON (with createP1mesh)";
#if 1//defined(FEELPP_HAS_VTK)
        myexporterType = M_exporter_ho->type();
        myexporterFreq = M_exporter_ho->freq();
#endif
    }
    else
    {
        hovisuMode = "OFF";
    }
    if ( M_exporter )
    {
        myexporterType = M_exporter->type();
        myexporterFreq = M_exporter->freq();
    }

    std::string doExport_str;
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Velocity ) )
        doExport_str=(doExport_str.empty())?"velocity":doExport_str+" - velocity";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pressure ) )
        doExport_str=(doExport_str.empty())?"pressure":doExport_str+" - pressure";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Displacement ) )
        doExport_str=(doExport_str.empty())?"displacement":doExport_str+" - displacement";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Vorticity ) )
        doExport_str=(doExport_str.empty())?"vorticity":doExport_str+" - vorticity";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) )
        doExport_str=(doExport_str.empty())?"normal stress":doExport_str+" - normal stress";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
        doExport_str=(doExport_str.empty())?"wall shear stress":doExport_str+" - wall shear stress";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Density ) )
        doExport_str=(doExport_str.empty())?"density":doExport_str+" - density";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Viscosity ) )
        doExport_str=(doExport_str.empty())?"viscosity":doExport_str+" - viscosity";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pid ) )
        doExport_str=(doExport_str.empty())?"pid":doExport_str+" - pid";
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::ALEMesh ) )
        doExport_str=(doExport_str.empty())?"alemesh":doExport_str+" - alemesh";
    for ( std::string const& userFieldName : M_postProcessUserFieldExported )
        doExport_str=(doExport_str.empty())?userFieldName:doExport_str+" - "+userFieldName;

    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : FluidMechanics---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository()
           << "\n   Physical Model"
           << "\n     -- pde name  : " << M_modelName
           << "\n     -- stress tensor law  : " << this->densityViscosityModel()->dynamicViscosityLaw()
           << "\n     -- time mode : " << StateTemporal
           << "\n     -- ale mode  : " << ALEmode
           << "\n     -- gravity  : " << std::boolalpha << M_useGravityForce;
    *_ostr << this->densityViscosityModel()->getInfoMaterialParameters()->str();
    // *_ostr << "\n   Physical Parameters"
    //        << "\n     -- rho : " << this->densityViscosityModel()->cstRho()
    //        << "\n     -- mu  : " << this->densityViscosityModel()->cstMu()
    //        << "\n     -- nu  : " << this->densityViscosityModel()->cstNu();
    // *_ostr << this->densityViscosityModel()->getInfo()->str();
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoPressureBC();
    for ( std::string typeOutlet : std::vector<std::string>({"free","windkessel"}) )
    {
        if ( this->hasFluidOutlet(typeOutlet) )
        {
            *_ostr << "\n       -- fluid outlets ["<<typeOutlet<<"] : ";
            bool hasDoneFirstElt = false;
            for (auto const& outletbc : M_fluidOutletsBCType )
            {
                if ( hasDoneFirstElt) *_ostr << " , ";
                *_ostr << std::get<0>( outletbc );
                if ( typeOutlet == "windkessel" )
                {
                    auto const& windkesselParam = std::get<2>( outletbc );
                    *_ostr << " (" << std::get<0>(windkesselParam)
                           << ",Rd=" << std::get<1>(windkesselParam)
                           << ",Rp=" << std::get<2>(windkesselParam)
                           << ",Cd=" << std::get<3>(windkesselParam) << ")";
                }
                hasDoneFirstElt=true;
            }
        }
    }
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
    *_ostr << this->getInfoALEMeshBC();
#endif
    *_ostr << "\n   Space Discretization";
    if ( this->hasGeofileStr() )
        *_ostr << "\n     -- geo file name   : " << this->geofileStr();
    *_ostr << "\n     -- mesh file name   : " << this->mshfileStr()
           << "\n     -- nb elt in mesh  : " << M_mesh->numGlobalElements()//numElements()
        // << "\n     -- nb elt in mesh  : " << M_mesh->numElements()
        // << "\n     -- nb face in mesh : " << M_mesh->numFaces()
           << "\n     -- hMin            : " << M_mesh->hMin()
           << "\n     -- hMax            : " << M_mesh->hMax()
           << "\n     -- hAverage        : " << M_mesh->hAverage()
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
           << "\n     -- freq save       : " << myexporterFreq
           << "\n     -- fields exported : " << doExport_str
           << "\n   Processors"
           << "\n     -- number of proc environnement : " << Environment::worldComm().globalSize()
           << "\n     -- environement rank : " << Environment::worldComm().rank()
           << "\n     -- global rank : " << this->worldComm().globalRank()
           << "\n     -- local rank : " << this->worldComm().localRank()
        //<< "\n   Matrix Index Start"
        //  << "\n     -- rowstart : " << this->rowStartInMatrix()
        //  << "\n     -- colstart : " << this->colStartInMatrix()
           << "\n   Numerical Solver"
           << "\n     -- solver : " << M_solverName;
    if ( M_useThermodynModel )
        *_ostr << M_thermodynModel->getInfo()->str();
    if ( M_algebraicFactory )
        *_ostr << M_algebraicFactory->getInfo()->str();
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
        *_ostr << this->meshALE()->getInfo()->str();
#endif
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::setModelName( std::string const& type )
{
    // if pde change -> force to rebuild all algebraic data at next solve
    if ( type != M_modelName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "Stokes" )
    {
        M_modelName="Stokes";
        M_solverName="LinearSystem";
    }
    else if ( type == "Oseen" ) // not realy a model but a solver for navier stokes
    {
        M_modelName="Navier-Stokes";
        M_solverName="Oseen";
    }
    else if ( type == "Navier-Stokes" )
    {
        M_modelName="Navier-Stokes";
        M_solverName="Newton";
    }
    else
        CHECK( false ) << "invalid modelName "<< type << "\n";
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::modelName() const
{
    return M_modelName;
}
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::setSolverName( std::string const& type )
{
    // if solver change -> force to rebuild all algebraic data at next solve
    if ( type != M_solverName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "LinearSystem" )
        M_solverName="LinearSystem";
    else if ( type == "Oseen" )
        M_solverName="Oseen";
    else if ( type == "Picard" || type == "FixPoint" )
        M_solverName="Picard";
    else if ( type == "Newton" )
        M_solverName="Newton";
    else
        CHECK( false ) << "invalid solver name " << type << "\n";
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::solverName() const
{
    return M_solverName;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::dynamicViscosityLaw() const
{
    return this->densityViscosityModel()->dynamicViscosityLaw();
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::setDynamicViscosityLaw( std::string const& type )
{
    // if viscosity model change -> force to rebuild all algebraic data at next solve
    if ( type != this->densityViscosityModel()->dynamicViscosityLaw() )
        this->setNeedToRebuildCstPart( true );

    this->densityViscosityModel()->setDynamicViscosityLaw( type );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::setDoExport(bool b)
{
    if ( M_exporter )
        M_exporter->setDoExport( b );
    if ( M_exporter_ho )
        M_exporter_ho->setDoExport( b );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("FluidMechanics","exportResults", (boost::format("start at time %1%")%time).str() );
    this->timerTool("PostProcessing").start();

    if ( this->isMoveDomain() && this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::ALEMesh ) )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        this->meshALE()->exportResults( time );
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

    if ( nOrderGeo == 1 )
    {
        this->exportResultsImpl( time );

        if ( this->hasMarkerPressureBC() && M_spaceLagrangeMultiplierPressureBC &&
             this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::LagrangeMultiplierPressureBC ) )
        {
            std::string geoExportType="static";//change_coords_only, change, static
            if ( !M_exporterLagrangeMultiplierPressureBC )
                M_exporterLagrangeMultiplierPressureBC = exporter( _mesh=M_spaceLagrangeMultiplierPressureBC->mesh(),
                                                                   _name="ExportLagrangeMultiplierPressureBC",
                                                                   _geo=geoExportType,
                                                                   _worldcomm=M_spaceLagrangeMultiplierPressureBC->worldComm(),
                                                                   _path=this->exporterPath() );
            M_exporterLagrangeMultiplierPressureBC->step( time )->add( prefixvm(this->prefix(),"pressurebc-lambda1"),
                                                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressurebc-lambda1")),
                                                                       *M_fieldLagrangeMultiplierPressureBC1 );
            if ( nDim == 3 )
                M_exporterLagrangeMultiplierPressureBC->step( time )->add( prefixvm(this->prefix(),"pressurebc-lambda2"),
                                                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressurebc-lambda2")),
                                                                           *M_fieldLagrangeMultiplierPressureBC2 );
            M_exporterLagrangeMultiplierPressureBC->save();
        }
    }

    if ( M_isHOVisu )
    {
        this->exportResultsImplHO( time );
    }

    if ( M_useThermodynModel && !M_thermodynModel->mesh()->isSameMesh( this->mesh() ) )
    {
        M_thermodynModel->exportResults( time );
    }


    this->exportMeasures( time );


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
    if ( !M_exporter ) return;
    if ( !M_exporter->doExport() ) return;


#if defined( FEELPP_MODELS_HAS_MESHALE )
    // because write geofile at each step ( TODO fix !!! )
    if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
        this->meshALE()->revertReferenceMesh();
#endif
    //M_exporter->step( time )->setMesh( M_mesh );
    bool hasFieldToExport = false;
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pid ) )
    {
        M_exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Velocity ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity")),
                                       M_Solution->template element<0>() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pressure ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"pressure"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure")),
                                       M_Solution->template element<1>() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Vorticity ) )
    {
        this->updateVorticity();
        M_exporter->step( time )->add( prefixvm(this->prefix(),"vorticity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"vorticity")),
                                       this->fieldVorticity() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) )
    {
        this->updateNormalStressOnCurrentMesh();
        M_exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")),
                                       this->fieldNormalStress() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
    {
        this->updateWallShearStress();
        M_exporter->step( time )->add( prefixvm(this->prefix(),"wallshearstress"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")),
                                       this->fieldWallShearStress() );
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Density ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"density"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"density")),
                                       this->densityViscosityModel()->fieldDensity() );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Viscosity ) )
    {
        if ( !M_XhNormalBoundaryStress ) this->createFunctionSpacesNormalStress();
        auto uCur = M_Solution->template element<0>();
        auto pCur = M_Solution->template element<1>();
        auto myViscosity = Feel::vf::FeelModels::fluidMecViscosity<2*nOrderVelocity>(uCur,pCur,*this->densityViscosityModel());
        auto viscosityField = M_XhNormalBoundaryStress->compSpace()->element(myViscosity);
        M_exporter->step( time )->add( prefixvm(this->prefix(),"viscosity"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"viscosity")),
                                       viscosityField );
        hasFieldToExport = true;
    }
    if ( this->isMoveDomain() )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )

        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Displacement ) )
        {
            auto drm = M_meshALE->dofRelationShipMap();
            auto thedisp = M_meshALE->functionSpace()->element();
            for (size_type i=0;i<thedisp.nLocalDof();++i)
                thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i);

            M_exporter->step( time )->add( prefixvm(this->prefix(),"displacement"),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacement")),
                                           thedisp );
            hasFieldToExport = true;
        }

        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::ALEMesh ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"displacementOnInterface"),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacementOnInterface")),
                                           this->meshDisplacementOnInterface() );
            M_exporter->step( time )->add( prefixvm(this->prefix(),"mesh-velocity"),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"mesh-velocity")),
                                           this->meshVelocity() );
            M_exporter->step( time )->add( prefixvm(this->prefix(),"mesh-velocity-interface"),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),"mesh-velocity-interface")),
                                           this->meshVelocity2() );
            hasFieldToExport = true;
        }
#endif
    }
    if ( M_useThermodynModel && M_thermodynModel->mesh()->isSameMesh( this->mesh() ) )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                       M_thermodynModel->fieldTemperature() );
        hasFieldToExport = true;
    }
    for ( std::string const& userFieldName : M_postProcessUserFieldExported )
    {
        if ( this->hasFieldUserScalar( userFieldName ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),userFieldName),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),userFieldName)),
                                           this->fieldUserScalar( userFieldName ) );
            hasFieldToExport = true;
        }
        else if ( this->hasFieldUserVectorial( userFieldName ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),userFieldName),
                                           prefixvm(this->prefix(),prefixvm(this->subPrefix(),userFieldName)),
                                           this->fieldUserVectorial( userFieldName ) );
            hasFieldToExport = true;
        }
    }

    //----------------------//
    if ( hasFieldToExport )
    {
        M_exporter->save();
        this->log("FluidMechanics","exportResults", "save done" );
    }

    //----------------------//
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
        this->meshALE()->revertMovingMesh();
#endif
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResultsImplHO( double time )
{
#if 1//defined(FEELPP_HAS_VTK)
    if ( !M_exporter_ho ) return;
    if ( !M_exporter_ho->doExport() ) return;

    // because write geofile at each step ( TODO fix !!! )
    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //    this->meshALE()->revertReferenceMesh();
    //M_exporter_ho->step( time )->setMesh( M_velocityVisuHO->mesh() );
    bool hasFieldToExport = false;
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pid ) )
    {
        M_exporter_ho->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Velocity ) )
    {
        M_opIvelocity->apply(M_Solution->template element<0>(),*M_velocityVisuHO);
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity_ho")), *M_velocityVisuHO );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Pressure ) )
    {
        M_opIpressure->apply(M_Solution->template element<1>(),*M_pressureVisuHO);
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure_ho")), *M_pressureVisuHO );
        hasFieldToExport = true;
    }

    if ( this->isMoveDomain() )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        auto drm = M_meshALE->dofRelationShipMap();
#if 0
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
#endif
        if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::Displacement ) )
        {
            //M_opImeshdisp->apply( thedisp , *M_meshdispVisuHO);
            M_opImeshdisp->apply( *M_meshALE->displacement() , *M_meshdispVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshdisp_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshdisp_ho")), *M_meshdispVisuHO );
            hasFieldToExport = true;
        }

        if ( false /*M_doExportAll*/ )
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

    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::NormalStress ) )
    {
        this->updateNormalStressOnCurrentMesh();
        M_opIstress->apply( this->fieldNormalStress(),*M_normalStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"normalstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")), *M_normalStressVisuHO );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( FluidMechanicsPostProcessFieldExported::WallShearStress ) )
    {
        this->updateWallShearStress();
        M_opIstress->apply( this->fieldWallShearStress(),*M_fieldWallShearStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"wallshearstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")), *M_fieldWallShearStressVisuHO );
        hasFieldToExport = true;
    }

    if ( hasFieldToExport )
        M_exporter_ho->save();

    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //   this->meshALE()->revertMovingMesh();

#endif
}


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    bool hasMeasure = false;

    // forces (lift,drag) measures
    for ( auto const& ppForces : M_postProcessMeasuresForces )
    {
        CHECK( ppForces.meshMarkers().size() == 1 ) << "TODO";
        auto measuredForce = this->computeForce( ppForces.meshMarkers().front() );
        std::string name = ppForces.name();
        this->postProcessMeasuresIO().setMeasure( "drag_"+name, measuredForce(0,0) );
        this->postProcessMeasuresIO().setMeasure( "lift_"+name, measuredForce(1,0) );
        hasMeasure = true;
    }
    // flow rate measures
    for ( auto const& ppFlowRate : M_postProcessMeasuresFlowRate )
    {
        double valFlowRate = this->computeFlowRate( ppFlowRate.meshMarkers(), ppFlowRate.useExteriorNormal() );
        this->postProcessMeasuresIO().setMeasure("flowrate_"+ppFlowRate.name(),valFlowRate);
        hasMeasure = true;
    }

    auto itFindMeasures = this->modelProperties().postProcess().find("Measures");
    if ( itFindMeasures != this->modelProperties().postProcess().end() )
    {
        bool hasMeasuresPressure = std::find( itFindMeasures->second.begin(), itFindMeasures->second.end(), "Pressure" ) != itFindMeasures->second.end();
        bool hasMeasuresVelocityDivergence = std::find( itFindMeasures->second.begin(), itFindMeasures->second.end(), "VelocityDivergence" ) != itFindMeasures->second.end();
        double area = 0;
        if ( hasMeasuresPressure || hasMeasuresVelocityDivergence )
            area = this->computeMeshArea();
        if ( hasMeasuresPressure )
        {
            double pressureSum = this->computePressureSum();
            double pressureMean = pressureSum/area;
            this->postProcessMeasuresIO().setMeasure("pressure_sum",pressureSum);
            this->postProcessMeasuresIO().setMeasure("pressure_mean",pressureMean);
            hasMeasure = true;
        }
        if ( hasMeasuresVelocityDivergence )
        {
            double velocityDivergenceSum = this->computeVelocityDivergenceSum();
            double velocityDivergenceMean = velocityDivergenceSum/area;
            double velocityDivergenceNormL2 = this->computeVelocityDivergenceNormL2();
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_sum",velocityDivergenceNormL2);
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_mean",velocityDivergenceMean);
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_normL2",velocityDivergenceNormL2);
            hasMeasure = true;
        }
    }


    // point measures
    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->modelProperties().postProcess().setParameterValues( paramValues );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint() )
    {
        auto const& ptPos = evalPoints.pointPosition();
        if ( !ptPos.hasExpression() )
            continue;
        node_type ptCoord(3);
        for ( int c=0;c<3;++c )
            ptCoord[c]=ptPos.value()(c);

        auto const& fields = evalPoints.fields();
        for ( std::string const& field : fields )
        {
            if ( field == "velocity" )
            {
                std::string ptNameExport = (boost::format("velocity_%1%")%ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId("velocity",ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextVelocity->replace( ptIdInCtx, ptCoord );
            }
            else if ( field == "pressure" )
            {
                std::string ptNameExport = (boost::format("pressure_%1%")%ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId("pressure",ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextPressure->replace( ptIdInCtx, ptCoord );
            }
        }
    }
    if ( M_postProcessMeasuresContextVelocity && this->postProcessMeasuresEvaluatorContext().has("velocity") )
    {
        auto evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextVelocity,
                                                _expr=idv(this->fieldVelocity()) );
        for ( int ctxId=0;ctxId<M_postProcessMeasuresContextVelocity->nPoints();++ctxId )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has( "velocity", ctxId ) ) continue;
            std::string const& ptNameExport = this->postProcessMeasuresEvaluatorContext().name( "velocity",ctxId );
            std::vector<double> vecValues = { evalAtNodes( ctxId*nDim ) };
            if ( nDim > 1 ) vecValues.push_back( evalAtNodes( ctxId*nDim+1 ) );
            if ( nDim > 2 ) vecValues.push_back( evalAtNodes( ctxId*nDim+2 ) );
            this->postProcessMeasuresIO().setMeasureComp( ptNameExport, vecValues );
            //std::cout << "export point " << ptNameExport << " with node " << M_postProcessMeasuresContextVelocity->node( ctxId ) << "\n";
            hasMeasure = true;
        }
    }
    if ( M_postProcessMeasuresContextPressure && this->postProcessMeasuresEvaluatorContext().has("pressure") )
    {
        auto evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextPressure,
                                                _expr=idv(this->fieldPressure()) );
        for ( int ctxId=0;ctxId<M_postProcessMeasuresContextPressure->nPoints();++ctxId )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has( "pressure", ctxId ) ) continue;
            std::string ptNameExport = this->postProcessMeasuresEvaluatorContext().name( "pressure",ctxId );
            this->postProcessMeasuresIO().setMeasure( ptNameExport, evalAtNodes( ctxId ) );
            //std::cout << "export point " << ptNameExport << " with node " << M_postProcessMeasuresContextPressure->node( ctxId ) << "\n";
            hasMeasure = true;
        }
    }

    if ( hasMeasure )
    {
        this->postProcessMeasuresIO().setParameter( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
    }
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

        std::string saveStressTensorLawType = this->dynamicViscosityLaw();
        std::string savePdeType = this->modelName();
        // prepare Stokes-stationary config
        this->setDynamicViscosityLaw( "newtonian" );
        this->setModelName( "Stokes" );
        this->setStationary( true );
        // possibility to config a specific time which appear in bc
        double timeUsed = doption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary.time-value-used-in-bc");
        this->updateTime( timeUsed );

        this->solve();
        this->hasSolveStokesStationaryAtKickOff( true );

        if ( boption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary.do-export") )
            this->exportResults(this->timeInitial());
        // revert parameters
        this->setDynamicViscosityLaw( saveStressTensorLawType );
        this->setModelName( savePdeType );
        this->setStationary( false );

        this->initTimeStep();
    }

    if ( this->startBySolveNewtonian() && this->densityViscosityModel()->dynamicViscosityLaw() != "newtonian" &&
         !this->hasSolveNewtonianAtKickOff() && !this->doRestart() )
    {
        this->log("FluidMechanics","solve", "start by solve newtonian" );

        std::string saveStressTensorLawType = this->dynamicViscosityLaw();
        this->setDynamicViscosityLaw("newtonian");
        this->solve();
        this->hasSolveNewtonianAtKickOff( true );
        this->setDynamicViscosityLaw( saveStressTensorLawType );
    }

    //--------------------------------------------------
    // run solver
    std::string algebraicSolver = M_solverName;
    if ( algebraicSolver == "Oseen" )
        algebraicSolver = "LinearSystem";
    M_algebraicFactory->solve( algebraicSolver, this->blockVectorSolution().vector() );
    // update sub vector
    M_blockVectorSolution.localize();

    //--------------------------------------------------
    // update windkessel solution ( todo put fluidOutletWindkesselPressureDistal and Proximal in blockVectorSolution )
    int cptBlock=1;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        ++cptBlock;
    if (this->hasMarkerDirichletBClm())
        ++cptBlock;
    if ( this->hasMarkerPressureBC() )
    {
        ++cptBlock;
        if ( nDim == 3 )
            ++cptBlock;
    }
    if (this->hasFluidOutletWindkesselImplicit() )
    {
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "implicit" )
                continue;

            if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
            {
                M_fluidOutletWindkesselPressureDistal[k] = M_blockVectorSolution(cptBlock)->operator()(0);
                M_fluidOutletWindkesselPressureProximal[k]= M_blockVectorSolution(cptBlock)->operator()(1);
            }
            ++cptBlock;
        }
    }
    //--------------------------------------------------
    // run thermodyn solver if not strong coupling
    if ( M_useThermodynModel && !M_useGravityForce )
    {
        M_thermodynModel->updateFieldVelocityConvection( idv(fieldVelocity()) );
        M_thermodynModel->solve();
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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    if ( this->definePressureCstMethod() == "algebraic" )
    {
        auto upSol = this->functionSpace()->element( sol, this->rowStartInVector() );
        auto pSol = upSol.template element<1>();
        CHECK( M_definePressureCstAlgebraicOperatorMeanPressure ) << "mean pressure operator does not init";
        double meanPressureCurrent = inner_product( *M_definePressureCstAlgebraicOperatorMeanPressure, pSol );
        pSol.add( -meanPressureCurrent );
    }
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

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


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateDefinePressureCst()
{
    if ( this->definePressureCstMethod() == "lagrange-multiplier" && !M_XhMeanPressureLM )
    {
        if ( M_densityViscosityModel->isDefinedOnWholeMesh() )
            M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
        else
            M_XhMeanPressureLM = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                                                 _range=M_rangeMeshElements );
    }
    else if ( this->definePressureCstMethod() == "algebraic" )
    {
        auto p = this->functionSpacePressure()->element();
        M_definePressureCstAlgebraicOperatorMeanPressure = form1_mean(_test=this->functionSpacePressure(),
                                                                      _range=M_rangeMeshElements,
                                                                      _expr=id(p) ).vectorPtr();
        M_definePressureCstAlgebraicOperatorMeanPressure->close();
    }
}


//---------------------------------------------------------------------------------------------------------//


namespace detail
{

template <typename VorticityFieldType, typename VelocityFieldType, typename RangeEltType >
void
updateVorticityImpl( VelocityFieldType /*const*/& fieldVelocity, VorticityFieldType & fieldVorticity, RangeEltType const& rangeElt, GeomapStrategyType geomapUsed, mpl::int_<2> /**/ )
{
    fieldVorticity.on(_range=rangeElt,
                      _expr=/*vf::abs*/(vf::dxv( fieldVelocity.template comp<ComponentType::Y>())
                                    -vf::dyv(fieldVelocity.template comp<ComponentType::X>())),
                      _geomap=geomapUsed );
}
template <typename VorticityFieldType, typename VelocityFieldType, typename RangeEltType >
void
updateVorticityImpl( VelocityFieldType /*const*/& fieldVelocity, VorticityFieldType & fieldVorticity, RangeEltType const& rangeElt, GeomapStrategyType geomapUsed, mpl::int_<3> /**/ )
{
    fieldVorticity.on(_range=rangeElt,
                      _expr=vf::curlv(fieldVelocity),
                      _geomap=geomapUsed );
}

} //  namesapce detail

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateVorticity()
{

    if (!M_XhVorticity) this->createFunctionSpacesVorticity();

    detail::updateVorticityImpl( this->fieldVelocity(),*M_fieldVorticity, M_rangeMeshElements, this->geomap(), mpl::int_<nDim>() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateTimeStepBDF()
{
    this->log("FluidMechanics","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // windkessel outlet
    if (this->hasFluidOutletWindkessel() )
    {
        bool doWriteOnDisk = this->worldComm().isMasterRank();
        // warning, in the implicit case, not all process have info
        // only process which have a face on outlet
        if ( this->hasFluidOutletWindkesselImplicit() )
            doWriteOnDisk = (bool)(M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0);

        if ( doWriteOnDisk )
        {
            std::string nameFile = this->rootRepository() + "/" + prefixvm(this->prefix(),"fluidoutletbc.windkessel.data");
            std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
            file.precision( 8 );
            file.setf( std::ios::scientific );
            file.width( 15 );
            file.setf( std::ios::left );
            file << M_bdf_fluid->iteration();
            file.width( 20 );
            file << this->time();

            for (int k=0;k<this->nFluidOutlet();++k)
            {
                if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;

                file.width( 20 );
                file << M_fluidOutletWindkesselPressureDistal[k];
                file.width( 20 );
                file << M_fluidOutletWindkesselPressureProximal[k];
            }
            file << "\n";
            file.close();
        }

        // update to next timestep
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" ) continue;

            const int sizeOld = M_fluidOutletWindkesselPressureDistal_old.find(k)->second.size();
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
    if ( M_useThermodynModel )
        M_thermodynModel->updateTimeStep();

    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf_fluid->time() );

    // update user functions which depend of time only
    this->updateUserFunctions(true);

    // maybe rebuild cst jacobian or linear
    if ( M_algebraicFactory &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (this->solverName() == "Newton" && !this->rebuildLinearPartInJacobian() )
        {
            this->log("FluidMechanics","updateTimeStepBDF", "do rebuildCstJacobian" );
            M_algebraicFactory->rebuildCstJacobian(M_Solution);
        }
        else if (this->solverName() == "LinearSystem" && !this->rebuildCstPartInLinearSystem())
        {
            this->log("FluidMechanics","updateTimeStepBDF", "do rebuildCstLinearPDE" );
            M_algebraicFactory->rebuildCstLinearPDE(M_Solution);
        }
    }


    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("FluidMechanics","updateTimeStepBDF", "finish" );
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
    this->meshALE()->revertReferenceMesh( false );

    M_fieldNormalStressRefMesh->zero();
    if ( listMarkers.empty() )
        M_fieldNormalStressRefMesh->on(_range=boundaryfaces(this->mesh()),
                                       _expr=val(Sigmav*trans(InvFa)*N()),
                                       _geomap=this->geomap() );
    else
        M_fieldNormalStressRefMesh->on(_range=markedfaces(this->mesh(),listMarkers/*this->markersNameMovingBoundary()*/),
                                       _expr=val(Sigmav*trans(InvFa)*N()),
                                       _geomap=this->geomap() );

    this->meshALE()->revertMovingMesh( false );

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

    this->meshALE()->revertReferenceMesh( false );


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

    this->meshALE()->revertMovingMesh( false );

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

    this->meshALE()->revertReferenceMesh();

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

    this->meshALE()->revertMovingMesh();

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
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updatePicardConvergence( vector_ptrtype const& Unew, vector_ptrtype const& Uold ) const
{
    using namespace Feel::vf;

    auto U1 = this->functionSpace()->element();U1 = *Unew;
    auto u1 = U1.template element<0>();
    auto p1 = U1.template element<1>();

    auto U2 = this->functionSpace()->element();U2 = *Uold;
    auto u2 = U2.template element<0>();
    auto p2 = U2.template element<1>();

    double err_u = integrate(_range=M_rangeMeshElements,
                             //_expr= trans(idv(u1)-idv(u2))*(idv(u1)-idv(u2)),
                             _expr= inner(idv(u1)-idv(u2)),
                             _geomap=this->geomap() ).evaluate()(0,0);
    double err_p = integrate(_range=M_rangeMeshElements,
                             //_expr= (idv(p1)-idv(p2))*(idv(p1)-idv(p2)),
                             _expr= inner(idv(p1)-idv(p2)),
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
    M_meshVelocityInterface->on( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                                 _expr=vf::idv(M_meshALE->velocity()),
                                 _geomap=this->geomap() );
    sync( *M_meshVelocityInterface, "=", M_dofsVelocityInterfaceOnMovingBoundary);

    //-------------------------------------------------------------------//
    // semi implicit optimisation
    if (this->useFSISemiImplicitScheme())
    {
        //this->meshALE()->revertReferenceMesh();
        this->updateNormalStressOnReferenceMeshOptPrecompute(this->markersNameMovingBoundary());
        //this->meshALE()->revertMovingMesh();
    }

    //-------------------------------------------------------------------//
    // move winkessel submesh
    if ( this->hasFluidOutletWindkesselImplicit() )
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
#if 0 //defined(FEELPP_HAS_VTK)
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
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMeshArea( std::string const& marker ) const
{
    return this->computeMeshArea( std::list<std::string>( { marker } ) );
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMeshArea( std::list<std::string> const& markers ) const
{
    double area = 0;
    if ( markers.empty() || markers.front().empty() )
        area = integrate(_range=M_rangeMeshElements,//elements(this->mesh()),
                         _expr=cst(1.),
                         _geomap=this->geomap() ).evaluate()(0,0);
    else
        area = integrate(_range=markedelements(this->mesh(),markers),
                         _expr=cst(1.),
                         _geomap=this->geomap() ).evaluate()(0,0);
    return area;
}


FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
Eigen::Matrix<typename FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::value_type,
              FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nDim,1>
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeForce(std::string const& markerName) const
{
    using namespace Feel::vf;

    auto solFluid = this->fieldVelocityPressurePtr();
    auto u = solFluid->template element<0>();
    auto p = solFluid->template element<1>();
#if 0
    //deformation tensor
    auto defv = sym(gradv(u));
    // Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Tenseur des contraintes
    auto Sigmav = val(-idv(p)*Id + 2*idv(this->densityViscosityModel()->fieldMu())*defv);
#endif
    auto Sigmav/*ViscousStressTensorExpr*/ = Feel::vf::FeelModels::fluidMecNewtonianStressTensor<2*nOrderVelocity>(u,p,*this->densityViscosityModel(),true);

    return integrate(_range=markedfaces(M_mesh,markerName),
                     _expr= Sigmav*N(),
                     _geomap=this->geomap() ).evaluate();

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeFlowRate( std::string const& marker, bool useExteriorNormal ) const
{
    return this->computeFlowRate( std::list<std::string>( { marker } ),useExteriorNormal );
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeFlowRate( std::list<std::string> const& markers, bool useExteriorNormal ) const
{
    using namespace Feel::vf;

    auto const& u = this->fieldVelocity();
    double res = integrate(_range=markedfaces(this->mesh(),markers),
                                 _expr= inner(idv(u),N()),
                                 _geomap=this->geomap() ).evaluate()(0,0);
    if ( !useExteriorNormal )
        res = -res;

    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computePressureSum() const
{
    auto const& p = this->fieldPressure();
    double res = integrate(_range=M_rangeMeshElements,
                           _expr= idv(p),
                           _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computePressureMean() const
{
    double area = this->computeMeshArea();
    double res = (1./area)*this->computePressureSum();
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceSum() const
{
    auto const& u = this->fieldVelocity();
    double res = integrate(_range=M_rangeMeshElements,
                           _expr= divv(u),
                           _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceMean() const
{
    double area = this->computeMeshArea();
    double res = (1./area)*this->computeVelocityDivergenceSum();
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceNormL2() const
{
    using namespace Feel::vf;

    auto const& u = this->fieldVelocity();

    double res = math::sqrt( integrate(_range=M_rangeMeshElements,
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
        this->meshALE()->revertReferenceMesh();
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
        this->meshALE()->revertMovingMesh();
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
        this->meshALE()->revertReferenceMesh();

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

        this->meshALE()->revertMovingMesh();
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

    bool doStabPSPG = M_stabilizationGLS &&
        ( this->stabilizationGLSType()== "pspg" ||
          this->stabilizationGLSType()== "supg-pspg" ||
          this->stabilizationGLSType() == "gls" );
    if ( doStabPSPG ||
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
    if ( this->hasMarkerPressureBC() )
    {
        ++nBlock;
        if ( nDim == 3 )
            ++nBlock;
    }
    if ( this->hasFluidOutletWindkesselImplicit() )
        nBlock += this->nFluidOutletWindkesselImplicit();
    if ( M_useThermodynModel && M_useGravityForce )
        nBlock += M_thermodynModel->nBlockMatrixGraph();
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
    if ( this->hasMarkerPressureBC() )
    {
        BlocksStencilPattern patCouplingLM(1,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCouplingLM(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(indexBlock,0) = stencil(_test=M_spaceLagrangeMultiplierPressureBC,_trial=this->functionSpace(),
                                             _pattern_block=patCouplingLM,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=this->functionSpace(),_trial=M_spaceLagrangeMultiplierPressureBC,
                                             _pattern_block=patCouplingLM.transpose(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
        if ( nDim == 3 )
        {
            myblockGraph(indexBlock,0) = myblockGraph(indexBlock-1,0);
            myblockGraph(0,indexBlock) = myblockGraph(0,indexBlock-1);
            ++indexBlock;
        }
    }
    if ( this->hasFluidOutletWindkesselImplicit() )
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
            if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "implicit" )
                continue;
            std::string markerOutlet = std::get<0>( M_fluidOutletsBCType[k] );
            // first column
            auto rangeFirstCol1 = stencilRange<0,0>( markedelements(this->fluidOutletWindkesselMesh(),markerOutlet) );
            auto rangeFirstCol2 = stencilRange<1,0>( markedelements(this->fluidOutletWindkesselMesh(),markerOutlet) );
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
        indexBlock += this->nFluidOutletWindkesselImplicit();//this->nFluidOutlet();
    }

    if ( M_useThermodynModel && M_useGravityForce )
    {
        myblockGraph(indexBlock,indexBlock) = M_thermodynModel->buildBlockMatrixGraph()(0,0);

        BlocksStencilPattern patCoupling1(1,space_fluid_type::nSpaces,size_type(Pattern::ZERO));
        patCoupling1(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(indexBlock,0) = stencil(_test=M_thermodynModel->spaceTemperature(),
                                             _trial=this->functionSpace(),
                                             _pattern_block=patCoupling1,
                                             _diag_is_nonzero=false,_close=false)->graph();

        BlocksStencilPattern patCoupling2(space_fluid_type::nSpaces,1,size_type(Pattern::ZERO));
        patCoupling2(0,0) = size_type(Pattern::COUPLED);
        myblockGraph(0,indexBlock) = stencil(_test=this->functionSpace(),
                                             _trial=M_thermodynModel->spaceTemperature(),
                                             _pattern_block=patCoupling2,
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
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

    if ( this->hasFluidOutletWindkesselImplicit() )
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
                    (*res)[k][i]+=2*this->nFluidOutletWindkesselImplicit();//nFluidOutlet();
            }
        }

        if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
        {
            const size_type velocitySplitSize = (*res)[0].size();
            (*res)[0].resize( velocitySplitSize + 2*this->nFluidOutletWindkesselImplicit/*nFluidOutlet*/(), 0 );
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
    if ( this->hasFluidOutletWindkesselImplicit() )
        res += 2*this->nFluidOutletWindkesselImplicit();
    return res;
}

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateBlockVectorSolution()
{
    // copy velocity/pressure in block
#if 0
    auto & vecAlgebraic = M_blockVectorSolution.vector();
    auto const& fieldVelPres = this->fieldVelocityPressure();
    for (int k=0;k< this->functionSpace()->nLocalDofWithGhost() ;++k)
        vecAlgebraic->set(k, fieldVelPres(k) );
#else
    //M_blockVectorSolution.vector()->close();
    M_blockVectorSolution.setVector( *M_blockVectorSolution.vector(), this->fieldVelocityPressure(), 0 );
#endif

    // do nothing for others block (fields define only in blockVectorSolution)
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {}
    if (this->hasMarkerDirichletBClm())
    {}
    if ( this->hasFluidOutletWindkesselImplicit() )
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

    if ( !M_couplingFSI_RNG_useInterfaceOperator )
    {
        if (this->doRestart())
            this->meshALE()->revertReferenceMesh();

        auto const& u = this->fieldVelocity();
        form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=M_couplingFSI_RNG_matrix,
               _rowstart=this->rowStartInMatrix(),
               _colstart=this->colStartInMatrix() ) +=
            integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                       _expr=inner(idt(u),id(u)),
                       _geomap=this->geomap() );
        M_couplingFSI_RNG_matrix->close();

        if (this->doRestart())
            this->meshALE()->revertMovingMesh();

        return;
    }

    //-----------------------------------------------------------------//
    // assembly : first version
    //-----------------------------------------------------------------//
#if 0 /////////////////////////
    this->meshALE()->revertReferenceMesh();
    form2( _test=Xh,_trial=Xh,_matrix=M_couplingFSI_RNG_matrix,
           _pattern=size_type(Pattern::COUPLED) ) +=
        integrate( _range=markedfaces(this->mesh(),this->markersNameMovingBoundary()),
                   _expr=this->couplingFSI_RNG_coeffForm2()*inner(idt(u),id(u)),
                   _geomap=this->geomap() );
    M_couplingFSI_RNG_matrix->close();
    this->meshALE()->revertMovingMesh();


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
        //        auto __face_it = faceMarked.template get<1>();
        // auto __face_en = faceMarked.template get<2>();
        // for( ; __face_it != __face_en; ++__face_it )
        // {
        //    auto const& face = boost::unwrap_ref( *__face_it );
        auto  const& face = boost::unwrap_ref( faceMarked );
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
            //}
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
    this->meshALE()->revertReferenceMesh();
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
        // auto __face_it = faceMarked.template get<1>();
        // auto __face_en = faceMarked.template get<2>();
        // for( ; __face_it != __face_en; ++__face_it )
        // {
            // auto const& face = *__face_it;
            auto const& face = boost::unwrap_ref( faceMarked );
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
            // }
        }
    }
#else
    auto myoperator = this->couplingFSI_RNG_interfaceOperator();
    std::vector<bool> dofdone(this->functionSpaceVelocity()->nLocalDofWithGhost(),false);
    for ( auto const& faceMarked : markedfaces(this->mesh(),this->markersNameMovingBoundary() ) )
    {
        // auto __face_it = faceMarked.template get<1>();
        // auto __face_en = faceMarked.template get<2>();
        // for( ; __face_it != __face_en; ++__face_it )
        // {
            auto const& face = boost::unwrap_ref( faceMarked );
            // auto const& face = boost::unwrap_ref( *__face_it );
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
        // }
    }
#endif

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateBoundaryConditionsForUse()
{
    auto XhVelocity = this->functionSpaceVelocity();
    auto XhCompVelocity = XhVelocity->compSpace();
    static const uint16_type nDofComponentsVelocity = XhVelocity->dof()->nDofComponents();
    auto mesh = this->mesh();

    std::set<std::string> velocityMarkers;
    std::map<ComponentType,std::set<std::string> > compVelocityMarkers;

    auto & dofsWithValueImposedVelocity = M_dofsWithValueImposed["velocity"];
    dofsWithValueImposedVelocity.clear();

    //-------------------------------------//
    // strong Dirichlet bc on velocity from fsi coupling
#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain() && this->couplingFSIcondition()=="dirichlet-neumann")
        velocityMarkers.insert( this->markersNameMovingBoundary().begin(), this->markersNameMovingBoundary().end() );
#endif
    // strong Dirichlet bc on velocity from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",marker(d) );
        velocityMarkers.insert( listMark.begin(), listMark.end() );
    }
    // strong Dirichlet bc on velocity component from expression
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto listMark = this->markerDirichletBCByNameId( "elimination",marker(d), comp );
            compVelocityMarkers[comp].insert( listMark.begin(), listMark.end() );
        }
    }
    // strong Dirichlet bc on velocity from inlet bc
    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        velocityMarkers.insert( marker );
    }

    //-------------------------------------//
    // distribute mesh markers by entity
    std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > meshMarkersVelocityByEntities;
    std::map<ComponentType, std::tuple< std::list<std::string>,std::list<std::string>,std::list<std::string>,std::list<std::string> > > meshMarkersCompVelocityByEntities;
    meshMarkersVelocityByEntities = detail::distributeMarkerListOnSubEntity( mesh, velocityMarkers );
    for ( auto const& compMarkerPair : compVelocityMarkers )
    {
        meshMarkersCompVelocityByEntities[compMarkerPair.first] = detail::distributeMarkerListOnSubEntity( mesh, compMarkerPair.second );
    }
    //-------------------------------------//
    // on topological faces
    auto const& listMarkedFacesVelocity = std::get<0>( meshMarkersVelocityByEntities );
    for ( auto const& faceWrap : markedfaces(mesh,listMarkedFacesVelocity ) )
    {
        auto const& face = unwrap_ref( faceWrap );
        auto facedof = XhVelocity->dof()->faceLocalDof( face.id() );
        for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
        {
            dofsWithValueImposedVelocity.insert( it->index() );
        }
    }
    // on marked edges (only 3d)
    auto const& listMarkedEdgesVelocity = std::get<1>( meshMarkersVelocityByEntities );
    for ( auto const& edgeWrap : markededges(mesh,listMarkedEdgesVelocity ) )
    {
        auto const& edge = unwrap_ref( edgeWrap );
        auto itEltInfo = edge.elements().begin();
        if ( itEltInfo == edge.elements().end() )
            continue;
        size_type eid = itEltInfo->first;
        uint16_type edgeid_in_element = itEltInfo->second;
        for( auto const& ldof : XhVelocity->dof()->edgeLocalDof( eid, edgeid_in_element ) )
        {
            dofsWithValueImposedVelocity.insert( ldof.index() );
        }
    }
    // on marked points
    auto const& listMarkedPointsVelocity = std::get<2>( meshMarkersVelocityByEntities );
    for ( auto const& pointWrap : markedpoints(mesh,listMarkedPointsVelocity ) )
    {
        auto const& point = unwrap_ref( pointWrap );
        auto itPointInfo = point.elements().begin();
        if ( itPointInfo == point.elements().end() )
            continue;
        size_type eid = itPointInfo->first;
        uint16_type ptid_in_element = itPointInfo->second;
        for( uint16_type c = 0; c < nDofComponentsVelocity; ++c )
        {
            size_type index = XhVelocity->dof()->localToGlobal( eid, ptid_in_element, c ).index();
            dofsWithValueImposedVelocity.insert( index );
        }
    }
    //-------------------------------------//
    // on velocity components
    for ( auto const& meshMarkersPair : meshMarkersCompVelocityByEntities )
    {
        ComponentType comp = meshMarkersPair.first;
        int compDofShift = ((int)comp);
        // topological faces
        auto const& listMarkedFacesCompVelocity = std::get<0>( meshMarkersPair.second );
        for ( auto const& faceWrap : markedfaces(mesh,listMarkedFacesCompVelocity ) )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = XhCompVelocity->dof()->faceLocalDof( face.id() );
            for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
            {
                size_type compdof = it->index();
                size_type thedof = compDofShift +  nDofComponentsVelocity*compdof;
                dofsWithValueImposedVelocity.insert( thedof );
            }
        }
        // edges (only 3d)
        auto const& listMarkedEdgesCompVelocity = std::get<1>( meshMarkersPair.second );
        for ( auto const& edgeWrap : markededges(mesh,listMarkedEdgesCompVelocity ) )
        {
            auto const& edge = unwrap_ref( edgeWrap );
            auto itEltInfo = edge.elements().begin();
            if ( itEltInfo == edge.elements().end() )
                continue;
            size_type eid = itEltInfo->first;
            uint16_type edgeid_in_element = itEltInfo->second;
            for( auto const& ldof : XhCompVelocity->dof()->edgeLocalDof( eid, edgeid_in_element ) )
            {
                size_type compdof = ldof.index();
                size_type thedof = compDofShift + nDofComponentsVelocity*compdof;
                dofsWithValueImposedVelocity.insert( ldof.index() );
            }
        }
        // points
        auto const& listMarkedPointsCompVelocity = std::get<2>( meshMarkersPair.second );
        for ( auto const& pointWrap : markedpoints(mesh,listMarkedPointsCompVelocity ) )
        {
            auto const& point = unwrap_ref( pointWrap );
            auto itPointInfo = point.elements().begin();
            if ( itPointInfo == point.elements().end() )
                continue;
            size_type eid = itPointInfo->first;
            uint16_type ptid_in_element = itPointInfo->second;
            size_type index = XhVelocity->dof()->localToGlobal( eid, ptid_in_element, compDofShift ).index();
            dofsWithValueImposedVelocity.insert( index );
        }
    }

    if ( this->hasMarkerPressureBC() && M_spaceLagrangeMultiplierPressureBC )
    {
        auto & dofsWithValueImposedPressureBC = M_dofsWithValueImposed["pressurebc-lm"];
        dofsWithValueImposedPressureBC.clear();
        for ( auto const& faceWrap : boundaryfaces(M_meshLagrangeMultiplierPressureBC) )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = M_spaceLagrangeMultiplierPressureBC->dof()->faceLocalDof( face.id() );
            for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                dofsWithValueImposedPressureBC.insert( it->index() );
        }
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if ( this->isMoveDomain() )
    {
        for ( auto const& faceWrap : markedfaces(mesh,this->markersNameMovingBoundary() ) )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = M_XhMeshVelocityInterface->dof()->faceLocalDof( face.id() );
            for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
            {
                M_dofsVelocityInterfaceOnMovingBoundary.insert( it->index() );
            }
        }
    }
#endif
}


} // namespace FeelModels
} // namespace Feel


