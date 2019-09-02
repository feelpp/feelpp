/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fluid/fluidmechanics.hpp>

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

#include <feel/feelpde/operatorpcd.hpp>

#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::useExtendedDofTable() const
{
    if ( this->worldComm().localSize() == 1 ) return false;
    bool useExtendedDofTable=false;
    for ( bool hasExt : M_Xh->extendedDofTableComposite() )
        useExtendedDofTable = useExtendedDofTable || hasExt;
    return useExtendedDofTable;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::getInfo() const
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
    for ( std::string const& fieldName : M_postProcessFieldExported )
        doExport_str=(doExport_str.empty())? fieldName : doExport_str + " - " + fieldName;
    for ( std::string const& fieldName : M_postProcessUserFieldExported )
        doExport_str=(doExport_str.empty())? fieldName : doExport_str + " - " + fieldName;

    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
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
        //<< "\n     -- stress tensor law  : " << this->materialProperties()->dynamicViscosityLaw()
           << "\n     -- time mode : " << StateTemporal
           << "\n     -- ale mode  : " << ALEmode
           << "\n     -- gravity  : " << std::boolalpha << M_useGravityForce;
    *_ostr << this->materialProperties()->getInfo()->str();
    // *_ostr << "\n   Physical Parameters"
    //        << "\n     -- rho : " << this->materialProperties()->cstRho()
    //        << "\n     -- mu  : " << this->materialProperties()->cstMu()
    //        << "\n     -- nu  : " << this->materialProperties()->cstNu();
    // *_ostr << this->materialProperties()->getInfo()->str();
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
    if ( this->hasGeoFile() )
        *_ostr << "\n     -- geo file name   : " << this->geoFile();
    *_ostr << "\n     -- mesh file name   : " << this->meshFile()
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
        if ( !M_definePressureCstMarkers.empty() )
        {
            *_ostr << "\n     -- define cst pressure on markers  : ";
            for ( auto const& markers : M_definePressureCstMarkers )
            {
                if ( markers.empty() ) continue;
                *_ostr << "[ ";
                for( auto it=markers.begin(),en=(--markers.end());it!=en;++it )
                    *_ostr << *it << " : ";
                *_ostr << *markers.rbegin() << " ]";
            }
        }
        *_ostr << "\n     -- define cst pressure with method  : " << this->definePressureCstMethod();
        if ( this->definePressureCstMethod() == "penalisation" )
            *_ostr << " ( beta=" << this->definePressureCstPenalisationBeta() << ")";
    }

    if ( !this->isStationary() )
    {
        *_ostr << "\n   Time Discretization"
               << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
               << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
               << "\n     -- time step    : " << this->timeStepBase()->timeStep()
               << "\n     -- type : " << M_timeStepping;
        if ( M_timeStepping == "BDF" )
            *_ostr << " ( order=" << this->timeStepBDF()->timeOrder() << " )";
        else if ( M_timeStepping == "Theta" )
            *_ostr << " ( theta=" << M_timeStepThetaValue << " )";
        *_ostr << "\n     -- restart mode : " << ResartMode
               << "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
        if ( this->timeStepBase()->saveFreq() )
            *_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
                   << "\n     -- file format save : " << this->timeStepBase()->fileFormat();
    }
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setModelName( std::string const& type )
{
    // if pde change -> force to rebuild all algebraic data at next solve
    if ( type != M_modelName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "Stokes" )
        M_modelName="Stokes";
    else if ( type == "StokesTransient" )
        M_modelName="StokesTransient";
    else if ( type == "Oseen" ) // not realy a model but a solver for navier stokes
    {
        M_modelName="Navier-Stokes";
        M_solverName="Oseen";
    }
    else if ( type == "Navier-Stokes" )
        M_modelName="Navier-Stokes";
    else
        CHECK( false ) << "invalid modelName "<< type << "\n";
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::modelName() const
{
    return M_modelName;
}
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setSolverName( std::string const& type )
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solverName() const
{
    return M_solverName;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::isStationaryModel() const
{
    if( this->modelName() == "Stokes" )
        return true;
    else
        return this->isStationary();
}

//---------------------------------------------------------------------------------------------------------//

#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::string const&
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::dynamicViscosityLaw() const
{
    return this->materialProperties()->dynamicViscosityLaw();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setDynamicViscosityLaw( std::string const& type )
{
    // if viscosity model change -> force to rebuild all algebraic data at next solve
    if ( type != this->materialProperties()->dynamicViscosityLaw() )
        this->setNeedToRebuildCstPart( true );

    this->materialProperties()->setDynamicViscosityLaw( type );
}
#endif
//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setDoExport(bool b)
{
    if ( M_exporter )
        M_exporter->setDoExport( b );
    if ( M_exporter_ho )
        M_exporter_ho->setDoExport( b );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("FluidMechanics","exportResults", (boost::format("start at time %1%")%time).str() );
    this->timerTool("PostProcessing").start();

    if ( this->isMoveDomain() && this->hasPostProcessFieldExported( "alemesh" ) )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        this->meshALE()->exportResults( time );
#endif
    }

    if ( nOrderGeo == 1 )
    {
        this->exportFields( time );

        if ( this->hasMarkerPressureBC() && M_spaceLagrangeMultiplierPressureBC &&
             this->hasPostProcessFieldExported( "pressurebc" ) )
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


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::exportFields( double time )
{
    bool hasFieldToExport = this->updateExportedFields( M_exporter, M_postProcessFieldExported, time );
    if ( hasFieldToExport )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        // TODO : only at the first export
        if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
            this->meshALE()->revertReferenceMesh( false );
#endif

        M_exporter->save();

#if defined( FEELPP_MODELS_HAS_MESHALE )
        if ( this->isMoveDomain() && M_exporter->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
            this->meshALE()->revertMovingMesh( false );
#endif

        this->upload( M_exporter->path() );
    }

    bool hasFieldOnTraceToExport = this->updateExportedFieldsOnTrace( M_exporterTrace, M_postProcessFieldOnTraceExported, time );
    if ( hasFieldOnTraceToExport )
    {
        M_exporterTrace->save();
    }

    if  ( hasFieldToExport || hasFieldOnTraceToExport )
    {
        if ( M_exporter )
            this->upload( M_exporter->path() );
        else
            this->upload( M_exporterTrace->path() );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    //exporter->step( time )->setMesh( M_mesh );
    bool hasFieldToExport = false;
    if ( fields.find( "pid" ) != fields.end() )
    {
        exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( fields.find( "velocity" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"velocity"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity")),
                                     M_Solution->template element<0>() );
        hasFieldToExport = true;
    }
    if ( fields.find( "pressure" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"pressure"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure")),
                                     M_Solution->template element<1>() );
        hasFieldToExport = true;
    }
    if ( fields.find( "vorticity" ) != fields.end() )
    {
        this->updateVorticity();
        exporter->step( time )->add( prefixvm(this->prefix(),"vorticity"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"vorticity")),
                                     this->fieldVorticity() );
        hasFieldToExport = true;
    }
    if ( fields.find( "density" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"density"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"density")),
                                     this->materialProperties()->fieldDensity() );
        hasFieldToExport = true;
    }
    if ( fields.find( "viscosity" ) != fields.end() )
    {
        //if ( !M_XhNormalBoundaryStress ) this->createFunctionSpacesNormalStress();
        auto uCur = M_Solution->template element<0>();
        auto pCur = M_Solution->template element<1>();
        auto viscosityField = this->materialProperties()->dynamicViscositySpace()->element();
        for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
        {
            std::string const& matName = rangeData.first;
            auto const& range = rangeData.second;
            //auto const& dynamicViscosity = this->materialProperties()->dynamicViscosity(matName);
            auto myViscosity = Feel::FeelModels::fluidMecViscosity(gradv(uCur),*this->materialProperties(),matName);
            viscosityField.on( _range=range,_expr=myViscosity );
        }
        exporter->step( time )->add( prefixvm(this->prefix(),"viscosity"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"viscosity")),
                                     viscosityField );
        hasFieldToExport = true;
    }
    if ( this->isMoveDomain() )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )

        if ( fields.find( "displacement" ) != fields.end() )
        {
            auto drm = M_meshALE->dofRelationShipMap();
            auto thedisp = M_meshALE->functionSpace()->element();
            for (size_type i=0;i<thedisp.nLocalDof();++i)
                thedisp(drm->dofRelMap()[i])=(*(M_meshALE->displacementInRef()))(i);

            exporter->step( time )->add( prefixvm(this->prefix(),"displacement"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacement")),
                                         thedisp );
            hasFieldToExport = true;
        }

        if ( fields.find( "alemesh" ) != fields.end() )
        {
            exporter->step( time )->add( prefixvm(this->prefix(),"displacementOnInterface"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"displacementOnInterface")),
                                         this->meshDisplacementOnInterface() );
            exporter->step( time )->add( prefixvm(this->prefix(),"mesh-velocity"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"mesh-velocity")),
                                         this->meshVelocity() );
#if 0
            exporter->step( time )->add( prefixvm(this->prefix(),"mesh-velocity-interface"),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),"mesh-velocity-interface")),
                                         this->meshVelocity2() );
#endif
            hasFieldToExport = true;
        }
#endif
    }

    for ( auto const& fieldUserScalar : this->fieldsUserScalar() )
    {
        std::string const& userFieldName = fieldUserScalar.first;
        if ( fields.find( userFieldName ) != fields.end() )
        {
            exporter->step( time )->add( prefixvm(this->prefix(),userFieldName),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),userFieldName)),
                                         this->fieldUserScalar( userFieldName ) );
            hasFieldToExport = true;
        }
    }
    for ( auto const& fieldUserVectorial : this->fieldsUserVectorial() )
    {
        std::string const& userFieldName = fieldUserVectorial.first;
        if ( fields.find( userFieldName ) != fields.end() )
        {
            exporter->step( time )->add( prefixvm(this->prefix(),userFieldName),
                                         prefixvm(this->prefix(),prefixvm(this->subPrefix(),userFieldName)),
                                         this->fieldUserVectorial( userFieldName ) );
            hasFieldToExport = true;
        }
    }

    //----------------------//
    return hasFieldToExport;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::exportResultsImplHO( double time )
{
#if 1//defined(FEELPP_HAS_VTK)
    if ( !M_exporter_ho ) return;
    if ( !M_exporter_ho->doExport() ) return;

    // because write geofile at each step ( TODO fix !!! )
    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //    this->meshALE()->revertReferenceMesh();
    //M_exporter_ho->step( time )->setMesh( M_velocityVisuHO->mesh() );
    bool hasFieldToExport = false;
    if ( this->hasPostProcessFieldExported( "pid" ) )
    {
        M_exporter_ho->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( "velocity" ) )
    {
        M_opIvelocity->apply(M_Solution->template element<0>(),*M_velocityVisuHO);
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity_ho")), *M_velocityVisuHO );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( "pressure" ) )
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
        if ( this->hasPostProcessFieldExported( "displacement" ) )
        {
            //M_opImeshdisp->apply( thedisp , *M_meshdispVisuHO);
            M_opImeshdisp->apply( *M_meshALE->displacement() , *M_meshdispVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"meshdisp_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"meshdisp_ho")), *M_meshdispVisuHO );
            hasFieldToExport = true;
        }

        if ( false /*M_doExportAll*/ )
        {
#if 0
            auto themeshvelocityOnInterfaceVisuHO = M_XhVectorialVisuHO->element();
            auto hola2 = vf::project(_space=this->functionSpaceVelocity(),_range=elements(M_Xh->mesh()),_expr=vf::idv(this->meshVelocity2Ptr()));

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
#endif
        }
#endif // HAS_MESHALE
    }

#if 0
    if ( this->hasPostProcessFieldExported( "normal-stress" ) )
    {
        this->updateNormalStressOnCurrentMesh();
        M_opIstress->apply( this->fieldNormalStress(),*M_normalStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"normalstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")), *M_normalStressVisuHO );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessFieldExported( "wall-shear-stress" ) )
    {
        this->updateWallShearStress();
        M_opIstress->apply( this->fieldWallShearStress(),*M_fieldWallShearStressVisuHO );
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"wallshearstress_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")), *M_fieldWallShearStressVisuHO );
        hasFieldToExport = true;
    }
#endif
    if ( hasFieldToExport )
        M_exporter_ho->save();

    //if (M_isMoveDomain && M_exporter_ho->exporterGeometry()==ExporterGeometry::EXPORTER_GEOMETRY_STATIC)
    //   this->meshALE()->revertMovingMesh();

#endif
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateExportedFieldsOnTrace( export_trace_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;
    bool hasFieldToExport = false;

    if ( fields.find( "normal-stress" ) != fields.end() )
    {
        this->updateNormalStressOnCurrentMesh( "trace_mesh", this->fieldNormalStressPtr() );
        exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"normalstress")),
                                     this->fieldNormalStress() );
        hasFieldToExport = true;
    }
    if ( fields.find( "wall-shear-stress" ) != fields.end() )
    {
        this->updateWallShearStress( "trace_mesh", this->fieldWallShearStressPtr() );
        exporter->step( time )->add( prefixvm(this->prefix(),"wallshearstress"),
                                     prefixvm(this->prefix(),prefixvm(this->subPrefix(),"wallshearstress")),
                                     this->fieldWallShearStress() );
        hasFieldToExport = true;
    }

    return hasFieldToExport;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
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

    
    if ( true )
    {
        bool hasMeasuresPressure = M_postProcessMeasuresFields.find("pressure") != M_postProcessMeasuresFields.end();
        bool hasMeasuresVelocityDivergence = M_postProcessMeasuresFields.find( "velocity-divergence" ) != M_postProcessMeasuresFields.end();
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
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( this->keyword() ) )
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

    auto fieldTuple = hana::make_tuple( std::make_pair( "velocity",this->fieldVelocity() ), std::make_pair( "pressure",this->fieldPressure() ) );
    for ( auto const& ppNorm : this->modelProperties().postProcess().measuresNorm( this->keyword() ) )
    {
        std::map<std::string,double> resPpNorms;
        measureNormEvaluation( this->mesh(), M_rangeMeshElements, ppNorm, resPpNorms, this->symbolsExpr(), fieldTuple );
        for ( auto const& resPpNorm : resPpNorms )
        {
            this->postProcessMeasuresIO().setMeasure( resPpNorm.first, resPpNorm.second );
            hasMeasure = true;
        }
    }
    for ( auto const& ppStat : this->modelProperties().postProcess().measuresStatistics( this->keyword() ) )
    {
        std::map<std::string,double> resPpStats;
        measureStatisticsEvaluation( this->mesh(), M_rangeMeshElements, ppStat, resPpStats, this->symbolsExpr(), fieldTuple );
        for ( auto const& resPpStat : resPpStats )
        {
            this->postProcessMeasuresIO().setMeasure( resPpStat.first, resPpStat.second );
            hasMeasure = true;
        }
    }

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    this->modelProperties().parameters().updateParameterValues();

    auto paramValues = this->modelProperties().parameters().toParameterValues();
    //this->modelProperties().materials().setParameterValues( paramValues );
    this->M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : this->M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    this->M_bcNeumannScalar.setParameterValues( paramValues );
    this->M_bcNeumannVectorial.setParameterValues( paramValues );
    this->M_bcNeumannTensor2.setParameterValues( paramValues );
    this->M_bcPressure.setParameterValues( paramValues );
    M_bcMovingBoundaryImposed.setParameterValues( paramValues );
    this->M_volumicForcesProperties.setParameterValues( paramValues );
    this->updateFluidInletVelocity();
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("FluidMechanics","solve", "start" );
    this->timerTool("Solve").start();

    this->updateParameterValues();

    // copy velocity/pressure in algebraic vector solution (maybe velocity/pressure has been changed externaly)
    this->updateBlockVectorSolution();

    if ( M_applyMovingMeshBeforeSolve && !M_bcMovingBoundaryImposed.empty() )
        this->updateALEmesh();

    if ( this->startBySolveStokesStationary() &&
         !this->hasSolveStokesStationaryAtKickOff() && !this->doRestart() )
    {
        this->log("FluidMechanics","solve", "start by solve stokes stationary" );
#if 0 // TODO
        std::string saveStressTensorLawType = this->dynamicViscosityLaw();
#endif
        std::string savePdeType = this->modelName();
        std::string saveSolverName = this->solverName();
        bool saveIsStationary = this->isStationary();
        // prepare Stokes-stationary config
#if 0 // TODO
        this->setDynamicViscosityLaw( "newtonian" );
#endif
        this->setModelName( "Stokes" );
        this->setStationary( true );
        //this->solve();
        M_algebraicFactory->solve( "LinearSystem", this->blockVectorSolution().vectorMonolithic() );
        this->hasSolveStokesStationaryAtKickOff( true );

        //if ( boption(_prefix=this->prefix(),_name="start-by-solve-stokes-stationary.do-export") )
        //   this->exportResults(this->timeInitial());
        // revert parameters
#if 0 // TODO
        this->setDynamicViscosityLaw( saveStressTensorLawType );
#endif
        this->setModelName( savePdeType );
        this->setSolverName( saveSolverName );
        this->setStationary( saveIsStationary );
    }
#if 0 // TODO
    if ( this->startBySolveNewtonian() && this->materialProperties()->dynamicViscosityLaw() != "newtonian" &&
         !this->hasSolveNewtonianAtKickOff() && !this->doRestart() )
    {
        this->log("FluidMechanics","solve", "start by solve newtonian" );

        std::string saveStressTensorLawType = this->dynamicViscosityLaw();
        this->setDynamicViscosityLaw("newtonian");
        this->solve();
        this->hasSolveNewtonianAtKickOff( true );
        this->setDynamicViscosityLaw( saveStressTensorLawType );
    }
#endif
    //--------------------------------------------------
    // run solver
    std::string algebraicSolver = M_solverName;
    if ( algebraicSolver == "Oseen" )
        algebraicSolver = "LinearSystem";
    M_algebraicFactory->solve( algebraicSolver, this->blockVectorSolution().vectorMonolithic() );
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    if ( this->definePressureCst() && this->definePressureCstMethod() == "algebraic" )
    {
        auto upSol = this->functionSpace()->element( sol, this->rowStartInVector() );
        auto pSol = upSol.template element<1>();
        CHECK( !M_definePressureCstAlgebraicOperatorMeanPressure.empty() ) << "mean pressure operator does not init";
        for ( int k=0;k<M_definePressureCstAlgebraicOperatorMeanPressure.size();++k )
        {
            double meanPressureImposed = 0;
            double meanPressureCurrent = inner_product( *M_definePressureCstAlgebraicOperatorMeanPressure[k].first, pSol );
            for ( size_type dofId : M_definePressureCstAlgebraicOperatorMeanPressure[k].second )
                pSol(dofId) += (meanPressureImposed - meanPressureCurrent);
        }
        sync( pSol, "=" );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->postSolveNewton( rhs,sol );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->postSolveNewton( rhs,sol );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype const& mat = data.matrix();
    vector_ptrtype const& vecSol = data.currentSolution();
    if ( M_preconditionerAttachPMM )
        this->updateInHousePreconditionerPMM( mat, vecSol );
    if ( M_preconditionerAttachPCD )
        this->updateInHousePreconditionerPCD( mat,vecSol, data );
}
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{
    sparse_matrix_ptrtype const& mat = data.jacobian();
    vector_ptrtype const& vecSol = data.currentSolution();
    if ( M_preconditionerAttachPMM )
        this->updateInHousePreconditionerPMM( mat, vecSol );
    if ( M_preconditionerAttachPCD )
        this->updateInHousePreconditionerPCD( mat,vecSol,data );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditionerPMM( sparse_matrix_ptrtype const& /*mat*/,vector_ptrtype const& vecSol ) const
{
    if ( !this->algebraicFactory() )
        return; // TODO : should be use with multiphysics toolboxes as heat-fluid

    bool hasAlreadyBuiltPMM = this->algebraicFactory()->preconditionerTool()->hasAuxiliarySparseMatrix( "pmm" );
    if ( hasAlreadyBuiltPMM && !M_pmmNeedUpdate )
        return;
    sparse_matrix_ptrtype pmmMat;
    if ( hasAlreadyBuiltPMM )
        pmmMat = this->algebraicFactory()->preconditionerTool()->auxiliarySparseMatrix( "pmm" );
    else
    {
        pmmMat = M_backend->newMatrix(_trial=this->functionSpacePressure(), _test=this->functionSpacePressure());
        this->algebraicFactory()->preconditionerTool()->attachAuxiliarySparseMatrix( "pmm", pmmMat );
    }
    CHECK( pmmMat ) << "pmmMat is not initialized";

    auto massbf = form2( _trial=this->functionSpacePressure(), _test=this->functionSpacePressure(),_matrix=pmmMat);
    auto const& p = this->fieldPressure();
    auto coeff = cst(1.)/idv(this->materialProperties()->fieldMu());
    massbf = integrate( _range=M_rangeMeshElements, _expr=coeff*inner( idt(p),id(p) ) );
    pmmMat->close();
    M_pmmNeedUpdate = false;
}
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditionerPCD( sparse_matrix_ptrtype const& mat,vector_ptrtype const& vecSol,  DataUpdateBase & data ) const
{
    this->log("FluidMechanics","updateInHousePreconditionerPCD", "start" );

    CHECK( this->hasOperatorPCD() ) << "operator PCD does not init";

    typedef Feel::Alternatives::OperatorPCD<space_fluid_velocity_type,space_fluid_pressure_type> op_pcd_type;
    std::shared_ptr<op_pcd_type> myOpPCD =
        std::dynamic_pointer_cast<op_pcd_type>( this->operatorPCD() );

    auto U = this->functionSpace()->element( vecSol, this->rowStartInVector() );
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    myOpPCD->updateStart();

    for ( auto const& rangeData : this->materialProperties()->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& therange = rangeData.second;
        auto const& dynamicViscosity = this->materialProperties()->dynamicViscosity(matName);
        auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u),*this->materialProperties(),matName);
        //auto rhoExpr = this->materialProperties()->density( matName ).expr();
        auto const& fieldRho = this->materialProperties()->fieldRho();
        auto rhoExpr = idv( fieldRho );
        if ( this->modelName() == "Stokes" || this->modelName() == "StokesTransient" )
        {
            if (this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, -rhoExpr*idv(this->meshVelocity()), true );
#endif
            }
            else
            {
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, zero<nDim,1>(), false );
            }
        }
        else if ( ( this->modelName() == "Navier-Stokes" && this->solverName() == "Oseen" ) || this->modelName() == "Oseen" )
        {
            auto BetaU = *M_fieldConvectionVelocityExtrapolated;//this->timeStepBDF()->poly();
            auto betaU = BetaU.template element<0>();
            if (this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*( idv(betaU)-idv(this->meshVelocity()) ), true );
#endif
            }
            else
            {
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*idv(betaU), true );
            }
        }
        else if ( this->modelName() == "Navier-Stokes" )
        {
            if (this->isMoveDomain() )
            {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*( idv(u)-idv(this->meshVelocity()) ), true );
#endif
            }
            else
            {
                myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*idv(u), true );
            }
        }

        if ( !this->isStationaryModel() )
        {
            myOpPCD->updateFpMass( therange, rhoExpr*this->timeStepBDF()->polyDerivCoefficient(0) );
        }
        if ( data.hasInfo( "use-pseudo-transient-continuation" ) )
        {
            //Feel::cout << "updsate PCD : use-pseudo-transient-continuation\n";
            //Warning : it's a copy past, should be improve : TODO!
            double pseudoTimeStepDelta = data.doubleInfo("pseudo-transient-continuation.delta");
            auto norm2_uu = this->materialProperties()->fieldRho().functionSpace()->element(); // TODO : improve this (maybe create an expression instead)
            //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(idv(u))/h());
            auto fieldNormu = u.functionSpace()->compSpace()->element( norm2(idv(u)) );
            auto maxu = fieldNormu.max( this->materialProperties()->fieldRho().functionSpace() );
            //auto maxux = u[ComponentType::X].max( this->materialProperties()->fieldRho().functionSpace() );
            //auto maxuy = u[ComponentType::Y].max( this->materialProperties()->fieldRho().functionSpace() );
            //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(vec(idv(maxux),idv(maxux)))/h());
            norm2_uu.on(_range=M_rangeMeshElements,_expr=idv(maxu)/h());

            myOpPCD->updateFpMass( therange, (1./pseudoTimeStepDelta)*idv(norm2_uu) );
        }
    }

    if ( !dynamic_cast<DataUpdateJacobian*>(&data) || !boption(_name="pcd.apply-homogeneous-dirichlet-in-newton",_prefix=this->prefix()) )
    {
        auto const& fieldRho = this->materialProperties()->fieldRho();
        auto rhoExpr = idv( fieldRho );
        for( auto const& d : M_bcDirichlet )
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, name(d), expression(d,this->symbolsExpr()) );
        for( auto const& d : M_bcMovingBoundaryImposed )
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, name(d), idv(M_meshALE->velocity()) );
        for ( auto const& inletbc : M_fluidInletDesc )
        {
            std::string const& marker = std::get<0>( inletbc );
            auto const& inletVel = std::get<0>( M_fluidInletVelocityInterpolated.find(marker)->second );
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, marker, -idv(inletVel)*N() );
        }
    }
    else
        Feel::cout << "PCD NOT UP BC \n";

    // updated from the outside
    for ( auto const& f : M_addUpdateInHousePreconditionerPCD )
        f.second.second( *myOpPCD, data );

    myOpPCD->updateFinish();

    this->log("FluidMechanics","updateInHousePreconditionerPCD", "finish" );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateDefinePressureCst()
{
    M_definePressureCstOnlyOneZoneAppliedOnWholeMesh = M_materialProperties->isDefinedOnWholeMesh();
    M_definePressureCstMeshRanges.clear();
    for ( auto const& markers : M_definePressureCstMarkers )
    {
        for ( std::string const& marker : markers )
            CHECK( M_mesh->hasElementMarker( marker ) ) << "marker " << marker << "does not found in mesh";
        M_definePressureCstMeshRanges.push_back( markedelements(M_mesh,markers) );
    }
    if ( M_definePressureCstMeshRanges.empty() )
        M_definePressureCstMeshRanges.push_back( M_rangeMeshElements );
    else
        M_definePressureCstOnlyOneZoneAppliedOnWholeMesh = false;



    if ( this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_XhMeanPressureLM.resize( M_definePressureCstMeshRanges.size() );
        if ( M_definePressureCstOnlyOneZoneAppliedOnWholeMesh )
            M_XhMeanPressureLM[0] = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
        else
        {
            for ( int k=0;k<M_definePressureCstMeshRanges.size();++k )
                M_XhMeanPressureLM[k] = space_meanpressurelm_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm(),
                                                                        _range=M_definePressureCstMeshRanges[k] );
        }
    }
    else if ( this->definePressureCstMethod() == "algebraic" )
    {
        auto p = this->functionSpacePressure()->element();

        M_definePressureCstAlgebraicOperatorMeanPressure.resize(M_definePressureCstMeshRanges.size());
        auto dofTablePressure = this->functionSpacePressure()->dof();
        for ( int k=0;k<M_definePressureCstMeshRanges.size();++k )
        {
            auto const& rangeElt = M_definePressureCstMeshRanges[k];
            M_definePressureCstAlgebraicOperatorMeanPressure[k].first = form1_mean(_test=this->functionSpacePressure(),
                                                                                   _range=rangeElt,
                                                                                   _expr=id(p) ).vectorPtr();
            M_definePressureCstAlgebraicOperatorMeanPressure[k].first->close();
            auto & dofsOnRange = M_definePressureCstAlgebraicOperatorMeanPressure[k].second;
            for ( auto const& elt : rangeElt )
            {
                for( auto const& ldof : dofTablePressure->localDof( unwrap_ref( elt ).id() ) )
                    dofsOnRange.insert( ldof.second.index() );
            }
        }
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateVorticity()
{

    if (!M_XhVorticity) this->createFunctionSpacesVorticity();

    detail::updateVorticityImpl( this->fieldVelocity(),*M_fieldVorticity, M_rangeMeshElements, this->geomap(), mpl::int_<nDim>() );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateConvectionVelocityExtrapolated()
{
    if ( this->isStationary() )
        return;

    if ( !M_fieldConvectionVelocityExtrapolated )
        M_fieldConvectionVelocityExtrapolated = this->functionSpace()->elementPtr();
    else
        M_fieldConvectionVelocityExtrapolated->zero();

    if ( this->currentTime() == this->timeInitial() )
    {
        if ( M_bdf_fluid->iteration() == 1 )
            M_fieldConvectionVelocityExtrapolated->add( 1, M_bdf_fluid->unknown(1) );
        else if ( M_bdf_fluid->iteration() > 1 )
        {
            M_fieldConvectionVelocityExtrapolated->add( 2, M_bdf_fluid->unknown(1) );
            M_fieldConvectionVelocityExtrapolated->add( -1, M_bdf_fluid->unknown(2) );
        }
    }
    else
    {
        if ( M_timeStepping == "BDF" )
        {
            *M_fieldConvectionVelocityExtrapolated = *M_bdf_fluid->polyPtr();
        }
        else if ( M_timeStepping == "Theta" )
        {
            if ( M_bdf_fluid->iteration() == 1 )
                M_fieldConvectionVelocityExtrapolated->add( 1, M_bdf_fluid->unknown(0) );
            else if ( M_bdf_fluid->iteration() > 1 )
            {
                M_fieldConvectionVelocityExtrapolated->add( 2 /*3./2.*/, M_bdf_fluid->unknown(0) );
                M_fieldConvectionVelocityExtrapolated->add( -1 /*-1./2.*/, M_bdf_fluid->unknown(1) );
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("FluidMechanics","startTimeStep", "start" );

    if ( this->solverName() == "Oseen" )
        this->updateConvectionVelocityExtrapolated();
    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    // start time step
    if ( !this->doRestart() )
        M_bdf_fluid->start(*M_Solution);
    // up current time
    this->updateTime( M_bdf_fluid->time() );

    if ( this->solverName() == "Oseen" )
        this->updateConvectionVelocityExtrapolated();

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    this->log("FluidMechanics","startTimeStep", "finish" );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("FluidMechanics","updateTimeStep", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

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

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->isMoveDomain())
        M_meshALE->updateTimeStep();
#endif

    bool rebuildCstAssembly = false;
    if ( M_timeStepping == "BDF" )
    {
        int previousTimeOrder = this->timeStepBDF()->timeOrder();
        M_bdf_fluid->next( *M_Solution );
        int currentTimeOrder = this->timeStepBDF()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
        this->updateTime( M_bdf_fluid->time() );
    }
    else if ( M_timeStepping == "Theta" )
    {
        M_bdf_fluid->next( *M_Solution );
        this->updateTime( M_bdf_fluid->time() );
    }

    if ( this->solverName() == "Oseen" )
        this->updateConvectionVelocityExtrapolated();

    // update user functions which depend of time only
    this->updateUserFunctions(true);
    // update all expressions in bc or in house prec
    this->updateParameterValues();

    // maybe rebuild cst jacobian or linear
    if ( M_algebraicFactory && rebuildCstAssembly )
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


    this->timerTool("TimeStepping").stop("updateTimeStep");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("FluidMechanics","updateTimeStep", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( this->isStationaryModel() )
        return;
    if ( !M_algebraicFactory )
        return;
    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        M_blockVectorSolution.updateVectorFromSubVectors();
        //this->updateBlockVectorSolution();
        ModelAlgebraic::DataUpdateResidual dataResidual( M_blockVectorSolution.vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );

        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        M_algebraicFactory->evaluateResidual( dataResidual );
        M_algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_stabilizationGLS )
        {
            auto & dataInfos = M_algebraicFactory->dataInfos();
            *dataInfos.vectorInfo( prefixvm( this->prefix(),"time-stepping.previous-solution") ) = *M_blockVectorSolution.vectorMonolithic();
            if ( this->solverName() == "Oseen" )
            {
                std::string convectionOseenEntry = prefixvm( this->prefix(),"time-stepping.previous-convection-field-extrapolated" );
                if ( !dataInfos.hasVectorInfo( convectionOseenEntry ) )
                    dataInfos.addVectorInfo( convectionOseenEntry, this->backend()->newVector( M_fieldConvectionVelocityExtrapolated->mapPtr() ) );
                *dataInfos.vectorInfo( convectionOseenEntry ) = *M_fieldConvectionVelocityExtrapolated;
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNormalStressOnCurrentMesh( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate )
{
    this->log("FluidMechanics","updateNormalStressOnCurrentMesh", "start" );
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    CHECK( M_rangeDistributionByMaterialName ) << "M_rangeDistributionByMaterialName is not init";
    for ( auto const& rangeFacesMat : M_rangeDistributionByMaterialName->rangeMeshFacesByMaterial( nameOfRange ) )
    {
        std::string matName = rangeFacesMat.first;
        auto const& rangeFaces = rangeFacesMat.second;
        auto const sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,true);
        fieldToUpdate->on(_range=rangeFaces,
                          _expr=sigmav*N(),
                          _geomap=this->geomap() );
    }
    this->log("FluidMechanics","updateNormalStressOnCurrentMesh", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateWallShearStress( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate )
{
    this->log("FluidMechanics","updateWallShearStress", "start" );
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    CHECK( M_rangeDistributionByMaterialName ) << "M_rangeDistributionByMaterialName is not init";
    for ( auto const& rangeFacesMat : M_rangeDistributionByMaterialName->rangeMeshFacesByMaterial( nameOfRange ) )
    {
        std::string matName = rangeFacesMat.first;
        auto const& rangeFaces = rangeFacesMat.second;
        auto const sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,true);
        fieldToUpdate->on(_range=rangeFaces,
                          _expr=sigmav*vf::N() - (trans(sigmav*vf::N())*vf::N())*vf::N(),
                          _geomap=this->geomap() );
    }
    this->log("FluidMechanics","updateWallShearStress", "finish" );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNormalStressOnReferenceMesh( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate )
{
    bool meshIsOnRefAtBegin = this->meshALE()->isOnReferenceMesh();
    if ( !meshIsOnRefAtBegin )
        this->meshALE()->revertReferenceMesh( false );

    // current solution
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());

    //auto itFindRange = M_rangeMeshFacesByMaterial.find( "moving-boundary" );
    //CHECK( itFindRange != M_rangeMeshFacesByMaterial.end() ) << "not find range moving-boundary";
    CHECK( M_rangeDistributionByMaterialName ) << "M_rangeDistributionByMaterialName is not init";
    for ( auto const& rangeFacesMat : M_rangeDistributionByMaterialName->rangeMeshFacesByMaterial( nameOfRange ) )
    {
        std::string matName = rangeFacesMat.first;
        auto const& rangeFaces = rangeFacesMat.second;
        // stress tensor : -p*Id + 2*mu*D(u)
        auto const sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u)*inv(Fa),idv(p),*this->materialProperties(),matName,true);
        fieldToUpdate->on(_range=rangeFaces,
                          _expr=sigmav*det(Fa)*trans(inv(Fa))*N(),
                          _geomap=this->geomap() );

#if 0
        //
        auto resRef = integrate(_range=rangeFaces,_expr=idv(fieldToUpdate)).evaluate();
        this->meshALE()->revertMovingMesh( false );
        // current solution
        auto const& u = this->fieldVelocity();
        auto const& p = this->fieldPressure();
        // identity Matrix
        auto const Id = eye<nDim,nDim>();
        // deformation tensor
        auto Fa = Id+gradv(*M_meshALE->displacement());
        auto InvFa = det(Fa)*inv(Fa);
        auto const Sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,true);
        auto resMove = integrate(_range=rangeFaces,_expr=Sigmav*N() ).evaluate();
        this->meshALE()->revertReferenceMesh( false );
        if ( this->worldComm().isMasterRank() )
        {
            std::cout << "resRef " << resRef << "\n";
            std::cout << "resMove " << resMove << "\n";
            std::cout << "diff " << (resRef-resMove).norm() << "\n";
        }
#endif
    }

    if ( !meshIsOnRefAtBegin )
        this->meshALE()->revertMovingMesh( false );
}

//---------------------------------------------------------------------------------------------------------//


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updatePicardConvergence( vector_ptrtype const& Unew, vector_ptrtype const& Uold ) const
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateALEmesh()
{
    this->log("FluidMechanics","updateALEmesh", "start");

    if ( !M_bcMovingBoundaryImposed.empty() )
    {
        // Warning : evaluate expression on reference mesh (maybe it will better to change the API in order to avoid tjese meshmoves)
        bool meshIsOnRefAtBegin = this->meshALE()->isOnReferenceMesh();
        if ( !meshIsOnRefAtBegin )
            this->meshALE()->revertReferenceMesh( false );
        for( auto const& d : M_bcMovingBoundaryImposed )
        {
            M_meshDisplacementOnInterface->on( _range=markedfaces(this->mesh(),markers(d)),
                                               _expr=expression(d,this->symbolsExpr()),
                                               _geomap=this->geomap() );
        }
        if ( !meshIsOnRefAtBegin )
            this->meshALE()->revertMovingMesh( false );
    }
    //-------------------------------------------------------------------//
    // compute ALE map
    //std::vector< mesh_ale_type::ale_map_element_type> polyBoundarySet = { *M_meshDisplacementOnInterface };
    M_meshALE->update(*M_meshDisplacementOnInterface/*polyBoundarySet*/);

    //-------------------------------------------------------------------//

    if ( this->doCIPStabConvection() )
    {
        M_fieldMeshVelocityUsedWithStabCIP->on(_range=M_rangeMeshElements,_expr=idv(M_meshALE->velocity()) );
        sync( *M_fieldMeshVelocityUsedWithStabCIP, "=" );
    }
    //-------------------------------------------------------------------//
    // update operator PCD
    if ( M_preconditionerAttachPCD && this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasOperatorPCD("pcd") )
    {
        this->log("FluidMechanics","updateALEmesh", "rebuild operatorPCD");
        CHECK( this->algebraicFactory()->preconditionerTool()->hasOperatorPCD("pcd") ) << "operator PCD does not init";
        typedef Feel::Alternatives::OperatorPCD<space_fluid_velocity_type,space_fluid_pressure_type> op_pcd_type;
        std::shared_ptr<op_pcd_type> myOpPCD =
            std::dynamic_pointer_cast<op_pcd_type>( this->algebraicFactory()->preconditionerTool()->operatorPCD( "pcd" ) );
        myOpPCD->assemble();
    }

    this->log("FluidMechanics","updateALEmesh", "finish");
}
#endif

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeMeshArea( std::string const& marker ) const
{
    return this->computeMeshArea( std::set<std::string>( { marker } ) );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeMeshArea( std::set<std::string> const& markers ) const
{
    double area = 0;
    if ( markers.empty() || markers.begin()->empty() )
        area = integrate(_range=M_rangeMeshElements,//elements(this->mesh()),
                         _expr=cst(1.),
                         _geomap=this->geomap() ).evaluate()(0,0);
    else
        area = integrate(_range=markedelements(this->mesh(),markers),
                         _expr=cst(1.),
                         _geomap=this->geomap() ).evaluate()(0,0);
    return area;
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::force_type
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeForce(std::string const& markerName) const
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
    auto Sigmav = val(-idv(p)*Id + 2*idv(this->materialProperties()->fieldMu())*defv);
#endif
    CHECK( this->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = this->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,true);

    return integrate(_range=markedfaces(M_mesh,markerName),
                     _expr= sigmav*N(),
                     _geomap=this->geomap() ).evaluate();

}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeFlowRate( std::string const& marker, bool useExteriorNormal ) const
{
    return this->computeFlowRate( std::list<std::string>( { marker } ),useExteriorNormal );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeFlowRate( std::list<std::string> const& markers, bool useExteriorNormal ) const
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computePressureSum() const
{
    auto const& p = this->fieldPressure();
    double res = integrate(_range=M_rangeMeshElements,
                           _expr= idv(p),
                           _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computePressureMean() const
{
    double area = this->computeMeshArea();
    double res = (1./area)*this->computePressureSum();
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceSum() const
{
    auto const& u = this->fieldVelocity();
    double res = integrate(_range=M_rangeMeshElements,
                           _expr= divv(u),
                           _geomap=this->geomap() ).evaluate()(0,0);
    return res;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceMean() const
{
    double area = this->computeMeshArea();
    double res = (1./area)*this->computeVelocityDivergenceSum();
    return res;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeVelocityDivergenceNormL2() const
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::vector<double>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeAveragedPreassure( std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
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
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::vector<double>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeFlowRate(std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::super_type::block_pattern_type
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::blockPattern() const
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
int
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        nBlock+=M_XhMeanPressureLM.size();
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
    return nBlock;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
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

        for ( int k=0;k<M_XhMeanPressureLM.size();++k )
        {
            myblockGraph(indexBlock,0) = stencil(_test=M_XhMeanPressureLM[k],_trial=this->functionSpace(),
                                                 _pattern_block=patCouplingLM,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(0,indexBlock) = stencil(_test=this->functionSpace(),_trial=M_XhMeanPressureLM[k],
                                                 _pattern_block=patCouplingLM.transpose(),
                                                 _diag_is_nonzero=false,_close=false)->graph();
            ++indexBlock;
        }
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

    this->log("FluidMechanics","buildBlockMatrixGraph", "finish" );

    return myblockGraph;
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::graph_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    auto blockGraph = this->buildBlockMatrixGraph();
    blockGraph.close();

    if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
        return blockGraph(0,0);
    else
        return graph_ptrtype( new graph_type( blockGraph ) );
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::indexsplit_ptrtype
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildIndexSplit() const
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
size_type
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    auto res = this->functionSpace()->nLocalDofWithGhost();
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        for ( int k=0;k<M_XhMeanPressureLM.size();++k )
            res += M_XhMeanPressureLM[k]->nLocalDofWithGhost();
    }
    if (this->hasMarkerDirichletBClm())
        res += this->XhDirichletLM()->nLocalDofWithGhost();
    if ( this->hasFluidOutletWindkesselImplicit() )
        res += 2*this->nFluidOutletWindkesselImplicit();
    return res;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBlockVectorSolution()
{
    // copy velocity/pressure in block
#if 0
    auto & vecAlgebraic = M_blockVectorSolution.vector();
    auto const& fieldVelPres = this->fieldVelocityPressure();
    for (int k=0;k< this->functionSpace()->nLocalDofWithGhost() ;++k)
        vecAlgebraic->set(k, fieldVelPres(k) );
#else
    //M_blockVectorSolution.vector()->close();
    M_blockVectorSolution.setVector( *M_blockVectorSolution.vectorMonolithic(), this->fieldVelocityPressure(), 0 );
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBoundaryConditionsForUse()
{
    auto XhVelocity = this->functionSpaceVelocity();
    auto XhCompVelocity = XhVelocity->compSpace();
    static const uint16_type nDofComponentsVelocity = XhVelocity->dof()->nDofComponents();
    auto mesh = this->mesh();

    std::set<std::string> velocityMarkers;
    std::map<ComponentType,std::set<std::string> > compVelocityMarkers;

    //-------------------------------------//
    // strong Dirichlet bc on velocity from expression
    for( auto const& d : M_bcDirichlet )
    {
        auto listMark = this->markerDirichletBCByNameId( "elimination",name(d) );
        velocityMarkers.insert( listMark.begin(), listMark.end() );
    }
    // strong Dirichlet bc on velocity component from expression
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto listMark = this->markerDirichletBCByNameId( "elimination",name(d), comp );
            compVelocityMarkers[comp].insert( listMark.begin(), listMark.end() );
        }
    }
    // strong Dirichlet bc on velocity from inlet bc
    for ( auto const& inletbc : M_fluidInletDesc )
    {
        std::string const& marker = std::get<0>( inletbc );
        velocityMarkers.insert( marker );
    }
    // strong  Dirichlet bc on velocity from moving boundary imposed
    for( auto const& d : M_bcMovingBoundaryImposed )
    {
        auto listMark = M_bcMarkersMovingBoundaryImposed.markerDirichletBCByNameId( "elimination",name(d) );
        velocityMarkers.insert( listMark.begin(), listMark.end() );
    }

    //-------------------------------------//
    // distribute mesh markers by entity
    std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > meshMarkersVelocityByEntities;
    std::map<ComponentType, std::tuple< std::set<std::string>,std::set<std::string>,std::set<std::string>,std::set<std::string> > > meshMarkersCompVelocityByEntities;
    meshMarkersVelocityByEntities = detail::distributeMarkerListOnSubEntity( mesh, velocityMarkers );
    for ( auto const& compMarkerPair : compVelocityMarkers )
    {
        meshMarkersCompVelocityByEntities[compMarkerPair.first] = detail::distributeMarkerListOnSubEntity( mesh, compMarkerPair.second );
    }
    //-------------------------------------//
    // on topological faces
    auto const& listMarkedFacesVelocity = std::get<0>( meshMarkersVelocityByEntities );
    if ( !listMarkedFacesVelocity.empty() )
    {
        auto therange = markedfaces( mesh,listMarkedFacesVelocity );
        auto dofsToAdd = XhVelocity->dofs( therange );
        XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("velocity",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("velocity",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
    // on marked edges (only 3d)
    auto const& listMarkedEdgesVelocity = std::get<1>( meshMarkersVelocityByEntities );
    if ( !listMarkedEdgesVelocity.empty() )
    {
        auto therange = markededges(mesh,listMarkedEdgesVelocity );
        auto dofsToAdd = XhVelocity->dofs( therange );
        XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("velocity",MESH_EDGES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("velocity",MESH_EDGES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
    // on marked points
    auto const& listMarkedPointsVelocity = std::get<2>( meshMarkersVelocityByEntities );
    if ( !listMarkedPointsVelocity.empty() )
    {
        auto therange = markedpoints(mesh,listMarkedPointsVelocity );
        auto dofsToAdd = XhVelocity->dofs( therange );
        XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("velocity",MESH_POINTS).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("velocity",MESH_POINTS).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
    //-------------------------------------//
    // on velocity components
    for ( auto const& meshMarkersPair : meshMarkersCompVelocityByEntities )
    {
        ComponentType comp = meshMarkersPair.first;
        auto const& listMarkedFacesCompVelocity = std::get<0>( meshMarkersPair.second );
        if ( !listMarkedFacesCompVelocity.empty() )
        {
            auto therange = markedfaces(mesh,listMarkedFacesCompVelocity );
            auto dofsToAdd = XhVelocity->dofs( therange, comp );
            XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("velocity",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("velocity",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
        // edges (only 3d)
        auto const& listMarkedEdgesCompVelocity = std::get<1>( meshMarkersPair.second );
        if ( !listMarkedEdgesCompVelocity.empty() )
        {
            auto therange = markededges(mesh,listMarkedEdgesCompVelocity );
            auto dofsToAdd = XhVelocity->dofs( therange, comp );
            XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("velocity",MESH_EDGES).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("velocity",MESH_EDGES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
        // points
        auto const& listMarkedPointsCompVelocity = std::get<2>( meshMarkersPair.second );
        if ( !listMarkedPointsCompVelocity.empty() )
        {
            auto therange = markedpoints(mesh,listMarkedPointsCompVelocity );
            auto dofsToAdd = XhVelocity->dofs( therange, comp );
            XhVelocity->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
            this->dofEliminationIdsAll("velocity",MESH_POINTS).insert( dofsToAdd.begin(), dofsToAdd.end() );
            auto dofsMultiProcessToAdd = XhVelocity->dofs( therange, comp, true );
            this->dofEliminationIdsMultiProcess("velocity",MESH_POINTS).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
        }
    }

    if ( this->hasMarkerPressureBC() && M_spaceLagrangeMultiplierPressureBC )
    {
        auto therange = boundaryfaces(M_meshLagrangeMultiplierPressureBC);
        auto dofsToAdd = M_spaceLagrangeMultiplierPressureBC->dofs( therange );
        M_spaceLagrangeMultiplierPressureBC->dof()->updateIndexSetWithParallelMissingDof( dofsToAdd );
        this->dofEliminationIdsAll("pressurebc-lm",MESH_FACES).insert( dofsToAdd.begin(), dofsToAdd.end() );
        auto dofsMultiProcessToAdd = M_spaceLagrangeMultiplierPressureBC->dofs( therange, ComponentType::NO_COMPONENT, true );
        this->dofEliminationIdsMultiProcess("pressurebc-lm",MESH_FACES).insert( dofsMultiProcessToAdd.begin(), dofsMultiProcessToAdd.end() );
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateRangeDistributionByMaterialName( std::string const& key, range_faces_type const& rangeFaces )
{
    if ( !M_rangeDistributionByMaterialName )
    {
        M_rangeDistributionByMaterialName = std::make_shared<RangeDistributionByMaterialName<mesh_type>>();
        M_rangeDistributionByMaterialName->init( this->materialProperties()->rangeMeshElementsByMaterial() );
    }
    M_rangeDistributionByMaterialName->update( key, rangeFaces );
}

} // namespace FeelModels
} // namespace Feel


