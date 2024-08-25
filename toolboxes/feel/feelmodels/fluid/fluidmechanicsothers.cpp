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

//#include <feel/feelpde/operatorpcd.hpp>

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
    if ( this->worldComm().localSize() == 1 )
        return false;
    return ( M_XhVelocity->extendedDofTable() || M_XhPressure->extendedDofTable() );
}

//---------------------------------------------------------------------------------------------------------//
#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::getInfo() const
{
    std::string ALEmode;if ( this->hasMeshMotion() ) ALEmode="Yes"; else ALEmode="No";
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
    //if (this->doStabDivDiv()) stabAll_str=(stabAll_str.empty())?"DivDiv":stabAll_str+" - DivDiv";
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
    for ( std::string const& fieldName : this->postProcessExportsFields() )
        doExport_str=(doExport_str.empty())? fieldName : doExport_str + " - " + fieldName;
    // for ( std::string const& fieldName : M_postProcessUserFieldExported )
    //     doExport_str=(doExport_str.empty())? fieldName : doExport_str + " - " + fieldName;

    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : FluidMechanics---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository();
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        *_ostr << "\n   Physical Model [" << physicName << "]"
               << "\n     -- equation : " << physicFluidData->equation()
               << "\n     -- time mode : " << StateTemporal
               << "\n     -- ale mode  : " << ALEmode
               << "\n     -- gravity force enabled  : " << physicFluidData->gravityForceEnabled();
        if ( physicFluidData->gravityForceEnabled() )
             *_ostr << "\n     -- gravity force expr  : " << str( physicFluidData->gravityForceExpr().expression() );
    }
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
#if 0
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
#endif
#if defined( FEELPP_MODELS_HAS_MESHALE )
    // if ( this->isMoveDomain() )
    // *_ostr << this->getInfoALEMeshBC();
#endif
    *_ostr << "\n   Space Discretization";
#if 0
    if ( this->hasGeoFile() )
        *_ostr << "\n     -- geo file name   : " << this->geoFile();
    *_ostr << "\n     -- mesh file name   : " << this->meshFile()
           << "\n     -- nb elt in mesh  : " << M_mesh->numGlobalElements()//numElements()
        // << "\n     -- nb elt in mesh  : " << M_mesh->numElements()
        // << "\n     -- nb face in mesh : " << M_mesh->numFaces()
           << "\n     -- hMin            : " << M_mesh->hMin()
           << "\n     -- hMax            : " << M_mesh->hMax()
           << "\n     -- hAverage        : " << M_mesh->hAverage()
           << "\n     -- geometry order  : " << nOrderGeo;
#endif
    *_ostr << "\n     -- velocity order  : " << nOrderVelocity
           << "\n     -- pressure order  : " << nOrderPressure
           << "\n     -- nb dof (u)      : " << M_XhVelocity->nDof() << " (" << M_XhVelocity->nLocalDof() << ")"
           << "\n     -- nb dof (p)      : " << M_XhPressure->nDof() << " (" << M_XhPressure->nLocalDof() << ")"
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
           << "\n     -- number of proc environment : " << Environment::worldComm().globalSize()
           << "\n     -- environment rank : " << Environment::worldComm().rank()
           << "\n     -- global rank : " << this->worldComm().globalRank()
           << "\n     -- local rank : " << this->worldComm().localRank()
        //<< "\n   Matrix Index Start"
        //  << "\n     -- rowstart : " << this->rowStartInMatrix()
        //  << "\n     -- colstart : " << this->colStartInMatrix()
           << "\n   Numerical Solver"
           << "\n     -- scheme : " << (this->useSemiImplicitTimeScheme()? std::string("semi-implicit"):std::string("implicit"))
           << "\n     -- solver : " << M_solverName;
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
// #if defined( FEELPP_MODELS_HAS_MESHALE )
//     if ( this->isMoveDomain() )
//         *_ostr << this->meshALE()->getInfo()->str();
// #endif
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    // if ( this->worldComm().isMasterRank() )
    //     std::cout << "symbolsExpr : \n "<<  this->symbolsExpr().names() << std::endl;

    return _ostr;
}
#endif

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;

    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    nl::json subPt;
    subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p["Physics2"] = subPt;

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    // Boundary Conditions
    M_boundaryConditions->updateInformationObject( p["Boundary Conditions"] );
#if 0 // VINCENT
    if ( !M_bodySetBC.empty() )
        M_bodySetBC.updateInformationObject( p["Boundary Conditions"] );
#endif

    // FunctionSpace
    subPt.clear();
    subPt["Velocity"] = this->functionSpaceVelocity()->journalSection().to_string();
    subPt["Pressure"] = this->functionSpacePressure()->journalSection().to_string();
    p["Function Spaces"] = subPt;

    this->modelFields().updateInformationObject( p["Fields"] );

    if ( this->algebraicFactory() )
        this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New( tabInfoProp );

    // Environment
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    // Physics
    if ( jsonInfo.contains("Physics2") )
    {
        Feel::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics2"), tabInfoProp );
        tabInfo->add( "Physics2", TabulateInformations::New( tabInfoPhysics, tabInfoProp ) );
    }

    // Materials Properties
    if ( this->materialsProperties() && jsonInfo.contains("Materials Properties") )
        tabInfo->add( "Materials Properties", this->materialsProperties()->tabulateInformations(jsonInfo.at("Materials Properties"), tabInfoProp ) );

    // Boundary conditions
    if ( jsonInfo.contains("Boundary Conditions") )
        tabInfo->add( "Boundary Conditions", boundary_conditions_type::tabulateInformations( jsonInfo.at("Boundary Conditions"), tabInfoProp ) );

    // Meshes
    if ( jsonInfo.contains("Meshes") )
        tabInfo->add( "Meshes", super_type::super_model_meshes_type::tabulateInformations( jsonInfo.at("Meshes"), tabInfoProp ) );

    // Function Spaces
    auto tabInfoFunctionSpaces = TabulateInformationsSections::New( tabInfoProp );
    tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        for ( std::string const& spaceName : std::vector<std::string>({"Velocity","Pressure"}) )
        {
            nl::json::json_pointer jsonPointerSpace( jsonInfoFunctionSpaces.at( spaceName ).template get<std::string>() );
            if ( JournalManager::journalData().contains( jsonPointerSpace ) )
                tabInfoFunctionSpaces->add( spaceName, TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpace ), tabInfoProp ) );
        }
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );


    std::map<std::string,uint16_type> jsonPtrFunctionSpacesToLevel;
#if 0 // VINCENT
    if ( jsonInfo.contains("Boundary Conditions") )
    {
        auto const& jsonInfoBoundaryConditions = jsonInfo.at("Boundary Conditions");
        auto tabInfoBoundaryConditions = TabulateInformationsSections::New( tabInfoProp );

        if ( !M_bodySetBC.empty() )
            M_bodySetBC.updateTabulateInformations( tabInfoBoundaryConditions, jsonInfoBoundaryConditions, tabInfoProp, jsonPtrFunctionSpacesToLevel );
        //tabInfoBoundaryConditions->add( "Body Set Boundary Condition", M_bodySetBC.tabulateInformations( jsonInfoBoundaryConditions, tabInfoProp ) );

        tabInfo->add( "Boundary Conditions", tabInfoBoundaryConditions );
    }
#endif

    for ( auto const& [jsonPtrStr,level] :jsonPtrFunctionSpacesToLevel )
    {
        nl::json::json_pointer jsonPointerSpace( jsonPtrStr );
        if ( JournalManager::journalData().contains( jsonPointerSpace ) )
            tabInfoFunctionSpaces->add( /*spaceName*/ jsonPtrStr, TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( JournalManager::journalData().at( jsonPointerSpace ), TabulateInformationProperties(level) ) );
    }

    // Algebraic Solver
    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    return tabInfo;
}


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setSolverName( std::string const& type )
{
    // if solver change -> force to rebuild all algebraic data at next solve
    if ( type != M_solverName )
        this->setNeedToRebuildCstPart(true);

    M_solverName = type;
    if ( type == "Oseen" )
    {
        M_useSemiImplicitTimeScheme = true;
        M_solverName="LinearSystem";
    }
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
#if 0
    if( this->modelName() == "Stokes" )
        return true;
    else
        return this->isStationary();
#else

// TODO : should be detected in loop assembly (only works with one physic)
bool isStokes = true;
for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
 {
      auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
      if ( physicFluidData->equation() != "Stokes" )
      {
          isStokes = false;
          break;
      }
}
if ( isStokes )
    return true;
 else
     return this->isStationary();
#endif
}

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
    this->exportResults( time, this->symbolsExpr() );
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
    if ( this->hasPostProcessExportsField( "pid" ) )
    {
        M_exporter_ho->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessExportsField( "velocity" ) )
    {
        M_opIvelocity->apply( this->fieldVelocity(),*M_velocityVisuHO);
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity_ho")), *M_velocityVisuHO );
        hasFieldToExport = true;
    }
    if ( this->hasPostProcessExportsField( "pressure" ) )
    {
        M_opIpressure->apply( this->fieldPressure(),*M_pressureVisuHO);
        M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure_ho"), prefixvm(this->prefix(),prefixvm(this->subPrefix(),"pressure_ho")), *M_pressureVisuHO );
        hasFieldToExport = true;
    }

    if ( this->hasMeshMotion() )
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
#if 0
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
#endif
        if ( this->hasPostProcessExportsField( "displacement" ) )
        {
            //M_opImeshdisp->apply( thedisp , *M_meshdispVisuHO);
            M_opImeshdisp->apply( *this->meshMotionTool()->displacement() , *M_meshdispVisuHO);
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

#if 0
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
#endif

#if 0
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::executePostProcessMeasures( double time )
{
    auto mfields = this->modelFields();
    this->executePostProcessMeasures( time, mfields, this->symbolsExpr( mfields ) );
}
#endif

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateParameterValues()
{
    if ( !this->manageParameterValues() )
        return;

    this->modelProperties().parameters().updateParameterValues();
    auto paramValues = this->modelProperties().parameters().toParameterValues();
    this->materialsProperties()->updateParameterValues( paramValues );
    for ( auto [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->updateParameterValues( paramValues );

    M_bodySetBC.updateParameterValues( paramValues );
#if 0
    std::cout << "parameter values\n";
    for ( auto const& [name,val] : paramValues )
        std::cout << name << " : " << val << std::endl;
#endif
    this->setParameterValues( paramValues );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    for ( auto const& [param,val] : paramValues )
        M_currentParameterValues[param] = val;

    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->modelProperties().initialConditions().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );


    super_type::super_model_meshes_type::setParameterValues( paramValues );

    M_boundaryConditions->setParameterValues( paramValues );
    M_bodySetBC.setParameterValues( paramValues );

#if 0 // VINCENT
    M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    M_bcNeumannScalar.setParameterValues( paramValues );
    M_bcNeumannVectorial.setParameterValues( paramValues );
    M_bcNeumannTensor2.setParameterValues( paramValues );
    M_bcPressure.setParameterValues( paramValues );
    M_bcMovingBoundaryImposed.setParameterValues( paramValues );
    M_volumicForcesProperties.setParameterValues( paramValues );
    this->updateFluidInletVelocity();
#endif

    if ( this->hasTurbulenceModel() )
        M_turbulenceModelType->setParameterValues( paramValues );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::solve()
{
    this->log("FluidMechanics","solve", "start" );
    this->timerTool("Solve").start();

    this->setStartBlockSpaceIndex( 0 );

    auto se = this->symbolsExpr();
    this->updateFluidInletVelocity( se ); // TODO VINCENT : create an updateBoundaryConditionForUse?

    // copy velocity/pressure in algebraic vector solution (maybe velocity/pressure has been changed externally)
    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
#if 0
    if ( this->worldComm().isMasterRank() )
        std::cout << "symbolsExpr : \n "<<  se.names() << std::endl;
#endif
    if ( this->hasMeshMotion() && M_applyMovingMeshBeforeSolve )
        this->updateALEmesh();

    M_bodySetBC.updateForUse( *this );
    M_bodySetBC.updateAlgebraicFactoryForUse( *this, this->algebraicFactory() );
#if 0 // TODO
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
#endif
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

#if 0
    if ( this->hasTurbulenceModel() && M_useSemiImplicitTurbulenceCoupling )
        M_turbulenceModelType->solve();
#endif
    //--------------------------------------------------
    // run solver
    std::string algebraicSolver = M_solverName;
    this->algebraicFactory()->solve( algebraicSolver, this->algebraicBlockVectorSolution()->vectorMonolithic() );
    // update sub vector
    this->algebraicBlockVectorSolution()->localize();

#if 0
    //--------------------------------------------------
    // update windkessel solution ( todo put fluidOutletWindkesselPressureDistal and Proximal in blockVectorSolution )
    int cptBlock=1;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        ++cptBlock;

    if ( M_boundaryConditions->hasVelocityImposedLagrangeMultiplier() )
        ++cptBlock;

    if ( !M_boundaryConditions->pressureImposed().empty() )
    {
        ++cptBlock;
        if ( nDim == 3 )
            ++cptBlock;
    }
#if 0 // VINCENT
    if (this->hasFluidOutletWindkesselImplicit() )
    {
        for (int k=0;k<this->nFluidOutlet();++k)
        {
            if ( std::get<1>( M_fluidOutletsBCType[k] ) != "windkessel" || std::get<0>( std::get<2>( M_fluidOutletsBCType[k] ) ) != "implicit" )
                continue;

            if ( M_fluidOutletWindkesselSpace->nLocalDofWithoutGhost() > 0 )
            {
                M_fluidOutletWindkesselPressureDistal[k] = this->algebraicBlockVectorSolution()->operator()(cptBlock)->operator()(0);
                M_fluidOutletWindkesselPressureProximal[k] = this->algebraicBlockVectorSolution()->operator()(cptBlock)->operator()(1);
            }
            ++cptBlock;
        }
    }
#endif
#endif
    //--------------------------------------------------
#if 0
    for ( auto const& [bpname,bbc] : M_bodySetBC )
    {
        if ( bbc.spaceTranslationalVelocity()->nLocalDofWithGhost() )
            std::cout << "["<< this->worldComm().localRank() << "] NAME="<<bbc.name()
                      << " : " << bbc.fieldTranslationalVelocityPtr()->operator()(0)
                      << " , " << bbc.fieldTranslationalVelocityPtr()->operator()(1)
                      << " , " << bbc.fieldTranslationalVelocityPtr()->operator()(2)
                      << std::endl;
    }
#endif

    if ( this->hasTurbulenceModel() && M_useSemiImplicitTurbulenceCoupling )
    {
        for (int k=0;k<ioption(_name="solver-inner-turbulence.nit");++k)
            M_turbulenceModelType->solve();
    }

    // if ( this->hasTurbulenceModel() )
    //     M_turbulenceModelType->solve();


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
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonIteration( int step, vector_ptrtype residual, vector_ptrtype sol,
                                                           typename backend_type::solvernonlinear_type::UpdateIterationData const& data ) const
{
    if ( this->hasTurbulenceModel() && !M_useSemiImplicitTurbulenceCoupling )
    {
        M_turbulenceModelType->algebraicFactory()->dataInfos().addVectorInfo( prefixvm(this->prefix(), "current-solution"), sol );
        M_turbulenceModelType->solve();
    }
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updatePicardIteration( int step, vector_ptrtype sol ) const
{
    if ( this->hasTurbulenceModel() && !M_useSemiImplicitTurbulenceCoupling )
        M_turbulenceModelType->solve();
}

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
        auto pSol = this->functionSpacePressure()->element( sol, this->rowStartInVector()+1 );
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
{
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const
{
    this->postSolveNewton( rhs,sol );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const
{

}

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
    if ( !M_preconditionerAttachPMM && !M_preconditionerAttachPCD )
        return;
    vector_ptrtype const& vecSol = data.currentSolution();
    this->updateInHousePreconditioner( data, this->modelContext( vecSol, this->rowStartInVector() ) );
}
FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInHousePreconditioner( DataUpdateJacobian & data ) const
{
    if ( !M_preconditionerAttachPMM && !M_preconditionerAttachPCD )
        return;
    vector_ptrtype const& vecSol = data.currentSolution();
    this->updateInHousePreconditioner( data, this->modelContext( vecSol, this->rowStartInVector() ) );
}


FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateDefinePressureCst()
{
    auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
    M_definePressureCstOnlyOneZoneAppliedOnWholeMesh = mom->isDefinedOnWholeMesh( this->physicsAvailableFromCurrentType() );
    M_definePressureCstMeshRanges.clear();
    for ( auto const& markers : M_definePressureCstMarkers )
    {
        for ( std::string const& marker : markers )
            CHECK( this->mesh()->hasElementMarker( marker ) ) << "marker " << marker << "does not found in mesh";
        M_definePressureCstMeshRanges.push_back( markedelements(this->mesh(),markers) );
    }
    if ( M_definePressureCstMeshRanges.empty() )
        M_definePressureCstMeshRanges.push_back( M_rangeMeshElements );
    else
        M_definePressureCstOnlyOneZoneAppliedOnWholeMesh = false;



    if ( this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        M_XhMeanPressureLM.resize( M_definePressureCstMeshRanges.size() );
        if ( M_definePressureCstOnlyOneZoneAppliedOnWholeMesh )
            M_XhMeanPressureLM[0] = space_meanpressurelm_type::New( _mesh=this->mesh(), _worldscomm=this->localNonCompositeWorldsComm() );
        else
        {
            for ( int k=0;k<M_definePressureCstMeshRanges.size();++k )
                M_XhMeanPressureLM[k] = space_meanpressurelm_type::New( _mesh=this->mesh(), _worldscomm=this->localNonCompositeWorldsComm(),
                                                                        _range=M_definePressureCstMeshRanges[k] );
        }
    }
    else if ( this->definePressureCstMethod() == "algebraic" )
    {
        auto p = this->functionSpacePressure()->element();

        M_definePressureCstAlgebraicOperatorMeanPressure.clear();
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

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateVelocityExtrapolated()
{
    if ( this->isStationary() )
        return;

    CHECK( M_fieldVelocityExtrapolated ) << "should be init";
    M_fieldVelocityExtrapolated->zero();

    if ( this->currentTime() == this->timeInitial() )
    {
        if ( M_bdfVelocity->iteration() == 0 )
            M_fieldVelocityExtrapolated->add( 1, *M_fieldVelocity );
        else if ( M_bdfVelocity->iteration() == 1 )
            M_fieldVelocityExtrapolated->add( 1, M_bdfVelocity->unknown(1) );
        else if ( M_bdfVelocity->iteration() > 1 )
        {
            M_fieldVelocityExtrapolated->add( 2, M_bdfVelocity->unknown(1) );
            M_fieldVelocityExtrapolated->add( -1, M_bdfVelocity->unknown(2) );
            // TODO : BDF with order > 2
        }
    }
    else
    {
        if ( M_timeStepping == "BDF" )
        {
            *M_fieldVelocityExtrapolated = *M_bdfVelocity->polyPtr();
        }
        else if ( M_timeStepping == "Theta" )
        {
            if ( M_bdfVelocity->iteration() == 1 )
                M_fieldVelocityExtrapolated->add( 1, M_bdfVelocity->unknown(0) );
            else if ( M_bdfVelocity->iteration() > 1 )
            {
                M_fieldVelocityExtrapolated->add( 2 /*3./2.*/, M_bdfVelocity->unknown(0) );
                M_fieldVelocityExtrapolated->add( -1 /*-1./2.*/, M_bdfVelocity->unknown(1) );
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::startTimeStepPreProcess()
{
    this->log("FluidMechanics","startTimeStepPreProcess", "start" );

    this->updateVelocityExtrapolated();

    if ( M_usePreviousSolution )
        *M_vectorPreviousSolution = *this->algebraicBlockVectorSolution()->vectorMonolithic();

    this->log("FluidMechanics","startTimeStepPreProcess", "finish" );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::startTimeStep( bool applyPreProcess )
{
    this->log("FluidMechanics","startTimeStep", "start" );

    if ( applyPreProcess )
        this->startTimeStepPreProcess();

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    // start time step
    if ( !this->doRestart() )
    {
        M_bdfVelocity->start( *M_fieldVelocity );
        M_savetsPressure->start( *M_fieldPressure );
        M_bodySetBC.startTimeStep();
        if ( this->hasTurbulenceModel() )
            M_turbulenceModelType->startTimeStep();

        // windkessel outlet
        for ( auto const& [bcName,bcData] : M_boundaryConditions->outletWindkessel() )
        {
            auto & windkesselData = M_fluidOutletWindkesselData.at(bcName);
            CHECK( std::get<0>( windkesselData ).size() == 2 && std::get<1>( windkesselData ).size() == 2 ) << "something wrong";
            for (int k=0;k<2;++k)
                std::get<1>( windkesselData )[k]->start( unwrap_ptr(std::get<0>( windkesselData )[k]) );
        }

    }
    // up current time
    this->updateTime( M_bdfVelocity->time() );

    if ( true )
    {
        *M_vectorPreviousVelocityExtrapolated = *M_vectorVelocityExtrapolated;
        this->updateVelocityExtrapolated();
    }

    if ( M_usePreviousSolution )
        *M_vectorPreviousSolution = *this->algebraicBlockVectorSolution()->vectorMonolithic();

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
    for ( auto const& [bcName,bcData] : M_boundaryConditions->outletWindkessel() )
    {
        auto & windkesselData = M_fluidOutletWindkesselData.at(bcName);
        CHECK( std::get<0>( windkesselData ).size() == 2 && std::get<1>( windkesselData ).size() == 2 ) << "something wrong";
        for (int k=0;k<2;++k)
            std::get<1>( windkesselData )[k]->next( unwrap_ptr(std::get<0>( windkesselData )[k]) );
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    if (this->hasMeshMotion())
        this->meshMotionTool()->updateTimeStep();
#endif

    M_bodySetBC.updateTimeStep();

    if ( this->hasTurbulenceModel() )
        M_turbulenceModelType->updateTimeStep();

    bool rebuildCstAssembly = false;
    if ( M_timeStepping == "BDF" )
    {
        int previousTimeOrder = this->timeStepBDF()->timeOrder();
        M_bdfVelocity->next( *M_fieldVelocity );
        M_savetsPressure->next( *M_fieldPressure );
        int currentTimeOrder = this->timeStepBDF()->timeOrder();
        rebuildCstAssembly = previousTimeOrder != currentTimeOrder && this->timeStepBase()->strategy() == TS_STRATEGY_DT_CONSTANT;
        this->updateTime( M_bdfVelocity->time() );
    }
    else if ( M_timeStepping == "Theta" )
    {
        M_bdfVelocity->next( *M_fieldVelocity );
        M_savetsPressure->next( *M_fieldPressure );
        this->updateTime( M_bdfVelocity->time() );
    }

    if ( rebuildCstAssembly )
        this->setNeedToRebuildCstPart(true);

    if ( true )
    {
        *M_vectorPreviousVelocityExtrapolated = *M_vectorVelocityExtrapolated;
        this->updateVelocityExtrapolated();
    }

    if ( M_usePreviousSolution )
        *M_vectorPreviousSolution = *this->algebraicBlockVectorSolution()->vectorMonolithic();

    // update all expressions in bc or in house prec
    this->updateParameterValues();

    using mesh_adaptation_type = typename super_type::super_model_meshes_type::mesh_adaptation_type;
    this->template updateMeshAdaptation<mesh_type>( this->keyword(),
                                                    mesh_adaptation_type::createEvent<mesh_adaptation_type::Event::Type::each_time_step>( this->time(),M_bdfVelocity->iteration() ),
                                                    this->symbolsExpr() );

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
    auto algebraicFactory = this->algebraicFactory();
    if ( !algebraicFactory )
        return;

    if ( M_timeStepping == "Theta" )
    {
        M_timeStepThetaSchemePreviousContrib->zero();
        this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
        ModelAlgebraic::DataUpdateResidual dataResidual( this->algebraicBlockVectorSolution()->vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, true, false );
        dataResidual.addInfo( prefixvm( this->prefix(), "time-stepping.evaluate-residual-without-time-derivative" ) );

        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
        algebraicFactory->evaluateResidual( dataResidual );
        algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );

        if ( M_stabilizationGLS )
        {
            auto & dataInfos = algebraicFactory->dataInfos();
            //*dataInfos.vectorInfo( "time-stepping.previous-solution" ) = *M_blockVectorSolution.vectorMonolithic();
            dataInfos.addParameterValuesInfo( "time-stepping.previous-parameter-values", M_currentParameterValues );
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

        auto mphysics = this->materialsProperties()->physicsFromMaterial( matName, this->physicsFromCurrentType() );
        CHECK( mphysics.size() <= 1 ) << "something wrong";
        if ( mphysics.empty() )
            continue;
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(mphysics.begin()->second);
        auto const& matProps = this->materialsProperties()->materialProperties( matName );

        auto const sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u),idv(p),*physicFluidData,matProps,true);
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

        auto mphysics = this->materialsProperties()->physicsFromMaterial( matName, this->physicsFromCurrentType() );
        CHECK( mphysics.size() <= 1 ) << "something wrong";
        if ( mphysics.empty() )
            continue;
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(mphysics.begin()->second);
        auto const& matProps = this->materialsProperties()->materialProperties( matName );

        auto const sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u),idv(p),*physicFluidData,matProps,true);
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
    auto mmt = this->meshMotionTool();
    bool meshIsOnRefAtBegin = mmt->isOnReferenceMesh();
    if ( !meshIsOnRefAtBegin )
        mmt->revertReferenceMesh( false );

    // current solution
    auto const& u = this->fieldVelocity();
    auto const& p = this->fieldPressure();
    // identity Matrix
    auto const Id = eye<nDim,nDim>();
    // deformation tensor
    auto Fa = Id+gradv(*mmt->displacement());

    //auto itFindRange = M_rangeMeshFacesByMaterial.find( "moving-boundary" );
    //CHECK( itFindRange != M_rangeMeshFacesByMaterial.end() ) << "not find range moving-boundary";
    CHECK( M_rangeDistributionByMaterialName ) << "M_rangeDistributionByMaterialName is not init";
    for ( auto const& rangeFacesMat : M_rangeDistributionByMaterialName->rangeMeshFacesByMaterial( nameOfRange ) )
    {
        std::string matName = rangeFacesMat.first;
        auto const& rangeFaces = rangeFacesMat.second;

        auto mphysics = this->materialsProperties()->physicsFromMaterial( matName, this->physicsFromCurrentType() );
        CHECK( mphysics.size() <= 1 ) << "something wrong";
        if ( mphysics.empty() )
            continue;
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(mphysics.begin()->second);
        auto const& matProps = this->materialsProperties()->materialProperties( matName );

        // stress tensor : -p*Id + 2*mu*D(u)
        auto const sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u)*inv(Fa),idv(p),*physicFluidData,matProps,true);
        fieldToUpdate->on(_range=rangeFaces,
                          _expr=sigmav*det(Fa)*trans(inv(Fa))*N(),
                          _geomap=this->geomap() );

#if 0
        //
        auto resRef = integrate(_range=rangeFaces,_expr=idv(fieldToUpdate)).evaluate();
        mmt->revertMovingMesh( false );
        // current solution
        auto const& u = this->fieldVelocity();
        auto const& p = this->fieldPressure();
        // identity Matrix
        auto const Id = eye<nDim,nDim>();
        // deformation tensor
        auto Fa = Id+gradv(*mmt->displacement());
        auto InvFa = det(Fa)*inv(Fa);
        auto const Sigmav = Feel::FeelModels::fluidMecStressTensor(gradv(u),idv(p),*this->materialProperties(),matName,true);
        auto resMove = integrate(_range=rangeFaces,_expr=Sigmav*N() ).evaluate();
        mmt->revertReferenceMesh( false );
        if ( this->worldComm().isMasterRank() )
        {
            std::cout << "resRef " << resRef << "\n";
            std::cout << "resMove " << resMove << "\n";
            std::cout << "diff " << (resRef-resMove).norm() << "\n";
        }
#endif
    }

    if ( !meshIsOnRefAtBegin )
        mmt->revertMovingMesh( false );
}

//---------------------------------------------------------------------------------------------------------//



//#if defined( FEELPP_MODELS_HAS_MESHALE )

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateALEmesh()
{
    this->updateALEmesh( this->symbolsExpr() );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateALEmeshImpl()
{
    // VINCENT TODO move this code
#if 0
    //-------------------------------------------------------------------//
    // compute ALE map
    //std::vector< mesh_ale_type::ale_map_element_type> polyBoundarySet = { *M_meshDisplacementOnInterface };
    //M_meshALE->update(*M_meshDisplacementOnInterface/*polyBoundarySet*/);
    M_meshALE->updateMovingMesh();
#endif
    //-------------------------------------------------------------------//

    if ( this->doCIPStabConvection() )
    {
        M_fieldMeshVelocityUsedWithStabCIP->on(_range=M_rangeMeshElements,_expr=idv(this->meshMotionTool()->velocity()) );
        sync( *M_fieldMeshVelocityUsedWithStabCIP, "=" );
    }
    //-------------------------------------------------------------------//
    // update operator PCD
    if ( M_preconditionerAttachPCD && this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasOperatorPCD("pcd") )
    {
        this->log("FluidMechanics","updateALEmesh", "rebuild operatorPCD");
        CHECK( this->algebraicFactory()->preconditionerTool()->hasOperatorPCD("pcd") ) << "operator PCD does not init";
        typedef Feel::Alternatives::OperatorPCD<space_velocity_type,space_pressure_type> op_pcd_type;
        std::shared_ptr<op_pcd_type> myOpPCD =
            std::dynamic_pointer_cast<op_pcd_type>( this->algebraicFactory()->preconditionerTool()->operatorPCD( "pcd" ) );
        myOpPCD->assemble();
    }
    //-------------------------------------------------------------------//
    // update operator PMM
    if ( M_preconditionerAttachPMM && this->algebraicFactory() && this->algebraicFactory()->preconditionerTool()->hasOperatorPMM("pmm") )
    {
        this->log("FluidMechanics","updateALEmesh", "rebuild operatorPMM");
        CHECK( this->algebraicFactory()->preconditionerTool()->hasOperatorPMM("pmm") ) << "operator PMM does not init";
        typedef Feel::Alternatives::OperatorPMM<space_pressure_type> op_pmm_type;
        std::shared_ptr<op_pmm_type> myOpPMM =
            std::dynamic_pointer_cast<op_pmm_type>( this->algebraicFactory()->preconditionerTool()->operatorPMM( "pmm" ) );
        myOpPMM->assemble();
    }
}

//#endif

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


//---------------------------------------------------------------------------------------------------------//

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeFlowRate( std::string const& marker, bool useExteriorNormal ) const
{
    return this->computeFlowRate( std::set<std::string>( { marker } ),useExteriorNormal );
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::computeFlowRate( std::set<std::string> const& markers, bool useExteriorNormal ) const
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
         (this->hasMeshMotion() && ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-neumann" || this->couplingFSIcondition() == "nitsche" ) )
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
    int nBlock = 2;
    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
        nBlock+=M_XhMeanPressureLM.size();

    if ( M_boundaryConditions->hasVelocityImposedLagrangeMultiplier() )
        ++nBlock;

    if ( !M_boundaryConditions->pressureImposed().empty() )
    {
        ++nBlock;
        if ( nDim == 3 )
            ++nBlock;
    }

    for ( auto const& [bcName,bcData] : M_boundaryConditions->outletWindkessel() )
        if ( bcData->useImplicitCoupling() )
            nBlock += 2;

    for ( auto const& [bpname,bbc] : M_bodySetBC )
    {
        ++nBlock; // translational velocity
        if ( !bbc.isInNBodyArticulated() )
            ++nBlock; // angular velocity
    }
    for ( auto const& nba : M_bodySetBC.nbodyArticulated() )
    {
        ++nBlock; // angular velocity
        if ( nba.articulationMethod() != "lm" )
            continue;
        nBlock += nba.articulations().size();
    }
    return nBlock;
}

FLUIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("FluidMechanics","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    auto XhV = this->functionSpaceVelocity();
    auto XhP = this->functionSpacePressure();

    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;
    // first field (velocity+pressure)
    auto ru = stencilRange<0,0>(marked3faces(this->mesh(),1));
    auto rp = stencilRange<1,1>(marked3faces(this->mesh(),1));
#if 0
    myblockGraph(indexBlock,indexBlock) =stencil(_test=this->functionSpace(),_trial=this->functionSpace(),
                                                 _pattern_block=this->blockPattern(),
                                                 //_diag_is_nonzero=false,_close=false,
                                                 _diag_is_nonzero=(nBlock==1),_close=(nBlock==1),
                                                 _range_extended=stencilRangeMap(ru,rp) )->graph();
#endif
    auto blockPat = blockPattern();
    myblockGraph(indexBlock,indexBlock) = stencil(_test=XhV,_trial=XhV,
                                                  _pattern=blockPat(0,0),
                                                  _diag_is_nonzero=false,_close=false,
                                                  _range_extended=stencilRangeMap(ru) )->graph();
    myblockGraph(indexBlock+1,indexBlock+1) = stencil(_test=XhP,_trial=XhP,
                                                      _pattern=blockPat(1,1),
                                                      _diag_is_nonzero=false,_close=false,
                                                      _range_extended=stencilRangeMap(ru) )->graph();
    myblockGraph(indexBlock,indexBlock+1) = stencil(_test=XhV,_trial=XhP,
                                                    _pattern=blockPat(0,1),
                                                    _diag_is_nonzero=false,_close=false )->graph();
    myblockGraph(indexBlock+1,indexBlock) =  myblockGraph(indexBlock,indexBlock+1)->transpose( false );

    indexBlock+=2;

    if ( this->definePressureCst() && this->definePressureCstMethod() == "lagrange-multiplier" )
    {
        for ( int k=0;k<M_XhMeanPressureLM.size();++k )
        {
            myblockGraph(indexBlock,1) = stencil(_test=M_XhMeanPressureLM[k],_trial=XhP,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(1,indexBlock) = stencil(_test=XhP,_trial=M_XhMeanPressureLM[k],
                                                 _diag_is_nonzero=false,_close=false)->graph();
            ++indexBlock;
        }
    }

    if ( M_boundaryConditions->hasVelocityImposedLagrangeMultiplier() )
    {
        myblockGraph(indexBlock,0) = stencil(_test=this->XhDirichletLM(),_trial=XhV,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=XhV,_trial=this->XhDirichletLM(),
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    if ( !M_boundaryConditions->pressureImposed().empty() )
    {
        myblockGraph(indexBlock,0) = stencil(_test=M_spaceLagrangeMultiplierPressureBC,_trial=XhV,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=XhV,_trial=M_spaceLagrangeMultiplierPressureBC,
                                             _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
        if ( nDim == 3 )
        {
            myblockGraph(indexBlock,0) = myblockGraph(indexBlock-1,0);
            myblockGraph(0,indexBlock) = myblockGraph(0,indexBlock-1);
            ++indexBlock;
        }
    }
    for ( auto const& [bcName,bcData] : M_boundaryConditions->outletWindkessel() )
    {
        if ( !bcData->useImplicitCoupling() )
            continue;

        auto range = markedelements(M_fluidOutletWindkesselMesh,bcData->markers());
        auto rangeFirstCol1 = stencilRange<0,0>( range );
         // column with velocity
        myblockGraph(indexBlock,0) = stencil(_test=M_fluidOutletWindkesselSpace,_trial=XhV,
                                             //_pattern_block=patCouplingFirstCol,
                                             _range=stencilRangeMap(rangeFirstCol1),
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(indexBlock+1,0) = myblockGraph(indexBlock,0) ;

        // row with velocity
        myblockGraph(0,indexBlock) = myblockGraph(indexBlock,0)->transpose( false );
        myblockGraph(0,indexBlock+1) = myblockGraph(0,indexBlock);

        // on diag
        myblockGraph(indexBlock,indexBlock) = stencil(_test=M_fluidOutletWindkesselSpace,_trial=M_fluidOutletWindkesselSpace,
                                                      //_pattern=(size_type)Pattern::COUPLED,
                                                      //_pattern_block=patCouplingDiag,
                                                      _diag_is_nonzero=false,_close=false)->graph();

        myblockGraph(indexBlock+1,indexBlock) = myblockGraph(indexBlock,indexBlock);
        myblockGraph(indexBlock+1,indexBlock+1) = myblockGraph(indexBlock,indexBlock);
        indexBlock += 2;
    }

    for ( auto const& [bpname,bbc] : M_bodySetBC )
    {
        size_type startBlockIndexTranslationalVelocity = this->startSubBlockSpaceIndex("body-bc."+bbc.name()+".translational-velocity");
        myblockGraph(startBlockIndexTranslationalVelocity,startBlockIndexTranslationalVelocity) = stencil(_test=bbc.spaceTranslationalVelocity(),_trial=bbc.spaceTranslationalVelocity(),
                                                                                                          _diag_is_nonzero=false,_close=false)->graph();
        if ( !bbc.isInNBodyArticulated() )
        {
            size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+bbc.name()+".angular-velocity");
            myblockGraph(startBlockIndexAngularVelocity,startBlockIndexAngularVelocity) = stencil(_test=bbc.spaceAngularVelocity(),_trial=bbc.spaceAngularVelocity(),
                                                                                                  _diag_is_nonzero=false,_close=false)->graph();
        }
        //indexBlock +=2;
    }

    for ( auto const& nba : M_bodySetBC.nbodyArticulated() )
    {
        size_type startBlockIndexAngularVelocity = this->startSubBlockSpaceIndex("body-bc."+nba.name()+".angular-velocity");
        myblockGraph(startBlockIndexAngularVelocity,startBlockIndexAngularVelocity) = stencil(_test=nba.spaceAngularVelocity(),_trial=nba.spaceAngularVelocity(),
                                                                                              _diag_is_nonzero=false,_close=false)->graph();
        if ( nba.articulationMethod() != "lm" )
            continue;
        for ( auto const& ba : nba.articulations() )
        {
            auto const& bbc1 = ba.body1();
            auto const& bbc2 = ba.body2();
            size_type startBlockIndexTranslationalVelocityBody1 = this->startSubBlockSpaceIndex("body-bc."+bbc1.name()+".translational-velocity");
            size_type startBlockIndexTranslationalVelocityBody2 = this->startSubBlockSpaceIndex("body-bc."+bbc2.name()+".translational-velocity");
            size_type startBlockIndexArticulationLMTranslationalVelocity = this->startSubBlockSpaceIndex( "body-bc.articulation-lm."+ba.name()+".translational-velocity");

            auto createArticulationGraph = [] ( auto const& map1, auto const& map2, bool closeGraph = true, size_type pattern = Pattern::COUPLED ) {
                Feel::Context graph_prop( pattern );
                auto sparsity_graph = std::make_shared<GraphCSR>( map1, map2 );
                if ( map2->nLocalDofWithGhost() > 0 )
                {
                    for ( size_type i=0;i<map1->nLocalDofWithGhost();++i )
                    {
                        const size_type ig1 = map1->mapGlobalProcessToGlobalCluster(i);
                        const rank_type theproc = map1->procOnGlobalCluster( ig1 );
                        // define row in graph
                        auto& row = sparsity_graph->row( ig1 );
                        row.template get<0>() = theproc;
                        const size_type il1 = ig1 - map1->firstDofGlobalCluster( theproc );
                        row.template get<1>() = il1;
                        //up the pattern graph
                        if ( graph_prop.test( Pattern::COUPLED ) )
                        {
                            for ( size_type j=0;j<map2->nLocalDofWithGhost();++j )
                                row.template get<2>().insert( map2->mapGlobalProcessToGlobalCluster(j) );
                        }
                        else if ( graph_prop.test( Pattern::DEFAULT ) )
                        {
                            row.template get<2>().insert( map2->mapGlobalProcessToGlobalCluster(i) );
                        }
                    }
                }
                if ( closeGraph )
                    sparsity_graph->close();
                return sparsity_graph;
            };

            myblockGraph(startBlockIndexTranslationalVelocityBody1,startBlockIndexArticulationLMTranslationalVelocity) = createArticulationGraph( bbc1.spaceTranslationalVelocity()->mapPtr(),
                                                                                                                                                  ba.dataMapLagrangeMultiplierTranslationalVelocity(),
                                                                                                                                                  false, Pattern::DEFAULT );
            myblockGraph(startBlockIndexTranslationalVelocityBody2,startBlockIndexArticulationLMTranslationalVelocity) = createArticulationGraph( bbc2.spaceTranslationalVelocity()->mapPtr(),
                                                                                                                                                  ba.dataMapLagrangeMultiplierTranslationalVelocity(),
                                                                                                                                                  false, Pattern::DEFAULT );
            myblockGraph(startBlockIndexArticulationLMTranslationalVelocity,startBlockIndexTranslationalVelocityBody1) = createArticulationGraph( ba.dataMapLagrangeMultiplierTranslationalVelocity(),
                                                                                                                                                  bbc1.spaceTranslationalVelocity()->mapPtr(),
                                                                                                                                                  false, Pattern::DEFAULT );
            myblockGraph(startBlockIndexArticulationLMTranslationalVelocity,startBlockIndexTranslationalVelocityBody2) = createArticulationGraph( ba.dataMapLagrangeMultiplierTranslationalVelocity(),
                                                                                                                                                  bbc2.spaceTranslationalVelocity()->mapPtr(),
                                                                                                                                                  false, Pattern::DEFAULT );
            //++indexBlock;

        }
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
void
FLUIDMECHANICS_CLASS_TEMPLATE_TYPE::updateRangeDistributionByMaterialName( std::string const& key, range_faces_type const& rangeFaces )
{
    if ( !M_rangeDistributionByMaterialName )
    {
        M_rangeDistributionByMaterialName = std::make_shared<RangeDistributionByMaterialName<mesh_type>>();
        auto mom = this->materialsProperties()->materialsOnMesh( this->mesh() );
        M_rangeDistributionByMaterialName->init( mom->rangeMeshElementsByMaterial() );
    }
    M_rangeDistributionByMaterialName->update( key, rangeFaces );
}

} // namespace FeelModels
} // namespace Feel


