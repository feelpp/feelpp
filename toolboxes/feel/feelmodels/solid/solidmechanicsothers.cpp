/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */


#include <feel/feelmodels/solid/solidmechanics.hpp>
#include <feel/feelmodels/modelvf/solidmecfirstpiolakirchhoff.hpp>
#include <feel/feelmodels/modelvf/solidmecincompressibility.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresnormevaluation.hpp>
#include <feel/feelmodels/modelcore/modelmeasuresstatisticsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
std::shared_ptr<std::ostringstream>
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::getInfo() const
{

    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
#if 0
    this->log("SolidMechanics","getInfo", "start" );
    std::string StateTemporal = (this->isStationary())? "Stationary" : "Transient";
    size_type nElt,nDof;
    if (this->hasSolidEquationStandard()) {/*nElt=M_mesh->numGlobalElements();*/ nDof=M_XhDisplacement->nDof();}
#if 0
    else {nElt=M_mesh_1dReduced->numGlobalElements(); nDof=M_Xh_1dReduced->nDof();}
#endif
    std::string ResartMode;
    if (this->doRestart()) ResartMode = "Yes";
    else ResartMode = "No";

    std::string hovisuMode,myexporterType;
    int myexporterFreq=1;
    std::string doExport_str;
    if ( this->hasSolidEquationStandard() )
    {
        if ( M_isHOVisu && M_exporter_ho )
        {
            std::string hovisuSpaceUsed = soption(_name="hovisu.space-used",_prefix=this->prefix());
            if ( hovisuSpaceUsed == "displacement" || hovisuSpaceUsed == "pressure" )
                hovisuMode = "ON (with OperatorLagrangeP1 on "+hovisuSpaceUsed+")";
            else if ( hovisuSpaceUsed == "p1" )
                hovisuMode = "ON (with createP1mesh)";
#if 1 //defined(FEELPP_HAS_VTK)
            myexporterType = M_exporter_ho->type();
            myexporterFreq = M_exporter_ho->freq();
#endif
        }
        else if ( M_exporter )
        {
            hovisuMode = "OFF";
            myexporterType = M_exporter->type();
            myexporterFreq = M_exporter->freq();
        }
        for ( std::string const& fieldName : this->postProcessExportsFields() )
            doExport_str=(doExport_str.empty())? fieldName : doExport_str + " - " + fieldName;
    }
    else
    {
#if 0
        if ( M_exporter_1dReduced )
        {
            hovisuMode = "OFF";
            myexporterType = M_exporter_1dReduced->type();
            myexporterFreq = M_exporter_1dReduced->freq();
        }
#endif
    }

    *_ostr << "\n||==============================================||"
           << "\n||----------Info : SolidMechanics---------------||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository()
           << "\n   Physical Model"
        //<< "\n     -- model : " << M_modelName
        //           << "\n     -- material law : " << this->mechanicalProperties()->materialLaw()
           << "\n     -- has displacement-pressure formulation : " << std::boolalpha << this->hasDisplacementPressureFormulation()
           << "\n     -- time mode : " << StateTemporal;
    *_ostr << this->materialsProperties()->getInfoMaterialParameters()->str();
    if ( this->hasSolidEquationStandard() )
    {
        *_ostr << "\n   Boundary conditions"
               << M_bcDirichletMarkerManagement.getInfoDirichletBC()
               << M_bcNeumannMarkerManagement.getInfoNeumannBC()
               << M_bcNeumannEulerianFrameMarkerManagement.getInfoNeumannEulerianFrameBC()
               << M_bcRobinMarkerManagement.getInfoRobinBC()
               << M_bcFSIMarkerManagement.getInfoFluidStructureInterfaceBC();
#if 0
        *_ostr << "\n   Space Discretization";
        if ( this->hasGeoFile() )
            *_ostr << "\n     -- geo file name   : " << this->geoFile();
        *_ostr << "\n     -- mesh file name   : " << this->meshFile()
               << "\n     -- nb elt in mesh : " << nElt
               << "\n     -- nb dof (displacement) : " << nDof
               << "\n     -- polynomial order : " << nOrder;
        if ( this->hasDisplacementPressureFormulation() )
            *_ostr << "\n     -- nb dof (pressure) : " << M_XhPressure->nDof();
#endif
    }
    if ( !this->isStationary() )
    {
        *_ostr << "\n   Time Discretization"
               << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
               << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
               << "\n     -- time step    : " << this->timeStepBase()->timeStep()
               << "\n     -- type : " << M_timeStepping;
        if ( M_timeStepping == "Newmark" )
        {
            if (this->hasSolidEquationStandard())
                *_ostr << " ( gamma="<< this->timeStepNewmark()->gamma() <<", beta="<< this->timeStepNewmark()->beta() << " )";
#if 0
            else// if (this->is1dReducedModel())
                *_ostr << " ( gamma="<< this->timeStepNewmark1dReduced()->gamma() <<", beta="<< this->timeStepNewmark1dReduced()->beta() << " )";
#endif
        }
        else if ( M_timeStepping == "BDF" )
        {
            if (this->hasSolidEquationStandard())
                *_ostr << " ( order=" << M_timeStepBdfDisplacement->timeOrder() << " )";
        }
        else if ( M_timeStepping == "Theta" )
        {
            *_ostr << " ( theta=" << M_timeStepThetaValue << " )";
        }
        *_ostr << "\n     -- restart mode : " << ResartMode
               << "\n     -- save on disk : " << std::boolalpha << this->timeStepBase()->saveInFile();
        if ( this->timeStepBase()->saveFreq() )
            *_ostr << "\n     -- freq save : " << this->timeStepBase()->saveFreq()
                   << "\n     -- file format save : " << this->timeStepBase()->fileFormat();
    }
    if ( !myexporterType.empty() )
        *_ostr << "\n   Exporter"
               << "\n     -- type            : " << myexporterType
               << "\n     -- high order visu : " << hovisuMode
               << "\n     -- freq save       : " << myexporterFreq
               << "\n     -- fields exported : " << doExport_str;
    *_ostr << "\n   Processors"
           << "\n     -- number of proc environnement : " << Environment::worldComm().globalSize()
           << "\n     -- environement rank : " << Environment::worldComm().rank()
           << "\n     -- global rank : " << this->worldComm().globalRank()
           << "\n     -- local rank : " << this->worldComm().localRank()
           << "\n   Numerical Solver"
           << "\n     -- solver : " << M_solverName;
    if ( this->algebraicFactory() )
        *_ostr << this->algebraicFactory()->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n";

    this->log("SolidMechanics","getInfo", "finish" );
#endif
    return _ostr;
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateInformationObject( nl::json & p ) const
{
    if ( !this->isUpdatedForUse() )
        return;
    if ( p.contains( "Environment" ) )
        return;

    super_type::super_model_base_type::updateInformationObject( p["Environment"] );

    super_physics_type::updateInformationObjectFromCurrentType( p["Physics"] );

    nl::json subPt;
    subPt.emplace( "time mode", std::string( (this->isStationary())?"Stationary":"Transient") );
    p["Physics2"] = subPt;

    // Materials properties
    if ( this->materialsProperties() )
        this->materialsProperties()->updateInformationObject( p["Materials Properties"] );

    if (this->hasSolidEquationStandard())
    {
        super_type::super_model_meshes_type::updateInformationObject( p["Meshes"] );

        // Boundary Conditions
        M_boundaryConditions.updateInformationObject( p["Boundary Conditions"] );

        subPt.clear();
        this->functionSpaceDisplacement()->updateInformationObject( subPt["Displacement"] );
        if ( this->hasDisplacementPressureFormulation() )
            this->functionSpacePressure()->updateInformationObject( subPt["Pressure"] );
        p["Function Spaces"] = subPt;

        this->modelFields().updateInformationObject( p["Fields"] );

        if ( this->algebraicFactory() )
            this->algebraicFactory()->updateInformationObject( p["Algebraic Solver"] );
    }

    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->updateInformationObject( p["Toolbox Solid 1d Reduced"] );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
tabulate_informations_ptr_t
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const
{
    auto tabInfo = TabulateInformationsSections::New();

    // Environment
    if ( jsonInfo.contains("Environment") )
        tabInfo->add( "Environment",  super_type::super_model_base_type::tabulateInformations( jsonInfo.at("Environment"), tabInfoProp ) );

    // Physics
    if ( jsonInfo.contains("Physics") )
        tabInfo->add( "Physics", super_physics_type::tabulateInformations( jsonInfo.at("Physics"), tabInfoProp ) );

    if ( jsonInfo.contains("Physics2") )
    {
        Feel::Table tabInfoPhysics;
        TabulateInformationTools::FromJSON::addAllKeyToValues( tabInfoPhysics, jsonInfo.at("Physics2"), tabInfoProp );
        tabInfo->add( "Physics2", TabulateInformations::New( tabInfoPhysics ) );
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
    if ( jsonInfo.contains("Function Spaces") )
    {
        auto const& jsonInfoFunctionSpaces = jsonInfo.at("Function Spaces");
        auto tabInfoFunctionSpaces = TabulateInformationsSections::New();
        for ( std::string const& spaceName : std::vector<std::string>({"Displacement","Pressure"}) )
        {
            if ( jsonInfoFunctionSpaces.contains( spaceName ) )
                tabInfoFunctionSpaces->add( spaceName, TabulateInformationTools::FromJSON::tabulateInformationsFunctionSpace( jsonInfoFunctionSpaces.at( spaceName ), tabInfoProp ) );
        }
        tabInfo->add( "Function Spaces", tabInfoFunctionSpaces );
    }

    // fields
    if ( jsonInfo.contains("Fields") )
        tabInfo->add( "Fields", TabulateInformationTools::FromJSON::tabulateInformationsModelFields( jsonInfo.at("Fields"), tabInfoProp ) );

    // Algebraic Solver
    if ( jsonInfo.contains( "Algebraic Solver" ) )
        tabInfo->add( "Algebraic Solver", model_algebraic_factory_type::tabulateInformations( jsonInfo.at("Algebraic Solver"), tabInfoProp ) );

    // Subtoolbox
    if ( this->hasSolidEquation1dReduced() && jsonInfo.contains( "Toolbox Solid 1d Reduced" ) )
        tabInfo->add( "Toolbox Solid 1d Reduced", M_solid1dReduced->tabulateInformations( jsonInfo.at("Toolbox Solid 1d Reduced"), tabInfoProp ) );

    return tabInfo;
}

//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::setModelName( std::string const& type )
{
    // if pde change -> force to rebuild all algebraic data at next solve
    if ( type != M_modelName )
        this->setNeedToRebuildCstPart(true);

    if ( type == "Elasticity" )
        M_modelName = "Elasticity";
    else if ( type == "Elasticity-Large-Deformation" )
        M_modelName="Elasticity-Large-Deformation";
    else if ( type == "Hyper-Elasticity" )
        M_modelName = "Hyper-Elasticity";
    else if ( type == "Generalised-String" )
        M_modelName = "Generalised-String";
    else
        CHECK( false ) << "invalid modelName "<< type << "\n";

    M_isStandardModel= (M_modelName != "Generalised-String" );
    M_is1dReducedModel=!M_isStandardModel;
}
//---------------------------------------------------------------------------------------------------//
#endif

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::setSolver( std::string const& type )
{
    // if solver change -> force to rebuild all algebraic data at next solve
    if ( type != M_solverName )
        this->setNeedToRebuildCstPart(true);
    if ( type == "automatic" )
        M_solverName = "automatic";
    else if ( type == "LinearSystem" )
        M_solverName="LinearSystem";
    else if ( type == "Newton" )
        M_solverName="Newton";
    else
        CHECK( false ) << "invalid solver name " << type << "\n";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
int
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 0;
    if ( this->hasSolidEquationStandard() )
    {
        ++nBlock;
        if ( this->hasDisplacementPressureFormulation() )
            ++nBlock;
        if ( M_timeSteppingUseMixedFormulation )
            ++nBlock;
    }
    if ( this->hasSolidEquation1dReduced() )
    {
        nBlock += M_solid1dReduced->algebraicBlockVectorSolution()->size();
    }
    return nBlock;
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("SolidMechanics","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;

    if ( this->hasSolidEquationStandard() )
    {
        myblockGraph(indexBlock,indexBlock) =stencil(_test=this->functionSpace(),_trial=this->functionSpace(),
                                                     //_pattern_block=this->blockPattern(),
                                                     _diag_is_nonzero=(nBlock==1),
                                                     _close=(nBlock==1) )->graph();
        ++indexBlock;

        if ( this->hasDisplacementPressureFormulation() )
        {
            myblockGraph(indexBlock,0) = stencil(_test=M_XhPressure,_trial=M_XhDisplacement,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(0,indexBlock) = stencil(_test=M_XhDisplacement,_trial=M_XhPressure,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(indexBlock,indexBlock) = stencil(_test=M_XhPressure,_trial=M_XhPressure,
                                                          _diag_is_nonzero=false,_close=false)->graph();
            ++indexBlock;
        }

        if ( M_timeSteppingUseMixedFormulation )
        {
            myblockGraph(indexBlock,0) = stencil(_test=M_XhDisplacement,_trial=M_XhDisplacement,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(0,indexBlock) = stencil(_test=M_XhDisplacement,_trial=M_XhDisplacement,
                                                 _diag_is_nonzero=false,_close=false)->graph();
            myblockGraph(indexBlock,indexBlock) = stencil(_test=M_XhDisplacement,_trial=M_XhDisplacement,
                                                          _diag_is_nonzero=false,_close=false)->graph();
            ++indexBlock;
        }
    }

    if ( this->hasSolidEquation1dReduced() )
    {
        int startIndexBlock1dReduced = indexBlock;
        auto blockMat1dReduced = M_solid1dReduced->buildBlockMatrixGraph();
        for (int tk1=0;tk1<blockMat1dReduced.nRow() ;++tk1 )
            for (int tk2=0;tk2<blockMat1dReduced.nCol() ;++tk2 )
                myblockGraph(startIndexBlock1dReduced+tk1,startIndexBlock1dReduced+tk2) = blockMat1dReduced(tk1,tk2);
    }

    myblockGraph.close();

    this->log("SolidMechanics","buildBlockMatrixGraph", "finish" );
    return myblockGraph;
}


//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    auto mfields = this->modelFields();
    auto se = this->symbolsExpr( mfields );
    if ( this->hasSolidEquationStandard() )
        this->exportResults( time, mfields, se, this->exprPostProcessExports( se ) );
    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->exportResults( time, se );

    // this->log("SolidMechanics","exportResults",(boost::format("start at time %1%")%time).str() );
    // this->timerTool("PostProcessing").start();

    // if (!M_isHOVisu)
    //     this->exportFields( time );
    // else
    //     this->exportFieldsImplHO( time );

    // this->exportMeasures( time );

    // this->timerTool("PostProcessing").stop("exportResults");
    // if ( this->scalabilitySave() )
    // {
    //     if ( !this->isStationary() )
    //         this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
    //     this->timerTool("PostProcessing").save();
    // }

    // this->log("SolidMechanics","exportResults", "finish" );

} // SolidMechanics::export

#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::exportFields( double time )
{
    if (this->isStandardModel())
    {
        bool hasFieldToExport = this->updateExportedFields( M_exporter, M_postProcessFieldExported, time );
        if ( hasFieldToExport )
        {
            M_exporter->save();
            this->upload( M_exporter->path() );
        }
    }
    else
    {
        bool hasFieldToExport = this->updateExportedFields1dReduced( M_exporter_1dReduced, M_postProcessFieldExported, time );
        if ( hasFieldToExport )
        {
            M_exporter_1dReduced->save();
            this->upload( M_exporter_1dReduced->path() );
        }
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateExportedFields( exporter_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    bool hasFieldToExport = false;
    if ( fields.find( "pid" ) != fields.end() )
    {
        exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
        hasFieldToExport = true;
    }
    if ( fields.find( "displacement" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"displacement"), this->fieldDisplacement() );
        hasFieldToExport = true;
    }
    if ( fields.find( "pressure" ) != fields.end() && M_useDisplacementPressureFormulation )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"pressure"), *M_fieldPressure );
        hasFieldToExport = true;
    }
    if ( fields.find( "velocity" ) != fields.end() && M_fieldVelocity )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"velocity"), this->fieldVelocity() );
        hasFieldToExport = true;
    }
    if ( fields.find( "acceleration" ) != fields.end() && M_fieldAcceleration )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"acceleration"), this->fieldAcceleration() );
        hasFieldToExport = true;
    }
    if ( fields.find( "normal-stress" ) != fields.end() )
    {
        this->updateNormalStressFromStruct();
        exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"), *M_fieldNormalStressFromStruct );
        hasFieldToExport = true;
    }
    if ( ( fields.find( "Von-Mises" ) != fields.end() ) ||
         ( fields.find( "Tresca" ) != fields.end() ) ||
         ( fields.find( "principal-stresses" ) != fields.end() ) )
    {
        this->updateStressCriterions();
        if ( fields.find( "Von-Mises" ) != fields.end() )
            exporter->step( time )->add( prefixvm(this->prefix(),"von-mises-criterions"), *M_fieldVonMisesCriterions );
        if ( fields.find( "Tresca" ) != fields.end() )
            exporter->step( time )->add( prefixvm(this->prefix(),"tresca-criterions"), *M_fieldTrescaCriterions );
        if ( fields.find( "principal-stresses" ) != fields.end() )
            for (int d=0;d<M_fieldsPrincipalStresses.size();++d)
                exporter->step( time )->add( prefixvm(this->prefix(),(boost::format("princial-stress-%1%")%d).str() ), *M_fieldsPrincipalStresses[d] );
#if 0
            exporter->step( time )->add( prefixvm(this->prefix(),"sigma_xx"), M_fieldStressTensor->comp( Component::X,Component::X ) );
            exporter->step( time )->add( prefixvm(this->prefix(),"sigma_xy"), M_fieldStressTensor->comp( Component::X,Component::Y ) );
            exporter->step( time )->add( prefixvm(this->prefix(),"sigma_yx"), M_fieldStressTensor->comp( Component::Y,Component::X ) );
            exporter->step( time )->add( prefixvm(this->prefix(),"sigma_yy"), M_fieldStressTensor->comp( Component::Y,Component::Y ) );
            if (nDim==3)
            {
                exporter->step( time )->add( prefixvm(this->prefix(),"sigma_xz"), M_fieldStressTensor->comp( Component::X,Component::Z ) );
                exporter->step( time )->add( prefixvm(this->prefix(),"sigma_yz"), M_fieldStressTensor->comp( Component::Y,Component::Z ) );
                exporter->step( time )->add( prefixvm(this->prefix(),"sigma_zx"), M_fieldStressTensor->comp( Component::Z,Component::X ) );
                exporter->step( time )->add( prefixvm(this->prefix(),"sigma_zy"), M_fieldStressTensor->comp( Component::Z,Component::Y ) );
                exporter->step( time )->add( prefixvm(this->prefix(),"sigma_zz"), M_fieldStressTensor->comp( Component::Z,Component::Z ) );
            }
#endif
            hasFieldToExport = true;
        }
    if ( fields.find( "material-properties" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"YoungModulus"), this->mechanicalProperties()->fieldYoungModulus() );
        exporter->step( time )->add( prefixvm(this->prefix(),"PoissonCoefficient"), this->mechanicalProperties()->fieldCoeffPoisson() );
        exporter->step( time )->add( prefixvm(this->prefix(),"Lame1"), this->mechanicalProperties()->fieldCoeffLame1() );
        exporter->step( time )->add( prefixvm(this->prefix(),"Lame2"), this->mechanicalProperties()->fieldCoeffLame2() );
        hasFieldToExport = true;
    }
    /*if ( fields.find( "fsi" ) != fields.end() )
     {
         exporter->step( time )->add( prefixvm(this->prefix(),"velocity-interface-from-fluid"), this->fieldVelocityInterfaceFromFluid() );
         hasFieldToExport = true;
     }*/

    return hasFieldToExport;
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
bool
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateExportedFields1dReduced( exporter_1dreduced_ptrtype exporter, std::set<std::string> const& fields, double time )
{
    if ( !exporter ) return false;
    if ( !exporter->doExport() ) return false;

    bool hasFieldToExport = false;
    if ( fields.find( "pid" ) != fields.end() )
    {
        exporter->step( time )->addRegions( prefixvm(this->prefix(),"1d-reduced"),
                                            prefixvm(this->prefix(),prefixvm(this->subPrefix(),"1d-reduced") ) );
        hasFieldToExport = true;
    }
    if ( fields.find( "displacement" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"1d-reduced-displacement"), *M_disp_1dReduced );
        hasFieldToExport = true;
    }
    if ( fields.find( "velocity" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"1d-reduced-velocity"), this->timeStepNewmark1dReduced()->currentVelocity()  /**M_velocity_1dReduced*/ );
        hasFieldToExport = true;
    }
    if ( fields.find( "acceleration" ) != fields.end() )
    {
        exporter->step( time )->add( prefixvm(this->prefix(),"1d-reduced-acceleration"), this->timeStepNewmark1dReduced()->currentAcceleration() );
        hasFieldToExport = true;
    }
    return hasFieldToExport;
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::exportFieldsImplHO( double time )
{
    std::set<std::string> const& fields = M_postProcessFieldExported;
    if (this->isStandardModel())
    {
#if 1 //defined(FEELPP_HAS_VTK)
        if ( !M_exporter_ho->doExport() ) return;
        //M_exporter_ho->step( time )->setMesh( M_displacementVisuHO->mesh() );

        bool hasFieldToExport = false;
        if ( fields.find( "pid" ) != fields.end() )
        {
            M_exporter_ho->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
            hasFieldToExport = true;
        }
        if ( fields.find( "displacement" ) != fields.end() )
        {
            M_opIdisplacement->apply(this->fieldDisplacement(),*M_displacementVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"displacement-ho"), *M_displacementVisuHO );
            hasFieldToExport = true;
        }
        if ( fields.find( "pressure" ) != fields.end() && M_useDisplacementPressureFormulation )
        {
            M_opIpressure->apply(*M_fieldPressure,*M_pressureVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure-ho"), *M_pressureVisuHO );
            hasFieldToExport = true;
        }
        if ( fields.find( "velocity" ) != fields.end() )
        {
            auto velocityVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentVelocity(),velocityVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity-ho"), velocityVisuHO );
            hasFieldToExport = true;
        }
        if ( fields.find( "acceleration" ) != fields.end() )
        {
            auto accelerationVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentAcceleration(),accelerationVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"acceleration-ho"), accelerationVisuHO );
            hasFieldToExport = true;
        }
        if ( fields.find( "normal-stress" ) != fields.end() )
        {
            auto normalstressVisuHO = M_XhVectorialVisuHO->element();
            this->updateNormalStressFromStruct();
            M_opInormalstress->apply( *M_fieldNormalStressFromStruct,normalstressVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"normalstress-ho"), normalstressVisuHO );
            hasFieldToExport = true;
        }
        if ( hasFieldToExport )
            M_exporter_ho->save();
#endif
    }

}
#endif

#if 0
namespace detail
{
template <typename ElementTensor2Type>
typename ElementTensor2Type::component_functionspace_type::element_ptrtype
componentFieldFromTensor2Field( ElementTensor2Type const uTensor2, uint16_type c1, uint16_type c2 )
{
    auto compSpace = uTensor2.compSpace();
    auto uComp = compSpace->elementPtr();
    static const uint16_type nComponents2 = ElementTensor2Type::functionspace_type::nComponents2;
    auto dofTensor2 = uTensor2.functionSpace()->dof();
    auto dofComp = compSpace->dof();
    int nLocDofPerComp = dofTensor2->nLocalDof( true );
    CHECK(  nLocDofPerComp == dofComp->nLocalDof() ) << "must be same size";

    for ( auto const& elt : elements(uTensor2.mesh()) )
    {
        for ( size_type j =0; j < nLocDofPerComp;++j )
        {
            uint16_type compTensor2 = c1*nComponents2+c2;
            size_type gdofTensor2 = dofTensor2->localToGlobal( elt, j, compTensor2 ).index();
            double val = uTensor2( gdofTensor2 );

            size_type gdofComp = dofComp->localToGlobal( elt, j, 0 ).index();
            uComp->set( gdofComp, val );
        }
    }
    return uComp;
}
} // namespace detail
#endif


#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    if (!this->isStandardModel())
        return;

    std::string modelName = "solid";
    bool hasMeasure = false;

    // volume variation
    for ( auto const& ppvv : M_postProcessVolumeVariation )
    {
        std::string const& vvname = ppvv.first;
        auto const& vvmarkers = ppvv.second;
        elements_reference_wrapper_t<mesh_type> vvrange = ( vvmarkers.size() == 1 && vvmarkers.begin()->empty() )?
            M_rangeMeshElements : markedelements( this->mesh(),vvmarkers );
        double volVar = this->computeVolumeVariation( vvrange );
        this->postProcessMeasuresIO().setMeasure( vvname, volVar );
        hasMeasure = true;
    }

    std::set<std::string> fieldNameStressScalar = { "Von-Mises","Tresca","princial-stress-1","princial-stress-2","princial-stress-3",
                                                    "stress_xx","stress_xy","stress_xz","stress_yx","stress_yy","stress_yz","stress_zx","stress_zy","stress_zz" };
    // points evaluation
    // this->modelProperties().parameters().updateParameterValues();
    // auto paramValues = this->modelProperties().parameters().toParameterValues();
    // this->modelProperties().postProcess().setParameterValues( paramValues );
    for ( auto const& evalPoints : this->modelProperties().postProcess().measuresPoint( modelName ) )
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
            if ( field == "displacement" || field == "velocity" || field == "acceleration" )
            {
                std::string ptNameExport = (boost::format("%1%_%2%")%field %ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId(field,ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextDisplacement->replace( ptIdInCtx, ptCoord );
            }
            else if ( field == "pressure" )
            {
                if ( !M_useDisplacementPressureFormulation )
                    continue;
                std::string ptNameExport = (boost::format("pressure_%1%")%ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId("pressure",ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextPressure->replace( ptIdInCtx, ptCoord );
            }
            else if ( fieldNameStressScalar.find( field ) != fieldNameStressScalar.end() )
            {
                std::string ptNameExport = (boost::format("%1%_%2%")%field %ptPos.name()).str();
                int ptIdInCtx = this->postProcessMeasuresEvaluatorContext().ctxId(field,ptNameExport);
                if ( ptIdInCtx >= 0 )
                    M_postProcessMeasuresContextStressScalar->replace( ptIdInCtx, ptCoord );
            }
        }
    }

    if ( M_postProcessMeasuresContextDisplacement )
    {
        Eigen::Matrix<value_type,Eigen::Dynamic,1> evalAtNodes;
        for ( std::string const& field : std::vector<std::string>( { "displacement","velocity","acceleration" } ) )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has(field) ) continue;

            if ( field == "displacement")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextDisplacement,
                                                   _expr=idv(this->fieldDisplacement()) );
            else if ( field == "velocity")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextDisplacement,
                                                   _expr=idv(this->fieldVelocity()) );
            else if ( field == "acceleration")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextDisplacement,
                                                   _expr=idv(this->fieldAcceleration()) );
            for ( int ctxId=0;ctxId<M_postProcessMeasuresContextDisplacement->nPoints();++ctxId )
            {
                if ( !this->postProcessMeasuresEvaluatorContext().has( field, ctxId ) ) continue;
                std::string const& ptNameExport = this->postProcessMeasuresEvaluatorContext().name( field,ctxId );
                std::vector<double> vecValues = { evalAtNodes( ctxId*nDim ) };
                if ( nDim > 1 ) vecValues.push_back( evalAtNodes( ctxId*nDim+1 ) );
                if ( nDim > 2 ) vecValues.push_back( evalAtNodes( ctxId*nDim+2 ) );
                this->postProcessMeasuresIO().setMeasureComp( ptNameExport, vecValues );
                //std::cout << "export point " << ptNameExport << " with node " << M_postProcessEvalPointDisplacement->node( ctxId ) << "\n";
                hasMeasure = true;
            }
        }
    }
    if ( M_useDisplacementPressureFormulation && M_postProcessMeasuresContextPressure && this->postProcessMeasuresEvaluatorContext().has("pressure") )
    {
        auto evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextPressure,
                                                _expr=idv(this->fieldPressure()) );
        for ( int ctxId=0;ctxId<M_postProcessMeasuresContextPressure->nPoints();++ctxId )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has( "pressure", ctxId ) ) continue;
            std::string ptNameExport = this->postProcessMeasuresEvaluatorContext().name( "pressure",ctxId );
            this->postProcessMeasuresIO().setMeasure( ptNameExport, evalAtNodes( ctxId ) );
            //std::cout << "export point " << ptNameExport << " with node " << M_postProcessEvalPointPressure->node( ctxId ) << "\n";
            hasMeasure = true;
        }
    }
    if ( M_postProcessMeasuresContextStressScalar )
    {
        this->updateStressCriterions();
        Eigen::Matrix<value_type,Eigen::Dynamic,1> evalAtNodes;
        for ( std::string const& field : fieldNameStressScalar )
        {
            if ( !this->postProcessMeasuresEvaluatorContext().has(field) ) continue;

            if ( field == "Von-Mises")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldVonMisesCriterions()) );
            else if ( field == "Tresca")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldTrescaCriterions()) );
            else if ( field == "princial-stress-1")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldPrincipalStresses(0)) );
            else if ( field == "princial-stress-2")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldPrincipalStresses(1)) );
            else if ( field == "princial-stress-3" && nDim == 3 )
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldPrincipalStresses(2)) );
            else if ( field == "stress_xx")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::X,Component::X)) );
            else if ( field == "stress_xy")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::X,Component::Y)) );
            else if ( field == "stress_xz" && nDim == 3 )
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::X,Component::Z)) );
            else if ( field == "stress_yx")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Y,Component::X)) );
            else if ( field == "stress_yy")
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Y,Component::Y)) );
            else if ( field == "stress_yz" && nDim == 3)
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Y,Component::Z)) );
            else if ( field == "stress_zx" && nDim == 3)
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Z,Component::X)) );
            else if ( field == "stress_zy" && nDim == 3)
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Z,Component::Y)) );
            else if ( field == "stress_zz" && nDim == 3)
                evalAtNodes = evaluateFromContext( _context=*M_postProcessMeasuresContextStressScalar,
                                                   _expr=idv(this->fieldStressTensor().comp(Component::Z,Component::Z)) );
            for ( int ctxId=0;ctxId<M_postProcessMeasuresContextStressScalar->nPoints();++ctxId )
            {
                if ( !this->postProcessMeasuresEvaluatorContext().has( field, ctxId ) ) continue;
                std::string const& ptNameExport = this->postProcessMeasuresEvaluatorContext().name( field,ctxId );
                this->postProcessMeasuresIO().setMeasure( ptNameExport, evalAtNodes( ctxId ) );
                hasMeasure = true;
            }

        }
    }


    auto mfields = this->modelFields();
    auto symbolsExpr = this->symbolsExpr( mfields );
    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, symbolsExpr, mfields );
    bool hasMeasurePoint = false;//this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, mfields );
    if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
        hasMeasure = true;

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }

}
#endif
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::startTimeStep()
{
    this->log("SolidMechanics","startTimeStep", "start" );

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    // go to next time step : ti + \Delta t
    if ( this->hasSolidEquationStandard() )
    {
        if ( M_timeStepping == "Newmark" )
        {
            // start time step
            if ( !this->doRestart() )
                M_timeStepNewmark->start(*M_fieldDisplacement);
        }
        else if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
        {
            // start time step
            if ( !this->doRestart() )
            {
                M_timeStepBdfDisplacement->start( *M_fieldDisplacement );
                M_timeStepBdfVelocity->start( *M_fieldVelocity );
            }
        }
        // start save pressure
        if ( this->hasDisplacementPressureFormulation() && !this->doRestart() )
            M_savetsPressure->start( *M_fieldPressure );
    }

    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->startTimeStep();

    // up current time
    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();

    this->log("SolidMechanics","startTimeStep", "finish" );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("SolidMechanics","updateTimeStep", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // some time stepping require to compute residual without time derivative
    this->updateTimeStepCurrentResidual();

    // go to next time step
    if (this->hasSolidEquationStandard())
    {
        if ( M_timeStepping == "Newmark" )
        {
            // next time step
            M_timeStepNewmark->next( *M_fieldDisplacement );
        }
        else if ( M_timeStepping == "BDF" || M_timeStepping == "Theta" )
        {
            M_timeStepBdfDisplacement->next( *M_fieldDisplacement );
            M_timeStepBdfVelocity->next( *M_fieldVelocity );
        }

        if ( this->hasDisplacementPressureFormulation() )
            M_savetsPressure->next(*M_fieldPressure);
    }

    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->updateTimeStep();

    // up current time
    this->updateTime( this->timeStepBase()->time() );

    this->updateParameterValues();

    this->timerTool("TimeStepping").stop("updateTimeStep");

    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("SolidMechanics","updateTimeStep", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateTimeStepCurrentResidual()
{
    if ( this->hasSolidEquationStandard() )
    {
        auto algebraicFactory = this->algebraicFactory();
        if ( !algebraicFactory || M_timeStepping == "BDF" || M_timeStepping == "Newmark" )
            return;

        if ( M_timeStepping == "Theta" )
        {
            M_timeStepThetaSchemePreviousContrib->zero();
            this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
            std::vector<std::string> infos = { prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") };
            algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", false );
            algebraicFactory->evaluateResidual( this->algebraicBlockVectorSolution()->vectorMonolithic(), M_timeStepThetaSchemePreviousContrib, infos, false );
            algebraicFactory->setActivationAddVectorResidualAssembly( "Theta-Time-Stepping-Previous-Contrib", true );
        }
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::predictorDispl()
{
    this->log("SolidMechanics","predictorDispl", "start" );

    if (this->hasSolidEquationStandard())
    {
        if ( M_timeStepping == "Newmark" )
        {
            //std::string mytype = soption(_name="predictor-disp-type",_prefix=this->prefix());
            double dt = M_timeStepNewmark->timeStep();
            if (M_timeStepNewmark->iteration() == 1) //order 1: \tilde{u} = u_n + dt* \dot{u}_n
            {
                //this->fieldDisplacement().add(M_timeStepNewmark->timeStep(),M_timeStepNewmark->currentVelocity());
                this->fieldDisplacement().add( dt, M_timeStepNewmark->previousVelocity(0));
            }
            else //order 2: \tilde{u} = u_n + dt*( (3/2)*\dot{u}_n-(1/2)*\dot{u}_{n-1}
            {
                //this->fieldDisplacement().add((3./2.)*M_timeStepNewmark->timeStep(), M_timeStepNewmark->currentVelocity() );
                //this->fieldDisplacement().add((-1./2.)*M_timeStepNewmark->timeStep(), M_timeStepNewmark->previousVelocity());
                this->fieldDisplacement().add( (3./2.)*dt, M_timeStepNewmark->previousVelocity(0) );
                this->fieldDisplacement().add( (-1./2.)*dt, M_timeStepNewmark->previousVelocity(1) );
            }
        }
        else
        {
            *M_fieldDisplacement = M_timeStepBdfDisplacement->poly();
            *M_fieldVelocity = M_timeStepBdfVelocity->poly();
        }
    }

    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->predictorDispl();

    this->log("SolidMechanics","predictorDispl", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateParameterValues()
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
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::setParameterValues( std::map<std::string,double> const& paramValues )
{
    if ( this->manageParameterValuesOfModelProperties() )
    {
        this->modelProperties().parameters().setParameterValues( paramValues );
        this->modelProperties().postProcess().setParameterValues( paramValues );
        this->materialsProperties()->setParameterValues( paramValues );
    }

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        physicData->setParameterValues( paramValues );

    M_boundaryConditions.setParameterValues( paramValues );


#if 0 // VINCENT
    this->M_bcDirichlet.setParameterValues( paramValues );
    for ( auto & bcDirComp : this->M_bcDirichletComponents )
        bcDirComp.second.setParameterValues( paramValues );
    this->M_bcNeumannScalar.setParameterValues( paramValues );
    this->M_bcNeumannVectorial.setParameterValues( paramValues );
    this->M_bcNeumannTensor2.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameScalar.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameVectorial.setParameterValues( paramValues );
    this->M_bcNeumannEulerianFrameTensor2.setParameterValues( paramValues );
    this->M_bcRobin.setParameterValues( paramValues );
    this->M_volumicForcesProperties.setParameterValues( paramValues );
#endif
    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->setParameterValues( paramValues );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::solve( bool upVelAcc )
{
    this->log("SolidMechanics","solve", "start" );
    this->timerTool("Solve").start();

    //this->updateParameterValues();

    this->algebraicBlockVectorSolution()->updateVectorFromSubVectors();
    this->algebraicFactory()->solve( M_solverName, this->algebraicBlockVectorSolution()->vectorMonolithic() );
    this->algebraicBlockVectorSolution()->localize();

    if ( upVelAcc && !this->isStationary() )
        this->updateVelocity();

    double tElapsed = this->timerTool("Solve").stop("solve");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("Solve").setAdditionalParameter("time",this->currentTime());
        this->timerTool("Solve").save();
    }

    this->log("SolidMechanics","solve", (boost::format("finish in %1% s")%tElapsed).str() );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateVelocity()
{
    if ( M_timeStepping != "Newmark" )
        return;

    this->log("SolidMechanics","updateVelocityAndAcceleration", "start" );

    if (this->hasSolidEquationStandard())
    {
        M_timeStepNewmark->updateFromDisp(*M_fieldDisplacement);
    }

    if ( this->hasSolidEquation1dReduced() )
    {
        M_solid1dReduced->updateVelocity();
    }

    this->log("SolidMechanics","updateVelocityAndAcceleration", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNormalStressFromStruct()
{
#if 0 // TODO VINCENT
    if ( !M_XhNormalStress )
        this->createAdditionalFunctionSpacesNormalStress();

    auto const& u = this->fieldDisplacement();
    auto range = boundaryfaces(this->mesh());
    //auto range = markedfaces(this->mesh(),this->markerNameFSI());


    if ( M_modelName == "Elasticity" )
    {
        auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
        auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

        auto const Id = eye<nDim,nDim>();
        auto eps = sym(gradv(u));//0.5*(gradv(u)+trans(gradv(u)));
        if ( !this->hasDisplacementPressureFormulation() )
        {
            auto sigma = idv(coeffLame1)*trace(eps)*Id + 2*idv(coeffLame2)*eps;
            M_fieldNormalStressFromStruct->on( _range=range,
                                               _expr=sigma*N(),
                                               _geomap=this->geomap() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigma = idv(p)*Id + 2*idv(coeffLame2)*eps;
            M_fieldNormalStressFromStruct->on( _range=range,
                                               _expr=sigma*N(),
                                               _geomap=this->geomap() );
        }
    }
    else if ( M_modelName=="Elasticity-Large-Deformation" )
    {
        CHECK( false ) << "TODO";
    }
    else if ( M_modelName == "Hyper-Elasticity" )
    {

        auto const sigma = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
        if ( !this->useDisplacementPressureFormulation() )
        {
            M_fieldNormalStressFromStruct->on( _range=range,
                                               _expr=sigma*N(),
                                               _geomap=this->geomap() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigmaWithPressure = Feel::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()) + sigma ;
            M_fieldNormalStressFromStruct->on( _range=range,
                                               _expr=sigmaWithPressure*N(),
                                               _geomap=this->geomap() );
        }
    }
#endif

}


//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateStressCriterions()
{
    this->log("SolidMechanics","updateStressCriterions", "start" );

    this->createAdditionalFunctionSpacesStressTensor();

    auto const& u = this->fieldDisplacement();
    if ( M_modelName == "Elasticity" )
    {
        auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
        auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

        auto const Id = eye<nDim,nDim>();
        auto eps = sym(gradv(u));//0.5*(gradv(u)+trans(gradv(u)));
        if ( !this->useDisplacementPressureFormulation() )
        {
            auto sigma = idv(coeffLame1)*trace(eps)*Id + 2*idv(coeffLame2)*eps;
            M_fieldStressTensor->on(_range=M_rangeMeshElements,_expr=sigma );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigma = idv(p)*Id + 2*idv(coeffLame2)*eps;
            M_fieldStressTensor->on(_range=M_rangeMeshElements,_expr=sigma );
        }
    }
    else if ( M_modelName == "Hyper-Elasticity" )
    {
#if 0 // TODO VINCENT
        auto sigma = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
        if ( !this->useDisplacementPressureFormulation() )
        {
            M_fieldStressTensor->on(_range=M_rangeMeshElements,_expr=sigma );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigmaWithPressure = Feel::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()) + sigma ;
            M_fieldStressTensor->on(_range=M_rangeMeshElements,_expr=sigmaWithPressure );
        }

#if 0
        auto Id = eye<nDim,nDim>();
        auto sigma_dev = sigma-(1./3.)*trace(sigma)*Id;
        M_fieldVonMisesCriterions->on(_range=M_rangeMeshElements,
                                      _expr=sqrt((3./2.)*inner(sigma_dev,sigma_dev, mpl::int_<InnerProperties::IS_SAME>() )) );
#endif
#endif
    }

    typedef Eigen::Matrix<double, nDim, nDim> matrixN_type;
    Eigen::EigenSolver< matrixN_type > eigenSolver;
    matrixN_type sigma_eigen_matrix;
    std::vector<double> eigenValuesSorted(nDim);

    auto dof = M_XhStressTensor->dof();
    for ( auto const& eltWrap : M_rangeMeshElements )
    {
        auto const& elt = boost::unwrap_ref( eltWrap );
        int nLocDofPerComp = dof->nLocalDof( true );
        for ( size_type j =0; j < nLocDofPerComp;++j )
        {
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for (uint16_type comp2=0; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                {
                    uint16_type comp = comp1*space_stress_tensor_type::nComponents2+comp2;
                    size_type gdof = dof->localToGlobal( elt, j, comp ).index();
                    double val = M_fieldStressTensor->operator()( gdof );
                    sigma_eigen_matrix(comp1,comp2) =val;
                }

            //compute eigenvalues
            eigenSolver.compute(sigma_eigen_matrix);

            for ( uint16_type comp=0; comp < space_stress_tensor_type::nComponents1;++comp )
                eigenValuesSorted[comp] = real(eigenSolver.eigenvalues()[comp]);
            std::sort(eigenValuesSorted.begin(), eigenValuesSorted.end(), std::greater/*less*/<double>() );

            // global dof id
            size_type gdofScal = M_XhStressTensor->compSpace()->dof()->localToGlobal( elt, j, 0 ).index();

            // update principal stress
            for ( uint16_type comp=0; comp < space_stress_tensor_type::nComponents1;++comp )
                M_fieldsPrincipalStresses[comp]->set( gdofScal, eigenValuesSorted[comp] );

#if 0
            // update invariants
            double invariant1=0;
            for ( uint16_type comp=0; comp < space_stress_tensor_type::nComponents1;++comp )
                invariant1 += eigenValuesSorted[comp];
            M_fieldsPrincipalStresses[0]->set( gdofScal, invariant1 );
            double invariant2=0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for ( uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    invariant2 += eigenValuesSorted[comp1]*eigenValuesSorted[comp2];
            M_fieldsPrincipalStresses[1]->set( gdofScal, invariant2 );
            double invariant3 = 1;
            for ( uint16_type comp=0; comp < space_stress_tensor_type::nComponents1;++comp )
                invariant3 *= eigenValuesSorted[comp];
            M_fieldsPrincipalStresses[2]->set( gdofScal, invariant3 );
#endif
            // update Tresca Criterions
            double resTresca = 0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for (uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    resTresca =  std::max( resTresca, std::abs( eigenValuesSorted[comp1] - eigenValuesSorted[comp2] ) );
            M_fieldTrescaCriterions->set( gdofScal, resTresca );

            // update Von Mises Criterions
            double resVonMises = 0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for (uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    resVonMises += (1./2.)*std::pow( eigenValuesSorted[comp1] - eigenValuesSorted[comp2], 2 );
            resVonMises = std::sqrt( resVonMises );
            M_fieldVonMisesCriterions->set( gdofScal, resVonMises );
        }
    }

    this->log("SolidMechanics","updateStressCriterions", "finish" );

}
#endif
//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::super_type::block_pattern_type
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::blockPattern() const
{
    if ( this->hasDisplacementPressureFormulation() )
        return BlocksStencilPattern(2,2) << size_type(Pattern::COUPLED) << size_type(Pattern::COUPLED)
                                         << size_type(Pattern::COUPLED) << size_type(Pattern::ZERO);
    else
        return BlocksStencilPattern(1,1) << size_type(Pattern::COUPLED);
}

//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::computeMaxDisp() const
{
    double res = 0;
    if(this->isStandardModel())
    {
        auto u = this->fieldDisplacement();
        auto dispMagnOnBoundary = vf::project(_space=functionSpaceDisplacement()->compSpace(),
                                              _range=M_rangeMeshElements,
                                              _expr=sqrt(trans(idv(u))*idv(u)) );
        res = dispMagnOnBoundary.max();
    }
    else if(this->is1dReducedModel())
    {
        res = this->fieldDisplacementScal1dReduced().max();
    }

    //if(this->isStandardModel()) res = this->fieldDisplacement().max();
    //else if(this->is1dReducedModel()) res = this->fieldDisplacementScal1dReduced().max();

    return res;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::computeMaxDispOnBoundary(std::string __marker) const
{
    double res = 0;
    auto u = this->fieldDisplacement();
    auto dispMagnOnBoundary = vf::project(_space=functionSpaceDisplacement()->compSpace(),
                                          _range=markedfaces(this->mesh(),__marker),
                                          _expr=sqrt(trans(idv(u))*idv(u)) );
    res = dispMagnOnBoundary.max();

    //if(this->isStandardModel()) res = this->fieldDisplacement().max();
    //else if(this->is1dReducedModel()) res = this->fieldDisplacementScal1dReduced().max();

    return res;
}
#endif
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::computeExtremumValue( std::string const& field, std::set<std::string> const& markers, std::string const& type ) const
{
    if(!this->hasSolidEquationStandard()) return 0.;
    CHECK( type == "max" || type == "min" ) << "invalid type " << type;

    if ( field == "displacement" || field == "velocity" || field == "acceleration" )
    {
        auto fieldMagnitude = functionSpaceDisplacement()->compSpace()->elementPtr();
        if ( markers.empty() )
        {
            if ( field == "displacement" )
                fieldMagnitude->on(_range=M_rangeMeshElements, _expr=norm2( idv(this->fieldDisplacement()) ) );
            else if ( field == "velocity" )
                fieldMagnitude->on(_range=M_rangeMeshElements, _expr=norm2( idv(this->fieldVelocity()) ) );
            else if ( field == "acceleration" )
                fieldMagnitude->on(_range=M_rangeMeshElements, _expr=norm2( idv(this->fieldAcceleration()) ) );
        }
        else
        {
            if ( this->mesh()->hasFaceMarker( *markers.begin() ) )
            {
                for ( std::string const& marker : markers )
                    CHECK( this->mesh()->hasFaceMarker( marker ) ) << "marker list must be same type (here marker : " << marker << " must be a face marker)";
                if ( field == "displacement" )
                    fieldMagnitude->on(_range=markedfaces(this->mesh(),markers), _expr=norm2( idv(this->fieldDisplacement()) ) );
                else if ( field == "velocity" )
                    fieldMagnitude->on(_range=markedfaces(this->mesh(),markers), _expr=norm2( idv(this->fieldVelocity()) ) );
                else if ( field == "acceleration" )
                    fieldMagnitude->on(_range=markedfaces(this->mesh(),markers), _expr=norm2( idv(this->fieldAcceleration()) ) );
            }
        }
        if ( type == "max" )
            return fieldMagnitude->max();
        else
            return fieldMagnitude->min();
    }
    else if ( field == "pressure" )
    {
        CHECK( this->hasDisplacementPressureFormulation() ) << "model does not take into account the pressure";
        CHECK( false ) << "TODO pressure max/min";
    }

    //if(this->isStandardModel()) res = this->fieldDisplacement().max();
    //else if(this->is1dReducedModel()) res = this->fieldDisplacementScal1dReduced().max();

    return 0.;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::computeVolumeVariation( elements_reference_wrapper_t<mesh_type> const& rangeElt ) const
{
    //using namespace Feel::vf;
    if(!this->hasSolidEquationStandard()) return 0.;

    auto const Id = eye<nDim,nDim>();
    auto const& u = this->fieldDisplacement();
    auto Fv = Id + gradv(u);
    auto detFv = det(Fv);
    double newArea = integrate(_range=rangeElt,
                               _expr=detFv,
                               _geomap=this->geomap() ).evaluate()(0,0);
    double refAera = integrate(_range=rangeElt,
                               _expr=cst(1.),
                               _geomap=this->geomap() ).evaluate()(0,0);

    return (newArea-refAera)/refAera;
}

//---------------------------------------------------------------------------------------------------//



SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateMassMatrixLumped()
{
    CHECK ( this->hasSolidEquationStandard() ) << "only compute when isStandardModel";
    auto Vh = this->functionSpaceDisplacement();
    auto const& u = this->fieldDisplacement();
    auto mesh = Vh->mesh();
    auto se = this->symbolsExpr();
    // mass matrix of Vh
    auto massMatrix = this->backend()->newMatrix(_test=Vh,_trial=Vh);
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& densityProp = this->materialsProperties()->density( matName );
            auto densityExpr = expr( densityProp.expr(), se );
            form2( _trial=Vh, _test=Vh,_matrix=massMatrix) =
                integrate(_range=M_rangeMeshElements,
                          _expr=densityExpr*inner(idt(u),id(u)) );
        }
    }
    massMatrix->close();

    // mass matrix lumped
    auto graph = std::make_shared<graph_type>( Vh->dof(),Vh->dof() );
    graph->addMissingZeroEntriesDiagonal();
    graph->close();
    M_massMatrixLumped = this->backend()->newMatrix( 0,0,0,0,graph );
    //M_massMatrixLumped = this->backend()->newMatrix(_test=Vh,_trial=Vh);
    if ( Vh->fe()->nOrder==1)
    {
        M_vecDiagMassMatrixLumped = this->backend()->newVector(Vh);
        auto unityVec = this->backend()->newVector(Vh);
        unityVec->setConstant(1);
        massMatrix->multVector( unityVec,M_vecDiagMassMatrixLumped );
        M_massMatrixLumped->setDiagonal( M_vecDiagMassMatrixLumped );
    }
    else
    {
        // ref for high order :
        // - https://www.sharcnet.ca/Software/Ansys/17.2/en-us/help/ans_thry/thy_et2.html
        // - https://scicomp.stackexchange.com/questions/19704/how-to-formulate-lumped-mass-matrix-in-fem
        // - section 16.2.4 of : Zhu, J., Z. R. L. Taylor, and O. C. Zienkiewicz. "The finite element method: its basis and fundamentals." (2005): 54-102.
        auto sumRow = this->backend()->newVector(Vh);
        auto unityVec = this->backend()->newVector(Vh);
        unityVec->setConstant(1);
        massMatrix->multVector( unityVec,sumRow );
        double sumMatrix = sumRow->sum();
        M_vecDiagMassMatrixLumped = massMatrix->diagonal();
        double sumDiag = M_vecDiagMassMatrixLumped->sum();
        M_vecDiagMassMatrixLumped->scale( sumMatrix/sumDiag );
        M_massMatrixLumped->setDiagonal( M_vecDiagMassMatrixLumped );
    }
}

} // FeelModels


} // Feel




