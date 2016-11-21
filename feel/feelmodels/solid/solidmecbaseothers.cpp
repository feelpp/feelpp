/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */


#include <feel/feelmodels/solid/solidmecbase.hpp>
#include <feel/feelmodels/modelvf/solidmecfirstpiolakirchhoff.hpp>
#include <feel/feelmodels/modelvf/solidmecincompressibility.hpp>


namespace Feel
{
namespace FeelModels
{

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::restartExporters( double time )
{
    // restart exporter
    if (this->doRestart() && this->restartPath().empty())
    {
        if ( this->isStandardModel() )
        {
            if (!M_isHOVisu)
            {
                if ( M_exporter && M_exporter->doExport() )
                    M_exporter->restart(this->timeInitial());
            }
            else
            {
#if 1 // defined(FEELPP_HAS_VTK)
                if ( M_exporter_ho && M_exporter_ho->doExport() )
                    M_exporter_ho->restart(this->timeInitial());
#endif
            }
        }
        else
        {
            if ( M_exporter_1dReduced && M_exporter_1dReduced->doExport() )
                M_exporter_1dReduced->restart(this->timeInitial());
        }
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    this->log("SolidMechanics","getInfo", "start" );

    std::string StateTemporal = (this->isStationary())? "Stationary" : "Transient";
    size_type nElt,nDof;
    if (isStandardModel()) {nElt=M_mesh->numGlobalElements(); nDof=M_XhDisplacement->nDof();}
    else {nElt=M_mesh_1dReduced->numGlobalElements(); nDof=M_Xh_1dReduced->nDof();}

    std::string SchemaTimeType = soption(_name="time-schema",_prefix=this->prefix());
    std::string ResartMode;
    if (this->doRestart()) ResartMode = "Yes";
    else ResartMode = "No";

    std::string hovisuMode,myexporterType;
    int myexporterFreq=1;
    if ( this->isStandardModel() )
    {
        if (M_isHOVisu)
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
        else
        {
            hovisuMode = "OFF";
            myexporterType = M_exporter->type();
            myexporterFreq = M_exporter->freq();
        }
    }
    else
    {
        hovisuMode = "OFF";
        myexporterType = M_exporter_1dReduced->type();
        myexporterFreq = M_exporter_1dReduced->freq();
    }

    std::string doExport_str;
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Displacement ) )
        doExport_str=(doExport_str.empty())?"displacement":doExport_str+" - displacement";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Velocity ) )
        doExport_str=(doExport_str.empty())?"velocity":doExport_str+" - velocity";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Acceleration ) )
        doExport_str=(doExport_str.empty())?"acceleration":doExport_str+" - acceleration";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pressure ) )
        doExport_str=(doExport_str.empty())?"pressure":doExport_str+" - pressure";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
        doExport_str=(doExport_str.empty())?"normal stress":doExport_str+" - normal stress";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::MaterialProperties ) )
        doExport_str=(doExport_str.empty())?"material properties":doExport_str+" - material properties";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::VonMises ) )
        doExport_str=(doExport_str.empty())?"fsi":doExport_str+" - Von Mises";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Tresca ) )
        doExport_str=(doExport_str.empty())?"fsi":doExport_str+" - Tresca";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::PrincipalStresses ) )
        doExport_str=(doExport_str.empty())?"fsi":doExport_str+" - principal stresses";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::FSI ) )
        doExport_str=(doExport_str.empty())?"fsi":doExport_str+" - fsi";
    if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pid ) )
        doExport_str=(doExport_str.empty())?"pid":doExport_str+" - pid";
    for ( std::string const& userFieldName : M_postProcessUserFieldExported )
        doExport_str=(doExport_str.empty())?userFieldName:doExport_str+" - "+userFieldName;

    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n||==============================================||"
           << "\n||----------Info : SolidMechanics---------------||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Root Repository : " << this->rootRepository()
           << "\n   Physical Model"
           << "\n     -- pde name : " << M_pdeType
           << "\n     -- material law : " << this->mechanicalProperties()->materialLaw()
           << "\n     -- use displacement-pressure formulation : " << std::boolalpha << M_useDisplacementPressureFormulation
           << "\n     -- time mode : " << StateTemporal;
    *_ostr << this->mechanicalProperties()->getInfoMaterialParameters()->str();
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoNeumannEulerianFrameBC()
           << this->getInfoRobinBC()
           << this->getInfoFluidStructureInterfaceBC();
    *_ostr << "\n   Space Discretization";
    if ( this->hasGeofileStr() )
        *_ostr << "\n     -- geo file name   : " << this->geofileStr();
    *_ostr << "\n     -- msh file name   : " << this->mshfileStr()
           << "\n     -- nb elt in mesh : " << nElt
           << "\n     -- nb dof (displacement) : " << nDof
           << "\n     -- polynomial order : " << nOrder;
    if ( M_useDisplacementPressureFormulation )
        *_ostr << "\n     -- nb dof (pressure) : " << M_XhPressure->nDof();
    *_ostr << "\n     -- mechanical properties : "
           << (boost::format("P%1%%2%")%mechanicalproperties_type::space_type::basis_type::nOrder %std::string( (use_continous_mechanical_properties)? "c":"d")).str();
    *_ostr << "\n   Time Discretization"
           << "\n     -- initial time : " << this->timeStepBase()->timeInitial()
           << "\n     -- final time   : " << this->timeStepBase()->timeFinal()
           << "\n     -- time step    : " << this->timeStepBase()->timeStep()
           << "\n     -- type : " << SchemaTimeType;
    if ( SchemaTimeType == "Newmark" )
    {
        if (this->isStandardModel())
            *_ostr << " ( gamma="<< this->timeStepNewmark()->gamma() <<", beta="<< this->timeStepNewmark()->beta() << " )";
        else if (this->is1dReducedModel())
            *_ostr << " ( gamma="<< this->timeStepNewmark1dReduced()->gamma() <<", beta="<< this->timeStepNewmark1dReduced()->beta() << " )";
    }
    *_ostr << "\n     -- restart mode : " << ResartMode
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
           << "\n   Numerical Solver"
           << "\n     -- solver : " << M_pdeSolver;
    if (this->isStandardModel() && M_algebraicFactory)
        *_ostr << M_algebraicFactory->getInfo()->str();
    else if (this->is1dReducedModel() && M_algebraicFactory_1dReduced)
        *_ostr << M_algebraicFactory_1dReduced->getInfo()->str();
    *_ostr << "\n||==============================================||"
           << "\n";

    this->log("SolidMechanics","getInfo", "finish" );

    return _ostr;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeType(std::string __type)
{
    if ( __type == "Elasticity" )                      { M_pdeType="Elasticity"; M_pdeSolver="LinearSystem"; }
    else if ( __type == "Elasticity-Large-Deformation" ) { M_pdeType="Elasticity-Large-Deformation"; M_pdeSolver="Newton"; }
    else if ( __type == "Hyper-Elasticity" )           { M_pdeType="Hyper-Elasticity"; M_pdeSolver="Newton"; }
    else if ( __type == "Generalised-String" )         { M_pdeType="Generalised-String"; M_pdeSolver="LinearSystem"; }
    else
        CHECK( false ) << "invalid pdeType "<< __type << "\n";

    //if ( M_useDisplacementPressureFormulation && M_pdeSolver != "Newton" )
    //    M_pdeSolver="Newton";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::pdeSolver(std::string __type)
{
    if ( __type == "LinearSystem" )
        M_pdeSolver="LinearSystem";
    else if ( __type == "Newton" )
        M_pdeSolver="Newton";
    else
        CHECK( false ) << "invalid pdeSolver " << __type << "\n";
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
int
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::nBlockMatrixGraph() const
{
    int nBlock = 1;
    if ( M_useDisplacementPressureFormulation )
        ++nBlock;
    return nBlock;
}

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    this->log("SolidMechanics","buildBlockMatrixGraph", "start" );
    int nBlock = this->nBlockMatrixGraph();

    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    int indexBlock=0;
    myblockGraph(indexBlock,indexBlock) =stencil(_test=this->functionSpace(),_trial=this->functionSpace(),
                                                 //_pattern_block=this->blockPattern(),
                                                 _diag_is_nonzero=(nBlock==1),
                                                 _close=(nBlock==1) )->graph();
    ++indexBlock;

    if ( M_useDisplacementPressureFormulation )
    {
        myblockGraph(indexBlock,0) = stencil(_test=M_XhPressure,_trial=M_XhDisplacement,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=M_XhDisplacement,_trial=M_XhPressure,
                                             _diag_is_nonzero=false,_close=false)->graph();
        //if ( M_pdeSolver=="LinearSystem" )
        myblockGraph(indexBlock,indexBlock) = stencil(_test=M_XhPressure,_trial=M_XhPressure,
                                                      _diag_is_nonzero=false,_close=false)->graph();
        ++indexBlock;
    }

    myblockGraph.close();

    this->log("SolidMechanics","buildBlockMatrixGraph", "finish" );
    return myblockGraph;
}

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::graph_ptrtype
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::buildMatrixGraph() const
{
    if (this->isStandardModel())
    {
        auto blockGraph = this->buildBlockMatrixGraph();
        if ( blockGraph.nRow() == 1 && blockGraph.nCol() == 1 )
            return blockGraph(0,0);
        else
            return graph_ptrtype( new graph_type( blockGraph ) );
    }
    else //if (this->is1dReducedModel())
    {
        return graph_ptrtype();
    }
}

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateBlockVectorSolution()
{
    if (this->isStandardModel())
    {
        M_blockVectorSolution.setVector( *M_blockVectorSolution.vector(), this->fieldDisplacement(), 0 );

        if ( M_useDisplacementPressureFormulation )
        {
            size_type blockIndexPressure = this->startBlockIndexFieldsInMatrix().find("pressure")->second;
            M_blockVectorSolution.setVector( *M_blockVectorSolution.vector(), this->fieldPressure(), blockIndexPressure );
        }
        M_blockVectorSolution.vector()->close();
    }
    else if (this->is1dReducedModel())
    {
        // nothing because do in solve (TODO mv here)
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    this->log("SolidMechanics","exportResults",(boost::format("start at time %1%")%time).str() );
    this->timerTool("PostProcessing").start();

    if (!M_isHOVisu)
        this->exportFieldsImpl( time );
    else
        this->exportFieldsImplHO( time );

    this->exportMeasures( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }

    this->log("SolidMechanics","exportResults", "finish" );

} // SolidMechanics::export
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportFieldsImpl( double time )
{
    if (this->isStandardModel())
    {
        if ( !M_exporter->doExport() ) return;
        //M_exporter->step( time )->setMesh( M_mesh );
        bool hasFieldToExport = false;
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pid ) )
        {
            M_exporter->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Displacement ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"displacement"), this->fieldDisplacement() );
            hasFieldToExport = true;
        }
        if ( M_useDisplacementPressureFormulation && this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pressure ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"pressure"), *M_fieldPressure );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Velocity ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity"), this->timeStepNewmark()->currentVelocity() );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Acceleration ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"acceleration"), this->timeStepNewmark()->currentAcceleration() );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
        {
            this->updateNormalStressFromStruct();
            M_exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"), *M_fieldNormalStressFromStruct );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::VonMises ) ||
             this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Tresca ) ||
             this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::PrincipalStresses ) )
        {
            this->updateStressCriterions();
            if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::VonMises ) )
                M_exporter->step( time )->add( prefixvm(this->prefix(),"von-mises-criterions"), *M_fieldVonMisesCriterions );
            if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Tresca ) )
                M_exporter->step( time )->add( prefixvm(this->prefix(),"tresca-criterions"), *M_fieldTrescaCriterions );
            if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::PrincipalStresses ) )
                for (int d=0;d<M_fieldsPrincipalStresses.size();++d)
                    M_exporter->step( time )->add( prefixvm(this->prefix(),(boost::format("princial-stress-%1%")%d).str() ), *M_fieldsPrincipalStresses[d] );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::MaterialProperties ) )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"YoungModulus"), this->mechanicalProperties()->fieldYoungModulus() );
            M_exporter->step( time )->add( prefixvm(this->prefix(),"PoissonCoefficient"), this->mechanicalProperties()->fieldCoeffPoisson() );
            M_exporter->step( time )->add( prefixvm(this->prefix(),"Lame1"), this->mechanicalProperties()->fieldCoeffLame1() );
            M_exporter->step( time )->add( prefixvm(this->prefix(),"Lame2"), this->mechanicalProperties()->fieldCoeffLame2() );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::FSI) && this->fieldVelocityInterfaceFromFluidPtr() )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity-interface-from-fluid"), this->fieldVelocityInterfaceFromFluid() );
            hasFieldToExport = true;
        }
        for ( std::string const& userFieldName : M_postProcessUserFieldExported )
        {
            if ( this->hasFieldUserScalar( userFieldName ) )
            {
                M_exporter->step( time )->add( prefixvm(this->prefix(),userFieldName), this->fieldUserScalar( userFieldName ) );
                hasFieldToExport = true;
            }
            else if ( this->hasFieldUserVectorial( userFieldName ) )
            {
                M_exporter->step( time )->add( prefixvm(this->prefix(),userFieldName), this->fieldUserVectorial( userFieldName ) );
                hasFieldToExport = true;
            }
        }

        //M_exporter->step( time )->add( "prestress", *U_displ_struct_prestress);
        if ( hasFieldToExport )
            M_exporter->save();
    }
    else
    {
        if ( !M_exporter_1dReduced->doExport() ) return;
        //M_exporter_1dReduced->step( time )->setMesh( M_mesh_1dReduced );
        bool hasFieldToExport = false;
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pid ) )
        {
            M_exporter_1dReduced->step( time )->addRegions( prefixvm(this->prefix(),"1d-reduced"),
                                                            prefixvm(this->prefix(),prefixvm(this->subPrefix(),"1d-reduced") ) );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Displacement ) )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-displacement"), *M_disp_1dReduced );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Velocity ) )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-velocity"), this->timeStepNewmark1dReduced()->currentVelocity()  /**M_velocity_1dReduced*/ );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Acceleration ) )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-acceleration"), this->timeStepNewmark1dReduced()->currentAcceleration() );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-stress"),*M_stress_1dReduced );
            hasFieldToExport = true;
        }
        if ( false )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-stressvec"),*M_stress_vect_1dReduced);
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-stress"),*M_stress_1dReduced);
            hasFieldToExport = true;
        }
        if ( hasFieldToExport )
            M_exporter_1dReduced->save();
    }
}

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportFieldsImplHO( double time )
{
    if (this->isStandardModel())
    {
#if 1 //defined(FEELPP_HAS_VTK)
        if ( !M_exporter_ho->doExport() ) return;
        //M_exporter_ho->step( time )->setMesh( M_displacementVisuHO->mesh() );

        bool hasFieldToExport = false;
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pid ) )
        {
            M_exporter_ho->step( time )->addRegions( this->prefix(), this->subPrefix().empty()? this->prefix() : prefixvm(this->prefix(),this->subPrefix()) );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Displacement ) )
        {
            M_opIdisplacement->apply(this->fieldDisplacement(),*M_displacementVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"displacement-ho"), *M_displacementVisuHO );
            hasFieldToExport = true;
        }
        if ( M_useDisplacementPressureFormulation && this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Pressure ) )
        {
            M_opIpressure->apply(*M_fieldPressure,*M_pressureVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure-ho"), *M_pressureVisuHO );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Velocity ) )
        {
            auto velocityVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentVelocity(),velocityVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity-ho"), velocityVisuHO );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::Acceleration ) )
        {
            auto accelerationVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentAcceleration(),accelerationVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"acceleration-ho"), accelerationVisuHO );
            hasFieldToExport = true;
        }
        if ( this->hasPostProcessFieldExported( SolidMechanicsPostProcessFieldExported::NormalStress ) )
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

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportMeasures( double time )
{
    if (!this->isStandardModel())
        return;
    bool hasMeasure = false;

    // volume variation
    auto itFindMeasures = this->modelProperties().postProcess().find("Measures");
    if ( itFindMeasures != this->modelProperties().postProcess().end() )
    {
        if ( std::find( itFindMeasures->second.begin(), itFindMeasures->second.end(), "VolumeVariation" ) != itFindMeasures->second.end() )
        {
            double volVar = this->computeVolumeVariation();
            this->postProcessMeasuresIO().setMeasure( "volume_variation", volVar );
            hasMeasure = true;
        }
    }

    std::set<std::string> fieldNameStressScalar = { "Von-Mises","Tresca","princial-stress-1","princial-stress-2","princial-stress-3",
                                                    "stress_xx","stress_xy","stress_xz","stress_yx","stress_yy","stress_yz","stress_zx","stress_zy","stress_zz" };
    // points evaluation
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


    // extremum evaluation
    for ( auto const& measureExtremum : this->modelProperties().postProcess().measuresExtremum() )
    {
        auto const& fields = measureExtremum.fields();
        std::string const& name = measureExtremum.extremum().name();
        std::string const& type = measureExtremum.extremum().type();
        auto const& meshMarkers = measureExtremum.extremum().meshMarkers();
        for ( std::string const& field : fields )
        {
            double val = this->computeExtremumValue( field, meshMarkers, type );
            if ( field == "displacement" || field == "velocity" || field == "acceleration" )
            {
                std::string nameExport = (boost::format("%1%_magnitude_%2%_%3%")%field %type %name).str();
                this->postProcessMeasuresIO().setMeasure( nameExport,val );
                hasMeasure = true;
            }
        }
    }


    if ( hasMeasure )
    {
        this->postProcessMeasuresIO().setParameter( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
    }

}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateTimeStep()
{
    this->log("SolidMechanics","updateTimeStep", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    if (this->isStandardModel())
    {
        // next time step
        M_timeStepNewmark->next( *M_fieldDisplacement );
        if ( M_useDisplacementPressureFormulation ) M_savetsPressure->next(*M_fieldPressure);
        // up current time
        this->updateTime( M_timeStepNewmark->time() );
    }
    else if (this->is1dReducedModel())
    {
        // next time step
        M_newmark_displ_1dReduced->next( *M_disp_1dReduced );
        // up current time
        this->updateTime( M_newmark_displ_1dReduced->time() );
    }

    // update user functions which depend of time only
    this->updateUserFunctions(true);

    this->timerTool("TimeStepping").stop("updateTimeStep");

    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("SolidMechanics","updateTimeStep", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::predictorDispl()
{
    this->log("SolidMechanics","predictorDispl", "start" );

    if (isStandardModel())
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
        //this->fieldDisplacementScal1dReduced().add(M_newmark_displ_1dReduced->timeStep(),M_newmark_displ_1dReduced->currentVelocity() );
        double dt = M_newmark_displ_1dReduced->timeStep();
        if (M_newmark_displ_1dReduced->iteration() == 1)
        {
            this->fieldDisplacementScal1dReduced().add( dt, M_newmark_displ_1dReduced->previousVelocity(0) );
        }
        else
        {
            this->fieldDisplacementScal1dReduced().add( (3./2.)*dt, M_newmark_displ_1dReduced->previousVelocity(0) );
            this->fieldDisplacementScal1dReduced().add( (-1./2.)*dt, M_newmark_displ_1dReduced->previousVelocity(1) );
        }
    }

    this->log("SolidMechanics","predictorDispl", "finish" );
}

//---------------------------------------------------------------------------------------------------//


SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::solve( bool upVelAcc )
{
    this->log("SolidMechanics","solve", "start" );
    this->timerTool("Solve").start();

    this->updateBlockVectorSolution();
    if (M_pdeType=="Elasticity")
    {
        if (M_pdeSolver=="LinearSystem")
        {
            M_algebraicFactory->linearSolver( M_blockVectorSolution.vector() );
        }
        else if (M_pdeSolver == "Newton")
        {
            M_algebraicFactory->AlgoNewton2( M_blockVectorSolution.vector() );
        }
    }
    else if ( M_pdeType=="Hyper-Elasticity" || M_pdeType=="Elasticity-Large-Deformation" )
    {
        if (M_pdeSolver=="LinearSystem")
        {
            CHECK(false) <<" M_pdeType " << M_pdeType << "can not be used with LinearSystem\n";
        }
        else if (M_pdeSolver == "Newton")
        {
            M_algebraicFactory->AlgoNewton2( M_blockVectorSolution.vector() );
        }
    }
    else if (M_pdeType=="Generalised-String")
    {
        auto Uvec = this->backend1dReduced()->newVector( M_Xh_1dReduced );
        *Uvec = *M_disp_1dReduced;
        //Uvec->close();
        if (M_pdeSolver=="LinearSystem")
        {
            M_algebraicFactory_1dReduced->linearSolver(Uvec);
        }
        else if (M_pdeSolver == "Newton")
        {
            CHECK(false) <<" M_pdeType " << M_pdeType << "can not be used with Newton\n";
        }
        Uvec->close();
        *M_disp_1dReduced=*Uvec;
    }

    if (this->isStandardModel())
        M_blockVectorSolution.localize();


    if ( upVelAcc )
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

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateVelocity()
{
    this->log("SolidMechanics","updateVelocityAndAcceleration", "start" );

    if (this->isStandardModel())
    {
        M_timeStepNewmark->updateFromDisp(*M_fieldDisplacement);
    }
    else if (this->is1dReducedModel())
    {
        M_newmark_displ_1dReduced->updateFromDisp(*M_disp_1dReduced);
    }

    this->log("SolidMechanics","updateVelocityAndAcceleration", "finish" );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateNormalStressFromStruct()
{
    using namespace Feel::vf;

    if ( !M_XhNormalStress )
        this->createAdditionalFunctionSpacesNormalStress();

    auto const& u = this->fieldDisplacement();
    if ( M_pdeType=="Elasticity" )
    {
        auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
        auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

        auto const Id = eye<nDim,nDim>();
        auto eps = sym(gradv(u));//0.5*(gradv(u)+trans(gradv(u)));
        if ( !this->useDisplacementPressureFormulation() )
        {
            auto sigma = idv(coeffLame1)*trace(eps)*Id + 2*idv(coeffLame2)*eps;
            M_fieldNormalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                               _expr=sigma*vf::N() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigma = idv(p)*Id + 2*idv(coeffLame2)*eps;
            M_fieldNormalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                               _expr=sigma*vf::N() );
        }
    }
    else if ( M_pdeType=="Elasticity-Large-Deformation" )
    {
        CHECK( false ) << "TODO";
    }
    else if ( M_pdeType=="Hyper-Elasticity" )
    {
        auto const sigma = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
        if ( !this->useDisplacementPressureFormulation() )
        {
            M_fieldNormalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                               _expr=sigma*vf::N() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigmaWithPressure = Feel::vf::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()) + sigma ;
            M_fieldNormalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                               _expr=sigmaWithPressure*vf::N() );
        }
    }
}


//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateStressCriterions()
{
    this->log("SolidMechanics","updateStressCriterions", "start" );

    this->createAdditionalFunctionSpacesStressTensor();

    auto const& u = this->fieldDisplacement();
    if ( M_pdeType=="Elasticity" )
    {
        auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
        auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

        auto const Id = eye<nDim,nDim>();
        auto eps = sym(gradv(u));//0.5*(gradv(u)+trans(gradv(u)));
        if ( !this->useDisplacementPressureFormulation() )
        {
            auto sigma = idv(coeffLame1)*trace(eps)*Id + 2*idv(coeffLame2)*eps;
            M_fieldStressTensor->on(_range=elements(this->mesh()),_expr=sigma );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigma = idv(p)*Id + 2*idv(coeffLame2)*eps;
            M_fieldStressTensor->on(_range=elements(this->mesh()),_expr=sigma );
        }
    }
    else if ( M_pdeType=="Hyper-Elasticity" )
    {
        auto sigma = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
        if ( !this->useDisplacementPressureFormulation() )
        {
            M_fieldStressTensor->on(_range=elements(this->mesh()),_expr=sigma );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigmaWithPressure = Feel::vf::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()) + sigma ;
            M_fieldStressTensor->on(_range=elements(this->mesh()),_expr=sigmaWithPressure );
        }

#if 0
        auto Id = eye<nDim,nDim>();
        auto sigma_dev = sigma-(1./3.)*trace(sigma)*Id;
        M_fieldVonMisesCriterions->on(_range=elements(this->mesh()),
                                      _expr=sqrt((3./2.)*inner(sigma_dev,sigma_dev, mpl::int_<InnerProperties::IS_SAME>() )) );
#endif

    }

    typedef Eigen::Matrix<double, nDim, nDim> matrixN_type;
    Eigen::EigenSolver< matrixN_type > eigenSolver;
    matrixN_type sigma_eigen_matrix;

    auto dof = M_XhStressTensor->dof();
    for ( auto const& elt : elements(this->mesh()) )
    {
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

            size_type gdofScal = M_XhStressTensor->compSpace()->dof()->localToGlobal( elt, j, 0 ).index();

            double resPrincipalStress1=0;
            for ( uint16_type comp=0; comp < space_stress_tensor_type::nComponents1;++comp )
                resPrincipalStress1 += real(eigenSolver.eigenvalues()[comp]);
            M_fieldsPrincipalStresses[0]->set( gdofScal, resPrincipalStress1 );
            double resPrincipalStress2=0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for ( uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    resPrincipalStress2 += real(eigenSolver.eigenvalues()[comp1])*real(eigenSolver.eigenvalues()[comp2]);
            M_fieldsPrincipalStresses[1]->set( gdofScal, resPrincipalStress2 );
            if ( space_stress_tensor_type::nComponents1 == 3 )
            {
                double resPrincipalStress3 = real(eigenSolver.eigenvalues()[0])*real(eigenSolver.eigenvalues()[1])*real(eigenSolver.eigenvalues()[2]);
                M_fieldsPrincipalStresses[2]->set( gdofScal, resPrincipalStress3 );
            }

            double resTresca = 0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for (uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    resTresca =  std::max( resTresca, std::abs( real(eigenSolver.eigenvalues()[comp1])- real(eigenSolver.eigenvalues()[comp2]) ) );

            double resVonMises = 0;
            for ( uint16_type comp1=0; comp1 < space_stress_tensor_type::nComponents1;++comp1 )
                for (uint16_type comp2=comp1+1; comp2 < space_stress_tensor_type::nComponents2;++comp2 )
                    resVonMises += (1./2.)*std::pow( real(eigenSolver.eigenvalues()[comp1])- real(eigenSolver.eigenvalues()[comp2]), 2 );
            resVonMises = std::sqrt( resVonMises );

            M_fieldTrescaCriterions->set( gdofScal, resTresca );
            M_fieldVonMisesCriterions->set( gdofScal, resVonMises );
        }
    }

    this->log("SolidMechanics","updateStressCriterions", "finish" );

}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::block_pattern_type
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::blockPattern() const
{
    if (M_useDisplacementPressureFormulation)
        return BlocksStencilPattern(2,2) << size_type(Pattern::COUPLED) << size_type(Pattern::COUPLED)
                                         << size_type(Pattern::COUPLED) << size_type(Pattern::ZERO);
    else
        return BlocksStencilPattern(1,1) << size_type(Pattern::COUPLED);
}

//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateStressTensorBis( element_stress_ptrtype stressN)
{

    if (this->verbose()) std::cout << "[SolidMechanics] : updateStressTensor start\n";

    boost::timer btime; btime.restart();

    // maybe just pointer??
    *M_normalStressFromFluid = *stressN;
#if 0
    if (M_is1dReduced)
    {
#if 0
        auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
        ForEachBC( bcDef,cl::paroi_mobile,
                   TransfertStress2dTo1d(PhysicalName); )
#else
            this->transfertNormalStress2dTo1dWithInterpolation();
#endif
    }
#endif
    if (this->verbose()) std::cout << "[SolidMechanics] : updateStressTensor finish in "<<btime.elapsed()<<"\n";

}
#endif
//---------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMaxDisp() const
{
    double res = 0;
    if(this->isStandardModel())
    {
        auto u = this->fieldDisplacement();
        auto dispMagnOnBoundary = vf::project(_space=functionSpaceDisplacement()->compSpace(),
                                              _range=elements(this->mesh()),
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

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeMaxDispOnBoundary(std::string __marker) const
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
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeExtremumValue( std::string const& field, std::list<std::string> const& markers, std::string const& type ) const
{
    if(!this->isStandardModel()) return 0.;
    CHECK( type == "max" || type == "min" ) << "invalid type " << type;

    if ( field == "displacement" || field == "velocity" || field == "acceleration" )
    {
        auto fieldMagnitude = functionSpaceDisplacement()->compSpace()->elementPtr();
        if ( markers.empty() )
        {
            if ( field == "displacement" )
                fieldMagnitude->on(_range=elements(this->mesh()), _expr=norm2( idv(this->fieldDisplacement()) ) );
            else if ( field == "velocity" )
                fieldMagnitude->on(_range=elements(this->mesh()), _expr=norm2( idv(this->fieldVelocity()) ) );
            else if ( field == "acceleration" )
                fieldMagnitude->on(_range=elements(this->mesh()), _expr=norm2( idv(this->fieldAcceleration()) ) );
        }
        else
        {
            if ( this->mesh()->hasFaceMarker( markers.front() ) )
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
        CHECK( M_useDisplacementPressureFormulation ) << "model does not take into account the pressure";
        CHECK( false ) << "TODO pressure max/min";
    }

    //if(this->isStandardModel()) res = this->fieldDisplacement().max();
    //else if(this->is1dReducedModel()) res = this->fieldDisplacementScal1dReduced().max();

    return 0.;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeVolumeVariation() const
{
    //using namespace Feel::vf;
    if(!this->isStandardModel()) return 0.;

    auto const Id = eye<nDim,nDim>();
    auto const& u = this->fieldDisplacement();
    auto Fv = Id + gradv(u);
    auto detFv = det(Fv);
    double newArea = integrate(_range=elements(u.mesh()),
                             _expr=detFv,
                             _geomap=this->geomap() ).evaluate()(0,0);
    double refAera = integrate(_range=elements(u.mesh()),
                               _expr=cst(1.),
                               _geomap=this->geomap() ).evaluate()(0,0);

    return (newArea-refAera)/refAera;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceDispFrom1dDisp()
{
    M_disp_vect_1dReduced->on( _range=elements(M_mesh_1dReduced),
                               _expr=idv(M_disp_1dReduced)*oneY() );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceScalStressDispFromVectStress()
{
    M_stress_1dReduced->on( _range=elements(M_mesh_1dReduced),
                            _expr=-inner(idv(M_stress_vect_1dReduced),oneY()) );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceVelocityFrom1dVelocity()
{
    M_velocity_vect_1dReduced->on( _range=elements(M_mesh_1dReduced),
                                   _expr=idv(M_newmark_displ_1dReduced->currentVelocity())*oneY() );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
typename SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::element_vect_1dreduced_ptrtype
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::extendVelocity1dReducedVectorial( element_1dreduced_type const& vel1d ) const
{
    auto res = M_Xh_vect_1dReduced->elementPtr( vf::idv(vel1d)*vf::oneY() );
    return res;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateSubMeshDispFSIFromPrevious()
{
    auto subfsimesh = M_fieldSubMeshDispFSI->mesh();
    M_fieldSubMeshDispFSI->on(_range=elements(subfsimesh),
                              _expr=idv(this->timeStepNewmark()->previousUnknown()) );
    M_meshMoverTrace.apply( subfsimesh,*M_fieldSubMeshDispFSI );
}


} // FeelModels

} // Feel




