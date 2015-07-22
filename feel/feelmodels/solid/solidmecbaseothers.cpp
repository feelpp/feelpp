/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */


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
    if (this->doRestart() && this->restartPath().empty())
    {
        if (!M_isHOVisu)
        {
            if ( M_exporter->doExport() ) M_exporter->restart(time);
        }
        else
        {
#if defined(FEELPP_HAS_VTK)
            if ( M_exporter_ho->doExport() ) M_exporter_ho->restart(time);
#endif
        }
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::restartExporters1dReduced( double time )
{
    if (this->doRestart() && this->restartPath().empty())
    {
        if ( M_exporter_1dReduced->doExport() ) M_exporter_1dReduced->restart(time);
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
    if (isStandardModel()) {nElt=M_mesh->numGlobalElements(); nDof=M_Xh->nDof();}
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
#if defined(FEELPP_HAS_VTK)
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
    if (M_doExportDisplacement) doExport_str=(doExport_str.empty())?"displacement":doExport_str+" - displacement";
    if (M_doExportVelocity) doExport_str=(doExport_str.empty())?"velocity":doExport_str+" - velocity";
    if (M_doExportAcceleration) doExport_str=(doExport_str.empty())?"acceleration":doExport_str+" - acceleration";
    if (M_doExportPressure) doExport_str=(doExport_str.empty())?"pressure":doExport_str+" - pressure";
    if (M_doExportNormalStress) doExport_str=(doExport_str.empty())?"normal stress":doExport_str+" - normal stress";
    if (M_doExportVelocityInterfaceFromFluid) doExport_str=(doExport_str.empty())?"velocty interface from fluid":doExport_str+" - velocty interface from fluid";

    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n||==============================================||"
           << "\n||----------Info : SolidMechanics---------------||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Appli Repository : " << this->appliRepository()
           << "\n   Physical Model"
           << "\n     -- pde name : " << M_pdeType
           << "\n     -- material law : " << this->mechanicalProperties()->materialLaw()
           << "\n     -- use displacement-pressure formulation : " << std::boolalpha << M_useDisplacementPressureFormulation
           << "\n     -- time mode : " << StateTemporal
           << "\n   Physical Parameters"
           << "\n     -- rho : " << this->mechanicalProperties()->cstRho()
           << "\n     -- young modulus : " << this->mechanicalProperties()->cstYoungModulus()
           << "\n     -- coeff poisson : " << this->mechanicalProperties()->cstCoeffPoisson()
           << "\n     -- coeff Lame 1 : " << this->mechanicalProperties()->cstCoeffLame1()
           << "\n     -- coeff Lame 2 : " << this->mechanicalProperties()->cstCoeffLame2()
           << "\n   Boundary conditions";
    *_ostr << this->getInfoDirichletBC()
           << this->getInfoNeumannBC()
           << this->getInfoRobinBC()
           << this->getInfoFluidStructureInterfaceBC();
    *_ostr << "\n   Space Discretization"
           << "\n     -- msh file name   : " << this->mshfileStr()
           << "\n     -- nb elt in mesh : " << nElt
           << "\n     -- nb dof (displacement) : " << nDof
           << "\n     -- polynomial order : " << nOrder;
    if ( M_useDisplacementPressureFormulation )
        *_ostr << "\n     -- nb dof (pressure) : " << M_XhPressure->nDof();
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
        myblockGraph(indexBlock,0) = stencil(_test=M_XhPressure,_trial=M_Xh,
                                             _diag_is_nonzero=false,_close=false)->graph();
        myblockGraph(0,indexBlock) = stencil(_test=M_Xh,_trial=M_XhPressure,
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
        auto & vecAlgebraic = M_blockVectorSolution.vector();
        auto const& fieldDisp = this->fieldDisplacement();
        for (int k=0;k< M_Xh->nLocalDofWithGhost() ;++k)
            vecAlgebraic->set(k, fieldDisp(k) );

        if ( M_useDisplacementPressureFormulation )
        {
            auto const& fieldPressure = this->fieldPressure();
            size_type startDofIndexPressure = this->startDofIndexFieldsInMatrix().find("pressure")->second;
            for (int k=0;k< M_XhPressure->nLocalDofWithGhost() ;++k)
                vecAlgebraic->set( startDofIndexPressure+k, fieldPressure(k) );
        }
        vecAlgebraic->close();
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

    this->log("SolidMechanics","exportResults", "finish" );

} // SolidMechanics::export
SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResultsImpl( double time )
{
    if (this->isStandardModel())
    {
        if ( !M_exporter->doExport() ) return;
        //M_exporter->step( time )->setMesh( M_mesh );
        bool hasFieldToExport = false;
        if ( M_doExportDisplacement )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"displacement"), this->fieldDisplacement() );
            hasFieldToExport = true;
        }
        if ( M_useDisplacementPressureFormulation && M_doExportPressure )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"pressure"), *M_fieldPressure );
            hasFieldToExport = true;
        }
        if ( M_doExportVelocity)
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity"), this->timeStepNewmark()->currentVelocity() );
            hasFieldToExport = true;
        }
        if ( M_doExportAcceleration )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"acceleration"), this->timeStepNewmark()->currentAcceleration() );
            hasFieldToExport = true;
        }
        if ( M_doExportNormalStress )
        {
            this->updateNormalStressFromStruct();
            M_exporter->step( time )->add( prefixvm(this->prefix(),"normalstress"), *M_normalStressFromStruct/*M_normalStressFromFluid*/ );
            hasFieldToExport = true;
        }
        if ( M_doExportVelocityInterfaceFromFluid && this->velocityInterfaceFromFluid() )
        {
            M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity-interface-from-fluid"), *this->velocityInterfaceFromFluid() );
            hasFieldToExport = true;
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
        if ( M_doExportDisplacement )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-displacement"), *M_disp_1dReduced );
            hasFieldToExport = true;
        }
        if ( M_doExportVelocity )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-velocity"), this->timeStepNewmark1dReduced()->currentVelocity()  /**M_velocity_1dReduced*/ );
            hasFieldToExport = true;
        }
        if ( M_doExportAcceleration )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-acceleration"), this->timeStepNewmark1dReduced()->currentAcceleration() );
            hasFieldToExport = true;
        }
        if ( M_doExportNormalStress )
        {
            M_exporter_1dReduced->step( time )->add( prefixvm(this->prefix(),"1d-reduced-stress"),*M_stress_1dReduced );
            hasFieldToExport = true;
        }
        //M_exporter_1dReduced->step( time )->add( "structure-1d-reduced-stressvec",*M_stress_vect_1dReduced);
        if ( hasFieldToExport )
            M_exporter_1dReduced->save();
    }
}

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::exportResultsImplHO( double time )
{
    if (this->isStandardModel())
    {
#if defined(FEELPP_HAS_VTK)
        if ( !M_exporter_ho->doExport() ) return;
        //M_exporter_ho->step( time )->setMesh( M_displacementVisuHO->mesh() );

        bool hasFieldToExport = false;
        if ( M_doExportDisplacement )
        {
            M_opIdisplacement->apply(this->fieldDisplacement(),*M_displacementVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"displacement-ho"), *M_displacementVisuHO );
            hasFieldToExport = true;
        }
        if ( M_useDisplacementPressureFormulation && M_doExportPressure )
        {
            M_opIpressure->apply(*M_fieldPressure,*M_pressureVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"pressure-ho"), *M_pressureVisuHO );
            hasFieldToExport = true;
        }
        if ( M_doExportVelocity )
        {
            auto velocityVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentVelocity(),velocityVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"velocity-ho"), velocityVisuHO );
            hasFieldToExport = true;
        }
        if ( M_doExportAcceleration )
        {
            auto accelerationVisuHO = M_XhVectorialVisuHO->element();
            M_opIdisplacement->apply( this->timeStepNewmark()->currentAcceleration(),accelerationVisuHO );
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"acceleration-ho"), accelerationVisuHO );
            hasFieldToExport = true;
        }
        if ( M_doExportNormalStress )
        {
            auto normalstressVisuHO = M_XhVectorialVisuHO->element();
            this->updateNormalStressFromStruct();
            M_opInormalstress->apply( *M_normalStressFromStruct/*M_normalStressFromFluid*/,normalstressVisuHO);
            M_exporter_ho->step( time )->add( prefixvm(this->prefix(),"normalstress-ho"), normalstressVisuHO );
            hasFieldToExport = true;
        }
        if ( hasFieldToExport )
            M_exporter_ho->save();
#endif
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
        M_newmark_displ_struct->next( *M_fieldDisplacement );
        if ( M_useDisplacementPressureFormulation ) M_savetsPressure->next(*M_fieldPressure);
        // up current time
        this->updateTime( M_newmark_displ_struct->time() );
    }
    else if (this->is1dReducedModel())
    {
        // next time step
        M_newmark_displ_1dReduced->next( *M_disp_1dReduced );
        // up current time
        this->updateTime( M_newmark_displ_1dReduced->time() );
    }

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
        double dt = M_newmark_displ_struct->timeStep();
        if (M_newmark_displ_struct->iteration() == 1) //order 1: \tilde{u} = u_n + dt* \dot{u}_n
        {
            //this->fieldDisplacement().add(M_newmark_displ_struct->timeStep(),M_newmark_displ_struct->currentVelocity());
            this->fieldDisplacement().add( dt, M_newmark_displ_struct->previousVelocity(0));
        }
        else //order 2: \tilde{u} = u_n + dt*( (3/2)*\dot{u}_n-(1/2)*\dot{u}_{n-1}
        {
            //this->fieldDisplacement().add((3./2.)*M_newmark_displ_struct->timeStep(), M_newmark_displ_struct->currentVelocity() );
            //this->fieldDisplacement().add((-1./2.)*M_newmark_displ_struct->timeStep(), M_newmark_displ_struct->previousVelocity());
            this->fieldDisplacement().add( (3./2.)*dt, M_newmark_displ_struct->previousVelocity(0) );
            this->fieldDisplacement().add( (-1./2.)*dt, M_newmark_displ_struct->previousVelocity(1) );
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
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::saveDataInPoint(const std::list<boost::tuple<std::string,typename mesh_type::node_type> > & __listPt, bool extrapolate)
{
    auto it = __listPt.begin();
    auto en = __listPt.end();

    for ( ; it!=en ; ++it)
    {
        auto nameFile = prefixvm(this->prefix(),"disp-")+boost::get<0>(*it)+".data";
        if (M_nameFilesData.find(nameFile)==M_nameFilesData.end())
        {
            if (!this->doRestart() && this->worldComm().isMasterRank())
            {
                //ATTENTION EFFACER LES DONNEES REDONTANTE EN CAS DE RESTART
                std::ofstream newfileDisp(nameFile.c_str(), std::ios::out | std::ios::trunc);
                newfileDisp.close();
                M_nameFilesData.insert(nameFile);
            }
        }
#if 1
        auto const dispOnPt = this->fieldDisplacement()(boost::get<1>(*it),extrapolate);
        if (this->worldComm().isMasterRank())
        {
            std::ofstream fileDisp(nameFile.c_str(), std::ios::out | std::ios::app);
            fileDisp << this->time();
            for (uint16_type i = 0; i<mesh_type::nRealDim;++i)
                fileDisp << " " << dispOnPt(i,0,0);
            //fileDisp << " " << this->fieldDisplacement()(boost::get<1>(*it))(i,0,0);
            fileDisp << "\n";
            fileDisp.close();
        }
#else
        auto ctx = this->fieldDisplacement()->functionSpace()->context();
        ctx.add( it->get<1>() );
        auto eval_at_node = evaluateFromContext ( ctx , _expr = idv (u) );
        double value1 = eval_at_node ( it->get<1>() );
#endif
    }
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateVelocity()
{
    this->log("SolidMechanics","updateVelocityAndAcceleration", "start" );

    if (this->isStandardModel())
    {
        M_newmark_displ_struct->updateFromDisp(*M_fieldDisplacement);
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

    if ( !M_XhStress )
        M_XhStress = space_stress_type::New( _mesh=M_mesh, _worldscomm=this->localNonCompositeWorldsComm() );
    if ( !M_normalStressFromStruct )
        M_normalStressFromStruct.reset(new element_stress_type( M_XhStress, "normalStressBoundaryFromStruct" ));

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
            M_normalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                          _expr=sigma*vf::N() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigma = idv(p)*Id + 2*idv(coeffLame2)*eps;
            M_normalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
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
            M_normalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                          _expr=sigma*vf::N() );
        }
        else
        {
            auto const& p = this->fieldPressure();
            auto sigmaWithPressure = Feel::vf::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()) + sigma ;
            M_normalStressFromStruct->on( _range=boundaryfaces(this->mesh()),
                                          _expr=sigmaWithPressure*vf::N() );
        }
    }
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

//---------------------------------------------------------------------------------------------------//

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

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
double
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::computeIncompressibility() const
{
    using namespace Feel::vf;

    auto const Id = eye<nDim,nDim>();

    auto u = this->fieldDisplacement();
    auto Fv = Id + gradv(u);
    auto detFv = det(Fv);
    auto res = integrate(_range=elements(u.mesh()),
                         _expr=detFv,
                         _geomap=this->geomap() ).evaluate()(0,0);
    auto aera = integrate(_range=elements(u.mesh()),
                          _expr=cst(1.),
                          _geomap=this->geomap() ).evaluate()(0,0);

    return res/aera;
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceDispFrom1dDisp()
{
    *M_disp_vect_1dReduced = vf::project( _space=M_Xh_vect_1dReduced,
                                           _range=elements(M_mesh_1dReduced),
                                           _expr=vf::idv(M_disp_1dReduced)*vf::oneY());
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceScalStressDispFromVectStress()
{
    *M_stress_1dReduced = vf::project( _space=M_stress_1dReduced->functionSpace(),
                                        _range=elements(M_mesh_1dReduced),
                                        _expr=-vf::trans(vf::idv(M_stress_vect_1dReduced))*vf::oneY() );
}

//---------------------------------------------------------------------------------------------------//

SOLIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateInterfaceVelocityFrom1dVelocity()
{
    *M_velocity_vect_1dReduced = vf::project(_space=M_Xh_vect_1dReduced,
                                              _range=elements(M_mesh_1dReduced),
                                              _expr=vf::idv(M_newmark_displ_1dReduced->currentVelocity())*vf::oneY());
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




