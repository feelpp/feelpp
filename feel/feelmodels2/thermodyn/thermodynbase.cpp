/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/thermodyn/thermodynbase.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/geotool.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels2/modelmesh/reloadmesh.hpp>

namespace Feel
{
namespace FeelModels
{

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::ThermoDynamicsBase( bool __isStationary,
                                                            std::string __prefix,
                                                            WorldComm const& __worldComm,
                                                            bool __buildMesh,
                                                            std::string __subPrefix,
                                                            std::string __appliShortRepository )
    :
    super_type( __isStationary,__prefix,__worldComm,__subPrefix,__appliShortRepository)
{}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::build()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","build", "start",
                                               this->worldComm(),this->verboseAllProc());

    boost::mpi::timer mpiTimer;
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    this->createMesh();
    double timeMesh = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    double timeFunctionSpaces = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // bdf time schema
    this->createTimeDiscretisation();
    double timeDiscretisation = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // ALE mode (maybe)
    //this->createALE();
    //double timeALE = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // physical parameters
    //this->createOthers();
    //double timeOthers = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    //export
    this->createExporters();
    double timeExporters = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // save timers
    //if ( this->scalabilitySave() ) this->saveTimerBuild(timeMesh,timeFunctionSpaces,timeDiscretisation,
    //                                                    timeOthers,timeALE,timeExporters);
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","build", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::loadMesh( mesh_ptrtype mesh )
{
    boost::mpi::timer mpiTimer;
    //-----------------------------------------------------------------------------//
    // create or reload mesh
    if (this->doRestart()) this->createMesh();
    else this->M_mesh = mesh;
    double timeMesh = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // functionSpaces and elements
    this->createFunctionSpaces();
    double timeFunctionSpaces = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // bdf time schema
    this->createTimeDiscretisation();
    double timeDiscretisation = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // ALE mode (maybe)
    //this->createALE();
    //double timeALE = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // physical parameters
    //this->createOthers();
    //double timeOthers = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    //export
    this->createExporters();
    double timeExporters = mpiTimer.elapsed();mpiTimer.restart();
    //-----------------------------------------------------------------------------//
    // save timers
    //if ( this->scalabilitySave() ) this->saveTimerBuild(timeMesh,timeFunctionSpaces,timeDiscretisation,
    //                                                    timeOthers,timeALE,timeExporters);
    //-----------------------------------------------------------------------------//

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","build", "finish",
                                               this->worldComm(),this->verboseAllProc());

}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::loadParameterFromOptionsVm()
{
    M_meshSize = doption(_name="hsize",_prefix=this->prefix());
    M_thermalConductivity = doption(_name="thermal-conductivity",_prefix=this->prefix()); // [ W/(m*K) ]
    M_rho = doption(_name="rho",_prefix=this->prefix()); // density [ kg/(m^3) ]
    M_heatCapacity = doption(_name="heat-capacity",_prefix=this->prefix()); // [ J/(kg*K) ]

    M_fieldVelocityConvectionIsUsed = boption(_name="use_velocity-convection",_prefix=this->prefix());
    M_fieldVelocityConvectionIsIncompressible = boption(_name="velocity-convection_is_incompressible",_prefix=this->prefix());

    M_doExportAll = boption(_name="do_export_all",_prefix=this->prefix());
    M_doExportVelocityConvection = boption(_name="do_export_velocity-convection",_prefix=this->prefix());
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createMesh()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh", "start",
                                               this->worldComm(),this->verboseAllProc());
    boost::timer thetimer;

    // save path of file mesh
    auto fmpath = this->fileNameMeshPath();
    if (this->doRestart())
    {
        if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh",
                                                   "restart with : "+fmpath,
                                                   this->worldComm(),this->verboseAllProc());

        if ( !this->restartPath().empty() )
        {
            this->M_mesh = reloadMesh<mesh_type>(this->restartPath()+"/"+fmpath,this->worldComm());
        }
        else
        {
            this->M_mesh = reloadMesh<mesh_type>(fmpath,this->worldComm());
        }
    }
    else
    {
        if (this->hasMshfileStr())
        {
            std::string path = this->appliRepository();
            std::string mshfileRebuildPartitions = path + "/" + this->prefix() + ".msh";

            if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh", "load msh file : " + this->mshfileStr(),
                                                       this->worldComm(),this->verboseAllProc());

            this->M_mesh = loadGMSHMesh(_mesh=new mesh_type,
                                        _filename=this->mshfileStr(),
                                        _worldcomm=this->worldComm(),
                                        _rebuild_partitions=this->rebuildMeshPartitions(),
                                        _rebuild_partitions_filename=mshfileRebuildPartitions,
                                        _partitions=this->worldComm().localSize(),
                                        _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK);

            if (this->rebuildMeshPartitions()) this->setMshfileStr(mshfileRebuildPartitions);
        }
        else if (this->hasGeofileStr())
        {
            std::string path = this->appliRepository();
            std::string mshfile = path + "/" + this->prefix() + ".msh";
            this->setMshfileStr(mshfile);
            this->M_mesh = GeoTool::createMeshFromGeoFile<mesh_type>(this->geofileStr(),this->prefix(),this->meshSize(),1,
                                                                     this->worldComm().localSize(),this->worldComm());
        }
        else
        {
            std::string geotoolSavePath;

            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh", "change rep -> "+ this->geotoolSaveDirectory(),
                                                           this->worldComm(),this->verboseAllProc());
                Environment::changeRepository( _directory=boost::format(this->geotoolSaveDirectory()), _subdir=false );

                geotoolSavePath = Environment::rootRepository()+"/"+ this->geotoolSaveDirectory();
            }
            else
            {
                geotoolSavePath = this->appliRepository();
            }
            std::string geotoolSaveName = this->geotoolSaveName();
            std::string geofilename = geotoolSavePath + "/" + geotoolSaveName;// without .geo
            std::string mshfilename = geotoolSavePath + "/" + geotoolSaveName + ".msh";
            this->setMshfileStr(mshfilename);

            if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh", "build mesh by using geotool desc",
                                                       this->worldComm(),this->verboseAllProc());

            this->loadConfigMeshFile(geofilename);

            if ( this->geotoolSaveDirectory()!=this->appliShortRepository() )
            {
                if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createMesh", "change rep -> " + this->appliRepository() ,
                                                           this->worldComm(),this->verboseAllProc());
                Environment::changeRepository( _directory=boost::format(this->appliShortRepository()), _subdir=true );
            }

        }
        this->saveMSHfilePath(fmpath);
    }

    double tElpased = thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log( this->prefix()+".ThermoDynamics","createMesh",
                                                (boost::format("finish in %1% s")%tElpased).str(),
                                                this->worldComm(),this->verboseAllProc() );

} // createMesh()

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createFunctionSpaces()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createFunctionSpaces", "start",
                                               this->worldComm(),this->verboseAllProc());

    // functionspace
    M_Xh = space_temperature_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );

    M_fieldTemperature.reset( new element_temperature_type(M_Xh,"U"));

    if ( this->fieldVelocityConvectionIsUsed() )
        this->updateForUseFunctionSpacesVelocityConvection();

    // backend : use worldComm of Xh
    M_backend = backend_type::build( soption( _name="backend" ), this->prefix(), M_Xh->worldComm() );

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createFunctionSpaces", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateForUseFunctionSpacesVelocityConvection()
{
    if ( !M_XhVelocityConvection )
        M_XhVelocityConvection = space_velocityconvection_type::New( _mesh=M_mesh, _worldscomm=this->worldsComm() );
    if ( !M_fieldVelocityConvection )
        M_fieldVelocityConvection.reset( new element_velocityconvection_type(M_XhVelocityConvection,"VelocityConvection"));
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createTimeDiscretisation()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createTimeDiscretisation", "start",
                                               this->worldComm(),this->verboseAllProc());

    std::string suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdfTemperature = bdf( _vm=Environment::vm(), _space=this->spaceTemperature(),
                            _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature"+suffixName)),
                            _prefix=this->prefix(),
                            // don't use the fluid.bdf {initial,final,step}time but the general bdf info, the order will be from fluid.bdf
                            _initial_time=this->timeInitial(),
                            _final_time=this->timeFinal(),
                            _time_step=this->timeStep(),
                            _restart=this->doRestart(),
                            _restart_path=this->restartPath(),
                            _restart_at_last_save=this->restartAtLastSave(),
                            _save=this->bdfSaveInFile(), _freq=this->bdfSaveFreq() );

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createTimeDiscretisation", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::createExporters()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createExporters", "start",
                                               this->worldComm(),this->verboseAllProc());

    std::string geoExportType="static";//change_coords_only, change, static
    M_exporter = exporter( _mesh=this->mesh(),
                           _name="Export",
                           _geo=geoExportType,
                           _path=this->exporterPath() );

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","createExporters", "finish",
                                               this->worldComm(),this->verboseAllProc());
}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
BlocksBaseGraphCSR
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::buildBlockMatrixGraph() const
{
    int nBlock = this->nBlockMatrixGraph();
    BlocksBaseGraphCSR myblockGraph(nBlock,nBlock);
    myblockGraph(0,0) = stencil(_test=this->spaceTemperature(),
                                _trial=this->spaceTemperature() )->graph();
    return myblockGraph;
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
size_type
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::nLocalDof() const
{
    size_type res = this->spaceTemperature()->nLocalDofWithGhost();
    return res;
}

THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::init( bool buildMethodNum, model_algebraic_factory_type::appli_ptrtype const& app )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","init", "start",
                                               this->worldComm(),this->verboseAllProc());

    if ( this->fieldVelocityConvectionIsUsed() )
        this->updateForUseFunctionSpacesVelocityConvection();

    //-------------------------------------------------//
    // vector solution
    M_blockVectorSolution.resize( 1 );
    M_blockVectorSolution(0) = this->fieldTemperature();
    // init petsc vector associated to the block
    M_blockVectorSolution.buildVector( this->backend() );
#if 0
    if ( buildMethodNum )
    {
        // matrix graph of non zero
        //auto graph = stencil(_test=this->spaceTemperature(),
        //                     _trial=this->spaceTemperature() )->graph();
        auto graph = this->buildBlockMatrixGraph()(0,0);
        // tool for assembly and solver
        M_methodNum.reset( new methodsnum_type(this->shared_from_this(),this->backend(),
                                               graph, graph->mapRow().indexSplit() ) );
    }
#endif
    //-------------------------------------------------//
    // start or restart time step scheme and exporter
    if (!this->doRestart())
    {
        // start time step
        M_bdfTemperature->start(*this->fieldTemperature());
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
    else
    {
        // start time step
        M_bdfTemperature->restart();
        // load a previous solution as current solution
        *this->fieldTemperature() = M_bdfTemperature->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdfTemperature->timeInitial() );
        // restart exporter
        this->restartExporters();
        // up current time
        this->updateTime( M_bdfTemperature->time() );
    }
    //-------------------------------------------------//
    if ( buildMethodNum )
    {
        // matrix graph of non zero
        auto graph = this->buildBlockMatrixGraph()(0,0);
        // tool for assembly and solver
        M_methodNum.reset( new model_algebraic_factory_type(app,this->backend(),
                                                            graph, graph->mapRow().indexSplit() ) );
    }

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","init", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::restartExporters()
{
    if (this->doRestart() && this->restartPath().empty() )
    {
        if ( M_exporter->doExport() )
            M_exporter->restart(this->timeInitial());
    }
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
boost::shared_ptr<std::ostringstream>
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||----------Info : ThermoDynamics---------------||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n   Prefix : " << this->prefix()
           << "\n   Appli Repository : " << this->appliRepository();
    *_ostr << "\n   Physical Model"
           << "\n     -- time mode           : " << std::string( (this->isStationary())?"Stationary":"Transient")
           << "\n     -- velocity-convection : " << std::string( (this->fieldVelocityConvectionIsUsedAndOperational())?"Yes":"No" );
    *_ostr << "\n   Boundary conditions"
           << this->getInfoDirichletBC()
           << this->getInfoNeumannBC();
    *_ostr << "\n   Physical Parameters"
           << "\n     -- thermal conductivity : " << this->thermalConductivity()
           << "\n     -- heat capacity        : " << this->heatCapacity()
           << "\n     -- density              : "  << this->rho();
    *_ostr << "\n   Mesh Discretization"
           << "\n     -- msh filename      : " << this->mshfileStr()
           << "\n     -- number of element : " << M_mesh->numGlobalElements()
           << "\n     -- order             : " << nOrderGeo;
    *_ostr << "\n   Space Temperature Discretization"
           << "\n     -- order         : " << nOrderPoly
           << "\n     -- number of dof : " << M_Xh->nDof() << " (" << M_Xh->nLocalDof() << ")";
    if ( this->fieldVelocityConvectionIsUsedAndOperational() )
        *_ostr << "\n   Space Velocity Convection Discretization"
               << "\n     -- number of dof : " << M_XhVelocityConvection->nDof() << " (" << M_XhVelocityConvection->nLocalDof() << ")";
    *_ostr << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n||==============================================||"
           << "\n";

    return _ostr;
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::solve()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","solve", "start",
                                               this->worldComm(),this->verboseAllProc());

    M_methodNum->linearSolver(this->blockVectorSolution().vector());
    M_blockVectorSolution.localize();

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","solve", "finish",
                                               this->worldComm(),this->verboseAllProc());

}



THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::exportResults( double time )
{
    if ( !M_exporter->doExport() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","exportResults", "start",
                                               this->worldComm(),this->verboseAllProc());

    M_exporter->step( time )->add( prefixvm(this->prefix(),"temperature"),
                                   prefixvm(this->prefix(),prefixvm(this->subPrefix(),"temperature")),
                                   *this->fieldTemperature() );
    if ( ( M_doExportVelocityConvection || M_doExportAll ) && this->fieldVelocityConvectionIsOperational() )
    {
        M_exporter->step( time )->add( prefixvm(this->prefix(),"velocity-convection"),
                                       prefixvm(this->prefix(),prefixvm(this->subPrefix(),"velocity-convection")),
                                       *this->fieldVelocityConvection() );
    }
    M_exporter->save();

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","exportResults", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateBdf()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBdf", "start",
                                               this->worldComm(),this->verboseAllProc());
    int previousTimeOrder = this->timeStepBdfTemperature()->timeOrder();

    M_bdfTemperature->next( *this->fieldTemperature() );

    int currentTimeOrder = this->timeStepBdfTemperature()->timeOrder();

    this->updateTime( this->timeStepBdfTemperature()->time() );

    // maybe rebuild cst jacobian or linear
    if ( M_methodNum &&
         previousTimeOrder!=currentTimeOrder &&
         this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        if (!this->rebuildCstPartInLinearSystem())
        {
            if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBdf", "do rebuildCstLinearPDE",
                                                       this->worldComm(),this->verboseAllProc());
            M_methodNum->rebuildCstLinearPDE(this->blockVectorSolution().vector());
        }
    }
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBdf", "finish",
                                               this->worldComm(),this->verboseAllProc());

}









THERMODYNAMICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
THERMODYNAMICSBASE_CLASS_TEMPLATE_TYPE::updateLinearPDE( const vector_ptrtype& X, sparse_matrix_ptrtype& A, vector_ptrtype& F, bool buildCstPart,
                                                         sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                                         bool _doClose, bool _doBCStrongDirichlet ) const
{
    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateLinearPDE", "start"+sc,
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer thetimer;

    //--------------------------------------------------------------------------------------------------//

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();

    auto const& u = *this->fieldTemperature();
    auto const& v = *this->fieldTemperature();

    //boundaries conditions
    //auto const bcDef = THERMODYNAMICS_BC(this->shared_from_this());

    //--------------------------------------------------------------------------------------------------//

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=rowStartInVector );


    //--------------------------------------------------------------------------------------------------//

    if ( buildCstPart )
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= this->thermalConductivity()*inner(gradt(u),grad(v)),
                       _geomap=this->geomap() );
    }

    if ( this->fieldVelocityConvectionIsUsedAndOperational() && !buildCstPart )
    {
        double thecoeff = this->rho()*this->heatCapacity();
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= thecoeff*(gradt(u)*idv(this->fieldVelocityConvection()))*id(v),
                       _geomap=this->geomap() );
        if ( !this->fieldVelocityConvectionIsIncompressible() )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= thecoeff*(idt(u)*divv(this->fieldVelocityConvection()))*id(v),
                           _geomap=this->geomap() );
        }
    }

    this->updateSourceTermLinearPDE(F, buildCstPart);

    if (!this->isStationary())
    {
        bool buildNonCstPart=!buildCstPart;
        bool BuildNonCstPart_Form2TransientTerm = buildNonCstPart;
        bool BuildNonCstPart_Form1TransientTerm = buildNonCstPart;
        if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
        {
            BuildNonCstPart_Form2TransientTerm = buildCstPart;
        }

        if (BuildNonCstPart_Form2TransientTerm)
        {
            double thecoeff = this->rho()*this->heatCapacity()*this->timeStepBdfTemperature()->polyDerivCoefficient(0);
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= thecoeff*idt(u)*id(v),
                           _geomap=this->geomap() );
        }

        if (BuildNonCstPart_Form1TransientTerm)
        {
            double thecoeff = this->rho()*this->heatCapacity();
            auto rhsTimeStep = this->timeStepBdfTemperature()->polyDeriv();
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= thecoeff*idv(rhsTimeStep)*id(v),
                           _geomap=this->geomap() );
        }
    }

    // update neumann bc
    this->updateWeakBCLinearPDE(A,F,buildCstPart);

    if (/*this->hasDirichletBCelimination()  &&*/ !buildCstPart && _doBCStrongDirichlet)
    {
        this->updateBCStrongDirichletLinearPDE(A,F);
    }


    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateLinearPDE",
                                               "finish in "+(boost::format("%1% s") % timeElapsed).str(),
                                               this->worldComm(),this->verboseAllProc());
}


} // end namespace FeelModels
} // end namespace Feel
