#include <feel/feelmodels/hdg/mixedelasticity.hpp>

namespace Feel
{
namespace FeelModels
{

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::MixedElasticity( std::string const& prefix,
                                                      worldcomm_ptr_t const& worldComm,
                                                      std::string const& subPrefix,
                                                      ModelBaseRepository const& modelRep )
    : super_type( prefix, worldComm, subPrefix, modelRep ),
      ModelBase( prefix, worldComm, subPrefix, modelRep ),
      M_tauCst(doption (prefixvm(this->prefix(), "tau_constant") )),
      M_tauOrder(ioption( prefixvm(this->prefix(), "tau_order") )),
      M_hFace(ioption( prefixvm(this->prefix(), "hface")) ),
      M_useSC(boption(prefixvm(this->prefix(), "use-sc"))),
      M_nullspace(boption(prefixvm(this->prefix(), "nullspace")) )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    if (M_useSC)
    {
        if ( this->prefix().empty())
            M_backend = Feel::backend( _name="sc", _rebuild=true);
        else
            M_backend = Feel::backend( _name=prefixvm(prefix,"sc"), _rebuild=true);
    }
    else
    {
        if ( this->prefix().empty())
            M_backend = Feel::backend( _rebuild=true);
        else
            M_backend = Feel::backend( _name=prefix, _rebuild=true);
    }

    M_useUserIBC = false;

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticityConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticitySolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void MixedElasticity<Dim,Order,G_Order, E_Order>::setIBCList( std::vector<std::string> markersIBC )
{
    M_useUserIBC = true;
    M_IBCList.clear();
    for( auto const& marker : markersIBC )
    {
        ExpressionStringAtMarker exAtMark(std::make_tuple("expression", marker, std::string(""), std::string(""), std::string("")));
        M_IBCList.push_back(exAtMark);
    }
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::initTimeStep()
{
    // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_nm_mixedelasticity -> start( M_pp );
        // up current time
        this->updateTime( M_nm_mixedelasticity -> time() );
    }
    else
    {
        // start time step
        M_nm_mixedelasticity->restart();
        // load a previous solution as current solution
        M_pp = M_nm_mixedelasticity->previousUnknown();
        // up initial time
        this->setTimeInitial( M_nm_mixedelasticity->timeInitial() );
        // restart exporter
        //this->restartPostProcess();
        // up current time
        this->updateTime( M_nm_mixedelasticity->time() );

        this->log("MixedElasticity","initTimeStep", "restart nm/exporter done" );
    }

}
MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::updateTimeStepNM()
{
    this->log("MixedElasticity","updateTimeStepNM", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    // int previousTimeOrder = this->timeStepNM()->timeOrder();

    M_nm_mixedelasticity->next( M_pp );
    this->timeStepNM()->updateFromDisp( M_pp );

    // int currentTimeOrder = this->timeStepNM()->timeOrder();

    this->updateTime( M_nm_mixedelasticity->time() );


    this->timerTool("TimeStepping").stop("updateNm");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("MixedElasticity","updateTimeStepNM", "finish" );
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::createTimeDiscretization()
{
    this->log("MixedElasticity","createTimeDiscretization", "start" );
    this->timerTool("Constructor").start();


    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    auto dt = this->timeStep();

    if ( myFileFormat == "binary" )
        suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();

    M_nm_mixedelasticity = newmark( _vm=Environment::vm(), _space=M_Wh,
                                    _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"newmark"+suffixName)),
                                    _prefix="",
                                    _initial_time=this->timeInitial(),
                                    _final_time=this->timeFinal(),
                                    _time_step=this->timeStep(),
                                    _restart=this->doRestart(), _restart_path=this->restartPath(),_restart_at_last_save=this->restartAtLastSave(),
                                    _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() );
    M_nm_mixedelasticity->setfileFormat( myFileFormat );
    M_nm_mixedelasticity->setPathSave( (fs::path(this->rootRepository()) /
                                        fs::path( prefixvm(this->prefix(), (boost::format("newmark_dt_%1%")%dt).str() ) ) ).string() );




    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("MixedElasticity","createTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
typename MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::self_ptrtype
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::New( std::string const& prefix,
                                          worldcomm_ptr_t const& worldComm, std::string const& subPrefix,
                                          ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type> ( prefix,worldComm,subPrefix,modelRep );
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::init( mesh_ptrtype mesh, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( _mesh=new mesh_type);
    else
        M_mesh = mesh;
    M_timers["mesh"].push_back(toc("initMesh", this->verbose() || FLAGS_v > 0));

    tic();
    this->initModel();
    M_timers["initModel"].push_back(toc("initModel", this->verbose() || FLAGS_v > 0));

    tic();
    this->initSpaces();
    M_timers["spaces"].push_back(toc("initSpaces", this->verbose() || FLAGS_v > 0));

    if (!isStationary())
    {
        tic();
        this->createTimeDiscretization();
        this->initTimeStep();
        toc("time_discretization", this->verbose() || FLAGS_v > 0);
    }

    tic();
    this->initExporter(meshVisu);
    M_timers["exporter"].push_back(toc("initExporter", this->verbose() || FLAGS_v > 0));


    tic();
    this->assemble();
    M_timers["asbMatrix"].push_back(toc("assembleMatrix", this->verbose() || FLAGS_v > 0));


}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::initModel()
{

    // initialize marker lists for each boundary condition type
    // Strain
    auto itField = modelProperties().boundaryConditions().find( "stress");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );


        if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann_scalar" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann scalar:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }

        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann computed from displacement:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Robin:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
    }

    // Displacement
    itField = modelProperties().boundaryConditions().find( "displacement");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );

        if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
    }


    if ( !M_IBCList.empty() )
    {
        M_IBCList.clear();
    }
    itField = modelProperties().boundaryConditions().find( "stress");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );

        if ( itType != mapField.end() )
        {
            Feel::cout << "Integral:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                M_IBCList.push_back(exAtMarker);
            }
            Feel::cout << std::endl;
        }
    }


    if ( M_IBCList.empty() )
        M_integralCondition = 0;
    else
        M_integralCondition = M_IBCList.size();


}


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::initSpaces()
{

    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement(faces(M_mesh),[this]( auto const& e ) {
                                                                for( auto exAtMarker : this->M_IBCList)
                                                                {
                                                                    if ( e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                                                                        return true;
                                                                }
                                                                return false;
                                                            });

    auto face_mesh = createSubmesh( _mesh=M_mesh, _range=complement_integral_bdy, _update=0 );


    M_Vh = Pdhms<Order>( M_mesh, true );
    M_Wh = Pdhv<Order>( M_mesh, true );
    M_Mh = Pdhv<Order>( face_mesh, true );
    M_M0h = Pdh<0>( face_mesh );

    std::vector<std::string> ibc_markers(M_integralCondition);
    for( int i = 0; i < M_integralCondition; i++)
    {
        ibc_markers.push_back(M_IBCList[i].marker());
    }

    auto ibc_mesh = createSubmesh( _mesh=M_mesh, _range=markedfaces(M_mesh, ibc_markers), _update=0 );
    M_Ch = Pchv<0>( ibc_mesh, true );
    // M_Ch = Pchv<0>( M_mesh, true );

    auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >( M_integralCondition, M_Ch);
    M_ps = std::make_shared<product2_space_type>(product2(ibcSpaces,M_Vh,M_Wh,M_Mh));

    // M_ps = std::make_shared<product_space_std>(product(M_Vh,M_Wh,M_Mh));

    M_up = M_Vh->element( "u" ); // Strain
    M_pp = M_Wh->element( "p" ); // Displacement

    Feel::cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
               << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
               << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        Feel::cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;
}

} // namespace Feel
} // namespace FeelModels
