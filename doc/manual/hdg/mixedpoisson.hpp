#ifndef _MIXEDPOISSON_HPP
#define _MIXEDPOISSON_HPP 

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/complement.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelts/bdf.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>


namespace Feel {

namespace FeelModels {

inline
po::options_description
makeMixedPoissonOptions( std::string prefix = "mixedpoisson" )
{
    po::options_description mpOptions( "Mixed Poisson HDG options");
    mpOptions.add_options()
        ( prefixvm( prefix, "gmsh.submesh").c_str(), po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix, "tau_order").c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "picard.itol").c_str(), po::value<double>()->default_value( 1e-4 ), "tolerance" )
        ( prefixvm( prefix, "picard.itmax").c_str(), po::value<int>()->default_value( 10 ), "iterations max" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( prefixvm( prefix, "conductivity_json").c_str(), po::value<std::string>()->default_value( "cond" ), "key for conductivity in json" )
	( prefixvm( prefix, "p_exact").c_str(), po::value<std::string>()->default_value( "1.0" ), "p exact" )
        ( prefixvm( prefix, "conductivityNL_json").c_str(), po::value<std::string>()->default_value( "condNL" ), "key for non linear conductivity in json (depends on potential p)" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ;
    mpOptions.add ( envfeelmodels_options( prefix ) ).add( modelnumerical_options( prefix ) );
    return mpOptions;
}

inline po::options_description
makeMixedPoissonLibOptions( std::string prefix = "mixedpoisson" )
{
    po::options_description mpLibOptions( "Mixed Poisson HDG Lib options");
    // if ( !prefix.empty() )
    //    mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order, int G_Order = 1>
class MixedPoisson    :	public ModelNumerical
{
public:
    typedef ModelNumerical super_type;

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedPoisson<Dim,Order,G_Order> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,G_Order,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;
    
    // Vh
    using Vh_t = Pdhv_type<mesh_type,Order>;
    using Vh_ptr_t = Pdhv_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
    // Wh
    using Wh_t = Pdh_type<mesh_type,Order>;
    using Wh_ptr_t = Pdh_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
    // Mh
    using Mh_t = Pdh_type<face_mesh_type,Order>;
    using Mh_ptr_t = Pdh_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
    // Ch
    using Ch_t = Pch_type<mesh_type,0>;
    using Ch_ptr_t = Pch_ptrtype<mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    // M0h
    using M0h_t = Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t = Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;

    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    using linearAssembly_function_type = boost::function<void ( sparse_matrix_ptrtype& A,vector_ptrtype& F )>;

    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef boost::shared_ptr <exporter_type> exporter_ptrtype;
   
    // time
    // typedef Lagrange<Order,Scalar,Discontinuous> basis_scalar_type;
    // typedef FunctionSpace<mesh_type, bases<basis_scalar_type>> space_mixedpoisson_type;
    
    // typedef Bdf<space_mixedpoisson_type>  bdf_type;
    typedef Bdf <Wh_t> bdf_type;
    // typedef Bdf<Pdh_type<mesh_type,Order>> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    
//private:
protected:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh; // potential trace
    Ch_ptr_t M_Ch; // Lagrange multiplier
    M0h_ptr_t M_M0h;

    backend_ptrtype M_backend;
    BlocksBaseGraphCSR M_hdg_graph;
    sparse_matrix_ptrtype M_A;
    sparse_matrix_ptrtype M_A_cst;
    BlocksBaseVector<double> M_hdg_vec;
    vector_ptrtype M_F;
    BlocksBaseVector<double> M_hdg_sol;
    vector_ptrtype M_U;

    Vh_element_ptr_t M_up; // flux solution
    Wh_element_ptr_t M_pp; // potential solution 
    Ch_element_ptr_t M_mup; // potential solution on the integral boundary condition

    // time discretization
    bdf_ptrtype M_bdf_mixedpoisson;
    
    // map_scalar_field<2> M_dirichlet;
    // map_vector_field<Dim,1,2> M_neumann;
    // map_scalar_field<2> M_source;    

    int M_tau_order;
    
    bool M_integralCondition;
    bool M_isPicard;

    std::list<std::string> M_integralMarkersList;
    
    
    void initGraphs();
    void initGraphsWithIntegralCond();
    void assembleACst();

public:
    linearAssembly_function_type M_updateAssembly;
    
    // constructor
    // MixedPoisson( std::string prefix = "" );
    MixedPoisson( std::string const& prefix = "mixedpoisson",                   
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    
    MixedPoisson( self_type const& MP ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedpoisson",
                             WorldComm const& worldComm = Environment::worldComm(),
                             std::string const& subPrefix = "",
                             std::string const& rootRepository = ModelBase::rootRepositoryByDefault() ); 

    using op_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Wh_t, Pdh_type<mesh_type,1>>>;
    using opv_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Vh_t, Pdhv_type<mesh_type,1>>>;
    void init( mesh_ptrtype mesh = nullptr, int extraRow = 0, int extraCol = 0, mesh_ptrtype meshVisu = nullptr);
               
    virtual void initModel();
    virtual void initSpaces();
    virtual void initGraphs(int extraRow, int extraCol);
    virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void assembleA();
    virtual void assembleF();
    void solve();
    void solveNL();
    template<typename ExprT>
    void updateConductivityTerm( Expr<ExprT> expr, std::string marker = "");
    void updateConductivityTerm( bool isNL = false);
    template<typename ExprT>
    void updatePotentialRHS( Expr<ExprT> expr, std::string marker = "");
    template<typename ExprT>
    void updateFluxRHS( Expr<ExprT> expr, std::string marker = "");
    void computeError();    

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    M0h_ptr_t traceSpaceOrder0() const { return M_M0h; }
    Ch_ptr_t constantSpace() const {return M_Ch;}
    
    Vh_element_ptr_t fluxField() const { return M_up; }
    Wh_element_ptr_t potentialField() const { return M_pp; }
    model_prop_type modelProperties() { return *M_modelProperties; }
    model_prop_type modelProperties() const { return *M_modelProperties; }
    std::list<std::string> integralMarkersList() const { return M_integralMarkersList; }
    bool integralCondition() const { return M_integralCondition; }
    int tau_order() const { return M_tau_order; }
    backend_ptrtype get_backend() { return M_backend; }
    
    // time step scheme
    virtual void createTimeDiscretization() ;
    bdf_ptrtype timeStepBDF() { return M_bdf_mixedpoisson; }
    bdf_ptrtype const& timeStepBDF() const { return M_bdf_mixedpoisson; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    virtual void updateTimeStepBDF();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepBDF(); }
   
    // Exporter
    // void exportResults() { this->exportResults (this->currentTime() ); M_exporter->save(); }
    // virtual void exportResults ( double Time) ;
    void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr )
        {
            this->exportResults (this->currentTime(), mesh, Idh, Idhv );
        }
    void exportResults ( double Time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  ) ;
    exporter_ptrtype M_exporter;
    exporter_ptrtype exporterMP() { return M_exporter; }
};



template<int Dim, int Order, int G_Order> 
void MixedPoisson<Dim, Order, G_Order>::initTimeStep()
{
        // start or restart time step scheme
    if (!this->doRestart())
    {
        // start time step
        M_bdf_mixedpoisson -> start( *M_pp );
        // up current time
        this->updateTime( M_bdf_mixedpoisson -> time() );
    }
    else
    {
        // start time step
        M_bdf_mixedpoisson->restart();
        // load a previous solution as current solution
        *M_pp = M_bdf_mixedpoisson->unknown(0);
        // up initial time
        this->setTimeInitial( M_bdf_mixedpoisson->timeInitial() );
        // restart exporter
        //this->restartPostProcess();
        // up current time
        this->updateTime( M_bdf_mixedpoisson->time() );

        this->log("MixedPoisson","initTimeStep", "restart bdf/exporter done" );
    }
    
}


template<int Dim, int Order, int G_Order>
void MixedPoisson<Dim, Order, G_Order>::updateTimeStepBDF()
{
    this->log("MixedPoisson","updateTimeStepBDF", "start" );
    this->timerTool("TimeStepping").setAdditionalParameter("time",this->currentTime());
    this->timerTool("TimeStepping").start();

    int previousTimeOrder = this->timeStepBDF()->timeOrder();

    M_bdf_mixedpoisson->next( *M_pp );

    int currentTimeOrder = this->timeStepBDF()->timeOrder();

    this->updateTime( M_bdf_mixedpoisson->time() );


    this->timerTool("TimeStepping").stop("updateBdf");
    if ( this->scalabilitySave() ) this->timerTool("TimeStepping").save();
    this->log("MixedPoisson","updateTimeStepBDF", "finish" );
}


template<int Dim, int Order, int G_Order>
MixedPoisson<Dim, Order, G_Order>::MixedPoisson( std::string const& prefix,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    : super_type( prefix, worldComm, subPrefix, rootRepository )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedPoisson","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());


    this->setFilenameSaveInfo( prefixvm(this->prefix(),"MixedPoisson.info") );


    M_prefix = prefix;
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedPoissonSolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedPoisson","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


template<int Dim, int Order, int G_Order>
typename MixedPoisson<Dim,Order, G_Order>::self_ptrtype
MixedPoisson<Dim,Order,G_Order>::New( std::string const& prefix,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type> ( prefix,worldComm,subPrefix,rootRepository );
}

/*
template<int Dim, int Order, int G_Order>
MixedPoisson<Dim, Order, G_Order>::MixedPoisson(std::string prefix )
{
    M_prefix = prefix;
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );
}
*/

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::init( mesh_ptrtype mesh, int extraRow, int extraCol, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;
    toc("mesh");

    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);


    tic();
    this->initModel();
    toc("model");

    tic();
    this->initSpaces();
    toc("spaces");

    if(!this->isStationary()){
        tic();
        this->createTimeDiscretization();
        this->initTimeStep();
        toc("timeDiscretization",true);
    }


    tic();
    this->initGraphs(extraRow, extraCol);
    M_A = this->get_backend()->newBlockMatrix(_block=M_hdg_graph);
    M_A_cst = this->get_backend()->newBlockMatrix(_block=M_hdg_graph);
    M_F = this->get_backend()->newBlockVector(_block=M_hdg_vec, _copy_values=false);
    M_U = this->get_backend()->newBlockVector(_block=M_hdg_sol, _copy_values=false);
    toc("graphs");
    tic();
    this->initExporter( meshVisu );
    toc("exporter");
    tic();
    this->assembleA();
    M_A_cst->close();
    toc("assemble");
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::initModel()
{
    // M_dirichlet = M_modelProperties -> boundaryConditions().getScalarFields( "potential", "Dirichlet");
    // M_neumann = M_modelProperties -> boundaryConditions().template getVectorFields<Dim> ( "flux", "Neumann" );
    // M_source = M_modelProperties -> boundaryConditions().getScalarFields( "potential", "SourceTerm");

    // initialize marker lists for each boundary condition type
    M_integralMarkersList.clear();
    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        /*
	if ( itType != mapField.end() )
        {
            cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    cout << " " << marker;
                else
                    cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            cout << std::endl;
        }
	*/ 
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            cout << "Neumann:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    cout << " " << marker;
                else
                    cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            cout << std::endl;
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            cout << "Robin:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    cout << " " << marker;
                else
                    cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            cout << std::endl;
        }
    }
    itField = M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Integral" );
        if ( itType != mapField.end() )
        {
            cout << "Integral:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    cout << " " << marker;
                else
                    cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
                M_integralMarkersList.push_back(marker);
            }
            cout << std::endl;
        }
    }

    if ( M_integralMarkersList.empty() )
        M_integralCondition = false;
    else
        M_integralCondition = true;
    if ( boost::icontains(M_modelProperties->model(),"picard") )
        M_isPicard = true;
    else
        M_isPicard = false;
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::initSpaces()
{
    // Mh only on the faces whitout integral condition
    auto complement_integral_bdy = complement(faces(M_mesh),[this]( auto const& e ) {
            for( auto marker : this->M_integralMarkersList)
            {
                if ( e.marker().value() == this->M_mesh->markerName( marker ) )
                    return true;
            }
            return false; });
    auto face_mesh = createSubmesh( M_mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );

    M_Vh = Pdhv<Order>( M_mesh, true );
    M_Wh = Pdh<Order>( M_mesh, true );
    M_Mh = Pdh<Order>( face_mesh, true );
    M_Ch = Pch<0>( M_mesh );
    M_M0h = Pdh<0>( face_mesh, true );

    M_up = M_Vh->elementPtr( "u" );
    M_pp = M_Wh->elementPtr( "p" );
    

    cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;
}


template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::initExporter( mesh_ptrtype meshVisu ) 
{
    std::string geoExportType="static"; //change_coords_only, change, static
    M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh() ,
                            _name="Export",
                            _geo=geoExportType );
        //_path=this->exporterPath() ); 
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim,Order,G_Order>::createTimeDiscretization()
{
    this->log("MixedPoisson","createTimeDiscretization", "start" );
    this->timerTool("Constructor").start();
    
    
    std::string myFileFormat = soption(_name="ts.file-format");// without prefix
    std::string suffixName = "";
    if ( myFileFormat == "binary" )
         suffixName = (boost::format("_rank%1%_%2%")%this->worldComm().rank()%this->worldComm().size() ).str();
    M_bdf_mixedpoisson = bdf( _vm=Environment::vm(), _space=M_Wh ,
                       _name=prefixvm(this->prefix(),prefixvm(this->subPrefix(),"p"+suffixName)) ,
                       _prefix="",
                       _initial_time=this->timeInitial(),
                       _final_time=this->timeFinal(),
                       _time_step=this->timeStep(),
                       _restart=this->doRestart(),
                       _restart_path=this->restartPath(),
                       _restart_at_last_save=this->restartAtLastSave(),
                       _save=this->tsSaveInFile(), _freq=this->tsSaveFreq() ); 
    M_bdf_mixedpoisson->setfileFormat( myFileFormat );
    M_bdf_mixedpoisson->setPathSave( (fs::path(this->rootRepository()) /
                               fs::path( prefixvm(this->prefix(), (boost::format("bdf_o_%1%_dt_%2%")%M_bdf_mixedpoisson->bdfOrder()%this->timeStep() ).str() ) ) ).string() );

    double tElapsed = this->timerTool("Constructor").stop("createTimeDiscr");
    this->log("MixedPoisson","createTimeDiscretization", (boost::format("finish in %1% s") %tElapsed).str() );
}

// when a derived class need to add extra block,
// redefine initGraphs by calling the super method
// and just specify the extra block
template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::initGraphs(int extraRow, int extraCol)
{
    int baseRow = M_integralCondition ? 4 : 3;
    int baseCol = M_integralCondition ? 4 : 3;
    M_hdg_graph = BlocksBaseGraphCSR(baseRow+extraRow,baseCol+extraCol);
    M_hdg_vec = BlocksBaseVector<double>(baseRow+extraRow);
    M_hdg_sol = BlocksBaseVector<double>(baseCol+extraCol);

    M_hdg_graph(0,0) = stencil( _test=M_Vh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(1,0) = stencil( _test=M_Wh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(2,0) = stencil( _test=M_Mh,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();

    M_hdg_graph(0,1) = stencil( _test=M_Vh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(1,1) = stencil( _test=M_Wh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(2,1) = stencil( _test=M_Mh,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();

    M_hdg_graph(0,2) = stencil( _test=M_Vh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(1,2) = stencil( _test=M_Wh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();
    M_hdg_graph(2,2) = stencil( _test=M_Mh,_trial=M_Mh, _diag_is_nonzero=false, _close=false)->graph();

    M_hdg_vec(0,0) = M_backend->newVector( M_Vh );
    M_hdg_vec(1,0) = M_backend->newVector( M_Wh );
    M_hdg_vec(2,0) = M_backend->newVector( M_Mh );

    auto phatp = M_Mh->elementPtr( "phat" );

    M_hdg_sol(0,0) = M_up;
    M_hdg_sol(1,0) = M_pp;
    M_hdg_sol(2,0) = phatp;

    if ( M_integralCondition )
    {
        M_hdg_graph(3,0) = stencil( _test=M_Ch,_trial=M_Vh, _diag_is_nonzero=false, _close=false)->graph();
        M_hdg_graph(3,1) = stencil( _test=M_Ch,_trial=M_Wh, _diag_is_nonzero=false, _close=false)->graph();
        M_hdg_graph(3,2) = stencil( _test=M_Ch,_trial=M_Mh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
        M_hdg_graph(0,3) = stencil( _test=M_Vh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
        M_hdg_graph(1,3) = stencil( _test=M_Wh,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();
        M_hdg_graph(2,3) = stencil( _test=M_Mh,_trial=M_Ch, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
        M_hdg_graph(3,3) = stencil( _test=M_Ch,_trial=M_Ch, _diag_is_nonzero=false, _close=false)->graph();

        M_hdg_vec(3,0) = M_backend->newVector( M_Ch );
        M_mup = M_Ch->elementPtr( "mup" );
        M_hdg_sol(3,0) = M_mup;
    }
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::solve()
{
    tic();
    M_modelProperties -> parameters().updateParameterValues();
    updateConductivityTerm();
    assembleF();
    if ( M_updateAssembly != NULL )
        this->M_updateAssembly( M_A, M_F );
    M_backend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U);
    M_hdg_sol.localize(M_U);
    toc("solve", true);
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::solveNL()
{
    tic();
    if ( M_updateAssembly != NULL )
        this->M_updateAssembly( M_A, M_F );
    M_backend->solve( _matrix=M_A, _rhs=M_F, _solution=M_U);
    M_hdg_sol.localize(M_U);
    toc("solve", true);
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::assembleA()
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    auto p = M_Wh->element( "p" );
    auto q = M_Wh->element( "p" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto phat = M_Mh->element( "phat" );
    auto l = M_Mh->element( "lambda" );
    auto H = M_M0h->element( "H" );
    if ( ioption(prefixvm(M_prefix, "hface") ) == 0 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 1 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 2 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
        H.on( _range=elements(M_M0h->mesh()), _expr=h() );
    // stabilisation parameter
    auto tau_constant = cst(doption(prefixvm(M_prefix, "tau_constant")));

    auto gammaMinusIntegral = complement(boundaryfaces(M_mesh),[this]( auto const& e ) {
            for( auto marker : this->M_integralMarkersList)
            {
                if ( e.marker().value() == this->M_mesh->markerName( marker ) )
                    return true;
            }
            return false; });

    auto a12 = form2( _trial=M_Wh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0, _colstart=1 );
    auto a13 = form2( _trial=M_Mh, _test=M_Vh,_matrix=M_A_cst,
                      _rowstart=0, _colstart=2);
    auto a21 = form2( _trial=M_Vh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=0);
    auto a22 = form2( _trial=M_Wh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=1 );
    auto a23 = form2( _trial=M_Mh, _test=M_Wh,_matrix=M_A_cst,
                      _rowstart=1, _colstart=2);
    auto a31 = form2( _trial=M_Vh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2, _colstart=0);
    auto a32 = form2( _trial=M_Wh, _test=M_Mh,_matrix=M_A_cst,
                      _rowstart=2, _colstart=1);
    auto a33 = form2(_trial=M_Mh, _test=M_Mh,_matrix=M_A_cst,
                     _rowstart=2, _colstart=2);

    // -(p,div(v))_Omega
    a12 += integrate(_range=elements(M_mesh),_expr=-(idt(p)*div(v)));

    // <phat,v.n>_Gamma\Gamma_I
    a13 += integrate(_range=internalfaces(M_mesh),
                     _expr=( idt(phat)*leftface(trans(id(v))*N())+
                             idt(phat)*rightface(trans(id(v))*N())) );
    a13 += integrate(_range=gammaMinusIntegral,
                     _expr=idt(phat)*trans(id(v))*N());


    // -(j, grad(w))
    a21 += integrate(_range=elements(M_mesh),_expr=(-grad(w)*idt(u)));
    // <j.n,w>_Gamma
    a21 += integrate(_range=internalfaces(M_mesh),
                     _expr=( leftface(id(w))*leftfacet(trans(idt(u))*N()) ) );
    a21 += integrate(_range=internalfaces(M_mesh),
                     _expr=(rightface(id(w))*rightfacet(trans(idt(u))*N())) );
    a21 += integrate(_range=boundaryfaces(M_mesh),
                     _expr=(id(w)*trans(idt(u))*N()));


    // <tau p, w>_Gamma
    a22 += integrate(_range=internalfaces(M_mesh),
                     _expr=tau_constant *
                     ( leftfacet( pow(idv(H),M_tau_order)*idt(p))*leftface(id(w)) +
                       rightfacet( pow(idv(H),M_tau_order)*idt(p))*rightface(id(w) )));
    a22 += integrate(_range=boundaryfaces(M_mesh),
                     _expr=(tau_constant * pow(idv(H),M_tau_order)*id(w)*idt(p)));
    
    // (1/delta_t p, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() ) {
        a22 += integrate(_range=elements(M_mesh),
                         _expr = (this->timeStepBDF()->polyDerivCoefficient(0)*idt(p)*id(w)) );
    }
    
    // <-tau phat, w>_Gamma\Gamma_I
    a23 += integrate(_range=internalfaces(M_mesh),
                     _expr=-tau_constant * idt(phat) *
                     ( leftface( pow(idv(H),M_tau_order)*id(w) )+
                       rightface( pow(idv(H),M_tau_order)*id(w) )));
    a23 += integrate(_range=gammaMinusIntegral,
                     _expr=-tau_constant * idt(phat) * pow(idv(H),M_tau_order)*id(w) );


    // <j.n,mu>_Omega/Gamma
    a31 += integrate(_range=internalfaces(M_mesh),
                     _expr=( id(l)*(leftfacet(trans(idt(u))*N())+
                                    rightfacet(trans(idt(u))*N())) ) );


    // <tau p, mu>_Gamma_N
    a32 += integrate(_range=internalfaces(M_mesh),
                     _expr=tau_constant * id(l) * ( leftfacet( pow(idv(H),M_tau_order)*idt(p) )+
                                                    rightfacet( pow(idv(H),M_tau_order)*idt(p) )));


    // <-tau phat, mu>_Omega/Gamma
    a33 += integrate(_range=internalfaces(M_mesh),
                     _expr=-tau_constant * idt(phat) * id(l) * ( leftface( pow(idv(H),M_tau_order) )+
                                                                 rightface( pow(idv(H),M_tau_order) )));

    // BC
    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <phat, mu>_Gamma_D
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=idt(phat) * id(l) );
            }
        }
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                // <j.n,mu>_Gamma_N
                a31 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));
                // <tau p, mu>_Gamma_N
                a32 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(idv(H),M_tau_order)*idt(p) ) );
                // <-tau phat, mu>_Gamma_N
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(idv(H),M_tau_order) ) );
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression1());
                // <j.n,mu>_Gamma_R
                a31 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=( id(l)*(trans(idt(u))*N()) ));
                // <tau p, mu>_Gamma_R
                a32 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=tau_constant * id(l) * ( pow(idv(H),M_tau_order)*idt(p) ) );
                // <-tau phat, mu>_Gamma_R
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=-tau_constant * idt(phat) * id(l) * ( pow(idv(H),M_tau_order) ) );
                // <g_R^1 phat, mu>_Gamma_R
                a33 += integrate(_range=markedfaces(M_mesh,marker),
                                 _expr=g*idt(phat) * id(l) );
            }
        }
    }

    if ( M_integralCondition )
    {
        auto a14 = form2(_trial=M_Ch, _test=M_Vh,_matrix=M_A_cst,
                         _rowstart=0,
                         _colstart=3);
        auto a34 = form2(_trial=M_Ch, _test=M_Mh,_matrix=M_A_cst,
                         _rowstart=2,
                         _colstart=3);
        auto a41 = form2(_trial=M_Vh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=0);
        auto a42 = form2(_trial=M_Wh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=1);
        auto a43 = form2(_trial=M_Mh, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=2);
        auto a24 = form2(_trial=M_Ch, _test=M_Wh,_matrix=M_A_cst,
                         _rowstart=1,
                         _colstart=3);
        auto a44 = form2(_trial=M_Ch, _test=M_Ch,_matrix=M_A_cst,
                         _rowstart=3,
                         _colstart=3);

        auto itField = M_modelProperties->boundaryConditions().find( "flux");
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    // <lambda, v.n>_Gamma_I
                    a14 += integrate( _range=markedfaces(M_mesh,marker),
                                      _expr=trans(id(u))*N()*idt(nu) );

                    // <lambda, tau w>_Gamma_I
                    a24 += integrate( _range=markedfaces(M_mesh,marker),
                                      _expr=-tau_constant * ( pow(idv(H),M_tau_order)*id(w) ) * idt(nu) );

                    // <j.n, m>_Gamma_I
                    a41 += integrate( _range=markedfaces(M_mesh,marker), _expr=trans(idt(u))*N()*id(nu) );

                    // <tau p, m>_Gamma_I
                    a42 += integrate( _range=markedfaces(M_mesh,marker), _expr=tau_constant *
                          ( pow(idv(H),M_tau_order)*idt(p) ) * id(nu) );

                    // -<lambda2, m>_Gamma_I
                    a44 += integrate( _range=markedfaces(M_mesh,marker), _expr=-pow(idv(H),M_tau_order)*id(nu)*idt(nu) );
                }
            }
        }
    }

    
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::assembleF()
{
    M_F->zero();
    auto nu = M_Ch->element( "nu" );
    auto l = M_Mh->element( "lambda" );
    auto w = M_Wh->element();
    auto v = M_Vh->element();

    // Building the RHS

    auto rhs1 = form1( _test=M_Wh, _vector=M_F, _rowstart=0);
    auto rhs2 = form1( _test=M_Wh, _vector=M_F, _rowstart=1);
    auto rhs3 = form1( _test=M_Mh, _vector=M_F, _rowstart=2);
    /*
    M_source.setParameterValues( M_modelProperties->parameters().toParameterValues() );
    if (!this->isStationary())
        M_source.setParameterValues( {{"t",M_bdf_mixedpoisson->time()}} ); 
    for ( auto const& s : M_source ){
	// (f, w)_Omega
	rhs2 += integrate( _range = markedelements(M_mesh,marker(s)), 
			   _expr = expression(s)*id(w) );
        // (p_old,w)_Omega
        if ( !this->isStationary() ){
            rhs2 += integrate( _range = markedelements(M_mesh,marker(s)),
                               _expr = idv(this->timeStepBDF()->polyDeriv()) * id(w));
        }
    }
    // Dirichlet potential BC
    M_dirichlet.setParameterValues( M_modelProperties->parameters().toParameterValues() );
    if (!this->isStationary() )
        M_dirichlet.setParameterValues( {{"t",M_bdf_mixedpoisson->time()}} );
    for ( auto const& d : M_dirichlet ){
        // <g_D, mu>_Gamma_D
        rhs3 += integrate( _range = markedfaces(M_mesh,marker(d)),
                           _expr = id(l)*expression(d));

    }*/
    auto itField = M_modelProperties->boundaryConditions().find( "potential");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // (f, w)_Omega
                rhs2 += integrate( _range=markedelements(M_mesh,marker),
                                   _expr=g*id(w));
                // (p_old,w)_Omega
                if ( !this->isStationary() ){
                    rhs2 += integrate( _range=markedelements(M_mesh,marker),
                                      _expr= idv(this->timeStepBDF()->polyDeriv()) * id(w));
                }
            }
        }
	
        itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( exAtMarker.isExpression() )
                {
                    auto g = expr(exAtMarker.expression());
		    // <g_D, mu>_Gamma_D
		    rhs3 += integrate(_range=markedfaces(M_mesh,marker),
                                  _expr=id(l)*g);                    
                }
                else if ( exAtMarker.isFile() )
                {
                    double g = 0;
                    if ( !this->isStationary() )
                    {
                    	std::cout << "use data file to set rhs for Dirichlet BC at time " << M_bdf_mixedpoisson->time() << std::endl;
                        LOG(INFO) << "use data file to set rhs for Dirichlet BC at time " << M_bdf_mixedpoisson->time() << std::endl;

                        // data may depend on time
                        g = exAtMarker.data(M_bdf_mixedpoisson->time());
                    }
                    else
                    	g = exAtMarker.data(0.1);

                    LOG(INFO) << "use g=" << g << std::endl;
                    std::cout << "g=" << g << std::endl;
 		    // <g_D, mu>_Gamma_D
		    rhs3 += integrate(_range=markedfaces(M_mesh,marker),
                                      _expr=id(l)*g);
                }

            }
        }

        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression());
                // <g_N,mu>_Gamma_N
                rhs3 += integrate( _range=markedfaces(M_mesh, marker),
                                   _expr=id(l)*g);
            }
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr(exAtMarker.expression2());
                // <g_R^2,mu>_Gamma_R
                rhs3 += integrate( _range=markedfaces(M_mesh, marker),
                                   _expr=id(l)*g);
            }
        }
    }

    itField = M_modelProperties->boundaryConditions().find( "flux");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "SourceTerm" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                auto g = expr<3,1>(exAtMarker.expression());
                // (g, v)_Omega
                rhs2 += integrate( _range=markedelements(M_mesh,marker),
                                   _expr=inner(g,id(v)));
            }
        }
    }

    if ( M_integralCondition )
    {
        auto rhs4 = form1( _test=M_Ch, _vector=M_F,
                           _rowstart=3);
        itField = M_modelProperties->boundaryConditions().find( "flux");
        if ( itField != M_modelProperties->boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Integral" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    double meas = integrate( _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
                    if ( exAtMarker.isExpression() )
                    {
                        auto g = expr(exAtMarker.expression());
                        // <I_target,m>_Gamma_I
                        rhs4 += integrate(_range=markedfaces(M_mesh,marker),
                                          _expr=g*id(nu)/meas);
                    }
                    else if ( exAtMarker.isFile() )
                    {
                        double g = 0;
                        if ( !this->isStationary() )
                        {
                            std::cout << "use data file to set rhs for IBC at time " << M_bdf_mixedpoisson->time() << std::endl;
                            LOG(INFO) << "use data file to set rhs for IBC at time " << M_bdf_mixedpoisson->time() << std::endl;
                            
                            // data may depend on time
                            g = exAtMarker.data(M_bdf_mixedpoisson->time());
                        }
                        else
                            g = exAtMarker.data(0.1);
                            
                            LOG(INFO) << "use g=" << g << std::endl;
                            std::cout << "g=" << g << std::endl;
                            rhs4 += integrate(_range=markedfaces(M_mesh,marker),
                                              _expr=g*id(nu)/meas);
                    }
                }
            }
        }
    }
}

template<int Dim, int Order, int G_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updateConductivityTerm( Expr<ExprT> expr, std::string marker)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    MatConvert(toPETSc(M_A_cst)->mat(), MATSAME, MAT_INITIAL_MATRIX, &(toPETSc(M_A)->mat()));
    auto a11 = form2( _trial=M_Vh, _test=M_Vh,_matrix=M_A );
    if ( marker.empty() )
        a11 += integrate( _range=elements(M_mesh), _expr=inner(idt(u),id(v))/expr);
    else
        a11 += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/expr);
}

template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::updateConductivityTerm( bool isNL)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    MatConvert(toPETSc(M_A_cst)->mat(), MATSAME, MAT_INITIAL_MATRIX, &(toPETSc(M_A)->mat()));

    auto a11 = form2( _trial=M_Vh, _test=M_Vh, _matrix=M_A );
    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( !isNL )
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivity_json")));
            // (sigma^-1 j, v)
            a11 += integrate(_range=markedelements(M_mesh,marker),
                             _expr=(trans(idt(u))*id(v))/cond );
        }
        else
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivityNL_json")), "p", idv(*M_pp));
	    // (sigma(p)^-1 j, v)
            a11 += integrate(_range=markedelements(M_mesh,marker),
                             _expr=(trans(idt(u))*id(v))/cond );
        }
    }
}

template<int Dim, int Order, int G_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updatePotentialRHS( Expr<ExprT> expr, std::string marker)
{
    auto rhs = form1( _test=M_Wh, _vector=M_F, _rowstart=1);
    auto w = M_Wh->element();
    if ( marker.empty() )
        rhs += integrate(_range=elements(M_mesh), _expr=expr*id(w) );
    else
        rhs += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(w)) );
}

template<int Dim, int Order, int G_Order>
template<typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updateFluxRHS( Expr<ExprT> expr, std::string marker)
{
    auto rhs = form1( _test=M_Vh, _vector=M_F, _rowstart=0);
    auto v = M_Vh->element();
    if ( marker.empty() )
        rhs += integrate(_range=elements(M_mesh), _expr=expr*id(v) );
    else
        rhs += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(v)) );
}

/*
template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::exportResults()
{
    auto e = exporter( M_mesh);

    auto postProcess = M_modelProperties->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "flux" )
                e->add(prefixvm(M_prefix, "flux"), *M_up);
            if ( field == "potential" )
                e->add(prefixvm(M_prefix, "potential"), *M_pp);
        }
    }

    e->save();
}*/

template <int Dim, int Order, int G_Order>
void
MixedPoisson<Dim,Order, G_Order>::exportResults( double time, mesh_ptrtype mesh, op_interp_ptrtype Idh, opv_interp_ptrtype Idhv  )
{

    this->log("MixedPoisson","exportResults", "start");
    this->timerTool("PostProcessing").start();

    if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC && mesh  )
    {
        LOG(INFO) << "exporting on visualisation mesh at time " << time;
        M_exporter->step( time )->setMesh( mesh );
    }
    else if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
    {
        LOG(INFO) << "exporting on computational mesh at time " << time;
        M_exporter->step( time )->setMesh( M_mesh );
    }
    
    // Export computed solutions
    auto postProcess = M_modelProperties->postProcess();
    auto itField = postProcess.find( "Fields");
    if ( itField != postProcess.end() )
    {
        for ( auto const& field : (*itField).second )
        {
            if ( field == "flux" )
            {
                LOG(INFO) << "exporting flux at time " << time;
                M_exporter->step( time )->add(prefixvm(M_prefix, "flux"), Idhv?(*Idhv)( *M_up):*M_up );
                if (M_integralCondition)
                {
                    double j_integral = 0;
                    for( auto marker : this->M_integralMarkersList)
                    {
                        LOG(INFO) << "exporting integral flux at time " << time << " on marker " << marker;
                        j_integral += integrate(_range=markedfaces(M_mesh,marker),_expr=trans(idv(M_up))*N()).evaluate()(0,0);
                    }
                    M_exporter->step( time )->add(prefixvm(M_prefix, "integralFlux"), j_integral);
                }
            }
            if ( field == "potential" )
            {
                LOG(INFO) << "exporting potential at time " << time;
                M_exporter->step( time )->add(prefixvm(M_prefix, "potential"), Idh?(*Idh)(*M_pp):*M_pp);
                if (M_integralCondition)
                {
                    LOG(INFO) << "exporting IBC potential at time " << time << " value " << (*M_mup)[0];
                    M_exporter->step( time )->add(prefixvm(M_prefix, "cstPotential"),(*M_mup)[0] );
                }
            }
        }
    }
    M_exporter->save();
    /*/ Export exact solutions
     if ( this->isStationary() ){
     auto K = 10;
     auto p_exact = expr(soption(prefixvm(M_prefix,"p_exact") ));
     auto gradp_exact = grad<Dim>(p_exact);
     auto u_exact = -K*trans(gradp_exact);
        
     *
     for( auto const& pairMat : M_modelProperties->materials() )
     {
     auto marker = pairMat.first;
     auto material = pairMat.second;
     auto K = material.getScalar(soption(prefixvm(M_prefix,"conductivity_json")));
     u_exact = -K*trans(gradp_exact) ;
     }*

     auto p_exact_proj = project( _space=M_Wh, _range=elements(M_mesh), _expr=p_exact);
     auto u_exact_proj = project( _space=M_Vh, _range=elements(M_mesh), _expr=u_exact);
     M_exporter -> step (0) -> add(prefixvm(M_prefix, "p_exact"), p_exact_proj);

     M_exporter -> step (0) -> add(prefixvm(M_prefix, "u_exact"), u_exact_proj);
     }*/

    

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("MixedPoisson","exportResults", "finish");
}


template<int Dim, int Order, int G_Order>
void
MixedPoisson<Dim, Order, G_Order>::computeError(){

    auto K = 10;
    auto p_exact = expr(soption(prefixvm(M_prefix,"p_exact") ));
    auto gradp_exact = grad<Dim>(p_exact);
    auto u_exact = -K*trans(gradp_exact);
    /*
    for( auto const& pairMat : M_modelProperties->materials() )
    {
         auto marker = pairMat.first;
         auto material = pairMat.second;
         auto K = material.getScalar(soption(prefixvm(M_prefix,"conductivity_json")));
	 u_exact = -K*trans(gradp_exact) ;
    } */   
    
    tic();

    bool has_dirichlet = nelements(markedfaces(M_mesh,"Dirichlet"),true) >= 1;

    auto l2err_u = normL2( _range=elements(M_mesh), _expr=u_exact - idv(*M_up) );
    auto l2norm_uex = normL2( _range=elements(M_mesh), _expr=u_exact );
    if (l2norm_uex < 1)
	l2norm_uex = 1.0;
    
    auto l2err_p = normL2( _range=elements(M_mesh), _expr=p_exact - idv(*M_pp) );
    auto l2norm_pex = normL2( _range=elements(M_mesh), _expr=p_exact );
    if (l2norm_pex < 1)
	l2norm_pex = 1.0;
    // Feel::cout << "Has Dirichlet: " << has_dirichlet << std::endl;
   
    if ( !has_dirichlet ){
	auto mean_p_exact = mean( elements(M_mesh), p_exact )(0,0);
	auto mean_p = mean( _range=elements(this->mesh()), _expr=idv(*M_pp) )(0,0);
	l2err_p = normL2( elements(M_mesh), (p_exact - cst(mean_p_exact)) - (idv(*M_pp) - cst(mean_p)) );
    }
    
    Feel::cout << "============================================" << std::endl;
    Feel::cout << "||p-p_ex||_L2=\t" << l2err_p/l2norm_pex << std::endl;
    Feel::cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
    Feel::cout << "============================================" << std::endl;
    toc("error");


}



} // Namespace FeelModels

} // Namespace Feel

#endif
