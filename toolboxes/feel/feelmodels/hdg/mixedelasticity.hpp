#ifndef _MIXEDELASTICITY_HPP
#define _MIXEDELASTICITY_HPP 

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
#include <feel/feelts/newmark.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/projector.hpp>

// #define USE_SAME_MATH 1

namespace Feel {


template <typename SpaceType>
NullSpace<double> hdgNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
    return userNullSpace;
}

template <typename SpaceType>
NullSpace<double> hdgNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    return userNullSpace;
}



namespace FeelModels {

inline
po::options_description
makeMixedElasticityOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpOptions( "Mixed Elasticity HDG options");
    mpOptions.add_options()
        ( prefixvm( prefix, "gmsh.submesh").c_str(), po::value<std::string>()->default_value( "" ), "submesh extraction" )
        // ( "gmsh.submesh2", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ( prefixvm( prefix,"tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix,"tau_order").c_str(), po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( prefixvm( prefix, "use-sc").c_str(), po::value<bool>()->default_value(true), "use static condensation")           
        ( prefixvm( prefix, "nullspace").c_str(), po::value<bool>()->default_value( false ), "add null space" )
        ;
    mpOptions.add( modelnumerical_options( prefix ) );
	mpOptions.add ( backend_options( prefix+".sc" ) );
    return mpOptions;
}

inline po::options_description
makeMixedElasticityLibOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpLibOptions( "Mixed Elasticity HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order, int G_Order = 1, int E_Order = 4>
class MixedElasticity    :	public ModelNumerical
{
public:
    typedef ModelNumerical super_type;


    static const uint16_type expr_order = Order + E_Order;
    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedElasticity<Dim,Order,G_Order,E_Order> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,G_Order,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef std::shared_ptr<face_mesh_type> face_mesh_ptrtype;
    

// ---- //
    using Vh_t =  Pdhms_type<mesh_type,Order>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
// ---- //
    using Wh_t =  Pdhv_type<mesh_type,Order>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
// ---- //
    using Mh_t =  Pdhv_type<face_mesh_type,Order>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
// ---- //
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;    
// ---- //
    using Ch_t = Pchv_type<face_mesh_type,0>;
    using Ch_ptr_t = Pchv_ptrtype<face_mesh_type,0>;
    using Ch_element_t = typename Ch_t::element_type;
    using Ch_element_ptr_t = typename Ch_t::element_ptrtype;
    using Ch_element_vector_type = std::vector<Ch_element_t>;



	using op_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Wh_t, Pdhv_type<mesh_type,Order>>>;
    using opv_interp_ptrtype = std::shared_ptr<OperatorInterpolation<Vh_t, Pdhms_type<mesh_type,Order>>>;
 
    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    /* 
    using product_space_std = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
	using product_space_ptrtype = std::shared_ptr<product_space_std>;
    using bilinear_block_std = BlockBilinearForm<product_space_std>;
    */

    using product2_space_type = ProductSpaces2<Ch_ptr_t,Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using product2_space_ptrtype = std::shared_ptr<product2_space_type>;
    using integral_boundary_list_type = std::vector<ExpressionStringAtMarker>;
    
    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef std::shared_ptr <exporter_type> exporter_ptrtype;
    
    // typedef Newmark<space_mixedelasticity_type>  newmark_type;
    typedef Newmark <Wh_t> newmark_type;
    typedef std::shared_ptr<newmark_type> newmark_ptrtype;
    
//private:
protected:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // stress
    Wh_ptr_t M_Wh; // displacement
    Mh_ptr_t M_Mh; // displacement trace 
    M0h_ptr_t M_M0h;  
    Ch_ptr_t M_Ch; // Lagrange multiplier for IBC

	product2_space_ptrtype M_ps;

    backend_ptrtype M_backend;
    condensed_matrix_ptr_t<value_type> M_A_cst;
    condensed_vector_ptr_t<value_type> M_F;

    Vh_element_t M_up; // stress solution
    Wh_element_t M_pp; // displacement solution 
    Ch_element_vector_type M_mup; // displacement solution on the IBC

    double M_tau_constant;
    int M_tau_order;
    
    int M_integralCondition;
    int M_useUserIBC;
    integral_boundary_list_type M_IBCList;

    // time discretization
    newmark_ptrtype M_nm_mixedelasticity;

    // save time
    std::map<std::string, std::vector<double> > M_timers;

public:
    
    // constructor
    MixedElasticity( std::string const& prefix = "mixedelasticity",
                     worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                     std::string const& subPrefix = "",
                     ModelBaseRepository const& modelRep = ModelBaseRepository() );
    
    MixedElasticity( self_type const& ME ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedelasticity",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }
    M0h_ptr_t traceSpaceOrder0() const { return M_M0h; }
    Ch_ptr_t constantSpace() const {return M_Ch;}

    Vh_element_t fluxField() const { return M_up; }
    Wh_element_t potentialField() const { return M_pp; }
    model_prop_ptrtype modelProperties() { return M_modelProperties; }
    model_prop_ptrtype modelProperties() const { return M_modelProperties; }
    
    integral_boundary_list_type integralBoundaryList() const { return M_IBCList; }
    int integralCondition() const { return M_integralCondition; }
    void setIBCList(std::vector<std::string> markersIbc);
    product2_space_ptrtype getPS() const { return M_ps; }

    int tau_order() const { return M_tau_order; }
    backend_ptrtype get_backend() { return M_backend; }
    condensed_vector_ptr_t<value_type> getF() {return M_F; }
    std::map<std::string, std::vector<double> > getTimers() {return M_timers; }

    // Exporter
    virtual void exportResults( mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  )
    {
    	this->exportResults (this->currentTime(), mesh , Idh, Idhv);
        M_exporter -> save();
    }
    
    void exportResults ( double Time, mesh_ptrtype mesh = nullptr, op_interp_ptrtype Idh = nullptr, opv_interp_ptrtype Idhv = nullptr  ) ;
    exporter_ptrtype M_exporter;
    exporter_ptrtype exporterME() { return M_exporter; }
	
    void init( mesh_ptrtype mesh = nullptr, mesh_ptrtype meshVisu = nullptr);

    virtual void initModel();
    virtual void initSpaces();
    virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void exportTimers(); 
    
    virtual void assemble();
    void assembleSTD();   
    void assembleF();
    void assembleMatrixIBC(int i, std::string markerOpt = "" ); 
    void assembleRhsIBC(int i, std::string marker = "", double intjn = 0);

	void assembleCst();
	void assembleNonCst();

    void geometricTest();
    
    void solve();


    // time step scheme
    virtual void createTimeDiscretization() ;
    newmark_ptrtype timeStepNM() { return M_nm_mixedelasticity; }
    newmark_ptrtype const& timeStepNM() const { return M_nm_mixedelasticity; }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepNM(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepNM(); }
    virtual void updateTimeStepNM();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepNM(); }

};

template<int Dim, int Order, int G_Order, int E_Order>
MixedElasticity<Dim, Order, G_Order, E_Order>::MixedElasticity( std::string const& prefix,
                                                                worldcomm_ptr_t const& worldComm,
                                                                std::string const& subPrefix,
                                                                ModelBaseRepository const& modelRep )
    : super_type( prefix, worldComm, subPrefix, modelRep )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "start",
            this->worldComm(),this->verboseAllProc());

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"MixedElasticity.info") );


    M_prefix = prefix;

    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
	if (boption(prefixvm(this->prefix(), "use-sc")))
	{
    	if ( M_prefix.empty())
        	M_backend = backend( _name="sc", _rebuild=true);
    	else
        	M_backend = backend( _name=prefixvm(prefix,"sc"), _rebuild=true);
	} 
	else	
	{
    	if ( M_prefix.empty())
        	M_backend = backend( _rebuild=true);
    	else
        	M_backend = backend( _name=prefix, _rebuild=true);
	}

    M_tau_constant = doption (prefixvm(M_prefix, "tau_constant") );
    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );
    M_useUserIBC = false;

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticityConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticitySolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

template<int Dim, int Order, int G_Order, int E_Order>
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

template<int Dim, int Order, int G_Order, int E_Order>
void MixedElasticity<Dim, Order, G_Order, E_Order>::initTimeStep()
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
template<int Dim, int Order, int G_Order, int E_Order>
void MixedElasticity<Dim, Order, G_Order, E_Order>::updateTimeStepNM()
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

template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim,Order,G_Order, E_Order>::createTimeDiscretization()
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

template<int Dim, int Order, int G_Order, int E_Order>
typename MixedElasticity<Dim,Order, G_Order, E_Order>::self_ptrtype
MixedElasticity<Dim,Order,G_Order,E_Order>::New( std::string const& prefix,
                                                 worldcomm_ptr_t const& worldComm, std::string const& subPrefix,
                                                 ModelBaseRepository const& modelRep )
{
    return std::make_shared<self_type> ( prefix,worldComm,subPrefix,modelRep );
}

template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::init( mesh_ptrtype mesh, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;
    M_timers["mesh"].push_back(toc("initMesh"));

    tic();
    this->initModel();
    M_timers["initModel"].push_back(toc("initModel"));

    tic();
    this->initSpaces();
    M_timers["spaces"].push_back(toc("initSpaces"));

	if (!isStationary())
	{
		tic();
		this->createTimeDiscretization();	
        this->initTimeStep();
		toc("time_discretization");
	}

    tic();
    this->initExporter(meshVisu);
    M_timers["exporter"].push_back(toc("initExporter"));

    	
    tic();
    this->assemble();
    M_timers["asbMatrix"].push_back(toc("assembleMatrix"));
	

}

template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::initModel()
{

    // initialize marker lists for each boundary condition type
    // Strain
    auto itField = M_modelProperties->boundaryConditions().find( "stress");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
       
 
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
        
        itType = mapField.find( "Neumann_scalar" );
        if ( itType != mapField.end() )
        {
            cout << "Neumann scalar:";
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
	 
        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            cout << "Neumann computed from displacement:";
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

    // Displacement
    itField = M_modelProperties->boundaryConditions().find( "displacement");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        
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
    }


    if ( !M_IBCList.empty() )
    {
        M_IBCList.clear();
    }
    itField = M_modelProperties->boundaryConditions().find( "stress");
    if ( itField != M_modelProperties->boundaryConditions().end() )
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


template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::initSpaces()
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

    auto face_mesh = createSubmesh( M_mesh, complement_integral_bdy, EXTRACTION_KEEP_MESH_RELATION, 0 );


    M_Vh = Pdhms<Order>( M_mesh, true );
    M_Wh = Pdhv<Order>( M_mesh, true );
    M_Mh = Pdhv<Order>( face_mesh, true );
    M_M0h = Pdh<0>( face_mesh );

	std::vector<std::string> ibc_markers(M_integralCondition);
	for( int i = 0; i < M_integralCondition; i++)
	{
		ibc_markers.push_back(M_IBCList[i].marker());
	}

    auto ibc_mesh = createSubmesh( M_mesh, markedfaces(M_mesh, ibc_markers), EXTRACTION_KEEP_MESH_RELATION, 0 );	
    M_Ch = Pchv<0>( ibc_mesh, true ); 
    // M_Ch = Pchv<0>( M_mesh, true );

    auto ibcSpaces = std::make_shared<ProductSpace<Ch_ptr_t,true> >( M_integralCondition, M_Ch);
    M_ps = std::make_shared<product2_space_type>(product2(ibcSpaces,M_Vh,M_Wh,M_Mh));

	// M_ps = std::make_shared<product_space_std>(product(M_Vh,M_Wh,M_Mh));

    M_up = M_Vh->element( "u" ); // Strain
    M_pp = M_Wh->element( "p" ); // Displacement

    cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
    if ( M_integralCondition )
        cout << "Ch<" << 0 << "> : " << M_Ch->nDof() << std::endl;
}

template<int Dim, int Order, int G_Order, int E_Order>
void MixedElasticity<Dim, Order, G_Order, E_Order>::assembleMatrixIBC( int i , std::string markerOpt)
{


    auto bbf = blockform2( *M_ps, M_A_cst);

    auto v = M_Vh->element( "v" );
    auto u = M_Wh->element( "u" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto uI = M_Ch->element( "uI" );
   

    auto H = M_M0h->element( "H" );

    if ( ioption(prefixvm(prefix(), "hface") ) == 0 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(prefix(), "hface") ) == 1 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(prefix(), "hface") ) == 2 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
        H.on( _range=elements(M_M0h->mesh()), _expr=h() );

    // stabilisation parameter
    auto tau_constant = cst(doption(prefixvm(prefix(), "tau_constant")));

    std::string marker;
    if ( !markerOpt.empty())
    {
        marker = markerOpt;
    }
    else
    {
        auto exAtMarker = M_IBCList[i];
        marker = exAtMarker.marker();
        Feel::cout << "Integral on: " << marker << std::endl;
    }

    
    // <lambda, v.n>_Gamma_I
    bbf( 0_c, 3_c, 0, i) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh, marker), _expr=-trans(idt(uI))*(id(v)*N()) );
    
    // <lambda, tau w>_Gamma_I
    bbf( 1_c, 3_c, 1, i ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr= tau_constant * trans(idt(uI)) * pow(idv(H),M_tau_order)*id(w) );
    
    // <sigma.n, m>_Gamma_I
    bbf( 3_c, 0_c, i, 0 ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr= inner(idt(v)*N(),id(nu)) );
     

    // <tau u, m>_Gamma_I
    bbf( 3_c, 1_c, i, 1 ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                        _expr= tau_constant * pow(idv(H),M_tau_order)* inner(idt(u),id(nu)) ),

    // -<lambda2, m>_Gamma_I
    bbf( 3_c, 3_c, i, i ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                        _expr=-tau_constant * pow(idv(H),M_tau_order) * inner(idt(uI),id(nu)) );
    


}

template<int Dim, int Order, int G_Order, int E_Order>
void MixedElasticity<Dim, Order, G_Order, E_Order>::assembleRhsIBC( int i, std::string markerOpt, double intjn )
{
    auto blf = blockform1( *M_ps, M_F );
    auto nu = M_Ch->element( "nu" );

    auto exAtMarker = M_IBCList[i];
    auto marker = exAtMarker.marker();
    
    auto g = expr<Dim,1,expr_order>(exAtMarker.expression());
    if ( !this->isStationary() )
        g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );

    Feel::cout << "IBC condition: " << g << std::endl; 

    double meas = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
    Feel::cout << "Measure of the ibc: " << meas << std::endl;

    // <F_target,m>_Gamma_I
    blf(3_c,i) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=inner(g,id(nu))/meas);

}




template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::initExporter( mesh_ptrtype meshVisu ) 
{
     std::string geoExportType="static"; //change_coords_only, change, static
     M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh(),
                             _name="Export",
                             _geo=geoExportType,
         		    		 _path=this->exporterPath() ); 
}


template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::assemble()
{

    tic();
	solve::strategy s = boption(prefixvm(prefix(), "use-sc"))?solve::strategy::static_condensation:solve::strategy::monolithic;

    auto U = M_ps -> element();
    M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps), *M_backend ); //M_backend->newBlockMatrix(_block=csrGraphBlocks(ps)); 
	//M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s==solve::strategy::static_condensation)?Pattern::COUPLED:pattern), *M_backend, (s==solve::strategy::static_condensation)?false:true);

    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *M_backend, false);//M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    //    M_A_cst = M_backend->newBlockMatrix(_block=csrGraphBlocks(*M_ps));
    //M_F = M_backend->newBlockVector(_block=blockVector(*M_ps), _copy_values=false);
    toc("creating matrices and vectors");
    

}

template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::assembleCst()
{
    // Assembling standard matrix
    tic();
	M_A_cst->zero();
    this->assembleSTD();
    M_timers["asbStd"].push_back(toc("assembleStandardMatrix"));

    // Assembling ibc part
    tic();
    for ( int i = 0; i < M_IBCList.size(); i++ )
        this->assembleMatrixIBC( i );
    M_timers["asbIbc"].push_back(toc("assembleIbcMatrix"));

    M_A_cst->close();
}


template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::assembleNonCst()
{
    tic();
	M_F->zero();
    this->assembleF( );

    for ( int i = 0; i < M_IBCList.size(); i++ )
        this->assembleRhsIBC( i );
    M_F->close();
    M_timers["asbRHS"].push_back(toc("assembleRHS"));

}


template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::solve()
{
    
    auto U = M_ps -> element();
	auto bbf = blockform2(*M_ps, M_A_cst);
    
	auto blf = blockform1(*M_ps, M_F);


	std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(get_backend(),hdgNullSpace(M_Wh,mpl::int_<FEELPP_DIM>())) );
	get_backend()->attachNearNullSpace( myNullSpace );
    if ( boption(_name=prefixvm( this->prefix(), "nullspace").c_str()) )
	    get_backend()->attachNearNullSpace( myNullSpace );


    std::string solver_string = "MixedElasticity : ";
    if( boption(prefixvm(this->prefix(), "use-sc")) )
        solver_string += "static condensation";
    else
        solver_string += "monolithic";
    
    tic();
    tic();
    bbf.solve(_solution=U, _rhs=blf, _rebuild=false, _condense=boption(prefixvm(this->prefix(), "use-sc")), _name= this->prefix());
    M_timers["solver"].push_back(toc("solver"));
    toc(solver_string);
    
    M_up = U(0_c);
    M_pp = U(1_c);

    if ( VLOG_IS_ON(2) )
    {
        cout << "u_hat=" << U(2_c) << std::endl;
        cout << "u=" << U(1_c) << std::endl;
        cout << "sigma=" << U(0_c) << std::endl;
    }

    for( int i = 0; i < M_integralCondition; i++ )
        M_mup.push_back(U(3_c,i));
 
}


template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order, E_Order>::assembleSTD()
{
    auto tau_constant = cst(M_tau_constant); 

    auto sigma = M_Vh->element( "sigma" ); 
    auto v     = M_Vh->element( "v" ); 
    auto u     = M_Wh->element( "u" ); 
    auto w     = M_Wh->element( "w" ); 
    auto uhat  = M_Mh->element( "uhat" ); 
    auto m     = M_Mh->element( "m" ); 
    auto H     = M_M0h->element( "H" ); 

    auto gammaMinusIntegral = complement(boundaryfaces(M_mesh),[this]( auto const& e ) {
        for( auto exAtMarker : this->M_IBCList)
        {
            if ( e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                return true;
        }
        return false; });



    if ( ioption(prefixvm(M_prefix, "hface") ) == 0 )
	    H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 1 )
	    H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 2 )
	    H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
	    H.on( _range=elements(M_M0h->mesh()), _expr=h() );

    auto sc_param = boption(prefixvm(prefix(), "use-sc")) ? 0.5 : 1.0;

    auto bbf = blockform2 ( *M_ps, M_A_cst );


    for( auto const& pairMat : M_modelProperties->materials() )
    {
		auto material = pairMat.second;
		auto lambda = material.getScalar("lambda");
		Feel::cout << "Lambda: " << lambda << std::endl;
		auto mu = material.getScalar("mu");
		Feel::cout << "Mu: " << mu << std::endl;
		auto c1 = cst(0.5)/mu; 
    	auto c2 = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu)); 
        Feel::cout << "c1: " << mean(_range=elements(M_mesh),_expr=c1) << std::endl;
        Feel::cout << "c2: " << mean(_range=elements(M_mesh),_expr=c2) << std::endl;

		bbf( 0_c, 0_c ) +=  integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    	bbf( 0_c, 0_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );
    }


    bbf( 0_c, 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(trans(idt(u))*div(v)));

    bbf( 0_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
							    _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
			    						 trans(idt(uhat))*rightface(id(v)*N())) );
    bbf( 0_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=gammaMinusIntegral,
		    					_expr=-trans(idt(uhat))*(id(v)*N()) );

    bbf( 1_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
								_expr=(trans(id(w))*divt(sigma)));

    // ( d^2u/dt^2, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() ) {
		auto dt = this->timeStep();
    	for( auto const& pairMat : M_modelProperties->materials() )
    	{
			auto material = pairMat.second;
			auto rho = material.getScalar("rho");
			// bbf( 1_c, 1_c ) += integrate(_range=elements(M_mesh),
        	//                              _expr = this->timeStepNM()->polySecondDerivCoefficient()*rho*inner(idt(u),id(w)) );
			bbf( 1_c, 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
            	                         _expr = -rho*inner(idt(u),id(w))/(dt*dt) );
		}
    }

    // begin dp: here we need to put the projection of u on the faces
    bbf( 1_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),_expr=-tau_constant *
		    ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
		      rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));

    bbf( 1_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=boundaryfaces(M_mesh),
		    _expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));

    bbf( 1_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
		    _expr=tau_constant *
		    ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
		      rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));

    bbf( 1_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=gammaMinusIntegral,
		    _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );


    bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
		    _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
				    rightfacet(idt(sigma)*N())) ) );


    // BC
    bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
			    _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
				    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));

    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
    _expr=sc_param*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
				    rightface( pow(idv(H),M_tau_order) )));

    auto itField = M_modelProperties->boundaryConditions().find( "displacement");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
	    auto mapField = (*itField).second;
	    auto itType = mapField.find( "Dirichlet" );
	    if ( itType != mapField.end() )
	    {
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();
		    	cout << "Dirichlet on " << marker << std::endl;
			    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
					    					_expr=trans(idt(uhat)) * id(m) );
		    }
	    }
	
    }
    
	itField = M_modelProperties->boundaryConditions().find( "stress");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {	
	    auto mapField = (*itField).second;
	    auto itType = mapField.find( "Neumann" );
	    if ( itType != mapField.end() )
	    {
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();
		    	cout << "Neumann on " << marker << std::endl;
			    bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker), 
					    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );
		    }
	    }
	    itType = mapField.find( "Neumann_scalar" );
	    if ( itType != mapField.end() )
	    {
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();
		    	cout << "Neumann on " << marker << std::endl;
			    bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker), 
					    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );
		    }
	    }
	    itType = mapField.find( "Neumann_exact" );
	    if ( itType != mapField.end() )
	    {
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();
		    	cout << "Neumann on " << marker << std::endl;
			    bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker), 
					    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );
		    }
	    }
    }

} // end assemble STD
  

template<int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim, Order, G_Order,E_Order>::assembleF()
{

    auto blf = blockform1( *M_ps, M_F );

    auto w     = M_Wh->element( "w" ); 
    auto m     = M_Mh->element( "m" ); 
    
    // Building the RHS

    auto itField = M_modelProperties->boundaryConditions().find("stress");
    if (itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("SourceTerm");
        
		if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto g = expr<Dim,1,expr_order> (exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
				blf( 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
                    			  	      _expr=trans(g)*id(w));
            }
        }
		itType = mapField.find("Neumann");
		if ( itType != mapField.end() )
		{
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto marker = exAtMarker.marker();
                auto g = expr<Dim,1,expr_order> (exAtMarker.expression());
				if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
			    cout << "Neumann condition on " << marker << ": " << g << std::endl;
			    blf( 2_c ) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
					    				_expr=trans(id(m))* g );
            }
		}
        	
		itType = mapField.find("Neumann_scalar");
		if ( itType != mapField.end() )
		{
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto marker = exAtMarker.marker();
                auto g = expr<expr_order> (exAtMarker.expression());
				if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
			    cout << "Neumann condition on " << marker << ": " << g << std::endl;
			    blf( 2_c ) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
				    					_expr= inner(expr(g)*N(), id(m)) );
            }
		}	
	    
        itType = mapField.find("Neumann_exact");
	    if ( itType != mapField.end() )
	    {
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    auto marker = exAtMarker.marker();
			    auto g = expr<Dim,1, expr_order>(exAtMarker.expression());
			    if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
                 
			    for( auto const& pairMat : M_modelProperties->materials() )
			    {
				    auto gradu_exact = grad( g );
				    auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
			        auto material = pairMat.second;
				    auto lambda = material.getScalar("lambda");
				    auto mu = material.getScalar("mu");
				    auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;

	    			cout << "Neumann condition computed from displacement on " << marker << std::endl;
		    		blf( 2_c ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(id(m)) * sigma_exact *N() );
                }
            }
	    }
    } 
    
    itField = M_modelProperties->boundaryConditions().find("displacement");
    if (itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("Dirichlet");
		if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto marker = exAtMarker.marker();
                auto g = expr<Dim,1, expr_order>(exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
				cout << "Dirichlet condition on " << marker << ": " << g << std::endl;
				blf( 2_c ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(id(m))*g);
            }
        }
    }


    // (u_old,w)_Omega
    if ( !this->isStationary() )
    {
    	for( auto const& pairMat : M_modelProperties->materials() )
    	{
			auto material = pairMat.second;
			auto rho = material.getScalar("rho");
			auto u = this->timeStepNM()->previousUnknown(0);
			auto u1 = this->timeStepNM()->previousUnknown(1);
			auto dt = this-> timeStep();
        	// blf(1_c) += integrate( _range=elements(M_mesh),
        	//                         _expr= rho*inner(idv(this->timeStepNM()->polyDeriv()),id(w)) );
        	blf(1_c) += integrate(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr= -rho*inner( 2*idv(u)-idv(u1) ,id(w))/(dt*dt) );
		}
    }

} // end assembleF



template <int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim,Order, G_Order, E_Order>::exportResults( double time, mesh_ptrtype mesh , op_interp_ptrtype Idh , opv_interp_ptrtype Idhv )
{
    this->log("MixedElasticity","exportResults", "start");
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
     {
         for ( auto const& field : M_modelProperties->postProcess().exports().fields() )
         {
            if ( field == "stress" )
            {
                LOG(INFO) << "exporting stress at time " << time;
		
                // Exporting the stress component by component
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::X)) );
                M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaXX"), SXX );

                if (Dim > 1)
                {
                    auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Y, Component::Y)) );
                    auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::Y)) );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaYY"), SYY );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaXY"), SXY );
                }
                if (Dim > 2)
                {            
                    auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Z, Component::Z)) );
                    auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Y, Component::Z)) );
                    auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::Z)) );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaZZ"), SZZ );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaYZ"), SYZ );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"sigmaXZ"), SXZ );
                }

         		// M_exporter->step(time)->add(prefixvm(M_prefix, "stress"), Idhv?(*Idhv)( M_up):M_up );
         
               if (M_integralCondition)
                {

                    for( auto exAtMarker : this->M_IBCList)
                    {
                        std::vector<double> force_integral(Dim);
                        auto marker = exAtMarker.marker();
                        LOG(INFO) << "exporting integral flux at time "
                                  << time << " on marker " << marker;
                        auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                            _expr=trans(idv(M_up))*N());
                    
                        Feel::cout << "Force computed: " << std::endl;
                        for( auto i=0;i < Dim;i++ )
                        {
                            std::string stringForce_help = (boost::format("integralForce_%1%")%i).str();
                            force_integral[i] = j_integral.evaluate()(i,0);
                            Feel::cout << force_integral[i] << std::endl;
                            M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                        }
                    }

                }

            }
            else if ( M_mesh->hasFaceMarker(field) )
            { 
                auto marker = field;
            	LOG(INFO) << "exporting computed force on " << marker << " at time " << time;
                std::vector<double> force_integral(Dim);

                auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                            _expr=trans(idv(M_up))*N());
                    
                Feel::cout << "Force computed: " << std::endl;
                for( auto i=0;i < Dim;i++ )
                {
                    std::string stringForce_help = (boost::format("integralForce_%1%")%i).str();
                    force_integral[i] = j_integral.evaluate()(i,0);
                    Feel::cout << force_integral[i] << std::endl;
                    M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                }
            }
            else if ( field == "scaled_displacement" )
            {
                auto scaled_displ = M_Wh->element("scaled_displacement");
                for( auto const& pairMat : modelProperties()->materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.getScalar( "scale_displacement" );

                    scaled_displ.on( _range=markedelements(M_mesh,marker) , _expr= kk*idv(M_pp));
                }

                M_exporter->step(time)->add(prefixvm(M_prefix, "displacement"),Idh?(*Idh)( scaled_displ):scaled_displ ) ;    	



            }
            else if ( field == "scaled_stress" )
            {
                auto scaled_stress = M_Vh->element("scaled_stress");
                for( auto const& pairMat : modelProperties()->materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.getScalar( "scale_stress" );

                    scaled_stress.on( _range=markedelements(M_mesh,marker) , _expr= kk*idv(M_up));
                }
                
                // Exporting the scaled stress component by component
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::X)) );
                M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaXX"), SXX );

                if (Dim > 1)
                {
                    auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Y, Component::Y)) );
                    auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::Y)) );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaYY"), SYY );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaXY"), SXY );
                }
                if (Dim > 2)
                {            
                    auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Z, Component::Z)) );
                    auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Y, Component::Z)) );
                    auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::Z)) );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaZZ"), SZZ );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaYZ"), SYZ );
                    M_exporter->step(time)->add(prefixvm(M_prefix,"scaled_sigmaXZ"), SXZ );
                }
               
                // Exporting scaled stress integral
                if (M_integralCondition)
                {
                    for( auto exAtMarker : this->M_IBCList)
                    {
                        std::vector<double> force_integral(Dim);
                        auto marker = exAtMarker.marker();
                        LOG(INFO) << "exporting scaled integral flux at time "
                                  << time << " on marker " << marker;
                        auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(idv(scaled_stress))*N());
                    
                        Feel::cout << "Force computed: " << std::endl;
                        for( auto i=0;i < Dim;i++ )
                        {
                            std::string stringForce_help = (boost::format("scaled_integralForce_%1%")%i).str();
                            force_integral[i] = j_integral.evaluate()(i,0);
                            Feel::cout << force_integral[i] << std::endl;
                            M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                        }
                    }

                }


            }
            else if ( field == "displacement" )
        	{
            	LOG(INFO) << "exporting displacement at time " << time;
                M_exporter->step(time)->add(prefixvm(M_prefix, "displacement"),Idh?(*Idh)( M_pp):M_pp ) ;    	
                // Projecting on L2 space for continuity.
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto UX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::X]) );     
                M_exporter->step(time)->add(prefixvm(M_prefix, "UX"),UX ) ;  
  	
                if (Dim > 1)
                {
                    auto UY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::Y]) );
                    M_exporter->step(time)->add(prefixvm(M_prefix, "UY"),UY ) ;    	
                }
                if (Dim > 2)
                {   
                    auto UZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::Z]) );
                    M_exporter->step(time)->add(prefixvm(M_prefix, "UZ"),UZ ) ;    	
	            }

				auto itField = M_modelProperties->boundaryConditions().find("ExactSolution");
 				if ( itField != M_modelProperties->boundaryConditions().end() )
 				{
 		    		auto mapField = (*itField).second;
 		    		auto itType = mapField.find( "u_exact" );
 				    if (itType != mapField.end() )
 				    {
 						for (auto const& exAtMarker : (*itType).second )
 						{
    		    			if (exAtMarker.isExpression() )
 			    			{		
								auto u_exact = expr<Dim,1,expr_order> (exAtMarker.expression());
								if ( !this->isStationary() )
								    u_exact.setParameterValues( { {"t", time } } );

								auto export_uEX = project(_quad=_Q<expr_order>(), _space=M_Wh, _range=elements( M_mesh ), _expr=u_exact);
								M_exporter->step(time)->add(prefixvm(M_prefix, "u_exact"), Idh?(*Idh)( export_uEX): export_uEX );		
                                
								auto l2err_u = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=u_exact - idv(M_pp) );
							
                                auto l2norm_uex = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=u_exact );
								if (l2norm_uex < 1)
									l2norm_uex = 1.0;	
                                
								cout << "----- Computed Errors -----" << std::endl;
								// cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
								cout << "||u-u_ex||_L2=\t" << l2err_u << std::endl;
								// Export the errors
								M_exporter -> step( time )->add(prefixvm(M_prefix, "u_error_L2"), l2err_u/l2norm_uex );
							    //------ Sigma 	------//
								auto gradu_exact = grad(u_exact);
								auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
								for( auto const& pairMat : M_modelProperties->materials() )
								{
									auto material = pairMat.second;
									auto lambda = material.getScalar("lambda");
									auto mu = material.getScalar("mu");
									auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;

									// EXPORT SIGMA EXACT
                                    auto export_sigmaEX = project(_quad=_Q<expr_order>(), _space=M_Vh, _range=elements(M_mesh), _expr=sigma_exact); 
				
                                    auto Sh = Pch<Order> (M_mesh);
                                    auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);
                                    
                                    auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::X)) );
                                    M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactXX"), SXX );
                                    
                                    if (Dim > 1)
                                    {
                                        auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Y, Component::Y)) );
                                        auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::Y)) );
                                        M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactYY"), SYY );
                                        M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactXY"), SXY );
                                    }
                                    if (Dim > 2)
                                    {
                                        auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Y, Component::Z)) );
                                        auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::Z)) );
                                        auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Z, Component::Z)) );
                                        M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactZZ"), SZZ );
                                        M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactYZ"), SYZ );
                                        M_exporter->step(time)->add(prefixvm(M_prefix,"s_exactXZ"), SXZ );
                                    }

                					// M_exporter->add(prefixvm(M_prefix, "sigma_exact"), Idhv?(*Idhv)( export_sigmaEX): export_sigmaEX );

									auto l2err_sigma = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=sigma_exact - idv(M_up) );
									auto l2norm_sigmaex = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=sigma_exact );
									if (l2norm_sigmaex < 1)
										l2norm_sigmaex = 1.0;		
									cout << "||sigma-sigma_ex||_L2=\t" << l2err_sigma/l2norm_sigmaex << std::endl;
									cout << "---------------------------" << std::endl;
									// Export the errors
									M_exporter -> step( time )->add(prefixvm(M_prefix, "sigma_error_L2"), l2err_sigma/l2norm_sigmaex );
								} 
 							}
 		    			}
 					}
				}
             }
         }
     }

     this->timerTool("PostProcessing").stop("exportResults");
     if ( this->scalabilitySave() )
     {
         if ( !this->isStationary() )
             this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
         this->timerTool("PostProcessing").save();
     }
     this->log("MixedElasticity","exportResults", "finish");
}

template <int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim,Order, G_Order, E_Order>::geometricTest( ){

    auto itField = M_modelProperties -> boundaryConditions().find("GeometricalTest");
    if ( itField != M_modelProperties -> boundaryConditions().end() )
    {
        auto mapField = itField -> second;
        auto itType = mapField.find( "force_F" );
        if (itType != mapField.end() )
        {
            for (auto const& exAtMarker : itType->second )
            {
                auto curvedForce = (integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,exAtMarker.marker()), _expr = idv(M_up)*N() )).evaluate();
                auto forceF = (expr<Dim,1,expr_order> (exAtMarker.expression() )).evaluate();
                /* 
                Feel::cout << "Force F uploaded:\t" << forceF << std::endl;
                Feel::cout << "Force F computed from M_up:\t" << curvedForce << std::endl; 
                */
                auto curveError = (curvedForce - forceF).cwiseAbs();
    
                Feel::cout << "Error for geometrical order:\t" << curveError << std::endl;  
            }
        }
        /* 
        itType = mapField.find( "force_F_2" );
        if (itType != mapField.end() )
        {
            for (auto const& exAtMarker : itType->second )
            {
                auto forceF_2 = expr<Dim,1,expr_order> (exAtMarker.expression() );
                auto forceIntegral = (integrate( _range=markedfaces(M_mesh,exAtMarker.marker()), _expr = forceF_2 * N() )).evaluate();
                Feel::cout << "Force F computed from input:\t" << forceIntegral << std::endl; 
    
            }
        }
        */
    }
}


// Time exporter
template <int Dim, int Order, int G_Order, int E_Order>
void
MixedElasticity<Dim,Order, G_Order, E_Order>::exportTimers()
{
    if( Environment::isMasterRank() )
    {
        std::ofstream timers( "timers.dat", std::ios::out | std::ios::trunc);
        std::string fmtS = "";
        for( int i = 1; i <= M_timers.size(); ++i )
            fmtS += "%" + std::to_string(i) + "% %|" + std::to_string(14*i) + "t|";
        boost::format fmt(fmtS);
        for( auto const& pair : M_timers )
            fmt % pair.first;
        timers << fmt << std::endl;
        if (isStationary())
        {
                for( auto const& pair : M_timers )
                    fmt % pair.second;
                timers << fmt << std::endl;
        }

        /* 
           //( for( int i = 0; this->timeStepBase()->isFinished() ; ++i )
            {
                for( auto const& pair : M_timers )
                    fmt % pair.second[i];
                timers << fmt << std::endl;
        }*/
        
        timers.close();
    }
}

} // Namespace FeelModels

} // Namespace Feel

#endif
