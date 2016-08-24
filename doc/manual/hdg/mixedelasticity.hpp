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

namespace Feel {

namespace FeelModels {

inline
po::options_description
makeMixedElasticityOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpOptions( "Mixed Elasticity HDG options");
    mpOptions.add_options()
        ( "gmsh.submesh", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( "gmsh.submesh2", po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ( prefixvm( prefix,"tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix,"tau_order").c_str(), po::value<int>()->default_value( -1 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ;
    mpOptions.add ( envfeelmodels_options( prefix ) ).add( modelnumerical_options( prefix ) );
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

template<int Dim, int Order, int G_Order = 1>
class MixedElasticity    :	public ModelNumerical
{
public:
    typedef ModelNumerical super_type;

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedElasticity<Dim,Order,G_Order> self_type;
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
    
	using op_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Wh_t, Pdhv_type<mesh_type,Order>>>;
    using opv_interp_ptrtype = boost::shared_ptr<OperatorInterpolation<Vh_t, Pdhms_type<mesh_type,Order>>>;
 
    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    using product_space_std = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using bilinear_block_std = BlockBilinearForm<product_space_std>;

    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef boost::shared_ptr <exporter_type> exporter_ptrtype;
    
    // typedef Newmark<space_mixedelasticity_type>  newmark_type;
    typedef Newmark <Wh_t> newmark_type;
    typedef boost::shared_ptr<newmark_type> newmark_ptrtype;
    
//private:
protected:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh; // potential trace

    backend_ptrtype M_backend;
    sparse_matrix_ptrtype M_A_cst;
    vector_ptrtype M_F;

    Vh_element_t M_up; // flux solution
    Wh_element_t M_pp; // potential solution 

    double M_tau_constant;
    int M_tau_order;
    
    // time discretization
    newmark_ptrtype M_nm_mixedelasticity;

public:
    
    // constructor
   MixedElasticity( std::string const& prefix = "mixedelasticity",                   
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    
    MixedElasticity( self_type const& ME ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedelasticity",
                             WorldComm const& worldComm = Environment::worldComm(),
                             std::string const& subPrefix = "",
                             std::string const& rootRepository = ModelBase::rootRepositoryByDefault() ); 

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }

    Vh_element_t fluxField() const { return M_up; }
    Wh_element_t potentialField() const { return M_pp; }
    model_prop_type modelProperties() { return *M_modelProperties; }
    model_prop_type modelProperties() const { return *M_modelProperties; }
    
    int tau_order() const { return M_tau_order; }
    backend_ptrtype get_backend() { return M_backend; }
    vector_ptrtype getF() {return M_F; }
 
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
    virtual void assemble();
    
    template<typename PS>
    void assembleSTD(PS&& ps);
    template<typename PS>
    void assembleF(PS&& ps);
    void solve();

    // time step scheme
    virtual void createTimeDiscretization() ;
    newmark_ptrtype timeStepNM() { return M_nm_mixedelasticity; }
    newmark_ptrtype const& timeStepNM() const { return M_nm_mixedelasticity; }
    boost::shared_ptr<TSBase> timeStepBase() { return this->timeStepNM(); }
    boost::shared_ptr<TSBase> timeStepBase() const { return this->timeStepNM(); }
    virtual void updateTimeStepNM();
    virtual void initTimeStep();
    void updateTimeStep() { this->updateTimeStepNM(); }

};

template<int Dim, int Order, int G_Order>
MixedElasticity<Dim, Order, G_Order>::MixedElasticity( std::string const& prefix,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    : super_type( prefix, worldComm, subPrefix, rootRepository )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"MixedElasticity.info") );

    M_prefix = prefix;
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_tau_constant = doption (prefixvm(M_prefix, "tau_constant") );
    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticityConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticitySolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


template<int Dim, int Order, int G_Order>
void MixedElasticity<Dim, Order, G_Order>::initTimeStep()
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
template<int Dim, int Order, int G_Order>
void MixedElasticity<Dim, Order, G_Order>::updateTimeStepNM()
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

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim,Order,G_Order>::createTimeDiscretization()
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

template<int Dim, int Order, int G_Order>
typename MixedElasticity<Dim,Order, G_Order>::self_ptrtype
MixedElasticity<Dim,Order,G_Order>::New( std::string const& prefix,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type> ( prefix,worldComm,subPrefix,rootRepository );
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::init( mesh_ptrtype mesh, mesh_ptrtype meshVisu )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;
    toc("mesh");

    tic();
    this->initModel();
    toc("model");

    tic();
    this->initSpaces();
    toc("spaces");

	if (!isStationary())
	{
		tic();
		this->createTimeDiscretization();	
        this->initTimeStep();
		toc("time_discretization");
	}

    tic();
    this->initExporter(meshVisu);
    toc("exporter");

    tic();
    this->assemble();
    toc("assemble");


}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::initModel()
{


    // initialize marker lists for each boundary condition type
    // Strain
    auto itField = M_modelProperties->boundaryConditions().find( "strain");
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

}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::initSpaces()
{


    auto face_mesh = createSubmesh( M_mesh, faces(M_mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );

    M_Vh = Pdhms<Order>( M_mesh, true );
    M_Wh = Pdhv<Order>( M_mesh, true );
    M_Mh = Pdhv<Order>( face_mesh, true );

    M_up = M_Vh->element( "u" ); // Strain
    M_pp = M_Wh->element( "p" ); // Displacement

    cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::initExporter( mesh_ptrtype meshVisu ) 
{
     std::string geoExportType="static"; //change_coords_only, change, static
     M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh(),
                             _name="Export",
                             _geo=geoExportType,
         		    		 _path=this->exporterPath() ); 
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::assemble()
{

    auto ps = product(M_Vh,M_Wh,M_Mh);

    tic();
    auto U = ps.element();
    M_A_cst = M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    M_F = M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    toc("creating matrices and vectors");
    
    tic();
    this->assembleSTD( ps );
    M_A_cst->close();
    toc("assemble A");
    /*
    tic();
    this->assembleF( ps );
    M_F->close();
    toc("assemble F");
    */
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::solve()
{
    
    auto ps = product(M_Vh,M_Wh,M_Mh);

    tic();
	// M_F->clear();
    this->assembleF( ps );
    M_F->close();
    toc("assemble F");
    
	tic();
    auto U = ps.element();
	auto M_U = M_backend->newBlockVector(_block=U,_copy_values= false);
	M_backend->solve(_matrix=M_A_cst, _rhs = M_F , _solution=M_U) ;
	U.localize(M_U);	
    M_up = U(0_c);
    M_pp = U(1_c);
   
    toc("solve");
}

template<int Dim, int Order, int G_Order>
template<typename PS>
void
MixedElasticity<Dim, Order, G_Order>::assembleSTD(PS&& ps)
{
    auto tau_constant = cst(M_tau_constant); 
    auto face_mesh = M_Mh->mesh();//createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 ); 
    M0h_ptr_t M0h = Pdh<0>( face_mesh,true ); 

    auto sigma = M_Vh->element( "sigma" ); 
    auto v     = M_Vh->element( "v" ); 
    auto u     = M_Wh->element( "u" ); 
    auto w     = M_Wh->element( "w" ); 
    auto uhat  = M_Mh->element( "uhat" ); 
    auto m     = M_Mh->element( "m" ); 
    auto H     = M0h->element( "H" ); 

    if ( ioption(prefixvm(M_prefix, "hface") ) == 0 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 1 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 2 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
	    H.on( _range=elements(M0h->mesh()), _expr=h() );

    auto bbf = blockform2 ( ps, M_A_cst );


    for( auto const& pairMat : M_modelProperties->materials() )
    {
		auto material = pairMat.second;
		auto lambda = expr(material.getScalar("lambda"));
		Feel::cout << "Lambda: " << lambda << std::endl;
		auto mu = expr(material.getScalar("mu"));
		auto c1 = cst(0.5)/mu; 
    	auto c2 = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu)); 
		bbf( 0_c, 0_c ) +=  integrate(_range=elements(M_mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    	bbf( 0_c, 0_c ) += integrate(_range=elements(M_mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );
    }


    bbf( 0_c, 1_c ) += integrate(_range=elements(M_mesh),_expr=(trans(idt(u))*div(v)));

    bbf( 0_c, 2_c) += integrate(_range=internalfaces(M_mesh),
							    _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
			    						 trans(idt(uhat))*rightface(id(v)*N())) );
    bbf( 0_c, 2_c) += integrate(_range=boundaryfaces(M_mesh),
		    					_expr=-trans(idt(uhat))*(id(v)*N()));

    bbf( 1_c, 0_c) += integrate(_range=elements(M_mesh),
								_expr=(trans(id(w))*divt(sigma)));

    // ( d^2u/dt^2, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() ) {
		auto dt = this->timeStep();
    	for( auto const& pairMat : M_modelProperties->materials() )
    	{
			auto material = pairMat.second;
			auto rho = expr(material.getScalar("rho"));
			// bbf( 1_c, 1_c ) += integrate(_range=elements(M_mesh),
        	//                              _expr = this->timeStepNM()->polySecondDerivCoefficient()*rho*inner(idt(u),id(w)) );
			bbf( 1_c, 1_c ) += integrate(_range=elements(M_mesh),
            	                         _expr = rho*inner(idt(u),id(w))/(dt*dt) );
		}
    }

    // begin dp: here we need to put the projection of u on the faces
    bbf( 1_c, 1_c) += integrate(_range=internalfaces(M_mesh),_expr=-tau_constant *
		    ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
		      rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));

    bbf( 1_c, 1_c) += integrate(_range=boundaryfaces(M_mesh),
		    _expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));

    bbf( 1_c, 2_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=tau_constant *
		    ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
		      rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));

    bbf( 1_c, 2_c) += integrate(_range=boundaryfaces(M_mesh),
		    _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );


    bbf( 2_c, 0_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
				    rightfacet(idt(sigma)*N())) ) );


    // BC
    bbf( 2_c, 1_c) += integrate(_range=internalfaces(M_mesh),
			    _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
				    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));

    bbf( 2_c, 2_c) += integrate(_range=internalfaces(M_mesh),
    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
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
			    bbf( 2_c, 2_c) += integrate(_range=markedfaces(M_mesh,marker),
					    					_expr=trans(idt(uhat)) * id(m) );
		    }
	    }
	
    }
    
	itField = M_modelProperties->boundaryConditions().find( "strain");
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
			    bbf( 2_c, 0_c) += integrate(_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_range=markedfaces(M_mesh,marker), 
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
			    bbf( 2_c, 0_c) += integrate(_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_range=markedfaces(M_mesh,marker), 
					    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );
		    }
	    }
    }

} // end assemble STD
  

template<int Dim, int Order, int G_Order>
template< typename PS>
void
MixedElasticity<Dim, Order, G_Order>::assembleF(PS&& ps)
{
    M_F->zero();
    auto blf = blockform1( ps, M_F );

    auto w     = M_Wh->element( "w" ); 
    auto m     = M_Mh->element( "m" ); 
    
    // Building the RHS

    auto itField = M_modelProperties->boundaryConditions().find("strain");
    if (itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("SourceTerm");
        
		if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto g = expr<Dim,G_Order> (exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
				blf( 1_c ) += integrate(_range=elements(M_mesh),
                    			  	      _expr=trans(g)*id(w));
            }
        }
		itType = mapField.find("Neumann");
		if ( itType != mapField.end() )
		{
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto marker = exAtMarker.marker();
                auto g = expr<Dim,Dim> (exAtMarker.expression());
				if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
				cout << "Neumann condition on " << marker << ": " << g << std::endl;
				blf( 2_c ) += integrate(_range=markedfaces(M_mesh,marker),
									    _expr=trans(id(m))*(g*N()));
            }
		}	
		itType = mapField.find("Neumann_scalar");
		if ( itType != mapField.end() )
		{
            for ( auto const& exAtMarker : (*itType).second )
            {
				auto marker = exAtMarker.marker();
				auto g = expr(exAtMarker.expression());
				if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
                // auto g = expr<Dim,Dim> ( scalar*eye<Dim,Dim>() );
				cout << "Neumann scalar condition on " << marker << ": " << g << std::endl;
				blf( 2_c ) += integrate(_range=markedfaces(M_mesh,marker),
									    _expr=trans(id(m))*(eye<Dim,Dim>()*g*N()));
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
                auto g = expr<Dim,G_Order>(exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
				cout << "Dirichlet condition on " << marker << ": " << g << std::endl;
				blf( 2_c ) += integrate(_range=markedfaces(M_mesh,marker),
                		      _expr=trans(id(m))*g);
            }
        }
    }

    // (u_old,w)_Omega
    if ( !this->isStationary() )
    {
    	for( auto const& pairMat : M_modelProperties->materials() )
    	{
			auto material = pairMat.second;
			auto rho = expr(material.getScalar("rho"));
			auto u = this->timeStepNM()->previousUnknown(0);
			auto u1 = this->timeStepNM()->previousUnknown(1);
			auto dt = this-> timeStep();
        	// blf(1_c) += integrate( _range=elements(M_mesh),
        	//                         _expr= rho*inner(idv(this->timeStepNM()->polyDeriv()),id(w)) );
        	blf(1_c) += integrate( _range=elements(M_mesh), _expr= rho*inner( 2*idv(u)-idv(u1) ,id(w))/(dt*dt) );
		}
    }

} // end assembleF


template <int Dim, int Order, int G_Order>
void
MixedElasticity<Dim,Order, G_Order>::exportResults( double time, mesh_ptrtype mesh , op_interp_ptrtype Idh , opv_interp_ptrtype Idhv )
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
     auto postProcess = M_modelProperties->postProcess();
     auto itField = postProcess.find( "Fields");
     
	 if ( itField != postProcess.end() )
     {
         for ( auto const& field : (*itField).second )
         {
            if ( field == "strain" )
            {
                LOG(INFO) << "exporting strain at time " << time;
		 		M_exporter->step(time)->add(prefixvm(M_prefix, "strain"), Idhv?(*Idhv)( M_up):M_up );
            }
            else if ( field == "displacement" )
        	{
            	LOG(INFO) << "exporting displacement at time " << time;
                M_exporter->step(time)->add(prefixvm(M_prefix, "displacement"),Idh?(*Idh)( M_pp):M_pp ) ;    	
		
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
								auto u_exact = expr<Dim,1> (exAtMarker.expression());
								if ( !this->isStationary() )
								    u_exact.setParameterValues( { {"t", time } } );

								auto export_uEX = project( _space=M_Wh, _range=elements( M_mesh ), _expr=u_exact);
								M_exporter->step(time)->add(prefixvm(M_prefix, "u_exact"), Idh?(*Idh)( export_uEX): export_uEX );		

								auto l2err_u = normL2( _range=elements(M_mesh), _expr=u_exact - idv(M_pp) );
								auto l2norm_uex = normL2( _range=elements(M_mesh), _expr=u_exact );
								if (l2norm_uex < 1)
									l2norm_uex = 1.0;	
								cout << "----- Computed Errors -----" << std::endl;
								cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
								// Export the errors
								M_exporter -> step( time )->add(prefixvm(M_prefix, "u_error_L2"), l2err_u/l2norm_uex );
							    //------ Sigma 	------//
								auto gradu_exact = grad(u_exact);
								auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
								for( auto const& pairMat : M_modelProperties->materials() )
								{
									auto material = pairMat.second;
									auto lambda = expr(material.getScalar("lambda"));
									auto mu = expr(material.getScalar("mu"));
									auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;

									auto export_sigmaEX = project( _space=M_Vh, _range=elements(M_mesh), _expr=sigma_exact); 
									M_exporter->add(prefixvm(M_prefix, "sigma_exact"), Idhv?(*Idhv)( export_sigmaEX): export_sigmaEX );

									auto l2err_sigma = normL2( _range=elements(M_mesh), _expr=sigma_exact - idv(M_up) );
									auto l2norm_sigmaex = normL2( _range=elements(M_mesh), _expr=sigma_exact );
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

} // Namespace FeelModels

} // Namespace Feel

#endif
