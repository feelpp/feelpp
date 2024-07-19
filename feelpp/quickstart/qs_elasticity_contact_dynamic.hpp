#pragma once
#include "qs_elasticity_contact.hpp"

namespace Feel
{



template <int Dim, int Order, int OrderGeo>
class ContactDynamic
{
public:
    using mesh_t = Mesh<Simplex<Dim,OrderGeo>>;
    using mesh_tP1 = Mesh<Simplex<Dim>>;
    using spacev_t = Pchv_type<mesh_t, Order>;
    using space_t = Pch_type<mesh_t, Order>;
    using spacev_ptr_t = Pchv_ptrtype<mesh_t, Order>; 
    using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
    using elementv_t = typename spacev_t::element_type;
    using element_t = typename space_t::element_type;
    using form2_type = form2_t<spacev_t,spacev_t>; 
    using form1_type = form1_t<spacev_t>; 
    using ts_ptrtype = std::shared_ptr<NewmarkContact<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_tP1>>; 

    // Constructors
    ContactDynamic() = default;
    ContactDynamic(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    elementv_t const& u() const { return u_; }
    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }

    // Methods
    void initialize();
    void initializeContact();
    void processLoading(form1_type& l);
    void processMaterials(form2_type &a);
    void processBoundaryConditions(form1_type& l, form2_type& a);
    void processContactPenalty(form1_type& l, form2_type& a, Range<mesh_t,MESH_FACES> const& elts, elementv_t const& u);
    void processContactNitsche(form1_type& l, form2_type& a, Range<mesh_t,MESH_FACES> const& elts, elementv_t const& u);
    void processContactPersistency(form1_type& l, form2_type& a, Range<mesh_t,MESH_FACES> const& elts, elementv_t const& u);
    void run();
    Range<mesh_t, MESH_FACES> getContactRegion( elementv_t const& u );
    void timeLoop();
    void timeLoopFixedPoint();
    void exportResults(double t);
    void writeResultsToFile(const std::string& filename) const;
    void initG();
    void storeData();
    void error();

private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    std::shared_ptr<mesh_tP1> meshP1;
    spacev_ptr_t Xhv_;
    space_ptr_t Xh_;

    elementv_t u_;
    element_t contactRegion_;
    element_t contactPressure_;
    element_t contactDisplacement_;
    element_t g_;
    int nbrFaces_;
    Range<mesh_t, MESH_FACES> myelts_;

    ts_ptrtype ts_;
    exporter_ptrtype e_;
    nl::json meas_;

    int quad_,quad1_;

    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    int withMarker_,save_, computeerror_;
    double fixedPointtol_;
    int fixedPoint_; 
    std::vector<double> pressurePoint_;
};

// Constructor
template <int Dim, int Order, int OrderGeo>
ContactDynamic<Dim, Order, OrderGeo>::ContactDynamic(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // Define Xhv
    Xhv_ = Pchv<Order>(mesh_); 

    // Quadratures
    std::string matQuad = fmt::format( "/Models/LinearElasticity/quad");
    quad_ = specs_[nl::json::json_pointer( matQuad )].get<int>();

    std::string matQuad1 = fmt::format( "/Models/LinearElasticity/quad1");
    quad1_ = specs_[nl::json::json_pointer( matQuad1 )].get<int>();

    // Get elastic structure parameters
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
    
    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
    
    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));
    
    // Initialize external forces
    if ( specs_["/Models/LinearElasticity"_json_pointer].contains("loading") )
    {
        for ( auto [key, loading] : specs_["/Models/LinearElasticity/loading"_json_pointer].items() )
        {
            std::string loadtype = fmt::format( "/Models/LinearElasticity/loading/{}/type", key );

            if ( specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Gravity" )
            {
                LOG( INFO ) << fmt::format( "Loading {}: Gravity found", key );
                std::string loadexpr = fmt::format( "/Models/LinearElasticity/loading/{}/parameters/expr", key );
                externalforce_ = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
            }
        }
    }

    // Initialize exporter

    meshP1 = loadMesh( _mesh = new mesh_tP1, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    e_ = Feel::exporter(_mesh = meshP1, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    

    // Initialize Newmark scheme
    bool steady = get_value(specs_, "/TimeStepping/LinearElasticity/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/LinearElasticity/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    double gamma = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);

    // Set initial conditions
    u_ = Xhv_->element();
    auto u0_ = Xhv_->element();

    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    auto init_displ = expr<Dim,1>(get_value(specs_, "/InitialConditions/LinearElasticity/displacement/expr", default_displ ));
    
    u0_.on(_range=elements(mesh_), _expr=init_displ);
    
    ts_ = newmarkContact(Xhv_, steady, initial_time, final_time, time_step, gamma, beta );
    
    ts_->start();
    ts_->initialize( u0_ );
    u_ = u0_;
    

    LOG(INFO) << "The step is  " << ts_->timeStep() << "\n"
              << "The initial time is " << ts_->timeInitial() << "\n"
              << "The final time is " << ts_->timeFinal() << "\n";

    ts_->updateFromDisp(u_);
}

// Initialization of the contact terms
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::initializeContact()
{
    // Define Xh
    Xh_ = Pch<Order>(mesh_);

    // Initialize contact field
    contactRegion_ =  project(_space=Xh_, _range=elements(mesh_), _expr = cst(0.));
    contactPressure_ =  project(_space=Xh_, _range=elements(mesh_), _expr = cst(0.));
    contactDisplacement_ = project(_space=Xh_, _range=elements(mesh_), _expr = cst(0.)); 
    nbrFaces_ = 0;
    
    // Get contact parameters
    std::string matMethod = fmt::format( "/Collision/LinearElasticity/method" );
    method_ = specs_[nl::json::json_pointer( matMethod )].get<std::string>();

    std::string matEpsilon = fmt::format( "/Collision/LinearElasticity/epsilon" );
    epsilon_ = specs_[nl::json::json_pointer( matEpsilon )].get<double>(); 
    
    std::string matDirection = fmt::format( "/Collision/LinearElasticity/direction");
    direction_ = specs_[nl::json::json_pointer( matDirection )].get<std::string>();

    std::string matDirectionD = fmt::format( "/Collision/LinearElasticity/ddirection");
    ddirection_ = specs_[nl::json::json_pointer( matDirectionD )].get<std::vector<double>>();

    std::string matTheta = fmt::format( "/Collision/LinearElasticity/theta" );
    theta_ = specs_[nl::json::json_pointer( matTheta )].get<double>(); 
    
    std::string matGamma0 = fmt::format( "/Collision/LinearElasticity/gamma0" );
    gamma0_ = specs_[nl::json::json_pointer( matGamma0 )].get<double>(); 
    gamma_ = gamma0_/H_;

    std::string matwithMarker = fmt::format( "/Collision/LinearElasticity/marker" );
    withMarker_ = specs_[nl::json::json_pointer( matwithMarker )].get<int>();  
    
    std::string mattolContactRegion = fmt::format( "/Collision/LinearElasticity/tolContactRegion" );
    tolContactRegion_ = specs_[nl::json::json_pointer( mattolContactRegion )].get<double>();  

    std::string mattolDistance = fmt::format("/Collision/LinearElasticity/tolDistance");
    tolDistance_ = specs_[nl::json::json_pointer( mattolDistance )].get<double>();  
    
    std::string matcomputeerror = fmt::format( "/Collision/LinearElasticity/error" );
    computeerror_ = specs_[nl::json::json_pointer( matcomputeerror )].get<int>(); 
    
    std::string matsave = fmt::format( "/Collision/LinearElasticity/save" );
    save_ = specs_[nl::json::json_pointer( matsave )].get<int>();    

    std::string matFixedPointtol = fmt::format( "/Collision/LinearElasticity/fixedPointTol" );
    fixedPointtol_ = specs_[nl::json::json_pointer( matFixedPointtol )].get<double>(); 
    
    std::string matFixedPoint = fmt::format( "/Collision/LinearElasticity/fixedPoint" );
    fixedPoint_ = specs_[nl::json::json_pointer( matFixedPoint )].get<int>(); 

    std::string matpressurePoint = fmt::format( "/Collision/LinearElasticity/pressurePoint" );
    pressurePoint_ = specs_[nl::json::json_pointer( matpressurePoint )].get<std::vector<double>>();  

}

// Process loading
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processLoading(form1_type& l)
{
    //l += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_),_quad=quad_,_quad1=quad1_ );
    l += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_));
}

// Process materials
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processMaterials( form2_type &a )
{
    auto deft = sym(gradt(u_));
    auto def = sym(grad(u_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    //a += integrate( _range = elements(mesh_), _expr = cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u_),id( u_ ) ) + inner(sigmat,def),_quad=quad_,_quad1=quad1_);
    a += integrate( _range = elements(mesh_), _expr = cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u_),id( u_ ) ) + inner(sigmat,def));
}

// Process contact conditions Penalty method
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processContactPenalty(form1_type& l, form2_type& a, Range<mesh_t, MESH_FACES> const& elts , elementv_t const& u )
{    
    if (withMarker_ == 0)
    {
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=elts, _expr = cst(1.));

        //a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces),_quad=quad_,_quad1=quad1_);
        //l += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces),_quad=quad_,_quad1=quad1_);     
        a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
        l += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));     
    
    }
    else 
    {
        //a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)),_quad=quad_,_quad1=quad1_);
        //l += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)),_quad=quad_,_quad1=quad1_);   
        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)));
        l += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)));   
    
    }
}


// Process persistency conditions
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processContactPersistency(form1_type& l, form2_type& a, Range<mesh_t, MESH_FACES> const& elts, elementv_t const& u)
{
    if (withMarker_ == 0)
    {
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=elts, _expr = cst(1.));

        a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*(ts_->polyFirstDerivCoefficient()*idt(u_)-idv(ts_->polyFirstDeriv())),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
    
        //a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*(ts_->polyFirstDerivCoefficient()*idt(u_)-idv(ts_->polyFirstDeriv())),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces),_quad=quad_,_quad1=quad1_);
    }
    else 
        //a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*(ts_->polyFirstDerivCoefficient()*idt(u_)-idv(ts_->polyFirstDeriv())),trans(expr<Dim,1>(direction_))*id(u)),_quad=quad_,_quad1=quad1_);
        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*(ts_->polyFirstDerivCoefficient()*idt(u_)-idv(ts_->polyFirstDeriv())),trans(expr<Dim,1>(direction_))*id(u)));

}

// Process contact conditions Nitsche method
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processContactNitsche(form1_type& l, form2_type& a, Range<mesh_t, MESH_FACES> const& elts , elementv_t const& u )
{    
    if (withMarker_ == 0)
    {
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=elts, _expr = cst(1.));

        auto const Id = eye<Dim,Dim>();
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
        auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

        a += integrate (_range=boundaryfaces(mesh_),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces)); 
        a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces));
        l += integrate (_range=boundaryfaces(mesh_),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces)); 
    

        //a += integrate (_range=boundaryfaces(mesh_),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces),_quad=quad_,_quad1=quad1_); 
        //a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces),_quad=quad_,_quad1=quad1_);
        //l += integrate (_range=boundaryfaces(mesh_),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces),_quad=quad_,_quad1=quad1_); 
    
    }
    else 
    {
        auto const Id = eye<Dim,Dim>();
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
        auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)); 
        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));
        l += integrate (_range=markedfaces(mesh_,"contact"),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));     
    

        //a += integrate (_range=markedfaces(mesh_,"contact"),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma),_quad=quad_,_quad1=quad1_); 
        //a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma),_quad=quad_,_quad1=quad1_);
        //l += integrate (_range=markedfaces(mesh_,"contact"),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma),_quad=quad_,_quad1=quad1_);     
    }
    
}


// Process boundary conditions
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processBoundaryConditions(form1_type& l, form2_type& a)
{
    // Boundary Condition Dirichlet
    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Dirichlet") )
    {
        for ( auto [key, bc] : specs_["/BoundaryConditions/LinearElasticity/Dirichlet"_json_pointer].items() )
        {
            std::cout << "Add Dirichlet conditions" << std::endl;
            LOG( INFO ) << fmt::format( "Dirichlet conditions found: {}", key );
            std::string e = fmt::format("/BoundaryConditions/LinearElasticity/Dirichlet/{}/g/expr",key);
            auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
            LOG(INFO) << "BoundaryCondition Dirichlet : " << bc_dir << std::endl;
            a+=on(_range=markedfaces(mesh_,key), _rhs=l, _element=u_, _expr=expr<Dim,1>( bc_dir ) );
        }
    }
}


// Time loop
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::timeLoop()
{
    // Initialize linear and bilinear forms
    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_ = form1( _test = Xhv_ );
    auto lt_ = form1( _test = Xhv_ );
    
    a_.zero();
    at_.zero();
    l_.zero();
    lt_.zero();

    std::cout << "***** Process loading *****" << std::endl;
    processLoading(l_);

    std::cout << "***** Process materials *****" << std::endl;
    processMaterials(a_);

    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        lt_ = l_;
        at_ = a_;

        //lt_ +=  integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ),_quad=quad_,_quad1=quad1_ );

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_ +=  integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ));
        
        std::cout << "***** Process contact *****" << std::endl;
        if (withMarker_ == 0)
        {
            myelts_ = getContactRegion(u_);
            std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        }
        else 
            contactRegion_.on( _range=markedfaces(mesh_,"contact"), _expr = cst(1.));
        
        if (method_.compare("penalty") == 0)
            processContactPenalty(lt_, at_, myelts_, u_);
        else if (method_.compare("persistency") == 0)
            processContactPersistency(lt_, at_, myelts_, u_);
        else if (method_.compare("nitsche") == 0)
            processContactNitsche(lt_, at_, myelts_, u_);
        
        std::cout << "***** Process boundary conditions *****" << std::endl;
        processBoundaryConditions(lt_, at_);

        std::cout << "***** Solve *****" << std::endl;
        at_.solve( _rhs = lt_, _solution = u_ );

        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

        // Reset
        at_.zero();
        lt_.zero();

    }    
}

template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::timeLoopFixedPoint()
{
    // Initialize linear and bilinear forms
    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_tmp =  form2( _test = Xhv_, _trial = Xhv_ );
    
    auto l_ = form1( _test = Xhv_ );
    auto lt_ = form1( _test = Xhv_ );
    auto lt_tmp = form1( _test = Xhv_ );
    
    a_.zero();
    at_.zero();
    at_tmp.zero();

    l_.zero();
    lt_.zero();
    lt_tmp.zero();

    std::cout << "***** Process loading *****" << std::endl;
    processLoading(l_);

    std::cout << "***** Process materials *****" << std::endl;
    processMaterials(a_);


    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        lt_ = l_;
        at_ = a_;

        //lt_ +=  integrate( _range=markedelements( mesh_, material.get<std::string>() ), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ),_quad=quad_,_quad1=quad1_ );
        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_ +=  integrate( _range=markedelements( mesh_, material.get<std::string>() ), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ));
        
        auto u_tmp =  Xhv_->element();
        u_tmp.on( _range=elements(mesh_), _expr = idv(u_));

        auto u_tmpNew = Xhv_->element();
        u_tmpNew.on( _range=elements(mesh_), _expr = idv(u_)); 

        int fixedPointIteration = 0;
        double fixedPointerror = 0.;

        while ((fixedPointerror > fixedPointtol_) || (fixedPointIteration < 1))
        {

            lt_tmp = lt_;
            at_tmp = at_;

            u_tmp.on( _range=elements(mesh_), _expr = idv(u_tmpNew)); ;

            std::cout << "***** Process contact *****" << std::endl;
            if (withMarker_ == 0)
            {
                myelts_ = getContactRegion(u_tmp);
                std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
            }
            else 
                contactRegion_.on(_range=markedfaces(mesh_,"contact"), _expr = cst(1.));
            
            if ((nbrFaces_ > 0) || (withMarker_ == 1) )
            {
                std::cout << "Fixed point iteration : " << fixedPointIteration << std::endl;
                std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
                std::cout << "Error : " << fixedPointerror << std::endl;

                if (method_.compare("penalty") == 0)
                    processContactPenalty(lt_tmp, at_tmp, myelts_, u_tmp);
                else if (method_.compare("persistency") == 0)
                    processContactPersistency(lt_tmp, at_tmp, myelts_, u_tmp);
                else if (method_.compare("nitsche") == 0)
                    processContactNitsche(lt_tmp, at_tmp, myelts_, u_tmp);
                
                std::cout << "***** Process boundary conditions *****" << std::endl;
                processBoundaryConditions(lt_tmp, at_tmp);
                
                at_tmp.solve(_rhs = lt_tmp, _solution = u_tmpNew);

                //                fixedPointerror = integrate(_range=elements(mesh_), _expr = norm2( idv(u_tmp)-idv(u_tmpNew)),_quad=quad_,_quad1=quad1_).evaluate()(0,0) / integrate(_range=elements(mesh_),_expr=norm2(idv(u_)),_quad=quad_,_quad1=quad1_).evaluate()(0,0); 
                fixedPointerror = integrate(_range=elements(mesh_), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(mesh_),_expr=norm2(idv(u_))).evaluate()(0,0); 
                fixedPointIteration++;
            }
            else if ((fixedPointIteration == 0) || (nbrFaces_ == 0))
                break;
            
                
            
            if (fixedPointIteration == 10)
                break;
            
            
            lt_tmp.zero();
            at_tmp.zero();
        }
        
        if (withMarker_ == 0)
        {
            myelts_ = getContactRegion(u_tmpNew);
            std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        }
        else 
            contactRegion_.on( _range=markedfaces(mesh_,"contact"), _expr = cst(1.));

        if (method_.compare("penalty") == 0)
            processContactPenalty(lt_, at_, myelts_, u_tmpNew);
        else if (method_.compare("persistency") == 0)
            processContactPersistency(lt_, at_, myelts_, u_tmpNew);
        else if (method_.compare("nitsche") == 0)
            processContactNitsche(lt_, at_, myelts_, u_tmpNew);
                
        processBoundaryConditions(lt_, at_);

        at_.solve( _rhs = lt_, _solution = u_ );
        
        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

        // Reset
        at_.zero();
        lt_.zero();
    }
}

// Run method 
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::run()
{
    std::cout << "***** Run dynamic elasticity with unilateral contact *****" << std::endl;

    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    std::cout <<  "***** Initialize distance g *****" << std::endl;
    initG();

    std::cout <<  "***** Start time loop *****" << std::endl;
    
    if ((method_.compare("penalty") == 0) || (method_.compare("nitsche") == 0) || (method_.compare("persistency") == 0) )
    {
        this->exportResults(0);

        if (fixedPoint_ == 1)
            timeLoopFixedPoint();
        else 
            timeLoop();
                
    }

}


template<int Dim, int Order, int OrderGeo>
Range<typename ContactDynamic<Dim, Order, OrderGeo>::mesh_t, MESH_FACES>
ContactDynamic<Dim, Order, OrderGeo>::getContactRegion( elementv_t const& u )
{   
    Range<mesh_t,MESH_FACES> myelts( mesh_ );

    auto defv = sym(gradv(u));
    auto Id = eye<Dim,Dim>();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
    
    contactRegion_.on( _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));    

    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=Xh_, _trial=Xh_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(mesh_) )
    {                
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : Xh_->dof()->faceLocalDof( face.id() ) )
        {
            index_type thedof = ldof.index();
            thedof = trialDofIdToContainerId[ thedof ];

            if (contactRegion_[thedof] >= tolContactRegion_)
                contactDof++;
                    
            if (Order == 1)
            {
                if (Dim == 2)
                {
                    if (contactDof == 2)
                    {
                        nbrFaces_++;
                        myelts.push_back( face );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 3)
                    {
                        nbrFaces_++;
                        myelts.push_back( face );
                    }
                }         
            }
            else if (Order == 2)
            {
                if (Dim == 2)
                {
                    if (contactDof == 3)
                    {
                        nbrFaces_++;
                        myelts.push_back( face );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 4)
                    {
                        nbrFaces_++;
                        myelts.push_back( face );
                    }
                } 
            }
        }
    }
    myelts.shrink_to_fit();   
  
    return myelts;
}

template <int Dim, int Order, int OrderGeo>
void 
ContactDynamic<Dim, Order, OrderGeo>::initG()
{   
    // Init the distance fields 
    g_ = Xh_->element();

    // Load and export rigid obstacles
    auto wall = loadMesh(_mesh=new mesh_t, _filename = "$cfgdir/wall.geo",_h = 0.4);

    auto expWall = Feel::exporter(_mesh = wall, _name = fmt::format("Wall"));
    expWall->addRegions();
    expWall->save();
    
    // Raytracing to compute distance
    using bvh_ray_type = BVHRay<Dim>;
    Eigen::VectorXd origin(Dim);
    Eigen::VectorXd dir(Dim);
    
    if constexpr(Dim == 2)
        dir << ddirection_[0], ddirection_[1];
    else if constexpr(Dim == 3)
        dir << ddirection_[0], ddirection_[1], ddirection_[2];
    
    for (auto const& theface : boundaryfaces(mesh_) )
    {                
        auto & face = boost::unwrap_ref( theface );

        auto &point = face.point(0);
        if (point.isOnBoundary())
        {    

            if constexpr(Dim == 2)
                origin << point.node()[0], point.node()[1];
            else if constexpr(Dim == 3)
                origin << point.node()[0], point.node()[1], point.node()[2];

            bvh_ray_type ray(origin,dir);
            auto bvh = boundingVolumeHierarchy(_range=boundaryfaces(wall));
            auto rayIntersection = bvh->intersect(_ray=ray) ;

            if (!rayIntersection.empty()) 
            {
                for ( auto const& rir : rayIntersection )
                {
                    for (auto const& ldof  : Xh_->dof()->faceLocalDof( face.id() ))
                        g_[ldof.index()] = rir.distance() - tolDistance_;
                }
            } 
        }    
    }
}


template <int Dim, int Order, int OrderGeo>
void 
ContactDynamic<Dim, Order, OrderGeo>::storeData()
{
    u_.saveHDF5("solution_ref.h5");
    saveGMSHMesh(_mesh=mesh_,_filename= "mesh_ref.msh" );
}


template <int Dim, int Order, int OrderGeo>
void 
ContactDynamic<Dim, Order, OrderGeo>::error()
{
    auto meshref = loadMesh(_mesh=new mesh_t, _filename = "/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/mesh_ref.msh");
    auto uref = spacev_t::New(meshref)->elementPtr() ;
    uref->loadHDF5("/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/solution_ref.h5");

    auto op_inter = opInterpolation(_domainSpace =  spacev_t::New(mesh_), _imageSpace = spacev_t::New(meshref) );
    auto uinter =  spacev_t::New(meshref)->element(); 
    op_inter->apply(u_, uinter);

    auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
    std::cout << "H1 error : " << h1err << std::endl;
}

// Export results
template <int Dim, int Order, int OrderGeo>
void 
ContactDynamic<Dim, Order, OrderGeo>::exportResults(double t)
{
    /*
    if (OrderGeo == 1)
    {
        // Define exports
        e_->step(t)->addRegions();
        e_->step(t)->add( "displacement", u_ );
        e_->step(t)->add( "velocity", ts_->currentVelocity() );

        // Compute new contact region 
        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for export : " << nbrFaces_ << std::endl;
    
        auto realcontactRegion = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_));  
        e_->step(t)->add( "realcontactRegion", realcontactRegion) ;

        auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts_->begin(), myelts_->end(), myelts_ );
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = project(_space=XhCFaces, _range=myfaces, _expr = cst(1.));

        auto const Id = eye<Dim,Dim>();
        auto defv = sym(gradv(u_));
        auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
        contactRegion_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = idv(contactFaces));
        contactPressure_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
        contactDisplacement_ = project(_space=Xh_, _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_))*idv(contactFaces));

        auto expression = Xh_->element();
        if (method_.compare("penalty") == 0)
            expression = project(_space=Xh_, _range=elements(mesh_), _expr =  (trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_))*idv(contactFaces)  );
        else if (method_.compare("nitsche") == 0)
            expression = project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = (cst(gamma_)*(trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_)) - trans(expr<Dim,1>(direction_))*sigmav)*idv(contactFaces)  );
        
        e_->step(t)->add( "contactRegion", contactRegion_) ;
        e_->step(t)->add( "contactPressure", contactPressure_);
        e_->step(t)->add( "contactDisplacement", contactDisplacement_ );
        e_->step(t)->add( "expression", expression);
        e_->step(t)->add("g", g_);
        e_->save();

        // Save values in json file
        meas_["time"].push_back(ts_->time());

        auto sig = lambda_*trace(defv)*Id + 2*mu_*defv;
        auto J = det(Id + gradv(u_));


        double E1 = 0.5*rho_*normL2Squared(_range=elements(mesh_),_expr=idv(ts_->currentVelocity()));
        double E2 = 0.5*integrate( _range= elements( mesh_ ), _expr= inner(sig,defv) ).evaluate()( 0,0 );
        
        meas_["Eh1"].push_back(E1);
        meas_["Eh2"].push_back(E2);
        meas_["Eh"].push_back(E1 + E2);

        double disp = integrate( _range= boundaryfaces(mesh_), _expr= (trans(expr<Dim,1>(direction_))*idv(u_)  - idv(g_)) * idv(contactFaces )).evaluate()( 0,0 );
        meas_["disp"].push_back(disp);

        double Lv = integrate( _range=elements(mesh_),_expr=  cst(rho_)*abs(trans(expr<Dim,1>( externalforce_ )))*idv(u_)).evaluate()( 0,0 );
        meas_["Lv"].push_back(Lv);

        double volume = integrate(_range=elements(mesh_), _expr = det(Id + gradv(u_))).evaluate()( 0, 0 );
        meas_["volume"].push_back(volume);


        auto ctx = Xh_->context();
        node_type t1(Dim);
       
        if (Dim == 2)
        {
            t1(0)=pressurePoint_[0]; t1(1)=pressurePoint_[1];
        }
        else 
        {
            t1(0)=pressurePoint_[0]; t1(1)=pressurePoint_[1]; t1(2)=pressurePoint_[2];
        }    
                
        ctx.add( t1 );

        auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure_) ); 
        auto evaluateDispExpr = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
            
        meas_["evaluateStress"].push_back(evaluateStress(0,0));
        meas_["evaluateDisp"].push_back(evaluateDisp(0,0));
    
        if (method_.compare("penalty") == 0)
            meas_["E"].push_back(E1 + E2 + Lv);
        else if (method_.compare("nitsche") == 0)
        {
            double R1 = normL2Squared(_range= boundaryfaces(mesh_), _expr= sqrt(cst(gamma0_)/cst(gamma_)) * trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
            double R2 = normL2Squared( _range= boundaryfaces(mesh_), _expr= sqrt(cst(gamma0_)/cst(gamma_)) * ( cst(gamma_) * ( trans(expr<Dim,1>(direction_)) *idv(u_)  - idv(g_) ) - trans(expr<Dim,1>(direction_))*sigmav ) * idv(contactFaces));

            meas_["R1"].push_back(R1);
            meas_["R2"].push_back(R2);
            meas_["R"].push_back((R1 - R2)/(2.*gamma0_));

            double Es = (E1 + E2) - theta_*(R1 - R2)/(2.*gamma0_);
            meas_["Es"].push_back(Es);      

            meas_["E"].push_back(Es + Lv);
        }
    
        this->writeResultsToFile("measures.json");
    }
    */
    // Interpolation
    e_->step(t)->addRegions();
        
    auto Xhv_P1 = Pchv<Order>(meshP1); 
    auto uinter =  Xhv_P1->element(); 

    auto op_inter = opInterpolation(_domainSpace =  Xhv_, _imageSpace = Xhv_P1 );
    op_inter->apply(u_, uinter);


    e_->step(t)->add( "displacement", uinter );

    auto const Id = eye<Dim,Dim>();
    auto defvinter = sym(gradv(uinter));
    auto sigmavinter = (lambda_*trace(defvinter)*Id + 2*mu_*defvinter)*N();
   
    auto contactPressureinter =  Pch<Order>(meshP1)->element();
    contactPressureinter.on( _range=boundaryfaces(meshP1), _expr = trans(expr<Dim,1>(direction_))*sigmavinter);

    auto contactDisplacementinter = Pch<Order>(meshP1)->element();
    contactDisplacementinter.on( _range=elements(meshP1), _expr = (trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)));

    e_->step(t)->add( "contactPressure", contactPressureinter);
    e_->step(t)->add( "contactDisplacement", contactDisplacementinter );
    e_->save();

    myelts_ = getContactRegion(u_);
    std::cout << "Nbr faces for export : " << nbrFaces_ << std::endl;
    
    auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
    auto XhCFaces = Pdh<0>(face_mesh);
    auto contactFaces =XhCFaces->element();
    contactFaces.on(_range=myelts_, _expr = cst(1.));

    auto defv = sym(gradv(u_));
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();

    contactPressure_.on(_range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));


    // Save values in json file
    meas_["time"].push_back(ts_->time());

    auto sig = lambda_*trace(defv)*Id + 2*mu_*defv;
    auto J = det(Id + gradv(u_));

    double E1 = 0.5*rho_*normL2Squared(_range=elements(mesh_),_expr=idv(ts_->currentVelocity()));
    //    double E2 = 0.5*integrate( _range= elements( mesh_ ), _expr= inner(sig,defv),_quad=quad_,_quad1=quad1_ ).evaluate()( 0,0 );

    double E2 = 0.5*integrate( _range= elements( mesh_ ), _expr= inner(sig,defv)).evaluate()( 0,0 );
        
    meas_["Eh1"].push_back(E1);
    meas_["Eh2"].push_back(E2);
    meas_["Eh"].push_back(E1 + E2);

    //    double disp = integrate( _range= boundaryfaces(mesh_), _expr= (trans(expr<Dim,1>(direction_))*idv(u_)  - idv(g_)) * idv(contactFaces ),_quad=quad_,_quad1=quad1_).evaluate()( 0,0 );

    double disp = integrate( _range= boundaryfaces(mesh_), _expr= (trans(expr<Dim,1>(direction_))*idv(u_)  - idv(g_)) * idv(contactFaces )).evaluate()( 0,0 );
    meas_["disp"].push_back(disp);

    //    double Lv = integrate( _range=elements(mesh_),_expr=  cst(rho_)*abs(trans(expr<Dim,1>( externalforce_ )))*idv(u_),_quad=quad_,_quad1=quad1_).evaluate()( 0,0 );

    double Lv = integrate( _range=elements(mesh_),_expr=  cst(rho_)*abs(trans(expr<Dim,1>( externalforce_ )))*idv(u_)).evaluate()( 0,0 );
    meas_["Lv"].push_back(Lv);

    //    double volume = integrate(_range=elements(mesh_), _expr = det(Id + gradv(u_)),_quad=quad_,_quad1=quad1_).evaluate()( 0, 0 );

    double volume = integrate(_range=elements(mesh_), _expr = det(Id + gradv(u_))).evaluate()( 0, 0 );
    meas_["volume"].push_back(volume);

    auto ctx = Xh_->context();
    node_type t1(Dim);
       
    if (Dim == 2)
    {
        t1(0)=pressurePoint_[0]; t1(1)=pressurePoint_[1];
    }
    else 
    {
        t1(0)=pressurePoint_[0]; t1(1)=pressurePoint_[1]; t1(2)=pressurePoint_[2];
    }    
                
    ctx.add( t1 );

    auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure_) ); 
    auto evaluateDispExpr = Xh_->element();
    evaluateDispExpr.on(_range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
    auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
            
    meas_["evaluateStress"].push_back(evaluateStress(0,0));
    meas_["evaluateDisp"].push_back(evaluateDisp(0,0));
    
    if ((method_.compare("penalty") == 0) || (method_.compare("persistency") == 0))
        meas_["E"].push_back(E1 + E2 + Lv);
    else if (method_.compare("nitsche") == 0)
    {
        double R1 = normL2Squared(_range= boundaryfaces(mesh_), _expr= sqrt(cst(gamma0_)/cst(gamma_)) * trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
        double R2 = normL2Squared( _range= boundaryfaces(mesh_), _expr= sqrt(cst(gamma0_)/cst(gamma_)) * ( cst(gamma_) * ( trans(expr<Dim,1>(direction_)) *idv(u_)  - idv(g_) ) - trans(expr<Dim,1>(direction_))*sigmav ) * idv(contactFaces));

        meas_["R1"].push_back(R1);
        meas_["R2"].push_back(R2);
        meas_["R"].push_back((R1 - R2)/(2.*gamma0_));

        double Es = (E1 + E2) - theta_*(R1 - R2)/(2.*gamma0_);
        meas_["Es"].push_back(Es);      

        meas_["E"].push_back(Es + Lv);
    }
    
    this->writeResultsToFile("measures.json");
}


template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::writeResultsToFile(const std::string& filename) const
{
    if ( Environment::isMasterRank() )
    {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << meas_.dump(4);  // Indent of 4 spaces for readability
            file.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }
}

} 