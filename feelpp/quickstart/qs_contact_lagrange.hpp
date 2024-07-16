#pragma once
#include <iostream>

#include <chrono>
#include <fmt/chrono.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>


namespace Feel
{

template <int Dim, int Order, int OrderGeo>
class ContactDynamicLagrange
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
    using ts_ptrtype = std::shared_ptr<Newmark<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_tP1>>; 

    // Constructors
    ContactDynamicLagrange() = default;
    ContactDynamicLagrange(nl::json const& specs);

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
    void processLoading(form1_type& l, elementv_t const& u);
    void processMaterials(form2_type &a, elementv_t const& u);
    void run();
    typename MeshTraits<Mesh<Simplex<Dim, OrderGeo>>>::faces_reference_wrapper_ptrtype getContactRegion(elementv_t const& u);
    void timeLoop();
    void timeLoopSimple();
    typename MeshTraits<Mesh<Simplex<Dim, OrderGeo>>>::faces_reference_wrapper_ptrtype getRandomFace();
    void writeResultsToFile(const std::string& filename) const;
    void initG();

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
    typename MeshTraits<Mesh<Simplex<Dim, OrderGeo>>>::faces_reference_wrapper_ptrtype myelts_;

    ts_ptrtype ts_;
    exporter_ptrtype e_;
    nl::json meas_;



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
ContactDynamicLagrange<Dim, Order, OrderGeo>::ContactDynamicLagrange(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order, int OrderGeo>
void ContactDynamicLagrange<Dim, Order, OrderGeo>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // Define Xhv
    Xhv_ = Pchv<Order>(mesh_); 

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
    
    ts_ = newmark( _space = Xhv_, _steady=steady, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _order=time_order );
    
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
void ContactDynamicLagrange<Dim, Order, OrderGeo>::initializeContact()
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
void ContactDynamicLagrange<Dim, Order, OrderGeo>::processLoading(form1_type& l, elementv_t const& u)
{
    l += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u) );
}

// Process materials
template <int Dim, int Order, int OrderGeo>
void ContactDynamicLagrange<Dim, Order, OrderGeo>::processMaterials( form2_type &a, elementv_t const& u )
{
    auto deft = sym(gradt(u));
    auto def = sym(grad(u));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a += integrate( _range = elements(mesh_), _expr = cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) + inner(sigmat,def));
}

// Time loop
template <int Dim, int Order, int OrderGeo>
void ContactDynamicLagrange<Dim, Order, OrderGeo>::timeLoop()
{   
    if (Order != 2)
        std::cout << "u has to be P2" << std::endl;
    
    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    // Init solution u
    auto uRes = Xhv_->element();
    uRes = u_; // u intial

    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_ = form1( _test = Xhv_ );

    a_.zero();
    l_.zero();

    processLoading(l_,uRes);
    processMaterials(a_,uRes);

    // export init time
    e_->step(0)->addRegions();

    auto Xhv_P1 = Pchv<Order>(meshP1); 
    auto uinter =  Xhv_P1->element(); 

    auto op_inter = opInterpolation(_domainSpace =  Xhv_, _imageSpace = Xhv_P1 );
    op_inter->apply(uRes, uinter);

    auto expression = Pch<Order>(meshP1)->element();

    e_->step(0)->add( "displacement", uinter );
    //e_->step(0)->add( "expression", expression );
    //e_->step(0)->add( "realcontactRegion", uinter) ;
    //e_->step(0)->add( "contactPressure", contactPressure_);
    //e_->step(0)->add( "contactDisplacement", contactDisplacement_ );


    e_->save();

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(uRes) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        // Compute faces in contact
        myelts_ = getContactRegion(uRes);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;

        auto Id = eye<Dim,Dim>();
        
        if (nbrFaces_ > 0)
        {
            std::cout << "Contact" << std::endl;
            // Contact
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts_->begin(), myelts_->end(), myelts_ );
        
            // Define lagrange multiplier
            auto submesh = createSubmesh(_mesh=mesh_,_range=myfaces);
            auto XhLambda = Pch<1>(submesh);
            auto lambda = XhLambda->element();

            auto ps = product( Xhv_, XhLambda );

            // Process loading and materials

            std::cout << "create a and rhs" << std::endl;
            auto a = blockform2( ps ,solve::strategy::monolithic,backend() );
            auto rhs = blockform1( ps,solve::strategy::monolithic,backend() );

            //a.zero();
            //rhs.zero();

            std::cout << "add rhs" << std::endl;
            //rhs(0_c) += l_;
            rhs(0_c) = integrate(_range=elements(mesh_), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(uRes));
            rhs(0_c) += integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( uRes ) ) );
            rhs(1_c) = integrate(_range=elements(submesh), _expr=idv(g_)*id(lambda));

            std::cout << "define a" << std::endl;
            auto deft = sym(gradt(uRes));
            auto def = sym(grad(uRes));
            auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

            std::cout << "add a" << std::endl;
            //a(0_c,0_c) = a_;
            a(0_c,0_c) = integrate(_range=elements(mesh_),_expr=inner(sigmat,def) );
            std::cout << "add a 1" << std::endl;
            a(0_c,1_c) += integrate( _range=elements(submesh),_expr=idt(lambda)*(trans(expr<Dim,1>(direction_))*id(uRes)));
            std::cout << "add a 2" << std::endl;
            a(1_c,0_c) += integrate( _range=elements(submesh),_expr=(trans(expr<Dim,1>(direction_))*idt(uRes))*id(lambda) );
            //a(1_c,1_c) += integrate( _range=elements(submesh),_expr=cst(1.)/cst(gamma_)*idt(lambda)*id(lambda) );

            std::cout << "***** Solve *****" << std::endl;
            auto U=ps.element();
            a.solve( _solution=U, _rhs=rhs, _rebuild=true);

            uRes = U(0_c);

            ts_->updateFromDisp(uRes);
        }
        else 
        {
            std::cout << "No Contact" << std::endl;
            // No contact
            
            /*
            auto elts = getRandomFace();
            auto myface =  boost::make_tuple( mpl::size_t<MESH_FACES>(), elts->begin(), elts->end(), elts );
            auto submesh = createSubmesh(_mesh=mesh_,_range=myface);
            auto XhLambda = Pch<1>(submesh);
            auto lambda = XhLambda->element();

            auto ps = product( Xhv_, XhLambda );

            // Process loading and materials
            auto a = blockform2( ps,solve::strategy::monolithic,backend() );
            auto rhs = blockform1( ps,solve::strategy::monolithic,backend() );
            */


            auto a = form2( _test = Xhv_, _trial = Xhv_ );
            auto rhs = form1( _test = Xhv_ ); 

            a.zero();
            rhs.zero();

            auto deft = sym(gradt(uRes));
            auto def = sym(grad(uRes));
            auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

            //a = integrate(_range=elements(mesh_),_expr=inner(sigmat,def) );
            a = a_;
            rhs = l_;
            //rhs = integrate(_range=elements(mesh_), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(uRes));
            rhs += integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( uRes ) ) );
            
            std::cout << "***** Solve *****" << std::endl;
            //auto U=ps.element();
            a.solve( _rhs=rhs, _solution=uRes, _rebuild=true);

            //uRes = U(0_c);
            ts_->updateFromDisp(uRes);
        }

        e_->step(ts_->time())->addRegions();

        op_inter->apply(uRes, uinter);

        e_->step(ts_->time())->add( "displacement", uinter );

        //expression = project(_space=Pch<Order>(meshP1), _range=boundaryfaces(mesh_), _expr = cst(gamma_)*(trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)) - ); 

        //auto realcontactRegion = project(_space=Pch<Order>(meshP1), _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_));  
        //e_->step(ts_->time())->add( "realcontactRegion", realcontactRegion) ;
   
        //auto defv = sym(gradv(uinter));
        //auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
        //auto contactPressure =  project(_space=Pch<Order>(meshP1), _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav);
        //auto contactDisplacement = project(_space=Pch<Order>(meshP1), _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)));

        //e_->step(ts_->time())->add( "contactPressure", contactPressure);
        //e_->step(ts_->time())->add( "contactDisplacement", contactDisplacement );
        e_->save();

        meas_["time"].push_back(ts_->time());

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
        auto evaluateDispExpr = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(uRes));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
            
        meas_["evaluateStress"].push_back(evaluateStress(0,0));
        meas_["evaluateDisp"].push_back(evaluateDisp(0,0));
    

        this->writeResultsToFile("measures.json");
    }
}

template <int Dim, int Order, int OrderGeo>
void ContactDynamicLagrange<Dim, Order, OrderGeo>::timeLoopSimple()
{
    if (Order != 2)
        std::cout << "u has to be P2" << std::endl;
    
    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    // Init solution u
    auto uRes = Xhv_->element();
    uRes = u_; // u intial

    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_ = form1( _test = Xhv_ );

    a_.zero();
    l_.zero();

    processLoading(l_,uRes);
    processMaterials(a_,uRes);

    // Define Lagrange multiplier
    auto submesh = createSubmesh(_mesh=mesh_,_range=markedfaces(mesh_, "contact"));
    auto XhLambda = Pch<1>(submesh);
    auto lambda = XhLambda->element();

    auto ps = product( Xhv_, XhLambda );

    // export init time
    e_->step(0)->addRegions();

    auto Xhv_P1 = Pchv<Order>(meshP1); 
    auto uinter =  Xhv_P1->element(); 

    auto op_inter = opInterpolation(_domainSpace =  Xhv_, _imageSpace = Xhv_P1 );
    op_inter->apply(uRes, uinter);

    auto expression = Pch<Order>(meshP1)->element();

    e_->step(0)->add( "displacement", uinter );
    //e_->step(0)->add( "expression", expression );
    //e_->step(0)->add( "realcontactRegion", uinter) ;
    //e_->step(0)->add( "contactPressure", contactPressure_);
    //e_->step(0)->add( "contactDisplacement", contactDisplacement_ );


    e_->save();

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(uRes) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        // Compute faces in contact
        myelts_ = getContactRegion(uRes);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;

        auto Id = eye<Dim,Dim>();
        
        if (nbrFaces_ > 0)
        {
            std::cout << "Contact" << std::endl;
            // Contact
            
            std::cout << "create a and rhs" << std::endl;
            auto a = blockform2( ps ,solve::strategy::monolithic,backend() );
            auto rhs = blockform1( ps,solve::strategy::monolithic,backend() );

            lambda = XhLambda->element();

            std::cout << "add rhs" << std::endl;

            rhs(0_c) = integrate(_range=elements(mesh_), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(uRes));
            rhs(0_c) += integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( uRes ) ) );
            rhs(1_c) = integrate(_range=elements(submesh), _expr=idv(g_)*id(lambda));

            std::cout << "define a" << std::endl;
            auto deft = sym(gradt(uRes));
            auto def = sym(grad(uRes));
            auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

            std::cout << "add a" << std::endl;
 
            a(0_c,0_c) = integrate(_range=elements(mesh_),_expr=inner(sigmat,def) );
            std::cout << "add a 1" << std::endl;
            a(0_c,1_c) += integrate( _range=elements(submesh),_expr=idt(lambda)*(trans(expr<Dim,1>(direction_))*id(uRes)));
            std::cout << "add a 2" << std::endl;
            a(1_c,0_c) += integrate( _range=elements(submesh),_expr=(trans(expr<Dim,1>(direction_))*idt(uRes))*id(lambda) );
            //a(1_c,1_c) += integrate( _range=elements(submesh),_expr=cst(1.)/cst(gamma_)*idt(lambda)*id(lambda) );

            std::cout << "***** Solve *****" << std::endl;
            auto U=ps.element();
            a.solve( _solution=U, _rhs=rhs, _rebuild=true);

            uRes = U(0_c);
            
            ts_->updateFromDisp(uRes);
        }
        else 
        {
            std::cout << "No Contact" << std::endl;
            
            std::cout << "create a and rhs" << std::endl;
            auto a = blockform2( ps ,solve::strategy::monolithic,backend() );
            auto rhs = blockform1( ps,solve::strategy::monolithic,backend() );

            lambda = XhLambda->element();

            std::cout << "add rhs" << std::endl;

            rhs(0_c) = integrate(_range=elements(mesh_), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(uRes));
            rhs(0_c) += integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( uRes ) ) );
            
            std::cout << "define a" << std::endl;
            auto deft = sym(gradt(uRes));
            auto def = sym(grad(uRes));
            auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

            std::cout << "add a" << std::endl;
 
            a(0_c,0_c) = integrate(_range=elements(mesh_),_expr=inner(sigmat,def) );
            std::cout << "add a 1" << std::endl;
            //a(0_c,1_c) += integrate( _range=elements(submesh),_expr=idt(lambda)*(trans(expr<Dim,1>(direction_))*id(uRes)));
            std::cout << "add a 2" << std::endl;
            //a(1_c,0_c) += integrate( _range=elements(submesh),_expr=(trans(expr<Dim,1>(direction_))*idt(uRes))*id(lambda) );
            a(1_c,1_c) += integrate( _range=elements(submesh),_expr=cst(1.)/cst(gamma_)*idt(lambda)*id(lambda) );

            std::cout << "***** Solve *****" << std::endl;
            auto U=ps.element();
            a.solve( _solution=U, _rhs=rhs, _rebuild=true);

            uRes = U(0_c);

            ts_->updateFromDisp(uRes);
        }

        e_->step(ts_->time())->addRegions();

        op_inter->apply(uRes, uinter);

        e_->step(ts_->time())->add( "displacement", uinter );

        //expression = project(_space=Pch<Order>(meshP1), _range=boundaryfaces(mesh_), _expr = cst(gamma_)*(trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)) - ); 

        //auto realcontactRegion = project(_space=Pch<Order>(meshP1), _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_));  
        //e_->step(ts_->time())->add( "realcontactRegion", realcontactRegion) ;
   
        //auto defv = sym(gradv(uinter));
        //auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
        //auto contactPressure =  project(_space=Pch<Order>(meshP1), _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav);
        //auto contactDisplacement = project(_space=Pch<Order>(meshP1), _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)));

        //e_->step(ts_->time())->add( "contactPressure", contactPressure);
        //e_->step(ts_->time())->add( "contactDisplacement", contactDisplacement );
        e_->save();

        meas_["time"].push_back(ts_->time());

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
        auto evaluateDispExpr = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(uRes));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
            
        meas_["evaluateStress"].push_back(evaluateStress(0,0));
        meas_["evaluateDisp"].push_back(evaluateDisp(0,0));
    

        this->writeResultsToFile("measures.json");
    }
}

// Run method 
template <int Dim, int Order, int OrderGeo>
void ContactDynamicLagrange<Dim, Order, OrderGeo>::run()
{
    std::cout << "***** Run dynamic elasticity with unilateral contact *****" << std::endl;

    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    std::cout <<  "***** Initialize distance g *****" << std::endl;
    initG();

    std::cout <<  "***** Start time loop *****" << std::endl;
    
    if (fixedPoint_ == 1)
        timeLoopSimple();
    else 
        timeLoop();
}



template<int Dim, int Order, int OrderGeo>
typename MeshTraits<Mesh<Simplex<Dim,OrderGeo>>>::faces_reference_wrapper_ptrtype
ContactDynamicLagrange<Dim, Order, OrderGeo>::getRandomFace()
{
    typename MeshTraits<Mesh<Simplex<Dim, OrderGeo>>>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<Mesh<Simplex<Dim,OrderGeo>>>::faces_reference_wrapper_type );

    auto const& trialDofIdToContainerId =  form2(_test=Xh_, _trial=Xh_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(mesh_) )
    {                
        auto & face = boost::unwrap_ref( theface );
        myelts->push_back( boost::cref( face ) );
        break;
    }

    myelts->shrink_to_fit();   
  
    return myelts;
}

template<int Dim, int Order, int OrderGeo>
typename MeshTraits<Mesh<Simplex<Dim,OrderGeo>>>::faces_reference_wrapper_ptrtype
ContactDynamicLagrange<Dim, Order, OrderGeo>::getContactRegion(elementv_t const& u)
{   
    typename MeshTraits<Mesh<Simplex<Dim, OrderGeo>>>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<Mesh<Simplex<Dim,OrderGeo>>>::faces_reference_wrapper_type );

    auto defv = sym(gradv(u));
    auto Id = eye<Dim,Dim>();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
    
    contactRegion_ = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));    

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
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 3)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }         
            }
            else if (Order == 2)
            {
                if (Dim == 2)
                {
                    if (contactDof == 2)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 4)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                } 
            }
        }
    }
    myelts->shrink_to_fit();   
  
    return myelts;
}

template <int Dim, int Order, int OrderGeo>
void 
ContactDynamicLagrange<Dim, Order, OrderGeo>::initG()
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
void ContactDynamicLagrange<Dim, Order, OrderGeo>::writeResultsToFile(const std::string& filename) const
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