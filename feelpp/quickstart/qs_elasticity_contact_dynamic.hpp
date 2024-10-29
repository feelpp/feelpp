#pragma once
#include "qs_elasticity_contact.hpp"

namespace Feel
{

template <int Dim, int Order, int OrderGeo>
class ContactDynamic
{
public:
    using mesh_t = Mesh<Simplex<Dim,OrderGeo>>;
    using spacev_t = Pchv_type<mesh_t, Order>;
    using space_t = Pch_type<mesh_t, Order>;
    using spacev_ptr_t = Pchv_ptrtype<mesh_t, Order>; 
    using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
    using elementv_t = typename spacev_t::element_type;
    using element_t = typename space_t::element_type;
    using form2_type = form2_t<spacev_t,spacev_t>; 
    using form1_type = form1_t<spacev_t>; 
    using ts_ptrtype = std::shared_ptr<Newmark<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

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
    void run();
    Range<mesh_t, MESH_FACES> getContactRegion( elementv_t const& u );
    void timeLoop();
    void timeLoopFixedPoint();
    void exportResults(double t);
    void writeResultsToFile(const std::string& filename) const;
    void initG();

private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
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

    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    double fixedPointtol_;
    int fixedPoint_; 
    std::vector<double> pressurePoint_;
    int nbrObs_;
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
    tic();
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    toc("mesh");
    // Define Xhv
    tic();
    Xhv_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );; 

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
    e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    

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
    
    u0_.on(_range=elements(support(Xhv_)), _expr=init_displ);
    
    ts_ = newmark(_space = Xhv_, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _gamma=gamma, _beta=beta );
    
    ts_->start();
    ts_->initialize( u0_ );
    u_ = u0_;
    
    ts_->updateFromDisp(u_);
    toc("init");

    LOG(INFO) << "The step is  " << ts_->timeStep() << "\n"
              << "The initial time is " << ts_->timeInitial() << "\n"
              << "The final time is " << ts_->timeFinal() << "\n";

    
}

// Initialization of the contact terms
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::initializeContact()
{
    // Define Xh
    tic();
    Xh_ = Pch<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );

    // Initialize contact field
    contactRegion_ =  project(_space=Xh_, _range=elements(support(Xh_)), _expr = cst(0.));
    contactPressure_ =  project(_space=Xh_, _range=elements(support(Xh_)), _expr = cst(0.));
    contactDisplacement_ = project(_space=Xh_, _range=elements(support(Xh_)), _expr = cst(0.)); 
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

    
    std::string mattolContactRegion = fmt::format( "/Collision/LinearElasticity/tolContactRegion" );
    tolContactRegion_ = specs_[nl::json::json_pointer( mattolContactRegion )].get<double>();  

    std::string mattolDistance = fmt::format("/Collision/LinearElasticity/tolDistance");
    tolDistance_ = specs_[nl::json::json_pointer( mattolDistance )].get<double>();  
    
    std::string matFixedPointtol = fmt::format( "/Collision/LinearElasticity/fixedPointTol" );
    fixedPointtol_ = specs_[nl::json::json_pointer( matFixedPointtol )].get<double>(); 
    
    std::string matFixedPoint = fmt::format( "/Collision/LinearElasticity/fixedPoint" );
    fixedPoint_ = specs_[nl::json::json_pointer( matFixedPoint )].get<int>(); 

    std::string matpressurePoint = fmt::format( "/Collision/LinearElasticity/pressurePoint" );
    pressurePoint_ = specs_[nl::json::json_pointer( matpressurePoint )].get<std::vector<double>>();  

    std::string matnbrObs = fmt::format( "/Collision/LinearElasticity/nbrObs");
    nbrObs_ = specs_[nl::json::json_pointer( matnbrObs )].get<int>();
    toc("init");
}

// Process loading
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processLoading(form1_type& l)
{
    l += integrate( _range = elements(support(Xhv_)), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_));
}

// Process materials
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processMaterials( form2_type &a )
{
    auto deft = sym(gradt(u_));
    auto def = sym(grad(u_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;
    
    a += integrate( _range = elements(support(Xhv_)), _expr = cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u_),id( u_ ) ) + inner(sigmat,def));
}

// Process contact conditions Penalty method
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processContactPenalty(form1_type& l, form2_type& a, Range<mesh_t, MESH_FACES> const& elts , elementv_t const& u )
{    
    a += integrate (_range=elts,_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)));
    l += integrate (_range=elts,_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)));     
}

// Process contact conditions Nitsche method
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::processContactNitsche(form1_type& l, form2_type& a, Range<mesh_t, MESH_FACES> const& elts , elementv_t const& u )
{    
    auto const Id = eye<Dim,Dim>();
    auto deft = sym(gradt(u));
    auto def = sym(grad(u));
    auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
    auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

    a += integrate (_range=elts,_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)); 
    a += integrate (_range=elts,_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));
    l += integrate (_range=elts,_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));     
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
            LOG( INFO ) << fmt::format( "Dirichlet conditions found: {}", key );
            std::string e = fmt::format("/BoundaryConditions/LinearElasticity/Dirichlet/{}/g/expr",key);
            auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
            LOG(INFO) << "BoundaryCondition Dirichlet : " << bc_dir << std::endl;
            a+=on(_range=markedfaces(support(Xhv_),key), _rhs=l, _element=u_, _expr=expr<Dim,1>( bc_dir ) );
            a+=on(_range=markedpoints(mesh_,key), _rhs=l, _element=u_, _expr=expr<Dim,1>( bc_dir ) );

        }
    }

    // Boundary Condition Neumann
    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Neumann") )
    {
        for ( auto [key, bc] : specs_["/BoundaryConditions/LinearElasticity/Neumann"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Neumann conditions found: {}", key );
            std::string e = fmt::format("/BoundaryConditions/LinearElasticity/Neumann/{}/h/expr",key);
            auto bc_neu = specs_[nl::json::json_pointer( e )].get<std::string>();
            LOG(INFO) << "BoundaryCondition Neumann : " << bc_neu << std::endl;
            l += integrate( _range = markedfaces(support(Xhv_),key), _expr = trans(expr<Dim,1>( bc_neu ))*id(u_));
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


    if (Environment::isMasterRank())
        std::cout << "***** Process loading *****" << std::endl;
    processLoading(l_);

    if (Environment::isMasterRank())
        std::cout << "***** Process materials *****" << std::endl;
    processMaterials(a_);

    if (Environment::isMasterRank())
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

    
        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_ +=  integrate( _range=elements( support(Xhv_)), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ));
        
        if (Environment::isMasterRank())
            std::cout << "***** Process contact *****" << std::endl;
        myelts_ = getContactRegion(u_);
        
        if (Environment::isMasterRank())
            std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        
        
        if (method_.compare("penalty") == 0)
            processContactPenalty(lt_, at_, myelts_, u_);
        else if (method_.compare("nitsche") == 0)
            processContactNitsche(lt_, at_, myelts_, u_);
        
        if (Environment::isMasterRank())
            std::cout << "***** Process boundary conditions *****" << std::endl;
        processBoundaryConditions(lt_, at_);

        if (Environment::isMasterRank())
            std::cout << "***** Solve *****" << std::endl;
        at_.solve( _rhs = lt_, _solution = u_ );

        if (Environment::isMasterRank())
            std::cout << "***** Export *****" << std::endl;
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
    tic();
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
    toc("init");

    /*
    if (Environment::isMasterRank())
        std::cout << "***** Process loading *****" << std::endl;
    */
   tic();
    processLoading(l_);

    /*
    if (Environment::isMasterRank())
        std::cout << "***** Process materials *****" << std::endl;
    */
    processMaterials(a_);
    toc("assSta");

    
    if (Environment::isMasterRank())
        std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        tic();
        lt_ = l_;
        at_ = a_;

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_ +=  integrate( _range=markedelements( support(Xhv_), material.get<std::string>() ), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ));
        toc("assInSta");

        tic();
        auto u_tmp =  Xhv_->element();
        u_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_));

        auto u_tmpNew = Xhv_->element();
        u_tmpNew.on( _range=elements(support(Xhv_)), _expr = idv(u_));
        

        int fixedPointIteration = 0;
        double fixedPointerror = 0.;
        toc("init"); 

        while ((fixedPointerror > fixedPointtol_) || (fixedPointIteration < 1))
        {
            tic();
            lt_tmp = lt_;
            at_tmp = at_;

            u_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_tmpNew)); ;

            /*
            if (Environment::isMasterRank())
                std::cout << "***** Process contact *****" << std::endl;
            */
            myelts_ = getContactRegion(u_tmp);
            
            /*
            if (Environment::isMasterRank())
                std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;


            if (Environment::isMasterRank())
            {
                std::cout << "Fixed point iteration : " << fixedPointIteration << std::endl;
                std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
                std::cout << "Error : " << fixedPointerror << std::endl;
            }
            */

            if (method_.compare("penalty") == 0)
                processContactPenalty(lt_tmp, at_tmp, myelts_, u_tmp);
            else if (method_.compare("nitsche") == 0)
                processContactNitsche(lt_tmp, at_tmp, myelts_, u_tmp);
            
            /*
            if (Environment::isMasterRank())
                std::cout << "***** Process boundary conditions *****" << std::endl;
            */
            processBoundaryConditions(lt_tmp, at_tmp);
            toc("AssP");

            tic();
            at_tmp.solve(_rhs = lt_tmp, _solution = u_tmpNew);
            toc("solve");

            tic();
            fixedPointerror = integrate(_range=elements(support(Xhv_)), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(support(Xhv_)),_expr=norm2(idv(u_))).evaluate()(0,0); 
            fixedPointIteration++;
            
            if (fixedPointIteration == 5)
                break;
            
            
            lt_tmp.zero();
            at_tmp.zero();
            toc("init");
        }
        
        myelts_ = getContactRegion(u_tmpNew);
        
        if (method_.compare("penalty") == 0)
            processContactPenalty(lt_, at_, myelts_, u_tmpNew);
        else if (method_.compare("nitsche") == 0)
            processContactNitsche(lt_, at_, myelts_, u_tmpNew);
                
        processBoundaryConditions(lt_, at_);

        at_.solve( _rhs = lt_, _solution = u_ );
        
        tic();
        ts_->updateFromDisp(u_);
        // Reset
        at_.zero();
        lt_.zero();
        toc("init");

        tic();
        this->exportResults(ts_->time());
        toc("export");

        
    }
}

// Run method 
template <int Dim, int Order, int OrderGeo>
void ContactDynamic<Dim, Order, OrderGeo>::run()
{
    if (Environment::isMasterRank())
        std::cout << "***** Run dynamic elasticity with unilateral contact *****" << std::endl;

    if (Environment::isMasterRank())
        std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    if (Environment::isMasterRank())
        std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    if (Environment::isMasterRank())
        std::cout <<  "***** Initialize distance g *****" << std::endl;
    tic();
    initG();
    toc("Raytracing");

    tic();
    if (Environment::isMasterRank())
        std::cout << "Init preconditioner" << std::endl;
    if constexpr(Dim == 3)
    {
        std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(backend(),qsNullSpace(Xhv_,mpl::int_<Dim>())) );
        backend()->attachNearNullSpace( myNullSpace );
    }
    toc("prec");

    if (Environment::isMasterRank())
        std::cout <<  "***** Start time loop *****" << std::endl;
    
    if ((method_.compare("penalty") == 0) || (method_.compare("nitsche") == 0) || (method_.compare("persistency") == 0) )
    {
        tic();
        this->exportResults(0);
        toc("export");
        
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
    
    contactRegion_.on( _range=elements(support(Xhv_)), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));    

    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=Xh_, _trial=Xh_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(support(Xh_)) )
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
    g_.on(_range=elements(support(Xh_)), _expr=cst(100000.));

    // Raytracing to compute distance
    using bvh_ray_type = BVHRay<Dim>;
    Eigen::VectorXd origin(Dim);
    Eigen::VectorXd dir(Dim);

    if constexpr(Dim == 2)
        dir << ddirection_[0], ddirection_[1];
    else if constexpr(Dim == 3)
        dir << ddirection_[0], ddirection_[1], ddirection_[2];

    std::string kind = (Dim==2)?"in-house":"third-party";

    auto bvh = boundingVolumeHierarchy(_range=markedfaces(mesh_, "Obs1"), _kind=kind);

    std::vector<std::size_t> faceIDs;

    BVHRaysDistributed<Dim> allrays;

    for ( auto const& theface : markedfaces( mesh_, "Wall" ) )
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
 
            
            if constexpr(Dim == 2)
            {
                auto rayIntersection = bvh->intersect( _ray = ray );
                if (!rayIntersection.empty())
                {
                    for ( auto const& rir : rayIntersection )
                    {
                        for (auto const& ldof  : Xh_->dof()->faceLocalDof( face.id() ))
                            g_[ldof.index()] = rir.distance() - tolDistance_;
                    }
                }
            }
            else 
            {
                allrays.push_back(std::move(ray));
                faceIDs.push_back(face.id());
            }            

        }
        
    }
    if constexpr(Dim == 3)
    {
            // compute intersections for allrays
        auto multiRayIntersectionResult = bvh->intersect(_ray=allrays);//,_parallel=false);
    
        for (auto const& [fid,rirs] : enumerate(multiRayIntersectionResult))
        {
            for (auto const& [rid,rir] : enumerate(rirs))
            {
                for (auto const& ldof : Xh_->dof()->faceLocalDof(faceIDs[fid]))
                    g_[ldof.index()] = rir.distance() - tolDistance_;
            }
            
        }
    }

}

// Export results
template <int Dim, int Order, int OrderGeo>
void 
ContactDynamic<Dim, Order, OrderGeo>::exportResults(double t)
{
    
    // Interpolation
    e_->step(t)->addRegions();
    myelts_ = getContactRegion(u_);

    if (t == 0)
    {
        if (Environment::isMasterRank())
            std::cout << "Export g" << std::endl;
        e_->step(t)->add( "g", g_ );
    }

    e_->step(t)->add( "displacement", u_ );
    
    auto const Id = eye<Dim,Dim>();
    auto defv = sym(gradv(u_));
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
    contactPressure_ =  project(_space=Xh_, _range=elements(support(Xh_)), _expr = cst(0.));
    contactPressure_.on(_range=myelts_, _expr = trans(expr<Dim,1>(direction_))*sigmav);
    
    auto contactDisplacement= Xh_->element();
    contactDisplacement.on( _range=elements(support(Xhv_)), _expr = (trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_)));

    e_->step(t)->add( "contactPressure", contactPressure_);
    e_->step(t)->add( "contactDisplacement", contactDisplacement );
    e_->save();

    
    
    // Save values in json file
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
    auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(contactDisplacement) );     
            
    meas_["evaluateStress"].push_back(evaluateStress(0,0));
    meas_["evaluateDisp"].push_back(evaluateDisp(0,0));

    this->writeResultsToFile("measures.json");

    /*
    int nbr = 1;
    std::ofstream ofs("outputs.csv");
    ofs << fmt::format("x, y, pressure") << std::endl;

    std::vector<double> press;
    double max;
        
    for (auto &bfaceC : markedfaces(mesh_,"Wall"))
    {
        auto & faceC = boost::unwrap_ref( bfaceC );

        auto ctx = Xh_->context();
        node_type t1(Dim);
        t1(0)=faceC.point(0).node()[0]; t1(1)=faceC.point(0).node()[1];
        ctx.add( t1 );

        auto evaluateStresstmp = evaluateFromContext( _context=ctx, _expr = idv(contactPressure_) );


        if (evaluateStresstmp(0,0)!=0)
        {
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}", faceC.point(0).node()[0], faceC.point(0).node()[1], evaluateStresstmp(0,0)) << std::endl;
            press.push_back(evaluateStresstmp(0,0));
        }

    
        nbr++;
    }

    if (press.size() != 0)
    {
        max = *min_element(press.begin(), press.end());
        std::cout << "Max : " << max << std::endl;
    }
    
    ofs.close();
    */
      
    

    
    
    

    /*
    auto sig = lambda_*trace(defv)*Id + 2*mu_*defv;
    auto J = det(Id + gradv(u_));

    double E1 = 0.5*rho_*normL2Squared(_range=elements(support(Xhv_)),_expr=idv(ts_->currentVelocity()));
    double E2 = 0.5*integrate( _range= elements( support(Xhv_) ), _expr= inner(sig,defv)).evaluate()( 0,0 );
        
    meas_["Eh1"].push_back(E1);
    meas_["Eh2"].push_back(E2);
    meas_["Eh"].push_back(E1 + E2);

    double disp = integrate( _range= myelts_, _expr= (trans(expr<Dim,1>(direction_))*idv(u_)  - idv(g_)) ).evaluate()( 0,0 );
    meas_["disp"].push_back(disp);

    
    double Lv = integrate( _range=elements(support(Xhv_)),_expr=  cst(rho_)*abs(trans(expr<Dim,1>( externalforce_ )))*idv(u_)).evaluate()( 0,0 );
    meas_["Lv"].push_back(Lv);

    
    double volume = integrate(_range=elements(support(Xhv_)), _expr = det(Id + gradv(u_))).evaluate()( 0, 0 );
    meas_["volume"].push_back(volume);
    */
    
    /*
    if ((method_.compare("penalty") == 0) || (method_.compare("persistency") == 0))
        meas_["E"].push_back(E1 + E2 + Lv);
    else if (method_.compare("nitsche") == 0)
    {
        double R1 = normL2Squared(_range= myelts_, _expr= sqrt(cst(gamma0_)/cst(gamma_)) * trans(expr<Dim,1>(direction_))*sigmav);
        double R2 = normL2Squared( _range= myelts_, _expr= sqrt(cst(gamma0_)/cst(gamma_)) * ( cst(gamma_) * ( trans(expr<Dim,1>(direction_)) *idv(u_)  - idv(g_) ) - trans(expr<Dim,1>(direction_))*sigmav ));

        meas_["R1"].push_back(R1);
        meas_["R2"].push_back(R2);
        meas_["R"].push_back((R1 - R2)/(2.*gamma0_));

        double Es = (E1 + E2) - theta_*(R1 - R2)/(2.*gamma0_);
        meas_["Es"].push_back(Es);      

        meas_["E"].push_back(Es + Lv);
    }
    */
    
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