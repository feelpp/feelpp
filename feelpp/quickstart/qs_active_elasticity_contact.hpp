#pragma once
#include "qs_elasticity_contact.hpp"

namespace Feel
{

template <int Dim, int Order>
class ActiveContact
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;
    using spacev_t = Pchv_type<mesh_t, Order>;
    using space_t = Pch_type<mesh_t, Order>;
    using spacev_ptr_t = Pchv_ptrtype<mesh_t, Order>; 
    using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
    using elementv_t = typename spacev_t::element_type;
    using element_t = typename space_t::element_type;
    using form2_type = form2_t<spacev_t,spacev_t>; 
    using form1_type = form1_t<spacev_t>; 
    using ts_ptrtype = std::shared_ptr<NewmarkContact<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

    // Constructors
    ActiveContact() = default;
    ActiveContact(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    elementv_t const& u() const { return u_; }
    exporter_ptrtype const& exporter() const { return e_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }

    // Methods
    void initialize();
    void processBoundaryConditions(form1_type& l, form2_type& a);
    void run();
    void timeLoop();
    void exportResults(double t);


private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t Xhv_;
    space_ptr_t Xh_;

    elementv_t u_;

    ts_ptrtype ts_;
    exporter_ptrtype e_;

    double Ca_, Lc_, rc_, fa_;
    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;
};

// Constructor
template <int Dim, int Order>
ActiveContact<Dim, Order>::ActiveContact(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order>
void ActiveContact<Dim, Order>::initialize()
{
    // Mesh
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    Xhv_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) ); 
    Xh_ = Pch<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) ); 

    // Structure parameters
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
    
    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
    
    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));

    std::string matCa = fmt::format("/Materials/Caoutchouc/parameters/Ca/value");
    Ca_ = std::stod(specs_[nl::json::json_pointer( matCa )].get<std::string>());

    std::string matLc = fmt::format("/Materials/Caoutchouc/parameters/Lc/value");
    Lc_ = std::stod(specs_[nl::json::json_pointer( matLc )].get<std::string>());

    std::string matRc = fmt::format("/Materials/Caoutchouc/parameters/rc/value");
    rc_ = std::stod(specs_[nl::json::json_pointer( matRc )].get<std::string>());

    std::string matFa = fmt::format("/Materials/Caoutchouc/parameters/fa/value");
    fa_ = std::stod(specs_[nl::json::json_pointer( matFa )].get<std::string>());
    
    // External force
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
    std::cout << "External force : " << externalforce_ << std::endl;

    // Exporter
    e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    
    // Newmark scheme
    bool steady = get_value(specs_, "/TimeStepping/LinearElasticity/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/LinearElasticity/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    double gamma = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);

    // Initial conditions
    u_ = Xhv_->element();
    auto u0_ = Xhv_->element();

    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    auto init_displ = expr<Dim,1>(get_value(specs_, "/InitialConditions/LinearElasticity/displacement/expr", default_displ ));
    u0_.on(_range=elements(support(Xhv_)), _expr=init_displ);
    
    ts_ = newmarkContact(Xhv_, steady, initial_time, final_time, time_step, gamma, beta );
    ts_->start();
    ts_->initialize( u0_ );
    u_ = u0_;
    
    ts_->updateFromDisp(u_);
}

// Process boundary conditions
template <int Dim, int Order>
void ActiveContact<Dim, Order>::processBoundaryConditions(form1_type& l, form2_type& a)
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
            a+=on(_range=markedfaces(support(Xhv_),key), _rhs=l, _element=u_, _expr=expr<Dim,1>( bc_dir ) );
        }
    }
}


// Time loop
template <int Dim, int Order>
void ActiveContact<Dim, Order>::timeLoop()
{
    // Initialize linear and bilinear forms
    auto Res = backend()->newVector(Xhv_);
    auto Jac = backend()->newMatrix( _test=Xhv_, _trial=Xhv_ );
    auto Id = eye<Dim,Dim>();


    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        std::cout << "time : " << ts_->time() << std::endl;
        auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
        std::cout << "sigma_a : " << sigma_a << std::endl;

        auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto u = Xhv_->element();
            u = *X;
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;
            auto dF = gradt(u);
            auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            auto dS = lambda_*trace(dE)*Id + 2*mu_*dE;
            
            auto ea = vec(cst(1.), cst(0.));
            auto eaea = ea*trans(ea);
            //auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*ts_->time())*(-Px());
            auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
            

            
            auto a = form2( _test=Xhv_, _trial=Xhv_, _matrix=J );

            a = integrate( _range=elements(support(Xhv_)),
                           _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
            a += integrate( _range=elements(support(Xhv_)),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );

            a += integrate(_range=elements( support(Xhv_) ),
                    _expr=inner( - dF*sigma_a*Px()*eaea, grad(u) ) );

            auto RR = backend()->newVector( Xhv_ );
            a += on( _range=markedfaces(mesh_,"dirichlet"),
                     _element=u, _rhs=RR,
                     _expr=zero<Dim,1>() );
        };
    
        auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Xhv_->element();
            u = *X;
            
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

            auto ea = vec(cst(1.), cst(0.));
            auto eaea = ea*trans(ea);
            auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
        

            auto r = form1( _test=Xhv_, _vector=R );
            r = integrate( _range=elements(support(Xhv_)),
                           _expr= inner( val(Fv*Sv) , grad(u) ) );
            r += integrate( _range=elements(support(Xhv_)),
                            _expr= - trans( expr<Dim, 1>( externalforce_ ) )*id( u )  );
            r += integrate( _range=elements(support(Xhv_)),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
            r += integrate(_range=elements( support(Xhv_)),
                    _expr= inner( val(Fv)*sigma_a*Px()*eaea,grad( u )));

            R->close();
            auto temp = Xhv_->element();
            temp = *R;
            temp.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
            *R = temp;
        };

        u_.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u_,_jacobian=Jac,_residual=Res );


        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

    }    
}

// Run method 
template <int Dim, int Order>
void ActiveContact<Dim, Order>::run()
{
    std::cout << "***** Run active elasticity with contact *****" << std::endl;

    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    this->initialize();

    std::cout <<  "***** Start time loop *****" << std::endl;
    this->exportResults(0);
    this->timeLoop();
}


// Export results
template <int Dim, int Order>
void 
ActiveContact<Dim, Order>::exportResults(double t)
{
    // Define exports
    e_->step(t)->addRegions();
    e_->step(t)->add( "displacement", u_ );
    e_->step(t)->add( "velocity", ts_->currentVelocity() );
    e_->save();
}

} 