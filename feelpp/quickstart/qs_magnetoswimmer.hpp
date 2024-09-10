#pragma once
#include "qs_active_elasticity.hpp"

namespace Feel
{

template <int Dim, int Order>
class MagnetoSwimmer
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
    MagnetoSwimmer() = default;
    MagnetoSwimmer(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    elementv_t const& u() const { return u_; }
    exporter_ptrtype const& exporter() const { return e_; }
    
    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }

    // Methods
    void initialize();
    void run();
    void timeLoop();
    void exportResults(double t);
    

private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t Xhv_, Xhv_F, Xhv_H;

    elementv_t u_, u_F, u_H;
    ts_ptrtype ts_F, ts_H;
    exporter_ptrtype e_;

    double H_;
    double E_F_, nu_F_, lambda_F_, mu_F_, rho_F_;
    double E_H_, nu_H_, lambda_H_, mu_H_, rho_H_;
    std::string disp_, source_;
    std::string activation_;
    std::vector<double> center_;
    double disp_x, disp_y;
    double Ca_, Lc_, rc_, fa_, va_;
    std::string type_;

};

// Constructor
template <int Dim, int Order>
MagnetoSwimmer<Dim, Order>::MagnetoSwimmer(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order>
void MagnetoSwimmer<Dim, Order>::initialize()
{
    // Get mesh parameters
    H_ = specs_["/Meshes/HyperElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/HyperElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    
    // Define space
    Xhv_ = Pchv<Order>( mesh_);
    Xhv_F = Pchv<Order>( mesh_, markedelements(mesh_,"flagellum") );
    Xhv_H = Pchv<Order>( mesh_, markedelements(mesh_,"head") );

    // Get flagellum parameters
    std::string matRho_F = fmt::format( "/Materials/Swimmer/parameters/rho_F/value");
    rho_F_ = std::stod(specs_[nl::json::json_pointer( matRho_F )].get<std::string>());

    std::string matE_F = fmt::format( "/Materials/Swimmer/parameters/E_F/value" );
    E_F_ = std::stod(specs_[nl::json::json_pointer( matE_F )].get<std::string>());
    
    std::string matNu_F = fmt::format( "/Materials/Swimmer/parameters/nu_F/value" );
    nu_F_ = std::stod(specs_[nl::json::json_pointer( matNu_F )].get<std::string>());
    
    lambda_F_ = E_F_*nu_F_/( (1+nu_F_)*(1-2*nu_F_) );
    mu_F_ = E_F_/(2*(1+nu_F_));

    std::string matCa = fmt::format("/Materials/Swimmer/parameters/Ca/value");
    Ca_ = std::stod(specs_[nl::json::json_pointer( matCa )].get<std::string>());

    std::string matLc = fmt::format("/Materials/Swimmer/parameters/Lc/value");
    Lc_ = std::stod(specs_[nl::json::json_pointer( matLc )].get<std::string>());

    std::string matRc = fmt::format("/Materials/Swimmer/parameters/rc/value");
    rc_ = std::stod(specs_[nl::json::json_pointer( matRc )].get<std::string>());

    std::string matFa = fmt::format("/Materials/Swimmer/parameters/fa/value");
    fa_ = std::stod(specs_[nl::json::json_pointer( matFa )].get<std::string>());

    std::string matVa = fmt::format("/Materials/Swimmer/parameters/va/value");
    va_ = std::stod(specs_[nl::json::json_pointer( matVa )].get<std::string>());
    
    std::string matType = fmt::format("/Materials/Swimmer/parameters/type/value");
    type_ = specs_[nl::json::json_pointer( matType )].get<std::string>();

    // Get head parameters
    std::string matRho_H = fmt::format( "/Materials/Swimmer/parameters/rho_H/value");
    rho_H_ = std::stod(specs_[nl::json::json_pointer( matRho_H )].get<std::string>());

    std::string matE_H = fmt::format( "/Materials/Swimmer/parameters/E_H/value" );
    E_H_ = std::stod(specs_[nl::json::json_pointer( matE_H )].get<std::string>());
    
    std::string matNu_H = fmt::format( "/Materials/Swimmer/parameters/nu_H/value" );
    nu_H_ = std::stod(specs_[nl::json::json_pointer( matNu_H )].get<std::string>());
    
    lambda_H_ = E_H_*nu_H_/( (1+nu_H_)*(1-2*nu_H_) );
    mu_H_ = E_H_/(2*(1+nu_H_));

    std::string matActivation = fmt::format( "/Models/HyperElasticity/activation" );
    activation_ = specs_[nl::json::json_pointer( matActivation )].get<std::string>();

    if (activation_.compare("disp") == 0)
    {
        std::string matDisp = fmt::format( "/Materials/Swimmer/parameters/disp/value" );
        disp_ = specs_[nl::json::json_pointer( matDisp )].get<std::string>();
    }
    else if (activation_.compare("source") == 0)
    {
        std::string matSource = fmt::format( "/Materials/Swimmer/parameters/source/value" );
        source_ = specs_[nl::json::json_pointer( matSource )].get<std::string>();

        std::string matCenter = fmt::format( "/Materials/Swimmer/parameters/center/value" );
        center_ = specs_[nl::json::json_pointer( matCenter )].get<std::vector<double>>();  

    }
    
    // Initialize exporter
    e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    
    // Initialize Newmark scheme
    bool steady = get_value(specs_, "/TimeStepping/HyperElasticity/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/HyperElasticity/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/HyperElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/HyperElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/HyperElasticity/step", std::string("0.1"))).evaluate()(0,0);
    double gamma = get_value(specs_, "/TimeStepping/HyperElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/HyperElasticity/beta", 0.25);

    // Set initial conditions
    u_ = Xhv_->element();
    u_ .on(_range=elements(support(Xhv_)), _expr = 0.*one());
    u_F = Xhv_F->element();
    u_F.on(_range=elements(support(Xhv_F)), _expr = 0.*one() );
    u_H = Xhv_H->element();
    u_H.on(_range=elements(support(Xhv_H)), _expr = 0.*one());

    ts_F = newmarkContact(Xhv_F, steady, initial_time, final_time, time_step, gamma, beta );
    ts_H = newmarkContact(Xhv_H, steady, initial_time, final_time, time_step, gamma, beta );
    
    ts_F->start();
    ts_F->initialize( u_F );

    ts_H->start();
    ts_H->initialize( u_H );

    ts_F->updateFromDisp(u_F);
    ts_H->updateFromDisp(u_H);
}

// Time loop
template <int Dim, int Order>
void MagnetoSwimmer<Dim, Order>::timeLoop()
{
    auto Id = eye<Dim,Dim>();
    
    if (activation_.compare("disp") == 0)  
    {
        std::size_t offset = 0;
        disp_x = std::stod(&disp_[1],&offset);
        disp_y = std::stod(&disp_[offset+2]);
    }


    // Initialize Flagellum
    
    
    auto a_F_ = form2( _test = Xhv_F, _trial = Xhv_F );
    auto at_F_ = form2( _test = Xhv_F, _trial = Xhv_F );
    auto l_F_ = form1( _test = Xhv_F );
    auto lt_F_ = form1( _test = Xhv_F );
    
    a_F_.zero();
    at_F_.zero();
    l_F_.zero();
    lt_F_.zero();

    auto deft_F = sym(gradt(u_F));
    auto def_F = sym(grad(u_F));

    auto sigmat_F_ = lambda_F_*trace(deft_F)*Id + 2*mu_F_*deft_F;
    a_F_ += integrate( _range = elements(support(Xhv_F)), _expr = cst(rho_F_)*inner( ts_F->polyDerivCoefficient()*idt(u_F),id( u_F ) ) + inner(sigmat_F_,def_F));
    

 
    //auto ResF = backend()->newVector(Xhv_F);
    //auto JacF = backend()->newMatrix( _test=Xhv_F, _trial=Xhv_F );
    
    // Initialize Head
    auto a_H_ = form2( _test = Xhv_H, _trial = Xhv_H );
    auto at_H_ = form2( _test = Xhv_H, _trial = Xhv_H );
    auto l_H_ = form1( _test = Xhv_H );
    auto lt_H_ = form1( _test = Xhv_H );
    
    a_H_.zero();
    at_H_.zero();
    l_H_.zero();
    lt_H_.zero();

    auto deft_H = sym(gradt(u_H));
    auto def_H = sym(grad(u_H));

    auto sigmat_H_ = lambda_H_*trace(deft_H)*Id + 2*mu_H_*deft_H;
    a_H_ += integrate( _range = elements(support(Xhv_H)), _expr = cst(rho_H_)*inner( ts_H->polyDerivCoefficient()*idt(u_H),id( u_H ) ) + inner(sigmat_H_,def_H));

    if (activation_.compare("source") == 0)
    {
        std::cout << "Activation : source" << std::endl;
        l_H_ += integrate( _range = elements( support(Xhv_H)), _expr = trans(expr<Dim,1>( source_ ))*id(u_H));
    }
        

    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_F->timeInitial(),ts_F->timeFinal(), ts_F->timeStep()) << std::endl;
    
    for ( ts_F->start(); ts_F->isFinished()==false; ts_F->next(u_F), ts_H->next(u_H) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_F->time(),ts_F->timeFinal()) << std::endl;

        

        // Flagellum
        /*
        auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto u = Xhv_F->element();
            u = *X;
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_F_*trace(Ev)*Id + 2*mu_F_*Ev;
            
            auto dF = gradt(u);
            auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            auto dS = lambda_F_*trace(dE)*Id + 2*mu_F_*dE;
            
            auto ea = vec(cst(0.), cst(1.));
            auto eaea = ea*trans(ea);
            
            auto a = form2( _test=Xhv_F, _trial=Xhv_F, _matrix=J );

            a = integrate( _range=elements(support(Xhv_F)),
                           _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
            
            a += integrate( _range=elements(support(Xhv_F)),
                            _expr= cst(rho_F_)*inner( ts_F->polyDerivCoefficient()*idt(u),id( u ) ) );

            
            if (type_.compare("bending") == 0)
            {
                std::cout << "Bending" << std::endl;
                auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_F->time());
                a += integrate(_range=elements( support(Xhv_F) ), _expr= - inner(dF*sigma_a*Px()*eaea, grad(u) ) );
            }    
            else if (type_.compare("flapping") == 0)
            {
                std::cout << "Flapping" << std::endl;
                auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_F->time()));
                a += integrate(_range=elements( support(Xhv_F) ), _expr= - inner(dF*sigma_a*Px()*eaea, grad(u) ) );
            }
            
            
            auto RR = backend()->newVector( Xhv_F );

            if (activation_.compare("disp") == 0)  
            {
                std::cout << "Activation : displacement" << std::endl;
                auto dx = disp_x*ts_F->time();
                auto dy = disp_y*ts_F->time();
                a += on(_range=markedfaces(mesh_,"motor"), _rhs=RR, _element=u, _expr= vec(cst(dx),cst(dy)));
            }  
            else if (activation_.compare("source") == 0)  
            {
                std::cout << "Activation : source" << std::endl;

                // On détermine le déplacement u_H
                auto ctx = Xhv_H->context();
                node_type t(Dim);
                t(0)=center_[0]; t(1)=center_[1];  
                ctx.add( t );

                auto ifv_uH = evaluateFromContext( _context=ctx, _expr= idv(u_H) ); 
                std::cout << "Disp_x : " << ifv_uH(0,0)  << " Disp_y : " << ifv_uH(1,0) << std::endl;

                a += on(_range=markedfaces(mesh_,"motor"), _rhs=RR, _element=u, _expr =  vec(cst(ifv_uH(0,0)),cst(ifv_uH(1,0))));
            }

        };

        auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Xhv_F->element();
            u = *X;
            
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_F_*trace(Ev)*Id + 2*mu_F_*Ev;

            auto ea = vec(cst(1.), cst(0.));
            auto eaea = ea*trans(ea);
            auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_F->time());
        

            auto r = form1( _test=Xhv_F, _vector=R );
            r = integrate( _range=elements(support(Xhv_F)),
                           _expr= inner( val(Fv*Sv) , grad(u) ) );

            r += integrate( _range=elements(support(Xhv_F)),
                            _expr= cst(rho_F_)*inner( ts_F->polyDerivCoefficient()*idv(u) -idv(ts_F->polyDeriv()),id( u ) ) );
            
            
            if (type_.compare("bending") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_F->time());
                r += integrate(_range=elements( support(Xhv_F)), _expr= - inner( val(Fv)*sigma_a*Px()*eaea,grad( u )));
            }    
            else if (type_.compare("flapping") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_F->time()));
                r += integrate(_range=elements( support(Xhv_F)), _expr= - inner( val(Fv)*sigma_a*Px()*eaea,grad( u )));
            }
            

            R->close();
            auto temp = Xhv_F->element();
            temp = *R;

            if (activation_.compare("disp") == 0)  
            {
                std::cout << "Activation : displacement" << std::endl;
                auto dx = disp_x*ts_F->time();
                auto dy = disp_y*ts_F->time();
                temp.on( _range=markedfaces(mesh_,"motor"),_expr=vec(cst(dx),cst(dy)) );
            }  
            else if (activation_.compare("source") == 0)  
            {
                std::cout << "Activation : source" << std::endl;

                // On détermine le déplacement u_H
                auto ctx = Xhv_H->context();
                node_type t(Dim);
                t(0)=center_[0]; t(1)=center_[1];  
                ctx.add( t );

                auto ifv_uH = evaluateFromContext( _context=ctx, _expr= idv(u_H) ); 
                std::cout << "Disp_x : " << ifv_uH(0,0)  << " Disp_y : " << ifv_uH(1,0) << std::endl;

                temp.on(_range=markedfaces(mesh_,"motor"),_expr =  vec(cst(ifv_uH(0,0)),cst(ifv_uH(1,0))));
            }

            *R = temp;
        };

        
        if (activation_.compare("disp") == 0)  
        {
            std::cout << "Activation : displacement" << std::endl;
            auto dx = disp_x*ts_F->time();
            auto dy = disp_y*ts_F->time();
            u_F.on( _range=markedfaces(mesh_,"motor"),_expr=vec(cst(dx),cst(dy)) );
        }  
        else if (activation_.compare("source") == 0)  
        {
            std::cout << "Activation : source" << std::endl;

            // On détermine le déplacement u_H
            auto ctx = Xhv_H->context();
            node_type t(Dim);
            t(0)=center_[0]; t(1)=center_[1];  
            ctx.add( t );

            auto ifv_uH = evaluateFromContext( _context=ctx, _expr= idv(u_H) ); 
            std::cout << "Disp_x : " << ifv_uH(0,0)  << " Disp_y : " << ifv_uH(1,0) << std::endl;

            u_F.on(_range=markedfaces(mesh_,"motor"),_expr =  vec(cst(ifv_uH(0,0)),cst(ifv_uH(1,0))));
        }
        
        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u_F,_jacobian=JacF,_residual=ResF);
        ts_F->updateFromDisp(u_F);
        */


        
        lt_F_ = l_F_;
        at_F_ = a_F_;
    
        std::cout << "Solve flagellum" << std::endl;
        lt_F_ +=  integrate( _range=elements( support(Xhv_F)), _expr= cst(rho_F_)*inner( idv(ts_F->polyDeriv()),id( u_F ) ));

        if (activation_.compare("disp") == 0)  
        {
            std::cout << "Activation : displacement" << std::endl;
            auto dx = disp_x*ts_F->time();
            auto dy = disp_y*ts_F->time();
            at_F_ += on(_range=markedfaces(mesh_,"motor"), _rhs=lt_F_, _element=u_F, _expr=  vec(cst(dx),cst(dy)));
        }  
        else if (activation_.compare("source") == 0)  
        {
            std::cout << "Activation : source" << std::endl;

            // On détermine le déplacement u_H
            auto ctx = Xhv_H->context();
            node_type t(Dim);
            t(0)=center_[0]; t(1)=center_[1];  
            ctx.add( t );

            auto ifv_uH = evaluateFromContext( _context=ctx, _expr= idv(u_H) ); 
            std::cout << "Disp_x : " << ifv_uH(0,0)  << " Disp_y : " << ifv_uH(1,0) << std::endl;

            at_F_ += on(_range=markedfaces(mesh_,"motor"), _rhs=lt_F_, _element=u_F, _expr =  vec(cst(ifv_uH(0,0)),cst(ifv_uH(1,0))));
        } 

        at_F_.solve( _rhs = lt_F_, _solution = u_F, _rebuild=true );
        
        

        // Head
        lt_H_ = l_H_;
        at_H_ = a_H_;

        std::cout << "Solve head" << std::endl;
        lt_H_ +=  integrate( _range=elements( support(Xhv_H)), _expr= cst(rho_H_)*inner( idv(ts_H->polyDeriv()),id( u_H ) ));

        if (activation_.compare("disp") == 0)  
        {
            std::cout << "Activation : displacement" << std::endl;  
            auto dx = disp_x*ts_F->time();
            auto dy = disp_y*ts_F->time();
            std::cout << "Disp_x : " << dx  << " Disp_y : " << dy << std::endl;
            //at_H_ += on(_range=markedfaces(mesh_,"headB"), _rhs=lt_H_, _element=u_H, _expr= expr<Dim,1>( disp_ ));
            at_H_ += on(_range=markedfaces(mesh_,"headB"), _rhs=lt_H_, _element=u_H, _expr= vec(cst(dx),cst(dy)));
        }

        at_H_.solve( _rhs = lt_H_, _solution = u_H, _rebuild=true );
  
        ts_H->updateFromDisp(u_H);

        

        // Total displacement
        std::cout << "Solve swimmer" << std::endl;
        u_.on(_range=elements(support(Xhv_)),_expr=idv(u_F) + idv(u_H));

        this->exportResults(ts_F->time());

        // Reset
        //at_F_.zero();
        //lt_F_.zero();
        at_H_.zero();
        lt_H_.zero();
    }    
}

// Run method 
template <int Dim, int Order>
void MagnetoSwimmer<Dim, Order>::run()
{
    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    std::cout <<  "***** Start time loop *****" << std::endl;
    this->exportResults(0);
    timeLoop();
}


// Export results
template <int Dim, int Order>
void 
MagnetoSwimmer<Dim, Order>::exportResults(double t)
{
    e_->step(t)->addRegions();
    e_->step(t)->add( "displacement", u_ );
    e_->step(t)->add( "displacementH", u_H );
    e_->step(t)->add( "displacementF", u_F );
    e_->save();
}

} 