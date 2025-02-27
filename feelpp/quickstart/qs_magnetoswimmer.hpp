#pragma once
#include "qs_active_elasticity.hpp"
#include <feel/feelts/bdf.hpp>
#include <cmath> 

namespace Feel
{

template <int Dim>
class MagnetoSwimmer
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;

    using spacer_t = Pch_type<mesh_t, 0>;
    using spacev_t = Pchv_type<mesh_t, 1>;

    using spacev_ptr_t = Pchv_ptrtype<mesh_t, 1>; 
    using spacer_ptr_t = Pch_ptrtype<mesh_t, 0>;

    using elementr_t = typename spacer_t::element_type;
    using elementv_t = typename spacev_t::element_type;

    using form2_type = form2_t<spacev_t,spacev_t>; 
    using form1_type = form1_t<spacev_t>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 
    using ts_ptrtype = std::shared_ptr<Newmark<spacev_t>>;


    std::vector<double> velocities; //Stores the velocities

    // Constructors
    MagnetoSwimmer() = default;
    MagnetoSwimmer(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    
    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }

    // Methods
    void initialize();
    void run_rigid();
    void run_elastic();
    void timeloop_rigid();
    void timeloop_elastic();
    void exportResults_rigid(double t);
    void exportResults_elastic(double t);
    

private:
    nl::json specs_; // json

    double H_; //Mesh step
    std::shared_ptr<mesh_t> mesh_; //Mesh 

    spacev_ptr_t Vh_disp_; //Displacement space
    spacev_ptr_t Vh_r_; 
    spacev_ptr_t Vh_;

    double rho_; // density
    double lambda_;
    double mu_;
    double lambda_tail;
    double mu_tail;
    double rho_tail;
    double rescale_;
    double freq_;
    double bx_;
    double by_;
    double initial_time, final_time, time_step;
    std::string externalforce_;

    exporter_ptrtype e_r; // exporter
    exporter_ptrtype e_e; // exporter
    ts_ptrtype ts_;

    double omega_; //Angular speed
    double theta_; //Angle
    elementv_t u_r; //Rotation displacement
    elementv_t u_e; //Rotation displacement
    elementv_t u_tot;

};

// Constructor
template <int Dim>
MagnetoSwimmer<Dim>::MagnetoSwimmer(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim>
void MagnetoSwimmer<Dim>::initialize()
{
    // Get mesh parameters
    H_ = specs_["/Meshes/Rigid/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/Rigid/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    // Define space
    Vh_disp_ = Pchv<1>( mesh_, markedelements( mesh_, "tail" ) );
    //Vh_r_ =  Pchv<1>( mesh_, markedelements( mesh_, "head" ) );
    Vh_ =  Pchv<1>( mesh_ );

    
    // Get density
    std::string matRho = fmt::format( "/Materials/MagnetoObject/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

    std::string matRho_tail = fmt::format( "/Materials/MagnetoObject/parameters/rho_tail/value");
    rho_tail = std::stod(specs_[nl::json::json_pointer( matRho_tail )].get<std::string>());

    std::string matExt = fmt::format( "/Materials/MagnetoObject/parameters/ext/value");
    externalforce_ = specs_[nl::json::json_pointer( matExt )].get<std::string>();


    std::string matE = fmt::format( "/Materials/MagnetoObject/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
    
    std::string matNu = fmt::format( "/Materials/MagnetoObject/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
    
    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));

    std::string matE_tail = fmt::format( "/Materials/MagnetoObject/parameters/E_tail/value" );
    double E_tail = std::stod(specs_[nl::json::json_pointer( matE_tail )].get<std::string>());
    
    std::string matNu_tail = fmt::format( "/Materials/MagnetoObject/parameters/nu_tail/value" );
    double nu_tail = std::stod(specs_[nl::json::json_pointer( matNu_tail )].get<std::string>());
    
    lambda_tail = E_tail*nu_tail/( (1+nu_tail)*(1-2*nu_tail) );
    mu_tail = E_tail/(2*(1+nu_tail));

    std::cout << "head : " << rho_ << ", " << lambda_ << ", " << mu_ << std::endl;
    std::cout << "tail : " << rho_tail << ", " << lambda_tail << ", " << mu_tail << std::endl;
    
    // Get param
    std::string matRescale = fmt::format( "/Materials/MagnetoObject/parameters/rescale/value");
    rescale_ = std::stod(specs_[nl::json::json_pointer( matRescale )].get<std::string>());

    std::string matFreq = fmt::format( "/Materials/MagnetoObject/parameters/freq/value");
    freq_ = std::stod(specs_[nl::json::json_pointer( matFreq )].get<std::string>());

    std::string matBx = fmt::format( "/Materials/MagnetoObject/parameters/bx/value");
    bx_ = std::stod(specs_[nl::json::json_pointer( matBx )].get<std::string>());

    std::string matBy = fmt::format( "/Materials/MagnetoObject/parameters/by/value");
    by_ = std::stod(specs_[nl::json::json_pointer( matBy )].get<std::string>());

    // Initialize exporter
    e_e = Feel::exporter(_mesh = mesh_, _name = "elastic", _geo="change" );
    
    // Time scheme
    initial_time = get_value(specs_, "/TimeStepping/Rigid/start", 0.0);
    final_time = get_value(specs_, "/TimeStepping/Rigid/end", 1.0);
    time_step = expr(get_value(specs_, "/TimeStepping/Rigid/step", std::string("0.1"))).evaluate()(0,0);
    
    // Set initial conditions
    omega_ = 0.;
    theta_ = 0.;

    //u_r = Vh_r_->element();
    //u_r.on(_range=elements(support(Vh_r_)), _expr=vec(cst(0.),cst(0.)));

    u_e =  Vh_disp_->element();
    u_e.on(_range=elements(support(Vh_disp_)), _expr=vec(cst(0.),cst(0.)));

    u_tot = Vh_->element();
    u_tot.on(_range=elements(mesh_), _expr=vec(cst(0.),cst(0.)));
    
    double beta = 0.25;
    double gamma = 0.5;

    ts_ =  newmark(_space = Vh_disp_, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _gamma=gamma, _beta=beta );
    ts_->start();
    ts_->initialize( u_e );    
    ts_->updateFromDisp(u_e);
}


template<int Dim>
void MagnetoSwimmer<Dim>::timeloop_rigid()
{
    // External magnetic field
    double bx = bx_/rescale_;
    double by = by_/rescale_;
    //double By = 0;
    
    // Magnetic momentum 
    double Mx = std::cos(0);
    double My = std::sin(0);

    // Compute moment of inertia J
    auto massCenter = mean( _range = elements(mesh_), _expr = P());
    auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
    std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

    auto momentOfInertia = integrate(_range=elements(mesh_),_expr=cst(rho_)*( inner(P()-massCenterVec) ) ).evaluate()(0,0);
    std::cout << "momentOfInertia : " << momentOfInertia << std::endl;
   
    // Time 
    double time = time_step;

    // Output
    std::ofstream ofs("res_rigid.csv");
    ofs << fmt::format("theta") << std::endl;
    ofs << fmt::format( "{:.6f}",theta_) << std::endl;

    while (time < final_time)
    {
        std::cout << "Time : " << time << std::endl;

        //By = by * std::sin(2*freq_*M_PI*time);

        //Compute torque
        //Mx = std::cos(theta_);
        //My = std::sin(theta_);


        //double Tm = (Mx*By - My*bx);
        double Tm = (Mx*by - My*bx);
        std::cout << "Torque : " << Tm << std::endl;

        //Update
        omega_ = omega_ + time_step * (Tm/momentOfInertia);
        //std::cout << "Angular velocity : " << omega_ << std::endl;
        theta_ = theta_ + omega_*time_step;
        //std::cout << "Rotation angle : " << theta_ << std::endl;
        
        ofs << fmt::format( "{:.6f}",theta_) << std::endl;

        auto rot = vec(
            cos(theta_) * (Px() - massCenter(0,0)) - sin(theta_) * (Py() - massCenter(1,0)) - Px() + massCenter(0,0),
            sin(theta_) * (Px() - massCenter(0,0)) + cos(theta_) * (Py() - massCenter(1,0)) - Py() + massCenter(1,0)
        );

        u_r.on(_range = elements(mesh_), _expr = rot);

        
        // Export results
        this->exportResults_rigid(time);
        time += time_step;
    }
}

// Time loop
/*
template <int Dim>
void MagnetoSwimmer<Dim>::timeloop_elastic()
{
    // External magnetic field
    double bx = bx_/rescale_;
    double by = by_/rescale_;
    double By = 0;
    double hx = 0;
    
    // Magnetic momentum 
    double Mx = std::cos(0);
    double My = std::sin(0);

    // Compute moment of inertia J
    auto massCenter = mean( _range = elements(mesh_), _expr = P());
    auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
    std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

    auto momentOfInertia = integrate(_range=elements(mesh_),_expr=cst(rho_)*( inner(P()-massCenterVec) ) ).evaluate()(0,0);
    std::cout << "momentOfInertia : " << momentOfInertia << std::endl;
   

    auto Id = eye<2,2>();
    auto Res = backend()->newVector(Vh_disp_);
    auto Jac = backend()->newMatrix( _test=Vh_disp_, _trial=Vh_disp_ );
    
    int iter = 0;

    std::ofstream ofs("res_elastic_hyper.csv");
    ofs << fmt::format("time,theta,error") << std::endl;
    

    while (time_step * iter < final_time)
    {
        std::cout << "time : " << time_step * iter << std::endl;
        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////

        //hx = 0.75 * (Mx*by - My*bx);
        hx = (Mx*by - My*bx)/0.75;
        std::cout << "hx : " << hx << std::endl;

        std::cout << "Solve elasticity" << std::endl;
        // Initialize linear and bilinear forms
    
        auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto u = Vh_disp_->element();
            u = *X;
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;
            auto Sv_tail = lambda_tail*trace(Ev)*Id + 2*mu_tail*Ev;

            
            auto dF = gradt(u);
            auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            auto dS = lambda_*trace(dE)*Id + 2*mu_*dE;
            auto dS_tail = lambda_tail*trace(dE)*Id + 2*mu_tail*dE;

            
            auto a = form2( _test=Vh_disp_, _trial=Vh_disp_, _matrix=J );

            a = integrate( _range=markedelements(mesh_,"head"),
                           _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
            a = integrate( _range=markedelements(mesh_,"tail"),
                           _expr= inner( dF*val(Sv_tail) + val(Fv)*dS_tail , grad(u) ) );
            a += integrate( _range=markedelements(mesh_,"head"),
                           _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );
            a += integrate( _range=markedelements(mesh_,"tail"),
                            _expr= cst(rho_tail)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );

        };
    
        auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Vh_disp_->element();
            u = *X;
            
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;
            auto Sv_tail = lambda_tail*trace(Ev)*Id + 2*mu_tail*Ev;


            auto r = form1( _test=Vh_disp_, _vector=R );
            r = integrate( _range=markedelements(mesh_,"head"),
                           _expr= inner( val(Fv*Sv) , grad(u) ) );
            r = integrate( _range=markedelements(mesh_,"tail"),
                           _expr= inner( val(Fv*Sv_tail) , grad(u) ) );

            r += integrate( _range=markedelements(mesh_,"head"),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
            r += integrate( _range=markedelements(mesh_,"tail"),
                            _expr= cst(rho_tail)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
            
            r += integrate( _range = markedfaces(mesh_, "Upper"), _expr = trans(vec(cst(0.),cst(-hx)))*id(u));
            r += integrate( _range = markedfaces(mesh_, "Lower"), _expr = trans(vec(cst(0.),cst(hx)))*id(u));

            R->close();
            
        };
        
        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u_e,_jacobian=Jac,_residual=Res );


        ts_->updateFromDisp(u_e);
        ts_->next(u_e);
        
        std::cout << "Export" << std::endl;
        this->exportResults_elastic(time_step*iter);

        //Update
        omega_ = omega_ + time_step * (0.75*hx/momentOfInertia);
        //std::cout << "Angular velocity : " << omega_ << std::endl;
        theta_ = theta_ + omega_*time_step;
        //std::cout << "Rotation angle : " << theta_ << std::endl;
        
        auto rot = vec(
            cos(theta_) * (Px() - massCenter(0,0)) - sin(theta_) * (Py() - massCenter(1,0)) - Px() + massCenter(0,0),
            sin(theta_) * (Px() - massCenter(0,0)) + cos(theta_) * (Py() - massCenter(1,0)) - Py() + massCenter(1,0)
        );

        u_r.on(_range = elements(mesh_), _expr = rot);

        
        // Export results
        this->exportResults_rigid(time_step*iter);

        // Check
        auto l2err = normL2( _range=elements(mesh_), _expr=idv(u_r) - idv(u_e) );
        std::cout << "Error : " << l2err << std::endl;

        ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",time_step*iter, theta_, l2err) << std::endl;
        iter++;
    } 
       
}
*/

template <int Dim>
void MagnetoSwimmer<Dim>::timeloop_elastic()
{
    // External magnetic field
    double bx = bx_/rescale_;
    double by = by_/rescale_;
    double By = 0;
    double hx = 0;
    
    // Magnetic momentum 
    double Mx = std::cos(0);
    double My = std::sin(0);

    // Compute moment of inertia J
    //auto massCenter = mean( _range = markedelements(mesh_,"head"), _expr = P());
    //auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
    //std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

    //auto momentOfInertia = integrate(_range = markedelements(mesh_,"head"),_expr=cst(rho_)*( inner(P()-massCenterVec) ) ).evaluate()(0,0);
    double momentOfInertia = 0.15625;
    std::cout << "momentOfInertia : " << momentOfInertia << std::endl;
    

    auto Id = eye<2,2>();
    auto Res = backend()->newVector(Vh_disp_);
    auto Jac = backend()->newMatrix( _test=Vh_disp_, _trial=Vh_disp_ );
    
    int iter = 0;

    std::ofstream ofs("res_elastic_hyper.csv");
    ofs << fmt::format("time,theta") << std::endl;
    

    while (time_step * iter < final_time)
    {
        std::cout << "time : " << time_step * iter << std::endl;
        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        hx = (Mx*by - My*bx);
        
        //Update
        omega_ = omega_ + time_step * (hx/momentOfInertia);
        std::cout << "Angular velocity : " << omega_ << std::endl;
        theta_ = theta_ + omega_*time_step;
        std::cout << "Rotation angle : " << theta_ << std::endl;
                
        auto rot = vec(
            cos(theta_) * (Px() -  cst(0.25)) - sin(theta_) * (Py() - cst(0.75)) - Px() + cst(0.25),
            sin(theta_) * (Px() -  cst(0.25)) + cos(theta_) * (Py() - cst(0.75)) - Py() + cst(0.75)
        );

        //u_r.on(_range = elements(support(Vh_r_)), _expr = rot);
    

        std::cout << "Solve elasticity" << std::endl;
        // Initialize linear and bilinear forms
    
        auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto u = Vh_disp_->element();
            u = *X;
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv_tail = lambda_tail*trace(Ev)*Id + 2*mu_tail*Ev;

            
            auto dF = gradt(u);
            auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
            auto dS_tail = lambda_tail*trace(dE)*Id + 2*mu_tail*dE;

            
            auto a = form2( _test=Vh_disp_, _trial=Vh_disp_, _matrix=J );


            a = integrate( _range=elements(support(Vh_disp_)),
                           _expr= inner( dF*val(Sv_tail) + val(Fv)*dS_tail , grad(u) ) );

            a += integrate( _range=elements(support(Vh_disp_)),
                            _expr= cst(rho_tail)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );
            
            auto RR = backend()->newVector( Vh_disp_ );
            a += on( _range=markedfaces(mesh_,"fsi"),
                        _element=u, _rhs=RR,
                        _expr=cst(0)*rot );

        };
    
        auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Vh_disp_->element();
            u = *X;
            
            auto Fv = Id + gradv(u);
            auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
            auto Sv_tail = lambda_tail*trace(Ev)*Id + 2*mu_tail*Ev;


            auto r = form1( _test=Vh_disp_, _vector=R );
  
            r = integrate( _range=elements(support(Vh_disp_)),
                           _expr= inner( val(Fv*Sv_tail) , grad(u) ) );

            r += integrate( _range=elements(support(Vh_disp_)),
                            _expr= cst(rho_tail)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
            r += integrate( _range = elements(support(Vh_disp_)), _expr = -trans(expr<Dim,1>( externalforce_ ))*id( u ));

            R->close();
            auto temp = Vh_disp_->element();
            temp = *R;
            temp.on( _range=markedfaces(mesh_,"fsi"),_expr=cst(0)*rot );
                            
            *R = temp;
            
        };
        u_e.zero();
        u_e.on( _range=markedfaces(mesh_,"fsi"),_expr=rot );

        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u_e,_jacobian=Jac,_residual=Res );


        ts_->updateFromDisp(u_e);
        ts_->next(u_e);

        //u_tot.on(_range= elements(support(Vh_r_)), _expr=idv(u_r));
        //u_tot.on(_range= elements(support(Vh_disp_)), _expr=idv(u_e));

        
        std::cout << "Export" << std::endl;
        this->exportResults_elastic(time_step*iter);

        ofs << fmt::format( "{:.6f}, {:.6f}",time_step*iter, theta_) << std::endl;
        iter++;
    } 
       
}

/*
// Time loop
template <int Dim>
void MagnetoSwimmer<Dim>::timeloop_elastic()
{
    // External magnetic field
    double bx = bx_/rescale_;
    double by = by_/rescale_;
    double By = 0;
    double hx = 0;
    
    // Magnetic momentum 
    double Mx = std::cos(0);
    double My = std::sin(0);

    // Compute moment of inertia J
    auto massCenter = mean( _range = elements(mesh_), _expr = P());
    auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
    std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

    auto momentOfInertia = integrate(_range=elements(mesh_),_expr=cst(rho_)*( inner(P()-massCenterVec) ) ).evaluate()(0,0);
    std::cout << "momentOfInertia : " << momentOfInertia << std::endl;
   

    // define fields
    auto utot = Vh_disp_->element(); 
    utot.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));  

    auto utot_e = Vh_disp_->element(); 
    utot_e.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));  

    auto dt_u_tot = Vh_disp_->element();
    dt_u_tot.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));   

    auto dtt_u_tot = Vh_disp_->element();
    dtt_u_tot.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));  

    auto dt_u_old = Vh_disp_->element(); 
    dt_u_old.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));  

    auto dtt_u_old = Vh_disp_->element();
    dtt_u_old.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));   

    auto dtt_u_old2 = Vh_disp_->element();
    dtt_u_old2.on(_range=elements( mesh_init ), _expr= vec(cst(0.),cst(0.)));  

    auto Id = eye<2,2>();
    double beta = 0.25;
    double gamma = 0.5;

    int iter = 0;

    std::ofstream ofs("res_elastic.csv");
    ofs << fmt::format("time,theta,error") << std::endl;
    

    while (time_step * iter < final_time)
    {
        std::cout << "time : " << time_step * iter << std::endl;
        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////

        //hx = 0.75 * (Mx*by - My*bx);
        hx = (Mx*by - My*bx)/0.75;
        std::cout << "hx : " << hx << std::endl;

        std::cout << "Solve elasticity" << std::endl;
        
        Vh_disp_current = Pchv<1>( mesh_);
        u_e =  Vh_disp_current->element();
        u_e.on(_range=elements(mesh_), _expr=vec(cst(0.),cst(0.)));
    
        auto eps_curr = sym(gradt(u_e));
        auto sigma_curr = lambda_*trace(eps_curr)*Id + 2*mu_*eps_curr;

        auto eps_tot = sym(gradv(utot)); 
        auto F = Id + gradv(utot); 
        auto J = det(F); 
        auto Js = J * sqrt(trans(N())*( trans(F) * F) * N()); 
        auto sigma_tot = lambda_*trace(eps_tot)*Id + 2*mu_*eps_tot;

        auto eps_mix = 0.5 * (gradt(u_e)*gradv(utot) + trans(gradv(utot))*trans(gradt(u_e)));
        auto sigma_mix = lambda_ * trace(eps_mix) * Id + 2 * mu_ * eps_mix;
       
        auto a = form2( _trial=Vh_disp_current, _test=Vh_disp_current);
        auto l = form1( _test=Vh_disp_current );

        a.zero();
        l.zero();

        a += integrate(_range = elements(mesh_), _expr = cst(1.)/J * inner( sigma_curr * trans(F), grad(u_e) ));
        a += integrate(_range = elements(mesh_), _expr = cst(1.)/J * inner( sigma_mix * trans(F), grad(u_e) ));
        a += integrate(_range = elements(mesh_), _expr= cst(rho_)/J * inner( cst(1.0)/(cst(beta)*std::pow(time_step,2)) * idt( u_e ),id( u_e ) ) );
           
        l += integrate( _range = elements(mesh_), _expr= - cst(1.)/J * inner( sigma_tot * trans(F), grad(u_e) ) );
        l += integrate( _range = markedfaces(mesh_, "Upper"), _expr = cst(1)/Js * trans(vec(cst(-hx),cst(0)))*id(u_e));
        l += integrate( _range = markedfaces(mesh_, "Lower"), _expr = cst(1)/Js * trans(vec(cst(hx),cst(0)))*id(u_e));
           
        l += integrate( _range = elements(mesh_), _expr = cst( rho_ )/J * inner( cst(1.0)/(cst(beta)*time_step) * idv( dt_u_old ), id( u_e ) ) );
        l += integrate( _range = elements(mesh_), _expr = cst( rho_ )/J * inner(  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_old ), id( u_e ) ) );

        a.solve(_rhs=l,_solution=u_e);

        // Update time scheme
        dtt_u_old2 = project(_space = Vh_disp_, _range = elements(mesh_init), _expr = idv(dtt_u_old));
        dtt_u_old = project(_space = Vh_disp_, _range = elements(mesh_init), _expr = cst(1.0) / (cst(beta)*std::pow(time_step,2)) * idv(u_e) - cst(1.0) / (cst(beta)*cst(time_step)) * idv( dt_u_old ) -  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_old ));
        dt_u_old = project(_space = Vh_disp_, _range = elements(mesh_init), _expr = idv(dt_u_old) + cst(1.0)*cst(time_step) * ((cst(1.0) - cst(gamma)) * idv(dtt_u_old2) + cst(gamma) * idv(dtt_u_old)) );
        utot  =  project(_space = Vh_disp_, _range = elements(mesh_init), _expr = idv(utot) + idv(u_e) );

        std::cout << "Export" << std::endl;
        this->exportResults_elastic(time_step*iter);

        std::cout << "Update current" << std::endl;
        meshMove( mesh_, u_e );

        //Update
        omega_ = omega_ + time_step * (0.75*hx/momentOfInertia);
        //std::cout << "Angular velocity : " << omega_ << std::endl;
        theta_ = theta_ + omega_*time_step;
        //std::cout << "Rotation angle : " << theta_ << std::endl;
        
        auto rot = vec(
            cos(theta_) * (Px() - massCenter(0,0)) - sin(theta_) * (Py() - massCenter(1,0)) - Px() + massCenter(0,0),
            sin(theta_) * (Px() - massCenter(0,0)) + cos(theta_) * (Py() - massCenter(1,0)) - Py() + massCenter(1,0)
        );

        u_r.on(_range = elements(mesh_), _expr = rot);

        
        // Export results
        this->exportResults_rigid(time_step*iter);

        // Check
        auto l2err = normL2( _range=elements(mesh_), _expr=idv(u_r) - idv(utot) );
        std::cout << "Error : " << l2err << std::endl;

        ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",time_step*iter, theta_, l2err) << std::endl;
        iter++;
    } 
       
}
*/

// Run method 
template <int Dim>
void MagnetoSwimmer<Dim>::run_rigid()
{
    std::cout << "***** Init *****" << std::endl;
    initialize();

    std::cout <<  "***** Export *****" << std::endl;
    this->exportResults_rigid(0);

    std::cout <<  "***** Start time loop *****" << std::endl;
    timeloop_rigid();
}

// Run method 
template <int Dim>
void MagnetoSwimmer<Dim>::run_elastic()
{
    std::cout << "***** Init *****" << std::endl;
    initialize();

    std::cout <<  "***** Export *****" << std::endl;
    //this->exportResults_elastic(0);

    std::cout <<  "***** Start time loop *****" << std::endl;
    timeloop_elastic();
}

// Export results
template <int Dim>
void 
MagnetoSwimmer<Dim>::exportResults_rigid(double t)
{
    e_r->step(t)->addRegions();
    e_r->step(t)->add( "disp", u_r );
    e_r->save();
}

template <int Dim>
void 
MagnetoSwimmer<Dim>::exportResults_elastic(double t)
{
    e_e->step(t)->setMesh(mesh_);
    e_e->step(t)->add( "disp", u_e );
    //e_e->step(t)->add( "u_r", u_r);
    //e_e->step(t)->add( "u_tot", u_tot);
    e_e->save();
}
}

