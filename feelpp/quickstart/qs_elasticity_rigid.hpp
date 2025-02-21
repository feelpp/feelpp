#pragma once

#include "qs_elasticity_contact.hpp"


#if FEELPP_DIM == 2
#define curl_op curlx
#define curlt_op curlxt
#define curlv_op curlxv
#else
#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv
#endif
namespace Feel
{
template <int Dim, int Order>
class ElasticRigid
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;

    // Rigid motion
    using spacev_t_R = Pchv_type<mesh_t, 0>;
    using space_t_R = Pch_type<mesh_t,0>;
    using spacev_ptr_t_R = Pchv_ptrtype<mesh_t, 0>;
    using space_ptr_t_R = Pch_ptrtype<mesh_t, 0>;
    using elementv_t_R = typename spacev_t_R::element_type;
    using element_t_R = typename space_t_R::element_type;

    // elastic motion
    using spacev_t_E = Pchv_type<mesh_t, Order>;
    using space_t_E = Pch_type<mesh_t, Order>;
    using spacev_ptr_t_E = Pchv_ptrtype<mesh_t, Order>;
    using space_ptr_t_E = Pch_ptrtype<mesh_t, Order>;
    using elementv_t_E = typename spacev_t_E::element_type;
    using element_t_E = typename space_t_E::element_type;

    // time schemes
    using ts_ptrtype_E = std::shared_ptr<Newmark<spacev_t_E>>;
    using ts_ptrtype_Translation = std::shared_ptr<Newmark<spacev_t_R>>;
    
    // exporter
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>;

    // Constructors
    ElasticRigid() = default;
    ElasticRigid(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }

    // Inits
    void initializeMesh();
    void initializeParam();
    void initializeFields();
    void initializeTs_Exp();

    // Contact
    void initializeContact();
    Range<mesh_t, MESH_FACES> getContactRegion(elementv_t_E const& u);
    void initG();

    // Run
    void rotationNeumann2DNew();
    void ball_tmp();
    void rigidSwimmer();

    // Export
    void exportResults(double t);
    

private:
    // Constructeur
    nl::json specs_;

    // Mesh
    double H_;
    std::shared_ptr<mesh_t> mesh_init;
    std::shared_ptr<mesh_t> mesh_current;

    // Spaces
    spacev_ptr_t_E Xv_init;
    space_ptr_t_E X_init;
    spacev_ptr_t_E Xv_current;
    space_ptr_t_E X_current;
    
    // Param Solid
    double density_, mass_;
    double E_, nu_, lambda_, mu_;
    std::string externalforce_, neumannUpper_, neumannLower_;
    double extFx_,extFy_;
    double Upper_x, Upper_y, Lower_x, Lower_y;

    // Fields
    elementv_t_E u_e_curr;
    elementv_t_E dt_u_e_tot;
    elementv_t_E dtt_u_e_tot;
    elementv_t_E dt_u_e_old;
    elementv_t_E dtt_u_e_old;
    elementv_t_E dtt_u_e_old2;
    elementv_t_E utotal;

    // Newmark schemes
    double initial_time_, final_time_, time_step_;
    double gamma_, beta_;

    // Exporter
    exporter_ptrtype e_init;
    exporter_ptrtype e_current;

    // Contact param
    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_contact;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    element_t_E contactRegion_;
    element_t_E g_;
    element_t_E contactFaces_;
    int nbrFaces_;
    Range<mesh_t, MESH_FACES> myelts_;

};

// Constructor
template <int Dim, int Order>
ElasticRigid<Dim, Order>::ElasticRigid(nl::json const& specs) : specs_(specs)
{
}


// Initialization Mesh
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeMesh()
{
    // Mesh
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_init = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    mesh_current = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    //mesh_current = createSubmesh( _mesh=mesh_current_total, _range=markedelements( mesh_current_total, "Solid" ) , _worldcomm=Environment::worldCommSeqPtr() );

    // Spaces
    Xv_init  = Pchv<Order>( mesh_init, markedelements( mesh_init, "Solid" ) );
    X_init = Pch<Order>(mesh_init, markedelements( mesh_init, "Solid" ) );
    Xv_current  = Pchv<Order>( mesh_current, markedelements( mesh_current, "Solid" ) );
    X_current = Pch<Order>( mesh_current, markedelements( mesh_current, "Solid" ) );
}


// Initialization body parameters
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeParam()
{
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    density_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    mass_ = integrate( _range = elements(support(Xv_init)), _expr = cst(density_) ).evaluate()(0,0);

    // Young modulus and poisson coefficient
    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());

    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());

    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));

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

                std::size_t offsetE = 0;
                extFx_ = std::stod(&externalforce_[1],&offsetE);
                extFy_ = std::stod(&externalforce_[offsetE+2]);
            }
            else 
                externalforce_ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
        }
    }
    // Neumann boundary conditions
    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Neumann") )
    {
        if ( specs_["/BoundaryConditions/LinearElasticity/Neumann"_json_pointer].contains("Upper") )
        {    
            neumannUpper_ = specs_[nl::json::json_pointer( "/BoundaryConditions/LinearElasticity/Neumann/Upper/expr" )].get<std::string>();
            std::size_t offsetU = 0;
            Upper_x = std::stod(&neumannUpper_[1],&offsetU);
            Upper_y = std::stod(&neumannUpper_[offsetU+2]);
        }
        else 
            neumannUpper_ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
        if ( specs_["/BoundaryConditions/LinearElasticity/Neumann"_json_pointer].contains("Lower") )
        {
            neumannLower_ = specs_[nl::json::json_pointer( "/BoundaryConditions/LinearElasticity/Neumann/Lower/expr" )].get<std::string>();
            std::size_t offsetL = 0;
            Lower_x = std::stod(&neumannLower_[1],&offsetL);
            Lower_y = std::stod(&neumannLower_[offsetL+2]);
        }
        else 
            neumannLower_ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    }
}


// Initialization exporter and newmark time scheme
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeTs_Exp()
{
    // Exporter
    e_init = Feel::exporter(_mesh = mesh_init, _name = "initial", _geo="change" );
    e_current =  Feel::exporter(_mesh = mesh_current, _name = "current", _geo="change" );

    // Initialize Newmark scheme
    initial_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    final_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    time_step_ = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    gamma_ = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    beta_ = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);
}

// Initialization displacement fields
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeFields()
{
    /*
        Elasticity
    */
    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");

    u_e_curr = Xv_current->element();
    u_e_curr.on(_range=elements(support(Xv_current)), _expr= expr<Dim,1>(default_displ));  

    dt_u_e_tot = Xv_init->element();
    dt_u_e_tot.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));   

    dtt_u_e_tot = Xv_init->element();
    dtt_u_e_tot.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));  

    dt_u_e_old = Xv_init->element(); 
    dt_u_e_old.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));  

    dtt_u_e_old = Xv_init->element();
    dtt_u_e_old.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));   

    dtt_u_e_old2 = Xv_init->element();
    dtt_u_e_old2.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));  
    
    /*
        Total displacement
    */
    utotal = Xv_init->element(); 
    utotal.on(_range=elements(support(Xv_init)), _expr= expr<Dim,1>(default_displ));    
}


// Export results
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::exportResults(double t)
{
    std::cout << "Export" << std::endl;
    e_init->step(t)->setMesh(mesh_init);
    e_init->step(t)->add( "utotal", utotal );
    e_init->step(t)->add("distance",g_);
    e_init->step(t)->add("contactRegion",contactRegion_);
    e_init->save(); 
    

    e_current->step(t)->setMesh(mesh_current);
    e_current->step(t)->add( "u_e_curr", u_e_curr );
    e_current->step(t)->add("contactFaces",contactFaces_);
    e_current->save();
}

template <int Dim, int Order>
void ElasticRigid<Dim, Order>::ball_tmp()
{
    if constexpr(Dim == 2)
    {
        // Init
        std::cout << "Init mesh" << std::endl;
        this->initializeMesh();
        std::cout << "Init param solid" << std::endl;
        this->initializeParam();
        std::cout << "Init param contact" << std::endl;
        this->initializeContact();
        std::cout << "Init exporter and time data" << std::endl;
        this->initializeTs_Exp();
        std::cout << "Init fields" << std::endl;
        this->initializeFields();
        
        
        // Get distance
        std::cout << "Distance" << std::endl;
        this->initG();

        // Export
        this->exportResults(0);

        // Translation
        std::cout << "Init translation" << std::endl;
        Eigen::Vector2d u_trans(0.,0.);
        Eigen::Vector2d u_trans_tot(0.,0.);
        Eigen::Vector2d dt_u_trans_old(0.,0.);
        Eigen::Vector2d dtt_u_trans_old(0.,0.);
        Eigen::Vector2d dtt_u_trans_old2(0.,0.);
        Eigen::Vector2d Force(extFx_,extFy_);
        Eigen::Vector2d ForceContact(0.,0.);

        // Rotation
        std::cout << "Init rotation" << std::endl;
        double velocity = 0;
        double angle = 0;


        // Starting time loop
        int iter = 0;
        std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
        auto Id = eye<Dim,Dim>();

        std::ofstream ofs("outputs_ball.csv");
        ofs << fmt::format("utot_x,utot_y,u_trans_tot_x,u_trans_tot_y, fc_x, fc_y") << std::endl;

        while (time_step_ * iter < final_time_)
        {
            iter++;
            if (Environment::isMasterRank())
                std::cout << "Time step : " << iter*time_step_ << std::endl;
    
            /*
                Outputs 
            */ 
            auto meanTot = mean(_range=elements(support(Xv_init)), _expr=idv(utotal));
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}",meanTot(0,0),meanTot(1,0),u_trans_tot[0],u_trans_tot[1],ForceContact[0],ForceContact[1]) << std::endl;
            
            /*
                Solve translation 
            */
            std::cout << "Solve translation" << std::endl;
            std::cout << "Print contactForce : " << ForceContact[0] << " " << ForceContact[1] << std::endl;
            u_trans = beta_ * std::pow(time_step_,2) * ((1./(beta_ * time_step_) * dt_u_trans_old + (1. - 2*beta_)/(2*beta_) * dtt_u_trans_old) + Force + ForceContact) ; 

            // Update Newmark scheme
            dtt_u_trans_old2 = dtt_u_trans_old;
            dtt_u_trans_old = 1. / (beta_ *std::pow(time_step_,2)) * u_trans - 1. / (beta_ * time_step_) *  dt_u_trans_old -  (1./(2.*beta_) - 1.) *  dtt_u_trans_old;
            dt_u_trans_old = dt_u_trans_old + time_step_ * ((1. - gamma_) * dtt_u_trans_old2 + gamma_ * dtt_u_trans_old) ;
            u_trans_tot  += u_trans;
            
            utotal = project(_space = Xv_init, _range = elements(support(Xv_init)), _expr = idv(utotal) + vec(cst(u_trans[0]),cst(u_trans[1])) );
            
            auto translational_move =  project(_space =  Pchv<Order>( mesh_current ),  _range = elements(mesh_current), _expr = vec(cst(0.),cst(0.) )); 
            translational_move = project(_space = Xv_current,  _range = elements(support(Xv_current)), _expr = vec(cst(u_trans[0]),cst(u_trans[1])) );
            meshMove( mesh_current, translational_move); 

            /*
                Solve rotation
            */

            /*
                Solve contact
            */
            
            std::cout << "Solve contact" << std::endl;
            
            Xv_current = Pchv<Order>( mesh_current, markedelements( mesh_current, "Solid" ) );
            X_current = Pch<Order>( mesh_current, markedelements( mesh_current, "Solid" ) );

            u_e_curr = Xv_current->element();
            u_e_curr.on( _range=elements(support(Xv_current)), _expr=expr<Dim,1>(default_displ));
                
            // Compute contact faces
            myelts_ = getContactRegion(utotal);
            std::cout << "Nombre of faces in contact : " << nbrFaces_ << std::endl;
             
            /*
                Solve elasticity
            */
            
            std::cout << "Solve elasticity" << std::endl;

            auto eps_curr = sym(gradt(u_e_curr));
            auto sigma_curr = lambda_*trace(eps_curr)*Id + 2*mu_*eps_curr;

            auto eps_tot = sym(gradv(utotal)); 
            auto F = Id + gradv(utotal); 
            auto J = det(F); 
            auto Js = J * sqrt(trans(N())*( trans(F) * F) * N()); 
            auto sigma_tot = lambda_*trace(eps_tot)*Id + 2*mu_*eps_tot;

            auto eps_mix = 0.5 * (gradt(u_e_curr)*gradv(utotal) + trans(gradv(utotal))*trans(gradt(u_e_curr)));
            auto sigma_mix = lambda_ * trace(eps_mix) * Id + 2 * mu_ * eps_mix;
        
            auto a = form2( _trial=Xv_current, _test=Xv_current);
            auto l = form1( _test=Xv_current );

            a.zero();
            l.zero();

            a += integrate(_range = elements(support(Xv_current)), _expr = cst(1.)/J * inner( sigma_curr * trans(F), grad(u_e_curr) ));
            a += integrate(_range = elements(support(Xv_current)), _expr = cst(1.)/J * inner( sigma_mix * trans(F), grad(u_e_curr) ));
            a += integrate(_range = elements(support(Xv_current)), _expr= cst(density_)/J * inner( cst(1.0)/(cst(beta_)*std::pow(time_step_,2)) * idt( u_e_curr ),id( u_e_curr ) ) );
            
                
            l += integrate( _range = elements(support(Xv_current)), _expr= - cst( density_ )/J * inner(vec(cst(dtt_u_trans_old[0]),cst(dtt_u_trans_old[1])), id( u_e_curr ) )); 
            l += integrate( _range = elements(support(Xv_current)), _expr= - cst(1.)/J * inner( sigma_tot * trans(F), grad(u_e_curr) ) );
            
            l += integrate( _range = elements(support(Xv_current)), _expr = cst( density_ )/J * inner( cst(1.0)/(cst(beta_)*time_step_) * idv( dt_u_e_old ), id( u_e_curr ) ) );
            l += integrate( _range = elements(support(Xv_current)), _expr = cst( density_ )/J * inner(  (cst(1.0)/(cst(2.0)*cst(beta_)) - cst(1.0)) * idv( dtt_u_e_old ), id( u_e_curr ) ) );
            l += integrate( _range = elements(support(Xv_current)), _expr= cst( density_ ) / J * trans(expr<Dim,1>(externalforce_))*id(u_e_curr) );

            // Contact 
            if (nbrFaces_ > 0)
            {
                //a += integrate(_range = myelts_, _expr= cst(1.)/(cst(epsilon_)*Js) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_curr),trans(expr<Dim,1>(direction_))*id(u_e_curr)));
                //l += integrate (_range = myelts_,_expr= - cst(1.)/(cst(epsilon_)*Js) * inner(trans(expr<Dim,1>(direction_))*idv(utotal) - idv(g_),trans(expr<Dim,1>(direction_))*id(u_e_curr)) );
                a += integrate(_range = myelts_, _expr = cst(1.)/(cst(epsilon_)*Js) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_curr),trans(expr<Dim,1>(direction_))*id(u_e_curr)));
                l += integrate (_range = myelts_,_expr = - cst(1.)/(cst(epsilon_)*Js) * inner(trans(expr<Dim,1>(direction_))*idv(utotal) - idv(g_),trans(expr<Dim,1>(direction_))*id(u_e_curr)) );
            }

            a.solve(_rhs=l,_solution=u_e_curr);
            

            // Update time scheme
            dtt_u_e_old2 = project(_space = Xv_init, _range = elements(support(Xv_init)), _expr = idv(dtt_u_e_old));
            dtt_u_e_old = project(_space = Xv_init, _range = elements(support(Xv_init)), _expr = cst(1.0) / (cst(beta_)*std::pow(time_step_,2)) * idv(u_e_curr) - cst(1.0) / (cst(beta_)*cst(time_step_)) * idv( dt_u_e_old ) -  (cst(1.0)/(cst(2.0)*cst(beta_)) - cst(1.0)) * idv( dtt_u_e_old ));
            dt_u_e_old = project(_space = Xv_init, _range = elements(support(Xv_init)), _expr = idv(dt_u_e_old) + cst(time_step_) * ((cst(1.0) - cst(gamma_)) * idv(dtt_u_e_old2) + cst(gamma_) * idv(dtt_u_e_old)) );
            utotal  =  project(_space = Xv_init, _range = elements(support(Xv_init)), _expr = idv(utotal) + idv(u_e_curr) );

            // Compute new contact force
            if (nbrFaces_ == 0)
            {
                ForceContact[0] = 0.;
                ForceContact[1] = 0.;
            } 
            else 
            {
                auto f_c = integrate(_range=myelts_,_expr= cst(1)/Js * trans(expr<Dim,1>(direction_))*(lambda_*trace(sym(gradv(u_e_curr)))*Id + 2*mu_*sym(gradv(u_e_curr)))*expr<Dim,1>(direction_)).evaluate();
                std::cout << f_c(0,0) << std::endl;
    
                auto f_c_3 = integrate( _range = myelts_, _expr = cst(1)/Js * (lambda_*trace(sym(gradv(u_e_curr)))*Id + 2*mu_*sym(gradv(u_e_curr)))*N() ).evaluate();
                std::cout << f_c_3(0,0) << ", " <<  f_c_3(1,0) << std::endl;

                ForceContact[0] = f_c_3(0,0); 
                ForceContact[1] = f_c_3(1,0);

            }

            // Checks
            auto meanDisp = mean(_range=elements(support(Xv_current)), _expr=idv(u_e_curr));
            std::cout << "Mean displacement x : " << meanDisp(0,0) << std::endl;
            std::cout << "Mean displacement y : " << meanDisp(1,0) << std::endl;
            auto meanCurl = mean(_range=elements(support(Xv_current)), _expr = curlv_op(u_e_curr));
            std::cout << "Mean curl x : " << meanCurl(0,0) << std::endl;
            std::cout << "Mean curl y : " << meanCurl(1,0) << std::endl;
            std::cout << "Mean curl z : " << meanCurl(2,0) << std::endl;

            // Export 
            this->exportResults(iter*time_step_);

            // Move
            auto elastic_move =  project(_space =  Pchv<Order>( mesh_current ),  _range = elements(mesh_current), _expr = vec(cst(0.),cst(0.) )); 
            elastic_move = project(_space =  Xv_current,  _range = elements(support(Xv_current)), _expr =idv(u_e_curr)); 
            meshMove( mesh_current, elastic_move );
 
        }
        ofs.close();
    }
}

template <int Dim, int Order>
void ElasticRigid<Dim, Order>::rotationNeumann2DNew()
{
    if constexpr(Dim == 2)
    {
        // Get parameters
        std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
        double density = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

        // Young modulus and poisson coefficient
        std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
        double E = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());

        std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
        double nu = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());

        double lambda = E*nu/( (1+nu)*(1-2*nu) );
        double mu = E/(2*(1+nu));

        std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
        auto Id = eye<Dim,Dim>();

        // Neumann boundary conditions
        std::string neumannLower = specs_[nl::json::json_pointer( "/BoundaryConditions/LinearElasticity/Neumann/Lower/expr" )].get<std::string>();
        std::string neumannUpper = specs_[nl::json::json_pointer( "/BoundaryConditions/LinearElasticity/Neumann/Upper/expr" )].get<std::string>();

        std::size_t offsetU = 0;
        double Upper_x = std::stod(&neumannUpper[1],&offsetU);
        double Upper_y = std::stod(&neumannUpper[offsetU+2]);

        std::cout << "Upper_x : " << Upper_x << " Upper_y : " << Upper_y << std::endl;
        
        std::size_t offsetL = 0;
        double Lower_x = std::stod(&neumannLower[1],&offsetL);
        double Lower_y = std::stod(&neumannLower[offsetL+2]);

        std::cout << "Lower_x : " << Lower_x << " Lower_y : " << Lower_y << std::endl;
        
        // External force
        std::string externalforce = specs_[nl::json::json_pointer( "/Models/LinearElasticity/loading/F1/parameters/expr" )].get<std::string>();

        // Time scheme parameters
        double initial_time = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
        double final_time = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
        double time_step = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    
        double gamma = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
        double beta = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);


        // Get init mesh
        double H = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
        auto mesh = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H);
        auto mesh_init = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H);

        // Exporter
        auto e = Feel::exporter(_mesh = mesh, _name = "ucurr", _geo="change" );
        auto e_init =  Feel::exporter(_mesh = mesh, _name = "utot", _geo="change" );

        // Initialisation
        auto Vhv = Pchv<Order>( mesh );
        auto Vhv_init = Pchv<Order>( mesh_init );

        // Elasticity
        auto u_curr = Vhv->element();
        u_curr.on(_range=elements( mesh ), _expr= expr<Dim,1>(default_displ));  

        auto u_tot_curr = Vhv->element(); 
        u_tot_curr.on(_range=elements( mesh ), _expr= expr<Dim,1>(default_displ));  

        auto utot = Vhv_init->element(); 
        utot.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto utot_e = Vhv_init->element(); 
        utot_e.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto dt_u_tot = Vhv_init->element();
        dt_u_tot.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));   

        auto dtt_u_tot = Vhv_init->element();
        dtt_u_tot.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto dt_u_old = Vhv_init->element(); 
        dt_u_old.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto dtt_u_old = Vhv_init->element();
        dtt_u_old.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));   

        auto dtt_u_old2 = Vhv_init->element();
        dtt_u_old2.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        // Rotation
        auto u_rot = Vhv->element();
        u_rot.on(_range=elements( mesh ), _expr= expr<Dim,1>(default_displ));  

        auto u_rot_tot = Vhv_init->element(); 
        u_rot_tot.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto dt_u_rot_old = Vhv_init->element(); 
        dt_u_rot_old.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        auto dtt_u_rot_old = Vhv_init->element();
        dtt_u_rot_old.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));   

        auto dtt_u_rot_old2 = Vhv_init->element();
        dtt_u_rot_old2.on(_range=elements( mesh_init ), _expr= expr<Dim,1>(default_displ));  

        std::cout << "Export" << std::endl;
        e_init->step(0)->setMesh(mesh_init);
        e_init->step(0)->add( "utot", utot );
        e_init->step(0)->add( "urot", u_rot_tot);
        e_init->step(0)->add( "utot_e", utot_e);
        e_init->save(); 
            
        e->step(0)->setMesh(mesh);
        e->step(0)->add( "u_curr", u_curr );
        e->step(0)->add( "u_rot", u_rot );
        e->save();

        double velocity = 0;

        // Starting time loop
        int iter = 0;

        std::ofstream ofs("outputs_new.csv");
        ofs << fmt::format("x_curr,y_curr,utot_x,utot_y,utot_e_x,utot_e_y") << std::endl;


        while (time_step * iter < final_time)
        {
            iter++;
            if (Environment::isMasterRank())
                std::cout << "Time step : " << iter*time_step << std::endl;
    
            // Exports
            auto mass = mean( _range = elements(mesh), _expr = P());
            auto meanE = mean(_range=elements(mesh_init), _expr=idv(utot_e));
            auto meanTot = mean(_range=elements(mesh_init), _expr=idv(utot));
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}",mass(0,0),mass(1,0),meanTot(0,0),meanTot(1,0),meanE(0,0),meanE(1,0)) << std::endl;
           
            
            /*
                Solve rotational 
            */
            std::cout << "Solve rotational" << std::endl;

            auto massCenter = mean( _range = elements(mesh), _expr = P());
            auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
            std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

            // Matrix of inertia
            auto momentOfInertia = integrate(_range=elements(mesh),_expr=cst(density)*( (Px()-massCenter(0,0))*(Px()-massCenter(0,0)) + (Py()-massCenter(1,0))*(Py()-massCenter(1,0)) ) ).evaluate()(0,0);
            std::cout << "momentOfInertia : " << momentOfInertia << std::endl;

            // Rotation
            double h = integrate(_range= markedfaces(mesh, "Upper"), _expr= -(Py()-massCenter(1,0))*cst(Upper_x)).evaluate()(0,0);
            h += integrate(_range= markedfaces(mesh, "Lower"), _expr= -(Py()-massCenter(1,0))*cst(Lower_x)).evaluate()(0,0);
            std::cout << "h : " << h << std::endl;
            velocity = h/momentOfInertia * time_step + velocity;
            double theta = velocity * time_step;

            std::cout << "Velocity : " << velocity << std::endl;
            std::cout << "theta : " << theta << std::endl;

            auto rot = vec(cos(theta)*(Px() - massCenter(0,0)) - sin(theta)*(Py() - massCenter(1,0)) - (Px() - massCenter(0,0)), sin(theta)*(Px() - massCenter(0,0)) + cos(theta)* (Py() - massCenter(1,0)) - (Py() - massCenter(1,0)));
            //u_rot.on(_range=elements(mesh),_expr = rot);
            u_rot = project(_space = Vhv, _range = elements(mesh), _expr = rot);

            // Update time scheme
            dtt_u_rot_old2 = project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(dtt_u_rot_old));
            dtt_u_rot_old = project(_space = Vhv_init, _range = elements(mesh_init), _expr = cst(1.0) / (cst(beta)*std::pow(time_step,2)) * idv(u_rot) - cst(1.0) / (cst(beta)*cst(time_step)) * idv( dt_u_rot_old ) -  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_rot_old ));
            dt_u_rot_old = project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(dt_u_rot_old) + cst(1.0)*cst(time_step) * ((cst(1.0) - cst(gamma)) * idv(dtt_u_rot_old2) + cst(gamma) * idv(dtt_u_rot_old)) );
            u_rot_tot  =  project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(u_rot_tot) + idv(u_rot)  );
            utot  =  project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(utot) + idv(u_rot) );

            auto test = mean( _range = elements(mesh_init), _expr = idv(dtt_u_rot_old));
            std::cout << "Test : " << test(0,0) << std::endl;

            meshMove( mesh, u_rot);

            /*
                Solve elasticity
            */
            std::cout << "Solve elasticity" << std::endl;

            auto eps_curr = sym(gradt(u_curr));
            auto sigma_curr = lambda*trace(eps_curr)*Id + 2*mu*eps_curr;

            auto eps_tot = sym(gradv(utot)); 
            auto F = Id + gradv(utot); 
            auto J = det(F); 
            auto Js = J * sqrt(trans(N())*( trans(F) * F) * N()); 
            auto sigma_tot = lambda*trace(eps_tot)*Id + 2*mu*eps_tot;

            auto eps_mix = 0.5 * (gradt(u_curr)*gradv(utot) + trans(gradv(utot))*trans(gradt(u_curr)));
            auto sigma_mix = lambda * trace(eps_mix) * Id + 2 * mu * eps_mix;
        
            auto a = form2( _trial=Vhv, _test=Vhv);
            auto l = form1( _test=Vhv );

            a.zero();
            l.zero();

            std::cout << "Assemble a" << std::endl;

            a += integrate(_range = elements(mesh), _expr = cst(1.)/J * inner( sigma_curr * trans(F), grad(u_curr) ));
            a += integrate(_range = elements(mesh), _expr = cst(1.)/J * inner( sigma_mix * trans(F), grad(u_curr) ));
            a += integrate(_range = elements(mesh), _expr= cst(density)/J * inner( cst(1.0)/(cst(beta)*std::pow(time_step,2)) * idt( u_curr ),id( u_curr ) ) );
            
            std::cout << "assemle l" << std::endl;
            l += integrate( _range = elements(mesh), _expr= - cst( density )/J * inner(idv(dtt_u_rot_old), id( u_curr ) )); 
            l += integrate( _range = elements(mesh), _expr= - cst(1.)/J * inner( sigma_tot * trans(F), grad(u_curr) ) );
            l += integrate( _range = markedfaces(mesh, "Upper"), _expr = cst(1)/Js * trans(expr<Dim,1>(neumannUpper))*id(u_curr));
            l += integrate( _range = markedfaces(mesh, "Lower"), _expr = cst(1)/Js * trans(expr<Dim,1>(neumannLower))*id(u_curr));
            
            l += integrate( _range = elements(mesh), _expr = cst( density )/J * inner( cst(1.0)/(cst(beta)*time_step) * idv( dt_u_old ), id( u_curr ) ) );
            l += integrate( _range = elements(mesh), _expr = cst( density )/J * inner(  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_old ), id( u_curr ) ) );
            l += integrate( _range = elements(mesh), _expr= cst( density ) / J * trans(expr<Dim,1>(externalforce))*id(u_curr) );

            a.solve(_rhs=l,_solution=u_curr, _rebuild = true);

            // Update time scheme
            dtt_u_old2 = project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(dtt_u_old));
            dtt_u_old = project(_space = Vhv_init, _range = elements(mesh_init), _expr = cst(1.0) / (cst(beta)*std::pow(time_step,2)) * idv(u_curr) - cst(1.0) / (cst(beta)*cst(time_step)) * idv( dt_u_old ) -  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_old ));
            dt_u_old = project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(dt_u_old) + cst(1.0)*cst(time_step) * ((cst(1.0) - cst(gamma)) * idv(dtt_u_old2) + cst(gamma) * idv(dtt_u_old)) );
            utot  =  project(_space = Vhv_init, _range = elements(mesh_init), _expr = idv(utot) + idv(u_curr) );
            utot_e = project(_space = Vhv_init,  _range = elements(mesh_init), _expr = idv(utot_e) + idv(u_curr) );

            // Checks
            auto meanDisp = mean(_range=elements(mesh), _expr=idv(u_curr));
            std::cout << "Mean displacement x : " << meanDisp(0,0) << std::endl;
            std::cout << "Mean displacement y : " << meanDisp(1,0) << std::endl;
        
            auto meanCurl = mean(_range=elements(mesh), _expr = curlv_op(u_curr));
            std::cout << "Mean curl x : " << meanCurl(0,0) << std::endl;
            std::cout << "Mean curl y : " << meanCurl(1,0) << std::endl;
            std::cout << "Mean curl z : " << meanCurl(2,0) << std::endl;

            std::cout << "Export" << std::endl;
            e_init->step(time_step*iter)->setMesh(mesh_init);
            e_init->step(time_step*iter)->add( "utot", utot );
            e_init->step(time_step*iter)->add( "urot", u_rot_tot);
            e_init->step(time_step*iter)->add( "utot_e", utot_e);
            e_init->save(); 
            
            e->step(time_step*iter)->setMesh(mesh);
            e->step(time_step*iter)->add( "u_curr", u_curr );
            e->step(time_step*iter)->add( "u_rot", u_rot );
            e->save();

            std::cout << "Update current" << std::endl;
            meshMove( mesh, u_curr );
            Vhv = Pchv<Order> (mesh);
            
            u_rot = Vhv->element();
            u_rot.on(_range=elements(mesh), _expr=expr<Dim,1>(default_displ)); 

            u_curr = Vhv->element();
            u_curr.on(_range=elements(mesh), _expr=expr<Dim,1>(default_displ)); 
            
        }
        ofs.close();

    }
    
}


// Initialization contact parameters 
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initializeContact()
{
    // Initialize contact field
    contactRegion_ =  project(_space=X_init, _range=elements(support(X_init)), _expr = cst(0.));
    contactFaces_ = project(_space=X_current, _range=elements(support(X_current)), _expr = cst(0.));
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
    gamma_contact = gamma0_/H_;

    std::string mattolContactRegion = fmt::format( "/Collision/LinearElasticity/tolContactRegion" );
    tolContactRegion_ = specs_[nl::json::json_pointer( mattolContactRegion )].get<double>();

    std::string mattolDistance = fmt::format("/Collision/LinearElasticity/tolDistance");
    tolDistance_ = specs_[nl::json::json_pointer( mattolDistance )].get<double>();

}

template <int Dim, int Order>
Range<typename ElasticRigid<Dim, Order>::mesh_t, MESH_FACES>
ElasticRigid<Dim, Order>::getContactRegion(elementv_t_E const& u)
{
    Range<mesh_t,MESH_FACES> myelts(mesh_init);
    
    contactRegion_ = project(_space = X_init,  _range = elements(support(X_init)), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));
    
    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=X_init, _trial=X_init).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(support(X_init)) )
    {
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : X_init->dof()->faceLocalDof( face.id() ) )
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
                        contactFaces_[thedof] = 1.;
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
        }
    }
    myelts.shrink_to_fit();

    return myelts;
}

/*
template <int Dim, int Order>
Range<typename ElasticRigid<Dim, Order>::mesh_t, MESH_FACES>
ElasticRigid<Dim, Order>::getContactRegion(elementv_t_E const& u)
{
    Range<mesh_t,MESH_FACES> myelts(mesh_init);
    contactRegion_.on( _range=elements(support(X_init)), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));
    
    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=X_init, _trial=X_init).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(support(X_init)) )
    {
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : X_init->dof()->faceLocalDof( face.id() ) )
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
                        contactFaces_[thedof] = 1.;
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
        }
    }
    myelts.shrink_to_fit();

    return myelts;
}
*/
/*
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initG()
{
    // Init the distance fields
    g_ = X_current->element();
    g_.on(_range=elements(support(X_current)), _expr=cst(1000.));

    // Raytracing to compute distance
    using bvh_ray_type = BVHRay<Dim>;
    Eigen::VectorXd origin(Dim);
    Eigen::VectorXd dir(Dim);

    if constexpr(Dim == 2)
        dir << ddirection_[0], ddirection_[1];
    else if constexpr(Dim == 3)
        dir << ddirection_[0], ddirection_[1], ddirection_[2];

    std::string kind = (Dim==2)?"in-house":"third-party";

    auto bvh = boundingVolumeHierarchy(_range=markedfaces(mesh_current, "Obs1"), _kind=kind);

    std::vector<std::size_t> faceIDs;

    BVHRaysDistributed<Dim> allrays;

    for ( auto const& theface : markedfaces( mesh_current, "Wall" ) )
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
                        for (auto const& ldof  : X_current->dof()->faceLocalDof( face.id() ))
                        {
                            g_[ldof.index()] = rir.distance() - tolDistance_ ;
                        }
                            
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
}
*/

template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initG()
{
    // Init the distance fields
    g_ = X_init->element();
    g_.on(_range=elements(support(X_init)), _expr=cst(10.));

    // Raytracing to compute distance
    using bvh_ray_type = BVHRay<Dim>;
    Eigen::VectorXd origin(Dim);
    Eigen::VectorXd dir(Dim);

    if constexpr(Dim == 2)
        dir << ddirection_[0], ddirection_[1];
    else if constexpr(Dim == 3)
        dir << ddirection_[0], ddirection_[1], ddirection_[2];

    std::string kind = (Dim==2)?"in-house":"third-party";

    auto bvh = boundingVolumeHierarchy(_range=markedfaces(mesh_init, "Obs1"), _kind=kind);

    std::vector<std::size_t> faceIDs;

    BVHRaysDistributed<Dim> allrays;

    for ( auto const& theface : markedfaces( mesh_init, "Wall" ) )
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
                        for (auto const& ldof  : X_init->dof()->faceLocalDof( face.id() ))
                        {
                            g_[ldof.index()] = rir.distance() - tolDistance_ ;
                        }
                            
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
}

}
/*
auto Res = backend()->newVector(Vhv);
            auto Jac = backend()->newMatrix( _test=Vhv, _trial=Vhv );

            auto eps_tot = sym(gradv(utot)); 
            auto F = Id + gradv(utot); 
            auto J = det(F); 
            auto Js = J * sqrt(trans(N())*( trans(F) * F) * N()); 
            auto sigma_tot = lambda*trace(eps_tot)*Id + 2*mu*eps_tot;

            auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J_)
            {
                auto u = Vhv->element();
                u = *X;

                auto deps = sym(gradt(u));
                auto dsigma = lambda*trace(deps)*Id + 2*mu*deps;

                auto deps_mix = 0.5 * (gradt(u)*gradv(utot) + trans(gradv(utot))*trans(gradt(u)));
                auto dsigma_mix = lambda * trace(deps_mix) * Id + 2 * mu * deps_mix;
        
                auto a = form2( _test=Vhv, _trial=Vhv, _matrix=J_ );
                
                a = integrate( _range = elements(mesh), _expr= cst(1.)/J * inner( dsigma * trans(F), grad(u) ) );
                a += integrate( _range = elements(mesh), _expr= cst(1.)/J * inner( dsigma_mix * trans(F), grad(u) ) );

                a += integrate( _range=elements(mesh), _expr=  cst( density )/J * inner( cst(1.0)/(cst(beta)*std::pow(time_step,2)) * idt( u ),id( u ) ) );
            
            };

            auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
            {
                auto u = Vhv->element();
                u = *X;
            
                auto eps_v = sym(gradv(u));
                auto sigma_v = lambda*trace(eps_v)*Id + 2*mu*eps_v;

                auto eps_mix = 0.5 * (gradv(u)*gradv(utot) + trans(gradv(utot))*trans(gradv(u)));
                auto sigma_mix = lambda * trace(eps_mix) * Id + 2 * mu * eps_mix;
        
                auto r = form1( _test=Vhv, _vector=R );

                r += integrate( _range = elements(mesh), _expr = -cst( density )/J * trans( expr<Dim, 1>( externalforce ) )* id( u )  ); // external forces
                r += integrate( _range = elements(mesh), _expr =  cst( density )/J * inner( cst(1.0)/(cst(beta)*std::pow(time_step,2)) * idv( u ),id( u ) ) ); // acceleration elastic
                r += integrate( _range = elements(mesh), _expr = -cst( density )/J * inner( cst(1.0)/(cst(beta)*time_step) * idv( dt_u_old ), id( u ) ) ); // acceleration elastic
                r += integrate( _range = elements(mesh), _expr = -cst( density )/J * inner(  (cst(1.0)/(cst(2.0)*cst(beta)) - cst(1.0)) * idv( dtt_u_old ), id( u ) ) ); // acceleration elastic
                r += integrate( _range = elements(mesh), _expr =  cst( density )/J * inner(vec(cst(dtt_u_trans_old[0]),cst(dtt_u_trans_old[1])), id( u ) )); // acceleration translation
                r += integrate( _range = elements(mesh), _expr = cst(1.)/J * inner( val(sigma_v * trans(F)), grad(u) )); // elastic term
                r += integrate( _range = elements(mesh), _expr = cst(1.)/J * inner( val(sigma_mix * trans(F)), grad(u) )); // elastic term
                r += integrate( _range = elements(mesh), _expr = cst(1.)/J * inner( sigma_tot * trans(F), grad(u) ) ); // elastic term
                
                R->close();
            };

            backend()->nlSolver()->residual = Residual;
            backend()->nlSolver()->jacobian = Jacobian;
            backend()->nlSolve( _solution=u_curr,_jacobian=Jac,_residual=Res );

*/