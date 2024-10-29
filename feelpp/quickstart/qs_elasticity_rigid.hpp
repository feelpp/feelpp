#pragma once

#include "qs_elasticity_contact.hpp"
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
    void rotation();
    void rotationNeumann();
    void translationContact();

    // Export
    void exportResults(double t);
    

private:
    // Constructeur
    nl::json specs_;

    // Mesh
    double H_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t_E XvE_;
    space_ptr_t_E XE_;
    spacev_ptr_t_R XvR_;
    space_ptr_t_R XR_;

    // Param
    double density_, mass_;
    double E_, nu_, lambda_, mu_;
    std::string externalforce_;
    double extFx_,extFy_;

    // Contact param
    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    std::vector<double> pressurePoint_;
    int fixedPoint_;
    double fixedPointtol_;
    element_t_E contactRegion_;
    element_t_E g_;
    int nbrFaces_;
    Range<mesh_t, MESH_FACES> myelts_;
    

    // Fields
    elementv_t_E u_total_; 
    elementv_t_E u_E_;
    elementv_t_E u_R_;
    elementv_t_R u_translation_; 
    elementv_t_E u_rotation_;

    // Newmark schemes
    double initial_time_, final_time_, time_step_;
    ts_ptrtype_E ts_E_;
    ts_ptrtype_Translation ts_Translation_;

    // Exporter
    exporter_ptrtype e_;

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
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    XvE_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );
    XvR_ = Pchv<0>(mesh_, markedelements(mesh_, "Caoutchouc"));
    XE_ = Pch<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );
    XR_ = Pch<0>(mesh_, markedelements(mesh_, "Caoutchouc"));

}

// Initialization body parameters
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeParam()
{
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    density_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    mass_ = integrate( _range = elements(support(XvE_)), _expr = cst(density_) ).evaluate()(0,0);

    // Young modulus and poisson coefficient
    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());

    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());

    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));

    // External force
    std::size_t offset = 0;
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
                extFx_ = std::stod(&externalforce_[1],&offset);
                extFy_ = std::stod(&externalforce_[offset+2]);
            }
        }
    }
}

// Initialization displacement fields
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeFields()
{
    // Displacement
    u_total_ = XvE_->element();
    u_E_ = XvE_->element();
    u_R_ = XvE_->element();
    u_translation_ = XvR_->element();
    u_rotation_ = XvE_->element();

    // Init to zero
    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    
    u_total_.on(_range=elements(support(XvE_)), _expr= expr<Dim,1>(default_displ));    
    u_E_.on(_range=elements(support(XvE_)), _expr=expr<Dim,1>(default_displ));    
    u_R_.on(_range=elements(support(XvE_)), _expr=expr<Dim,1>(default_displ));    
    u_translation_.on(_range=elements(support(XvR_)), _expr=expr<Dim,1>(default_displ));    
    u_rotation_.on(_range=elements(support(XvE_)), _expr=expr<Dim,1>(default_displ));    
}

// Initialization exporter and newmark time scheme
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeTs_Exp()
{
    // Exporter
    e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    this->exportResults(0);

    // Initialize Newmark scheme
    initial_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    final_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    time_step_ = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    
    double gamma = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);

    ts_E_ = newmark(_space = XvE_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_, _gamma=gamma, _beta=beta );
    ts_E_->start();
    ts_E_->initialize( u_E_ );
    ts_E_->updateFromDisp(u_E_);

    ts_Translation_ = newmark(_space = XvR_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_, _gamma=gamma, _beta=beta );
    ts_Translation_->start();
    ts_Translation_->initialize( u_translation_ );
    ts_Translation_->updateFromDisp( u_translation_);
}

// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::rotationNeumann()
{
    
    auto Id = eye<Dim,Dim>();
    auto iter = -1;
    // Inits
    this->initializeMesh();
    this->initializeParam();
    this->initializeFields();
    this->initializeTs_Exp();

    
    /*
        Assemblage : elasticity
    */
    auto a_e = form2( _test = XvE_, _trial = XvE_ );
    auto at_e = form2( _test = XvE_, _trial = XvE_ );
    auto l_e = form1( _test = XvE_ );
    auto lt_e = form1( _test = XvE_ );

    a_e.zero();
    at_e.zero();
    l_e.zero();
    lt_e.zero();

    auto deft = sym(gradt(u_E_));
    auto def = sym(grad(u_E_));
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a_e += integrate( _range = elements(support(XvE_)), _expr = cst(density_) * inner( ts_E_->polyDerivCoefficient()*idt(u_E_),id( u_E_ ) ) + inner(sigmat,def));
    l_e = integrate( _range = elements(support(XvE_)), _expr = cst(density_) * trans( expr<Dim, 1>( externalforce_ ) ) * id( u_E_ ) );
    l_e += integrate( _range = markedfaces(mesh_, "Neumann2"), _expr = trans(vec(cst(1.),cst(0.),cst(0.)))*id(u_E_));
    l_e += integrate( _range = markedfaces(mesh_, "Neumann1"), _expr = trans(vec(cst(-1.),cst(0.),cst(0.)))*id(u_E_));

    // Starting time loop
    for ( ts_E_->start(); ts_E_->isFinished() == false; ts_E_->next( u_E_ ))
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_E_->time(),ts_E_->timeFinal()) << std::endl;

        /*
            Solve elasticity
        */
        lt_e = l_e;
        at_e = a_e;

        lt_e += integrate( _range=elements(support(XvE_)), _expr= cst(density_)*inner( idv(ts_E_->polyDeriv()),id( u_E_ ) ));
                
        // Delete rigid motion 
        
        at_e+=on(_range=markedpoints(mesh_,"a0"), _rhs=lt_e, _element=u_E_, _expr=0.*one()); // translation
        
        // Rotation 
        
        if constexpr(Dim == 3)
        {
            std::cout << "Delete rotation" << std::endl;
            at_e+=on(_range=markedpoints(mesh_,"a1"), _rhs=lt_e, _element=u_E_, _expr=vec(cst(0.),cst(1.),cst(1.))*(idv(u_E_)));
            at_e+=on(_range=markedpoints(mesh_,"a2"), _rhs=lt_e, _element=u_E_, _expr=vec(cst(1.),cst(0.),cst(1.))*(idv(u_E_)));
            at_e+=on(_range=markedpoints(mesh_,"a3"), _rhs=lt_e, _element=u_E_, _expr=vec(cst(1.),cst(1.),cst(0.))*(idv(u_E_)));
            //at_e+=on(_range=markedpoints(mesh_,"a1"), _rhs=lt_e, _element=u_E_[Component::X], _expr=cst(0.));
            //at_e+=on(_range=markedpoints(mesh_,"a2"), _rhs=lt_e, _element=u_E_[Component::Y], _expr=cst(0.));
            //at_e+=on(_range=markedpoints(mesh_,"a3"), _rhs=lt_e, _element=u_E_[Component::Z], _expr=cst(0.));
        }

        at_e.solve( _rhs = lt_e, _solution = u_E_);

        auto dispE = mean(_range=elements(support(XvE_)), _expr=idv(u_E_));
        std::cout << "Mean displacement x : " << dispE(0,0) << std::endl;
        std::cout << "Mean displacement y : " << dispE(1,0) << std::endl;
        std::cout << "Mean displacement z : " << dispE(2,0) << std::endl;
        auto rotE = mean(_range=elements(support(XvE_)), _expr = curlv(u_E_));
        std::cout << "Mean curl x : " << rotE(0,0) << std::endl;
        std::cout << "Mean curl y : " << rotE(1,0) << std::endl;
        std::cout << "Mean curl z : " << rotE(2,0) << std::endl;



        at_e.matrixPtr()->printMatlab(fmt::format("A{}.m",ts_E_->iteration()));
        lt_e.vectorPtr()->printMatlab(fmt::format("l{}.m",ts_E_->iteration()));

        ts_E_->updateFromDisp(u_E_);

        // Export
        this->exportResults(ts_E_->time());

        // Set to zero
        at_e.zero();
        lt_e.zero();

    }
}


// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::rotation()
{
    auto Id = eye<Dim,Dim>();
    auto iter = -1;
    // Inits
    this->initializeMesh();
    this->initializeParam();
    this->initializeFields();
    this->initializeTs_Exp();

    // Initialize rotation
    auto elem_const = [&]() 
    { 
        if constexpr(Dim == 2) 
            return XR_->element(); 
        else if constexpr(Dim == 3)
            return XvR_->element(); 
    };

    auto ts_const = [&]() 
    { 
        if constexpr(Dim == 2) 
            return bdf(_space = XR_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_ );
        else if constexpr(Dim == 3)
            return bdf(_space = XvR_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_ );

    };

    // Angular velocity
    std::string default_displ = std::string("{0.,0.,0.}");
    auto v_angular_ = elem_const();
    auto dt_v_angular_ = elem_const();

    if constexpr(Dim == 2)
    {
        v_angular_.on(_range=elements(support(XR_)), _expr= cst(0.));
        dt_v_angular_.on(_range=elements(support(XR_)), _expr= cst(0.));
    }
    else if constexpr(Dim == 3)
    {
        v_angular_.on(_range=elements(support(XvR_)), _expr= expr<Dim,1>(default_displ));
        dt_v_angular_.on(_range=elements(support(XvR_)), _expr= expr<Dim,1>(default_displ));
    }
    

    auto ts_V_angular_ = ts_const();
    ts_V_angular_->start();

    // Rotation Angle 
    auto theta_ = elem_const();

    if constexpr(Dim == 2)
        theta_.on(_range=elements(support(XR_)), _expr= cst(0.));
    else if constexpr(Dim == 3)
        theta_.on(_range=elements(support(XvR_)), _expr= expr<Dim,1>(default_displ));

    auto ts_Theta_ = ts_const();
    ts_Theta_->start();

    // Mass center
    auto massCenter = mean( _range = elements(support(XvR_)), _expr = P());
    auto massCenter_const = [&]()
    {
        if constexpr(Dim == 2) 
            return vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
        else if constexpr(Dim == 3)
            return vec(cst(massCenter(0,0)),cst(massCenter(1,0)),cst(massCenter(2,0)));
    };
    auto massCenterVec = massCenter_const();
    std::cout << "massCenter : " << massCenter(0,0) << ", " << massCenter(1,0) << std::endl;

    // Matrix of inertia
    auto momentOfInertia_const = [&]() 
    { 
        if constexpr(Dim == 2) 
            return integrate(_range=elements(support(XvR_)),_expr=cst(density_)*( inner(P()-massCenterVec) ) ).evaluate()(0,0);
        else if constexpr(Dim == 3)
        {
            auto rvec = P()-massCenterVec;
            return integrate(_range=elements(support(XvR_)),_expr=cst(density_)*( inner(rvec)*Id - rvec*trans(rvec) ) ).evaluate();
        } 
    };
    auto momentOfInertia = momentOfInertia_const();
    std::cout << "momentOfInertia : " << momentOfInertia << std::endl;

    // Bilinear and linear forms

    /*
        Assemblage : translation
    */
    auto a_translation_ = form2( _test = XvR_, _trial = XvR_ );
    auto l_translation_ = form1( _test = XvR_ );
    auto lt_translation_ = form1( _test = XvR_ );
    
    a_translation_.zero();
    l_translation_.zero();
    lt_translation_.zero();

    l_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst(density_)*trans(expr<Dim,1>( externalforce_ ))*id(u_translation_));
    a_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst(density_)*inner( ts_Translation_->polyDerivCoefficient()*idt(u_translation_),id( u_translation_ ) ) );

    /*
        Assemblage : rotation
    */
    auto a_const = [&]() 
    { 
        if constexpr(Dim == 2) 
            return form2( _test = XR_, _trial = XR_ );
        else if constexpr(Dim == 3)
            return form2( _test = XvR_, _trial = XvR_ );
    };

    auto l_const = [&]() 
    { 
        if constexpr(Dim == 2) 
            return form1( _test = XR_ );
        else if constexpr(Dim == 3)
            return form1( _test = XvR_ );
    };

    auto a_v_angular_ = a_const();
    auto l_v_angular_ = l_const();
    auto lt_v_angular_ = l_const();
    a_v_angular_.zero();
    l_v_angular_.zero();

    
    if constexpr(Dim == 2)
    {
        a_v_angular_ += integrate( _range = elements(support(XR_)), _expr = cst(momentOfInertia)*inner( ts_V_angular_->polyDerivCoefficient(0)*idt(v_angular_),id( v_angular_ ) ) );
        l_v_angular_ += integrate( _range = elements(support(XR_)), _expr = cst(momentOfInertia) * (cst(extFy_)*(Px()-massCenter(0,0)) - cst(extFx_)*(Py() - massCenter(1,0))) * id(v_angular_));
    }

    auto a_theta_ = a_const();
    auto l_theta_ = l_const();
    a_theta_.zero();
    l_theta_.zero();

    if constexpr(Dim == 2)
    {
        a_theta_ += integrate( _range = elements(support(XR_)), _expr =  inner( ts_Theta_->polyDerivCoefficient(0)*idt(theta_),id( theta_ ) ) );
    }
    
    /*
        Assemblage : elasticity
    */
    auto a_e = form2( _test = XvE_, _trial = XvE_ );
    auto at_e = form2( _test = XvE_, _trial = XvE_ );
    auto l_e = form1( _test = XvE_ );
    auto lt_e = form1( _test = XvE_ );

    a_e.zero();
    at_e.zero();
    l_e.zero();
    lt_e.zero();

    auto deft = sym(gradt(u_E_));
    auto def = sym(grad(u_E_));
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a_e += integrate( _range = elements(support(XvE_)), _expr = cst(density_) * inner( ts_E_->polyDerivCoefficient()*idt(u_E_),id( u_E_ ) ) + inner(sigmat,def));
    l_e = integrate( _range = elements(support(XvE_)), _expr = cst(density_) * trans( expr<Dim, 1>( externalforce_ ) ) * id( u_E_ ) );

    // Starting time loop
    for ( ts_Translation_->start(); ts_Translation_->isFinished() == false; ts_Translation_->next( u_translation_ ))
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_Translation_->time(),ts_Translation_->timeFinal()) << std::endl;

        /*
            Solve translation
        */
        lt_translation_ = l_translation_;
        lt_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst( density_ ) * inner( idv( ts_Translation_->polyDeriv() ), id( u_translation_ ) ) );
        a_translation_.solve( _rhs = lt_translation_, _solution = u_translation_, _rebuild = true );
        ts_Translation_->updateFromDisp(u_translation_);

        /*
            Solve rotation
        */
        if constexpr(Dim == 2)
        {
            lt_v_angular_ = l_v_angular_;
            lt_v_angular_ += integrate( _range = elements(support(XR_)), _expr = cst(momentOfInertia) * inner( idv( ts_V_angular_->polyDeriv()) , id( v_angular_ ) ) );
            // on donne une vitesse angulaire
            a_v_angular_+=on(_range= elements(support(XE_)), _rhs=lt_v_angular_, _element=v_angular_,_expr=cst(0.1));
        }
        a_v_angular_.solve( _rhs = lt_v_angular_, _solution = v_angular_, _rebuild = true);
        ts_V_angular_->next(v_angular_);

        double velocity = 0;
        if constexpr(Dim == 2)
        {
            velocity = mean(_range = elements(support(XR_)), _expr = idv(v_angular_))(0,0);     
            std::cout << "Angular velocity : " << velocity << std::endl; 
        } 

        if constexpr(Dim == 2)
        {
            l_theta_ += integrate( _range = elements(support(XR_)), _expr = idv(v_angular_) * id(theta_));
            l_theta_ += integrate( _range = elements(support(XR_)), _expr =  inner( idv( ts_Theta_->polyDeriv()) , id( theta_ ) ) );
        }
        a_theta_.solve( _rhs = l_theta_, _solution = theta_, _rebuild = true);
        ts_Theta_->next(theta_);

        double angle = 0;
        if constexpr(Dim == 2)
        {
            angle = mean(_range = elements(support(XR_)), _expr = idv(theta_))(0,0);    
            std::cout << "Rotation angle : " << angle << std::endl;
            auto rot = vec(cos(angle)*(Px() - massCenter(0,0)) - sin(angle)*(Py() - massCenter(1,0)) - Px() + massCenter(0,0), sin(angle)*(Px() - massCenter(0,0)) + cos(angle)*(Py() - massCenter(1,0)) - Py() + massCenter(1,0));
            u_rotation_.on(_range=elements(support(XvE_)),_expr =  rot);
        }

        /*
            Compute rigid motion
        */
        u_R_.on(_range=elements(support(XE_)),_expr =  idv(u_translation_) + idv(u_rotation_));

        /*
            Solve elasticity
        */
        lt_e = l_e;
        at_e = a_e;

        lt_e += integrate( _range=elements(support(XvE_)), _expr= cst(density_)*inner( idv(ts_E_->polyDeriv()),id( u_E_ ) ));
        lt_e += integrate( _range=elements(support(XvE_)), _expr= -cst(density_)*inner( idv(ts_Translation_->currentAcceleration()),id( u_E_ ) )); // translational acceleration
        

        if constexpr(Dim == 2)
        {
            ts_V_angular_->updateDerivative(v_angular_,dt_v_angular_);
            double acceleration = mean(_range = elements(support(XR_)), _expr = idv(dt_v_angular_))(0,0); 
            std::cout << "Rotation acceleration : " << acceleration << std::endl;

            auto coef_1 = - cst(acceleration) * sin(angle) - cst(velocity)*cst(velocity) * cos(angle);
            auto coef_2 = - cst(acceleration) * cos(angle) + cst(velocity)*cst(velocity) * sin(angle);
            auto coef_3 = cst(acceleration) * cos(angle) - cst(velocity)*cst(velocity) * sin(angle);
            auto coef_4 = - cst(acceleration) * sin(angle) - cst(velocity)*cst(velocity) * cos(angle);

            auto rot = vec(coef_1*(Px() - massCenter(0,0)) - coef_2*(Py() - massCenter(1,0)), coef_3*(Px() - massCenter(0,0)) - coef_4*(Py() - massCenter(1,0)));
            lt_e += integrate( _range=elements(support(XvE_)), _expr= - inner(rot, id( u_E_ ) ) ); // rotational acceleration
        }
                
        // Delete rigid motion 
        at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_E_, _expr=0.*one()); // translation
        
        // Rotation 
        if constexpr(Dim == 2)
        {
            at_e+=on(_range=markedpoints(mesh_,"rx"), _rhs=lt_e, _element=u_E_[Component::X], _expr=cst(0.));
            //at_e+=on(_range=markedpoints(mesh_,"rx2"), _rhs=lt_e, _element=u_E_[Component::X], _expr=cst(0.));
            at_e+=on(_range=markedpoints(mesh_,"ry"), _rhs=lt_e, _element=u_E_[Component::Y], _expr=cst(0.));
            //at_e+=on(_range=markedpoints(mesh_,"ry2"), _rhs=lt_e, _element=u_E_[Component::Y], _expr=cst(0.));
        }
        
        at_e.solve( _rhs = lt_e, _solution = u_E_ , _rebuild = true);

        at_e.matrixPtr()->printMatlab(fmt::format("A{}.m",ts_E_->iteration()));
        lt_e.vectorPtr()->printMatlab(fmt::format("l{}.m",ts_E_->iteration()));

        ts_E_->updateFromDisp(u_E_);
        ts_E_->next(u_E_);

        /*
            Compute total motion
        */
        u_total_.on(_range=elements(support(XvE_)), _expr=idv(u_R_) + idv(u_E_));

        // Export
        this->exportResults(ts_Translation_->time());

        // Set to zero
        lt_translation_.zero();
        lt_v_angular_.zero();
        l_theta_.zero();
        at_e.zero();
        lt_e.zero();

    }
}



// Export results
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::exportResults(double t)
{
    e_->step(t)->addRegions();
    e_->step(t)->add( "u_total_", u_total_ );
    e_->step(t)->add( "u_E_", u_E_ );
    e_->step(t)->add( "u_R_", u_R_ );
    e_->step(t)->add( "u_translation_", u_translation_ );
    e_->step(t)->add( "u_rotation_", u_rotation_);
    e_->save();
}


/*
// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::rotation()
{
    
    for ( ts_trans_->start(); ts_trans_->isFinished() == false; ts_trans_->next( u_trans ), ts_omega_->next( omega_ ), ts_theta_->next(theta_), ts_elas_->next())
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_trans_->time(),ts_trans_->timeFinal()) << std::endl;

        if constexpr(Dim == 2)
        {
            
            // Compute Xcm in current domain
            double angle = mean(_range = elements(support(VTheta_)), _expr = idv(theta_))(0,0);  
            std::cout << "Current angle : " << angle << std::endl; 

            auto xC = cos(angle) * (Px() - xcmI) - sin(angle) * (Py() - ycmI) + xcmI ;
            auto yC = sin(angle) * (Px() - xcmI) + cos(angle) * (Py() - ycmI) + ycmI ; 

            auto xcmCurr = mean( _range =  elements(support(VUTheta_)), _expr = vec( xC, yC ));
            auto xcmC = xcmCurr(0,0);
            auto ycmC = xcmCurr(1,0);
            std::cout << "Mass center current domain: " << xcmC << ", " << ycmC << std::endl;
            
            // Compute current J
            auto J_expr = cst(rho_) * ((xC-xcmC)*(xC-xcmC) + (yC-ycmC)*(yC-ycmC));
            auto J = integrate( _range = elements(support(VTheta_)), _expr = J_expr).evaluate()(0,0);
            std::cout << "J : " << J << std::endl;

            // Compute T
            auto T_expr = cst(rho_) * (cst(extTy)*(xC-xcmC) - cst(extTx)*(yC - ycmC));
            auto T = integrate( _range = elements(support(VTheta_)), _expr = T_expr).evaluate()(0,0); 
            std::cout << "T : " << T << std::endl;
            
            // Compute omega
            auto bdf_poly_omega = ts_omega_->polyDeriv();
            auto extrap_omega = ts_omega_->poly();

            a_omega_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * inner( ts_omega_->polyDerivCoefficient(0)*idt(omega_),id( omega_ ) ) );
            //a_omega_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * trans(gradt(omega_)*idv(extrap_omega))*id( omega_ )  );
            l_omega_ += integrate( _range = elements(support(VTheta_)), _expr = T_expr * id(omega_));
            l_omega_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * inner( idv( bdf_poly_omega) , id( omega_ ) ) );
            a_omega_+=on(_range= elements(support(VTheta_)), _rhs=l_omega_, _element=omega_,_expr=cst(avel));
            a_omega_.solve( _rhs = l_omega_, _solution = omega_, _rebuild = true);

            double velAng = mean(_range = elements(support(VTheta_)), _expr = idv(omega_))(0,0);     
            std::cout << "velAng : " << velAng << std::endl;  

            // Compute theta
            auto bdf_poly_theta = ts_theta_->polyDeriv();
            auto extrap_theta = ts_theta_->poly();

            a_theta_ += integrate( _range = elements(support(VTheta_)), _expr =  inner( ts_theta_->polyDerivCoefficient(0)*idt(theta_),id( theta_ ) ) );
            //a_theta_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * trans(gradt(theta_)*idv(extrap_theta))*id( theta_ )  );
            l_theta_ += integrate( _range = elements(support(VTheta_)), _expr = idv(omega_) * id(theta_));
            l_theta_ += integrate( _range = elements(support(VTheta_)), _expr =  inner( idv( bdf_poly_theta) , id( theta_ ) ) );
            a_theta_.solve( _rhs = l_theta_, _solution = theta_, _rebuild = true);

            // Compute rotation matrix
            double angleN = mean(_range = elements(support(VTheta_)), _expr = idv(theta_))(0,0);     
            std::cout << "AngleN : " << angleN << std::endl;  
            auto rot = vec(cos(angleN)*(Px() - xcmI) - sin(angleN)*(Py() - ycmI) - Px() + xcmI, sin(angleN)*(Px() - xcmI) + cos(angleN)*(Py() - ycmI) - Py() + ycmI);

            // Compute new displacement
            u_t.on(_range=elements(support(VUTheta_)),_expr =  rot);

            // Compute translation
            lt_trans_ = l_trans_;
            lt_trans_ += integrate( _range = elements(support(VUTrans_)), _expr = cst( rho_ ) * inner( idv( ts_trans_->polyDeriv() ), id( u_trans ) ) );
            a_trans_.solve( _rhs = lt_trans_, _solution = u_trans, _rebuild = true );
            
            // Compute rigid displacement
            u_rigid.on(_range=elements(support(VUTheta_)),_expr =  idv(u_trans) + idv(u_t));

            lt_e = l_e;
            at_e = a_e;

            lt_e += integrate( _range=elements(support(VUElas_)), _expr= cst(rho_)*inner( idv(ts_elas_->polyDeriv()),id( u_elas ) ));
            lt_e += integrate( _range=elements(support(VUElas_)), _expr= -cst(rho_)*inner( idv(ts_trans_->currentAcceleration()),id( u_elas ) ));
            //lt_e += integrate( _range=elements(support(VUElas_)), _expr= J_expr * idv( ts_omega_->polyDeriv()) * id( u_elas ) ) ;
            
            at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_elas, _expr=0.*one());
            at_e.solve( _rhs = lt_e, _solution = u_elas , _rebuild = true);

            u_total.on(_range=elements(support(VUTheta_)), _expr=idv(u_rigid) + idv(u_elas));

            e_->step(ts_trans_->time())->addRegions();
            e_->step(ts_trans_->time())->add( "utheta", u_t );
            e_->step(ts_trans_->time())->add( "u_trans", u_trans );
            e_->step(ts_trans_->time())->add( "u_rigid", u_rigid );
            e_->step(ts_trans_->time())->add( "u_elas", u_elas );
            e_->step(ts_trans_->time())->add( "u_total", u_total );
            e_->save();

            ts_trans_->updateFromDisp(u_trans);
            ts_elas_->updateFromDisp(u_elas);

            a_omega_.zero();
            l_omega_.zero();
            a_theta_.zero();
            l_theta_.zero();
            lt_trans_.zero();
            lt_e.zero();
            at_e.zero();
        }
    
    }
}
*/

/*
// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::rotation()
{
    // Init parameters
    std::cout << "init parameters" << std::endl;
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

    std::size_t offset = 0;
    double extFx = 0;
    double extFy = 0;

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
                extFx = std::stod(&externalforce_[1],&offset);
                extFy = std::stod(&externalforce_[offset+2]);
            }
        }
    }

    // External torque
    double extTx = 0.;
    double extTy = 0.;
    if ( specs_["/Models/LinearElasticity"_json_pointer].contains("torque") )
    {
        for ( auto [key, loading] : specs_["/Models/LinearElasticity/torque"_json_pointer].items() )
        {
            std::string loadtype = fmt::format( "/Models/LinearElasticity/torque/{}/type", key );

            if ( specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Gravity" )
            {
                LOG( INFO ) << fmt::format( "Loading {}: Gravity found", key );
                std::string loadexpr = fmt::format( "/Models/LinearElasticity/torque/{}/parameters/expr", key );
                auto torque = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
                extTx = std::stod(&torque[1],&offset);
                extTy = std::stod(&torque[offset+2]);
            }
        }
    }

    // Space
    std::cout << "init spaces" << std::endl;
    auto VTheta_ = Pch<0>(mesh_, markedelements(mesh_, "Caoutchouc"));
    auto VUTheta_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );
    auto VUTrans_ = Pchv<0>(mesh_, markedelements(mesh_, "Caoutchouc"));
   
    // Exporter
    std::cout << "init exporter" << std::endl;
    auto e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );

    // Initialize Newmark scheme
    std::cout << "init newmark" << std::endl;
    double initial_time = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    double gamma = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);

    // Rotation
    auto theta_ = VTheta_->element();
    theta_.on(_range=elements(support(VTheta_)), _expr= cst(0.));

    tst_ptrtype ts_t_ = newmark(_space = VTheta_, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _gamma=gamma, _beta=beta );
    ts_t_->start();
    ts_t_->initialize( theta_ );
    ts_t_->updateFromDisp(theta_);

    // Rotational displacement
    auto u_t = VUTheta_->element();
    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    u_t.on(_range=elements(support(VUTheta_)), _expr=expr<Dim,1>(default_displ));

    // Translation displacement
    auto u_trans = VUTrans_->element();
    u_trans.on(_range=elements(support(VUTrans_)), _expr=expr<Dim,1>(default_displ));
    
    tsr_ptrtype ts_trans_ =  newmark(_space = VUTrans_, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _gamma=gamma, _beta=beta );
    ts_trans_->start();
    ts_trans_->initialize( u_trans );
    ts_trans_->updateFromDisp(u_trans);
    
    // Rigid displacement
    auto u_rigid = VUTheta_->element();
    u_rigid.on(_range=elements(support(VUTheta_)), _expr=idv(u_trans) + idv(u_t));

    // Export at init time
    e_->step(0)->addRegions();
    e_->step(0)->add( "utheta", u_t );
    e_->step(0)->add( "u_trans", u_trans );
    e_->step(0)->add( "u_rigid", u_rigid );
    e_->save();

    // Init forms
    std::cout << "init forms" << std::endl;
    auto a_theta_ = form2( _test = VTheta_, _trial = VTheta_ );
    auto l_theta_ = form1( _test = VTheta_ );
    a_theta_.zero();
    l_theta_.zero();

    auto a_trans_ = form2( _test = VUTrans_, _trial = VUTrans_ );
    auto l_trans_ = form1( _test = VUTrans_ );
    auto lt_trans_ = form1( _test = VUTrans_ );
    a_trans_.zero();
    l_trans_.zero();
    lt_trans_.zero();

    l_trans_ += integrate( _range = elements(support(VUTrans_)), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_trans));
    a_trans_ += integrate( _range = elements(support(VUTrans_)), _expr = cst(rho_)*inner( ts_trans_->polyDerivCoefficient()*idt(u_trans),id( u_trans ) ) );

    auto xcmInit = mean( _range =  elements(support(VUTheta_)), _expr = P());
    auto xcmI = xcmInit(0,0);
    auto ycmI = xcmInit(1,0);
    std::cout << "Mass center reference domain : " << xcmI << ", " << ycmI << std::endl;

    std::cout << "start time loop" << std::endl;
    for ( ts_t_->start(); ts_t_->isFinished() == false; ts_t_->next( theta_ ), ts_trans_->next( u_trans ))
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_t_->time(),ts_t_->timeFinal()) << std::endl;

        if constexpr(Dim == 2)
        {
            
            // Compute Xcm in current domain
            double angle = mean(_range = elements(support(VTheta_)), _expr = idv(theta_))(0,0);  
            std::cout << "Current angle : " << angle << std::endl; 

            auto xC = cos(angle) * (Px() - xcmI) - sin(angle) * (Py() - ycmI) + xcmI ;
            auto yC = sin(angle) * (Px() - xcmI) + cos(angle) * (Py() - ycmI) + ycmI ; 

            auto xcmCurr = mean( _range =  elements(support(VUTheta_)), _expr = vec( xC, yC ));
            auto xcmC = xcmCurr(0,0);
            auto ycmC = xcmCurr(1,0);
            std::cout << "Mass center current domain: " << xcmC << ", " << ycmC << std::endl;
            
            // Compute current J
            auto J_expr = cst(rho_) * ((xC-xcmC)*(xC-xcmC) + (yC-ycmC)*(yC-ycmC));
            auto J = integrate( _range = elements(support(VTheta_)), _expr = J_expr).evaluate()(0,0);
            std::cout << "J : " << J << std::endl;

            // Compute T
            auto T_expr = cst(rho_) * (cst(extTy)*(xC-xcmC) - cst(extTx)*(yC - ycmC));
            auto T = integrate( _range = elements(support(VTheta_)), _expr = T_expr).evaluate()(0,0); 
            std::cout << "T : " << T << std::endl;
            
            // Compute theta
            a_theta_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * inner( ts_t_->polyDerivCoefficient()*idt(theta_),id( theta_ ) ) );
            l_theta_ += integrate( _range = elements(support(VTheta_)), _expr = T_expr * id(theta_));
            l_theta_ += integrate( _range = elements(support(VTheta_)), _expr = J_expr * inner( idv( ts_t_->polyDeriv() ), id( theta_ ) ) );
            a_theta_.solve( _rhs = l_theta_, _solution = theta_, _rebuild = true);

            // Compute rotation matrix
            double angleN = mean(_range = elements(support(VTheta_)), _expr = idv(theta_))(0,0);     
            std::cout << "AngleN : " << angleN << std::endl;  
            auto rot = vec(cos(angleN)*(Px() - xcmI) - sin(angleN)*(Py() - ycmI) - Px() + xcmI, sin(angleN)*(Px() - xcmI) + cos(angleN)*(Py() - ycmI) - Py() + ycmI);

            // Compute new displacement
            u_t.on(_range=elements(support(VUTheta_)),_expr =  rot);

            // Compute translation
            lt_trans_ = l_trans_;
            lt_trans_ += integrate( _range = elements(support(VUTrans_)), _expr = cst( rho_ ) * inner( idv( ts_trans_->polyDeriv() ), id( u_trans ) ) );
            a_trans_.solve( _rhs = lt_trans_, _solution = u_trans, _rebuild = true );
            
            // Compute rigid displacement
            u_rigid.on(_range=elements(support(VUTheta_)),_expr =  idv(u_trans) + idv(u_t));

            e_->step(ts_t_->time())->addRegions();
            e_->step(ts_t_->time())->add( "utheta", u_t );
            e_->step(ts_t_->time())->add( "u_trans", u_trans );
            e_->step(ts_t_->time())->add( "u_rigid", u_rigid );
            e_->save();

            ts_t_->updateFromDisp(theta_);
            ts_trans_->updateFromDisp(u_trans);

            a_theta_.zero();
            l_theta_.zero();
            lt_trans_.zero();
        }
    
    }
}
*/



// Initialization contact parameters 
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initializeContact()
{
    // Initialize contact field
    contactRegion_ =  project(_space=XE_, _range=elements(support(XE_)), _expr = cst(0.));
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
}


template <int Dim, int Order>
Range<typename ElasticRigid<Dim, Order>::mesh_t, MESH_FACES>
ElasticRigid<Dim, Order>::getContactRegion(elementv_t_E const& u)
{
    Range<mesh_t,MESH_FACES> myelts(mesh_ );

    auto defv = sym(gradv(u));
    auto Id = eye<Dim,Dim>();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();

    contactRegion_.on( _range=elements(support(XvE_)), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));

    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=XE_, _trial=XE_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(support(XE_)) )
    {
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : XE_->dof()->faceLocalDof( face.id() ) )
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



template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initG()
{
    // Init the distance fields
    g_ = XE_->element();
    g_.on(_range=elements(support(XE_)), _expr=cst(100000.));

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
                        for (auto const& ldof  : XE_->dof()->faceLocalDof( face.id() ))
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
                for (auto const& ldof : XE_->dof()->faceLocalDof(faceIDs[fid]))
                    g_[ldof.index()] = rir.distance() - tolDistance_;
            }
            
        }
    }

}


// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::translationContact()
{
    auto Id = eye<Dim,Dim>();
    auto dir = oneZ();

    // Inits
    this->initializeMesh();
    this->initializeParam();
    this->initializeFields();
    this->initializeTs_Exp();

    // Contact inits
    this->initializeContact();
    this->initG();

    // Linear and bilinear forms

    /*
        Assemblage tranlsation
    */

    auto a_translation_ = form2( _test = XvR_, _trial = XvR_ );
    auto l_translation_ = form1( _test = XvR_ );
    auto lt_translation_ = form1( _test = XvR_ );
    
    a_translation_.zero();
    l_translation_.zero();
    lt_translation_.zero();

    l_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst(density_)*trans(expr<Dim,1>( externalforce_ ))*id(u_translation_));
    a_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst(density_)*inner( ts_Translation_->polyDerivCoefficient()*idt(u_translation_),id( u_translation_ ) ) );

    /*
        Assemblage elasticity
    */

    auto a_e = form2( _test = XvE_, _trial = XvE_ );
    auto at_e = form2( _test = XvE_, _trial = XvE_ );
    auto l_e = form1( _test = XvE_ );
    auto lt_e = form1( _test = XvE_ );

    a_e.zero();
    at_e.zero();
    l_e.zero();
    lt_e.zero();

    auto deft = sym(gradt(u_E_));
    auto def = sym(grad(u_E_));
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a_e += integrate( _range = elements(support(XvE_)), _expr = cst(density_) * inner( ts_E_->polyDerivCoefficient()*idt(u_E_),id( u_E_ ) ) + inner(sigmat,def));
    l_e = integrate( _range = elements(support(XvE_)), _expr = cst(density_) * trans( expr<Dim, 1>( externalforce_ ) ) * id( u_E_ ) );

    // output 
    std::ofstream ofs("outputs.csv");
    ofs << fmt::format( "time, evaluateStress, evaluateDisp") << std::endl;

    // Start time loop
    for ( ts_Translation_->start(); ts_Translation_->isFinished() == false; ts_Translation_->next( u_translation_ ))
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_Translation_->time(),ts_Translation_->timeFinal()) << std::endl;

        if (Environment::isMasterRank())
            std::cout << "***** Compute new contact region *****" << std::endl;

        myelts_ = getContactRegion(u_total_);
        
        if (Environment::isMasterRank())
            std::cout << " Nbr faces for processContact : " << nbrFaces_ << std::endl;
                
        /*
            Solve translation
        */

        lt_translation_ = l_translation_;
        lt_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst( density_ ) * inner( idv( ts_Translation_->polyDeriv() ), id( u_translation_ ) ) );
        if (nbrFaces_ > 0)//Add contact terms
        {
            auto epsv = sym(gradv(u_E_));
            auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
            lt_translation_ += integrate( _range = elements(support(XvR_)), _expr = inner(sigmav, id( u_translation_ ) ) );
        }
        a_translation_.solve( _rhs = lt_translation_, _solution = u_translation_, _rebuild = true );
        ts_Translation_->updateFromDisp(u_translation_);

    
        /*
            Solve elasticity
        */

        lt_e = l_e;
        at_e = a_e;

        lt_e += integrate( _range=elements(support(XvE_)), _expr= cst(density_)*inner( idv(ts_E_->polyDeriv()),id( u_E_ ) ));
        lt_e += integrate( _range=elements(support(XvE_)), _expr= -cst(density_)*inner( idv(ts_Translation_->currentAcceleration()),id( u_E_ ) )); // translational acceleration
        
        if (nbrFaces_ > 0) // add contact terms
        {
            
            if (method_.compare("penalty") == 0)
            {
                at_e += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u_E_),trans(expr<Dim,1>(direction_))*id(u_E_)) );
                lt_e += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_translation_), trans(expr<Dim,1>(direction_))*id(u_E_)) );
            }
            else if (method_.compare("nitsche") == 0)
            {
                auto deft = sym(gradt(u_E_));
                auto def = sym(grad(u_E_));
                auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
                auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

                at_e += integrate (_range=myelts_,_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)); 
                at_e += integrate (_range=myelts_,_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u_E_) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u_E_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));
                lt_e += integrate (_range=myelts_,_expr= inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_translation_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u_E_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));     
            }
        }

        // deleting translation
        at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_E_, _expr=0.*one());

        at_e.solve( _rhs = lt_e, _solution = u_E_ , _rebuild = true);
        ts_E_->updateFromDisp(u_E_);
        ts_E_->next(u_E_);

        /*
            Solve total displacement
        */
        u_total_.on(_range=elements(support(XvE_)), _expr=idv(u_translation_) + idv(u_E_));

        // Reset
        lt_translation_.zero();
        at_e.zero();
        lt_e.zero();

        // Exports
        this->exportResults(ts_Translation_->time());

        
        auto ctx = XE_->context();
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

        auto epsv = sym(gradv(u_E_));
        auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
        auto contactPressure = XE_->element();
        contactPressure.on(_range=myelts_, _expr = trans(expr<Dim,1>(direction_))*sigmav);
        auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure) ); 

        auto evaluateDispExpr = XE_->element();
        evaluateDispExpr.on(_range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_total_));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
    

        ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",ts_E_->time(),evaluateStress(0,0), evaluateDisp(0,0)) << std::endl;

    }
    ofs.close();
}


/*
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::timeLoopFixedPoint()
{
    // Elasticity equations   

    auto a_e = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_e = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_tmp_e = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_e = form1( _test = Xhv_ );
    auto lt_e = form1( _test = Xhv_ );
    auto lt_tmp_e = form1( _test = Xhv_ );

    a_e.zero();
    at_e.zero();
    at_tmp_e.zero();
    l_e.zero();
    lt_e.zero();
    lt_tmp_e.zero();

    auto deft = sym(gradt(u_e_));
    auto def = sym(grad(u_e_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a_e += integrate( _range = elements(support(Xhv_)), _expr = cst(rho_)*inner( ts_e_->polyDerivCoefficient()*idt(u_e_),id( u_e_ ) ) + inner(sigmat,def));
    l_e = integrate( _range = elements(support(Xhv_)), _expr = cst( rho_ ) * trans( expr<Dim, 1>( externalforce_ ) ) * id( u_e_ ) );

    // Rigid equations
    auto a_r_ = form2( _test = VhvC_, _trial = VhvC_ );
    auto at_r_ = form2( _test = VhvC_, _trial = VhvC_ );
    auto l_r_ = form1( _test = VhvC_ );
    auto lt_r_ = form1( _test = VhvC_ );
    auto lt_tmp_r_ = form1( _test = VhvC_ );

    a_r_.zero();
    at_r_.zero();
    l_r_.zero();
    lt_r_.zero();
    lt_tmp_r_.zero();

    auto dir = oneZ();//-expr<Dim,1>(direction_);

    l_r_ += integrate( _range = elements(support(VhvC_)), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_r_));
    a_r_ += integrate( _range = elements(support(VhvC_)), _expr = cst(rho_)*inner( ts_r_->polyDerivCoefficient()*idt(u_r_),id( u_r_ ) ) );

    // Outputs
    std::ofstream ofs("outputs.csv");
    ofs << fmt::format( "time, evaluateStress, evaluateDisp") << std::endl;

    int iteration = 0;
    for ( ts_e_->start(); ts_e_->isFinished() == false; ts_e_->next( u_e_ ), ts_r_->next( u_r_ ) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_e_->time(),ts_e_->timeFinal()) << std::endl;

        if (Environment::isMasterRank())
            std::cout << "***** Compute new contact region *****" << std::endl;
        myelts_ = getContactRegion(u_);
        
        if (Environment::isMasterRank())
            std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;

        if (nbrFaces_ == 0)
        {
            
            lt_r_ = l_r_;
            at_r_ = a_r_;

            lt_r_ += integrate( _range = markedelements( support( VhvC_ ), "Caoutchouc" ), _expr = cst( rho_ ) * inner( idv( ts_r_->polyDeriv() ), id( u_r_ ) ) );

            at_r_.solve( _rhs = lt_r_, _solution = u_r_, _rebuild = true );
            ts_r_->updateFromDisp(u_r_);

            lt_e = l_e;
            at_e = a_e;

            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= cst(rho_)*inner( idv(ts_e_->polyDeriv()),id( u_e_ ) ));
            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= -cst(rho_)*inner( idv(ts_r_->currentAcceleration()),id( u_e_ ) ));
            
            at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_e_, _expr=0.*one());
            at_e.solve( _rhs = lt_e, _solution = u_e_ , _rebuild = true);
            ts_e_->updateFromDisp(u_e_);

            u_.on(_range=elements(support(Xhv_)),_expr=idv(u_e_) + idv(u_r_));
        }
        else
        {
            
            int fixedPointIteration = 0;
            double fixedPointerror = 0.;

            lt_r_ = l_r_;
            at_r_ = a_r_;

            lt_r_ += integrate( _range = markedelements( support( VhvC_ ), "Caoutchouc" ), _expr = cst( rho_ ) * inner( idv( ts_r_->polyDeriv() ), id( u_r_ ) ) );

            lt_e = l_e;
            at_e = a_e;

            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= cst(rho_)*inner( idv(ts_e_->polyDeriv()),id( u_e_ ) ));
            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= -cst(rho_)*inner( idv(ts_r_->currentAcceleration()),id( u_e_ ) ));

            auto u_e_tmp =  Xhv_->element();
            u_e_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_e_));

            auto u_e_tmpNew = Xhv_->element();
            u_e_tmpNew.on( _range=elements(support(Xhv_)), _expr = idv(u_e_)); 

            auto u_r_tmp =  VhvC_->element();
            u_r_tmp.on( _range=elements(support(VhvC_)), _expr = idv(u_r_));

            auto u_r_tmpNew = VhvC_->element();
            u_r_tmpNew.on( _range=elements(support(VhvC_)), _expr = idv(u_r_)); 

            auto u_tmp =  Xhv_->element();
            u_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_));

            auto u_tmpNew = Xhv_->element();
            u_tmpNew.on( _range=elements(support(Xhv_)), _expr = idv(u_)); 
            
        
            while ((fixedPointerror > fixedPointtol_) || (fixedPointIteration < 1))
            {
            
                if (Environment::isMasterRank())
                    std::cout << "Fixed point iteration : " << fixedPointIteration << std::endl;
                u_e_tmp.on( _range=elements(mesh_), _expr = idv(u_e_tmpNew)); ;
                u_r_tmp.on( _range=elements(mesh_), _expr = idv(u_r_tmpNew)); ;
                u_tmp.on( _range=elements(mesh_), _expr = idv(u_tmpNew)); ;
            
                if (Environment::isMasterRank())
                    std::cout << "***** Compute new contact region *****" << std::endl;
                myelts_ = getContactRegion(u_);
                
                if (Environment::isMasterRank())
                    std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
 
                lt_tmp_r_ = lt_r_;

                //auto epsv = sym(gradv(u_e_tmp));
                auto epsv = sym(gradv(u_tmp));
                auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();

                lt_tmp_r_ +=  integrate( _range=myelts_, _expr= inner( sigmav,id( u_r_tmp )));
                at_r_.solve( _rhs = lt_tmp_r_, _solution = u_r_tmpNew, _rebuild = true );

                
                at_tmp_e = at_e;
                lt_tmp_e = lt_e;

                at_tmp_e += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_tmp),trans(expr<Dim,1>(direction_))*id(u_e_tmp)) );
                lt_tmp_e += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_r_tmp), trans(expr<Dim,1>(direction_))*id(u_e_tmp)) );
                        
                at_tmp_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_tmp_e, _element=u_e_tmp, _expr=0.*one());
                at_tmp_e.solve( _rhs = lt_tmp_e, _solution = u_e_tmpNew , _rebuild = true);

                u_tmpNew.on(_range=elements(support(Xhv_)),_expr=idv(u_e_tmpNew) + idv(u_r_tmpNew));

                fixedPointerror = integrate(_range=elements(support(Xhv_)), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(support(Xhv_)),_expr=norm2(idv(u_))).evaluate()(0,0); 

                if (Environment::isMasterRank())
                    std::cout << "Error fixed point : " << fixedPointerror << std::endl;
                fixedPointIteration++;

                if (fixedPointIteration == 10)
                    break;
            
                // Reset
                lt_tmp_r_.zero();
                at_tmp_e.zero();
                lt_tmp_e.zero();
            
            }

            u_e_.on( _range=elements(support(Xhv_)), _expr = idv(u_e_tmpNew));
            u_r_.on( _range=elements(support(VhvC_)), _expr = idv(u_r_tmpNew));
            u_.on( _range=elements(support(Xhv_)), _expr = idv(u_tmpNew)); 


            ts_r_->updateFromDisp(u_r_);
            ts_e_->updateFromDisp(u_e_);
        }

        // Reset
        at_r_.zero();
        lt_r_.zero();
        at_e.zero();
        lt_e.zero();


        // Exports
        myelts_ = getContactRegion(u_);

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

        auto epsv = sym(gradv(u_e_));
        auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();

        auto contactPressure = Xh_->element();
        contactPressure.on(_range=myelts_, _expr = trans(expr<Dim,1>(direction_))*sigmav);
        auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure) ); 

        auto evaluateDispExpr = Xh_->element();
        evaluateDispExpr.on(_range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
    
        if ( Dim == 2 )
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",
                                ts_e_->time(), evaluateStress(0,0), evaluateDisp(0,0))
                << std::endl;
        else
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",
                                ts_e_->time(), evaluateStress(0,0), evaluateDisp(0,0) )
                << std::endl;

        // Export
        this->exportResults( ts_e_->time() );

        iteration++;
    }
    ofs.close();
    
}
*/

/*
// Run method
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::run()
{
    if (Environment::isMasterRank())
        std::cout << "***** Initialize parameters *****" << std::endl;
    initialize();

    if (Environment::isMasterRank())
        std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    if (Environment::isMasterRank())
        std::cout <<  "***** Initialize distance g *****" << std::endl;
    initG();

    this->exportResults(0.);

    if (Environment::isMasterRank())
        std::cout <<  "***** Start time loop *****" << std::endl;
    if (fixedPoint_ == 1)
        timeLoopFixedPoint();
    else 
        timeLoop();
}
*/

}