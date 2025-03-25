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
    void initializeContact();
    void initG();
    Range<mesh_t, MESH_FACES> getContactRegion(elementv_t_E const& u);

    // Run
    void run();

    // Export
    void exportResults(double t);
    

private:
    // Constructeur
    nl::json specs_;

    // Mesh
    double H_;
    std::shared_ptr<mesh_t> mesh_;

    // Spaces
    spacev_ptr_t_E Xv_;
    space_ptr_t_E X_;
    spacev_ptr_t_R XvR_;
    
    // Param Solid
    double density_, mass_;
    double E_, nu_, lambda_, mu_;
    std::string externalforce_, neumannUpper_, neumannLower_;
    double extFx_,extFy_;
    Eigen::Vector2d Force;
    Eigen::Vector2d ForceContact;
    double Upper_x, Upper_y, Lower_x, Lower_y;

    // Translation
    Eigen::Vector2d u_trans;
    Eigen::Vector2d dt_u_trans;
    Eigen::Vector2d dtt_u_trans;
    Eigen::Vector2d dtt_u_trans_old;
    elementv_t_R u_trans_true;


    // Rotation
    double omega;
    double theta;
    double acceleration;

    // Fields
    elementv_t_E u_elastic;
    elementv_t_E u_theta;
    elementv_t_E dt_u_theta;
    elementv_t_E dtt_u_theta;
    elementv_t_E rot_a;
    elementv_t_E dtt_u_theta_old;
    elementv_t_E u_rigid;
    elementv_t_E u_total;
    elementv_t_E u_true;

    // Newmark schemes
    double initial_time_, final_time_, time_step_;
    double gamma_, beta_;

    // Exporter
    exporter_ptrtype e_;
    ts_ptrtype_E ts_;
    ts_ptrtype_E ts_true;
    ts_ptrtype_Translation ts_translation;
    


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
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    // Spaces
    Xv_  = Pchv<Order>( mesh_, markedelements( mesh_, "Solid" ) );
    X_ = Pch<Order>(mesh_, markedelements( mesh_, "Solid" ) );
    XvR_ = Pchv<0>(mesh_, markedelements( mesh_, "Solid" ) );
}


// Initialization body parameters
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeParam()
{
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    density_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    mass_ = integrate( _range = elements(support(X_)), _expr = cst(density_) ).evaluate()(0,0);

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
    Force.setZero();
    Force[0] = extFx_;
    Force[1] = extFy_;
    

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

// Initialization contact parameters 
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initializeContact()
{
    // Initialize contact field
    contactRegion_ =  project(_space=X_, _range=elements(support(X_)), _expr = cst(0.));
    contactFaces_ = project(_space=X_, _range=elements(support(X_)), _expr = cst(0.));
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

    ForceContact.setZero();

}


// Initialization exporter and newmark time scheme
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeTs_Exp()
{
    // Exporter
    e_ = Feel::exporter(_mesh = mesh_, _name = "initial", _geo="change" );

    // Initialize Newmark scheme
    initial_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    final_time_ = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    time_step_ = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    gamma_ = get_value(specs_, "/TimeStepping/LinearElasticity/gamma", 0.5);
    beta_ = get_value(specs_, "/TimeStepping/LinearElasticity/beta", 0.25);

    ts_ =  newmark(_space = Xv_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_, _gamma=gamma_, _beta=beta_ );
    ts_->start();
    ts_->initialize( u_elastic );    
    ts_->updateFromDisp(u_elastic);

    ts_true =  newmark(_space = Xv_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_, _gamma=gamma_, _beta=beta_ );
    ts_true->start();
    ts_true->initialize( u_true );    
    ts_true->updateFromDisp(u_true);

    ts_translation = newmark(_space = XvR_, _initial_time=initial_time_, _final_time=final_time_, _time_step=time_step_, _gamma=gamma_, _beta=beta_ );
    ts_translation->start();
    ts_translation->initialize( u_trans_true );    
    ts_translation->updateFromDisp(u_trans_true);

    
}

// Initialization displacement fields
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initializeFields()
{
    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");


    u_trans.setZero();
    dt_u_trans.setZero();
    dtt_u_trans.setZero();
    dtt_u_trans_old.setZero();

    u_trans_true = XvR_->element();
    u_trans_true.on(_range=elements(support(XvR_)), _expr= expr<Dim,1>(default_displ));  

    omega = 0.;
    theta = 0;
    acceleration = 0;

    u_elastic = Xv_->element();
    u_elastic.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));  

    u_theta = Xv_->element();
    u_theta.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));   

    dt_u_theta = Xv_->element();
    dt_u_theta.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));  

    dtt_u_theta = Xv_->element(); 
    dtt_u_theta.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));  

    rot_a = Xv_ -> element();
    rot_a.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));  

    dtt_u_theta_old = Xv_->element();
    dtt_u_theta_old.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));   

    u_rigid = Xv_->element();
    u_rigid.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));  

    u_total = Xv_->element(); 
    u_total.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));    

    u_true = Xv_->element();
    u_true.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>(default_displ));
}

template <int Dim, int Order>
Range<typename ElasticRigid<Dim, Order>::mesh_t, MESH_FACES>
ElasticRigid<Dim, Order>::getContactRegion(elementv_t_E const& u)
{
    Range<mesh_t,MESH_FACES> myelts(mesh_);

    contactFaces_ = project(_space=X_, _range=elements(support(X_)), _expr = cst(0.));
    contactRegion_ = project(_space = X_,  _range = elements(support(X_)), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));
    
    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=X_, _trial=X_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(support(X_)) )
    {
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : X_->dof()->faceLocalDof( face.id() ) )
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

template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initG()
{
    // Init the distance fields
    g_ = X_->element();
    g_.on(_range=elements(support(X_)), _expr=cst(100.));


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
                        for (auto const& ldof  : X_->dof()->faceLocalDof( face.id() ))
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
}


// Export results
template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::exportResults(double t)
{
    std::cout << "Export" << std::endl;
    e_->step(t)->setMesh(mesh_);
    e_->step(t)->add( "u_total", idv(u_total) );
    e_->step(t)->add("u_rigid", idv(u_rigid) );
    e_->step(t)->add("u_theta", idv(u_theta) );
    e_->step(t)->add("u_elastic",u_elastic);
    e_->step(t)->add("u_true",u_true);
    e_->step(t)->add("u_trans_true", u_trans_true);
    e_->step(t)->add( "dtt_u_theta", idv(dtt_u_theta) );
    e_->step(t)->add("rot_a", rot_a);
    // Contact
    e_->step(t)->add("g_",g_);
    e_->step(t)->add("contactRegion_",contactRegion_);
    e_->step(t)->add("contactFaces_",contactFaces_);
    e_->save(); 
}

template <int Dim, int Order>
void ElasticRigid<Dim, Order>::run()
{
    
    if constexpr(Dim == 2)
    {
        // Init
        std::cout << "Init mesh" << std::endl;
        this->initializeMesh();
        std::cout << "Init param solid" << std::endl;
        this->initializeParam();
        std::cout << "Init fields" << std::endl;
        this->initializeFields();
        std::cout << "Initialize contact" << std::endl;
        this->initializeContact();
        std::cout << "Init distance" << std::endl;
        this->initG();
        std::cout << "Init exporter and time data" << std::endl;
        this->initializeTs_Exp();
        
        
        // Export
        this->exportResults(0);

        // Starting time loop
        int iter = 1;
        auto Id = eye<Dim,Dim>();
        
        // Construct backend
        auto backendElastic = backend( _name="elastic", _worldcomm = Xv_->worldCommPtr() );
        auto backendRigid = backend( _name="rigid", _worldcomm = Xv_->worldCommPtr() );
        
        std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(backendElastic,qsNullSpace(Xv_,mpl::int_<Dim>())) );
        backendElastic->attachNearNullSpace( myNullSpace );
        
        auto Res = backendElastic->newVector(Xv_);
        auto Jac = backendElastic->newMatrix( _test=Xv_, _trial=Xv_ );

        auto ResTrue = backendRigid->newVector(Xv_);
        auto JacTrue = backendRigid->newMatrix( _test=Xv_, _trial=Xv_ );

        std::ofstream ofs("res.csv");
        ofs << fmt::format("time, errorL2, errorH1") << std::endl;

        // Translation
        auto a_translation_ = form2( _test = XvR_, _trial = XvR_ );
        auto l_translation_ = form1( _test = XvR_ );
        auto lt_translation_ = form1( _test = XvR_ );
    
        a_translation_.zero();
        l_translation_.zero();
        lt_translation_.zero();

        l_translation_ = integrate( _range = elements(support(XvR_)), _expr = cst(density_)*trans(expr<Dim,1>( externalforce_ ))*id(u_trans_true));
        a_translation_ = integrate( _range = elements(support(XvR_)), _expr = cst(density_)*inner( ts_translation->polyDerivCoefficient()*idt(u_trans_true),id( u_trans_true ) ) );
        
        while (time_step_ * iter < final_time_)
        {
            if (Environment::isMasterRank())
                std::cout << "Time step : " << iter*time_step_ << std::endl;
            
            // Solve contact
            //u_total = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = idv(u_rigid) + idv(u_elastic));
            myelts_ = getContactRegion(u_total);
            //myelts_ = getContactRegion(u_true);
            std::cout << "Faces in contact : " << nbrFaces_ << std::endl;

            
            /*
                Solve translation 
            */
            std::cout << "Solve translation" << std::endl;
            std::cout << "ForceContact : " << ForceContact[0] << ", " <<  ForceContact[1] << std::endl;
            auto u_trans_iter = time_step_ * dt_u_trans + time_step_*time_step_*(1.-2*beta_)/2. * dtt_u_trans + beta_*time_step_*time_step_ * (Force + ForceContact);
            u_trans += u_trans_iter;

            // Update Newmark scheme
            dtt_u_trans_old = dtt_u_trans;
            dtt_u_trans = 1. / (beta_ *std::pow(time_step_,2)) * u_trans_iter - 1. / (beta_ * time_step_) *  dt_u_trans -  (1./(2.*beta_) - 1.) *  dtt_u_trans;
            dt_u_trans = dt_u_trans + time_step_ * ((1. - gamma_) * dtt_u_trans_old + gamma_ * dtt_u_trans) ;
            
            //std::cout << "utrans : " << u_trans << std::endl;
            //std::cout << "dtt_u_trans : " << dtt_u_trans << std::endl;
            //std::cout << "dt_u_trans : " << dt_u_trans << std::endl;
            
            lt_translation_.zero();
            lt_translation_ = l_translation_;
            lt_translation_ += integrate( _range = elements(support(XvR_)), _expr = cst( density_ ) * inner( idv( ts_translation->polyDeriv() ), id( u_trans_true ) ) );
            if (nbrFaces_ > 0)//Add contact terms
            {
                auto F = Id + gradv(u_elastic);
                auto epsv = sym(gradv(u_elastic)) + 0.5*trans(gradv(u_elastic))*gradv(u_elastic);
                auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv);
                
                //lt_translation_ += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idv(u_elastic),  id(u_trans_true))) ;
                //lt_translation_ += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * (trans(expr<Dim,1>(direction_))*idv(u_elastic)) * id(u_trans_true)) ;
                //lt_translation_ += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * (trans(expr<Dim,1>(direction_))*idv(u_elastic)) * id(u_trans_true)) ;
                
                
	            lt_translation_ += integrate( _range = myelts_, _expr = cst(1)/det(F) * inner(val(F*sigmav*trans(F))*N(), id( u_trans_true ) ) );
            }
        
            a_translation_.solve( _rhs = lt_translation_, _solution = u_trans_true, _name="translation" );
            ts_translation->updateFromDisp(u_trans_true);
            ts_translation->next(u_trans_true);

            /*
                Solve rotation
            */
            auto massCenter = mean( _range = elements(support(Xv_)), _expr = P());
            auto massCenterVec = vec(cst(massCenter(0,0)),cst(massCenter(1,0)));
            auto momentOfInertia = integrate(_range=elements(support(Xv_)),_expr=cst(density_)*( (Px()-massCenter(0,0))*(Px()-massCenter(0,0)) + (Py()-massCenter(1,0))*(Py()-massCenter(1,0)) ) ).evaluate()(0,0);

            double T = 0.;
            T += integrate(_range= markedfaces(mesh_, "Upper"), _expr= -(Py()-massCenter(1,0))*cst(Upper_x)).evaluate()(0,0);
            T += integrate(_range= markedfaces(mesh_, "Lower"), _expr= -(Py()-massCenter(1,0))*cst(Lower_x)).evaluate()(0,0);

            //std::cout << "T : " << T << std::endl;
            acceleration = (T/momentOfInertia);
            omega = omega + time_step_ * (T/momentOfInertia);
            theta += time_step_ * omega;

            auto theta_iter = time_step_ * omega;

            auto rot = vec(
                cos(theta) * (Px() - massCenter(0,0)) - sin(theta) * (Py() - massCenter(1,0)) - Px() + massCenter(0,0),
                sin(theta) * (Px() - massCenter(0,0)) + cos(theta) * (Py() - massCenter(1,0)) - Py() + massCenter(1,0)
            );

            auto rot_iter = vec (
                cos(theta_iter) * (Px() - massCenter(0,0)) - sin(theta_iter) * (Py() - massCenter(1,0)) - Px() + massCenter(0,0),
                sin(theta_iter) * (Px() - massCenter(0,0)) + cos(theta_iter) * (Py() - massCenter(1,0)) - Py() + massCenter(1,0)
            );
    
            u_theta = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = rot);
            auto u_theta_iter = project(_space = Xv_, _range = elements(support(Xv_)), _expr = rot_iter);
            dtt_u_theta_old = project(_space = Xv_, _range =  elements(support(Xv_)), _expr= idv(dtt_u_theta));   
            dtt_u_theta = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = cst(1.0) / (cst(beta_)*std::pow(time_step_,2)) * idv(u_theta_iter) - cst(1.0) / (cst(beta_)*cst(time_step_)) * idv( dt_u_theta ) -  (cst(1.0)/(cst(2.0)*cst(beta_)) - cst(1.0)) * idv( dtt_u_theta ));
            dt_u_theta = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = idv(dt_u_theta) + cst(1.0)*cst(time_step_) * ((cst(1.0) - cst(gamma_)) * idv(dtt_u_theta_old) + cst(gamma_) * idv(dtt_u_theta)) ); 
            
            auto coef_1 = (- cst(acceleration) * sin(theta) - cst(omega)*cst(omega) * cos(theta));
            auto coef_2 = (- cst(acceleration) * cos(theta) + cst(omega)*cst(omega) * sin(theta));
            auto coef_3 = (cst(acceleration) * cos(theta) - cst(omega)*cst(omega) * sin(theta));
            auto coef_4 = (- cst(acceleration) * sin(theta) - cst(omega)*cst(omega) * cos(theta));
            rot_a = project(_space = Xv_, _range =  elements(support(Xv_)), _expr =  vec(coef_1*(Px() - massCenter(0,0)) + coef_2*(Py() - massCenter(1,0)), coef_3*(Px() - massCenter(0,0)) + coef_4*(Py() - massCenter(1,0))));
            
            u_rigid.on(_range=elements(support(Xv_)), _expr=  idv(u_theta) + vec(cst(u_trans[0]),cst(u_trans[1])));
            //u_rigid.on(_range=elements(support(Xv_)), _expr=  idv(u_theta) + idv(u_trans_true));
            u_total = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = idv(u_rigid) + idv(u_elastic));

            /*
                Solve elasticity
            */

            auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
            {
                auto u = Xv_->element();
                u = *X;

                //auto Fv = Id + gradv(u) + gradv(u_rigid);
                //auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*(gradv(u));
                //auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

                //auto dF = gradt(u);
                //auto dE = sym(gradt(u)) + 0.5*trans(gradt(u))*(gradv(u)) + 0.5*trans(gradv(u))*(gradt(u));
                //auto dS = lambda_*trace(dE)*Id + 2*mu_*dE;

                //auto Emixte = 0.5*trans(gradv(u))*gradv(u_rigid) + 0.5*trans(gradv(u_rigid))*gradv(u);
                //auto Smixte = lambda_*trace(Emixte)*Id + 2*mu_*Emixte;
                //auto dEmixte = 0.5*trans(gradt(u))*gradv(u_rigid) + 0.5*trans(gradv(u_rigid))*gradt(u);
                //auto dSmixte = lambda_*trace(dEmixte)*Id + 2*mu_*dEmixte;

                //auto Erigid = 0.5*trans(gradv(u_rigid))*gradv(u_rigid);
                //auto Srigid = lambda_*trace(Erigid)*Id + 2*mu_*Erigid;
            
            
                auto a = form2( _test=Xv_, _trial=Xv_, _matrix=J );

                auto Fv = Id + gradv(u);
                auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
                auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

                auto dF = gradt(u);
                auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
                auto dS = lambda_*trace(dE)*Id + 2*mu_*dE;

   
                a = integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner( ts_true->polyDerivCoefficient()*idt(u),id( u ) ) );
                a += integrate( _range=elements(support(Xv_)), _expr = inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );

   
                //a = integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );
                //a += integrate( _range=elements(support(Xv_)), _expr = inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
                //a += integrate( _range=elements(support(Xv_)), _expr = inner( dF*val(Smixte) + val(Fv)*dSmixte , grad(u) ) );
                //a += integrate( _range=elements(support(Xv_)), _expr =  inner( dF*val(Srigid),grad(u) ) );

                //if (nbrFaces_ > 0)
                //    a += integrate(_range=myelts_, _expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)));

                if (nbrFaces_ > 0)
                {
                    a += integrate (_range=myelts_,_expr= cst(gamma_contact) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) );
                    a += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*dF*val(Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)));
                    a += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv)*dS*N(),trans(expr<Dim,1>(direction_))*id(u)));
                }
            
            };
       
            auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
            {
                auto u = Xv_->element();
                u = *X;
               
                //auto Fv = Id + gradv(u) + gradv(u_rigid);
                //auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*(gradv(u));
                //auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

                //auto Emixte = 0.5*trans(gradv(u))*gradv(u_rigid) + 0.5*trans(gradv(u_rigid))*gradv(u);
                //auto Smixte = lambda_*trace(Emixte)*Id + 2*mu_*Emixte;

                //auto Erigid = 0.5*trans(gradv(u_rigid))*gradv(u_rigid);
                //auto Srigid = lambda_*trace(Erigid)*Id + 2*mu_*Erigid;

                auto Fv = Id + gradv(u);
                auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
                auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;
   
                auto r = form1( _test=Xv_, _vector=R );

                r = integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
                //r += integrate( _range=elements(support(Xv_)), _expr = inner( val(Fv*Sv) , grad(u) ) );
                //r += integrate( _range=elements(support(Xv_)), _expr = inner( val(Fv*Smixte) , grad(u) ) );
                //r += integrate( _range=elements(support(Xv_)), _expr = inner( val(Fv*Srigid) , grad(u) ) );
                //r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner(idv(dtt_u_theta), id(u)) ) ;
               
               
                r += integrate( _range=elements(support(Xv_)), _expr = inner( val(Fv*Sv) , grad(u) ) );
                r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)* trans(vec(-cst(acceleration)*(Py() - massCenter(1,0))-cst(omega)*cst(omega)*(Px() - massCenter(0,0)), cst(acceleration)*(Px() - massCenter(0,0)) -cst(omega)*cst(omega)*(Py() - massCenter(1,0)) ))*id(u) ) ;
                //r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)* inner(idv(dtt_u_theta),id(u)) ) ;

                r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)*(trans(vec(cst(dtt_u_trans[0]),cst(dtt_u_trans[1]))) - trans(expr<Dim,1>( externalforce_ )))*id( u ) ) ;
                
                //r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner(idv(ts_translation->currentAcceleration()),id( u ) ));
                //r += integrate( _range=elements(support(Xv_)), _expr= cst(density_)*(inner(idv(ts_translation->currentAcceleration()),id( u ) ) - trans( expr<Dim, 1>( externalforce_ ) )*id( u ))  );
                r += integrate( _range = markedfaces(mesh_, "Upper"), _expr = - trans(expr<Dim,1>( neumannUpper_ ))*id(u));
                r += integrate( _range = markedfaces(mesh_, "Lower"), _expr = - trans(expr<Dim,1>( neumannLower_ ))*id(u));

                //if (nbrFaces_ > 0)
                //{
                //    r += integrate (_range=myelts_,_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) );
                //    r += integrate (_range=myelts_,_expr= - cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_rigid),trans(expr<Dim,1>(direction_))*id(u)) );
                //}   

                if (nbrFaces_ > 0)
                {
                    r += integrate (_range=myelts_,_expr= cst(gamma_contact) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) );
                    r += integrate (_range=myelts_,_expr= - cst(gamma_contact) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_rigid),trans(expr<Dim,1>(direction_))*id(u)) );
                    r += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv*Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)) );
                } 
            
                R->close();
            };
            std::cout << "solve elasticity" << std::endl;
            backendElastic->nlSolver()->residual = Residual;
            backendElastic->nlSolver()->jacobian = Jacobian;
            backendElastic->nlSolve( _solution=u_elastic,_jacobian=Jac,_residual=Res );
   
            ts_->updateFromDisp(u_elastic);
            ts_->next(u_elastic);

            // Update contact force 
            if (nbrFaces_ > 0)//Add contact terms
            {
                auto F = Id + gradv(u_elastic);
                auto epsv = sym(gradv(u_elastic)) + 0.5*trans(gradv(u_elastic))*gradv(u_elastic);
                auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
                auto J = det(F);
                
                auto trac = Xv_->element();
                trac.on(_range=elements(support(Xv_)), _expr= expr<Dim,1>("{0.,0.}"));
                trac.on(_range=myelts_, _expr = cst(density_) * val(F*sigmav));

                auto force_C = integrate( _range = elements(support(Xv_)), _expr = divv(trac)/det(F) ).evaluate();
                //auto force_C = integrate( _range = myelts_, _expr = trans(expr<Dim,1>(direction_))*val(F*sigmav) ).evaluate();
                //auto force_C = integrate( _range = myelts_, _expr = - cst(1.)/cst(epsilon_) * (trans(expr<Dim,1>(direction_))*idv(u_total) - idv(g_) )).evaluate();
                //ForceContact[0] = ddirection_[0]*force_C(0,0);
                //ForceContact[1] = ddirection_[1]*force_C(0,0);
                std::cout << "force_C : " << force_C(0,0) << std::endl;
                
                
                auto volume = integrate( _range = elements(support(X_)), _expr = J ).evaluate()(0,0);
                mass_ = volume*density_;
                
                std::cout << "mass _ " << mass_ << std::endl;
                ForceContact[0] = ddirection_[0]*force_C(0,0)/mass_;
                ForceContact[1] = ddirection_[1]*force_C(0,0)/mass_;

            } 
            else  // on remet la force à 0
            {
                ForceContact[0] = 0;
                ForceContact[1] = 0;
            }
            
            u_total = project(_space = Xv_, _range =  elements(support(Xv_)), _expr = idv(u_rigid) + idv(u_elastic));

            //myelts_ = getContactRegion(u_total);
            //std::cout << "Faces in contact : " << nbrFaces_ << std::endl;

            // Checks
            auto meanDisp = mean(_range=elements(support(Xv_)), _expr=idv(u_elastic));
            std::cout << "Mean displacement x : " << meanDisp(0,0) << std::endl;
            std::cout << "Mean displacement y : " << meanDisp(1,0) << std::endl;
            auto meanCurl = mean(_range=elements(support(Xv_)), _expr = curlv_op(u_elastic));
            std::cout << "Mean curl x : " << meanCurl(0,0) << std::endl;
            std::cout << "Mean curl y : " << meanCurl(1,0) << std::endl;
            std::cout << "Mean curl z : " << meanCurl(2,0) << std::endl;

            // True solution
            
            auto JacobianTrue = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
            {
                auto u = Xv_->element();
                u = *X;

                auto Fv = Id + gradv(u);
                auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
                auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

                auto dF = gradt(u);
                auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
                auto dS = lambda_*trace(dE)*Id + 2*mu_*dE;

                auto a = form2( _test=Xv_, _trial=Xv_, _matrix=J );
   
                a = integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner( ts_true->polyDerivCoefficient()*idt(u),id( u ) ) );
                a += integrate( _range=elements(support(Xv_)), _expr = inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );

                if (nbrFaces_ > 0)
                {
                    a += integrate (_range=myelts_,_expr= cst(gamma_contact) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) );
                    a += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*dF*val(Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)));
                    a += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv)*dS*N(),trans(expr<Dim,1>(direction_))*id(u)));
                }
                

            };
       
            auto ResidualTrue = [=](const vector_ptrtype& X, vector_ptrtype& R)
            {
               auto u = Xv_->element();
               u = *X;
               
               auto Fv = Id + gradv(u);
               auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
               auto Sv = lambda_*trace(Ev)*Id + 2*mu_*Ev;

               auto r = form1( _test=Xv_, _vector=R );

               r = integrate( _range=elements(support(Xv_)), _expr = cst(density_)*inner( ts_true->polyDerivCoefficient()*idv(u) -idv(ts_true->polyDeriv()),id( u ) ) );
               r += integrate( _range=elements(support(Xv_)), _expr = inner( val(Fv*Sv) , grad(u) ) );
               r += integrate( _range = markedfaces(mesh_, "Upper"), _expr = - trans(expr<Dim,1>( neumannUpper_ ))*id(u));
               r += integrate( _range = markedfaces(mesh_, "Lower"), _expr = - trans(expr<Dim,1>( neumannLower_ ))*id(u));
               r += integrate( _range=elements(support(Xv_)), _expr = cst(density_)*(- trans(expr<Dim,1>( externalforce_ )))*id( u ) ) ;

                if (nbrFaces_ > 0)
                {
                    r += integrate (_range=myelts_,_expr= cst(gamma_contact) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) );
                    r += integrate (_range=myelts_,_expr= - cst(gamma_contact) * inner(idv(g_),trans(expr<Dim,1>(direction_))*id(u)) );
                    r += integrate (_range=myelts_,_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv*Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)) );
                }   
            
               R->close();
            };
            std::cout << "solve elasticity true" << std::endl;
            backendRigid->nlSolver()->residual = ResidualTrue;
            backendRigid->nlSolver()->jacobian = JacobianTrue;
            backendRigid->nlSolve( _solution=u_true,_jacobian=JacTrue,_residual=ResTrue );
   
            ts_true->updateFromDisp(u_true);
            ts_true->next(u_true);
            
            
            // Compute errors
            auto l2err_u = normL2( _range=elements(support(Xv_)), _expr=idv(u_true) - idv(u_total) );
            auto h1err_u = normH1( _range=elements(support(Xv_)), _expr=idv(u_true) - idv(u_total) , _grad_expr=gradv(u_true) - gradv(u_total) );
            //ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}, {:.6f}", time_step_*iter, l2err_u, h1err_u) << std::endl;

            if ( Environment::isMasterRank() )
            {
                std::cout << " l2err_u = " << l2err_u << std::endl;
                std::cout << " h1err_u = " << h1err_u << std::endl;
            }
            
            // Export 
            this->exportResults(iter*time_step_);
            iter++;
        }   
        ofs.close();
    }
    
}
}