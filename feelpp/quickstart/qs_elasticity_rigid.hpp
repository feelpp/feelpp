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

namespace Feel
{


template <int Dim, int Order>
class ElasticRigid
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;
    using spacerv_t = Pchv_type<mesh_t, 0>;
    using spacev_t = Pchv_type<mesh_t, Order>;
    using space_t = Pch_type<mesh_t, Order>;
    using spacev_ptr_t = Pchv_ptrtype<mesh_t, Order>; 
    using space_ptr_t = Pch_ptrtype<mesh_t, Order>;
    using spacerv_ptr_t = Pchv_ptrtype<mesh_t, 0>;
    using elementv_t = typename spacev_t::element_type;
    using elementrv_t = typename spacerv_t::element_type;
    using element_t = typename space_t::element_type;
    using form2_type = form2_t<spacev_t,spacev_t>; 
    using form1_type = form1_t<spacev_t>; 
    using ts_ptrtype = std::shared_ptr<NewmarkContact<spacev_t>>;
    using tsr_ptrtype = std::shared_ptr<NewmarkContact<spacerv_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

    // Constructors
    ElasticRigid() = default;
    ElasticRigid(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    
    elementv_t const& u() const { return u_; }
    elementv_t const& u_e() const { return u_e_; }
    elementv_t const& u_r() const { return u_r_; }


    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }
    void setU_e(elementv_t const& u_e) { u_e_ = u_e; }
    void setU_r(elementv_t const& u_r) { u_r_ = u_r; }

    // Methods
    void initialize();
    void initializeContact();
    Range<mesh_t, MESH_FACES> getContactRegion(elementv_t const& u);
    void run();
    void timeLoop();
    void exportResults(double t);
    void initG();

private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t Xhv_;
    space_ptr_t Xh_;
    spacerv_ptr_t VhvC_;

    elementv_t u_;
    elementv_t u_e_;
    elementrv_t u_r_;

    element_t contactRegion_;
    element_t contactPressure_;
    element_t contactDisplacement_;
    element_t g_;
    int nbrFaces_;
    Range<mesh_t, MESH_FACES> myelts_;

    ts_ptrtype ts_e_;
    tsr_ptrtype ts_r_;
    exporter_ptrtype e_;
    nl::json meas_;

    double H_,E_, nu_, lambda_, mu_, rho_,m_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    int withMarker_,save_, computeerror_;
};


// Constructor
template <int Dim, int Order>
ElasticRigid<Dim, Order>::ElasticRigid(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // Define Xhv
    Xhv_ = Pchv<Order>(mesh_); 
    VhvC_ = Pchv<0>(mesh_);

    // Get elastic structure parameters
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    m_ = integrate( _range = elements(mesh_), _expr = cst(rho_) ).evaluate()(0,0);

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
    u_e_ = Xhv_->element();
    u_r_ = VhvC_->element();

    auto u0_ = Xhv_->element();
    auto ur0_ =  VhvC_->element();

    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    auto init_displ = expr<Dim,1>(get_value(specs_, "/InitialConditions/LinearElasticity/displacement/expr", default_displ ));
    u0_.on(_range=elements(mesh_), _expr=init_displ);
    
    ts_e_ = newmarkContact(Xhv_, steady, initial_time, final_time, time_step, gamma, beta );
    ts_e_->start();
    ts_e_->initialize( u0_ );

    ts_r_ = newmarkContact(VhvC_, steady, initial_time, final_time, time_step, gamma, beta );
    ts_r_->start();
    ts_r_->initialize( ur0_ );

    u_e_ = u0_;
    u_r_ = ur0_;
    
    ts_e_->updateFromDisp(u_e_);
    ts_r_->updateFromDisp(u_r_);
}


template <int Dim, int Order>
void 
ElasticRigid<Dim, Order>::initializeContact()
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

    //std::string matFixedPointtol = fmt::format( "/Collision/LinearElasticity/fixedPointTol" );
    //fixedPointtol_ = specs_[nl::json::json_pointer( matFixedPointtol )].get<double>(); 
    
    //std::string matFixedPoint = fmt::format( "/Collision/LinearElasticity/fixedPoint" );
    //fixedPoint_ = specs_[nl::json::json_pointer( matFixedPoint )].get<int>(); 

    //std::string matpressurePoint = fmt::format( "/Collision/LinearElasticity/pressurePoint" );
    //pressurePoint_ = specs_[nl::json::json_pointer( matpressurePoint )].get<std::vector<double>>();  
}


template <int Dim, int Order>
Range<typename ElasticRigid<Dim, Order>::mesh_t, MESH_FACES>
ElasticRigid<Dim, Order>::getContactRegion(elementv_t const& u)
{   
    Range<mesh_t,MESH_FACES> myelts(mesh_ );

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


template <int Dim, int Order>
void 
ElasticRigid<Dim, Order>::initG()
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

// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::timeLoop()
{
    /*
    
    Elasticity equations

    */

    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_ = form1( _test = Xhv_ );
    auto lt_ = form1( _test = Xhv_ );
    
    a_.zero();
    at_.zero();
    l_.zero();
    lt_.zero();

    //l_ += integrate( _range = elements(mesh_), _expr = cst(0.));

    auto deft = sym(gradt(u_e_));
    auto def = sym(grad(u_e_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;
    a_ += integrate( _range = elements(mesh_), _expr = cst(rho_)*inner( ts_e_->polyDerivCoefficient()*idt(u_e_),id( u_e_ ) ) + inner(sigmat,def));


    /*
    
    Rigid equations

    */


    auto a_r_ = form2( _test = VhvC_, _trial = VhvC_ );
    auto at_r_ = form2( _test = VhvC_, _trial = VhvC_ );
    auto l_r_ = form1( _test = VhvC_ );
    auto lt_r_ = form1( _test = VhvC_ );
    
    a_r_.zero();
    at_r_.zero();
    l_r_.zero();
    lt_r_.zero();

 
    l_r_ += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_r_));
    a_r_ += integrate( _range = elements(mesh_), _expr = cst(rho_)*inner( ts_r_->polyDerivCoefficient()*idt(u_r_),id( u_r_ ) ) );


    for ( ts_e_->start(); ts_e_->isFinished()==false; ts_e_->next(u_e_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_e_->time(),ts_e_->timeFinal()) << std::endl;

        

        std::cout << "***** Compute new contact region *****" << std::endl;

        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=myelts_, _expr = cst(1.));

        std::cout << "***** Solve rigid *****" << std::endl;

        lt_r_ = l_r_;
        at_r_ = a_r_;

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_r_ +=  integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_r_->polyDeriv()),id( u_r_ ) ));

        
        if (nbrFaces_ > 0)
        {
            auto epsv = sym(gradv(u_));
            auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
            //std::cout << " Contact pressure : " << sigmav*vec(cst(0.),cst(-1.)) << std::endl;
    
            lt_r_ += integrate( _range=boundaryfaces(mesh_), _expr= inner( vec(cst(0.), -trans(vec(cst(0.),cst(-1.)))*sigmav*vec(cst(0.),cst(-1.)) ),id( u_r_ )) * idv(contactFaces));
        }

        at_r_.solve( _rhs = lt_r_, _solution = u_r_, _rebuild = true );
        ts_r_->updateFromDisp(u_r_);
        
        // Reset
        at_r_.zero();
        lt_r_.zero();

        ts_r_->next(u_r_);


        std::cout << "***** Solve elasticity *****" << std::endl;
        
        lt_ = l_;
        at_ = a_;

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
            lt_ +=  integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_e_->polyDeriv()),id( u_e_ ) ));
        
        if (nbrFaces_ > 0)
        {
            auto epsv = sym(gradv(u_));
            auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
            
            at_ += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_),trans(expr<Dim,1>(direction_))*id(u_e_)) * idv(contactFaces));
            lt_ += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_r_), trans(expr<Dim,1>(direction_))*id(u_e_)) * idv(contactFaces));     
            
            //lt_ += integrate( _range=boundaryfaces(mesh_), _expr= inner( vec(cst(0.), trans(vec(cst(0.),cst(-1.)))*sigmav*vec(cst(0.),cst(-1.)) ),id( u_r_ )) * idv(contactFaces));

        }
        
        at_.solve( _rhs = lt_, _solution = u_e_ , _rebuild = true);

        ts_e_->updateFromDisp(u_e_);

        // Reset
        at_.zero();
        lt_.zero();

        std::cout << "***** Solve displacement *****" << std::endl;
        u_.on(_range=elements(mesh_),_expr=idv(u_e_) + idv(u_r_));

        // Export
        this->exportResults(ts_e_->time());

    }    

}

// Run method 
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::run()
{
    std::cout << "***** Initialize parameters *****" << std::endl;
    initialize();

    std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    std::cout <<  "***** Initialize distance g *****" << std::endl;
    initG();

    this->exportResults(0.);

    std::cout <<  "***** Start time loop *****" << std::endl;
    timeLoop();
                
}

// Export results
template <int Dim, int Order>
void 
ElasticRigid<Dim, Order>::exportResults(double t)
{

    e_->step(t)->addRegions();

    e_->step(t)->add( "displacement", u_ );
    e_->step(t)->add( "displacement_elastic", u_e_ );
    e_->step(t)->add( "displacement_rigid", u_r_ );

    //e_->step(*it)->add( "velocity_elastic", ts_e_->currentVelocity() );
    //e_->step(t)->add( "velocity_rigid", ts_r_->currentVelocity() );

    e_->save();
}

} 