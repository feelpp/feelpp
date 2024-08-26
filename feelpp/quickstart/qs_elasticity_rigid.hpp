#pragma once

#include "qs_elasticity_contact.hpp"
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
    elementrv_t const& u_r() const { return u_r_; }


    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }
    void setU_e(elementv_t const& u_e) { u_e_ = u_e; }
    void setU_r(elementrv_t const& u_r) { u_r_ = u_r; }

    // Methods
    void initialize();
    void initializeContact();
    Range<mesh_t, MESH_FACES> getContactRegion(elementv_t const& u);
    void run();
    void timeLoop();
    void timeLoopFixedPoint();
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

    double H_,E_, nu_, lambda_, mu_, rho_, m_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    int withMarker_,save_, computeerror_;
    std::vector<double> pressurePoint_;
    int fixedPoint_;
    double fixedPointtol_;
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
    // Mesh and spaces
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

    Xhv_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );
    VhvC_ = Pchv<0>(mesh_, markedelements(mesh_, "Caoutchouc"));

    // Solid parameters
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    m_ = integrate( _range = elements(support(Xhv_)), _expr = cst(rho_) ).evaluate()(0,0);

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
            }
        }
    }

    // Exporter
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
    u0_.on(_range=elements(support(Xhv_)), _expr=init_displ);

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
    Xh_ = Pch<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );

    // Initialize contact field
    contactRegion_ =  project(_space=Xh_, _range=elements(support(Xhv_)), _expr = cst(0.));
    contactPressure_ =  project(_space=Xh_, _range=elements(support(Xhv_)), _expr = cst(0.));
    contactDisplacement_ = project(_space=Xh_, _range=elements(support(Xhv_)), _expr = cst(0.));
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
ElasticRigid<Dim, Order>::getContactRegion(elementv_t const& u)
{
    Range<mesh_t,MESH_FACES> myelts(mesh_ );

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




template <int Dim, int Order>
void
ElasticRigid<Dim, Order>::initG()
{
    // Init the distance fields
    g_ = Xh_->element();

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
            //allrays.push_back(ray);
#if 1            
            auto rayIntersection = bvh->intersect( _ray = ray );
            if (!rayIntersection.empty())
            {
                for ( auto const& rir : rayIntersection )
                {
                    for (auto const& ldof  : Xh_->dof()->faceLocalDof( face.id() ))
                        g_[ldof.index()] = rir.distance() - tolDistance_;
                }
            }
#endif            
        }
        
    }
#if 0    
    // compute intersections for allrays
    auto multiRayIntersectionResult = bvh->intersect(_ray=allrays);//,_parallel=false);
    for (auto const& rir : multiRayIntersectionResult)
    {
        for (auto const& ldof : Xh_->dof()->faceLocalDof(rir.face().id()))
            g_[ldof.index()] = rir.distance() - tolDistance_;
    }
#endif    

}

// Time loop
template <int Dim, int Order>
void ElasticRigid<Dim, Order>::timeLoop()
{
    
    // Elasticity equations

    auto a_e = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_e = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_e = form1( _test = Xhv_ );
    auto lt_e = form1( _test = Xhv_ );

    a_e.zero();
    at_e.zero();
    l_e.zero();
    lt_e.zero();

   
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
    auto lt_r_f_ = form1( _test = VhvC_ );

    a_r_.zero();
    at_r_.zero();
    l_r_.zero();
    lt_r_.zero();
    lt_r_f_.zero();

    auto dir = oneZ(); //-expr<Dim,1>(direction_);

    l_r_ += integrate( _range = elements(support(VhvC_)), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_r_));
    a_r_ += integrate( _range = elements(support(VhvC_)), _expr = cst(rho_)*inner( ts_r_->polyDerivCoefficient()*idt(u_r_),id( u_r_ ) ) );

    // output 
    std::ofstream ofs("outputs.csv");
    ofs << fmt::format( "time, f_s, f_C, f_t, mean_displ_ex, mean_displ_ey, mean_displ_rx, mean_displ_ry, mean_disp_x, mean_displ_y, evaluateStress, evaluateDisp") << std::endl;


    // Time iterations
    int iteration = 1;
    for ( ts_e_->start(); ts_e_->isFinished() == false; ts_e_->next( u_e_ ), ts_r_->next( u_r_ ) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.3f}/{}", fmt::localtime(std::time(nullptr)), ts_e_->time(),ts_e_->timeFinal()) << std::endl;

        std::cout << "***** Compute new contact region *****" << std::endl;

        myelts_ = getContactRegion(u_);
        std::cout << "At iteration : " << iteration << " Nbr faces for processContact : " << nbrFaces_ << std::endl;

        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xh_) ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=myelts_, _expr = cst(1.));

        
        
        std::cout << "***** Solve rigid *****" << std::endl;
        lt_r_f_ = l_r_;
        at_r_ = a_r_;

        lt_r_ += integrate( _range = markedelements( support( VhvC_ ), "Caoutchouc" ), _expr = cst( rho_ ) * inner( idv( ts_r_->polyDeriv() ), id( u_r_ ) ) );

        if (nbrFaces_ > 0)
        {
            auto epsv = sym(gradv(u_e_));
            auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
            
            lt_r_f_ +=  integrate( _range=boundaryfaces(support(VhvC_)), _expr= inner( sigmav,id( u_r_ )) * idv(contactFaces));            
        }

        lt_r_ += lt_r_f_;
        at_r_.solve( _rhs = lt_r_, _solution = u_r_, _rebuild = true );
        ts_r_->updateFromDisp(u_r_);

        std::cout << "***** Solve elasticity *****" << std::endl;

        lt_e = l_e;
        at_e = a_e;


        lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= cst(rho_)*inner( idv(ts_e_->polyDeriv()),id( u_e_ ) ));
        lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= -cst(rho_)*inner( idv(ts_r_->currentAcceleration()),id( u_e_ ) ));
        
        if (nbrFaces_ > 0)
        {
            
            if (method_.compare("penalty") == 0)
            {
                at_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_),trans(expr<Dim,1>(direction_))*id(u_e_)) * idv(contactFaces));
                lt_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_r_), trans(expr<Dim,1>(direction_))*id(u_e_)) * idv(contactFaces));
            }

            else if (method_.compare("nitsche") == 0)
            {
                auto const Id = eye<Dim,Dim>();
                auto deft = sym(gradt(u_e_));
                auto def = sym(grad(u_e_));
                auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
                auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

                at_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces)); 
                at_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u_e_) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u_e_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces));
                lt_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_r_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u_e_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces));     
            }
        }

        // deleting translation
        at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_e_, _expr=0.*one());

        at_e.solve( _rhs = lt_e, _solution = u_e_ , _rebuild = true);
        ts_e_->updateFromDisp(u_e_);

        std::cout << "***** Solve displacement *****" << std::endl;
        u_.on(_range=elements(support(Xhv_)),_expr=idv(u_e_) + idv(u_r_));

        // Exports
        auto vr = VhvC_->element();
        vr.on(_range=elements(support(VhvC_)),_expr=dir);
        
        auto f_s = l_r_(vr);
        auto f_t = lt_r_f_(vr);
        auto epsv = sym(gradv(u_e_));
        auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();
        auto fc_density_expr = -trans(expr<Dim,1>(direction_))*sigmav;
        auto f_C = integrate(_range=boundaryfaces(support(Xh_)), _expr = fc_density_expr * idv(contactFaces)).evaluate();

        std::cout << fmt::format("forces : gravity: [{}] contact [{}] sum: [{}]", f_s, f_C, f_t) << std::endl;
        
        double meas = integrate(_range=elements(support(Xhv_)),_expr=cst(1.)).evaluate()(0,0);
        auto dist2wall_r = integrate(_range=elements(support(VhvC_)),_expr=(idv(g_)-trans(expr<Dim,1>(direction_))*idv(u_r_))).evaluate()(0,0)/meas;
        auto dist2wall_ = integrate(_range=elements(support(Xhv_)),_expr=(idv(g_)-trans(expr<Dim,1>(direction_))*idv(u_))).evaluate()(0,0)/meas;

        std::cout << fmt::format("distance to wall, rigid: [{}] total: [{}]", dist2wall_r, dist2wall_) << std::endl;
       
        auto mean_displ_e = integrate(_range=elements(support(Xhv_)),_expr=idv(u_e_)).evaluate()/meas;
        auto mean_displ_r = integrate(_range=elements(support(VhvC_)),_expr=idv(u_r_)).evaluate()/meas;
        auto mean_displ = integrate(_range=elements(support(Xhv_)),_expr=idv(u_)).evaluate()/meas;

        std::cout << fmt::format("mean displacement, elastic: [{}] rigid: [{}] total: [{}]", mean_displ_e, mean_displ_r, mean_displ) << std::endl;

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

        auto contactPressure = Xh_->element();
        contactPressure.on(_range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
        auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure) ); 

        auto evaluateDispExpr = Xh_->element();
        evaluateDispExpr.on(_range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
    
        if ( Dim == 2 )
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}",
                                ts_e_->time(), f_s, f_C( 0, 0 ), f_t,
                                mean_displ_e( 0, 0 ), mean_displ_e( 1, 0 ), mean_displ_r( 0, 0 ), mean_displ_r( 1, 0 ), mean_displ( 0, 0 ), mean_displ( 1, 0 ), evaluateStress(0,0), evaluateDisp(0,0))
                << std::endl;
        else
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}, {:.6f}",
                                ts_e_->time(), f_s, f_C( 0, 0 ), f_t,
                                mean_displ_e( 0, 0 ), mean_displ_e( 1, 0 ), mean_displ_e( 2, 0 ),
                                mean_displ_r( 0, 0 ), mean_displ_r( 1, 0 ), mean_displ_r( 2, 0 ), 
                                mean_displ( 0, 0 ), mean_displ( 1, 0 ), mean_displ( 2, 0 ), evaluateStress(0,0), evaluateDisp(0,0) )
                << std::endl;

        // Export
        this->exportResults( ts_e_->time() );

        // Reset
        at_r_.zero();
        lt_r_.zero();
        lt_r_f_.zero();
        at_e.zero();
        lt_e.zero();

        iteration++;

    }
    ofs.close();
}



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

        std::cout << "***** Compute new contact region *****" << std::endl;
        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;

        if (nbrFaces_ == 0)
        {
            std::cout << "At iteration : " << iteration << " no contact" << std::endl;

            std::cout << "***** Solve rigid *****" << std::endl;
            lt_r_ = l_r_;
            at_r_ = a_r_;

            lt_r_ += integrate( _range = markedelements( support( VhvC_ ), "Caoutchouc" ), _expr = cst( rho_ ) * inner( idv( ts_r_->polyDeriv() ), id( u_r_ ) ) );

            at_r_.solve( _rhs = lt_r_, _solution = u_r_, _rebuild = true );
            ts_r_->updateFromDisp(u_r_);

            std::cout << "***** Solve elasticity *****" << std::endl;
            lt_e = l_e;
            at_e = a_e;

            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= cst(rho_)*inner( idv(ts_e_->polyDeriv()),id( u_e_ ) ));
            lt_e += integrate( _range=markedelements(support(Xh_),"Caoutchouc"), _expr= -cst(rho_)*inner( idv(ts_r_->currentAcceleration()),id( u_e_ ) ));
            
            at_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_e, _element=u_e_, _expr=0.*one());
            at_e.solve( _rhs = lt_e, _solution = u_e_ , _rebuild = true);
            ts_e_->updateFromDisp(u_e_);

            std::cout << "***** Solve displacement *****" << std::endl;
            u_.on(_range=elements(support(Xhv_)),_expr=idv(u_e_) + idv(u_r_));
        }
        else
        {
            std::cout << "At iteration : " << iteration << " contact" << std::endl;

            int fixedPointIteration = 0;
            double fixedPointerror = 0.;

            std::cout << "***** Solve rigid *****" << std::endl;
            lt_r_ = l_r_;
            at_r_ = a_r_;

            lt_r_ += integrate( _range = markedelements( support( VhvC_ ), "Caoutchouc" ), _expr = cst( rho_ ) * inner( idv( ts_r_->polyDeriv() ), id( u_r_ ) ) );

            std::cout << "***** Solve elasticity *****" << std::endl;
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
            
                std::cout << "Fixed point iteration : " << fixedPointIteration << std::endl;
                u_e_tmp.on( _range=elements(mesh_), _expr = idv(u_e_tmpNew)); ;
                u_r_tmp.on( _range=elements(mesh_), _expr = idv(u_r_tmpNew)); ;
                u_tmp.on( _range=elements(mesh_), _expr = idv(u_tmpNew)); ;
            
                std::cout << "***** Compute new contact region *****" << std::endl;
                myelts_ = getContactRegion(u_);
                std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
                auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xh_) ), _update=0 );
                auto XhCFaces = Pdh<0>(face_mesh);
                auto contactFaces = XhCFaces->element();
                contactFaces.on( _range=myelts_, _expr = cst(1.));

                lt_tmp_r_ = lt_r_;

                //auto epsv = sym(gradv(u_e_tmp));
                auto epsv = sym(gradv(u_tmp));
                auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();

                lt_tmp_r_ +=  integrate( _range=boundaryfaces(support(VhvC_)), _expr= inner( sigmav,id( u_r_tmp )) * idv(contactFaces));
                at_r_.solve( _rhs = lt_tmp_r_, _solution = u_r_tmpNew, _rebuild = true );

                
                at_tmp_e = at_e;
                lt_tmp_e = lt_e;

                at_tmp_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u_e_tmp),trans(expr<Dim,1>(direction_))*id(u_e_tmp)) * idv(contactFaces));
                lt_tmp_e += integrate (_range=boundaryfaces(support(Xh_)),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_) - trans(expr<Dim,1>(direction_))*idv(u_r_tmp), trans(expr<Dim,1>(direction_))*id(u_e_tmp)) * idv(contactFaces));
                        
                at_tmp_e+=on(_range=markedpoints(mesh_,"CM"), _rhs=lt_tmp_e, _element=u_e_tmp, _expr=0.*one());
                at_tmp_e.solve( _rhs = lt_tmp_e, _solution = u_e_tmpNew , _rebuild = true);

                std::cout << "***** Solve displacement *****" << std::endl;
                u_tmpNew.on(_range=elements(support(Xhv_)),_expr=idv(u_e_tmpNew) + idv(u_r_tmpNew));

                std::cout << "***** Compute error *****" << std::endl;
                fixedPointerror = integrate(_range=elements(support(Xhv_)), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(support(Xhv_)),_expr=norm2(idv(u_))).evaluate()(0,0); 
            
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
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xh_) ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces.on( _range=myelts_, _expr = cst(1.));

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
        contactPressure.on(_range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
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
    if (fixedPoint_ == 1)
        timeLoopFixedPoint();
    else 
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
    e_->step(t)->add( "g", g_ );
    e_->save();
}

}