#pragma once
#include "qs_elasticity_contact.hpp"

namespace Feel
{

template <int Dim, int Order, int OrderGeo>
class ContactLagrange
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
    using ts_ptrtype = std::shared_ptr<NewmarkContact<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_tP1>>; 

    // Constructors
    ContactLagrange() = default;
    ContactLagrange(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    elementv_t const& u() const { return u_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }

    // Methods
    void runDynamic();
    void runStatic();
    Range<mesh_t, MESH_FACES> getContactRegion( elementv_t const& u );
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
    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string direction_;
    std::vector<double> ddirection_;
    int withMarker_,save_, computeerror_;
    double fixedPointtol_;
    int fixedPoint_; 
    std::vector<double> pressurePoint_;
};

// Constructor
template <int Dim, int Order, int OrderGeo>
ContactLagrange<Dim, Order, OrderGeo>::ContactLagrange(nl::json const& specs) : specs_(specs)
{
}


// Run method 
template <int Dim, int Order, int OrderGeo>
void ContactLagrange<Dim, Order, OrderGeo>::runDynamic()
{
    if (Order != 2)
        std::cout << "u has to be P2, to satisfy inf-sup condition" << std::endl;

    std::cout << " ***** Initialize parameters ***** " << std::endl;
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    
    Xhv_ = Pchv<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );
    Xh_ = Pch<Order>( mesh_, markedelements( mesh_, "Caoutchouc" ) );

    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
    
    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));

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

    auto e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );

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
    
    ts_ = newmarkContact(Xhv_, steady, initial_time, final_time, time_step, gamma, beta );
    
    ts_->start();
    ts_->initialize( u0_ );
    u_ = u0_;

    // Export initial time
    e_->step(0)->addRegions();
    e_->step(0)->add("displacement",u_);
    e_->save();
    

    LOG(INFO) << "The step is  " << ts_->timeStep() << "\n"
              << "The initial time is " << ts_->timeInitial() << "\n"
              << "The final time is " << ts_->timeFinal() << "\n";

    ts_->updateFromDisp(u_);

    nbrFaces_ = 0;
    contactRegion_ =  project(_space=Xh_, _range=elements(support(Xhv_)), _expr = cst(0.));
    
    // Get contact parameters
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

    std::string matpressurePoint = fmt::format( "/Collision/LinearElasticity/pressurePoint" );
    pressurePoint_ = specs_[nl::json::json_pointer( matpressurePoint )].get<std::vector<double>>();  

    initG();
    
    // Outputs
    std::ofstream ofs("outputs.csv");
    ofs << fmt::format( "time, evaluateStress, evaluateDisp") << std::endl;


    std::cout << " ***** Initialize elements  and linear, bilinear forms***** " << std::endl;
    auto a_ = form2(_test=Xhv_, _trial=Xhv_);
    auto l_ = form1(_test=Xhv_);

    a_.zero();
    l_.zero();

    l_+=integrate(_range=elements(support(Xhv_)), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_));

    auto deft = sym(gradt(u_));
    auto def = sym(grad(u_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a_+= integrate(_range=elements(support(Xhv_)), _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u_),id( u_ ) ) + inner(sigmat,def));

    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;
    
        // Compute faces in contact
        std::cout << "Compute contact region" << std::endl;
        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;

        if (nbrFaces_ > 0)
        {
            std::cout << "Contact" << std::endl;

            double fixedPointerror = 0.;
            int fixedPointIteration = 0;
            fixedPointtol_ = 0.0001;

            auto u_tmp =  Xhv_->element();
            u_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_));

            auto u_tmpNew = Xhv_->element();
            u_tmpNew.on( _range=elements(support(Xhv_)), _expr = idv(u_)); 


            while ((fixedPointerror > fixedPointtol_) || (fixedPointIteration < 1))
            {
                std::cout << "Fixed Point iteration : " <<  fixedPointIteration << std::endl;
                std::cout << "Error : " << fixedPointerror << std::endl;


                u_tmp.on( _range=elements(support(Xhv_)), _expr = idv(u_tmpNew)); ;


                std::cout << "Compute contact region" << std::endl;
                myelts_ = getContactRegion(u_tmp);
                std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
    
                // Define lagrange multiplier
                auto submesh = createSubmesh( _mesh=mesh_, _range=myelts_, _context = EXTRACTION_KEEP_MESH_RELATION, _update = 0 );
                auto XhLambda = Pch<1>(submesh);
                auto lambda = XhLambda->element();

                auto ps = product( Xhv_, XhLambda );

                // Process loading and materials
                std::cout << "create a and rhs" << std::endl;
                auto at_ = blockform2( ps ,solve::strategy::monolithic,backend() );
                auto lt_ = blockform1( ps,solve::strategy::monolithic,backend() );


                std::cout << "add RHS" << std::endl;

                std::cout << "add RHS_1" << std::endl;
                lt_(0_c) += integrate(_range=elements(support(Xhv_)), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_tmp));
            
                std::cout << "add RHS_2" << std::endl;
                for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
                    lt_(0_c) +=  integrate( _range=elements( support(Xhv_)), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_tmp ) ));

                std::cout << "add RHS_3" << std::endl;
                lt_(1_c) += integrate(_range=elements(submesh), _expr=idv(g_)*id(lambda));

                std::cout << "add A" << std::endl;
                auto deft = sym(gradt(u_tmp));
                auto def = sym(grad(u_tmp));
                auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

                std::cout << "add A_1" << std::endl;
                at_(0_c,0_c) += integrate(_range=elements(support(Xhv_)), _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u_tmp),id( u_tmp ) ) + inner(sigmat,def));
                std::cout << "add A_2" << std::endl;
                at_(0_c,1_c) += integrate( _range=elements(submesh),_expr=idt(lambda)*(trans(expr<Dim,1>(direction_))*id(u_tmp)));
                std::cout << "add A_3" << std::endl;
                at_(1_c,0_c) += integrate( _range=elements(submesh),_expr=(trans(expr<Dim,1>(direction_))*idt(u_tmp))*id(lambda) );

            
                std::cout << "***** Solve *****" << std::endl;
                auto U=ps.element();
                at_.solve( _solution=U, _rhs=lt_, _rebuild=true);

                u_tmpNew = U(0_c);
                fixedPointerror = integrate(_range=elements(support(Xhv_)), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(support(Xhv_)),_expr=norm2(idv(u_))).evaluate()(0,0); 
                fixedPointIteration++;

                if (fixedPointIteration == 10)
                    break;

            }

            u_.on( _range=elements(support(Xhv_)), _expr = idv(u_tmpNew)); 
            ts_->updateFromDisp(u_);
            
        }
        else 
        {
            std::cout << "No Contact" << std::endl;
            
            auto at_ = form2( _test = Xhv_, _trial = Xhv_ );
            auto lt_ = form1( _test = Xhv_ ); 

            at_.zero();
            lt_.zero();

            at_ = a_;
            lt_ = l_;

            for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
                lt_ +=  integrate( _range=elements( support(Xhv_)), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ));
        

            std::cout << "***** Solve *****" << std::endl;
            at_.solve( _rhs=lt_, _solution=u_, _rebuild=true);

            ts_->updateFromDisp(u_);
        }

        // Exports
        std::cout << "***** Exports *****" << std::endl;
        
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

        auto epsv = sym(gradv(u_));
        auto sigmav = (lambda_*trace(epsv)*Id + 2*mu_*epsv)*N();

        std::cout << "***** Exports - contactPressure *****" << std::endl;
        auto contactPressure = Xh_->element();
        contactPressure.on(_range=boundaryfaces(support(Xhv_)), _expr = trans(expr<Dim,1>(direction_))*sigmav);
        auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(contactPressure) ); 

        std::cout << "***** Exports - evaluateDispExpr *****" << std::endl;
        auto evaluateDispExpr = Xh_->element();
        evaluateDispExpr.on(_range=elements(support(Xhv_)), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= idv(evaluateDispExpr) );     
    
        if ( Dim == 2 )
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",
                                ts_->time(), evaluateStress(0,0), evaluateDisp(0,0))
                << std::endl;
        else
            ofs << fmt::format( "{:.6f}, {:.6f}, {:.6f}",
                                ts_->time(), evaluateStress(0,0), evaluateDisp(0,0) )
                << std::endl;
        
        e_->step(ts_->time())->addRegions();
        e_->step(ts_->time())->add("displacement",u_);
        e_->save();
    }

    ofs.close();
   
}

template <int Dim, int Order, int OrderGeo>
Range<typename ContactLagrange<Dim, Order, OrderGeo>::mesh_t, MESH_FACES>
ContactLagrange<Dim, Order, OrderGeo>::getContactRegion( elementv_t const& u )
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
ContactLagrange<Dim, Order, OrderGeo>::initG()
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


    auto e = Feel::exporter(_mesh = mesh_, _name = "InitialDistance" );
    e->addRegions();
    e->add( "g", g_ );
    e->save();

}
   
template <int Dim, int Order, int OrderGeo>
void ContactLagrange<Dim, Order, OrderGeo>::runStatic()
{
    if (Order != 2)
        std::cout << "u has to be P2, to satisfy inf-sup condition" << std::endl;

    std::cout << " ***** Initialize parameters ***** " << std::endl;

    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    Xhv_ = Pchv<Order>(mesh_); 
    Xh_ = Pch<Order>(mesh_);

    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
    std::string matE = fmt::format( "/Materials/Caoutchouc/parameters/E/value" );
    double E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
    std::string matNu = fmt::format( "/Materials/Caoutchouc/parameters/nu/value" );
    double nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
    
    lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
    mu_ = E_/(2*(1+nu_));
    
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

    auto e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );

    nbrFaces_ = 0;
    
    // Get contact parameters
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
    
    std::string matcomputeerror = fmt::format( "/Collision/LinearElasticity/error" );
    computeerror_ = specs_[nl::json::json_pointer( matcomputeerror )].get<int>(); 
    
    std::string matsave = fmt::format( "/Collision/LinearElasticity/save" );
    save_ = specs_[nl::json::json_pointer( matsave )].get<int>();    
    

    std::cout << " ***** Initialize elements  and linear, bilinear forms***** " << std::endl;
    // Define Lagrange multiplier
    auto submesh = createSubmesh(_mesh=mesh_,_range=markedfaces(mesh_,"contact"));
    auto XhLambda = Pch<1>(submesh);

    //  Initialize elements
    auto u = Xhv_->elementPtr();
    auto lambda = XhLambda->elementPtr();

    BlocksBaseGraphCSR myblockGraph(2,2);
    myblockGraph(0,0) = stencil( _test=Xhv_,_trial=Xhv_, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(1,0) = stencil( _test=XhLambda,_trial=Xhv_, _diag_is_nonzero=false, _close=false)->graph();
    myblockGraph(0,1) = stencil( _test=Xhv_,_trial=XhLambda, _diag_is_nonzero=false, _close=false)->graph();
    auto A = backend()->newBlockMatrix(_block=myblockGraph);


    BlocksBaseVector<double> myblockVec(2);
    myblockVec(0,0) = backend()->newVector( Xhv_ );
    myblockVec(1,0) = backend()->newVector( XhLambda );
    auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

    BlocksBaseVector<double> myblockVecSol(2);
    myblockVecSol(0,0) = u;
    myblockVecSol(1,0) = lambda;
    auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);

    std::cout << "***** Process loading *****" << std::endl;
    form1( _test=Xhv_, _vector=F ) = integrate(_range=elements(mesh_), _expr=cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u));

    std::cout << "***** Process materials *****" << std::endl;
    auto deft = sym(gradt(u));
    auto def = sym(grad(u));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    auto a = form2(_trial=Xhv_,_test=Xhv_,_matrix=A);
    a += integrate(_range=elements(mesh_),_expr=inner(sigmat,def) );

    form2( _trial=XhLambda, _test=Xhv_ ,_matrix=A, _rowstart=0, _colstart=1 ) += integrate( _range=elements(submesh),_expr=idt(lambda)*(trans(expr<Dim,1>(direction_))*id(u)));
    form2( _trial=Xhv_, _test=XhLambda ,_matrix=A,_rowstart=1, _colstart=0 ) += integrate( _range=elements(submesh),_expr=(trans(expr<Dim,1>(direction_))*idt(u))*id(lambda) );

    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Dirichlet") )
    {
        for ( auto [key, bc] : specs_["/BoundaryConditions/LinearElasticity/Dirichlet"_json_pointer].items() )
        {
            std::string e = fmt::format("/BoundaryConditions/LinearElasticity/Dirichlet/{}/g/expr",key);
            auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
            a+=on(_range=markedfaces(mesh_,key), _rhs=F, _element=*u, _expr=expr<Dim,1>( bc_dir ) );
        }
    }

    std::cout << "***** Solve *****" << std::endl;
    backend(_rebuild=true)->solve( _matrix=A, _rhs=F, _solution=U );
    myblockVecSol.localize(U);

    e_->addRegions();
 
    auto Xhv_P1 = Pchv<OrderGeo>(mesh_); 
    auto Xh_P1 = Pch<OrderGeo>(mesh_);

    auto uinter =  Xhv_P1->element(); 
    auto op_inter = opInterpolation(_domainSpace =  Xhv_, _imageSpace = Xhv_P1 );
    op_inter->apply(*u, uinter);

    e_->add( "displacement", uinter );
    e_->save();

    
    Range<mesh_t,MESH_FACES> myelts( mesh_ );
    auto realcontactRegion = project(_space=Xh_P1, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(uinter));  
       
    auto const& trialDofIdToContainerId =  form2(_test=Xh_P1, _trial=Xh_P1).dofIdToContainerIdTest();
    for (auto const& theface : markedfaces(mesh_,"contact") )
    {
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : Xh_P1->dof()->faceLocalDof( face.id() ) )
        {
            index_type thedof = ldof.index();
            thedof = trialDofIdToContainerId[ thedof ];

            if (realcontactRegion[thedof] >= tolContactRegion_)
                contactDof++;

            if (Dim == 2)
            {
                if (contactDof == 2)
                {
                    nbrFaces_++;
                    myelts.push_back( face  );
                }    

            }
            else if (Dim == 3)
            {
                if (contactDof == 3)
                {
                    nbrFaces_++;
                    myelts.push_back( face  );
                }
                    
            }
        }
    }
    myelts.shrink_to_fit();

    std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
    
    if (save_ == 1)
    {
        std::cout << "***** Save results *****" << std::endl;
        (*u).saveHDF5("solution_ref.h5");
        saveGMSHMesh(_mesh=mesh_,_filename= "mesh_ref.msh" );
    }

    if (computeerror_ == 1)
    {
        std::cout << "***** Compute error *****" << std::endl;
        auto meshref = loadMesh(_mesh=new mesh_t, _filename = "/data/scratch/vanlandeghem/feel/qs_elasticity_contact/square/np_1/mesh_ref.msh");
        auto uref = spacev_t::New(meshref)->elementPtr() ;
        uref->loadHDF5("/data/scratch/vanlandeghem/feel/qs_elasticity_contact/square/np_1/solution_ref.h5");

        auto op_inter = opInterpolation(_domainSpace =  spacev_t::New(mesh_), _imageSpace = spacev_t::New(meshref) );
        auto uinter =  spacev_t::New(meshref)->element();
        op_inter->apply(*u, uinter);

        auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
        std::cout << "H1 relative error : " << h1err << std::endl;

    }
}

} 