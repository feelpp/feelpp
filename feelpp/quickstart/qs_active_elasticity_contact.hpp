#pragma once
#include "qs_active_elasticity.hpp"

namespace Feel
{

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "active elasticity contact options" );
    options.add_options()

        // mesh parameters
        ( "specs", Feel::po::value<std::string>(),
          "json spec file for rht" )

        ( "steady", Feel::po::value<bool>()->default_value( 1 ),
          "if 1: steady else unsteady" );

    return options.add( Feel::feel_options() );
}

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
    void run();
    void timeLoopActive();
    void timeLoopHyper();
    void exportResults(double t);
    void initG();
    Range<mesh_t, MESH_FACES> getContactRegion( elementv_t const& u );


private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t Xhv_;
    space_ptr_t Xh_;

    elementv_t u_;
    element_t contactRegion_;
    element_t g_;
    int nbrFaces_;
    Range<mesh_t, MESH_FACES> myelts_;

    ts_ptrtype ts_;
    exporter_ptrtype e_;

    int active_;
    double Ca_, Lc_, rc_, fa_, va_;
    std::string type_;

    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;

    double epsilon_, tolContactRegion_, tolDistance_, gamma0_, gamma_;
    std::string direction_, method_;
    std::vector<double> ddirection_;
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
    H_ = specs_["/Meshes/HyperElasticity/Import/h"_json_pointer].get<double>();
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/HyperElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);

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

    std::string matactive = fmt::format("/Models/HyperElasticity/active");
    active_ = specs_[nl::json::json_pointer( matactive )].get<int>();

    // External force
    if ( specs_["/Models/HyperElasticity"_json_pointer].contains("loading") )
    {
        for ( auto [key, loading] : specs_["/Models/HyperElasticity/loading"_json_pointer].items() )
        {
            std::string loadtype = fmt::format( "/Models/HyperElasticity/loading/{}/type", key );

            if ( specs_[nl::json::json_pointer( loadtype )].get<std::string>() == "Gravity" )
            {
                LOG( INFO ) << fmt::format( "Loading {}: Gravity found", key );
                std::string loadexpr = fmt::format( "/Models/HyperElasticity/loading/{}/parameters/expr", key );
                externalforce_ = specs_[nl::json::json_pointer( loadexpr )].get<std::string>();
            }
        }
    }
    std::cout << "External force : " << externalforce_ << std::endl;

    // Exporter
    e_ = Feel::exporter(_mesh = mesh_, _name = specs_["/ShortName"_json_pointer].get<std::string>() );
    
    // Newmark scheme
    bool steady = get_value(specs_, "/TimeStepping/HyperElasticity/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/HyperElasticity/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/HyperElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/HyperElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/HyperElasticity/step", std::string("0.1"))).evaluate()(0,0);
    double gamma = get_value(specs_, "/TimeStepping/HyperElasticity/gamma", 0.5);
    double beta = get_value(specs_, "/TimeStepping/HyperElasticity/beta", 0.25);

    // Initial conditions
    u_ = Xhv_->element();
    auto u0_ = Xhv_->element();

    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    auto init_displ = expr<Dim,1>(get_value(specs_, "/InitialConditions/HyperElasticity/displacement/expr", default_displ ));
    u0_.on(_range=elements(support(Xhv_)), _expr=init_displ);
    
    ts_ = newmarkContact(Xhv_, steady, initial_time, final_time, time_step, gamma, beta );
    ts_->start();
    ts_->initialize( u0_ );
    u_ = u0_;
    
    ts_->updateFromDisp(u_);

    nbrFaces_ = 0;
    contactRegion_ =  project(_space=Xh_, _range=elements(support(Xhv_)), _expr = cst(0.));
    
    std::string matEpsilon = fmt::format( "/Collision/HyperElasticity/epsilon" );
    epsilon_ = specs_[nl::json::json_pointer( matEpsilon )].get<double>(); 

    std::string matGamma0 = fmt::format( "/Collision/HyperElasticity/gamma0" );
    gamma0_ = specs_[nl::json::json_pointer( matGamma0 )].get<double>(); 
    gamma_ = gamma0_/H_;
    
    std::string matDirection = fmt::format( "/Collision/HyperElasticity/direction");
    direction_ = specs_[nl::json::json_pointer( matDirection )].get<std::string>();

    std::string matDirectionD = fmt::format( "/Collision/HyperElasticity/ddirection");
    ddirection_ = specs_[nl::json::json_pointer( matDirectionD )].get<std::vector<double>>();

    std::string mattolContactRegion = fmt::format( "/Collision/HyperElasticity/tolContactRegion" );
    tolContactRegion_ = specs_[nl::json::json_pointer( mattolContactRegion )].get<double>();  

    std::string mattolDistance = fmt::format("/Collision/HyperElasticity/tolDistance");
    tolDistance_ = specs_[nl::json::json_pointer( mattolDistance )].get<double>();  

    std::string matMethod = fmt::format( "/Collision/HyperElasticity/method");
    method_ = specs_[nl::json::json_pointer( matMethod )].get<std::string>();
    
}

template <int Dim, int Order>
void ActiveContact<Dim, Order>::initG()
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

template <int Dim, int Order>
Range<typename ActiveContact<Dim, Order>::mesh_t, MESH_FACES>
ActiveContact<Dim, Order>::getContactRegion( elementv_t const& u )
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

// Time loop active elasticity
template <int Dim, int Order>
void ActiveContact<Dim, Order>::timeLoopActive()
{
    // Initialize parameters
    std::string matCa = fmt::format("/Materials/Caoutchouc/parameters/Ca/value");
    Ca_ = std::stod(specs_[nl::json::json_pointer( matCa )].get<std::string>());

    std::string matLc = fmt::format("/Materials/Caoutchouc/parameters/Lc/value");
    Lc_ = std::stod(specs_[nl::json::json_pointer( matLc )].get<std::string>());

    std::string matRc = fmt::format("/Materials/Caoutchouc/parameters/rc/value");
    rc_ = std::stod(specs_[nl::json::json_pointer( matRc )].get<std::string>());

    std::string matFa = fmt::format("/Materials/Caoutchouc/parameters/fa/value");
    fa_ = std::stod(specs_[nl::json::json_pointer( matFa )].get<std::string>());

    std::string matVa = fmt::format("/Materials/Caoutchouc/parameters/va/value");
    va_ = std::stod(specs_[nl::json::json_pointer( matVa )].get<std::string>());
    
    std::string matType = fmt::format("/Materials/Caoutchouc/parameters/type/value");
    type_ = specs_[nl::json::json_pointer( matType )].get<std::string>();


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

        std::cout << "***** Process contact *****" << std::endl;
        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        
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
            
            auto ea = vec(cst(0.), cst(1.));
            auto eaea = ea*trans(ea);
            
            
            auto a = form2( _test=Xhv_, _trial=Xhv_, _matrix=J );

            a = integrate( _range=elements(support(Xhv_)),
                           _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
            
            a += integrate( _range=elements(support(Xhv_)),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );

            if (type_.compare("bending") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
                a += integrate(_range=elements( support(Xhv_) ), _expr= - inner(dF*sigma_a*Px()*eaea, grad(u) ) );
            }    
            else if (type_.compare("flapping") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_->time()));
                a += integrate(_range=elements( support(Xhv_) ), _expr= - inner(dF*sigma_a*Px()*eaea, grad(u) ) );
            }
            
            if (nbrFaces_ > 0)
            {
                std::cout << "Add contact terms to Jacobian" << std::endl;
                auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xhv_) ), _update=0 );
                auto XhCFaces = Pdh<0>(face_mesh);
                auto contactFaces = XhCFaces->element();
                contactFaces.on( _range=myelts_, _expr = cst(1.));

                if (method_.compare("penalty") == 0)
                {
                    a += integrate(_range=boundaryfaces(support(Xhv_) ), _expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }    
                else if (method_.compare("nitsche") == 0)
                {
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(gamma_)  * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(gamma_)  * inner(trans(expr<Dim,1>(direction_))*dF*val(Sv)*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(gamma_)  * inner(trans(expr<Dim,1>(direction_))*val(Fv)*dS*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                
                    
                    if (type_.compare("bending") == 0)
                    {
                        auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
                        a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*dF*sigma_a*Px()*eaea*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    }    
                    else if (type_.compare("flapping") == 0)
                    {
                        auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_->time()));
                        a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*dF*sigma_a*Px()*eaea*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    }
                    
                }
            }

            auto RR = backend()->newVector( Xhv_ );

            a += on( _range=markedpoints(mesh_,"dirichlet"),
                     _element=u, _rhs=RR,
                     _expr=zero<Dim,1>() );
            
            
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
            
            if (type_.compare("bending") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
                r += integrate(_range=elements( support(Xhv_)), _expr= - inner( val(Fv)*sigma_a*Px()*eaea,grad( u )));
            }    
            else if (type_.compare("flapping") == 0)
            {
                auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_->time()));
                r += integrate(_range=elements( support(Xhv_)), _expr= - inner( val(Fv)*sigma_a*Px()*eaea,grad( u )));
            }

            if (nbrFaces_ > 0)
            {
                std::cout << "Add contact terms to Residual" << std::endl;
                auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xhv_) ), _update=0 );
                auto XhCFaces = Pdh<0>(face_mesh);
                auto contactFaces = XhCFaces->element();
                contactFaces.on( _range=myelts_, _expr = cst(1.));


                if (method_.compare("penalty") == 0)
                {
                    std::cout << "Penalty method" << std::endl;
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(epsilon_) * inner(idv(g_),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }    
                else if (method_.compare("nitsche") == 0)
                {
                    std::cout << "Nitsche method" << std::endl;
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(gamma_)  *  inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idv(u),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(gamma_)  * inner(cst(gamma_) * idv(g_),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(gamma_)  * inner(trans(expr<Dim,1>(direction_))*val(Fv*Sv)*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    
                    
                    if (type_.compare("bending") == 0)
                    {
                        auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*ts_->time());
                        r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr=  cst(1.)/cst(gamma_)  * inner(trans(expr<Dim,1>(direction_))*val(Fv)*sigma_a*Px()*eaea*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    }    
                    else if (type_.compare("flapping") == 0)
                    {
                        auto sigma_a = Ca_/(Lc_*rc_)*sin(2*pi*fa_*(Py() - va_*ts_->time()));
                        r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr=  cst(1.)/cst(gamma_)  * inner(trans(expr<Dim,1>(direction_))*val(Fv)*sigma_a*Px()*eaea*N(),cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    }
                    
                }
            }

            R->close();
            auto temp = Xhv_->element();
            temp = *R;
            temp.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
            temp.on( _range=markedpoints(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
            *R = temp;
        };

        u_.on(_range=markedpoints(mesh_,"dirichlet"), _expr=zero<Dim,1>());
        u_.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
        backend()->nlSolver()->residual = Residual;
        backend()->nlSolver()->jacobian = Jacobian;
        backend()->nlSolve( _solution=u_,_jacobian=Jac,_residual=Res );


        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

    }    
}

// Time loop hyper elasticity
template <int Dim, int Order>
void ActiveContact<Dim, Order>::timeLoopHyper()
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

        
        std::cout << "***** Process contact *****" << std::endl;
        myelts_ = getContactRegion(u_);
        std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        
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
            
            auto a = form2( _test=Xhv_, _trial=Xhv_, _matrix=J );

            a = integrate( _range=elements(support(Xhv_)),
                           _expr= inner( dF*val(Sv) + val(Fv)*dS , grad(u) ) );
            
            a += integrate( _range=elements(support(Xhv_)),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idt(u),id( u ) ) );

            if (nbrFaces_ > 0)
            {
                std::cout << "Add contact terms to Jacobian" << std::endl;
                auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xhv_) ), _update=0 );
                auto XhCFaces = Pdh<0>(face_mesh);
                auto contactFaces = XhCFaces->element();
                contactFaces.on( _range=myelts_, _expr = cst(1.));
                
                if (method_.compare("penalty") == 0)
                {
                    a += integrate(_range=boundaryfaces(support(Xhv_) ), _expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }    
                else if (method_.compare("nitsche") == 0)
                {
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - inner(trans(expr<Dim,1>(direction_))*dF*val(Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    a += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv)*dS*N(),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }
            }

            auto RR = backend()->newVector( Xhv_ );

            a += on( _range=markedpoints(mesh_,"dirichlet"),
                     _element=u, _rhs=RR,
                     _expr=zero<Dim,1>() );
            
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

            auto r = form1( _test=Xhv_, _vector=R );
            r = integrate( _range=elements(support(Xhv_)),
                           _expr= inner( val(Fv*Sv) , grad(u) ) );
            r += integrate( _range=elements(support(Xhv_)),
                            _expr= - trans( expr<Dim, 1>( externalforce_ ) )*id( u )  );
            r += integrate( _range=elements(support(Xhv_)),
                            _expr= cst(rho_)*inner( ts_->polyDerivCoefficient()*idv(u) -idv(ts_->polyDeriv()),id( u ) ) );
            
            if (nbrFaces_ > 0)
            {
                std::cout << "Add contact terms to Residual" << std::endl;
                auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xhv_) ), _update=0 );
                auto XhCFaces = Pdh<0>(face_mesh);
                auto contactFaces = XhCFaces->element();
                contactFaces.on( _range=myelts_, _expr = cst(1.));

                if (method_.compare("penalty") == 0)
                {
                    std::cout << "Penalty method" << std::endl;
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(1.)/cst(epsilon_) * inner(idv(g_),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }    
                else if (method_.compare("nitsche") == 0)
                {
                    std::cout << "Nitsche method" << std::endl;
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*idv(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - cst(gamma_) * inner(idv(g_),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                    r += integrate (_range=boundaryfaces(support(Xhv_) ),_expr= - inner(trans(expr<Dim,1>(direction_))*val(Fv*Sv)*N(),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
                }
            }

            R->close();
            auto temp = Xhv_->element();
            temp = *R;
            temp.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
            temp.on( _range=markedpoints(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
            
            *R = temp;
        };

        u_.on( _range=markedfaces(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
        u_.on( _range=markedpoints(mesh_,"dirichlet"),_expr=zero<Dim,1>() );
        
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

    std::cout << "***** Fix initial distance *****" << std::endl;
    this->initG();

    std::cout <<  "***** Start time loop *****" << std::endl;
    this->exportResults(0);

    if (active_ == 1)
        this->timeLoopActive();
    else 
        this->timeLoopHyper();
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

    auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(support(Xhv_) ), _update=0 );
    auto XhCFaces = Pdh<0>(face_mesh);
    auto contactFaces = XhCFaces->element();
    contactFaces.on( _range=myelts_, _expr = cst(1.));
    auto region = project(_space=Xh_,_range=boundaryfaces(support(Xhv_)), _expr = cst(1)*idv(contactFaces));

    e_->step(t)->add("contactFaces",region);

    if (active_ == 1)
    {
        auto sigma_a = Ca_/(Lc_*rc_)*std::sin(2*pi*fa_*t);
        auto activity = project(_space=Xh_, _range=elements(support(Xh_)), _expr = sigma_a*(-Px()));
        e_->step(t)->add("activity",activity);
    }
    e_->save();
}

} 