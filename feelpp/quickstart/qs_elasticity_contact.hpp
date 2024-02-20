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
#include <feel/feelts/newmark.hpp>


namespace Feel
{
inline const int FEELPP_DIM=2;
inline const int FEELPP_ORDER=2;


inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "elasticity contact options" );
    options.add_options()

        // mesh parameters
        ( "specs", Feel::po::value<std::string>(),
          "json spec file for rht" )

        ( "steady", Feel::po::value<bool>()->default_value( 1 ),
          "if 1: steady else unsteady" );

    return options.add( Feel::feel_options() );
}

template<typename T>
T get_value(const nl::json& specs, const std::string& path, const T& default_value)
{
    auto json_pointer = nl::json::json_pointer(path);
    return specs.contains(json_pointer) ? specs[json_pointer].get<T>() : default_value;
}

template <int Dim, int Order>
class ElasticContact
{
public:
    using mesh_t = Mesh<Simplex<Dim>>;
    using space_t = Pchv_type<mesh_t, Order>;
    using spaceC_t = Pch_type<mesh_t, Order>;
    using space_ptr_t = Pchv_ptrtype<mesh_t, Order>; 
    using spaceC_ptr_t = Pch_ptrtype<mesh_t, Order>;
    using element_t = typename space_t::element_type;
    using elementC_t = typename spaceC_t::element_type;
    using form2_type = form2_t<space_t,space_t>; 
    using form1_type = form1_t<space_t>; 
    using ts_ptrtype = std::shared_ptr<Newmark<space_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

    // Constructors
    ElasticContact() = default;
    ElasticContact(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    space_ptr_t const& Xh() const { return Xh_; }
    spaceC_ptr_t const XhC() const { return XhC_; } 
    element_t const& u() const { return u_; }
    element_t const& v() const { return v_; }
    elementC_t const& contactRegion() const { return contactRegion_; }
    int const& nbrFaces() const { return nbrFaces_; }
    form2_type const& a() const { return a_; }
    form2_type const& at() const { return at_; }
    form1_type const& l() const { return l_; }
    form1_type const& lt() const { return lt_; }
    ts_ptrtype const& bdf() const { return ts_; }
    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }


    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(element_t const& u) { u_ = u; }
    void setContactRegion(elementC_t const&c) { contactRegion_ = c; }

    // Methods
    void initialize();
    void initializeContact();
    void processLoading(form1_type& l);
    void processMaterials(form2_type &a);
    void processBoundaryConditions(form1_type& l, form2_type& a);
    void processContact(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts, element_t const& u, element_t const& v);
    void run();
    typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype getContactRegion(element_t const& u);
    void timeLoop();
    void timeLoopfixedPoint();
    void exportResults( double t );
    void writeResultsToFile(const std::string& filename) const;


private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    space_ptr_t Xh_;
    spaceC_ptr_t XhC_;
    element_t u_, v_;
    elementC_t contactRegion_;
    elementC_t pressureRegion_;
    elementC_t contactPressure_;
    elementC_t contactDisplacement_;
    elementC_t contactFaces_;
    int nbrFaces_;
    form2_type a_, at_, at_tmp;
    form1_type l_, lt_, lt_tmp;
    ts_ptrtype ts_;
    exporter_ptrtype e_;
    nl::json meas_;
    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;
    int fixedPoint_;
    double fixedPointtolerance_;
    double distance_, theta_, gamma0_, gamma_;
    std::string direction_;
    // double beta, gamma;
    //std::string F, G;
};

// Constructor
template <int Dim, int Order>
ElasticContact<Dim, Order>::ElasticContact(nl::json const& specs) : specs_(specs)
{
}

// Initialization
template <int Dim, int Order>
void ElasticContact<Dim, Order>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // define Xh
    Xh_ = Pchv<Order>(mesh_); 

    // Initialize elements
    u_ = Xh_->element();
    v_ = Xh_->element();

    // Initialize linear and bilinear forms
    a_ = form2( _test = Xh_, _trial = Xh_ );
    at_ = form2( _test = Xh_, _trial = Xh_ );
    l_ = form1( _test = Xh_ );
    lt_ = form1( _test = Xh_ );
    at_tmp = form2( _test = Xh_, _trial = Xh_ );
    lt_tmp = form1( _test = Xh_ );
    // Initialize time stepping
    bool steady = get_value(specs_, "/TimeStepping/LinearElasticity/steady", true);
    int time_order = get_value(specs_, "/TimeStepping/LinearElasticity/order", 2);
    double initial_time = get_value(specs_, "/TimeStepping/LinearElasticity/start", 0.0);
    double final_time = get_value(specs_, "/TimeStepping/LinearElasticity/end", 1.0);
    double time_step = expr(get_value(specs_, "/TimeStepping/LinearElasticity/step", std::string("0.1"))).evaluate()(0,0);
    // Set initial conditions
    auto u0_ = Xh_->element();
    std::string default_displ = (Dim==2)?std::string("{0.,0.}"):std::string("{0.,0.,0.}");
    auto init_displ = expr<Dim,1>(get_value(specs_, "/InitialConditions/LinearElasticity/displacement/expr", default_displ ));
    u0_.on(_range=elements(mesh_), _expr=init_displ);

    LOG(INFO) << "**** elasticity contact parameters initialized **** \n" << std::endl;

    // Initialize external forces
    if ( specs_["/Models/LinearElasticity"_json_pointer].contains("loading") )
    {
        for ( auto [key, loading] : specs_["/Models/LinearElasticity/loading"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Loading {} found", key );
            std::string loadtype = fmt::format( "/Models/LinearElasticity/loading/{}/type", key );
            
            // External forces
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
    ts_ = newmark( _space = Xh_, _steady=steady, _initial_time=initial_time, _final_time=final_time, _time_step=time_step, _order=time_order );

    ////////////////////////////////////////////////////
    //          Newmark beta-model for dttun          //
    ////////////////////////////////////////////////////
    l_.zero();
    a_.zero();
    
    ts_->start();
    if ( steady )
        ts_->setSteady();

    ts_->initialize( u0_ ); // set u0_
    u_ = u0_;

    if ( steady )
        LOG(INFO) << "\n***** Compute Steady state *****" << std::endl;
    else
    {
        LOG(INFO) << "\n***** Compute Transient state *****" << std::endl;
        LOG(INFO) << "The step is  " << ts_->timeStep() << "\n"
                  << "The initial time is " << ts_->timeInitial() << "\n"
                  << "The final time is " << ts_->timeFinal() << "\n";
    }

    at_.zero();
    lt_.zero();

    at_tmp.zero();
    lt_tmp.zero();

    ts_->updateFromDisp(u_);
}

// Initialization contact terms
template <int Dim, int Order>
void ElasticContact<Dim, Order>::initializeContact()
{
    XhC_ = Pch<Order>(mesh_);

    contactRegion_ = XhC_->element();
    pressureRegion_ = XhC_->element();
    contactPressure_ = XhC_->element();
    contactDisplacement_ = XhC_->element();
    contactFaces_ = XhC_->element();

    nbrFaces_ = 0;

    std::string matDistance = fmt::format( "/Collision/LinearElasticity/distance" );
    distance_ = specs_[nl::json::json_pointer( matDistance )].get<double>(); 
    std::string matTheta = fmt::format( "/Collision/LinearElasticity/theta" );
    theta_ = specs_[nl::json::json_pointer( matTheta )].get<double>(); 
    std::string matGamma0 = fmt::format( "/Collision/LinearElasticity/gamma0" );
    gamma0_ = specs_[nl::json::json_pointer( matGamma0 )].get<double>(); 
    std::string matDirection = fmt::format( "/Collision/LinearElasticity/direction");
    direction_ = specs_[nl::json::json_pointer( matDirection )].get<std::string>();
    gamma_ = gamma0_/H_;
    std::string matFixedPoint = fmt::format( "/Collision/LinearElasticity/fixedPoint" );
    fixedPoint_ = specs_[nl::json::json_pointer( matFixedPoint )].get<int>(); 
    std::string matFixedPointtolerance = fmt::format( "/Collision/LinearElasticity/tolerance" );
    fixedPointtolerance_ = specs_[nl::json::json_pointer( matFixedPointtolerance )].get<double>(); 


    this->exportResults(0);
}

// Process loading
template <int Dim, int Order>
void ElasticContact<Dim, Order>::processLoading(form1_type& l)
{
    LOG(INFO) << fmt::format("Loading expr : {} ", externalforce_) << std::endl;
    l += integrate( _range = elements(mesh_), _expr = trans(expr<Dim,1>( externalforce_ ))*id(v_) );
}

// Process materials
template <int Dim, int Order>
void ElasticContact<Dim, Order>::processMaterials( form2_type &a )
{
    for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
    {
        LOG( INFO ) << fmt::format( "Material {} found", material.get<std::string>() );

        std::string matE = fmt::format( "/Materials/{}/parameters/E/value", material.get<std::string>() );
        E_ = std::stod(specs_[nl::json::json_pointer( matE )].get<std::string>());
        std::string matNu = fmt::format( "/Materials/{}/parameters/nu/value", material.get<std::string>() );
        nu_ = std::stod(specs_[nl::json::json_pointer( matNu )].get<std::string>());
        std::string matRho = fmt::format( "/Materials/{}/parameters/rho/value", material.get<std::string>() );
        rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
        lambda_ = E_*nu_/( (1+nu_)*(1-2*nu_) );
        mu_ = E_/(2*(1+nu_));

        std::cout << "Lambda : " << lambda_ << " mu : " << mu_ << std::endl;

        auto deft = sym(gradt(u_));
        auto def = sym(grad(v_));
        auto Id = eye<Dim,Dim>();
        auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

        a += integrate( _range = elements(mesh_), _expr= rho_*inner( ts_->polyDerivCoefficient()*idt(u_),id( v_ ) ) + inner(sigmat,def));
    }
}

// Process contact conditions
template <int Dim, int Order>
void ElasticContact<Dim, Order>::processContact(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts , element_t const& u, element_t const& v )
{
    
    auto const Id = eye<Dim,Dim>();

    auto deft = sym(gradt(u));
    auto def = sym(grad(v));
    auto defv = sym(gradv(u));

    auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();
    auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();

    auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), elts->begin(), elts->end(), elts );

    contactFaces_ = project(_space=XhC_, _range=myfaces, _expr = cst(1.));
                        
    //a += integrate (_range=myfaces,_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma) ); 
    //a += integrate (_range=myfaces,_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u_) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(v_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));
    //l += integrate (_range=myfaces,_expr= inner(abs(cst(distance_) - Py()), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(v_) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));     

    a += integrate (_range=boundaryfaces(mesh_),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma) * idv(contactFaces_)); 
    a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(v) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma) * idv(contactFaces_));
    l += integrate (_range=boundaryfaces(mesh_),_expr= inner(abs(cst(distance_) - Py()), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(v) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma) * idv(contactFaces_));     

}

    


// Process boundary conditions
template <int Dim, int Order>
void ElasticContact<Dim, Order>::processBoundaryConditions(form1_type& l, form2_type& a)
{
    // Boundary Condition Dirichlet
    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Dirichlet") )
    {
        LOG( INFO ) << fmt::format( "Dirichlet conditions found" );
        for ( auto [key, bc] : specs_["/BoundaryConditions/LinearElasticity/Dirichlet"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Dirichlet conditions found: {}", key );
            std::string e = fmt::format("/BoundaryConditions/LinearElasticity/Dirichlet/{}/g/expr",key);
            auto bc_dir = specs_[nl::json::json_pointer( e )].get<std::string>();
            LOG(INFO) << "BoundaryCondition Dirichlet : " << bc_dir << std::endl;
            a+=on(_range=markedfaces(mesh_,key), _rhs=l, _element=u_, _expr=expr<Dim,1>( bc_dir ) );
        }
    }
}

// Run method 
template <int Dim, int Order>
void ElasticContact<Dim, Order>::run()
{
    LOG(INFO) << "\n***** Initialize *****" << std::endl;
    initialize();
    initializeContact();
    LOG(INFO) << "\n***** Time loop *****" << std::endl;
    if (fixedPoint_ == 1)
        timeLoopfixedPoint();
    else 
        timeLoop();
}

template<int Dim, int Order>
typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype
ElasticContact<Dim, Order>::getContactRegion(element_t const& u)
{   
    typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_type );

    auto defv = sym(gradv(u));
    auto Id = eye<Dim,Dim>();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
    
    contactRegion_ = project(_space=XhC_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u) - abs(cst(distance_) - Py()));    
    //contactRegion_ = project(_space=XhC_, _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u) - abs(cst(distance_) - Py()) - (trans(expr<Dim,1>(direction_))*sigmav));


    nbrFaces_ = 0;
    auto const& trialDofIdToContainerId =  form2(_test=XhC_, _trial=XhC_).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(mesh_) )
    {                
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : XhC_->dof()->faceLocalDof( face.id() ) )
        {
            index_type thedof = ldof.index();
            thedof = trialDofIdToContainerId[ thedof ];

            if (contactRegion_[thedof] >= -H_/10.)
                contactDof++;
                    
            if (Order == 1)
            {
                if (Dim == 2)
                {
                    if (contactDof == 2)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                    std::cout << "TODO" << std::endl;
                       
            }
            else if (Order == 2)
            {
                if (Dim == 2)
                {
                    if (contactDof == 2)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                    std::cout << "TODO" << std::endl;
            }
        }
    }
    myelts->shrink_to_fit();   

    contactFaces_ = project(_space=XhC_, _range=elements(mesh_), _expr = cst(0.)); 
    contactPressure_ =  project(_space=XhC_, _range=elements(mesh_), _expr = cst(0.));
    contactDisplacement_ = project(_space=XhC_, _range=elements(mesh_), _expr = cst(0.));

    return myelts;
}


// Time loop
template <int Dim, int Order>
void ElasticContact<Dim, Order>::timeLoop()
{
    int it = 0; 
    processLoading(l_);
    processMaterials(a_);

    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        lt_ = l_;
        at_ = a_;

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Material {} found", material );

            std::string matRho = fmt::format( "/Materials/{}/parameters/rho/value", material.get<std::string>() );
            auto rho = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
            lt_ +=  integrate( _range=markedelements( mesh_, material.get<std::string>() ), _expr= rho*inner( idv(ts_->polyDeriv()),id( v_ ) ) );
        }

        auto myelts = getContactRegion(u_);

        if (nbrFaces_ > 0)
        {

            std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
            processContact(lt_,at_, myelts, u_, v_);
        }

        processBoundaryConditions(lt_, at_);

        at_.solve( _rhs = lt_, _solution = u_ );

        
        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

        // Reset
        at_.zero();
        lt_.zero();

        it += 1;
    }
}

// Time loop with fixed point iterations
template <int Dim, int Order>
void ElasticContact<Dim, Order>::timeLoopfixedPoint()
{
    int it = 0; 
    processLoading(l_);
    processMaterials(a_);

    std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] start time stepping start: {}, stop: {}, step: {}", 
                                fmt::localtime(std::time(nullptr)), ts_->timeInitial(),ts_->timeFinal(), ts_->timeStep()) << std::endl;
    
    for ( ts_->start(); ts_->isFinished()==false; ts_->next(u_) )
    {
        if (Environment::isMasterRank())
            std::cout << fmt::format( "[{:%Y-%m-%d :%H:%M:%S}] time {:.6f}/{}", fmt::localtime(std::time(nullptr)), ts_->time(),ts_->timeFinal()) << std::endl;

        ////////////////////////////////////////////////////
        //          Newmark beta-model for dttun          //
        ////////////////////////////////////////////////////
        lt_ = l_;
        at_ = a_;

        for ( auto [key, material] : specs_["/Models/LinearElasticity/Materials"_json_pointer].items() )
        {
            LOG( INFO ) << fmt::format( "Material {} found", material );

            std::string matRho = fmt::format( "/Materials/{}/parameters/rho/value", material.get<std::string>() );
            auto rho = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());
            lt_ +=  integrate( _range=markedelements( mesh_, material.get<std::string>() ), _expr= rho*inner( idv(ts_->polyDeriv()),id( v_ ) ) );
        }

        auto u_tmp =  project(_space=Xh_, _range=elements(mesh_), _expr = idv(u_));
        auto u_tmpNew = project(_space=Xh_, _range=elements(mesh_), _expr = idv(u_));

        int fixedPointIteration = 0;
        double fixedPointerror = 0.;
        while ((fixedPointerror > fixedPointtolerance_) || (fixedPointIteration < 1))
        {

            lt_tmp = lt_;
            at_tmp = at_;

            u_tmp = project(_space=Xh_, _range=elements(mesh_), _expr = idv(u_tmpNew)); ;

            auto myelts = getContactRegion(u_tmp);
            if (nbrFaces_ > 0)
            {
                std::cout << "Fixed point iteration : " << fixedPointIteration << std::endl;
                std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
                std::cout << "Error : " << fixedPointerror << std::endl;

                processContact(lt_tmp,at_tmp, myelts, u_tmp, u_tmp);
                at_tmp.solve(_rhs = lt_tmp, _solution = u_tmpNew);

                fixedPointerror = integrate(_range=elements(mesh_), _expr = norm2( idv(u_tmp)-idv(u_tmpNew))).evaluate()(0,0) / integrate(_range=elements(mesh_),_expr=norm2(idv(u_))).evaluate()(0,0); 
                fixedPointIteration++;
            }

            else if (fixedPointIteration == 0)
                break;
            
            if (fixedPointIteration == 10)
                break;
            
            
            lt_tmp.zero();
            at_tmp.zero();
        }
        
        auto myelts = getContactRegion(u_tmpNew);
        if (nbrFaces_ > 0)
        {
            std::cout << "Faces in contact : " << nbrFaces_ << std::endl;
            processContact(lt_,at_, myelts, u_tmpNew, u_tmpNew);
        }
        
        //processBoundaryConditions(lt_, at_);

        at_.solve( _rhs = lt_, _solution = u_ );
        
        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

        // Reset
        at_.zero();
        lt_.zero();

        it += 1;
    }
}


// Export results
template <int Dim, int Order>
void ElasticContact<Dim, Order>::exportResults(double t)
{
    e_->step(t)->addRegions();
    e_->step(t)->add( "displacement", u_ );
    e_->step(t)->add( "velocity", ts_->currentVelocity() );
    e_->step(t)->add( "acceleration", ts_->currentAcceleration() );
    e_->step(t)->add( "contactRegion", contactRegion_) ;

    if (nbrFaces_ > 0)
    {
        auto const Id = eye<Dim,Dim>();
        auto defv = sym(gradv(u_));
        auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
        contactPressure_ =  project(_space=XhC_, _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav);
        contactDisplacement_ = project(_space=XhC_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_));
    }


    e_->step(t)->add( "contactPressure", contactPressure_);
    e_->step(t)->add( "contactDisplacement", contactDisplacement_ );
    e_->step(t)->add( "contactFaces", contactFaces_ );
    e_->save();

    meas_["time"].push_back(ts_->time());

    

    auto compute = [this]( auto const& u, std::string const& name )
    {
        auto totalQuantity = integrate(_range=elements(mesh_), _expr=idv(u)).evaluate();
        auto totalFlux = integrate(_range=boundaryfaces(mesh_), _expr=gradv(u)*N()).evaluate();
        double meas=measure(_range=elements(mesh_), _expr=cst(1.0));
        
        if constexpr ( space_t::is_scalar )
        {
            meas_[fmt::format("{}_totalQuantity",name)].push_back(totalQuantity(0,0));
            meas_[fmt::format("{}_totalFlux",name)].push_back(totalFlux(0,0));
            meas_[fmt::format("{}_mean",name)].push_back(totalQuantity(0,0)/meas);
            meas_[fmt::format("{}_min",name)].push_back(u.min());
            meas_[fmt::format("{}_max",name)].push_back(u.max());
        }
        else
        {
            for (int i = 0; i < Dim; i++)
            {
                meas_[fmt::format("{}_totalQuantity_{}",name,i)].push_back(totalQuantity(i,0));
                meas_[fmt::format("{}_totalFlux_{}",name,i)].push_back(totalFlux(i,0));
                meas_[fmt::format("{}_mean_{}",name,i)].push_back(totalQuantity(i,0)/meas);
                //meas_[fmt::format("min_{}",i)].push_back(u.min()[i]);
                //meas_[fmt::format("max_{}",i)].push_back(u.max()[i]);
            }
        }
        for( auto [key,values] : mesh_->markerNames())
        {
            if ( values[1] == Dim )
            {
                double meas=measure(_range=markedelements(mesh_,key), _expr=cst(1.0));
                auto quantity = integrate(_range=markedelements(mesh_,key), _expr=idv(u)).evaluate();
                if constexpr ( space_t::is_scalar )
                {
                    meas_[fmt::format("{}_quantity_{}",name,key)].push_back(quantity(0,0));
                    meas_[fmt::format("{}_mean_{}",name,key)].push_back(quantity(0,0)/meas);
                }
                else
                {
                    for (int i = 0; i < Dim; i++)
                    {
                        meas_[fmt::format("{}_quantity_{}_{}",name,key,i)].push_back(quantity(i,0));
                        meas_[fmt::format("{}_mean_{}_{}",name,key,i)].push_back(quantity(i,0)/meas);
                    }
                
                }
            }
            else if ( values[1] == Dim-1 )
            {
                double meas=measure(_range=markedfaces(mesh_,key), _expr=cst(1.0));
                auto quantity = integrate(_range=markedfaces(mesh_,key), _expr=idv(u)).evaluate();
                auto flux = integrate(_range=markedfaces(mesh_,key), _expr=gradv(u)*N()).evaluate();

                if constexpr ( space_t::is_scalar )
                {
                    meas_[fmt::format("{}_quantity_{}",name,key)].push_back(quantity(0,0));
                    meas_[fmt::format("{}_mean_{}",name,key)].push_back(quantity/meas);
                    meas_[fmt::format("{}_flux_{}",name,key)].push_back(flux);
                }
                else
                {
                    for (int i = 0; i < Dim; i++)
                    {
                        meas_[fmt::format("{}_quantity_{}_{}",name,key,i)].push_back(quantity(i,0));
                        meas_[fmt::format("{}_mean_{}_{}",name,key,i)].push_back(quantity(i,0)/meas);
                        meas_[fmt::format("{}_flux_{}_{}",name,key,i)].push_back(flux(i,0));
                    }
                }
            }

        }
    };
    
    compute(u_,"displacement");
    compute(ts_->currentVelocity(),"velocity");
    compute(ts_->currentAcceleration(),"acceleration");

    auto Id = eye<Dim,Dim>();
    auto epsv = sym(gradv(u_));
    auto sig = lambda_*trace(epsv)*Id + 2*mu_*epsv;
    auto sigmav = (cst(lambda_)*trace(epsv)*Id + cst(2*mu_)*epsv)*N();

    double E1 = 0.5*normL2Squared(_range=elements(mesh_),_expr=idv(ts_->currentVelocity()));
    double E2 = 0.5*integrate( _range= elements( mesh_ ), _expr= inner(sig,epsv) ).evaluate()( 0,0 );
    meas_["Eh1"].push_back(E1);
    meas_["Eh2"].push_back(E2);
    meas_["Eh"].push_back(E1 + E2);

    
    double R1 = normL2Squared(_range= boundaryfaces(mesh_), _expr= sqrt(h()) * trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces_));
    double R2 = normL2Squared( _range= boundaryfaces(mesh_), _expr= sqrt(h()) * ( cst(gamma_) * ( trans(expr<Dim,1>(direction_)) *idv(u_)  - abs(cst(distance_) - Py()) ) - trans(expr<Dim,1>(direction_))*sigmav ) * idv(contactFaces_));

    double disp = integrate( _range= boundaryfaces(mesh_), _expr= (trans(expr<Dim,1>(direction_))*idv(u_)  - abs(cst(distance_) - Py())) * idv(contactFaces_ )).evaluate()( 0,0 );

    meas_["R1"].push_back(R1);
    meas_["R2"].push_back(R2);
    meas_["R"].push_back((R1 - R2)/(2.*gamma0_));
    meas_["disp"].push_back(disp);
    double Es = (E1 + E2) - theta_*(R1 - R2)/(2.*gamma0_);
    meas_["Es"].push_back(Es);      
    
    double Lv = integrate( _range=elements(mesh_),_expr=  trans(expr<Dim,1>( externalforce_ ))*idv(u_)).evaluate()( 0,0 );
    meas_["Lv"].push_back(Lv);
    meas_["E"].push_back(Es - Lv);

    this->writeResultsToFile("measures.json");
    
}

template <int Dim, int Order>
void ElasticContact<Dim, Order>::writeResultsToFile(const std::string& filename) const
{
    if ( Environment::isMasterRank() )
    {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << meas_.dump(4);  // Indent of 4 spaces for readability
            file.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }
}



} 