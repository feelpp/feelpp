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
class ContactStatic
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
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

    // Constructors
    ContactStatic() = default;
    ContactStatic(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    elementv_t const& u() const { return u_; }
    exporter_ptrtype const& exporter() const { return e_; }
    nl::json measures() const { return meas_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU(elementv_t const& u) { u_ = u; }

    // Methods
    void initialize();
    void initializeContact();
    void processLoading(form1_type& l);
    void processMaterials(form2_type &a);
    void processBoundaryConditions(form1_type& l, form2_type& a);
    void processContactPenalty(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts, elementv_t const& u);
    void processContactNitsche(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts, elementv_t const& u);
    void run();
    typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype getContactRegion(elementv_t const& u);
    void exportResults();
    void initG();
    void storeData();
    void error();

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
    typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype myelts_;

    exporter_ptrtype e_;
    nl::json meas_;

    double H_,E_, nu_, lambda_, mu_, rho_;
    std::string externalforce_;

    double epsilon_,tolContactRegion_,tolDistance_;
    double theta_, gamma0_, gamma_;
    std::string method_, direction_;
    std::vector<double> ddirection_;
    int withMarker_,save_, computeerror_;
};

// Constructor
template <int Dim, int Order>
ContactStatic<Dim, Order>::ContactStatic(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order>
void ContactStatic<Dim, Order>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // Define Xhv
    Xhv_ = Pchv<Order>(mesh_); 

    // Get elastic structure parameters
    std::string matRho = fmt::format( "/Materials/Caoutchouc/parameters/rho/value");
    rho_ = std::stod(specs_[nl::json::json_pointer( matRho )].get<std::string>());

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
    
}

// Initialization of the contact terms
template <int Dim, int Order>
void ContactStatic<Dim, Order>::initializeContact()
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
}

// Process loading
template <int Dim, int Order>
void ContactStatic<Dim, Order>::processLoading(form1_type& l)
{
    l += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_) );
}

// Process materials
template <int Dim, int Order>
void ContactStatic<Dim, Order>::processMaterials( form2_type &a )
{
    auto deft = sym(gradt(u_));
    auto def = sym(grad(u_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a += integrate( _range = elements(mesh_), _expr = inner(sigmat,def));
}

// Process contact conditions Penalty method
template <int Dim, int Order>
void ContactStatic<Dim, Order>::processContactPenalty(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts , elementv_t const& u )
{    
    if (withMarker_ == 0)
    {
        auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), elts->begin(), elts->end(), elts );
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces = project(_space=XhCFaces, _range=myfaces, _expr = cst(1.));

        a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));
        l += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)) * idv(contactFaces));     
    }
    else 
    {
        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(trans(expr<Dim,1>(direction_))*idt(u),trans(expr<Dim,1>(direction_))*id(u)));
        l += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(epsilon_) * inner(idv(g_), trans(expr<Dim,1>(direction_))*id(u)));   
    }
}

// Process contact conditions Nitsche method
template <int Dim, int Order>
void ContactStatic<Dim, Order>::processContactNitsche(form1_type& l, form2_type& a, typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype elts , elementv_t const& u )
{    
    if (withMarker_ == 0)
    {
        auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), elts->begin(), elts->end(), elts );
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = XhCFaces->element();
        contactFaces = project(_space=XhCFaces, _range=myfaces, _expr = cst(1.));

        auto const Id = eye<Dim,Dim>();
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
        auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

        a += integrate (_range=boundaryfaces(mesh_),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces)); 
        a += integrate (_range=boundaryfaces(mesh_),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces));
        l += integrate (_range=boundaryfaces(mesh_),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma)*idv(contactFaces)); 
    
    }
    else 
    {
        auto const Id = eye<Dim,Dim>();
        auto deft = sym(gradt(u));
        auto def = sym(grad(u));
        auto sigma = (lambda_*trace(def)*Id + 2*mu_*def)*N();
        auto sigmat = (lambda_*trace(deft)*Id + 2*mu_*deft)*N();

        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= - cst(theta_)/cst(gamma_) * inner(trans(expr<Dim,1>(direction_))*sigmat, trans(expr<Dim,1>(direction_))*sigma)); 
        a += integrate (_range=markedfaces(mesh_,"contact"),_expr= cst(1.)/cst(gamma_) * inner(cst(gamma_) * trans(expr<Dim,1>(direction_))*idt(u) - trans(expr<Dim,1>(direction_))*sigmat, cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));
        l += integrate (_range=markedfaces(mesh_,"contact"),_expr= inner(idv(g_), cst(gamma_) * trans(expr<Dim,1>(direction_))*id(u) - cst(theta_)*trans(expr<Dim,1>(direction_))*sigma));     
    }
    
}


// Process boundary conditions
template <int Dim, int Order>
void ContactStatic<Dim, Order>::processBoundaryConditions(form1_type& l, form2_type& a)
{
    // Boundary Condition Dirichlet
    if ( specs_["/BoundaryConditions/LinearElasticity"_json_pointer].contains("Dirichlet") )
    {
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
void ContactStatic<Dim, Order>::run()
{
    std::cout << "***** Run static elasticity with unilateral contact *****" << std::endl;

    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    std::cout << "***** Initialize contact parameters *****" << std::endl;
    initializeContact();

    std::cout <<  "***** Initialize distance g *****" << std::endl;
    initG();

    if ((method_.compare("penalty") == 0) || (method_.compare("nitsche") == 0))
    {
        // Initialize elements
        u_ = Xhv_->element();

        // Initialize linear and bilinear forms
        auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
        auto l_ = form1( _test = Xhv_ );

        l_.zero();
        a_.zero();

        std::cout << "***** Process loading *****" << std::endl;
        processLoading(l_);

        std::cout << "***** Process materials *****" << std::endl;
        processMaterials(a_);

        std::cout << "***** Process contact *****" << std::endl;
        if (withMarker_ == 0)
        {
            myelts_ = getContactRegion(u_);
            std::cout << "Nbr faces for processContact : " << nbrFaces_ << std::endl;
        }
        else 
        {
            contactRegion_ = project(_space=Xh_, _range=markedfaces(mesh_,"contact"), _expr = cst(1.));
        }

        if (method_.compare("penalty") == 0)
            processContactPenalty(l_, a_, myelts_, u_);
        else if (method_.compare("nitsche") == 0)
            processContactNitsche(l_, a_, myelts_, u_);
        
        std::cout << "***** Process boundary conditions *****" << std::endl;
        processBoundaryConditions(l_, a_);

        std::cout << "***** Solve *****" << std::endl;
        a_.solve( _rhs = l_, _solution = u_ );

        std::cout << "***** Export results *****" << std::endl;
        exportResults();

        if (save_ == 1)
        {
            std::cout << "***** Save results *****" << std::endl;
            storeData();
        }

        if (computeerror_ == 1)
        {
            std::cout << "***** Compute error *****" << std::endl;
            error();
        }
    
    }
    else if (method_.compare("al") == 0)
    {
        if (Order != 2)
            std::cout << "u has to be P2" << std::endl;


        // Define lagrange multiplier
        auto submesh = createSubmesh(_mesh=mesh_,_range=markedfaces(mesh_,"contact"));
        auto XhLambda = Pch<1>(submesh);
        
        // Initialize elements
        auto u = Xhv_->elementPtr();
        auto lambda = XhLambda->elementPtr();

        // Initialize linear and bilinear forms
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
        auto displacement = project(_space=Xhv_, _range=elements(mesh_), _expr=idv(*u));
        e_->add( "displacement", displacement );

        auto realcontactRegion = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(*u) - idv(g_));  
        e_->add( "realcontactRegion", realcontactRegion) ;


        int nbrFacesC = 0;
        typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_type );

        auto const& trialDofIdToContainerId =  form2(_test=Xh_, _trial=Xh_).dofIdToContainerIdTest();
        for (auto const& theface : boundaryfaces(mesh_) )
        {                
            auto & face = boost::unwrap_ref( theface );
            int contactDof = 0;
            for( auto const& ldof : Xh_->dof()->faceLocalDof( face.id() ) )
            {
                index_type thedof = ldof.index();
                thedof = trialDofIdToContainerId[ thedof ];

                if (realcontactRegion[thedof] >= tolContactRegion_)
                    contactDof++;
                    
                if (Dim == 2)
                {
                    if (contactDof == 2)
                        nbrFacesC++;
                      
                }      
                else if (Dim == 3)
                {
                    if (contactDof == 4)
                        nbrFacesC++;
                }
            }
        }
        std::cout << "Faces in contact export : " << nbrFacesC << std::endl;
        myelts->shrink_to_fit();  
        
        auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
        auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
        auto XhCFaces = Pdh<0>(face_mesh);
        auto contactFaces = project(_space=XhCFaces, _range=myfaces, _expr = cst(1.));

        auto defv = sym(gradv(u));
        auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
        contactRegion_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = idv(contactFaces));
        contactPressure_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
        contactDisplacement_ = project(_space=Xh_, _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(*u) - idv(g_))*idv(contactFaces));

        auto expression = project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = (cst(gamma_)*(trans(expr<Dim,1>(direction_))*idv(*u) - idv(g_)) - trans(expr<Dim,1>(direction_))*sigmav)*idv(contactFaces)  );

        e_->add( "contactRegion", contactRegion_) ;
        e_->add( "contactPressure", contactPressure_);
        e_->add( "contactDisplacement", contactDisplacement_ );
        e_->add( "expression", expression);
        e_->add("g", g_);
        e_->save();

        int nbr = 1;
    
        for (auto &bfaceC : markedfaces(mesh_,"contact"))
        {
            auto & faceC = boost::unwrap_ref( bfaceC );

            auto ctx = Xh_->context();
            node_type t1(Dim);
            t1(0)=faceC.point(0).node()[0]; t1(1)=faceC.point(0).node()[1];
            ctx.add( t1 );

            auto evaluateStress = evaluateFromContext( _context=ctx, _expr = idv(contactPressure_) ); 
            auto evaluateDisp = evaluateFromContext( _context=ctx, _expr = idv(contactDisplacement_) );     
        
            std::cout << "Evaluate Stress : " << evaluateStress(0,0) << " for node : " <<  nbr << " with coordinates : (" << faceC.point(0).node()[0] << "," << faceC.point(0).node()[1]  << ")" << std::endl;
            std::cout << "Evaluate Disp : " << evaluateDisp(0,0) << " for node : " <<  nbr << " with coordinates : (" << faceC.point(0).node()[0] << "," << faceC.point(0).node()[1]  << ")" << std::endl;
        
            nbr++;   
        } 

            
        if (save_ == 1)
        {
            std::cout << "***** Save results *****" << std::endl;
            (*u).saveHDF5("solution_ref.h5");
            saveGMSHMesh(_mesh=mesh_,_filename= "mesh_ref.msh" );
        }

        if (computeerror_ == 1)
        {
            std::cout << "***** Compute error *****" << std::endl;
            auto meshref = loadMesh(_mesh=new mesh_t, _filename = "/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/mesh_ref.msh");
            auto uref = spacev_t::New(meshref)->elementPtr() ;
            uref->loadHDF5("/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/solution_ref.h5");

            auto op_inter = opInterpolation(_domainSpace =  spacev_t::New(mesh_), _imageSpace = spacev_t::New(meshref) );
            auto uinter =  spacev_t::New(meshref)->element(); 
            op_inter->apply(*u, uinter);

            auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
            std::cout << "H1 relative error : " << h1err << std::endl;

        }    
    }

}

template<int Dim, int Order>
typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype
ContactStatic<Dim, Order>::getContactRegion(elementv_t const& u)
{   
    typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<Mesh<Simplex<Dim>>>::faces_reference_wrapper_type );

    auto defv = sym(gradv(u));
    auto Id = eye<Dim,Dim>();
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
    
    contactRegion_ = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u) - idv(g_));    

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
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 3)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
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
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (Dim == 3)
                {
                    if (contactDof == 4)
                    {
                        nbrFaces_++;
                        myelts->push_back( boost::cref( face ) );
                    }
                } 
            }
        }
    }
    myelts->shrink_to_fit();   
  
    return myelts;
}

template <int Dim, int Order>
void 
ContactStatic<Dim, Order>::initG()
{   
    // Init the distance fields 
    g_ = Xh_->element();

    // Load and export rigid obstacles
    auto wall = loadMesh(_mesh=new mesh_t, _filename = "$cfgdir/wall.geo");

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
            auto rayIntersection = bvh->intersect(ray) ;

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


template <int Dim, int Order>
void 
ContactStatic<Dim, Order>::storeData()
{
    u_.saveHDF5("solution_ref.h5");
    saveGMSHMesh(_mesh=mesh_,_filename= "mesh_ref.msh" );
}


template <int Dim, int Order>
void 
ContactStatic<Dim, Order>::error()
{
    auto meshref = loadMesh(_mesh=new mesh_t, _filename = "/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/mesh_ref.msh");
    auto uref = spacev_t::New(meshref)->elementPtr() ;
    uref->loadHDF5("/data/scratch/vanlandeghem/feel/qs_elasticity_contact/benchmark/np_1/solution_ref.h5");

    auto op_inter = opInterpolation(_domainSpace =  spacev_t::New(mesh_), _imageSpace = spacev_t::New(meshref) );
    auto uinter =  spacev_t::New(meshref)->element(); 
    op_inter->apply(u_, uinter);

    auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
    std::cout << "H1 error : " << h1err << std::endl;
}

// Export results
template <int Dim, int Order>
void 
ContactStatic<Dim, Order>::exportResults()
{

    // Define exports
    e_->addRegions();
    e_->add( "displacement", u_ );

    // Compute new contact region 
    myelts_ = getContactRegion(u_);
    std::cout << "Nbr faces for export : " << nbrFaces_ << std::endl;
    
    auto realcontactRegion = project(_space=Xh_, _range=elements(mesh_), _expr = trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_));  
    e_->add( "realcontactRegion", realcontactRegion) ;

    auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts_->begin(), myelts_->end(), myelts_ );
    auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
    auto XhCFaces = Pdh<0>(face_mesh);
    auto contactFaces = project(_space=XhCFaces, _range=myfaces, _expr = cst(1.));

    auto const Id = eye<Dim,Dim>();
    auto defv = sym(gradv(u_));
    auto sigmav = (lambda_*trace(defv)*Id + 2*mu_*defv)*N();
   
    contactRegion_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = idv(contactFaces));
    contactPressure_ =  project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmav*idv(contactFaces));
    contactDisplacement_ = project(_space=Xh_, _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_))*idv(contactFaces));

    auto expression = Xh_->element();
    if (method_.compare("penalty") == 0)
        expression = project(_space=Xh_, _range=elements(mesh_), _expr =  (trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_))*idv(contactFaces)  );
    else if (method_.compare("nitsche") == 0)
        expression = project(_space=Xh_, _range=boundaryfaces(mesh_), _expr = (cst(gamma_)*(trans(expr<Dim,1>(direction_))*idv(u_) - idv(g_)) - trans(expr<Dim,1>(direction_))*sigmav)*idv(contactFaces)  );
        
    e_->add( "contactRegion", contactRegion_) ;
    e_->add( "contactPressure", contactPressure_);
    e_->add( "contactDisplacement", contactDisplacement_ );
    e_->add( "expression", expression);
    e_->add("g", g_);
    e_->save();

    // Print values on contact boundary

    int nbr = 1;
    for (auto &bfaceC : markedfaces(mesh_,"contact"))
    {
        auto & faceC = boost::unwrap_ref( bfaceC );

        auto ctx = Xh_->context();
        node_type t1(Dim);
        t1(0)=faceC.point(0).node()[0]; t1(1)=faceC.point(0).node()[1];
        ctx.add( t1 );

        auto evaluateStress = evaluateFromContext( _context=ctx, _expr = idv(contactPressure_) ); 
        auto evaluateDisp = evaluateFromContext( _context=ctx, _expr = idv(contactDisplacement_) );     
        
        std::cout << "Evaluate Stress : " << evaluateStress(0,0) << " for node : " <<  nbr << " with coordinates : (" << faceC.point(0).node()[0] << "," << faceC.point(0).node()[1]  << ")" << std::endl;
        std::cout << "Evaluate Disp : " << evaluateDisp(0,0) << " for node : " <<  nbr << " with coordinates : (" << faceC.point(0).node()[0] << "," << faceC.point(0).node()[1]  << ")" << std::endl;
        
        nbr++;
    }
}

} 