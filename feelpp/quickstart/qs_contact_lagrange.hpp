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
    /*
    auto face_mesh = createSubmesh( _mesh=mesh_, _range=boundaryfaces(mesh_ ), _update=0 );
    auto XhCFaces = Pdh<0>(face_mesh);
    auto contactFaces = XhCFaces->element();
    contactFaces.on(_range=myelts, _expr = cst(1.));

    auto contactRegion = Xh_P1->element();
    contactRegion.on( _range=boundaryfaces(mesh_), _expr = idv(contactFaces));
    
    auto defvinter = sym(gradv(uinter));
    auto sigmavinter = (lambda_*trace(defvinter)*Id + 2*mu_*defvinter)*N();
   
    auto contactPressureinter =  Xh_P1->element();
    contactPressureinter.on( _range=boundaryfaces(mesh_), _expr = trans(expr<Dim,1>(direction_))*sigmavinter);

    auto contactDisplacementinter =Xh_P1->element();
    contactDisplacementinter.on( _range=elements(mesh_), _expr = (trans(expr<Dim,1>(direction_))*idv(uinter) - idv(g_)));

    e_->add( "contactRegion", contactRegion);
    e_->add( "contactPressure", contactPressureinter);
    e_->add( "contactDisplacement", contactDisplacementinter );
    
    */


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