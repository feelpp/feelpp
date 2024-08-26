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
class BilateralContact
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
    using ts_ptrtype = std::shared_ptr<Newmark<spacev_t>>;
    using exporter_ptrtype = std::shared_ptr<Exporter<mesh_t>>; 

    // Constructors
    BilateralContact() = default;
    BilateralContact(nl::json const& specs);

    // Accessors
    nl::json const& specs() const { return specs_; }
    std::shared_ptr<mesh_t> const& mesh() const { return mesh_; }
    spacev_ptr_t const& Xhv() const { return Xhv_; }
    space_ptr_t const Xh() const { return Xh_; } 
    elementv_t const& u1() const { return u1_; }
    elementv_t const& u2() const { return u2_; }
    exporter_ptrtype const& exporter() const { return e_; }

    // Mutators
    void setSpecs(nl::json const& specs) { specs_ = specs; }
    void setMesh(std::shared_ptr<mesh_t> const& mesh) { mesh_ = mesh; }
    void setU1(elementv_t const& u1) { u1_ = u1; }
    void setU2(elementv_t const& u2) { u2_ = u2; }

    // Methods
    void initialize();
    void processLoading(form1_type& l);
    void processMaterials(form2_type &a);
    void processBoundaryConditions(form1_type& l, form2_type& a);
    void run();
    void exportResults();


private:
    nl::json specs_;
    std::shared_ptr<mesh_t> mesh_;
    spacev_ptr_t Xhv_;
    space_ptr_t Xh_;

    elementv_t u1_;
    elementv_t u2_;
    exporter_ptrtype e_;


    double H_;
    double E1_, nu1_, lambda1_, mu1_, rho1_;
    double E2_, nu2_, lambda2_, mu2_, rho2_;
    std::string externalforce_;

};

// Constructor
template <int Dim, int Order>
BilateralContact<Dim, Order>::BilateralContact(nl::json const& specs) : specs_(specs)
{
}

// Initialization 
template <int Dim, int Order>
void BilateralContact<Dim, Order>::initialize()
{
    // Get mesh size
    H_ = specs_["/Meshes/LinearElasticity/Import/h"_json_pointer].get<double>();
    // Load mesh
    mesh_ = loadMesh( _mesh = new mesh_t, _filename = specs_["/Meshes/LinearElasticity/Import/filename"_json_pointer].get<std::string>(), _h = H_);
    // Define Xhv
    Xhv_ = Pchv<Order>(mesh_); 

    // Get elastic structure parameters
    std::string matRho1 = fmt::format( "/Materials/Solid1/parameters/rho/value");
    rho1_ = std::stod(specs_[nl::json::json_pointer( matRho1 )].get<std::string>());

    std::string matE1 = fmt::format( "/Materials/Solid1/parameters/E/value" );
    double E1_ = std::stod(specs_[nl::json::json_pointer( matE1 )].get<std::string>());
    
    std::string matNu1 = fmt::format( "/Materials/Solid1/parameters/nu/value" );
    double nu1_ = std::stod(specs_[nl::json::json_pointer( matNu1 )].get<std::string>());
    
    lambda1_ = E1_*nu1_/( (1+nu1_)*(1-2*nu1_) );
    mu1_ = E1_/(2*(1+nu1_));

    std::string matRho2 = fmt::format( "/Materials/Solid1/parameters/rho/value");
    rho2_ = std::stod(specs_[nl::json::json_pointer( matRho2 )].get<std::string>());

    std::string matE2 = fmt::format( "/Materials/Solid1/parameters/E/value" );
    double E2_ = std::stod(specs_[nl::json::json_pointer( matE2 )].get<std::string>());
    
    std::string matNu2 = fmt::format( "/Materials/Solid1/parameters/nu/value" );
    double nu2_ = std::stod(specs_[nl::json::json_pointer( matNu2 )].get<std::string>());
    
    lambda2_ = E2_*nu2_/( (1+nu2_)*(1-2*nu2_) );
    mu2_ = E2_/(2*(1+nu2_));


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



// Process loading
template <int Dim, int Order>
void BilateralContact<Dim, Order>::processLoading(form1_type& l)
{
    l += integrate( _range = elements(mesh_), _expr = cst(rho_)*trans(expr<Dim,1>( externalforce_ ))*id(u_) );
}

// Process materials
template <int Dim, int Order>
void BilateralContact<Dim, Order>::processMaterials( form2_type &a )
{
    auto deft = sym(gradt(u_));
    auto def = sym(grad(u_));
    auto Id = eye<Dim,Dim>();
    auto sigmat = lambda_*trace(deft)*Id + 2*mu_*deft;

    a += integrate( _range = elements(mesh_), _expr = inner(sigmat,def));
}






// Process boundary conditions
template <int Dim, int Order>
void BilateralContact<Dim, Order>::processBoundaryConditions(form1_type& l, form2_type& a)
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
void BilateralContact<Dim, Order>::run()
{
    std::cout << "***** Run static elasticity with unilateral contact *****" << std::endl;

    std::cout << "***** Initialize elasticity parameters *****" << std::endl;
    initialize();

    u1_ = Xhv_->element();
    u2_ = Xhv_->element();

    this->exportResults(0);

    auto a_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto at_ = form2( _test = Xhv_, _trial = Xhv_ );
    auto l_ = form1( _test = Xhv_ );
    auto lt_ = form1( _test = Xhv_ );

    a_.zero();
    at_.zero();
    l_.zero();
    lt_.zero();

    
    std::cout << "***** Process loading *****" << std::endl;
    processLoading(l_);

    std::cout << "***** Process materials *****" << std::endl;
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
            lt_ +=  integrate( _range=elements( mesh_), _expr= cst(rho_)*inner( idv(ts_->polyDeriv()),id( u_ ) ) );
        
        
        std::cout << "***** Process boundary conditions *****" << std::endl;
        processBoundaryConditions(lt_, at_);

        std::cout << "***** Solve *****" << std::endl;
        at_.solve( _rhs = lt_, _solution = u_ );

        ts_->updateFromDisp(u_);
        this->exportResults(ts_->time());

        // Reset
        at_.zero();
        lt_.zero();

    }  

    

}








// Export results
template <int Dim, int Order>
void 
BilateralContact<Dim, Order>::exportResults()
{
    // Define exports
    e_->addRegions();
    e_->add( "displacement", u_ );
    e_->save();

    
}

} 