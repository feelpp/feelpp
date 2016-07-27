#ifndef _MIXEDELASTICITY_HPP
#define _MIXEDELASTICITY_HPP 

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/stencil.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/complement.hpp>
#include <feel/feelalg/topetsc.hpp>
// #include <feel/feelts/bdf.hpp>
#include <boost/algorithm/string.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelpoly/raviartthomas.hpp>
#include <feel/feeldiscr/pdhm.hpp>

namespace Feel {

namespace FeelModels {

inline
po::options_description
makeMixedElasticityOptions( std::string prefix = "mixedelasticy" )
{
    po::options_description mpOptions( "Mixed Elasticity HDG options");
    mpOptions.add_options()
        ( prefixvm( prefix, "gmsh.submesh").c_str(), po::value<std::string>()->default_value( "" ), "submesh extraction" )
        ( prefixvm( prefix, "hface").c_str(), po::value<int>()->default_value( 0 ), "hface" )
        ( "lambda", po::value<std::string>()->default_value( "1" ), "lambda" )
        ( "mu", po::value<std::string>()->default_value( "1" ), "mu" )
        ( prefixvm( prefix, "model_json").c_str(), po::value<std::string>()->default_value("model.json"), "json file for the model")
        ( prefixvm( prefix,"tau_constant").c_str(), po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( prefixvm( prefix,"tau_order").c_str(), po::value<int>()->default_value( -1 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ;
    mpOptions.add ( envfeelmodels_options( prefix ) ).add( modelnumerical_options( prefix ) );
    return mpOptions;
}

inline po::options_description
makeMixedElasticityLibOptions( std::string prefix = "mixedelasticity" )
{
    po::options_description mpLibOptions( "Mixed Elasticity HDG Lib options");
    // if ( !prefix.empty() )
    //     mpLibOptions.add( backend_options( prefix ) );
    return mpLibOptions;
}

template<int Dim, int Order, int G_Order = 1>
class MixedElasticity    :	public ModelNumerical
{
public:
    typedef ModelNumerical super_type;

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    typedef MixedElasticity<Dim,Order,G_Order> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;
    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,G_Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // The Lagrange multiplier lives in R^n-1
    typedef Simplex<Dim-1,G_Order,Dim> face_convex_type;
    typedef Mesh<face_convex_type> face_mesh_type;
    typedef boost::shared_ptr<face_mesh_type> face_mesh_ptrtype;
    

// ---- //
    using Vh_t =  Pdhms_type<mesh_type,Order>;
    using Vh_ptr_t =  Pdhms_ptrtype<mesh_type,Order>;
    using Vh_element_t = typename Vh_t::element_type;
    using Vh_element_ptr_t = typename Vh_t::element_ptrtype;
// ---- //
    using Wh_t =  Pdhv_type<mesh_type,Order>;
    using Wh_ptr_t =  Pdhv_ptrtype<mesh_type,Order>;
    using Wh_element_t = typename Wh_t::element_type;
    using Wh_element_ptr_t = typename Wh_t::element_ptrtype;
// ---- //
    using Mh_t =  Pdhv_type<face_mesh_type,Order>;
    using Mh_ptr_t =  Pdhv_ptrtype<face_mesh_type,Order>;
    using Mh_element_t = typename Mh_t::element_type;
    using Mh_element_ptr_t = typename Mh_t::element_ptrtype;
// ---- //
    using M0h_t =  Pdh_type<face_mesh_type,0>;
    using M0h_ptr_t =  Pdh_ptrtype<face_mesh_type,0>;
    using M0h_element_t = typename M0h_t::element_type;
    using M0h_element_ptr_t = typename M0h_t::element_ptrtype;
    

 
    //! Model properties type
    using model_prop_type = ModelProperties;
    using model_prop_ptrtype = std::shared_ptr<model_prop_type>;

    using product_space_std = ProductSpaces<Vh_ptr_t,Wh_ptr_t,Mh_ptr_t>;
    using bilinear_block_std = BlockBilinearForm<product_space_std>;

    typedef Exporter<mesh_type,G_Order> exporter_type;
    typedef boost::shared_ptr <exporter_type> exporter_ptrtype;
    
    // typedef Bdf<space_mixedpoisson_type>  bdf_type;
    // typedef Bdf <Wh_t> bdf_type;
    // typedef Bdf<Pdh_type<mesh_type,Order>> bdf_type;
    // typedef boost::shared_ptr<bdf_type> bdf_ptrtype;
    
//private:
protected:
    model_prop_ptrtype M_modelProperties;
    std::string M_prefix;

    mesh_ptrtype M_mesh;

    Vh_ptr_t M_Vh; // flux
    Wh_ptr_t M_Wh; // potential
    Mh_ptr_t M_Mh; // potential trace

    backend_ptrtype M_backend;
    sparse_matrix_ptrtype M_A_cst;
    sparse_matrix_ptrtype M_A;
    vector_ptrtype M_F;
    vector_ptrtype M_U;

    Vh_element_ptr_t M_up; // flux solution
    Wh_element_ptr_t M_pp; // potential solution 

    double M_tau_constant;
    int M_tau_order;
    

public:
    
    // constructor
   MixedElasticity( std::string const& prefix = "mixedelasticity",                   
                    WorldComm const& _worldComm = Environment::worldComm(),
                    std::string const& subPrefix = "",
                    std::string const& rootRepository = ModelBase::rootRepositoryByDefault() );
    
    MixedElasticity( self_type const& ME ) = default;
    static self_ptrtype New( std::string const& prefix = "mixedelasticity",
                             WorldComm const& worldComm = Environment::worldComm(),
                             std::string const& subPrefix = "",
                             std::string const& rootRepository = ModelBase::rootRepositoryByDefault() ); 

    // Get Methods
    mesh_ptrtype mesh() const { return M_mesh; }
    Vh_ptr_t fluxSpace() const { return M_Vh; }
    Wh_ptr_t potentialSpace() const { return M_Wh; }
    Mh_ptr_t traceSpace() const { return M_Mh; }

    Vh_element_ptr_t fluxField() const { return M_up; }
    Wh_element_ptr_t potentialField() const { return M_pp; }
    model_prop_type modelProperties() { return *M_modelProperties; }
    model_prop_type modelProperties() const { return *M_modelProperties; }
    
    int tau_order() const { return M_tau_order; }
    backend_ptrtype get_backend() { return M_backend; }
    
    // Exporter
    // virtual void exportResults( mesh_ptrtype mesh = nullptr )
    //     {
    //         this->exportResults (this->currentTime(), mesh );
    //         M_exporter -> save();
    //     }
    // void exportResults ( double Time, mesh_ptrtype mesh = nullptr  ) ;
    // exporter_ptrtype M_exporter;
    // exporter_ptrtype exporterMP() { return M_exporter; }

    void init( mesh_ptrtype mesh = nullptr);

    virtual void initModel();
    virtual void initSpaces();
    // virtual void initExporter( mesh_ptrtype meshVisu = nullptr );
    virtual void assemble();
    
    template<typename PS>
    void assembleSTD(PS&& ps);
    template<typename PS>
    void assembleF(PS&& ps);
    /* 
    template<typename PS, typename ExprT>
    void updateConductivityTerm( PS&& ps, Expr<ExprT> expr, std::string marker = "");
    template<typename PS>
    void updateConductivityTerm( PS&& ps, bool isNL = false);
    template<typename PS, typename ExprT>
    void updatePotentialRHS( PS&& ps, Expr<ExprT> expr, std::string marker = "");
    template<typename PS, typename ExprT>
    void updateFluxRHS( PS&& ps, Expr<ExprT> expr, std::string marker = "");
    */
};

template<int Dim, int Order, int G_Order>
MixedElasticity<Dim, Order, G_Order>::MixedElasticity( std::string const& prefix,
                                                    WorldComm const& worldComm,
                                                    std::string const& subPrefix,
                                                    std::string const& rootRepository )
    : super_type( prefix, worldComm, subPrefix, rootRepository )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "start",
                                               this->worldComm(),this->verboseAllProc());


    this->setFilenameSaveInfo( prefixvm(this->prefix(),"MixedElasticity.info") );


    M_prefix = prefix;
    M_modelProperties = std::make_shared<model_prop_type>( Environment::expand( soption( prefixvm(M_prefix, "model_json") ) ) );
    if ( M_prefix.empty())
        M_backend = backend( _rebuild=true);
    else
        M_backend = backend( _name=M_prefix, _rebuild=true);

    M_tau_constant = doption (prefixvm(M_prefix, "tau_constant") );
    M_tau_order = ioption( prefixvm(M_prefix, "tau_order") );

    //-----------------------------------------------------------------------------//

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticityConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".MixedElasticitySolve.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MixedElasticity","constructor", "finish",
                                               this->worldComm(),this->verboseAllProc());
}


template<int Dim, int Order, int G_Order>
typename MixedElasticity<Dim,Order, G_Order>::self_ptrtype
MixedElasticity<Dim,Order,G_Order>::New( std::string const& prefix,
                                         WorldComm const& worldComm, std::string const& subPrefix,
                                         std::string const& rootRepository )
{
    return boost::make_shared<self_type> ( prefix,worldComm,subPrefix,rootRepository );
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::init( mesh_ptrtype mesh )
{
    tic();
    if ( !mesh )
        M_mesh = loadMesh( new mesh_type);
    else
        M_mesh = mesh;
    toc("mesh");

    tic();
    this->initModel();
    toc("model");

    tic();
    this->initSpaces();
    toc("spaces");

    // tic();
    // this->initExporter( );
    // toc("exporter");

    tic();
    this->assemble();
    toc("assemble");
}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::initModel()
{

    // initialize marker lists for each boundary condition type
    // Strain
    auto itField = M_modelProperties->boundaryConditions().find( "strain");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        
	if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
	 
        itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Neumann:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
        itType = mapField.find( "Robin" );
        if ( itType != mapField.end() )
        {
            Feel::cout << "Robin:";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
        }
    }

    // Displacement
    itField = M_modelProperties->boundaryConditions().find( "displacement");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        
	if ( itType != mapField.end() )
        {
            Feel::cout << "Dirichlet: ";
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                if ( M_mesh->hasFaceMarker(marker) )
                    Feel::cout << " " << marker;
                else
                    Feel::cout << std::endl << "WARNING!! marker " << marker << "does not exist!" << std::endl;
            }
            Feel::cout << std::endl;
	}
    }

}

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::initSpaces()
{


    auto face_mesh = createSubmesh( M_mesh, faces(M_mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );

    M_Vh = Pdhms<Order>( M_mesh, true );
    M_Wh = Pdhv<Order>( M_mesh, true );
    M_Mh = Pdhv<Order>( face_mesh, true );

    M_up = M_Vh->elementPtr( "u" ); // Strain
    M_pp = M_Wh->elementPtr( "p" ); // Displacement

    Feel::cout << "Vh<" << Order << "> : " << M_Vh->nDof() << std::endl
         << "Wh<" << Order << "> : " << M_Wh->nDof() << std::endl
         << "Mh<" << Order << "> : " << M_Mh->nDof() << std::endl;
}

// template<int Dim, int Order, int G_Order>
// void
// MixedPoisson<Dim, Order, G_Order>::initExporter( mesh_ptrtype meshVisu ) 
// {
//     std::string geoExportType="static"; //change_coords_only, change, static
//     M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh() ,
//                             _name="Export",
//                             _geo=geoExportType,
//         		    _path=this->exporterPath() ); 
// }

template<int Dim, int Order, int G_Order>
void
MixedElasticity<Dim, Order, G_Order>::assemble()
{

    auto ps = product(M_Vh,M_Wh,M_Mh);

    tic();
    M_A = M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    M_A_cst = M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    M_U = M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    M_F = M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    toc("creating matrices and vectors");
    
    tic();
    this->assembleSTD( ps );
    M_A_cst->close();
    // this->updateConductivityTerm( ps );
    M_A->close();
    toc("assemble A");
    
    tic();
    this->assembleF( ps );
    M_F->close();
    toc("assemble F");
    
    tic();
    M_backend->solve(_matrix=M_A,_rhs=M_F,_solution=M_U);
    toc("solve");
    
    // M_up = M_U[0];
    // M_pp = M_U[1];
}

template<int Dim, int Order, int G_Order>
template<typename PS>
void
MixedElasticity<Dim, Order, G_Order>::assembleSTD(PS&& ps)
{

    auto lambda = expr(soption("lambda")); 
    auto mu     = expr(soption("mu")); 
    auto c1     = cst(0.5)/mu; 
    auto c2     = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu)); 
    auto tau_constant = cst(M_tau_constant); 
    auto face_mesh = M_Mh->mesh();//createSubmesh( mesh, faces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 ); 
    M0h_ptr_t M0h = Pdh<0>( face_mesh,true ); 

    auto sigma = M_Vh->element( "sigma" ); 
    auto v     = M_Vh->element( "v" ); 
    auto u     = M_Wh->element( "u" ); 
    auto w     = M_Wh->element( "w" ); 
    auto uhat  = M_Mh->element( "uhat" ); 
    auto m     = M_Mh->element( "m" ); 
    auto H     = M0h->element( "H" ); 

    if ( ioption(prefixvm(M_prefix, "hface") ) == 0 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 1 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( ioption(prefixvm(M_prefix, "hface") ) == 2 )
	    H.on( _range=elements(M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
	    H.on( _range=elements(M0h->mesh()), _expr=h() );

    auto bbf = blockform2 ( ps, M_A_cst );


    bbf( 0_c, 0_c ) =  integrate(_range=elements(M_mesh),_expr=(c1*inner(idt(sigma),id(v))) );
    bbf( 0_c, 0_c ) += integrate(_range=elements(M_mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );

    cout << "a11 works fine" << std::endl;

    bbf( 0_c, 1_c ) += integrate(_range=elements(M_mesh),_expr=(trans(idt(u))*div(v)));

    cout << "a12 works fine" << std::endl;

    bbf( 0_c, 2_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
			    trans(idt(uhat))*rightface(id(v)*N())) );
    bbf( 0_c, 2_c) += integrate(_range=boundaryfaces(M_mesh),
		    _expr=-trans(idt(uhat))*(id(v)*N()));

    cout << "a13 works fine" << std::endl;

    bbf( 1_c, 0_c) += integrate(_range=elements(M_mesh),_expr=(trans(id(w))*divt(sigma)));

    cout << "a21 works fine" << std::endl;

    // begin dp: here we need to put the projection of u on the faces

    bbf( 1_c, 1_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=-tau_constant *
		    ( leftfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*leftface(id(w)) +
		      rightfacet( pow(idv(H),M_tau_order)*trans(idt(u)))*rightface(id(w) )));
    bbf( 1_c, 1_c) += integrate(_range=boundaryfaces(M_mesh),
		    _expr=-(tau_constant * pow(idv(H),M_tau_order)*trans(idt(u))*id(w)));

    cout << "a22 works fine" << std::endl;

    bbf( 1_c, 2_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=tau_constant *
		    ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tau_order)*id(w))+
		      rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tau_order)*id(w) )));

    bbf( 1_c, 2_c) += integrate(_range=boundaryfaces(M_mesh),
		    _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tau_order)*id(w) );

    cout << "a23 works fine" << std::endl;

    bbf( 2_c, 0_c) += integrate(_range=internalfaces(M_mesh),
		    _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
				    rightfacet(idt(sigma)*N())) ) );

    cout << "a31 works fine" << std::endl;

    // BC

    auto itField = M_modelProperties->boundaryConditions().find( "strain");
    if ( itField != M_modelProperties->boundaryConditions().end() )
    {
	    auto mapField = (*itField).second;
	    auto itType = mapField.find( "Dirichlet" );

	    if ( itType != mapField.end() )
	    {
		    Feel::cout << "Dirichlet: ";
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();
			    bbf( 2_c, 2_c) += integrate(_range=markedfaces(M_mesh,marker),
					    _expr=trans(idt(uhat)) * id(m) );
		    }
	    }

	    itType = mapField.find( "Neumann" );
	    if ( itType != mapField.end() )
	    {
		    Feel::cout << "Neumann:";
		    for ( auto const& exAtMarker : (*itType).second )
		    {
			    std::string marker = exAtMarker.marker();

			    bbf( 2_c, 0_c) += integrate(_range=markedfaces(M_mesh,marker ),
					    _expr=( trans(id(m))*(idt(sigma)*N()) ));

			    bbf( 2_c, 1_c) += integrate(_range=markedfaces(M_mesh,marker),
					    _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tau_order)*idt(u) ) );

			    bbf( 2_c, 2_c) += integrate(_range=markedfaces(M_mesh,marker), 
					    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tau_order) ) );

		    }
	    }
    }
    bbf( 2_c, 1_c) += integrate(_range=internalfaces(M_mesh),
			    _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tau_order)*idt(u) )+
				    rightfacet( pow(idv(H),M_tau_order)*idt(u) )));


    bbf( 2_c, 2_c) += integrate(_range=internalfaces(M_mesh),
    _expr=tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tau_order) )+
				    rightface( pow(idv(H),M_tau_order) )));


    
    cout << "BC works fine" << std::endl;


}


template<int Dim, int Order, int G_Order>
template< typename PS>
void
MixedElasticity<Dim, Order, G_Order>::assembleF(PS&& ps)
{
    M_F->zero();
    auto blf = blockform1( ps, M_F );

    auto w     = M_Wh->element( "w" ); 
    auto m     = M_Mh->element( "m" ); 
    
    // Building the RHS
    auto itField = M_modelProperties->boundaryConditions().find("strain");
    if (itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("SourceTerm");
        
	if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto g = expr<Dim,G_Order> (exAtMarker.expression());
		blf( 1_c ) += integrate(_range=elements(M_mesh),
                      	      _expr=trans(g)*id(w));
            }
        }
	itType = mapField.find("Neumann");
	if ( itType != mapField.end() )
	{
            for ( auto const& exAtMarker : (*itType).second )
            {
		auto marker = exAtMarker.marker();
                auto g = expr<2,2> (exAtMarker.expression());
		blf( 2_c ) += integrate(_range=markedfaces(M_mesh,marker),
				    _expr=trans(id(m))*(g*N()));
            }
	}
    } 

    itField = M_modelProperties->boundaryConditions().find("displacement");
    if (itField != M_modelProperties->boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("Dirichlet");
	if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
		auto marker = exAtMarker.marker();
                auto g = expr<Dim,G_Order>(exAtMarker.expression());
		blf( 2_c ) += integrate(_range=markedfaces(M_mesh,marker),
                		      _expr=trans(id(m))*g);
            }
        }
    }

    cout << "rhs works fine" << std::endl;
}

/*
template<int Dim, int Order, int G_Order>
template<typename PS, typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updateConductivityTerm( PS&& ps, Expr<ExprT> expr, std::string marker)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    MatConvert(toPETSc(M_A_cst)->mat(), MATSAME, MAT_INITIAL_MATRIX, &(toPETSc(M_A)->mat()));
    auto bbf = blockform2( ps, M_A);
    if ( marker.empty() )
        bbf(0_c,0_c) += integrate( _range=elements(M_mesh), _expr=inner(idt(u),id(v))/expr);
    else
        bbf(0_c,0_c) += integrate( _range=markedelements(M_mesh, marker), _expr=inner(idt(u),id(v))/expr);

}*/

/*
template<int Dim, int Order, int G_Order>
template<typename PS>
void
MixedPoisson<Dim, Order, G_Order>::updateConductivityTerm( PS&& ps, bool isNL)
{
    auto u = M_Vh->element( "u" );
    auto v = M_Vh->element( "v" );
    MatConvert(toPETSc(M_A_cst)->mat(), MATSAME, MAT_INITIAL_MATRIX, &(toPETSc(M_A)->mat()));

    auto bbf = blockform2( ps, M_A);
    for( auto const& pairMat : M_modelProperties->materials() )
    {
        auto marker = pairMat.first;
        auto material = pairMat.second;
        if ( !isNL )
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivity_json")));
            // (sigma^-1 j, v)
            bbf(0_c,0_c) += integrate(_range=markedelements(M_mesh,marker),
                                      _expr=(trans(idt(u))*id(v))/cond );
        }
        else
        {
            auto cond = material.getScalar(soption(prefixvm(M_prefix,"conductivityNL_json")), "p", idv(*M_pp));
	    // (sigma(p)^-1 j, v)
            bbf(0_c,0_c) += integrate(_range=markedelements(M_mesh,marker),
                                      _expr=(trans(idt(u))*id(v))/cond );
        }
    }
}*/

/*
template<int Dim, int Order, int G_Order>
template<typename PS, typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updatePotentialRHS( PS&& ps, Expr<ExprT> expr, std::string marker)
{
    auto blf = blockform1( ps, M_F);
    auto w = M_Wh->element();
    if ( marker.empty() )
        blf(1_c) += integrate(_range=elements(M_mesh), _expr=expr*id(w) );
    else
        blf(1_c) += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(w)) );
}
*/
/*
template<int Dim, int Order, int G_Order>
template<typename PS, typename ExprT>
void
MixedPoisson<Dim, Order, G_Order>::updateFluxRHS( PS&& ps, Expr<ExprT> expr, std::string marker)
{
    auto blf = blockform1( ps, M_F);
    auto v = M_Vh->element();
    if ( marker.empty() )
        blf(0_c) += integrate(_range=elements(M_mesh), _expr=expr*id(v) );
    else
        blf(0_c) += integrate( _range=markedelements(M_mesh, marker),
                          _expr=inner(expr,id(v)) );
}*/


// template <int Dim, int Order, int G_Order>
// void
// MixedPoisson<Dim,Order, G_Order>::exportResults( double time, mesh_ptrtype mesh  )
// {
//     this->log("MixedPoisson","exportResults", "start");
//     this->timerTool("PostProcessing").start();

//     if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC && mesh  )
//     {
//         LOG(INFO) << "exporting on visualisation mesh at time " << time;
//         M_exporter->step( time )->setMesh( mesh );
//     }
//     else if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
//     {
//         LOG(INFO) << "exporting on computational mesh at time " << time;
//         M_exporter->step( time )->setMesh( M_mesh );
//     }
    
//     // Export computed solutions
//     auto postProcess = M_modelProperties->postProcess();
//     auto itField = postProcess.find( "Fields");
//     if ( itField != postProcess.end() )
//     {
//         for ( auto const& field : (*itField).second )
//         {
//             if ( field == "flux" )
//             {
//                 LOG(INFO) << "exporting flux at time " << time;
//                 // M_exporter->step( time )->add(prefixvm(M_prefix, "flux"), Idhv?(*Idhv)( *M_up):*M_up );
//                 if (M_integralCondition)
//                 {
//                     double meas = 0.0;
//                     double j_integral = 0;
//                     for( auto marker : this->M_integralMarkersList)
//                     {
//                         LOG(INFO) << "exporting integral flux at time " << time << " on marker " << marker;
//                         j_integral += integrate(_range=markedfaces(M_mesh,marker),_expr=trans(idv(M_up))*N()).evaluate()(0,0);
//                         meas += integrate(_range=markedfaces(M_mesh,marker),_expr=cst(1.0)).evaluate()(0,0);
//                     }
//                     M_exporter->step( time )->add(prefixvm(M_prefix, "integralFlux"), j_integral);
// 		    M_exporter->step( time )->add(prefixvm(M_prefix, "integralVelocity"), j_integral/meas);
//                 }
//             }
//             else if ( field == "potential" )
//             {
//                 LOG(INFO) << "exporting potential at time " << time;
//                 M_exporter->step( time )->add(prefixvm(M_prefix, "potential"), Idh?(*Idh)(*M_pp):*M_pp);
//                 if (M_integralCondition)
//                 {
//                     LOG(INFO) << "exporting IBC potential at time " << time << " value " << (*M_mup)[0];
//                     M_exporter->step( time )->add(prefixvm(M_prefix, "cstPotential_1"),(*M_mup)[0] );
// 		    Feel::cout << "Integral value of potential(mup) on " << M_integralMarkersList.front() << " : \t " << (*M_mup)[0] << std::endl;
//                 }
// 		if (M_integralCondition == 2)
// 		{
//                     LOG(INFO) << "exporting IBC_2 potential at time " << time << " value " << (*M_mup2)[0];
//                     M_exporter->step( time )->add(prefixvm(M_prefix, "cstPotential_2"),(*M_mup2)[0] );
// 		    Feel::cout << "Integral value of potential(mup) on " << M_integralMarkersList.back() << " : \t " << (*M_mup2)[0] << std::endl;
// 		}
// 		auto itField = M_modelProperties->boundaryConditions().find("Exact solution");
// 		if ( itField != M_modelProperties->boundaryConditions().end() )
// 		{
// 		    auto mapField = (*itField).second;
// 		    auto itType = mapField.find( "p_exact" );
// 		    if (itType != mapField.end() )
// 		    {
// 			for (auto const& exAtMarker : (*itType).second )
// 			{
//    			    if (exAtMarker.isExpression() )
// 			    {
// 				auto p_exact = expr(exAtMarker.expression() );
// 				if ( !this->isStationary() )
//     			      	    p_exact.setParameterValues( { {"t", time } } );
//     				double K = 1; 
// 				for( auto const& pairMat : M_modelProperties->materials() )
//                 		{
//                     		   auto material = pairMat.second;
//                     		   K = material.getDouble( "k" );
//                 		}
//     				auto gradp_exact = grad<Dim>(p_exact);
//    				auto u_exact = expr(-K*trans(gradp_exact));
// 				M_exporter->step( time )->add(prefixvm(M_prefix, "p_exact"), project( _space=M_Wh, _range=elements(M_mesh), _expr=p_exact) );	
// 				M_exporter->step( time )->add(prefixvm(M_prefix, "u_exact"), project( _space=M_Vh, _range=elements(M_mesh), _expr=u_exact) );
				
// 				auto l2err_u = normL2( _range=elements(M_mesh), _expr=u_exact - idv(*M_up) );
// 				auto l2norm_uex = normL2( _range=elements(M_mesh), _expr=u_exact );
// 				if (l2norm_uex < 1)
// 				    l2norm_uex = 1.0;

//     				auto l2err_p = normL2( _range=elements(M_mesh), _expr=p_exact - idv(*M_pp) );
//     				auto l2norm_pex = normL2( _range=elements(M_mesh), _expr=p_exact );
//     				if (l2norm_pex < 1)
// 				    l2norm_pex = 1.0;
						
//     				Feel::cout << "----- Computed Errors -----" << std::endl;
//     				Feel::cout << "||p-p_ex||_L2=\t" << l2err_p/l2norm_pex << std::endl;
// 				Feel::cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
// 				Feel::cout << "---------------------------" << std::endl;
				
//     				// Export the errors
// 				M_exporter -> step( time )->add(prefixvm(M_prefix, "p_error_L2"), l2err_p/l2norm_pex );
// 				M_exporter -> step( time )->add(prefixvm(M_prefix, "u_error_L2"), l2err_u/l2norm_uex );
// 			    } 
// 			}
// 		    }

// 		}
			
//             } else
//             {
//        		// Import data
// 		LOG(INFO) << "importing " << field << " at time " << time;
//                 double extra_export = 0.0;
//                 auto itField = M_modelProperties->boundaryConditions().find( "Other quantities");
//   		if ( itField != M_modelProperties->boundaryConditions().end() )
//     		{
// 		    auto mapField = (*itField).second;
// 		    auto itType = mapField.find( field );
// 		    if ( itType != mapField.end() )
//         	    {		
//             	    	for ( auto const& exAtMarker : (*itType).second )
//             		{
// 			    if ( exAtMarker.isExpression() )
//                             {
// 				 LOG(INFO) << "WARNING: you are trying to export a single expression";
//                 	    }
//                 	    else if ( exAtMarker.isFile() )
//                 	    {
//                     		if ( !this->isStationary() )
//                     		{
//                         	    extra_export = exAtMarker.data(M_bdf_mixedpoisson->time());
//                     		}
//                    	    	else
//                         	    extra_export = exAtMarker.data(0.1);
// 			    }
//              		}
//             	    }
//                 }

// 		// Transform data if necessary
// 		LOG(INFO) << "transforming " << field << "at time " << time;
// 		std::string field_k = field;
// 		field_k += "_k";
// 		double kk = 0.0;
// 		for( auto const& pairMat : M_modelProperties->materials() )
//                 {
//                     auto material = pairMat.second;
//                     kk = material.getDouble( field_k );
//                 }
// 		if (std::abs(kk) > 1e-10)
// 		    extra_export *= kk;

// 		// Export data
//                 LOG(INFO) << "exporting " << field << " at time " << time;
// 		M_exporter->step( time )->add(prefixvm(M_prefix, field), extra_export);
//             }
//         }
//     }

//     /*
//     double Ui_mean = 0;
//     double meas = 0;
//     for( auto marker : this->M_integralMarkersList)
//     {
//         Ui_mean += integrate(_range=markedfaces(this->mesh(),marker),_expr=idv(*M_pp) ).evaluate()(0,0);
// 	meas += integrate(_range=markedfaces(M_mesh,marker),_expr=cst(1.0)).evaluate()(0,0);
//     }
//     if (M_integralCondition)
// 	Feel::cout << "Integral value of potential(mup) on " << M_integralMarkersList.front() << " : \t " << (*M_mup)[0] << std::endl;
//     if ( M_integralCondition == 2)
// 	Feel::cout << "Integral value of potential(mup) on " << M_integralMarkersList.back() << " : \t " << (*M_mup2)[0] << std::endl;
//     // Feel::cout << "Integral value of potential(mean u): \t " << Ui_mean/meas << std::endl;
//     */

//     this->timerTool("PostProcessing").stop("exportResults");
//     if ( this->scalabilitySave() )
//     {
//         if ( !this->isStationary() )
//             this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
//         this->timerTool("PostProcessing").save();
//     }
//     this->log("MixedPoisson","exportResults", "finish");
// }

} // Namespace FeelModels

} // Namespace Feel

#endif
