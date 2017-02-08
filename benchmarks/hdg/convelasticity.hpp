#include "feel/feelmodels/hdg/mixedelasticity.hpp"

namespace Feel {

// namespace FeelModels {

inline
po::options_description
makeConvOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "cvg.refine-nb", po::value<int>()->default_value(3), "number of refinement" )
        ( "cvg.refine-factor", po::value<double>()->default_value(2), "factor of refinement" )
        ( "cvg.u_exact", po::value<std::string>()->default_value( "{0,0,0}:x:y:z" ), "exact displacement" )
        ( "cvg.lambda", po::value<double>()->default_value( 1.), "lambda" )
		( "cvg.mu", po::value<double>()->default_value( 1.), "mu" )
        ;
    options.add( FeelModels::makeMixedElasticityOptions() );
    return options;
}

inline
po::options_description
makeConvLibOptions()
{
    po::options_description options ( "Mixed Elasticity lib options");
    options.add( FeelModels::makeMixedElasticityLibOptions() );
    return options ;
}

template<int Dim, int Order,int G_Order = 1,int E_Order = 2>
class ConvergenceElasticityTest
{
public:
    using convex_type = Simplex<Dim,G_Order>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    static const uint16_type expr_order = Order+E_Order;
    using expr_scalar_type = scalar_field_expression<expr_order>;
    using expr_vectorial_type = vector_field_expression<Dim,G_Order,expr_order>;

    using mixed_elasticity_type = FeelModels::MixedElasticity<Dim,Order,G_Order,E_Order>;
    using mixed_elasticity_ptrtype = boost::shared_ptr<mixed_elasticity_type> ;

    using flux_type = typename mixed_elasticity_type::Vh_element_t;
    using potential_type = typename mixed_elasticity_type::Wh_element_t;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    using markers_type = std::vector<std::string>;

private:
    mesh_ptrtype M_mesh;

    mixed_elasticity_ptrtype M_model;

    flux_type M_sigma;
    potential_type M_u;
    expr_scalar_type M_u_exact;

    markers_type M_markersDirichlet;
    markers_type M_markersNeumann;

public:
    ConvergenceElasticityTest();
    void assembleExact();
    void run();
};

template<int Dim, int Order, int G_Order, int E_Order>
ConvergenceElasticityTest<Dim,Order,G_Order,E_Order>::ConvergenceElasticityTest()
{
    M_model = mixed_elasticity_type::New("mixedelasticity");
    // M_u_exact = expr<Dim,1,expr_order>(soption("cvg.u_exact"));

	/*
    Feel::cout << "using lambda: " << doption("cvg.lambda") << std::endl;
    Feel::cout << "using mu: " << doption("cvg.mu") << std::endl;
    Feel::cout << "using u exact     : " << M_u_exact << std::endl;
    // Feel::cout << "grad_u_exact      : " << grad<Dim,expr_order>(M_p_exact) << std::endl;
	*/	

	
    /*auto repo = boost::format( "conv_mixedpoisson/%1%D/P%2%/N%3%/E%4%/Dir%5%Neu%6%Ibc%7%/" )
        % Dim
        % Order
        % G_Order
        % E_Order
        % M_markersDirichlet.size()
        % M_markersNeumann.size()
    Environment::changeRepository( repo );
	
    Feel::cout << "Results are now stored in " << tc::red << repo << tc::reset << std::endl;
	*/
}


template<int Dim, int Order, int G_Order, int E_Order>
void
ConvergenceElasticityTest<Dim,Order,G_Order,E_Order>::run()
{
    double h = doption("gmsh.hsize");
/*
    auto lambda = cst(doption("cvg.lambda"));
	auto mu = cst(doption("cvg.mu"));
    auto gradu_exact = grad(M_u_exact);
    auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
	auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;
*/
/*   
	std::ofstream cvg_sigma, cvg_u;
	cvg_sigma.open( "convergence_sigma.dat", std::ios::out | std::ios::trunc);
    cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    boost::format fmterOut("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    cvg_sigma << fmter % "#h" % "nDof" % "l2err";
    cvg_u << fmter % "#h" % "nDof" % "l2err";
*/	

    export_ptrtype e( export_type::New( "convergence") );

    for ( int i = 0; i < ioption("cvg.refine-nb"); i++)
    {
        M_mesh = loadMesh( _mesh=new mesh_type, _h=h);
        M_model -> init(M_mesh); // inside here assembleSTD

/*
        auto nDofSigma = M_model->fluxSpace()->nDof();
        auto nDofU = M_model->potentialSpace()->nDof();
        auto sigma_ex = M_model->fluxSpace()->element();
        auto u_ex = M_model->potentialSpace()->element();
        auto u_ex_mean = M_model->potentialSpace()->element();
        auto u_mean = M_model->potentialSpace()->element();
*/
        
		M_model->solve();	// inside here assembleF
		
		M_model->exportResults(M_mesh);

        // M_sigma = M_model->fluxField();
        // M_u = M_model->potentialField();


/*
        double errSigma = normL2(_range=elements(M_mesh), _expr=idv(M_sigma)-sigma_exact, _quad=_Q<expr_order>());
        double mean_u_exact = mean( elements(M_mesh), M_u_exact)(0,0);
        double mean_u = mean( elements(M_mesh), idv(M_u))(0,0);
        double errU = normL2(_range=elements(M_mesh), _expr=idv(M_u)-cst(mean_u)-M_u_exact+cst(mean_u_exact), _quad=_Q<expr_order>());

  		
	    cout << fmterOut % "h" % "nDofSigma" % "errSigma" % "nDofU" % "errU";
        cout << fmterOut % h % nDofSigma % errSigma % nDofU % errU;
        cvg_sigma << fmter % h % nDofSigma % errSigma;
        cvg_u << fmter % h % nDofU % errU;
		

        sigma_ex.on( elements(M_mesh), sigma_exact);
        u_ex.on( elements(M_mesh), M_u_exact);
        u_ex_mean.on( elements(M_mesh), M_u_exact-cst(mean_u_exact));
        u_mean.on( elements(M_mesh), idv(M_u)-cst(mean_u));
*/
        e->step(i)->setMesh(M_mesh);
/*
        e->step(i)->add("sigma", M_sigma);
        e->step(i)->add("u", M_u);
        e->step(i)->add("u_mean", u_mean);
        e->step(i)->add("sigma_ex", sigma_ex);
        e->step(i)->add("u_ex", u_ex);
*/
        e->save();

        h /= doption("cvg.refine-factor");
    }


#if 0
	// COMPUTE THE ERROR IN 2D FOR A CURVED LINE TO VERIFY GEOMETRICAL ORDER 

    auto itField = M_model->modelProperties().boundaryConditions().find("GeometricalTest");
    if ( itField != M_model->modelProperties().boundaryConditions().end() )
    {
	    auto mapField = itField -> second;
        auto itType = mapField.find( "force_F" );
        if (itType != mapField.end() )
        {
			for (auto const& exAtMarker : itType->second )
            {
				auto forceF = expr<Dim,1> (exAtMarker.expression() );
				auto curveError = forceF; //curvedForce - forceF);
				auto curvedForce = integrate(_range=markedfaces(M_mesh,exAtMarker.marker()), _expr = id(M_model->fluxField()) * N() );
				
				Feel::cout << "Error for geometrical order:\t" << curveError << std::endl;	
			}
		}
	}
	
#endif
/*
    cvg_sigma.close();
    cvg_u.close();
*/
}

// } // end namespace FeelModels
} // end namespace Feel
