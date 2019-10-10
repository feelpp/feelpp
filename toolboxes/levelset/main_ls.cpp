
#include <feel/feelmodels/levelset/levelset.hpp>

namespace Feel
{

template <uint16_type OrderLevelset, uint16_type OrderLevelsetPN = LEVELSET_PN_ORDER >
void
runLevelsetApplication()
{
    using namespace Feel;

    typedef Simplex<FEELPP_DIM,1> convex_type;
    typedef Lagrange<OrderLevelset, Scalar, Continuous, PointSetFekete> basis_type;
    //typedef FunctionSpace<Mesh<convex_type>, bases<typename FeelModels::detail::ChangeBasisPolySet<Vectorial, basis_type>::type>, Periodicity<NoPeriodicity>> space_advection_velocity_type;
    typedef FunctionSpace<Mesh<convex_type>, Feel::detail::bases<Lagrange<OrderLevelset, Vectorial, Continuous, PointSetFekete>>, double, Periodicity<NoPeriodicity>, mortars<NoMortar>> space_advection_velocity_type;
    typedef Lagrange<OrderLevelsetPN, Scalar, Continuous, PointSetFekete> basis_PN_type;

    typedef FeelModels::LevelSet<
        convex_type,
        basis_type,
        NoPeriodicity,
        space_advection_velocity_type,
        basis_PN_type
        > model_type;
    
    auto LS = model_type::New("levelset");

    double ls_x0 = 0.5, ls_y0 = 0.5, ls_z0 = 0.5;
    double ls_radius = 0.25;

    //auto phi_init = sqrt( 
            //(vf::Px()-ls_x0)*(vf::Px()-ls_x0) 
            //+ (vf::Py()-ls_y0) * (vf::Py()-ls_y0) 
            ////+ (vf::Pz()-ls_z0) * (vf::Pz()-ls_z0) 
            //) - ls_radius;
    //LS->setInitialValue( phi_init );

    Feel::cout << "============================================================\n";
    Feel::cout << "Levelset toolbox with errors measures\n";
    Feel::cout << "============================================================\n";
    LS->init();
    LS->printAndSaveInfo();
    if( !LS->doRestart() )
        LS->exportResults(0.);

    // Errors measures :
    // LS is initialized, we save the phi, heaviside and dirac functions' initial values
    auto phi_0 = LS->functionSpace()->element();
    phi_0 = *LS->phi();
    auto H_0 = LS->functionSpace()->element();
    H_0 = *LS->heaviside();
    auto dirac_0 = LS->functionSpace()->element();
    dirac_0 = *LS->dirac();

    // Quadrature order used for integrals computed in errors measures:
    int error_quad_order = ioption( _name="levelset.quad.order" );
    Feel::cout << "Error quadrature order = " << error_quad_order << std::endl;

    
    bool exportDistToBoundary = boption( _name="export-dist-to-boundary" );

    std::shared_ptr<Exporter<typename model_type::mesh_type>> myExporter;

    if( exportDistToBoundary )
    {
        myExporter = exporter( 
                _mesh=LS->mesh(),
                _name="DistExport",
                _path=LS->exporterPath()
                );
    }

    if ( LS->isStationary() )
    {
        LS->solve();
        LS->exportResults();
    }
    else
    {
        int reinit_every = ioption( _name="levelset.reinit-every" );

        for ( int iter = 1; !LS->timeStepBase()->isFinished(); LS->updateTimeStep(), ++iter )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << LS->time() << "s \n";
            Feel::cout << "============================================================\n";

            Feel::cout << "Iter since reinit: " << LS->iterSinceReinit() << std::endl;
            Feel::cout << "Levelset BDF order: " << LS->timeStepBDF()->timeOrder() << std::endl;

            LS->solve();
            if( reinit_every > 0 && iter%reinit_every == 0 )
            {
                Feel::cout << "Reinitializing... ";
                LS->reinitialize();
                Feel::cout << "done\n";
            }

            LS->exportResults();
            if( exportDistToBoundary )
            {
                auto distToBoundary = LS->distToBoundary();
                myExporter->step(iter)->add("distToBoundary", *distToBoundary );
                myExporter->save();
            }

	    	    
	    // Error measures :
	    // we compute the L2 norm error integrals :

	    auto  chi_of_dirac0_positive = chi( idv(dirac_0) > 0 );
	    double first_integral = integrate(
					      _range=elements(LS->mesh()),
					      _expr=chi_of_dirac0_positive
					      ).evaluate()(0,0);
	    double second_integral = integrate(
					       _range=elements(LS->mesh()),
					       _expr=( pow( idv(phi_0)-idv(LS->phi()), 2.0 ) * chi_of_dirac0_positive ),
					       _quad=error_quad_order
					       ).evaluate()(0,0);
	    double l2_norm_error = std::sqrt( 1/first_integral * second_integral );

	    //Feel::cout << "First integral = " << first_integral << std::endl;
	    //Feel::cout << "Second integral = " << second_integral << std::endl;
	    Feel::cout << "L2 norm error = " << l2_norm_error << std::endl;

	    // Sign change error :
	    double e_sc_integral = integrate(
					     _range=elements(LS->mesh()),
					     _expr=pow(
						       ( 1-idv(H_0) ) - ( 1-idv(LS->heaviside()) ),
						       2.0
						       ),
					     _quad=error_quad_order
					     ).evaluate()(0,0);
	    double sign_change_error_h = std::sqrt(
						 e_sc_integral
						 );

	    Feel::cout << "Sign change error (h) = " << sign_change_error_h << std::endl;

	    double sign_change_error_chi = std::sqrt( integrate(
								_range=elements(LS->mesh()),
								_expr=chi( idv( LS->phi() ) * idv(phi_0) < 0 ),
								_quad=error_quad_order
								).evaluate()(0,0) );
	    Feel::cout << "Sign change error (chi) = " << sign_change_error_chi << std::endl;
	    
	    // Mass error
	    auto chi_of_phi0_negative = chi(  idv(phi_0) < 0 );
	    auto chi_of_phi_negative = chi(  idv(LS->phi()) < 0 );

	    double em_phi0_integral = integrate(
						_range=elements(LS->mesh()),
						_expr=chi_of_phi0_negative,
						_quad=error_quad_order
						).evaluate()(0,0);
	    double em_phi_integral = integrate(
					       _range=elements(LS->mesh()),
					       _expr=chi_of_phi_negative,
					       _quad=error_quad_order
					       ).evaluate()(0,0);
	    double mass_error = std::abs( em_phi_integral - em_phi0_integral ) / em_phi0_integral;

	    //Feel::cout << "phi0 integral = " << em_phi0_integral << std::endl;
	    //Feel::cout << "phi integral = " << em_phi_integral << std::endl;
	    Feel::cout << "Mass error = " << mass_error << std::endl;
        }
    }
}

}  // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description levelsetoptions( "application levelset options" );
    levelsetoptions.add( toolboxes_options( "levelset" ) );
    levelsetoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2, P1" )
        ("levelset.reinit-every", Feel::po::value<int>()->default_value( -1 ), "reinitialize levelset every n iterations" )
        ("export-dist-to-boundary", Feel::po::value<bool>()->default_value( false ), "compute and export the distance to the boundary" )
      ("levelset.quad.order", Feel::po::value<int>()->default_value( 1 ), "Quadrature order used for integrals computed in errors measures")
        ;

    Environment env( 
            _argc=argc, _argv=argv,
            _desc=levelsetoptions,
            _about=about(_name="application_levelset",
                _author="Feel++ Consortium",
                _email="feelpp-devel@feelpp.org")
            );

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P2" )
        runLevelsetApplication<2>();
    else if ( feapprox == "P1" )
        runLevelsetApplication<1>();
    
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
