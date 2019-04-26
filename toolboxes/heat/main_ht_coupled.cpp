/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/



#include <feel/feelmodels/heat/heat.hpp>
#include <boost/process.hpp>
template <int nDim,int OrderT>
void
runApplicationHeat()
{
    namespace bp = boost::process;
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                                      Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;
    std::shared_ptr<model_type> heat( new model_type("heat") );
    heat->init();
    heat->printAndSaveInfo();

    if ( heat->isStationary() )
    {
        heat->solve();
        heat->exportResults();
    }
    else
    {
        if ( !heat->doRestart() )
            heat->exportResults(heat->timeInitial());

        double the_value_computed_in_heat_sh = 0.;
        
        for ( heat->startTimeStep() ; !heat->timeStepBase()->isFinished(); heat->updateTimeStep() )
        {
            if (heat->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heat->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heat->addParameterInModelProperties( "Q", the_value_computed_in_heat_sh );
            
            heat->solve();

            // work only on process rank 0
            if ( heat->worldComm().isMasterRank() )
            {
                // Seed with a real random value, if available
                std::random_device r;
 
                // Choose a random mean between 1 and 6
                std::default_random_engine e1(r());
                std::uniform_int_distribution<int> uniform_dist(1e6, 1.1e6);
                int mean = uniform_dist(e1);
                std::cout << "Randomly-chosen mean: " << mean << '\n';
                
                boost::filesystem::path p = bp::search_path("heat.sh"); 
                int result = bp::system(p, std::to_string(mean));
                std::ifstream ifs("q.dat");
                std::string q;
                ifs >> q >> the_value_computed_in_heat_sh;
                    
            }
            mpi::broadcast( heat->worldComm(), the_value_computed_in_heat_sh, heat->worldComm().rank() );
            heat->exportResults();
        }
    }
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );
    heatoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            runApplicationHeat<_dim,_torder>();
                    } );
    return 0;
}
