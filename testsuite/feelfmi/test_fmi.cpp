#define BOOST_TEST_MODULE test_fmi

#include <feel/feelcore/testsuite.hpp>
#if FEELPP_HAS_FMILIB
#include <feel/feelfmi/fmu.hpp>
#elif FEELPP_HAS_FMI4CPP
#include <feel/feelfmi/fmi4cpp.hpp>
#endif
using namespace Feel;

inline AboutData makeAbout()
{
    AboutData about( "test_fmi" ,
                     "test_fmi" ,
                     "0.1",
                     "FMI test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( test_fmi )

#if FEELPP_HAS_FMI4CPP
BOOST_AUTO_TEST_CASE( test_fmi2cs_fmi4cpp )
{
    Environment::setOptionValue<std::string>( "fmu.filename", Environment::expand("${cfgdir}/fmi2cs/fmi2cs.fmu") );
    auto fmu = fmi2::fmu( soption("fmu.filename") ).as_cs_fmu();
    auto md = fmu->get_model_description();
    size_t numOutputs = 0;
    for ( const auto& v : *fmu->get_model_description()->model_variables )
    {
        if ( v.causality == fmi2::causality::output )
        {
            numOutputs++;
        }
    }
    BOOST_CHECK_EQUAL( 0, numOutputs );

    auto slave = fmu->new_instance();
    std::cout << "model_identifier=" << slave->get_model_description()->model_identifier << std::endl;
    BOOST_CHECK( slave->setup_experiment() );
    BOOST_CHECK( slave->enter_initialization_mode() );
    BOOST_CHECK( slave->exit_initialization_mode() );
    std::vector<fmi2Real> ref( 2 );
    std::vector<fmi2ValueReference> vr = { md->get_variable_by_name( "y" ).value_reference, md->get_variable_by_name( "lambda" ).value_reference  };
    auto y = md->model_variables->getByValueReference( vr[0] ).as_real();
    auto lambda = md->model_variables->getByValueReference( vr[1] ).as_real();
    double step_size = 1e-4;
    const double stop = 5;

    double t;
    while ( ( t = slave->get_simulation_time() ) <= stop )
    {
        BOOST_CHECK( slave->step( step_size ) );
        BOOST_CHECK( slave->read_real( vr, ref ) );
//        BOOST_CHECK_CLOSE( ref[0], 0.1602522260970522, 1e-2 );
//        BOOST_CHECK_CLOSE( ref[1], 3.068004452918325, 1e-2 );
//        std::cout << "t=" << t << ", y=" << ref[0] << std::endl;
        BOOST_CHECK_CLOSE( ref[0], exp(-ref[1]*t), 1 );
    }

    BOOST_CHECK( slave->terminate() );
} 
#endif
#if FEELPP_HAS_FMILIB
BOOST_AUTO_TEST_CASE( test_fmi2cs_fmilib )
{
    Environment::setOptionValue<std::string>("fmu.filename", "${cfgdir}/fmi2cs/fmi2cs.fmu");

  
    FMU my_fmu;
    my_fmu.load();
    my_fmu.initialize(0,1); // initialize the fmu with t_init=0 and t_final=1
    my_fmu.doSteps(0.5); // interate until t_stop=0.5
    double h = my_fmu.getValue<double>( "h" );
    double v = my_fmu.getValue<double>( "v" );
    BOOST_CHECK_SMALL( math::abs(h-0.1602522260970522), 1e-2 );
    BOOST_CHECK_SMALL( math::abs(v-3.068004452918325), 1e-2 );

    my_fmu.setValue( "h", 0.1602522260970522 );// change the value of h
    my_fmu.setValue( "v", 3.068004452918325 ); // change the value of v
    my_fmu.doSteps(1);                         // do steps until t=1
    h = my_fmu.getValue<double>( "h" );        // get value of h
    v = my_fmu.getValue<double>( "v" );        // get value of v
    BOOST_CHECK_SMALL( math::abs(h-0.4680044525539696), 1e-2 );
    BOOST_CHECK_SMALL( math::abs(v+1.836995547081675), 1e-2 );


    my_fmu.reset();    // reset the model
    my_fmu.simulate(); // perform simulation with defaut t_init and t_final. no need to call initalize with this method
    h = my_fmu.getValue<double>( "h" );
    v = my_fmu.getValue<double>( "v" );
    BOOST_CHECK_SMALL( math::abs(h-0.4680044525539696), 1e-2 );
    BOOST_CHECK_SMALL( math::abs(v+1.836995547081675), 1e-2 );
   
}
#endif 
BOOST_AUTO_TEST_SUITE_END()
