#define BOOST_TEST_MODULE test_fmi

#include <testsuite.hpp>

#include <feel/feelfmi/fmu.hpp>

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

BOOST_AUTO_TEST_CASE( test_fmi2cs )
{
    Environment::setOptionValue<std::string>("fmu.filename", "${cfgdir}/fmi2cs/fmi2cs.fmu");
    FMU my_fmu;
    my_fmu.load();
    my_fmu.simulate(0,0.5);
    double h = my_fmu.getValue<double>( "h" );
    double v = my_fmu.getValue<double>( "v" );
    BOOST_CHECK_SMALL( math::abs(h-0.1602522260970522), 1e-2 );
    BOOST_CHECK_SMALL( math::abs(v-3.068004452918325), 1e-2 );

    my_fmu.reset();
    my_fmu.simulate();
    h = my_fmu.getValue<double>( "h" );
    v = my_fmu.getValue<double>( "v" );
    BOOST_CHECK_SMALL( math::abs(h-0.4680044525539696), 1e-2 );
    BOOST_CHECK_SMALL( math::abs(v+1.836995547081675), 1e-2 );
}

BOOST_AUTO_TEST_SUITE_END()
