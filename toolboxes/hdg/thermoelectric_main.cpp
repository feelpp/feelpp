#include <feel/feelmodels/hdg/thermoelectric.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "thermoelectric-hdg-model" ,
                     "thermoelectric-hdg-model" ,
                     "0.1",
                     "ThermoElectricHDG-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;
}

template<int nDim, int OrderT, int OrderV, int OrderG>
int
runApplicationThermoElectric()
{
    using namespace Feel;

    using thermoelectric_type = ThermoElectricHDG<nDim, OrderT, OrderV, OrderG>;
    using thermoelectric_ptrtype = std::shared_ptr<thermoelectric_type>;

    auto te = std::make_shared<thermoelectric_type>();
    te->init();
    te->printAndSaveInfo();
    tic();
    te->solve();
    toc("solve picard");
    te->exportResults();
    return !te->checkResults();
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description teoptions( "thermoelectric options" );
    teoptions.add( makeThermoElectricHDGOptions() );
    teoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=teoptions,
                           _desc_lib=feel_options()
                           );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(// hana::int_c<2>,
                                 hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( // hana::make_tuple("P3", hana::int_c<3> ),
                                             // hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P1", hana::int_c<1> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P1", hana::int_c<1> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    int status = -1;
    std::vector<std::string> combinations;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                    [&discretization,&dimension,&status,&combinations]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretizationt = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            combinations.push_back(std::to_string(_dim)+","+_discretizationt);
                            if ( dimension == _dim && discretization == _discretizationt )
                            {
                                status = runApplicationThermoElectric<_dim,_torder,_torder,1>();
                            }
                        } );
    if( status < 0 )
    {
        Feel::cout << tc::red << "Wrong dimension (" << dimension << ") or discretization (" << discretization
                   << ") Possible combination:" << tc::reset << std::endl;
        for( auto const& s : combinations )
            Feel::cout << "\t("<< s << ")" << std::endl;
        status = 1;
    }

    std::ofstream os ( "timers.md" );
    Environment::saveTimersMD(os);
    os.close();

    return status;
}
