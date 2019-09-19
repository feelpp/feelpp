#include "thermoelectric.hpp"

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
void
runApplicationThermoElectric()
{
    using namespace Feel;

    using thermoelectric_type = ThermoElectricHDG<nDim, OrderT, OrderV, OrderG>;
    using thermoelectric_ptrtype = std::shared_ptr<thermoelectric_type>;

    auto te = std::make_shared<thermoelectric_type>();
    te->run();
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

    auto dimt = hana::make_tuple(/*hana::int_c<2>,*/hana::int_c<3>);
    auto discretizationg = hana::make_tuple(hana::int_c<1>,hana::int_c<2>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( //hana::make_tuple("T3", hana::int_c<3> ),
                                             hana::make_tuple("T2", hana::int_c<2> ),
                                             hana::make_tuple("T1", hana::int_c<1> ) );
    auto discretizationv = hana::make_tuple( //hana::make_tuple("V3", hana::int_c<3> ),
                                             hana::make_tuple("V2", hana::int_c<2> ),
                                             hana::make_tuple("V1", hana::int_c<1> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("T2", hana::int_c<2> ),
                                             hana::make_tuple("T1", hana::int_c<1> ) );
    auto discretizationv = hana::make_tuple( hana::make_tuple("V2", hana::int_c<2> ),
                                             hana::make_tuple("V1", hana::int_c<1> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("T1", hana::int_c<1> ) );
    auto discretizationv = hana::make_tuple( hana::make_tuple("V1", hana::int_c<1> ) );
#endif

    bool hasRun = false;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationg,hana::cartesian_product(discretizationt,discretizationv))),
                    [&discretization,&dimension,&hasRun]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            constexpr int _gorder = std::decay_t<decltype(hana::at_c<1>(d))>::value;
                            std::string const& _discretizationt = hana::at_c<0>( hana::at_c<0>( hana::at_c<2>(d) ) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<0>( hana::at_c<2>(d) )))>::value;
                            std::string const& _discretizationv = hana::at_c<0>( hana::at_c<1>( hana::at_c<2>(d) ));
                            constexpr int _vorder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>( hana::at_c<2>(d) )))>::value;
                            if ( dimension == _dim && discretization == _discretizationt+_discretizationv+_gorder )
                            {
                                hasRun = true;
                                runApplicationThermoElectric<_dim,_torder,_vorder,_gorder>();
                            }
                        } );
    if( !hasRun )
    {
        Feel::cout << tc::red << "Wrong dimension (" << dimension << ") or discretization (" << discretization
                   << ") Possible combination:" << tc::reset << std::endl;
        hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationg,hana::cartesian_product(discretizationt,discretizationv))),
                        []( auto const& d )
                            {
                                Feel::cout << "\t(" << std::decay_t<decltype(hana::at_c<0>(d))>::value << ","
                                           << hana::at_c<0>( hana::at_c<1>(d) ) << hana::at_c<0>( hana::at_c<2>(d) )
                                           << std::decay_t<decltype(hana::at_c<1>(d))>::value <<  ")" << std::endl;
                            }
                        );
    }
    return 0;
}
