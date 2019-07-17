#include "cvg_mixedelasticity.hpp"

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "conv-mixed-elasticity",
                     "conv-mixed-elasticity",
                     "0.1",
                     "Convergence Mixed Elasticity Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;

}

template<int nDim, int OrderP, int OrderG, int OrderE>
int
runApplicationCvgElasticity()
{
    using namespace Feel;
    typedef ConvergenceElasticityTest<nDim,OrderP,OrderG,OrderE> cv_type;

    Environment::changeRepository( boost::format("conv-mixed-elasticity_%1%DP%2%G%3%E%4%") % nDim % OrderP % OrderG % OrderE );

    cv_type CV;
    return CV.run();
}

int main(int argc, char *argv[])
{
    using namespace Feel;
    po::options_description meoptions( "mixedelasticity options" );
    meoptions.add( makeConvOptions() );
    meoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=meoptions,
                           _desc_lib=makeConvLibOptions().add(feel_options()) );

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
                            runApplicationCvgElasticity<_dim,_torder,1,4>();
                    } );

    return 0;
}
