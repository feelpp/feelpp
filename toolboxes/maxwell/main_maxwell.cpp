#include <feel/feelmodels/maxwell/maxwell.hpp>
#include <feel/feelmodels/modelcore/options.hpp>

// template<int Order>
void
runApplicationMaxwell()
{
    using namespace Feel;

    using model_type = FeelModels::Maxwell< Simplex<FEELPP_DIM,1> >;

    boost::shared_ptr<model_type> maxwell( new model_type("maxwell") );
    maxwell->init();
    maxwell->printAndSaveInfo();
    maxwell->solve();
    maxwell->exportResults();
}

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description maxwelloptions( "application maxwell options" );
    maxwelloptions.add( toolboxes_options("maxwell") );
    maxwelloptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2,P3 ")
        ;

    Environment env( _argc=argc, _argv=argv,
                     _desc=maxwelloptions,
                     _about=about(_name="application_maxwell",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    runApplicationMaxwell();
//     if ( feapprox == "P1" )
//         runApplicationMaxwell<1>();
// #if FEELPP_INSTANTIATION_ORDER_MAX >= 2
//     else if ( feapprox == "P2" )
//         runApplicationMaxwell<2>();
// #endif
// #if FEELPP_INSTANTIATION_ORDER_MAX >= 3
//     else if ( feapprox == "P3" )
//         runApplicationMaxwell<3>();
// #endif
//     else
//         CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;

}
