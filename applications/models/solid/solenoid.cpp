#include <feel/feelmodels/solid/solidmechanics.hpp>

namespace Feel
{

template <uint16_type OrderDisp>
void
runApplicationSolid()
{
    typedef FeelModels::SolidMechanics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> > model_type;
    auto SM = model_type::New("solid");

    SM->init();
    SM->printAndSaveInfo();
    SM->solve();
    SM->updateStressCriterions();
    auto u = SM->fieldDisplacement();
    auto s = SM->fieldsPrincipalStresses();
    auto vm = SM->fieldVonMisesCriterions();
    auto tr = SM->fieldTrescaCriterions();

    auto e = exporter( _mesh=mesh, _name="solidMec" );
    e->add( "displacement", u );
    e->add( "vonMises", vm );
    e->add( "tresca", tr );
    for( auto i : range(FEELPP_DIM))
        e->add( (boost::format( "principal-stress-%1%") %i).str(), s[i] );
    e->save();
}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description solidmecoptions( "application solid-mechanics options" );
    solidmecoptions.add( feelmodels_options("solid") );
    solidmecoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P1,P2 ")
        ("solve-quasi-static", Feel::po::value<bool>()->default_value(false), "solve-quasi-static ")
        ("solve-quasi-static.variable-initial", Feel::po::value<double>()->default_value(0.0), "solve-quasi-static.variable-initial")
        ("solve-quasi-static.variable-step", Feel::po::value<double>()->default_value(0.1), "solve-quasi-static.variable-step")
        ("solve-quasi-static.variable-symbol", Feel::po::value<std::string>()->default_value(""), "solve-quasi-static.variable-symbol")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=solidmecoptions,
                     _about=about(_name="application_solenoid",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P1" )
        runApplicationSolid<1>();
    else if ( feapprox == "P2" )
        runApplicationSolid<2>();
    else CHECK( false ) << "invalid feapprox " << feapprox;


    return 0;
}
