#include <feel/feelmodels/hdg/hdg_darcy.hpp>

namespace Feel
{

namespace HDG
{

template <uint16_type OrderPot>
void
runApplicationHDGdarcy()
{
    Darcy<??> D(??);
    D.solve();
    D.exportResults();
}

} // namespace HDG

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::HDG;

    po::options_description hdgoptions( "application HDG options");
    hdgoptions.add( darcy_options( "hdg-darcy" ) );

	Environment env( _argc=argc, _argv=argv,
                     _desc=hdgoptions,
                     _about=about(_name="application_hdg_darcy",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int feapprox = ioption("order-potential");
    runApplicationHDGdarcy<feapprox>();

    return 0;
}
