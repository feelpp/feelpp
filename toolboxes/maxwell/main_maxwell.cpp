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

    Environment env( _argc=argc, _argv=argv,
                     _desc=maxwelloptions,
                     _about=about(_name="application_maxwell",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runApplicationMaxwell();

    return 0;

}
