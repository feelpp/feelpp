#include <feel/feelmodels/maxwell/maxwell.hpp>
#include <feel/feelmodels/modelcore/options.hpp>

// template<int Order>
template<int Dim>
void
runApplicationMaxwell()
{
    using namespace Feel;

    using model_type = FeelModels::Maxwell< Simplex<Dim,1> >;

    std::shared_ptr<model_type> maxwell( new model_type("maxwell") );
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
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ;

    Environment env( _argc=argc, _argv=argv,
                     _desc=maxwelloptions,
                     _about=about(_name="application_maxwell",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    if( dimension == 2 )
        runApplicationMaxwell<2>();
    else
        runApplicationMaxwell<3>();

    return 0;

}
