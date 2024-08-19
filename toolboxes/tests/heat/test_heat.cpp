#define BOOST_TEST_MODULE heat testsuite
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelmodels/execute.hpp>
#include <feel/feelmodels/heat/heat.hpp>

using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Heat Options" );

    options.add( feel_options() );
    options.add( toolboxes_options("heat") );
    options.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ("ensemble.split", Feel::po::value<int>()->default_value( 2 ), "split ensemble" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_heat" ,
                     "test_heat" ,
                     "0.1",
                     "Heat test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2022-2024 Feel++ Consortium" );

    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( heatsuite )

BOOST_AUTO_TEST_CASE( test_heat_ensemble )
{
    int nsplit = ioption("ensemble.split");
    if (  Environment::numberOfProcessors() % nsplit == 0 )
    {
        auto [color, w, wglob] = Environment::worldCommPtr()->split( nsplit );
        using toolbox_t = Feel::FeelModels::Heat< Simplex<2,1>, Lagrange<1, Scalar,Continuous,PointSetFekete> >;
        fs::path ensembleRepo = fmt::format( "{}/ensemble/{}_{}/{}/", Environment::appRepository(), Environment::numberOfProcessors(), nsplit, color );
        fs::path ensembleRepoExpr = fmt::format( "{}/ensemble/{}_{}/{}/exprs/", Environment::appRepository(), Environment::numberOfProcessors(), nsplit, color );
        LOG(INFO) << "ensembleRepo = " << ensembleRepo.string();
        w->print( fmt::format( "[rank {}] ensemble repo : {}\n", w->rank(), ensembleRepo ), FLAGS_v > 0, FLAGS_v > 0, FLAGS_v > 0 );
        auto toolbox = toolbox_t::New( _prefix = "heat", _worldcomm = w,
                                       _repository = Feel::FeelModels::ModelBaseRepository( ensembleRepo.string(),false, ensembleRepoExpr ) );
        if ( toolbox->isStationary() )
            toolbox->setTimeFinal( 10 * toolbox->timeStep() );
        execute( toolbox, _verbose=1, _save=false );
    }
}

BOOST_AUTO_TEST_SUITE_END()
