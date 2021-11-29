#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <feel/feelmodels/modelcore/convergencemode.hpp>

using namespace Feel;

template <typename ToolboxType>
int
runToolboxSimulation( std::shared_ptr<ToolboxType> toolbox )
{
    toolbox->init();
    toolbox->printAndSaveInfo();

    if ( toolbox->isStationary() )
    {
        toolbox->solve();
        if( boption("use-postprocess") )
            toolbox->solvePostProcess();
        toolbox->exportResults();
    }
    else
    {
        if ( !toolbox->doRestart() )
            toolbox->exportResults(toolbox->timeInitial());

        for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
        {
            if (toolbox->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << toolbox->time() << "s \n";
                std::cout << "============================================================\n";
            }

            toolbox->solve();
            toolbox->exportResults();
        }
    }
    return !toolbox->checkResults();
}

template <typename ToolboxType>
int
runMixedPoissonSimulation()
{
    using namespace Feel;
    auto ME = ToolboxType::New("hdg.elasticity");
    return runToolboxSimulation( ME );
}

int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description mpoptions( "hdg.elasticity options" );
    mpoptions.add( toolboxes_options("mixedpoisson", "hdg.elasticity") );
    mpoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ("case.mode", Feel::po::value<std::string>()->default_value( "simulation" ), "mode : simulation, h-convergence")
        ("case.mode.h-convergence.hsize", po::value<std::vector<double> >()->multitoken(), "mesh hsize used in h-convergence" )
        ("use-postprocess", po::value<bool>()->default_value(false), "post process potential" )
        ;
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=mpoptions,
                           _about=about(_name="toolboxes_mixedpoisson",
                                        _author="Feel++ Consortium",
                                        _email="feelpp-devel@feelpp.org")
                           );

    std::string mode = soption(_name="case.mode");
    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ));
    int status = 1;

    if ( mode == "simulation" || mode == "h-convergence" )
    {
        hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status,&mode]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                            {
                                typedef Feel::FeelModels::MixedPoisson< Simplex<_dim,1>,_torder,Vectorial> model_type;
                                if ( mode == "simulation" )
                                    status = runMixedPoissonSimulation<model_type>();
                                else
                                    status = runHConvergence<model_type>("", runToolboxSimulation<model_type>);
                            }
                        } );
    }
    else if ( mode == "p-convergence" )
    {
        // TODO
    }
    return status;
}
