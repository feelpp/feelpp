
#include <feel/feelmodels/heat/heat.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::Heat< Simplex<FEELPP_DIM,1>,
                                      Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> heat( new model_type("heat") );
    heat->init();
    heat->printAndSaveInfo();

    if ( heat->isStationary() )
    {
        heat->solve();
        heat->exportResults();
    }
    else
    {
        for ( ; !heat->timeStepBase()->isFinished(); heat->updateTimeStep() )
        {
            if (heat->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heat->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heat->solve();
            heat->exportResults();
        }
    }

    return 0;
}
