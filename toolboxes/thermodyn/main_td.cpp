
#include <feel/feelmodels/thermodyn/thermodynamics.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description thermodynoptions( "application thermo-dynamics options" );
    thermodynoptions.add( feelmodels_options("thermo-dynamics") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=thermodynoptions,
                     _about=about(_name="application_thermodyn",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::ThermoDynamics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> thermoDyn( new model_type("thermo") );
    thermoDyn->init();
    thermoDyn->printAndSaveInfo();

    if (thermoDyn->isStationary() )
    {
        thermoDyn->solve();
        thermoDyn->exportResults();
    }
    else
    {
        for ( ; !thermoDyn->timeStepBase()->isFinished(); thermoDyn->updateTimeStep() )
        {
            if (thermoDyn->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << thermoDyn->time() << "s \n";
                std::cout << "============================================================\n";
            }

            thermoDyn->solve();
            thermoDyn->exportResults();
        }
    }

    return 0;
}
