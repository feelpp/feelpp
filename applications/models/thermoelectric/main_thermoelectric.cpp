
#include <feel/feelmodels/thermoelectric/thermoelectric.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description thermoelectricoptions( "application thermo-dynamics options" );
    thermoelectricoptions.add( feelmodels_options("thermo-electric") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=thermoelectricoptions,
                     _about=about(_name="application_thermoelectric",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::ThermoElectric< Simplex<FEELPP_DIM,1>,
                                        Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> thermoElectric( new model_type("thermo-electric") );
    thermoElectric->init();
    thermoElectric->printAndSaveInfo();

    if (thermoElectric->isStationary() )
    {
        thermoElectric->solve();
        thermoElectric->exportResults();
    }
    else
    {
        for ( ; !thermoElectric->timeStepBase()->isFinished(); thermoElectric->updateTimeStep() )
        {
            if (thermoElectric->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << thermoElectric->time() << "s \n";
                std::cout << "============================================================\n";
            }

            thermoElectric->solve();
            thermoElectric->exportResults();
        }
    }

    return 0;
}
