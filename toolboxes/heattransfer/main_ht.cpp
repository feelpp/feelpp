
#include <feel/feelmodels/heattransfer/heattransfer.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description heattransferoptions( "heat-transfer options" );
    heattransferoptions.add( toolboxes_options("heat-transfer") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=heattransferoptions,
                     _about=about(_name="toolboxes_heattransfer",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::HeatTransfer< Simplex<FEELPP_DIM,1>,
                                      Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> heatTransfer( new model_type("heat-transfer") );
    heatTransfer->init();
    heatTransfer->printAndSaveInfo();

    if ( heatTransfer->isStationary() )
    {
        heatTransfer->solve();
        heatTransfer->exportResults();
    }
    else
    {
        for ( ; !heatTransfer->timeStepBase()->isFinished(); heatTransfer->updateTimeStep() )
        {
            if (heatTransfer->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heatTransfer->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heatTransfer->solve();
            heatTransfer->exportResults();
        }
    }

    return 0;
}
