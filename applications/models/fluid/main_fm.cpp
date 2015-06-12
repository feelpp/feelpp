/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/fluid/fluidmechanics.hpp>

int
main( int argc, char** argv )
{

    using namespace Feel;
	po::options_description fluidmecoptions( "application fluid-mechanics options" );
    fluidmecoptions.add( feelmodels_options("fluid") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=fluidmecoptions,
                   _about=about(_name="application_fluid",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::FluidMechanics< Simplex<2,1,2>,
                                        Lagrange<2, Vectorial,Continuous,PointSetFekete>/*,
                                        Lagrange<1, Scalar,Continuous,PointSetFekete>*/ > model_type;
    boost::shared_ptr<model_type> FM( new model_type("fluid") );

    FM->init();
    FM->printAndSaveInfo();

    if ( FM->isStationary() )
    {
        FM->solve();
        FM->exportResults();
    }
    else
    {
        for ( ; !FM->timeStepBase()->isFinished(); FM->updateTimeStep() )
        {
            if (FM->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << FM->time() << "s \n";
                std::cout << "============================================================\n";
            }

            FM->solve();
            FM->exportResults();
        }
    }

    return 0;
}
