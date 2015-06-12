/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/solid/solidmechanics.hpp>



int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description solidmecoptions( "application solid-mechanics options" );
    solidmecoptions.add( feelmodels_options("solid") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=solidmecoptions,
                   _about=about(_name="application_solid",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    typedef FeelModels::SolidMechanics< Simplex<2,1,2>,
                                        Lagrange<1, Vectorial,Continuous,PointSetFekete> > model_type;
    boost::shared_ptr<model_type> SM( new model_type("solid") );

    SM->init();
    SM->printAndSaveInfo();

    if ( SM->isStationary() )
    {
        SM->solve();
        SM->exportResults();
    }
    else
    {
        for ( ; !SM->timeStepBase()->isFinished(); SM->updateTimeStep() )
        {
            if (SM->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << SM->time() << "s \n";
                std::cout << "============================================================\n";
            }

            SM->solve();
            SM->exportResults();
        }
    }

    return 0;
}
