/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/fsi/fsi.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;

	po::options_description fsioptions( "application fsi options" );
    fsioptions.add( Feel::feelmodels_options("fsi") );

	Environment env( _argc=argc, _argv=argv,
                     _desc=fsioptions,
                     _about=about(_name="application_fsi",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    typedef FeelModels::FluidMechanics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<2, Vectorial,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::SolidMechanics< Simplex<FEELPP_DIM,1>,
                                        Lagrange<1, Vectorial,Continuous,PointSetFekete> > model_solid_type;
    typedef FeelModels::FSI< model_fluid_type,model_solid_type> model_fsi_type;
    boost::shared_ptr<model_fsi_type> FSImodel( new model_fsi_type("fsi") );

    FSImodel->init();
    FSImodel->printAndSaveInfo();

    for ( ; !FSImodel->timeStepBase()->isFinished(); FSImodel->updateTimeStep() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "\n====================================================================================="
                      << "\n current time : " << std::setprecision( 5 ) << std::fixed << FSImodel->currentTime()
                      << "\n=====================================================================================\n";
        FSImodel->solve();
        FSImodel->exportResults();
    }
    return 0;

}
