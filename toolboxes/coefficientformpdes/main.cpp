/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes_registered_type.hpp>

template <int nDim,int nOrderGeo>
int
runApplicationCoefficientFormPDEs()
{
    using namespace Feel;

    using model_type = FeelModels::coefficient_form_PDEs_t< Simplex<nDim,nOrderGeo> >;

    std::shared_ptr<model_type> cfpdes( new model_type("cfpdes") );
    cfpdes->init();
    cfpdes->printAndSaveInfo();

    if ( cfpdes->isStationary() )
    {
        cfpdes->solve();
        cfpdes->exportResults();
    }
    else
    {
        if ( !cfpdes->doRestart() )
            cfpdes->exportResults(cfpdes->timeInitial());

        for ( cfpdes->startTimeStep() ; !cfpdes->timeStepBase()->isFinished(); cfpdes->updateTimeStep() )
        {
            if (cfpdes->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << cfpdes->time() << "s \n";
                std::cout << "============================================================\n";
            }

            cfpdes->solve();
            cfpdes->exportResults();
        }
    }

    return !cfpdes->checkResults();
}

int
main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description cfpdesoptions( "coefficient-form-pdes options" );
    cfpdesoptions.add( toolboxes_options( "coefficient-form-pdes", "cfpdes" ) );
    cfpdesoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=cfpdesoptions,
                     _about=about(_name="toolboxes_cfpdes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    Feel::FeelModels::printToolboxApplication( "cfpdes" );

    int dimension = ioption(_name="case.dimension");
    //std::string discretization = soption(_name="case.discretization");
    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
    //auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
    int status = 0;
    hana::for_each( dimt/*hana::cartesian_product(hana::make_tuple(dimt,discretizationt))*/, [&dimension,&status]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(d)>::value;
                        if ( dimension == _dim  )
                            status = runApplicationCoefficientFormPDEs<_dim,1>();
                    } );
    return status;
}
