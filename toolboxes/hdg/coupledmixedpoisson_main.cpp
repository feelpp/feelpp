#include <feel/feelmodels/hdg/coupledmixedpoisson.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-model" ,
                     "mixed-poisson-model" ,
                     "0.1",
                     "Mixed-Poisson-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016-2019 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;
}

inline po::options_description
makeCoupledMixedPoissonOptions()
{
    po::options_description cmpOptions( "Coupled Mixed Poisson HDG options" );
    cmpOptions.add( toolboxes_options("mixedpoisson", "hdg.poisson") );
    return cmpOptions;
}



template<int nDim, int nOrder>
int
runApplicationCoupledMixedPoisson()
{
    using namespace Feel;

    typedef FeelModels::CoupledMixedPoisson<Simplex<nDim>,nOrder> cmp_type;

    auto toolbox = std::make_shared<cmp_type>();
    toolbox->init();
    toolbox->printAndSaveInfo();
    if ( toolbox->isStationary() )
    {
        Feel::cout << "Model should not be stationary !" << std::endl;
        return 1;
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


int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description mpoptions( "hdg.poisson options" );
    mpoptions.add( makeCoupledMixedPoissonOptions() );
    mpoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=mpoptions
                           );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    int status = 1;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                    [&discretization,&dimension,&status]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                            std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                            constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                            if ( dimension == _dim && discretization == _discretization )
                                status = runApplicationCoupledMixedPoisson<_dim,_torder>();
                        } );

    return status;
}
