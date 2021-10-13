/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/heatfluid/heatfluid.hpp>

template <int nDim, int OrderT, int OrderV, int OrderP>
void
runApplicationHeatFluid()
{
    using namespace Feel;

    typedef FeelModels::Heat< Simplex<nDim,1>,
                              Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_heat_type;
    typedef FeelModels::FluidMechanics< Simplex<nDim,1>,
                                        Lagrange<OrderV, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderP, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::HeatFluid< model_heat_type,model_fluid_type> model_type;
    std::shared_ptr<model_type> heatFluid( new model_type("heat-fluid") );
    heatFluid->init();
    heatFluid->printAndSaveInfo();

    if (heatFluid->isStationary() )
    {
        heatFluid->solve();
        heatFluid->exportResults();
    }
    else
    {
        if ( !heatFluid->doRestart() )
            heatFluid->exportResults();

        for ( heatFluid->startTimeStep() ; !heatFluid->timeStepBase()->isFinished(); heatFluid->updateTimeStep() )
        {
            if (heatFluid->worldComm().isMasterRank())
            {
                std::cout << "============================================================\n";
                std::cout << "time simulation: " << heatFluid->time() << "s \n";
                std::cout << "============================================================\n";
            }

            heatFluid->solve();
            heatFluid->exportResults();
        }
    }
}

int
main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description heatfluidoptions( "application heat-fluid options" );
    heatfluidoptions.add( toolboxes_options("heat-fluid") );
    heatfluidoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1-P2P1" ), "discretization : P1-P2P1")
     ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=heatfluidoptions,
                     _about=about(_name="feelpp_toolbox_heatfluid",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);

    auto discretizationt = hana::make_tuple( hana::make_tuple("P1-P1P1", hana::make_tuple( hana::int_c<1>,hana::int_c<1>,hana::int_c<1>) ),
                                             hana::make_tuple("P1-P2P1", hana::make_tuple( hana::int_c<1>,hana::int_c<2>,hana::int_c<1>) ) );

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _uorder = std::decay_t<decltype(hana::at_c<1>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _porder = std::decay_t<decltype(hana::at_c<2>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            runApplicationHeatFluid<_dim,_torder,_uorder,_porder>();
                    } );

    return 0;
}
