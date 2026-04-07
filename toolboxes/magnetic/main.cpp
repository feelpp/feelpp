/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/magnetic/magnetic.hpp>

template <int nDim,typename FeBasisType>
int
runApplicationMagnetic()
{
    using namespace Feel;

    typedef FeelModels::Magnetic< Simplex<nDim,1>, FeBasisType > model_type;
    std::shared_ptr<model_type> magnetic( new model_type("magnetic") );
    magnetic->init();
    magnetic->printAndSaveInfo();
    magnetic->solve();
    magnetic->exportResults();

    return !magnetic->checkResults();
}

int
main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description magneticoptions( "magnetic options" );
        magneticoptions.add( toolboxes_options("magnetic") );
        magneticoptions.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "Ned1h0" ), "discretization : Ned1h0 ")
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=magneticoptions,
                        _about=about(_name="toolboxes_magnetic",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");

        auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
        auto discretizationt = hana::make_tuple( hana::make_tuple("Ned1h0", std::type_identity<Nedelec<0,NedelecKind::NED1>>{} ) );

        int status = 0;
        hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status]( auto const& d )
                                                                                             {
                                                                                                 constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                                                                                                 std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                                                                                                 using _feBasisType = typename std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ) )>::type;
                                                                                                 if ( dimension == _dim && discretization == _discretization )
                                                                                                     status = runApplicationMagnetic<_dim,_feBasisType>();
                                                                                             } );



        // hana::for_each( Pc_t<2,3,1,2>, [&discretization, &dimension]( auto const& d )
        //                 {
        //                     constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
        //                     constexpr int _torder = std::decay_t<decltype( hana::at_c<1>( d ) )>::value;
        //                     std::string const& _discretization = hana::at_c<2>( d );
        //                     if ( dimension == _dim && discretization == _discretization )
        //                         runApplicationMagnetic<_dim,_torder>(); } );
        return status;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
