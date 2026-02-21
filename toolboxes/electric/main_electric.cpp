/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/electric/electric.hpp>

template <int nDim,int OrderT>
int
runApplicationElectric()
{
    using namespace Feel;

    typedef FeelModels::Electric< Simplex<nDim,1>,
                                  Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;
    auto electric = std::make_shared<model_type>("electric");
    electric->init();
    electric->printAndSaveInfo();
    electric->solve();
    electric->exportResults();
    return !electric->checkResults();
}

int
main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description electricoptions( "electric options" );
        electricoptions.add( toolboxes_options("electric") );
        electricoptions.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
            ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=electricoptions,
                        _about=about(_name="toolboxes_electric",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org"));

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");
        int status = 0;

        hana::for_each( Pc_t<2,3,1,2>, [&discretization, &dimension, &status]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
                            constexpr int _torder = std::decay_t<decltype( hana::at_c<1>( d ) )>::value;
                            std::string const& _discretization = hana::at_c<2>( d );
                            if ( dimension == _dim && discretization == _discretization )
                                status = runApplicationElectric<_dim,_torder>(); } );
        return status;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
