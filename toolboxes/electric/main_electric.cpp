/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/electric/electric.hpp>

template <int nDim,int OrderT>
void
runApplicationElectric()
{
    using namespace Feel;

    typedef FeelModels::Electric< Simplex<nDim,1>,
                                  Lagrange<OrderT, Scalar,Continuous,PointSetFekete> > model_type;
    std::shared_ptr<model_type> electric( new model_type("electric") );
    electric->init();
    electric->printAndSaveInfo();
    electric->solve();
    electric->exportResults();
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

        hana::for_each( Pc_t<>, [&discretization, &dimension]( auto const& d )
                        {
                            constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
                            constexpr int _torder = std::decay_t<decltype( hana::at_c<1>( d ) )>::value;
                            std::string const& _discretization = hana::at_c<2>( d );
                            if ( dimension == _dim && discretization == _discretization )
                                runApplicationElectric<_dim,_torder>(); } );
        return 0;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
