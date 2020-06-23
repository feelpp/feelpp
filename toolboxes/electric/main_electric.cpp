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

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> )
#if FEELPP_INSTANTIATION_ORDER_MAX >= 2
                                             ,hana::make_tuple("P2", hana::int_c<2> )
#endif
// #if FEELPP_INSTANTIATION_ORDER_MAX >= 3
//                                              ,hana::make_tuple("P3", hana::int_c<3> )
// #endif
                                             );

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            runApplicationElectric<_dim,_torder>();
                    } );
    return 0;
}
