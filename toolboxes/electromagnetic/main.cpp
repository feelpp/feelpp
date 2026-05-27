/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/electromagnetic/electromagnetic.hpp>
#include <feel/feelmodels/modelcore/convergencemode.hpp>

int
main(int argc, char**argv )
{
    using namespace Feel;
    try
    {
        po::options_description electromagneticoptions( "electromagnetic options" );
        electromagneticoptions.add( toolboxes_options("electromagnetic") );// TODO: to remove (missing use of custom vm in alg solver)
        electromagneticoptions.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "Pch1_Ned1h0" ), "discretization : Pch1_Ned1h0, Pch1_Pchv1" )
            ("case.mode", Feel::po::value<std::string>()->default_value( "simulation" ), "mode : simulation, h-convergence")
            ("case.mode.h-convergence.hsize", po::value<std::vector<double> >()->multitoken(), "mesh hsize used in h-convergence" )
            ("case.mode.h-convergence.measures", po::value<std::vector<std::string> >()->multitoken(), "measures names used in fit checker" )
            ("case.mode.h-convergence.slopes", po::value<std::vector<double> >()->multitoken(), "reference slope for the fit checker" )
            ;

        Environment env( _argc=argc, _argv=argv,
                         _desc=electromagneticoptions,
                         _about=about(_name="toolboxes_electromagnetic",
                                      _author="Feel++ Consortium",
                                      _email="feelpp-devel@feelpp.org"));

        std::string mode = soption(_name="case.mode");
        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");

        auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
        auto discretizationt = hana::make_tuple( hana::make_tuple( "Pch1_Pchv1",
                                                                   std::type_identity<Lagrange<1,Scalar,Continuous,PointSetFekete>>{},
                                                                   std::type_identity<Lagrange<1,Vectorial,Continuous,PointSetFekete>>{} ),
                                                 hana::make_tuple( "Pch1_Ned1h0",
                                                                   std::type_identity<Lagrange<1,Scalar,Continuous,PointSetFekete>>{},
                                                                   std::type_identity<Nedelec<0,NedelecKind::NED1>>{} )
                                                 );
        int status = 0;
        hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)),
                        [&discretization,&dimension,&status,&mode]( auto const& d )
                            {
                                constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                                std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                                using _feBasisElectricType = typename std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ) )>::type;
                                using _feBasisMagneticType = typename std::decay_t<decltype(hana::at_c<2>( hana::at_c<1>(d) ) )>::type;
                                using model_electric_type = FeelModels::Electric< Simplex<_dim,1>, _feBasisElectricType >;
                                using model_magnetic_type = FeelModels::Magnetic< Simplex<_dim,1>, _feBasisMagneticType >;

                                using model_type = FeelModels::Electromagnetic< model_electric_type, model_magnetic_type >;
                                if ( dimension == _dim && discretization == _discretization )
                                {
                                    if ( mode == "simulation" )
                                        status = Toolboxes::executeSingleRun<model_type>( "electromagnetic", "electromagnetic" );
                                    else if ( mode == "h-convergence" )
                                        status = Toolboxes::executeHConvergence<model_type>( "electromagnetic", "electromagnetic" );
                                }
                            } );
        return status;
    }
    catch(...)
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
