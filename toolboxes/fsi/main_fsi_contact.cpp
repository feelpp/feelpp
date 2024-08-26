/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/fsi/fsi.hpp>
#include "contactforce.hpp"

namespace Feel
{

template <uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderDisp=FEELPP_GEO_ORDER>
void
runApplicationFSI()
{
    using namespace Feel;

    typedef FeelModels::FluidMechanics< Simplex<FEELPP_DIM,FEELPP_GEO_ORDER>,
                                        Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    typedef FeelModels::SolidMechanics< Simplex<FEELPP_DIM,FEELPP_GEO_ORDER>,
                                        Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> > model_solid_type;
    typedef FeelModels::FSI< model_fluid_type,model_solid_type> model_fsi_type;
    std::shared_ptr<model_fsi_type> FSImodel( new model_fsi_type("fsi") );

    FSImodel->init();
    FSImodel->printAndSaveInfo();

    auto add_force_term = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateLinear & data)
    {
        auto const& t = unwrap_ptr(FSImodel->solidModel());
        std::string wall_normal = (FEELPP_DIM==2)?std::string("{0.,1.}"):std::string("{0.,1.,0.}");
        //double sigmaf = trans(expr<FEELPP_DIM,1>(wall_normal))*FSImodel->fieldNormalStressFromFluidPtr_solid();
        //auto test = inner(expr<FEELPP_DIM,1>(wall_normal), idv(FSImodel->fieldNormalStressFromFluidPtr_solid()));
        typename model_fsi_type::element_solid_normalstressfromfluid_ptrtype sigmaf = FSImodel->fieldNormalStressFromFluidPtr_solid();
        contactForceModelsFSI<FEELPP_DIM,0,model_fsi_type>(t, sigmaf, data);
    };

    /*
    auto add_force_termR = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateResidual & data)
    {
        auto const& t = unwrap_ptr(FSImodel->solidModel());
        std::string wall_normal = (FEELPP_DIM==2)?std::string("{0.,1.}"):std::string("{0.,1.,0.}");
        //double sigmaf = trans(expr<FEELPP_DIM,1>(wall_normal))*FSImodel->fieldNormalStressFromFluidPtr_solid();
        //auto test = inner(expr<FEELPP_DIM,1>(wall_normal), idv(FSImodel->fieldNormalStressFromFluidPtr_solid()));
        typename model_fsi_type::element_solid_normalstressfromfluid_ptrtype sigmaf = FSImodel->fieldNormalStressFromFluidPtr_solid();
        contactForceModelsFSI<FEELPP_DIM,1,model_fsi_type>(t, sigmaf, data);
    };
    */


    FSImodel->solidModel()->algebraicFactory()->addFunctionLinearAssembly( add_force_term );
    //FSImodel->solidModel()->algebraicFactory()->addFunctionResidualAssembly( add_force_termR );
    // FSImodel->solidModel()->algebraicFactory()->addFunctionJacobianAssembly( std::bind( &self_type::updateLinearPDE_Solid,
    //                                                                             std::ref( *this ), std::placeholders::_1 ) );
    // FSImodel->solidModel()->algebraicFactory()->addFunctionResidualAssembly( std::bind( &self_type::updateLinearPDE_Solid,
    //                                                                             std::ref( *this ), std::placeholders::_1 ) );

    for ( FSImodel->startTimeStep() ; !FSImodel->timeStepBase()->isFinished(); FSImodel->updateTimeStep() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "\n====================================================================================="
                      << "\n current time : " << std::setprecision( 5 ) << std::fixed << FSImodel->currentTime()
                      << "\n=====================================================================================\n";
        FSImodel->solve();
        FSImodel->exportResults();
    }

}

} // namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

	po::options_description fsioptions( "application fsi options" );
    fsioptions.add( Feel::toolboxes_options("fsi") );
    fsioptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P2P1" ), "fe-approximation : P2P1,P2P1-P2 ")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=fsioptions,
                     _about=about(_name="application_fsi",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    std::string feapprox = soption(_name="fe-approximation");
#if FEELPP_GEO_ORDER == 1
    if ( feapprox == "P2P1" )
        runApplicationFSI<2,1>();
    else if ( feapprox == "P2P1-P2" )
        runApplicationFSI<2,1,2>();
#if 0//FEELPP_DIM == 2
    else if ( feapprox == "P1P1" )
        runApplicationFSI<1,1>();
#endif
#elif FEELPP_GEO_ORDER == 2
    if ( true )
        runApplicationFSI<2,1,2>();
#endif
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
