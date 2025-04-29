#include <feel/feelmodels/fsi/fsi.hpp>
#include "swimmer.hpp"

namespace Feel
{
 
template <uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderDisp=FEELPP_GEO_ORDER>
void
runApplicationFSI_magneto()
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

    // Add magneto torque to fluid-rigid interaction
    auto add_torque = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateLinear & data)
    {
        auto const& t = unwrap_ptr(FSImodel->fluidModel());
        magnetoTorqueModelFSI<FEELPP_DIM,0,model_fsi_type>(t, data);
    };
    // add the lambda function to the algebraic factory
    FSImodel->fluidModel()->algebraicFactory()->addFunctionLinearAssembly( add_torque );
  
    for ( FSImodel->startTimeStep() ; !FSImodel->timeStepBase()->isFinished(); FSImodel->updateTimeStep() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "\n====================================================================================="
                      << "\n current time : " << std::setprecision( 5 ) << std::fixed << FSImodel->currentTime()
                      << "\n=====================================================================================\n";
        FSImodel->solveMagneto();
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

    runApplicationFSI_magneto<2,1>(); 
    return 0;
}
 