#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct MagnetoParam
    {
        std::string trajectory;
        double freq;
        double amp;
        double mx;
        double my;
        double mz;
        double bx;
        double by;
        double bz;
    };
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MagnetoParam,trajectory,freq,amp,mx,my,mz,bx,by,bz);
}

template<int nDim, std::size_t residualType, typename FSIModel, typename FluidMechanics, typename DataType>
void
magnetoTorqueModelFSI(FluidMechanics const& t, DataType & data)
{
    //std::cout << "apply fsi magneto torque" << std::endl;

    // Get parameters from JSON
    fs::path path (Environment::expand(soption(_name="fsi.filename")));
    json jsonMagneto;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i,nullptr,true,true);
        jsonMagneto = j["MagnetoTorque"]["body"]["setup"];
    }

    ns::MagnetoParam torqueParam = jsonMagneto["torqueParam"].get<ns::MagnetoParam>();
    std::string traj = torqueParam.trajectory;
    double freq = torqueParam.freq;
    double amp = torqueParam.amp;
    double bx = torqueParam.bx;
    double by = torqueParam.by;
    double bz = torqueParam.bz;
    double mx = torqueParam.mx;
    double my = torqueParam.my;
    double mz = torqueParam.mz;


    
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Get current orientation
    double orientation = 0;
    double x_curr = 0;
    
    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        auto angle = bpbc.body().rigidRotationAngles();

        if constexpr(nDim == 2)  
            orientation = angle(0,0);
        else if constexpr(nDim == 3)
            orientation = angle(2,0);

        auto massCenter = bpbc.body().massCenter();
        x_curr = massCenter(0,0);
        //std::cout << "Current orientation : " << orientation << std::endl;
    }
    
    // Compute true theta
    
    //double theta_true = amp * sin( 2 * M_PI * freq * t.currentTime() ); 
    //double T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * cst(by) * std::sin(theta_true) - cst(my) * std::sin(orientation) * cst(bx) * std::cos(theta_true)).evaluate()(0,0);
    double T_head = 0; 
    //std::cout << "Applied torue to head : " << T_head << std::endl;
    
    if (traj.compare("droite") == 0)
    {
        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) - cst(my) * std::sin(orientation) * cst(bx)).evaluate()(0,0);
    }
    else if (traj.compare("pipe") == 0)
    {
        T_head = 0;
    }
    else if (traj.compare("cos") == 0)
    {
        // n = (sin(x), 1), t = (1, - sin(x))

        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * (cst(bx) - cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) * std::sin(x_curr)) - cst(my) * std::sin(orientation) * (cst(bx)*std::sin(x_curr) + cst(by) * std::sin(2 * M_PI * freq * t.currentTime()))).evaluate()(0,0);
    }
    else if (traj.compare("circular") == 0)
    {
        double thetaTraj = (2 * M_PI * t.currentTime())/40.;
        T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * (cst(bx)*std::sin(thetaTraj) + cst(by) * std::sin(2 * M_PI * freq * t.currentTime()) * std::cos(thetaTraj)) - cst(my) * std::sin(orientation) * (cst(bx)*std::cos(thetaTraj) - cst(by) * std::sin(2 * M_PI * freq * t.currentTime())*std::sin(thetaTraj))).evaluate()(0,0);
    }

    // Add torque to newton eq
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };
    
    auto rowStartInVector = t.rowStartInVector();
    
    // Add torque
    for (auto const& [bpname,bpbc] : t.bodySetBC())
    {
        size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
        int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();
    
        if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
        {
            auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);
    
            if (residualType == 1) 
            {
                if constexpr(nDim == 2)  
                    r()->add(basisToContainerGpAngularVelocityVector[0],-T_head);
                else if constexpr(nDim == 3)  
                    r()->add(basisToContainerGpAngularVelocityVector[2],-T_head);
            }
                  
            else if (residualType == 0)
            {
                if constexpr(nDim == 2)  
                    r()->add(basisToContainerGpAngularVelocityVector[0],T_head);
                else if constexpr(nDim == 3)  
                    r()->add(basisToContainerGpAngularVelocityVector[2],T_head);
            }
        }       
                     
    }    
}




// namespace Feel
// {
 
// template <typename FluidStructureInteraction>
// void
// runApplicationFSI_magneto(int &nDim, int OrderT, int &OrderVelocity, int &OrderPressure, int &OrderGeo)
// {
//     using namespace Feel;
 
//     typedef FeelModels::FluidMechanics< Simplex<nDim,OrderGeo>,
//                                         Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
//                                         Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_fluid_type;
    
//     typedef FeelModels::SolidMechanics< Simplex<nDim,OrderGeo>,
//                                         Lagrange<OrderDisp, Vectorial,Continuous,PointSetFekete> > model_solid_type;
    
//     typedef FeelModels::FSI< model_fluid_type,model_solid_type> model_fsi_type;
    
//     std::shared_ptr<model_fsi_type> FSImodel( new model_fsi_type("fsi") );
 
//     FSImodel->init();
//     FSImodel->printAndSaveInfo();

//     // Add magneto torque to fluid-rigid interaction
//     auto add_torque = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateLinear & data)
//     {
//         auto const& t = unwrap_ptr(FSImodel->fluidModel());
//         magnetoTorqueModelFSI<nDim,0,model_fsi_type>(t, data);
//     };
//     // add the lambda function to the algebraic factory
//     FSImodel->fluidModel()->algebraicFactory()->addFunctionLinearAssembly( add_torque );
//     auto add_torque_residual = [&FSImodel](FeelModels::ModelAlgebraic::DataUpdateResidual & data)
//                                    {
//                                        auto const& t = unwrap_ptr(FSImodel->fluidModel());
//                                        magnetoTorqueModelFSI<nDim,1,model_fsi_type>(t, data);
//                                    };
//     // add the lambda function to the algebraic factory
//     FSImodel->fluidModel()->algebraicFactory()->addFunctionResidualAssembly( add_torque_residual );
  
// #if 0
//     auto eFluid = exporter( _mesh=FSImodel->fieldNormalStressRefMeshPtr_fluid()->mesh(), _name="ExportInterfaceFSI_fluid", _geo="static" );
//     auto eSolid = exporter( _mesh=FSImodel->fieldNormalStressFromFluidPtr_solid()->mesh(), _name="ExportInterfaceFSI_solid", _geo="static" );

//     auto const& bbc = FSImodel->fluidModel()->bodySetBC().begin()->second;
//     auto const& body = bbc.body();
//     auto eBody = exporter( _mesh=body.fieldDisplacement().mesh() , _name="ExportBody", _geo="static" );
// #endif
//     for ( FSImodel->startTimeStep() ; !FSImodel->timeStepBase()->isFinished(); FSImodel->updateTimeStep() )
//     {
//         if ( Environment::isMasterRank() )
//             std::cout << "\n====================================================================================="
//                       << "\n current time : " << std::setprecision( 5 ) << std::fixed << FSImodel->currentTime()
//                       << "\n=====================================================================================\n";
//         FSImodel->solve();
//         FSImodel->exportResults();
// #if 0
//         eFluid->step(FSImodel->currentTime())->add( "normal_stress_fluid", *FSImodel->fieldNormalStressRefMeshPtr_fluid() );
//         eFluid->save();
//         eSolid->step(FSImodel->currentTime())->add( "normal_stress_solid", *FSImodel->fieldNormalStressFromFluidPtr_solid() );
//         eSolid->save();
// #endif
// #if 0
//         auto rangeBody = elements(support(bbc.body().fieldDisplacement().functionSpace()));
//         //eBody->step(FSImodel->currentTime())->add( "disp", rangeBody, idv(bbc.body().fieldDisplacement()) );
//         eBody->step(FSImodel->currentTime())->add( "disp", bbc.body().fieldDisplacement() );
//         eBody->step(FSImodel->currentTime())->add( "elast-disp", bbc.body().fieldElasticDisplacement() );
//         eBody->step(FSImodel->currentTime())->add( "elast-vel", bbc.body().fieldElasticVelocity() );
//         eBody->step(FSImodel->currentTime())->add( "rigidVelocity", bbc.rigidVelocityExprFromFields() );
//         eBody->save();
// #endif
//         //t->updateDisplacementImposedOnInitialDomain( this->keyword()+"_body", idv(bbc.body().fieldDisplacement()), elements(support(bbc.body().fieldDisplacement().functionSpace())) );
//     }
// }
 
// } // namespace Feel

