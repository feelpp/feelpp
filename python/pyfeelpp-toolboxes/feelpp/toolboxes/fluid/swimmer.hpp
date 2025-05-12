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
        double freq;
        double amp;
        double mx;
        double my;
        double mz;
        double bx;
        double by;
        double bz;
    };
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MagnetoParam,freq,amp,mx,my,mz,bx,by,bz);
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
    
    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        auto angle = bpbc.body().rigidRotationAngles();
        orientation = angle(0,0);
        //std::cout << "Current orientation : " << orientation << std::endl;
    }
    
    // Compute true theta
    double theta_true = amp * sin( 2 * M_PI * freq * t.currentTime() );
       
    double T_head = integrate( _range = markedelements( t.mesh(), "Head" ), _expr = cst(mx) * std::cos(orientation) * cst(by) * std::sin(theta_true) - cst(my) * std::sin(orientation) * cst(bx) * std::cos(theta_true)).evaluate()(0,0);
    //std::cout << "Applied torue to head : " << T_head << std::endl;
    
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
                r()->add(basisToContainerGpAngularVelocityVector[0],-T_head);  
            else if (residualType == 0)
                r()->add(basisToContainerGpAngularVelocityVector[0],T_head); 
        }       
                     
    }    
}
