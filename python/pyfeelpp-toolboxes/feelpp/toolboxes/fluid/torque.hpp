
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct RigidTorqueParam
    {
        double T1;
        double T2;
        double T3;
    };

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RigidTorqueParam,T1, T2, T3);
}


template<std::size_t residualType,typename FluidMechanics, typename DataType>
void
TorqueRigid(FluidMechanics &t, DataType & data, ns::RigidTorqueParam rigidTorqueParam)
{
    int const dim = t.nDim;



    // Compute torque
    double T1 = rigidTorqueParam.T1;
    double T2 = rigidTorqueParam.T2;
    double T3;

    if (dim == 3)
    {
        T3 = rigidTorqueParam.T3;
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

              if (dim == 2)
            {
                std::cout << "Torque : " << T1 << " " << T2 << std::endl;
            }
            else if (dim == 3)
            {
                std::cout << "Torque : " << T1 << " " << T2 << " " << T3 << std::endl;
            }


            if (residualType == 1) 
            {
                r()->add(basisToContainerGpAngularVelocityVector[0],-T1);
                r()->add(basisToContainerGpAngularVelocityVector[1],-T2);
                if (dim == 3)
                    r()->add(basisToContainerGpAngularVelocityVector[2],-T3);
            }
            else if (residualType == 0)
            {
                r()->add(basisToContainerGpAngularVelocityVector[0],T1);
                r()->add(basisToContainerGpAngularVelocityVector[1],T2);
                if (dim == 3)
                    r()->add(basisToContainerGpAngularVelocityVector[2],T3);
            }
        }       
                
    }
    
}




template<std::size_t residualType, typename FluidMechanics, typename DataType>
void
rigidTorqueModel(FluidMechanics &t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Read json
    fs::path path (Environment::expand(soption(_name="pfluid.filename")));
    json jsonRigid;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonRigid = j["RigidTorque"]["body"]["setup"];
    }

    ns::RigidTorqueParam rigidTorqueParam = jsonRigid["torqueParam"].get<ns::RigidTorqueParam>();

    TorqueRigid<residualType>(t, data, rigidTorqueParam);

}
