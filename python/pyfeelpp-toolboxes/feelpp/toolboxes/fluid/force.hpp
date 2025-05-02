
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct RigidForceParam
    {
        double F1;
        double F2;
        double F3;
    };

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(RigidForceParam,F1, F2, F3);
}


template<std::size_t residualType,typename FluidMechanics, typename DataType>
void
ForceRigid(FluidMechanics &t, DataType & data, ns::RigidForceParam rigidForceParam)
{
    int const dim = t.nDim;
        

    // Compute force
    double F1 = rigidForceParam.F1;
    double F2 = rigidForceParam.F2;
    double F3;

    if (dim == 3)
    {
        F3 = rigidForceParam.F3;
    }
  

    // Add force to newton eq
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };

    auto rowStartInVector = t.rowStartInVector();


    // Add force
    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
        r()->setIsClosed(false);
        if ( bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0 )
        {
            auto const& basisToContainerGpTranslationalVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexTranslationalVelocity);   

            if (dim == 2)
            {
                std::cout << "Force : " << F1 << " " << F2 << std::endl;
            }
            else if (dim == 3)
            {
                std::cout << "Force : " << F1 << " " << F2 << " " << F3 << std::endl;
            }


            if (residualType == 1) 
            {
                r()->add(basisToContainerGpTranslationalVelocityVector[0],-F1);
                r()->add(basisToContainerGpTranslationalVelocityVector[1],-F2);
                if (dim == 3)
                    r()->add(basisToContainerGpTranslationalVelocityVector[2],-F3);
            }
            else if (residualType == 0)
            {
                r()->add(basisToContainerGpTranslationalVelocityVector[0],F1);
                r()->add(basisToContainerGpTranslationalVelocityVector[1],F2);
                if (dim == 3)
                    r()->add(basisToContainerGpTranslationalVelocityVector[2],F3);
            } 
    
          
        }              
    } 
    
}




template<std::size_t residualType, typename FluidMechanics, typename DataType>
void
rigidForceModel(FluidMechanics &t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Read json
    fs::path path (Environment::expand(soption(_name="dfluid.filename")));
    json jsonRigid;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonRigid = j["RigidForce"]["body"]["setup"];
    }

    ns::RigidForceParam rigidForceParam = jsonRigid["forceParam"].get<ns::RigidForceParam>();

    ForceRigid<residualType>(t, data, rigidForceParam);

}
