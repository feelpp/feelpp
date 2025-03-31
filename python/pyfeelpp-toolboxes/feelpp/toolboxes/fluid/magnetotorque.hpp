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
        double mx;
        double my;
        double mz;
        double bx;
        double by;
        double bz;
    };

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MagnetoParam,trajectory,mx,my,mz,bx,by,bz);
}

class Data
{
    public:
        static int nbr;
        static int iter;
        static std::vector<double> thetas;
        static std::vector<double> torques;
};

int Data::nbr = 0;
int Data::iter = 1;
std::vector<double> Data::thetas = {0.0};  
std::vector<double> Data::torques = {0.0}; 

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
MagnetoRigid(FluidMechanics &t, DataType & data, ns::MagnetoParam torqueParam)
{
    int const dim = t.nDim;

    if (dim == 2)
    {
        // Get collision parameters    
        std::string trajectory = torqueParam.trajectory;
        double bx = torqueParam.bx;
        double by = torqueParam.by;
        double bz = torqueParam.bz;
        double mx = torqueParam.mx;
        double my = torqueParam.my;
        double mz = torqueParam.mz;

        // Get current orientation
        double orientation;

        for ( auto const& [bpname,bpbc] : t.bodySetBC() )
        {
            auto angle = bpbc.body().rigidRotationAngles();
            orientation = angle(0,0);
            std::cout << "Orientation : " << orientation << std::endl;
            
            t.addParameterInModelProperties("theta", orientation);  
            
        }

        // Compute true theta
        double theta_true = 0;
        if (trajectory.compare("droite") == 0)
        { 
            theta_true = 0; 
        }

        if (trajectory.compare("cos") == 0)
        {
            double timestep = t.timeStep();
            double x0 = t.time();
            double x1 = t.time() + timestep;

            double y0 = std::cos(2 * M_PI * x0); 
            double y1 = std::cos(2 * M_PI * x1); 

            double vel = std::sqrt(timestep*timestep + (y1-y0)*(y1-y0))/timestep;
            std::cout << "vel : " << vel << std::endl;
            t.addParameterInModelProperties("vel", vel);  

            double dy_dx = - 2 * M_PI * std::sin(2 * M_PI * t.time());  
            theta_true = std::atan(dy_dx); 
        }
        t.updateParameterValues();
        std::cout << "Theta true : " << theta_true << std::endl;
        // Compute torque
        double T = mx * std::cos(orientation) * by * std::sin(theta_true) - my * std::sin(orientation) * bx * std::cos(theta_true);

        // Add torque to newton eq
        auto r = [&data]() 
        { 
            if constexpr(residualType == 1) 
                return data.residual(); 
            else if constexpr(residualType == 0)
                return data.rhs();
        };

        auto rowStartInVector = t.rowStartInVector();

        for (auto const& [bpname,bpbc] : t.bodySetBC())
        {
            size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();

            if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
            {
                auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);

                if (residualType == 1) 
                    r()->add(basisToContainerGpAngularVelocityVector[0],-T);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpAngularVelocityVector[0],T); 
            }       
                 
        }

        // Calcul du nombre de Reynlds
        auto Xh = Pch<1>(t.mesh());
        auto Re = Xh->element();
        
        fs::path path (Environment::expand(soption(_name="fluid.filename")));
        double mu;
    
        if (fs::exists(path))
        {
            std::ifstream i(path.string().c_str());
            json j = json::parse(i);
            mu = std::stod(j["Materials"]["Fluid"]["mu"].get<std::string>());
        }

        Re = project(_space = Xh, _range = elements(t.mesh()), _expr = norm2( idv(t.fieldVelocity()) )/cst(mu));
        t.modelMesh().template updateField<typename FluidMechanics::mesh_type>( "Re", idv(Re), elements(t.mesh()), "Pch1" );
    

        // add data
        Data::nbr += 1;
        Data::thetas.push_back(orientation);
        Data::torques.push_back(T);

    }
}

template<typename FluidMechanics>
void
reset_Data(FluidMechanics &t)
{
    Data::nbr = 0;
    Data::iter = 1;
    Data::thetas = {0.0};
    Data::torques = {0.0};
}

template<typename FluidMechanics>
void
write_Data(FluidMechanics &t)
{
    std::ofstream file("data.csv"); 

    // Écriture des en-têtes
    file << "theta,torque\n";
   
    for (int i = 1; i < Data::nbr; i ++)
    {
        file << Data::thetas[i] << "," << Data::torques[i] << "\n";
    }
    file.close();
}


template<std::size_t residualType, typename FluidMechanics, typename DataType>
void
magnetoTorqueModel(FluidMechanics &t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Read json
    fs::path path (Environment::expand(soption(_name="fluid.filename")));
    json jsonMagneto;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonMagneto = j["MagnetoTorque"]["body"]["setup"];
    }

    ns::MagnetoParam torqueParam = jsonMagneto["torqueParam"].get<ns::MagnetoParam>();

    MagnetoRigid<residualType>(t, data, torqueParam);

}