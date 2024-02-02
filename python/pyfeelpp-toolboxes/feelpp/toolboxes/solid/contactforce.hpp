#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/solid/solidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/matvec.hpp>
#include <mpi.h>
#include <feel/feelmesh/meshsupportbase.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/localization.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/detail/ginacmatrix.hpp>
#include <fmt/ostream.h>
#include <feel/feelmesh/bvh.hpp>

#include "contactforceStatic.hpp"
#include "contactforceDynamic.hpp"

using namespace Feel;
using namespace Feel::FeelModels;
using index_type = uint32_type;

template<int nDim, std::size_t residualType, typename SolidMechanics, typename DataType>
void
contactForceModels(SolidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    // Get parameters form json
    std::string model = jsonCollisionForce["model"].get<std::string>(); 
    std::string type = jsonCollisionForce["type"].get<std::string>();
    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
    double theta =  jsonCollisionForce["theta"].get<double>();
    int setA = jsonCollisionForce["updateA"].get<int>();
    std::vector<double> pressurePoint = jsonCollisionForce["pressurePoint"].get<std::vector<double>>();
    int withMarker = jsonCollisionForce["withMarker"].get<int>();
    double tolerance = jsonCollisionForce["tolerance"].get<double>();
    int fixedPoint = jsonCollisionForce["fixedPoint"].get<int>();            
    int dispCondition = jsonCollisionForce["dispCondition"].get<int>();
    std::vector<double> initCondition = jsonCollisionForce["initCondition"].get<std::vector<double>>();
    int raytracing = jsonCollisionForce["raytracing"].get<int>();
    double distance = 0;
    distance = jsonCollisionForce["walls"]["distances"].get<double>(); 
    int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

    // Init distance
    if (t.timeStepBase()->iteration() == 1)
        initG(t, direction, initCondition);

    // Model : steady
    if (model.compare("steady") == 0)
    {
        // Type case : unilateral
        if (type.compare("unilateral") == 0)
        {

            if (nbr_walls == 1)
            {                
                std::cout << " Steady collision algorithm " << std::endl;
                SingleUnilateralSteady<nDim, residualType>(t, data, direction, gamma0, theta, setA, pressurePoint, dispCondition, withMarker, raytracing, distance);                 
            }   
            else if (nbr_walls > 1)
            {
                // TODO
            }        
        }
        else if (type.compare("bilateral") == 0)
        {
            // TODO
        }
    }    

    // Model : dynamic
    if (model.compare("dynamic") == 0)
    {
        if (type.compare("unilateral") == 0)
        {
            if (nbr_walls == 1)
            {                
                if (fixedPoint == 1)
                {
                    std::cout << " Dynamic Fixed point collision algorithm " << std::endl;
                    SingleUnilateralDynamicFixedPoint<nDim, residualType>(t, data, direction, gamma0, theta, tolerance, setA,  dispCondition, withMarker, raytracing, distance);
                }
                else 
                {
                    std::cout << " Dynamic collision algorithm " << std::endl;
                    SingleUnilateralDynamic<nDim, residualType>(t, data, direction, gamma0, theta,setA, dispCondition, withMarker, raytracing, distance);
                }
                
            }   
            else if (nbr_walls > 1)
            {
                // TODO
            }

        }
        else if (type.compare("bilateral") == 0)
        {
            // TODO
        }
    }
}

