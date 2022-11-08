#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 
#include <feel/feelfilters/exporter.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feelvf/vf.hpp>
#include <mpi.h>

using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct CollisionForceParam
    {
        double forceRange;
        double epsBody;
        double epsWall;
    };

    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CollisionForceParam,forceRange,epsBody,epsWall);
}

struct Distance_param
{
    int id_A;
    int id_B;
    RowVectord coord_A;
    RowVectord coord_B;
    double dist;        
};

/*
    Implementation of ContactAvoidance model for spherical shaped bodies
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceSphericalType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time, double &time_pre, double &time_detection, double &time_com)
{
    tic();
    
    // Define force parameters

    int const dim = t.nDim;
    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;  
    
    // Pre-process phase

    auto start_pre = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> massCenters;
    std::vector<std::string> marker_bodies;
    std::vector<double> radii;
    int nbr_bodies = 0;
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        std::string marker = bpbc.name();
        marker_bodies.push_back(marker);
        massCenters.push_back(bpbc.body().massCenter());
        double meas = integrate( _range=markedfaces(t.mesh(),marker), _expr=cst(1.0)).evaluate()(0,0);
        if (dim == 2)
            radii.push_back(meas/(2*Pi));
        else if (dim == 3)
            radii.push_back(sqrt(meas/(4*Pi)));
        nbr_bodies ++;
    }

    for (auto m : t.mesh()->markerNames())
    {
        if (dim == 2)
        {
            if (m.second[1] == 1 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
        else if (dim == 3)
        {
            if (m.second[1] == 2 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
    }

    auto end_pre = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_pre = end_pre - start_pre;
    time_pre = float_pre.count();

    // Collision detection phase

    auto start_detection = std::chrono::high_resolution_clock::now();

    auto Vh = Pch<1>(t.mesh()); 
    std::vector<decltype( Vh->element() )> distToFluid;

    for (int b_w = 0;b_w < nbr_fluid; b_w++)
    {
        distToFluid.push_back(distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange));
    }

    // Define collision map
    std::map<std::string,Distance_param> DistanceParam;

    for (int b_i = 0;b_i < nbr_bodies; b_i++)
    {
        for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
        {
            double dist_ij = sqrt((massCenters[b_j]-massCenters[b_i]).squaredNorm()) -radii[b_i] - radii[b_j];
            if (dist_ij <= forceRange)
            {
                std::string id = marker_bodies[b_i] + "_" + marker_bodies[b_j];
                Distance_param d;
                d.id_A = b_i;
                d.id_B = b_j;
                d.coord_A = massCenters[b_i];
                d.coord_B = massCenters[b_j];
                d.dist = dist_ij;
                DistanceParam[id] = d;
            }
               
        }

        for (int b_w = 0;b_w < nbr_fluid; b_w++)
        {
            auto [minToBi,arg_minToBi] = minelt(_range=markedfaces(t.mesh(),marker_bodies[b_i]), _element=distToFluid[b_w]);

            if (minToBi <= forceRange)
            {    
                std::string id = marker_bodies[b_i] + "_" + marker_fluid[b_w];

                Distance_param d;
                d.id_A = b_i;
                d.id_B = -1; // Use id -1 for wall boundary
                d.coord_A = arg_minToBi;
                d.dist = minToBi;

                auto distToBi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange);                
                auto [minToBw,arg_minToBw] = minelt(_range=markedfaces(t.mesh(),marker_fluid[b_w]), _element=distToBi);
                d.coord_B = arg_minToBw;

                DistanceParam[id] = d;
            }
        }
    }


    auto end_detection = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_detection = end_detection - start_detection;
    time_detection =  float_detection.count();

    // Compute and insert lubrication force
    
    auto start_com = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));
           
    for (const auto& [key, value] : DistanceParam)
    {
        if (value.id_B != -1)
        {          
            RowVectord F_ij = 1/epsBody * std::pow(forceRange - value.dist,2)*(value.coord_A - value.coord_B);
            repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;
            repulsion_forces[value.id_B] = repulsion_forces[value.id_B] - F_ij;           
        }

        if (value.id_B == -1)
        {
            RowVectord F_ij = 1/epsWall * std::pow(forceRange - value.dist,2)*(value.coord_A - value.coord_B);
            repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;
        }
    }
    
    // Define rhs or residual                       
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };
    
    // Add the repulsion force 
    int B = 0;
    auto rowStartInVector = t.rowStartInVector();

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
        r()->setIsClosed(false);
        if ( bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0 )
        {
            auto const& basisToContainerGpTranslationalVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexTranslationalVelocity);    
        
            for (int d=0;d<dim;++d)
            {   
                if (residualType == 1) 
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],-repulsion_forces[B][d]);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);      
            }
        }
        B++;               
    }    

    auto end_com = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_com = end_com - start_com;
    time_com += float_com.count();

    time = toc("ContactAvoidanceSphericalType");
}


/*
    Implementation of ContactAvoidance model for complex shaped bodies
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceComplexType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time, double &time_pre, double &time_detection, double &time_com)
{
    tic();

    // Get collision parameters    

    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;
    int const dim = t.nDim;
    
    // Pre-process phase

    auto start_pre = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> massCenters;
    std::vector<std::string> marker_bodies;
    int nbr_bodies = 0;
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        marker_bodies.push_back(bpbc.name());
        massCenters.push_back(bpbc.body().massCenter());
        nbr_bodies ++;

    }

    for (auto m : t.mesh()->markerNames())
    {
        if (dim == 2)
        {
            if (m.second[1] == 1 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
        else if (dim == 3)
        {
            if (m.second[1] == 2 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
    }

    auto end_pre = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_pre = end_pre - start_pre;
    time_pre = float_pre.count();
   
    // Collision detection phase

    auto start_detection = std::chrono::high_resolution_clock::now();

    auto Vh = Pch<1>(t.mesh()); 
    std::vector<decltype( Vh->element() )> distToBodies;
    std::vector<decltype( Vh->element() )> distToFluid;

    for (int b_i = 0; b_i < nbr_bodies; b_i++)
    {
        distToBodies.push_back(distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange));
    }

    for (int b_w = 0;b_w < nbr_fluid; b_w++)
    {
        distToFluid.push_back(distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange));
    }

    // Define collision map
    std::map<std::string,Distance_param> DistanceParam;

    for (int b_i = 0; b_i < nbr_bodies; b_i++)
    {
        for (int b_j = b_i+1;b_j < nbr_bodies; b_j++)
        {
            auto [minToBj,arg_minToBj] = minelt(_range=markedfaces(t.mesh(),marker_bodies[b_j]), _element=distToBodies[b_i]);
                
            if (minToBj <= forceRange)
            {

                std::string id = marker_bodies[b_i] + "_" + marker_bodies[b_j];

                Distance_param d;
                d.id_A = b_i;
                d.id_B = b_j;
                d.coord_B = arg_minToBj;
                d.dist = minToBj;

                auto [minToBi,arg_minToBi] = minelt(_range=markedfaces(t.mesh(),marker_bodies[b_i]), _element=distToBodies[b_j]);
                d.coord_A = arg_minToBi;

                DistanceParam[id] = d;
            }
               
        }

        for (int b_w = 0;b_w < nbr_fluid; b_w++)
        {
            auto [minToBw,arg_minToBw] = minelt(_range=markedfaces(t.mesh(),marker_fluid[b_w]), _element=distToBodies[b_i]);

            if (minToBw <= forceRange)
            {    
                std::string id = marker_bodies[b_i] + "_" + marker_fluid[b_w];

                Distance_param d;
                d.id_A = b_i;
                d.id_B = -1; // Use id -1 for wall boundary
                d.coord_B = arg_minToBw;
                d.dist = minToBw;

                auto [minToBi,arg_minToBi] = minelt(_range=markedfaces(t.mesh(),marker_bodies[b_i]), _element=distToFluid[b_w]);
                d.coord_A = arg_minToBi;

                DistanceParam[id] = d;
            }
        }
    }

    auto end_detection = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_detection = end_detection - start_detection;
    time_detection =  float_detection.count();
    
    // Compute and insert lubrication force
    
    auto start_com = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));

    std::vector<Eigen::Vector3d> rot(nbr_bodies);
    std::fill(rot.begin(), rot.end(), Eigen::Vector3d::Zero());

    for (const auto& [key, value] : DistanceParam)
    {
        if (value.id_B != -1)
        {          
            RowVectord F_ij = 1/epsBody * std::pow(forceRange - value.dist,2)*(value.coord_A - value.coord_B);
            repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;
            repulsion_forces[value.id_B] = repulsion_forces[value.id_B] - F_ij; 
            
            if (dim == 2)
            {
                // We define Eigen::Vector3d types to use the function cross()    
                Eigen::Vector3d x_diff_F_i(value.coord_A[0] - massCenters[value.id_A][0],value.coord_A[1] - massCenters[value.id_A][1],0.);
                Eigen::Vector3d x_diff_F_j(value.coord_B[0] - massCenters[value.id_B][0],value.coord_B[1] - massCenters[value.id_B][1],0.);
                Eigen::Vector3d repulsion_force_i(F_ij[0],F_ij[1],0.);
                Eigen::Vector3d repulsion_force_j(-F_ij[0],-F_ij[1],0.);
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);
                rot[value.id_B] = rot[value.id_B] - repulsion_force_j.cross(x_diff_F_j);
            }    
            else if (dim == 3)
            {
                Eigen::Vector3d x_diff_F_i = value.coord_A - massCenters[value.id_A]; 
                Eigen::Vector3d x_diff_F_j = value.coord_B - massCenters[value.id_B]; 
                Eigen::Vector3d repulsion_force_i = F_ij;
                Eigen::Vector3d repulsion_force_j = -F_ij;
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);
                rot[value.id_B] = rot[value.id_B] - repulsion_force_j.cross(x_diff_F_j);                     
            }     
        }

        if (value.id_B == -1)
        {
            RowVectord F_ij = 1/epsWall * std::pow(forceRange - value.dist,2)*(value.coord_A - value.coord_B);
            repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;

            if (dim == 2)
            {   
                Eigen::Vector3d x_diff_F_i(value.coord_A[0] - massCenters[value.id_A][0],value.coord_A[1] - massCenters[value.id_A][1],0.);
                Eigen::Vector3d repulsion_force_i(F_ij[0],F_ij[1],0.);
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);
            }    
            else if (dim == 3)
            {
                Eigen::Vector3d x_diff_F_i = value.coord_A - massCenters[value.id_A]; 
                Eigen::Vector3d repulsion_force_i = F_ij;
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);                  
            } 
        }
    }
                          
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };

    int B = 0;
    auto rowStartInVector = t.rowStartInVector();

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
        r()->setIsClosed(false);
        if ( bpbc.spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0 )
        {
            auto const& basisToContainerGpTranslationalVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexTranslationalVelocity);    
        
            for (int d=0;d<dim;++d)
            {   
    
                if (residualType == 1) 
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],-repulsion_forces[B][d]);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);  
                       
            }
        }
        B++;               
    }
    
    B = 0;
    for (auto const& [bpname,bpbc] : t.bodySetBC())
    {
        size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
        int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();

        if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
        {
            auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);


            if (dim == 2)
            {
                if (residualType == 1) 
                    r()->add(basisToContainerGpAngularVelocityVector[0],-rot[B][2]);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpAngularVelocityVector[0],rot[B][2]); 
            }        
            else if (dim == 3)
            {
                for (int d=0;d<nLocalDofAngularVelocity;++d)
                {       
                    if (residualType == 1) 
                        r()->add(basisToContainerGpAngularVelocityVector[d],-rot[B][d]);  
                    else if (residualType == 0)
                        r()->add(basisToContainerGpAngularVelocityVector[d],rot[B][d]);     
                }         
            }     
        }
        B++;
    }

    auto end_com = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_com = end_com - start_com;
    time_com += float_com.count();

    time = toc("ContactAvoidanceComplexType");
}

/*
  Implementation of ContactAvoidance for articulated bodies
  Todo : interaction between articulated bodies
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceArticulatedType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time, double &time_pre, double &time_detection, double &time_com)
{
    tic();

    // Define collision parameters 
    int const dim = t.nDim;
    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;
    
    // Pre-process phase

    auto start_pre = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> massCenters;
    std::vector<RowVectord> G;
    std::vector<std::string> marker_bodies;
    std::vector<RowVectord> GravityForce;
    int nbr_bodies = 0;
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    for (auto const&nba : t.bodySetBC().nbodyArticulated())
    {
        auto gravityCenter = nba.massCenter();
        for (auto bbc : nba.bodyList())
        {
            G.push_back(gravityCenter);
            marker_bodies.push_back(bbc->name());
            massCenters.push_back(bbc->body().massCenter());
         
            if (bbc->gravityForceEnabled())
            {
                GravityForce.push_back(bbc->gravityForceWithMass());
            }
            else 
            {
                GravityForce.push_back(RowVectord::Zero(1));
            }

            nbr_bodies ++;
        }
    }

    for (auto m : t.mesh()->markerNames())
    {
        if (dim == 2)
        {
            if (m.second[1] == 1 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
        else if (dim == 3)
        {
            if (m.second[1] == 2 && std::find(marker_bodies.begin(), marker_bodies.end(), m.first) == marker_bodies.end())
            {
                marker_fluid.push_back(m.first);
                nbr_fluid += 1;
            }    
        }
    }

    auto end_pre = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_pre = end_pre - start_pre;
    time_pre = float_pre.count();
   
    // Collision detection phase

    auto start_detection = std::chrono::high_resolution_clock::now();
  
    auto Vh = Pch<1>(t.mesh()); 
    std::vector<decltype( Vh->element() )> distToFluid;

    for (int b_w = 0; b_w < nbr_fluid; b_w++)
    {
        distToFluid.push_back(distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange));
    }

    // Define collision map
    std::map<std::string,Distance_param> DistanceParam;

    for (int b_w = 0; b_w < nbr_fluid; b_w++)
    {
        for (int b_i = 0; b_i < nbr_bodies; b_i++)
        {
            auto [minToBi,arg_minToBi] = minelt(_range=markedfaces(t.mesh(),marker_bodies[b_i]), _element=distToFluid[b_w]);
                
            if (minToBi <= forceRange)
            {

                std::string id = marker_bodies[b_i] + "_" + marker_fluid[b_w];

                Distance_param d;
                d.id_A = b_i;
                d.id_B = -1; // Use id -1 for wall boundary
                d.coord_A = arg_minToBi;
                d.dist = minToBi;

                auto distToBi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange);
                auto [minToBw,arg_minToBw] = minelt(_range=markedfaces(t.mesh(),marker_fluid[b_w]), _element=distToBi);
                d.coord_B = arg_minToBw;

                DistanceParam[id] = d;
            }
               
        }

    }

    auto end_detection = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_detection = end_detection - start_detection;
    time_detection =  float_detection.count();

    // Compute and insert lubrication force
    
    auto start_com = std::chrono::high_resolution_clock::now();

    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));

    std::vector<Eigen::Vector3d> rot(nbr_bodies);
    std::fill(rot.begin(), rot.end(), Eigen::Vector3d::Zero());

    for (const auto& [key, value] : DistanceParam)
    {
        if (value.id_B == -1)
        {
            RowVectord F_ij = 1/epsWall * std::pow(forceRange - value.dist,2)*(value.coord_A - value.coord_B);
            repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;

            if (dim == 2)
            {   
                Eigen::Vector3d x_diff_F_i(value.coord_A[0] - G[value.id_A][0],value.coord_A[1] - G[value.id_A][1],0.);
                Eigen::Vector3d repulsion_force_i(F_ij[0],F_ij[1],0.);
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);
            }    
            else if (dim == 3)
            {
                Eigen::Vector3d x_diff_F_i = value.coord_A - G[value.id_A]; 
                Eigen::Vector3d repulsion_force_i = F_ij;
                rot[value.id_A] = rot[value.id_A] - repulsion_force_i.cross(x_diff_F_i);                  
            } 
        }
    }

    for (int b_i = 0; b_i < nbr_bodies; b_i++)
    {
        if (dim == 2)
        {
            Eigen::Vector3d x_diff_G_i(massCenters[b_i][0] - G[b_i][0],massCenters[b_i][1] - G[b_i][1],0.);
            Eigen::Vector3d gravity_force_i(GravityForce[b_i][0],GravityForce[b_i][1],0.);
            rot[b_i] = rot[b_i] + gravity_force_i.cross(x_diff_G_i);

        }        
        else if (dim == 3)
        {
            Eigen::Vector3d x_diff_G_i = massCenters[b_i] - G[b_i];
            Eigen::Vector3d gravity_force_i = GravityForce[b_i];
            rot[b_i] = rot[b_i] + gravity_force_i.cross(x_diff_G_i);                  
        }     
    }


    // Define rhs or residual                       
    auto r = [&data]() 
    { 
        if constexpr(residualType == 1) 
            return data.residual(); 
        else if constexpr(residualType == 0)
            return data.rhs();
    };

    // Add repulsion force
    int B = 0;
    auto rowStartInVector = t.rowStartInVector();

    for (auto const&nba : t.bodySetBC().nbodyArticulated())
    {
        for (auto bbc : nba.bodyList())
        {
            size_type startBlockIndexTranslationalVelocity = t.startSubBlockSpaceIndex("body-bc."+bbc->name()+".translational-velocity");
            r()->setIsClosed(false);

            if ( bbc->spaceTranslationalVelocity()->nLocalDofWithoutGhost() > 0 )
            {
                auto const& basisToContainerGpTranslationalVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexTranslationalVelocity);    
        
                for (int d=0;d<dim;++d)
                {   
                    if (residualType == 1) 
                        r()->add(basisToContainerGpTranslationalVelocityVector[d],-repulsion_forces[B][d]);  
                    else if (residualType == 0)
                        r()->add(basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);          
                }
            }
            B++; 
        }              
    }   

    // Add external forces on rotation
    for (auto const& [bpname,bpbc] : t.bodySetBC())
    {
        if ((bpbc.getNBodyArticulated().masterBodyBC().name() == bpbc.name()))
        {
            size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();

            if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
            {
                auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);
                
                for (int B=0;B<nbr_bodies;B++)
                {  
                    if (dim == 2)
                    {
                        if (residualType == 1) 
                            r()->add(basisToContainerGpAngularVelocityVector[0],-rot[B][2]);
                        else if (residualType == 0)
                            r()->add(basisToContainerGpAngularVelocityVector[0],rot[B][2]);
                        
                    }
                    else if (dim == 3)
                    {                    
                        for (int d=0;d<nLocalDofAngularVelocity;++d)
                        {   
                            if (residualType == 1) 
                                r()->add(basisToContainerGpAngularVelocityVector[d],-rot[B][d]);  
                            else if (residualType == 0)
                                r()->add(basisToContainerGpAngularVelocityVector[d],rot[B][d]);  
                       
                        }
                    }
                }
            }
        }
    }

    auto end_com = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_com = end_com - start_com;
    time_com += float_com.count();
 
    time = toc("ContactAvoidanceArticulatedType");  
}

/*
    This function reads the json file to get the collision force parameters and 
    then executes the corresponding function. It computes the execution time
*/


class Execution_time
{
    public:
        static int nbr;
        static double total_time;
        static double total_time_pre;
        static double total_time_detection;
        static double total_time_com;
};

int Execution_time::nbr = 0;
double Execution_time::total_time = 0.0;
double  Execution_time::total_time_pre = 0.0;
double  Execution_time::total_time_detection = 0.0;
double  Execution_time::total_time_com = 0.0;

template<typename FluidMechanics>
void
reset_executionTime(FluidMechanics const& t)
{
    Execution_time::nbr = 0;
    Execution_time::total_time = 0.0;
    Execution_time::total_time_pre = 0.0;
    Execution_time::total_time_detection = 0.0;
    Execution_time::total_time_com = 0.0;
}

template<std::size_t residualType, typename FluidMechanics, typename DataType>
void
contactForceModels(FluidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    // Read json
    fs::path path (Environment::expand(soption(_name="fluid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    // Get collision force model and type
    std::string model = jsonCollisionForce["model"].get<std::string>();
    std::string type = jsonCollisionForce["type"].get<std::string>();
    
    double time = 0.;
    double time_pre = 0;
    double time_detection = 0.0;
    double time_com = 0.;

    // Model case : contactAvoidance 
    if (model.compare("contactAvoidance") == 0)
    {

        // Type case : sphericalShapedBody
        if (type.compare("sphericalShapedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceSphericalType<residualType>(t, data, forceParam, time, time_pre, time_detection, time_com);
        }

        // Type case : complexShapedBody
        if (type.compare("complexShapedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceComplexType<residualType>(t, data, forceParam, time, time_pre, time_detection, time_com);
        }
            
        // Type case : articulatedBody
        if (type.compare("articulatedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceArticulatedType<residualType>(t, data, forceParam, time, time_pre, time_detection, time_com);
        }
            
    }

    // One can add other models

    // Execution time
    Execution_time::nbr += 1;        
    Execution_time::total_time += time;
    Execution_time::total_time_pre += time_pre;
    Execution_time::total_time_detection += time_detection;
    Execution_time::total_time_com += time_com;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        std::cout << "Total execution time for " << type << " function after " << Execution_time::nbr << " executions is " << Execution_time::total_time << std::endl;
        std::cout << "Total execution time for pre-process phase : " << Execution_time::total_time_pre << std::endl;
        std::cout << "Total execution time for contact detection : " << Execution_time::total_time_detection << std::endl;
        std::cout << "Total execution time for force comutation : " << Execution_time::total_time_com << std::endl; 
    }
  
    
}