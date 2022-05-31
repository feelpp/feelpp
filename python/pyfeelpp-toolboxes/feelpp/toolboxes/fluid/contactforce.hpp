#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 
#include <feel/feelfilters/exporter.hpp>
#include <feel/feells/distancetorange.hpp>
#include <feel/feelcore/ptreetools.hpp>

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

/*
    Implementation of ContactAvoidance model for spherical shaped bodies
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceSphericalType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time)
{
    tic();
    // Define force parameters
    int const dim = t.nDim;
    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;    
    
    // Get the number of bodies, boundary markers, radii and mass centers
    std::vector<RowVectord> massCenters;
    std::vector<std::string> marker_bodies;
    std::vector<double> radii;
    int nbr_bodies = 0;

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

    // Associate remaining boundary to fluid
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    BOOST_FOREACH(auto m, t.mesh()->markerNames())
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

    // Print mass centers, markers and radii
    fmt::print( "Mass centers : {}\n", massCenters);
    fmt::print( "Radii : {}\n", radii);
    fmt::print( "Boundary markers : {}\n", marker_bodies); 
    fmt::print( "Boundary markers fluid : {}\n", marker_fluid); 

    // Define repulsive force 
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(nbr_bodies));
    
    for (int b_i = 0;b_i < nbr_bodies; b_i++)
    {
        // Body-Body collision
        for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
        {
            double dist_ij = sqrt((massCenters[b_j]-massCenters[b_i]).squaredNorm());
            double activation = -(dist_ij-radii[b_i]-radii[b_j]-forceRange);

            if (activation > 0 && dist_ij >= radii[b_i]+radii[b_j])
            {
                RowVectord G_ij = massCenters[b_i]- massCenters[b_j];
                RowVectord F_ij = 1/epsBody * std::pow(activation,2)*G_ij;

                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
                repulsion_forces[b_j] = repulsion_forces[b_j] - F_ij; 
            }
        }
    }
    // Body-Wall collision
    auto Vh = Pch<1>(t.mesh());
    auto exp = exporter(_mesh = t.mesh(),_name="distance");
    exp->addRegions();

    for (int b_w=0;b_w < nbr_fluid;b_w++)
    {
        auto distToBw = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange);
        exp->add("dist_" + marker_fluid[b_w], distToBw);

        for (int b_i=0;b_i<nbr_bodies;b_i++)
        {
            auto minToBi = minmax( _range=markedfaces(t.mesh(),marker_bodies[b_i]), _pset=_Q<2>(), _expr=idv(distToBw));
            if (minToBi.min() <= forceRange)
            {
                std::cout << " Interaction body-wall " << marker_bodies[b_i] << " - " << marker_fluid[b_w] << std::endl;

                RowVectord coord_A = minToBi.argmin();
                
                auto distToBi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange);
                exp->add("dist_" + marker_bodies[b_i], distToBi);
                auto minToBw = minmax( _range=markedfaces(t.mesh(),marker_fluid[b_w]), _pset=_Q<2>(), _expr=idv(distToBi));
                RowVectord coord_B = minToBw.argmin();

                RowVectord G_ij = coord_A - coord_B;
                double activation = forceRange - minToBi.min();
                RowVectord F_ij = 1/epsWall * std::pow(activation,2)*G_ij;

                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
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
                // Print repulsion force
                std::cout << "Repulsion force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                
                if (residualType == 1) 
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],-repulsion_forces[B][d]);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);      
            }
        }
        B++;               
    }
    exp->save();
    time = toc("ContactAvoidanceSphericalType");
}

/*
    Implementation of ContactAvoidance model for spherical shaped bodies with fixed parameters
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
contactforce(FluidMechanics const& t, DataType & data, double &time)
{
    tic();
    // Define dimension
    int const dim = t.nDim;
    int nb_boundaries = 4;
    double radii = 0.125;
    double forceRange = 0.03;
    double epsBody = 100000000000000000;
    double epsWall = 0.00001;
    double height = 4;
    double width = 2;

    //Get the number of bodies, mass centers
    std::vector<RowVectord> massCenters;
    int nbr_bodies = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        massCenters.push_back(bpbc.body().massCenter());
        nbr_bodies ++;
    }

    // Print the centers of mass
    fmt::print( "Mass centers : {}\n", massCenters);

    // Define repulsive force 
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(nbr_bodies));
    
    for (int b_i = 0;b_i < nbr_bodies; b_i++)
    {
        // Body-Wall collision
        std::vector<RowVectord> Boundary;
        RowVectord Top{{massCenters[b_i][0],height+radii}};
        RowVectord Bottom{{massCenters[b_i][0],-radii}};
        RowVectord Left{{-radii,massCenters[b_i][1]}};
        RowVectord Right{{width+radii,massCenters[b_i][1]}};

        Boundary.push_back(Top);
        Boundary.push_back(Bottom);
        Boundary.push_back(Left);
        Boundary.push_back(Right);

        for (int b_j = 0;b_j < nb_boundaries;b_j++)
        {
            double dist_ij = sqrt((Boundary[b_j]-massCenters[b_i]).squaredNorm());
            double activation = -(dist_ij-2*radii-forceRange);

            if (activation > 0 && dist_ij >= 2*radii)
            {
                RowVectord G_ij = (massCenters[b_i]- Boundary[b_j]);
                RowVectord F_ij = 1/epsWall * std::pow(activation,2)*G_ij;
                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
        }
        
        // Body-Body collision
        for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
        {
            double dist_ij = sqrt((massCenters[b_j]-massCenters[b_i]).squaredNorm());
            double activation = -(dist_ij-2*radii-forceRange);

            if (activation > 0 && dist_ij >= 2*radii)
            {
                // Compute repulsion force
                RowVectord G_ij = massCenters[b_i]- massCenters[b_j];
                RowVectord F_ij = 1/epsBody * std::pow(activation,2)*G_ij;
                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
                repulsion_forces[b_j] = repulsion_forces[b_j] - F_ij;
            }
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
                // Print repulsion force
                std::cout << "Repulsion force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                        
                if (residualType == 1) 
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],-repulsion_forces[B][d]);  
                else if (residualType == 0)
                    r()->add(basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]); 
            }
        }
        B++;               
    }
    time = toc("contactforce");
}

/*
    Implementation of ContactAvoidance model for complex shaped bodies
*/

struct Distance_param
{
    int id_A;
    int id_B;
    RowVectord coord_A;
    RowVectord coord_B;
    double dist;        
};

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceComplexType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time)
{
    tic();
    // Get collision parameters    
    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;
    int const dim = t.nDim;

    // Associate boundary marker and mass center to each moving body
    std::vector<RowVectord> massCenters;
    std::vector<std::string> marker_bodies;
    int nbr_bodies = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        marker_bodies.push_back(bpbc.name());
        massCenters.push_back(bpbc.body().massCenter());
        nbr_bodies ++;
    }

    // Associate remaining boundary to fluid
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    BOOST_FOREACH(auto m, t.mesh()->markerNames())
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

    // Print mass centers and boundary markers
    fmt::print( "Mass centers : {}\n", massCenters);
    fmt::print( "Boundary markers : {}\n", marker_bodies); 
    fmt::print( "Boundary markers fluid : {}\n", marker_fluid); 

    // Define export for distance function
    auto Vh = Pch<1>(t.mesh());
    auto exp = exporter(_mesh = t.mesh(),_name="distance");
    exp->addRegions();

    // Define map containing distance parameters for body-body and body-wall interaction
    std::map<std::string,Distance_param> DistanceParam;

    for (int b_i = 0; b_i < nbr_bodies; b_i++)
    {
        // Compute an export distance function from b_i
        auto distToBi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange);
        exp->add("dist_" + marker_bodies[b_i], distToBi);

        // Body-body interaction
        for (int b_j = 0;b_j < nbr_bodies; b_j++)
        {
            if (b_i != b_j)
            {
                auto minToBj = minmax( _range=markedfaces(t.mesh(),marker_bodies[b_j]), _pset=_Q<2>(), _expr=idv(distToBi));
                
                // Verify if contact force is applied
                if (minToBj.min() <= forceRange)
                {
                    if (b_i < b_j)
                    {
                        // Define new key
                        std::string id = marker_bodies[b_i] + "_" + marker_bodies[b_j];

                        // Store distance and coordinates 
                        Distance_param d;
                        d.id_A = b_i;
                        d.id_B = b_j;
                        d.coord_B = minToBj.argmin();
                        d.dist = minToBj.min();

                        // Insert structure in map
                        DistanceParam[id] = d;
                    }
                    else
                    {
                        // Key already exists
                        std::string id = marker_bodies[b_j] + "_" + marker_bodies[b_i];
                        DistanceParam[id].coord_A = minToBj.argmin();
                    }
                }
            }   
        }
        // Body-Wall interaction
        for (int b_w = 0;b_w < nbr_fluid; b_w++)
        {
            auto minToBw = minmax( _range=markedfaces(t.mesh(),marker_fluid[b_w]), _pset=_Q<2>(), _expr=idv(distToBi));
            
            // Verify if contact force is applied
            if (minToBw.min() <= forceRange)
            {    
                // Define key
                std::string id = marker_bodies[b_i] + "_" + marker_fluid[b_w];

                // Store distance and coordinates 
                Distance_param d;
                d.id_A = b_i;
                d.id_B = -1; // Use id -1 for wall boundary
                d.coord_B = minToBw.argmin();
                d.dist = minToBw.min();

                // Compute coord_B
                auto distToBw = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange);
                exp->add("dist_" + marker_fluid[b_w], distToBw);
                auto minToBi = minmax( _range=markedfaces(t.mesh(),marker_bodies[b_i]), _pset=_Q<2>(), _expr=idv(distToBw));
                d.coord_A = minToBi.argmin();

                // Insert structure in map
                DistanceParam[id] = d;
            }
        }
    }
    
    // Print map
    for (const auto& [key, value] : DistanceParam) 
    {
        std::cout << '[' << key << "] : " <<std::endl;
        std::cout << "      id A : " << value.id_A << std::endl;
        std::cout << "      id B : " << value.id_B << std::endl;
        std::cout << "      coordinate A : " << value.coord_A << std::endl;
        std::cout << "      coordinate B : " << value.coord_B << std::endl;
        std::cout << "      Distance : " << value.dist << std::endl;
    }

    // Compute contact force Fi for each body
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));

    // Compute contact points for each body
    std::vector<RowVectord> contact_points(nbr_bodies);
    std::fill(contact_points.begin(), contact_points.end(), RowVectord::Zero(dim));
    
    for (const auto& [key, value] : DistanceParam)
    {
        // Body-body interaction
        if (value.id_B != -1)
        {
            double activation = forceRange - value.dist;
            if (activation >= 0 )
            {
                // Compute repulsion force
                RowVectord G_ij = value.coord_A - value.coord_B;
                RowVectord F_ij = 1/epsBody * std::pow(activation,2)*G_ij;

                repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;
                repulsion_forces[value.id_B] = repulsion_forces[value.id_B] - F_ij; 

                contact_points[value.id_A] = value.coord_A;
                contact_points[value.id_B] = value.coord_B; 
            }
        }
        else // Body-Wall interaction
        {            
            // Compute activation term
            double activation = forceRange - value.dist;

            if (activation > 0)
            {
                // Compute repulsion force
                RowVectord G_ij = value.coord_A - value.coord_B;
                RowVectord F_ij = 1/epsWall * std::pow(activation,2)*G_ij;
                    
                repulsion_forces[value.id_A] = repulsion_forces[value.id_A] + F_ij;
                contact_points[value.id_A] = value.coord_A;  
            }
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
                // Print repulsion force
                std::cout << "Repulsion force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                        
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
                // We define Eigen::Vector3d types to use the function cross()    
                Eigen::Vector3d x_diff_F(0.,0.,0.);
                x_diff_F[0] = contact_points[B][0] - massCenters[B][0];
                x_diff_F[1] = contact_points[B][1] - massCenters[B][1];

                Eigen::Vector3d repulsion_force_B(0.,0.,0.);
                repulsion_force_B[0] = repulsion_forces[B][0];
                repulsion_force_B[1] = repulsion_forces[B][1];

                double rot = - repulsion_force_B.cross(x_diff_F)[2];
                std::cout << "Rotational force : " << rot << " body : " << B << std::endl;
                    
                if (residualType == 1) 
                    r()->add(basisToContainerGpAngularVelocityVector[0],-rot);
                else if (residualType == 0)
                    r()->add(basisToContainerGpAngularVelocityVector[0],rot);
            }        
            else if (dim == 3)
            {
                Eigen::Vector3d x_diff_F = contact_points[B] - massCenters[B]; 
                Eigen::Vector3d repulsion_force_B = repulsion_forces[B];
                Eigen::Vector3d rot = - repulsion_force_B.cross(x_diff_F);

                for (int d=0;d<nLocalDofAngularVelocity;++d)
                {   
                    // Print moment
                    std::cout << "Rotation force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                        
                    if (residualType == 1) 
                        r()->add(basisToContainerGpAngularVelocityVector[d],-rot[d]);  
                    else if (residualType == 0)
                        r()->add(basisToContainerGpAngularVelocityVector[d],rot[d]);  
                       
                }
                
            }     
        }
        B++;
    }
    exp->save();
    time = toc("ContactAvoidanceComplexType");
}

/*
  Implementation of ContactAvoidance for articulated bodies
  Todo : interaction between articulated bodies
*/

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
ContactAvoidanceArticulatedType(FluidMechanics const& t, DataType & data, ns::CollisionForceParam forceParam, double &time)
{
    tic();
    // Define parameters dimension
    int const dim = t.nDim;
    double forceRange = forceParam.forceRange;
    double epsBody = forceParam.epsBody;
    double epsWall = forceParam.epsWall;
  
    
    // Get the number of bodies, boundary markers, mass centers, radii and the gravity force with mass 
    std::vector<RowVectord> massCenters;
    std::vector<std::string> marker_bodies;
    std::vector<RowVectord> GravityForce;
    std::vector<double> radii;
    int nbr_bodies = 0;

    for (auto const&nba : t.bodySetBC().nbodyArticulated())
    {
        for (auto bbc : nba.bodyList())
        {
            std::string marker = bbc->name();
            marker_bodies.push_back(marker);
            massCenters.push_back(bbc->body().massCenter());
         
            double meas = integrate(_range=markedfaces(t.mesh(),marker),_expr=cst(1.0)).evaluate()(0,0);
            if (dim == 2)
              radii.push_back(meas/(2*Pi));
            else if (dim == 3)
              radii.push_back(sqrt(meas/(4*Pi)));

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
 
    // Associate remaining boundary to fluid
    std::vector<std::string> marker_fluid;
    int nbr_fluid = 0;

    BOOST_FOREACH(auto m, t.mesh()->markerNames())
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
 
    // Print mass centers, markers, gravityforce and radii
    fmt::print( "Mass centers : {}\n", massCenters);
    fmt::print( "Radii : {}\n", radii);
    fmt::print( "Boundary markers : {}\n", marker_bodies); 
    fmt::print( "Boundary markers fluid : {}\n", marker_fluid); 
    fmt::print( "Gravity force with mass : {}\n", GravityForce);

    // Define repulsive force and contact points for body-wall interaction
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));
 
    std::vector<RowVectord> contact_points(nbr_bodies);
    std::fill(contact_points.begin(), contact_points.end(), RowVectord::Zero(dim));
 
    auto Vh = Pch<1>(t.mesh());
    auto exp = exporter(_mesh = t.mesh(),_name="distance");
    exp->addRegions();
 
    for (int b_w=0;b_w < nbr_fluid;b_w++)
    {
        auto distToBw = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_fluid[b_w]),_max_distance=2*forceRange);
        exp->add("dist_" + marker_fluid[b_w], distToBw);

        for (int b_i=0;b_i<nbr_bodies;b_i++)
        {
            auto minToBi = minmax( _range=markedfaces(t.mesh(),marker_bodies[b_i]), _pset=_Q<2>(), _expr=idv(distToBw));
            if (minToBi.min() <= forceRange)
            {
                RowVectord coord_A = minToBi.argmin();
                contact_points[b_i] = coord_A;
                
                auto distToBi = distanceToRange(_space=Vh, _range=markedfaces(t.mesh(),marker_bodies[b_i]),_max_distance=2*forceRange);
                exp->add("dist_" + marker_bodies[b_i], distToBi);
                auto minToBw = minmax( _range=markedfaces(t.mesh(),marker_fluid[b_w]), _pset=_Q<2>(), _expr=idv(distToBi));
                RowVectord coord_B = minToBw.argmin();

                RowVectord G_ij = coord_A - coord_B;
                double activation = forceRange - minToBi.min();
                RowVectord F_ij = 1/epsWall * std::pow(activation,2)*G_ij;

                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
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
                    // Print repulsion force
                    std::cout << "Repulsion force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                    
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
            RowVectord G = bpbc.getNBodyArticulated().massCenter();

            size_type startBlockIndexAngularVelocity = t.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
            int nLocalDofAngularVelocity = bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost();

            if (bpbc.spaceAngularVelocity()->nLocalDofWithoutGhost() > 0)
            {
                auto const& basisToContainerGpAngularVelocityVector = r()->map().dofIdToContainerId(rowStartInVector+startBlockIndexAngularVelocity);
                
                for (int B=0;B<nbr_bodies;B++)
                {  
                    if (dim == 2)
                    {
                        
                        // We define Eigen::Vector3d types to use the function cross()
                        
                       
                        Eigen::Vector3d x_diff_F(0.,0.,0.);
                        x_diff_F[0] = contact_points[B][0] - G[0];
                        x_diff_F[1] = contact_points[B][1] - G[1];

                        Eigen::Vector3d x_diff_G(0.,0.,0.);
                        x_diff_G[0] = massCenters[B][0] - G[0];
                        x_diff_G[1] = massCenters[B][1] - G[1];

                        Eigen::Vector3d repulsion_force_B(0.,0.,0.);
                        repulsion_force_B[0] = repulsion_forces[B][0];
                        repulsion_force_B[1] = repulsion_forces[B][1];

                        Eigen::Vector3d gravity_force_B(0.,0.,0.);
                        gravity_force_B[0] = GravityForce[B][0];
                        gravity_force_B[1] = GravityForce[B][1];

                        double rot = - repulsion_force_B.cross(x_diff_F)[2] + gravity_force_B.cross(x_diff_G)[2];
                        std::cout << "Rotational force : " << rot << " body : " << B << std::endl;
                    
                        if (residualType == 1) 
                            r()->add(basisToContainerGpAngularVelocityVector[0],-rot);
                        else if (residualType == 0)
                            r()->add(basisToContainerGpAngularVelocityVector[0],rot);
                        
                    }
                    else if (dim == 3)
                    {
                        Eigen::Vector3d x_diff_F = contact_points[B] - G;
                        Eigen::Vector3d x_diff_G = massCenters[B] - G;
                        Eigen::Vector3d repulsion_force_B = repulsion_forces[B];
                        Eigen::Vector3d gravity_force_B = GravityForce[B];
                 
                        Eigen::Vector3d rot = - repulsion_force_B.cross(x_diff_F) + gravity_force_B.cross(x_diff_G);
                    
                        for (int d=0;d<nLocalDofAngularVelocity;++d)
                        {   
                            // Print moment
                            std::cout << "Rotation force : " << repulsion_forces[B][d] << " dim : " << d << " body : " << B << std::endl;
                        
                            if (residualType == 1) 
                                r()->add(basisToContainerGpAngularVelocityVector[d],-rot[d]);  
                            else if (residualType == 0)
                                r()->add(basisToContainerGpAngularVelocityVector[d],rot[d]);  
                       
                        }
                    }
                }
            }
        }
    }
    exp->save();    
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
};

int Execution_time::nbr = 0;
double Execution_time::total_time = 0.0;

template<typename FluidMechanics>
void
reset_executionTime(FluidMechanics const& t)
{
    Execution_time::nbr = 0;
    Execution_time::total_time = 0.0;
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

    // Model case : contactAvoidance 
    if (model.compare("contactAvoidance") == 0)
    {
        // Type case : sphericalShapedBody
        if (type.compare("FixesphericalShapedBody") == 0)
            contactforce<residualType>(t, data, time);
        
        // Type case : sphericalShapedBody
        if (type.compare("sphericalShapedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceSphericalType<residualType>(t, data, forceParam, time);
        }

        // Type case : complexShapedBody
        if (type.compare("complexShapedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceComplexType<residualType>(t, data, forceParam, time);
        }
            
        // Type case : articulatedBody
        if (type.compare("articulatedBody") == 0)
        {
            ns::CollisionForceParam forceParam = jsonCollisionForce["forceParam"].get<ns::CollisionForceParam>();
            ContactAvoidanceArticulatedType<residualType>(t, data, forceParam, time);
        }
            
    }

    // One can add other models
    Execution_time::nbr += 1;        
    Execution_time::total_time += time;
    std::cout << "Total execution time for " << type << " function after " << Execution_time::nbr << " executions is " << Execution_time::total_time << std::endl;
}