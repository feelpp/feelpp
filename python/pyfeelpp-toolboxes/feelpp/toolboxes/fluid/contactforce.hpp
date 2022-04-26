#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 

using namespace Feel;
using namespace Feel::FeelModels;
using json = nl::json;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

namespace ns {
    struct Collisionforce 
    {
        double radius;
        double rho;
        double eps;
        double eps_;
        double epsW;
        double epsW_;
        double width;
        double height;
    };
        
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Collisionforce,radius,rho,eps,eps_,epsW,epsW_,width,height);
}


template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
contactforce(FluidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;

    // Define dimension
    int const dim = t.nDim;
    int const nb_boundaries = 4; // We consider rectangular fluid domains for the moment
    
    // Read json
    fs::path path (Environment::expand(soption(_name="fluid.filename")));
    json j_Contact;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        j_Contact = j["Contact"]["body"];
    }

    /*  
        Get the number of bodies, the centers of mass and the collision force parameters
    */

    std::vector<RowVectord> massCenters;
    std::vector<ns::Collisionforce> CollisionParam;
    int nbr_bodies = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {
        std::string marker = bpbc.name();
        CollisionParam.push_back(j_Contact[marker].get<ns::Collisionforce>());
        massCenters.push_back(bpbc.body().massCenter());
        nbr_bodies ++;
    }

    // Print the centers of mass
    fmt::print( "Mass centers : {}", massCenters);

    // Define repulsive force 
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(nbr_bodies));
    
    for (int b_i = 0;b_i < nbr_bodies; b_i++)
    {
        // Body-Wall collision
        std::vector<RowVectord> Boundary;
        RowVectord Top{{massCenters[b_i][0],CollisionParam[b_i].height+CollisionParam[b_i].radius}};
        RowVectord Bottom{{massCenters[b_i][0],-CollisionParam[b_i].radius}};
        RowVectord Left{{-CollisionParam[b_i].radius,massCenters[b_i][1]}};
        RowVectord Right{{CollisionParam[b_i].width+CollisionParam[b_i].radius,massCenters[b_i][1]}};

        Boundary.push_back(Top);
        Boundary.push_back(Bottom);
        Boundary.push_back(Left);
        Boundary.push_back(Right);

        for (int b_j = 0;b_j < nb_boundaries;b_j++)
        {
            // Compute the distance between the centers of mass Gi Gj
            double dist_ij = sqrt((Boundary[b_j]-massCenters[b_i]).squaredNorm());
            
            // Compute activation term
            double activation = -(dist_ij-CollisionParam[b_i].radius-CollisionParam[b_i].radius-CollisionParam[b_i].rho);

            if (activation > 0 && dist_ij >= 2*CollisionParam[b_i].radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- Boundary[b_j]);
                RowVectord F_ij = 1/CollisionParam[b_i].epsW * std::pow(activation,2)*G_ij;

                if (residualType == 1)
                    repulsion_forces[b_i] = repulsion_forces[b_i] - F_ij;
                else
                    repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
            else if (dist_ij < 2*CollisionParam[b_i].radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- Boundary[b_j]);
                RowVectord F_ij = 1/CollisionParam[b_i].epsW_ * (2*CollisionParam[b_i].radius - dist_ij) * G_ij;

                if (residualType == 1)
                    repulsion_forces[b_i] = repulsion_forces[b_i] - F_ij;
                else 
                    repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
        }
        

        // Body-Body collision
        for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
        {
            // Compute the distance between the centers of mass Gi Gj
            double dist_ij = sqrt((massCenters[b_j]-massCenters[b_i]).squaredNorm());
            
            // Compute activation term
            double activation = -(dist_ij-CollisionParam[b_i].radius-CollisionParam[b_j].radius-CollisionParam[b_i].rho);

            if (activation > 0 && dist_ij >= CollisionParam[b_i].radius + CollisionParam[b_j].radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- massCenters[b_j]);
                RowVectord F_ij = 1/CollisionParam[b_i].eps * std::pow(activation,2)*G_ij;

                if (residualType == 1)
                {
                    repulsion_forces[b_i] = repulsion_forces[b_i] - F_ij;
                    repulsion_forces[b_j] = repulsion_forces[b_j] + F_ij; 
                }
                else 
                {
                    repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
                    repulsion_forces[b_j] = repulsion_forces[b_j] - F_ij; 
                }
                
            }
            else if (dist_ij < CollisionParam[b_i].radius + CollisionParam[b_j].radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- massCenters[b_j]);
                RowVectord F_ij = 1/CollisionParam[b_i].eps_ * (CollisionParam[b_i].radius + CollisionParam[b_j].radius - dist_ij) * G_ij;

                if (residualType == 1)
                {
                    repulsion_forces[b_i] = repulsion_forces[b_i] - F_ij;
                    repulsion_forces[b_j] = repulsion_forces[b_j] + F_ij;
                }
                else 
                {
                    repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
                    repulsion_forces[b_j] = repulsion_forces[b_j] - F_ij;
                }
                
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
                        
                r()->add( basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);       
            }
        }
        B++;               
    }
}

/*
Contact force for articulated body - bottom interaction
Todo : Interaction with other boundaries
*/
template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
contactforceArticulatedBody(FluidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;

    // Define parameters dimension
    int const dim = t.nDim;

    // Read json
    fs::path path (Environment::expand(soption(_name="fluid.filename")));
    json j_Contact;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        j_Contact = j["Contact"]["body"];
    }
    
    /*  
        Get the number of bodies, the centers of mass, the gravity force with mass 
        and the collision force parameters
    */

    std::vector<RowVectord> massCenters;
    std::vector<RowVectord> GravityForce;
    std::vector<ns::Collisionforce> CollisionParam;
    int nbr_bodies = 0;

    for (auto const&nba : t.bodySetBC().nbodyArticulated())
    {
        for (auto bbc : nba.bodyList())
        {
            std::string marker = bbc->name();
            CollisionParam.push_back(j_Contact[marker].get<ns::Collisionforce>());
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

    // Print the centers of mass
    fmt::print( "Mass centers : {}\n", massCenters);

    // Define repulsive force F_contact
    std::vector<RowVectord> repulsion_forces(nbr_bodies);
    std::fill(repulsion_forces.begin(), repulsion_forces.end(), RowVectord::Zero(dim));

    for (int b_i = 0;b_i < nbr_bodies; b_i++)
    {
        // Body-Wall collision
        RowVectord Bottom{{massCenters[b_i][0],-CollisionParam[b_i].radius}};

        // Compute the distance between the centers of mass Gi Gj
        double dist_ij = sqrt((Bottom-massCenters[b_i]).squaredNorm());
            
        // Compute activation term
        double activation = -(dist_ij-CollisionParam[b_i].radius-CollisionParam[b_i].radius-CollisionParam[b_i].rho);

        if (activation > 0 && dist_ij >= 2*CollisionParam[b_i].radius)
        {
            // Compute repulsion force
            RowVectord G_ij = (massCenters[b_i]- Bottom);
            RowVectord F_ij = 1/CollisionParam[b_i].epsW * std::pow(activation,2)*G_ij;
            repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
        }
        else if (dist_ij < 2*CollisionParam[b_i].radius)
        {
            // Compute repulsion force
            RowVectord G_ij = (massCenters[b_i]- Bottom);
            RowVectord F_ij = 1/CollisionParam[b_i].epsW_ * (2*CollisionParam[b_i].radius - dist_ij) * G_ij;
            repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
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
                        /*
                            We define Eigen::Vector3d types to use the function cross()
                        */
                       
                        Eigen::Vector3d x_diff_F(0.,0.,0.);
                        x_diff_F[0] = massCenters[B][0] - G[0];
                        x_diff_F[1] = massCenters[B][1] - CollisionParam[B].radius - G[1];

                        Eigen::Vector3d x_diff_G(0.,0.,0.);
                        x_diff_G[0] = massCenters[B][0] - G[0];
                        x_diff_G[1] = massCenters[B][1] - G[1];

                        for (int d=0;d<nLocalDofAngularVelocity;++d)
                        {
                            Eigen::Vector3d repulsion_force_B(0.,0.,0.);
                            repulsion_force_B[0] = repulsion_forces[B][0];
                            repulsion_force_B[1] = repulsion_forces[B][1];

                            Eigen::Vector3d gravity_force_B(0.,0.,0.);
                            gravity_force_B[0] = GravityForce[B][0];
                            gravity_force_B[1] = GravityForce[B][1];

                            double rot = - repulsion_force_B.cross(x_diff_F)[2] + gravity_force_B.cross(x_diff_G)[2];
                            std::cout << "Rotational force : " << rot << " body : " << B << std::endl;
                    
                            if (residualType == 1) 
                                r()->add(basisToContainerGpAngularVelocityVector[d],-rot);
                            else if (residualType == 0)
                                r()->add(basisToContainerGpAngularVelocityVector[d],rot);
                        }
                    }
                    else if (dim == 3)
                    {
                        //TODO
                    }
                }
            }
        }
    }    
}
