#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

typedef Eigen::Matrix<double, 1, Eigen::Dynamic> RowVectord;

template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
contactforce(FluidMechanics const& t, DataType & data, nl::json const& j)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;

    // Define parameters dimension
    int const dim = t.nDim;
    double radius =  j["radius"].get<double>();
    double rho =  j["rho"].get<double>();
    double eps =  j["eps"].get<double>();
    
    // Get the centers of mass and the number of bodies
    std::vector<RowVectord> massCenters;
    int nbr_bodies = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {  
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
        for (int b_j=b_i+1; b_j < nbr_bodies; b_j++)
        {
            // Compute the distance between the centers of mass Gi Gj
            double dist_ij = sqrt((massCenters[b_j]-massCenters[b_i]).squaredNorm());
            
            // Print the distance
            std::cout << "Distance : " << dist_ij << std::endl;

            // Compute activation term
            // double activation = -(dist_ij-radius-radius-rho)/rho;
            double activation = -(dist_ij-radius-radius-rho);

            if (activation > 0 && dist_ij > 2*radius)
            {
                // Compute repulsion force
                //RowVectord G_ij = (massCenters[b_j]- massCenters[b_i])/dist_ij;
                RowVectord G_ij = (massCenters[b_i]- massCenters[b_j]);

                //RowVectord F_ij = c/eps * std::pow(activation,2)*G_ij;
                RowVectord F_ij = 1/eps * std::pow(activation,2)*G_ij;

                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
                repulsion_forces[b_j] = repulsion_forces[b_j] - F_ij; 
            }
            else if (dist_ij <= 2*radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- massCenters[b_j]);
                RowVectord F_ij = 1/eps * (2*radius - dist_ij) * G_ij;
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
                        
                r()->add( basisToContainerGpTranslationalVelocityVector[d],repulsion_forces[B][d]);       
            }
        }
        B++;               
    }   
}



template<std::size_t residualType,typename FluidMechanics,typename DataType>
void
contactforceBoundary(FluidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;

    // Define parameters dimension
    int const dim = t.nDim;
    double radius = 0.125;
    double rho = 0.03;
    double eps = 0.0000001/2;
    
    // Get the centers of mass and the number of bodies
    std::vector<RowVectord> massCenters;
    int nbr_bodies = 0;

    for ( auto const& [bpname,bpbc] : t.bodySetBC() )
    {  
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
            RowVectord Boundary{{massCenters[b_i][0],-radius}};
            // Compute the distance between the centers of mass Gi Gj
            double dist_ij = sqrt((Boundary-massCenters[b_i]).squaredNorm());
            
            // Print the distance
            std::cout << "Distance : " << dist_ij << std::endl;

            // Compute activation term
            // double activation = -(dist_ij-radius-radius-rho)/rho;
            double activation = -(dist_ij-radius-radius-rho);

            if (activation > 0 && dist_ij > 2*radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- Boundary);
                RowVectord F_ij = 1/eps * std::pow(activation,2)*G_ij;
                repulsion_forces[b_i] = repulsion_forces[b_i] + F_ij;
            }
            else if (dist_ij <= 2*radius)
            {
                // Compute repulsion force
                RowVectord G_ij = (massCenters[b_i]- Boundary);
                RowVectord F_ij = 1/eps * (2*radius - dist_ij) * G_ij;
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
