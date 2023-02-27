#include <feel/feelviewfactor/raytracingviewfactor.hpp>

namespace Feel {

// template<typename MeshType>
// Eigen::VectorXd 
// RayTracingViewFactor<MeshType>::get_random_point(matrix_node_type const& element_points)
// {            
//     double dimension;

//     dimension = column(element_points, 0).size();
    
//     // Choose points in a parallelogram, uniformly
//     Eigen::VectorXd p1(dimension),p2(dimension),p3(dimension),v(dimension),u(dimension),p(dimension);
//     for(int i=0;i<3;i++)
//     {
//         p1(i)=column(element_points, 0)[i];
//         p2(i)=column(element_points, 1)[i];
//         p3(i)=column(element_points, 2)[i];                
//     }
//     v = p2-p1;
//     u = p3-p1;
//     while(true)
//     {
//         unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();             
//         unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();             
//         std::default_random_engine generator1(seed),generator2(seed1);
//         std::uniform_real_distribution<double> xi1(0,1),xi2(0,1);
//         double s = xi1(generator1);
//         double t = xi2(generator2);
//         // If the point is on the left of the diagonal, keep it, else take the symmetric one
//         bool in_triangle = (s + t <= 1);
//         if(in_triangle)
//             p = p1 + s * u + t * v;  
//         else 
//             p= p1 + (1 - s) * u + (1 - t) * v;

//         if (isOnSurface(p,p1,p2,p3))
//             return p;
//     }

// }

// template<typename MeshType>
// void 
// RayTracingViewFactor<MeshType>::compute()
// {

// }

}
