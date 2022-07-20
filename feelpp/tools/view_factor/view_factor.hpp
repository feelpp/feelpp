#include <feel/feel.hpp>
#include <contrib/nanoflann/nanoflann.hpp>
using namespace nanoflann;
#include <limits.h>
#include "KDTreeVectorOfVectorsAdaptor.h"
#include <feel/feeldiscr/createsubmesh.hpp>
using namespace Feel;
#include <random>

// Compute the matrix of view factors F_ij associated to a Feelpp mesh
// F_ij is computed for surface markers i and j
template <typename MeshType>
class view_factor
{
   
    struct Ray 
    { 
        public: 
            //using vec_t = Eigen::Vector3d;
            using vec_t = std::vector<double>;
            Ray(const vec_t &orig, const vec_t &dir) : 
                origin(orig),
                dir(dir) 
            {}     
            vec_t origin, dir;  // ray origin and dir 
    };
    typedef typename MeshType::ptrtype mesh_ptrtype;
    typedef typename MeshType::trace_mesh_ptrtype tr_mesh_ptrtype;
    typedef typename MeshType::face_type face_type;
    typedef std::vector<std::vector<double> > my_vector_of_vectors_t;   
    typedef typename matrix_node<double>::type matrix_node_type;

    typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  KdTree_type;
    public:
        view_factor(mesh_ptrtype mesh,int Nrays)
        {
            if (mesh->dimension()==3)
            {
                // Extract the boundary points from the mesh
                M_mesh = mesh;
                M_Nrays = Nrays;
                M_submesh = createSubmesh(_mesh=mesh,_range=boundaryfaces(mesh),_update=MESH_ADD_ELEMENTS_INFO);   
                my_vector_of_vectors_t  boundary_points; 
                int boundary_pts_index = 0;

                boundary_points.resize(M_submesh->numPoints());
                M_point_indices.resize(M_submesh->numPoints());
                for(auto point : M_submesh->points())
                {   
                    boundary_points[boundary_pts_index].resize(int(M_mesh->dimension()));		

                    boundary_points[boundary_pts_index][0] = point.second.node()[0];
                    boundary_points[boundary_pts_index][1] = point.second.node()[1];
                    boundary_points[boundary_pts_index][2] = point.second.node()[2];
                    
                    M_point_indices[boundary_pts_index] = point.first;

                    boundary_pts_index++;
                    
                }    
                std::cout <<   "number of points" << M_submesh->numPoints() << std::endl;
                //M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 1000 /* max leaf */ );                      
                //M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 500 /* max leaf */ );                      
                //M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 200 /* max leaf */ );                      
                M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 10 /* max leaf */ );                      
                M_view_factors_matrix.resize(M_submesh->markerNames().size(),M_submesh->markerNames().size());
                M_view_factor_row.resize(M_submesh->markerNames().size());
                for(auto &m : M_submesh->markerNames())
                {        
                    //std::cout << m.first << m.second[0] << m.second[1] << m.second[2] << m.second[3] << std::endl;
                    if(m.second[1]==M_submesh->dimension())
                    {    
                        M_markers_string.push_back(m.first);
                        M_markers_int.push_back(m.second[0]);
                    }
                }
            }
            else{
                std::cout << "Only 3d meshes are supported" << std::endl;
            }
            // else if(mesh->dimension()==2)
            // {
            //     M_mesh = mesh;
            //     M_Nrays = Nrays;               
            //     my_vector_of_vectors_t  boundary_points; 
            //     int boundary_pts_index = 0;

            //     boundary_points.resize(M_mesh->numPoints());
            //     M_point_indices.resize(M_mesh->numPoints());
            //     for(auto point : M_mesh->points())
            //     {   
            //         boundary_points[boundary_pts_index].resize(int(M_mesh->dimension()));		

            //         boundary_points[boundary_pts_index][0] = point.second.node()[0];
            //         boundary_points[boundary_pts_index][1] = point.second.node()[1];
            //         boundary_points[boundary_pts_index][2] = point.second.node()[2];
                    
            //         M_point_indices[boundary_pts_index] = point.first;

            //         boundary_pts_index++;
                    
            //     }      
            //     M_kdtree_boundariesptr = new KDTreeVectorOfVectorsAdaptor(M_mesh->dimension() /*dim*/, boundary_points, 10 /* max leaf */ );                      
            //     M_view_factors_matrix.resize(M_mesh->markerNames().size(),M_mesh->markerNames().size());
            //     M_view_factor_row.resize(M_mesh->markerNames().size());
            //     for(auto &m : M_mesh->markerNames())
            //     {                            
            //         if(m.second[1]==M_mesh->dimension())
            //         {    
            //             M_markers_string.push_back(m.first);
            //             M_markers_int.push_back(m.second[0]);
            //         }
            //     }

            // }
        }        

        // Compute random direction, uniformly distributed on the sphere
        std::vector<double> get_random_direction(std::vector<double> &random_direction)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();              

            std::default_random_engine generator1(seed);
            std::default_random_engine generator2(seed1);
            std::uniform_real_distribution<double> xi2(0,1);
            std::uniform_real_distribution<double> xi1(0,1);

            double phi = 2.*M_PI*xi1(generator1);
            double theta = asin(sqrt(xi2(generator2)));
            random_direction[0]=sin(theta)*cos(phi);
            random_direction[1]=sin(theta)*sin(phi);
            random_direction[2]=cos(theta);

            return random_direction;
        }       
        
        // Compute the sum of the areas of three triangles
        double element_area(Eigen::Vector3d const& point,Eigen::Vector3d const& el_p1,Eigen::Vector3d const& el_p2,Eigen::Vector3d const& el_p3)
        {
            double area = 0.;            
            auto v1 = (point-el_p1).cross(point-el_p2);
            auto v2 = (point-el_p2).cross(point-el_p3);
            auto v3 = (point-el_p3).cross(point-el_p1);
            area = sqrt(v1.dot(v1))/2. + sqrt(v2.dot(v2))/2. + sqrt(v3.dot(v3))/2. ;            
            return area;
        }

        // Compare the area of the 2d simplex V1V2V3 (as sum of 3 subtriangles V_iV_jB) and the 2d triangle
        // created by the intersection P of the ray with the plane of  V1V2V3 (as sum of V_iV_jP)
        bool isOnSurface(Eigen::Vector3d &point,Eigen::Vector3d &el_p1,Eigen::Vector3d &el_p2,Eigen::Vector3d &el_p3)
        {
            auto c = (el_p1+el_p2+el_p3)/3.;
  
            auto elem_area = element_area(c, el_p1,el_p2,el_p3); 
            auto area = element_area(point, el_p1,el_p2,el_p3);             
            if ((area-elem_area)<1e-4)
                return true;
            else
                return false;
        }

        // Choose a random point in a triangle
        //Eigen::Vector3d get_random_point(tr_mesh_ptrtype &mesh, node_type const& element_barycenter,matrix_node_type const& element_points)
        Eigen::Vector3d get_random_point(tr_mesh_ptrtype &mesh, matrix_node_type const& element_points)
        {            
            // Choose points in a parallelogram, uniformly
            Eigen::Vector3d p1(3),p2(3),p3(3),v(3),u(3),p(3);
            for(int i=0;i<3;i++)
            {
                p1(i)=column(element_points, 0)[i];
                p2(i)=column(element_points, 1)[i];
                p3(i)=column(element_points, 2)[i];                
            }
            v = p2-p1;
            u = p3-p1;
            while(true)
            {
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();             
                unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();             
                std::default_random_engine generator1(seed),generator2(seed1);
                std::uniform_real_distribution<double> xi1(0,1),xi2(0,1);
                double s = xi1(generator1);
                double t = xi2(generator2);
                // If the point is on the left of the diagonal, keep it, else take the symmetric one
                bool in_triangle = (s + t <= 1);
                if(in_triangle)
                    p = p1 + s * u + t * v;  
                else 
                    p= p1 + (1 - s) * u + (1 - t) * v;

                if (isOnSurface(p,p1,p2,p3))
                    return p;
            }

        }

        // // Visit the kdtree and look for the leaves closer to the intersection point of the ray
        // void VisitNodes(KdTree_type::index_t::NodePtr &node, std::vector<double> &point,Ray const& rayon, double tmax,std::vector<int> &leaf_indices)
        // {               
        //     if(node == nullptr)
        //         return;
        //     /* If "node" is a leaf node, collect the indices of the associated points. */
        //     if ((node->child1 == nullptr) && (node->child2 == nullptr))
        //     {                                
        //         for(int j=node->node_type.lr.left; j<node->node_type.lr.right;j++)        
        //             leaf_indices.push_back(j);
        //         // std::cout << "size leaf " << size(leaf_indices) << std::endl;
        //         return;
        //     }
                
        //     /* Choose the child branch to be taken first */
        //     int idx = node->node_type.sub.divfeat; // The splitting dimension x=0, y=1, z=2, etc.
            
        //     double val  = point[idx]; // The [idx] coordinate of the intersection point "point"    
        //     //double splitValue = (node->node_type.sub.divlow + node->node_type.sub.divhigh) / 2.;

        //     std::vector<KdTree_type::index_t::NodePtr>    bestChild;
        //     KdTree_type::index_t::NodePtr    otherChild;


        //    /* if(val <node->node_type.sub.divlow && abs(node->node_type.sub.divhigh-node->node_type.sub.divlow)<1e-2) // without this, the algorithm has a hard time to find the 
        //     {                                                                                                       // intersection points of some rays whose origin is near the splitting plane
        //         bestChild = {node->child1,node->child2};       
        //         otherChild = nullptr;
        //     }
        //     else if(val >node->node_type.sub.divhigh && abs(node->node_type.sub.divlow-node->node_type.sub.divhigh)<1e-2) // without this, the algorithm has a hard time to find the 
        //     {                                                                                                             // intersection points of some rays whose origin is near the splitting plane
        //         bestChild = {node->child1,node->child2};       
        //         otherChild = nullptr;
        //     }
        //     else */                      

        //     if(val <node->node_type.sub.divlow) // if this is true, choose the left child node
        //     {

        //         bestChild  = {node->child1};
        //         otherChild = node->child2;  
                
        //     }
        //     else if(val >node->node_type.sub.divhigh) // if this is true, choose the right child node
        //     {

        //         bestChild  = {node->child2};
        //         otherChild = node->child1;        

        //     }
        //     else if (val >=node->node_type.sub.divlow && val <=node->node_type.sub.divhigh) // if this is true, choose both the left and right child nodes
        //     {

        //         bestChild = {node->child1,node->child2};       
        //         otherChild = nullptr;

        //     }    

        //     if (std::abs(rayon.dir[idx]) <= 1e-6) {

        //         for(int j=0;j<bestChild.size();j++)
        //             VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices);          
                            
        //     }
        //     else{

        //         double t; 
        //         if(point[idx] < node->node_type.sub.divlow)
        //         {             
        //             t = (node->node_type.sub.divlow - point[idx]) / rayon.dir[idx];
        //         }
        //         else if(point[idx] >node->node_type.sub.divhigh)
        //         {
        //             t = -(node->node_type.sub.divhigh - point[idx]) / rayon.dir[idx];            
        //         }
        //         else
        //         {
        //             auto t_low = -(node->node_type.sub.divlow - point[idx]) / rayon.dir[idx];
        //             auto t_high = -(node->node_type.sub.divhigh - point[idx]) / rayon.dir[idx];
        //             t = t_low*(0. <= t_low && t_low < tmax) + t_high*(0. <= t_high && t_high < tmax);
        //             // std::cout << node->node_type.sub.divlow - point[idx] << std::endl;
        //             // std::cout << node->node_type.sub.divhigh - point[idx] << std::endl;
        //             std::cout << rayon.dir[idx] << std::endl;
        //             std::cout << t_low << " " << t_high << std::endl;            
        //         }
        //         // Test if line segment straddles splitting plane
        //         if (0. <= t && t < tmax) {
        //             // Yes, traverse near side first, then far side
        //             for(int i=0;i<bestChild.size();i++)
        //                 VisitNodes(bestChild[i], point,rayon, t,leaf_indices);            


        //             std::vector<double> first_intersection_point(3);
        //             for(int i=0;i<3;i++)
        //                 first_intersection_point[i]= point[i] + rayon.dir[i] * t;
                    
        //             VisitNodes(otherChild, first_intersection_point, rayon, tmax - t,leaf_indices);             


        //         } 
        //         else {
        //             // No, so just traverse near side           
        //             for(int i=0;i<bestChild.size();i++)
        //                 VisitNodes(bestChild[i], point, rayon, tmax,leaf_indices); 
                    
        //         }
        //     }

        // }
        // Visit the kdtree and look for the leaves closer to the intersection point of the ray
        // FROM "Real time collision detection", Ericson C, p.323
        void VisitNodes(KdTree_type::index_t::NodePtr &node, std::vector<double> &point,Ray const& rayon, double tmax,std::vector<int> &leaf_indices)
        {               
            if(node == nullptr)
                return;
            /* If "node" is a leaf node, collect the indices of the associated points. */
            if ((node->child1 == nullptr) && (node->child2 == nullptr))
            {                                
                for(int j=node->node_type.lr.left; j<node->node_type.lr.right;j++)        
                    leaf_indices.push_back(j);
                return;
            }
                
            /* Choose the child branch to be taken first */
            int idx = node->node_type.sub.divfeat; // The splitting dimension x=0, y=1, z=2
            
            double val  = point[idx]; // The [idx] coordinate of the intersection point "point"                

            std::vector<KdTree_type::index_t::NodePtr>    bestChild;
            KdTree_type::index_t::NodePtr    otherChild;

            int first = (val > node->node_type.sub.divlow) && (val > node->node_type.sub.divhigh);
            if(first) // if this is true, choose the right child node
            {

                bestChild  = {node->child2};
                otherChild = node->child1;  
                
            }
            else
            {
                bestChild  = {node->child1};
                otherChild = {node->child2};  
            }

            if (std::abs(rayon.dir[idx]) <= 1e-6) {

                for(int j=0;j<bestChild.size();j++)
                    VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices);          
                            
            }
             else {
                // Find t value for intersection between segment and rightmost split plane
                // (nanoflann uses two split planes)
                double t = (node->node_type.sub.divhigh - val) / rayon.dir[idx];
                if (t>=tmax) // if the segment cuts the right split plane outside the bounding box, try with the left plane
                {
                    t = (node->node_type.sub.divlow - val) / rayon.dir[idx];
                    bestChild  = {node->child1,node->child2};
                    otherChild = {};
                }
                // Test if line segment intersects splitting plane
                if (0.0 <= t && t < tmax) {
                    // If it does, traverse the near child first, then the far one
                    for(int j=0;j<bestChild.size();j++)
                        VisitNodes(bestChild[j], point, rayon, t,leaf_indices);

                std::vector<double> first_intersection_point(3);
                for(int i=0;i<3;i++)
                        first_intersection_point[i]= point[i] + rayon.dir[i] * t;

                VisitNodes(otherChild, first_intersection_point, rayon, tmax - t,leaf_indices); } 
                else {
                    // If not, just traverse the near child
                for(int j=0;j<bestChild.size();j++)
                    VisitNodes(bestChild[j], point, rayon, tmax,leaf_indices); 
                }
                     
            }

        }

        // Checks intersection for the bounding-box of the whole kdtree
        // If an intersection is detected, the tree itself is inspected
        std::vector<int> get_intersection_indices(const Ray &rayon, std::vector<int> &leaf_indices)
        {
            // FROM "Real time collision detection", Ericson C, p.180
            double tmin = 0.;
            double tmax = INT_MAX*1;
            for (int i = 0; i < 3; i++) {       
                if (std::abs(rayon.dir[i]) < 1e-6) 
                {
                    // Ray is parallel to slab. No hit if origin not within slab
                    if (rayon.origin[i] < M_kdtree_boundariesptr->index->root_bbox[i].low || rayon.origin[i] > M_kdtree_boundariesptr->index->root_bbox[i].high) 
                        return leaf_indices;
                }
                else 
                {        
                    // Compute intersection t value of ray with near and far plane of slab                
                    double t1 = (M_kdtree_boundariesptr->index->root_bbox[i].low - rayon.origin[i]) * 1. / rayon.dir[i];
                    double t2 = (M_kdtree_boundariesptr->index->root_bbox[i].high - rayon.origin[i]) * 1. / rayon.dir[i];
                    // Make t1 be the intersection with the near plane, t2 with the far plane
                    if (t1 > t2) std::swap(t1, t2);
                    // Compute the intersection of slab intersection intervals
                    if (t1 > tmin) tmin = t1;
                    if (t2 < tmax) tmax = t2; 
                    // Exit with no collision as soon as slab intersection becomes empty if (tmin > tmax) return 0;
                    if (tmin > tmax) return leaf_indices;
                }
            }     
            std::vector<double> first_intersection_point(3); // the intersection with the first slab of the bounding box
            for(int i=0;i<3;i++)
            {
                first_intersection_point[i]= rayon.origin[i] + rayon.dir[i] * tmin;
            }            
            // Starting from the root_node, find the leaves where an intersection is still possible                  
            VisitNodes(M_kdtree_boundariesptr->index->root_node, first_intersection_point, rayon, tmax-tmin,leaf_indices); 
            return leaf_indices;
        }
        
        // Verify if the ray intersects the element
        bool check_intersection_with_triangle(std::pair<unsigned int,uint16_type> const& element, Ray &ray)
        {
            Eigen::Vector3d p1(3),p2(3),p3(3),n1(3),w(3),w_(3);
            Eigen::Matrix3d m(3,3);

            auto nodes =  M_submesh->element(element.first).vertices();

            for ( int i = 0; i < nodes.size1(); ++i )
            {p1( i ) = nodes(i,0);
                p2( i ) = nodes(i,1);
                p3( i ) = nodes(i,2);}
            Eigen::Vector3d origin(3),direction(3);
            direction << ray.dir[0],ray.dir[1],ray.dir[2];
            origin << ray.origin[0],ray.origin[1],ray.origin[2];

            
            // // normal vector
            n1 = (p2-p1).cross(p3-p1);
            n1 = n1/n1.norm();

            // intersection point
            w = origin + direction* ( (p1-origin).dot(n1))/(direction.dot(n1));            
            if(w==origin)
                return false;
            m.col(0) = p2-p1;
            m.col(1) = p3-p1;
            m.col(2) = n1;    
            w_ = m.inverse()*(w-p1);            

            return ((w_(0)>=0) && (w_(1)>=0) && (w_(0) +  w_(1)<=1));
            
        }

        // Check intersection with elements of the mesh, by looking at those connected to each point
        std::pair<bool,int> 
        check_intersection_with_elements(Ray &ray, std::vector<int> &leaf_indices, std::vector<int> &point_indices)
        {
            bool intersection_ = false;

            for(int j=0;j<leaf_indices.size();j++)
            {
                auto current_mesh_node = M_kdtree_boundariesptr->index->vAcc[leaf_indices.data()[j]];        
                // find current_mesh_node in mesh
                auto mesh_node = M_submesh->point(point_indices[current_mesh_node]);         

                for(auto &el : mesh_node.elements())
                {                               
                    if (check_intersection_with_triangle(el,ray))
                    {
                        intersection_ =  true;                        
                        auto marker_id = M_submesh->element(el.first).marker().value();                        
                        return std::make_pair(intersection_,marker_id);
                    }             
                }  
            }
            std::cout << "no intersection for Ray " << ray.origin.data()[0] << " "<< ray.origin.data()[1] << " "<< ray.origin.data()[2] 
             << " " << ray.dir.data()[0] <<" " << ray.dir.data()[1] << " " << ray.dir.data()[2] << std::endl;
            return std::make_pair(intersection_,0);
        }


        // Computes the whole view factor matrix (Nmarkers \times Nmarkers), row by row
        // Exploit the relationship Fij = Fji*Ai/Aj to compute only half of the coefficients
        Eigen::MatrixXd computeViewFactors(){
            for(int i=0;i<M_markers_string.size();i++)
            {
                // Surface with marker i launches rays to surface of marker j
                M_view_factors_matrix.row(i)=computeViewFactor(i,M_markers_string[i]);      
            }
            
            return M_view_factors_matrix;
        }

        // Computes the view factors of Fi with all Fj, j>=i
        Eigen::VectorXd computeViewFactor(int marker_index,std::string marker){
        
            M_view_factor_row.setZero();

            std::vector<double> random_direction(3);            
            auto ray_submesh = createSubmesh(_mesh=M_submesh,_range=markedelements(M_submesh,M_markers_string[marker_index]));               
            for(auto const &el: ray_submesh->elements())
            {
                for(int i=0;i<M_Nrays;i++)
                {
                    //auto pt = ray_submesh->point( point_bottom_indices[get_random_index(imax)]).node();   
                    //auto pt = get_random_point(ray_submesh,el.second.barycenter(),el.second.vertices());                 
                    auto pt = get_random_point(ray_submesh,el.second.vertices()); 
                    std::vector<double> random_origin = {pt(0),pt(1),pt(2)};     

                    get_random_direction(random_direction);  

                    Eigen::Vector3d rand_dir(3); 
                    Eigen::Vector3d p1(3),p2(3),p3(3),c(3);
                    for(int i=0;i<3;i++)
                    {
                        p1(i)=column(el.second.vertices(), 0)[i];
                        p2(i)=column(el.second.vertices(), 1)[i];
                        p3(i)=column(el.second.vertices(), 2)[i];
                        c(i) = el.second.barycenter()[i];
                        rand_dir(i) = random_direction[i];
                    }                               
                    auto element_normal = (p2-p1).cross(p3-p1);
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        random_direction={-rand_dir(0),-rand_dir(1),-rand_dir(2)};

                    Ray ray(random_origin,random_direction);
                    std::vector<int> leaf_indices={}; 
                    get_intersection_indices(ray,leaf_indices);                                            
                    
                    auto intersection= check_intersection_with_elements(ray,leaf_indices,M_point_indices);
                    //std::cout << M_view_factor_row.size() <<intersection.second <<  std::endl;
                    auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), intersection.second);
                    //std::cout <<std::distance(M_markers_int.begin(),index_view_factor) <<intersection.second << std::endl;
                    if(intersection.first)
                        M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
                }  
            }              
            M_view_factor_row /=(1.*M_Nrays*ray_submesh->numElements());

            std::cout << M_Nrays*ray_submesh->numElements() << std::endl;
            
            return M_view_factor_row;
        }

        // Computes the view factors of Fi with all Fj, j>=i
        Eigen::VectorXd computeViewFactorPatchMesh(Eigen::Vector3d patchOrigin,Eigen::Vector3d patchNormal){
        
            M_view_factor_row.setZero();

            std::vector<double> random_direction(3);            
            
            for(int i=0;i<M_Nrays;i++)
            {            
                auto pt = patchOrigin;
                std::vector<double> random_origin = {pt(0),pt(1),pt(2)};     

                get_random_direction(random_direction);  

                Eigen::Vector3d rand_dir(3); 
                for(int i=0;i<3;i++)
                {                    
                    rand_dir(i) = random_direction[i];
                }                               
                auto element_normal = patchNormal;
                element_normal.normalize();
                if(rand_dir.dot(element_normal)>0.)
                    random_direction={-rand_dir(0),-rand_dir(1),-rand_dir(2)};

                Ray ray(random_origin,random_direction);
                std::vector<int> leaf_indices={}; 
                get_intersection_indices(ray,leaf_indices);                                            
                
                auto intersection= check_intersection_with_elements(ray,leaf_indices,M_point_indices);
                
                std::cout << "intersection point" << random_direction[0]/random_direction[2]*2 << " " << random_direction[1]/random_direction[2]*2 <<" " <<2 <<std::endl;                
                auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), intersection.second);
                if(intersection.first)
                    M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
            }         
            M_view_factor_row /=(1.*M_Nrays);
            // }
            std::cout << M_Nrays << std::endl;
            
            return M_view_factor_row;
        }

        Eigen::VectorXd computeViewFactorPatchMeshControlledRays(Eigen::Vector3d patchOrigin,Eigen::Vector3d patchNormal){
        
            M_view_factor_row.setZero();

            std::vector<double> random_direction(3);
            std::vector<double> problematic_theta,problematic_phi;
            std::vector<double> theta_vec(M_Nrays),phi_vec(M_Nrays);  
            double theta,phi;          
            for(int j=0;j<M_Nrays;j++)
            {
                theta_vec[j] = asin(sqrt(j*1./M_Nrays*1.));
                for(int i=0;i<M_Nrays;i++)
                {                    
                    phi_vec[i] = 2.*M_PI*i*1./M_Nrays*1.;            
                    theta = theta_vec[j];
                    phi = phi_vec[i];
                    //auto pt = ray_submesh->point( point_bottom_indices[get_random_index(imax)]).node();   
                    auto pt = patchOrigin;
                    std::vector<double> random_origin = {pt(0),pt(1),pt(2)};     
                    random_direction[0]=sin(theta)*cos(phi);
                    random_direction[1]=sin(theta)*sin(phi);
                    random_direction[2]=cos(theta);
                    get_random_direction(random_direction);  

                    Eigen::Vector3d rand_dir(3); 
                    for(int k=0;k<3;k++)
                    {                    
                        rand_dir(k) = random_direction[k];
                    }                               
                    auto element_normal = patchNormal;
                    element_normal.normalize();
                    if(rand_dir.dot(element_normal)>0.)
                        random_direction={-rand_dir(0),-rand_dir(1),-rand_dir(2)};

                    Ray ray(random_origin,random_direction);
                    std::vector<int> leaf_indices={}; 
                    get_intersection_indices(ray,leaf_indices);  
                    
                    auto intersection= check_intersection_with_elements(ray,leaf_indices,M_point_indices);                
                    auto index_view_factor = std::find(M_markers_int.begin(), M_markers_int.end(), intersection.second);
                    if(intersection.first)
                        M_view_factor_row(std::distance(M_markers_int.begin(),index_view_factor))++;
                    if(!intersection.first)
                        {
                            std::cout << "theta " << theta << "phi " << phi << std::endl;
                            problematic_theta.push_back(theta);
                            problematic_phi.push_back(phi);
                        }
                }
            }         
            M_view_factor_row /=(1.*M_Nrays*M_Nrays);            
            // std::ofstream outFile("problematic_angles.csv");
            // // the important part
            // outFile << "Theta, Phi \n" ;
            // for (int i=0;i<problematic_theta.size();i++) outFile << problematic_theta[i] <<","<<problematic_phi[i] << "\n";
            // outFile.close();
            return M_view_factor_row;
        }


        mesh_ptrtype mesh() {return M_mesh;}
        std::vector<std::string> markerNames(){return M_markers_string;}
        
    mesh_ptrtype M_mesh;
    Eigen::MatrixXd M_view_factors_matrix;
    std::vector<std::string> M_markers_string;
    std::vector<int> M_markers_int;
    int M_Nrays;
    tr_mesh_ptrtype M_submesh;
    std::vector<int> M_point_indices;
    KdTree_type* M_kdtree_boundariesptr;
    Eigen::VectorXd M_view_factor_row;

};